import java.awt.*;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

class World{
    static Random diceRoller = new Random();//random number generator

    //cells
    static ArrayList cells = new ArrayList();//contains Cell objects
    private static int[] startNumCells = new int[3]; //from input
    static ArrayList<Integer> qList = new ArrayList<>();//list of quiescent cells
    private static ArrayList<Integer> aList = new ArrayList<>();//list of activated cells
    private static ArrayList<Integer> dList = new ArrayList<>();//list of dead cells
    private static float midX;//center pt of cell mass
    private static float midY;//center pt of cell mass
    private static int numLab;//number of labeled cells
    private static int numMass;//number of mass cells
    private static int numInf;//number infected cells
    private static int numRec;//number recruited cells
    private static float Ki67;

    //initial timing
    static boolean killIt = false;//flag to kill simulation
    static int datTimeInd = (Pars.txType==0) ? 0 : 2;//index for tracking data collection time
    private static int indexB = 0;//index for tracking diff time points of bulk scale data
    private static int indexI = 0;//index for tracking diff time points of individual scale data
    private static int killTime;//time to kill simulation

    //tracks
    static int numTracks = 200;//maximum number of tracks
	private static boolean tracksOn;//flag to indicate tracks are on
    private static int trackFrame0 = 0;//start time for tracks
    private static int tracksTime;//total time to track cells
	private static ArrayList<Integer> trackableList = new ArrayList<>();//list of possible cells to track
	static ArrayList<Integer> currentlyTracked = new ArrayList<>();//list of cells actually tracked
    static ArrayList tracks = new ArrayList();//contains Track objects

    /*************************
	*
    *      FRAME UPDATE
    *
    *************************/ 

	static void frameUpdate(int frameNum) {
        //record traits
        int tempInd = frameNum/(60*24);
        if( Pars.phenos && ((Pars.txType!=0 && frameNum==Pars.collectTimesBTx[indexI]) ||
                (frameNum==Pars.collectTimesI[indexI] && Pars.txType==0))){
            Data.writeAllData(tempInd,cells);
        }

		//apply treatment
		if(frameNum>=Pars.txSt && frameNum<=Pars.txSt+Pars.txDur){
			Treatment.get(Pars.txType,frameNum,cells);
		}
		if(frameNum>=killTime) {killIt=true;System.out.println("kill now");}

		reset_density(frameNum);//reset density;define qList,aList popR - record Tracks

        //display cell populations
        if(frameNum%(Pars.movTime)==0 && frameNum>0){
            System.out.println("   "+tempInd+"             "+aList.size()+"            "+Data.meanR);
        }

        setupTracksCollection(frameNum);//for tracks and metrics, collect=True
        if(frameNum%(Pars.movTime)==0) finalizeData(frameNum);//find & write rad, pop vals

        cellLoop(0,aList.size());//loop through activated cells to move, divide, secrete, consume

        Functions.fieldDiffusionBounded(Pars.numPointsW, Pars.numPointsH, Field.conc,Pars.hexDiag,
                Pars.Dc, Pars.frameTime, Pars.pDecay,Field.noMatter);

        killCells();
        metricsCollection(frameNum);
	}

    /*********************
     *
     * cell loop
     *
     **********************/
    private static void cellLoop(int start, int end){
        midX=0;
        midY=0;
        for(int i = start; i < end; i++) {
            int aInd = aList.get(i);//loop through activated cells only
            Cell cell = (Cell) cells.get(aInd);
            if (!cell.toDie) {
                midX += cell.x;
                midY += cell.y;
                if (cell.xDiv >= 1 && !cell.quiescent) {
                    if (!cell.toDie) division(cell);
                }
                else {
                    if (cell.walkTime <= 0) {
                        cell.reset();
                    } else {
                        cell.move();
                    }
                    cell.divUpdate();
                }
                Field.update(cell);//update cells' field with sec & cons
            }
        }
        midX=midX/aList.size();
        midY=midY/aList.size();
    }

	/*****************************
	*
    *    Environment Interaction
    *
    ******************************/

    private static void reset_density(int frameNum){
        resetForMetsCollection(frameNum);//timeT,divRate,meanSpeed,stdSpeed,timeM,spList,delR,timeS,numDivs
        Data.clearArrays(frameNum);//reset density meshes
        aList.clear();
        dList.clear();
        qList.clear();
        numInf=0;
        numRec=0;
        numLab=0;
        numMass=0;

        //find mesh points with cells and add to qList; find activated cells and add to aList
        for (int i = 0; i <cells.size(); i++) {
            Cell cell = (Cell) cells.get(i);

            if(cell.toDie){dList.add(i);}//add dead cells to list
            else {
                if (frameNum % 3 == 0 && tracksOn && Pars.tracks) {recordTracks(i, cell, frameNum);}//record tracks
                if(cell.mass && (frameNum % (Pars.movTime)) == 0){
                    collectData(cell);}//collect data from mass occasionally

                //if time to get metrics, initialize
                if((frameNum==Pars.collectTimesI[indexI] && Pars.txType==0) || (frameNum==Pars.collectTimesBTx[indexI] && Pars.txType!=0)){
                    initializeMets(cell,i);
                }

                cell.findEnvironment();//find the environment currently in, and set values to cell
                if (!cell.activated) {//if cell not currently activated, check if should be
                    cell.activated = !(cell.conc < Pars.concCutoff && cell.pop == 1);
                    if (cell.activated) {
                        cell.reset();
                    }
                }
                else{cell.activated = !(cell.conc < Pars.concCutoff && cell.pop == 1);}

                //find numbers of cell types
                if(frameNum%(Pars.movTime)==0) {
                    if (cell.labeled) numLab++;//count labeled cells
                    if (cell.mass) {
                        numMass++;//count cells in mass
                        if (cell.pop == 0) numInf++;//count inf cells
                        if (cell.pop == 1) numRec++;//count rec cells
                    }
                }
                if (cell.activated) {
                    aList.add(i);
                }

                cell.quiescent = false;//reset quiescent state

                if (cell.mass) {
                    cell.getKi67();//defines Ki67 in each cell in mass
                    //populate mesh defines qList, ids, dens's, and total Ki67
                    Data.populateMesh(cell.x, cell.y, i, cell.pop,cell.ki67);
                }
            }
        }
        setQui();//set quiescence from density - from qList
        if(frameNum%(Pars.movTime)==0) {//find the region for max activity for Ki67 measure
            findMaxKi67();
        }
    }

    private static void resetForMetsCollection(int frameNum){
        //reset lists and arrays for recording metrics
        if((frameNum==Pars.collectTimesI[indexI] && Pars.txType==0) || (frameNum==Pars.collectTimesBTx[indexI] && Pars.txType!=0)) {
            trackableList.clear();
            Data.resetMetricArrays();
        }
    }

    private static void setQui(){//check each populated mesh point (in qList) and set to quiescent if >maxGW
        int mX, mY;
        for (Integer integer : qList) {
            mY = (int) (Math.floor(1.f * integer / Data.meshRNx));
            mX = integer - Data.meshRNx * mY;
            if (Data.densCells[mX][mY] > Data.gwCap[mX][mY] * .75) {
                for (int j = 0; j < Data.densCells[mX][mY]; j++) {
                    Cell cell = (Cell) cells.get(Data.idCellsMesh[mX][mY][j]);
                    cell.vState = true;
                    cell.quiescent = Data.densCells[mX][mY] >= Data.gwCap[mX][mY];
                }
            }
        }
    }

   	/******************
   	*
    *    DIVISION
    *
    *******************/

    private static float[] pos = new float[2];
	private static int[] angs = new int[2];
	private static void division(Cell cell){
		pos=findNewPosition(cell);

 		Cell child = new Cell(cell.pop,pos[0],pos[1],cell.labeled);
		cells.add(child);
        child.whiteMatter = cell.whiteMatter;

		//reset cell
		cell.angleInDegree = angs[0];
		cell.divided=true;
		cell.divNewParams(cell.prevDiv,cell.prevSp);

		//set child
		child.angleInDegree = angs[1];
        child.divided=false;
		child.divNewParams(cell.prevDiv,cell.prevSp);
	}

	private static float[] findNewPosition(Cell cell){
	    //define daughter cell with respect to mother cell
        int cA = cell.angleInDegree+90;
        cA = (cA>360) ? cA-360 : (cA<0) ? cA+360 : cA;
        angs[0] = cA;
        angs[1] = (cA+180>360) ? cA+180-360 : (cA+180<0) ? cA+180+360 : cA+180;

        //define position of new cell based on angle
        float angRad = (float) Math.toRadians(angs[1]);
        pos[0]=(float) (cell.x+2*Pars.rad*(Math.cos(angRad)));
        pos[1]=(float) (cell.y+2*Pars.rad*(Math.sin(angRad)));

        //check that cell in not outside of tissue
        int hX=Functions.getHexCoordinates('x',pos[0], pos[1],Pars.numPointsW,Pars.numPointsH,Pars.hexSide,Pars.hexDiag);
        int hY=Functions.getHexCoordinates('y',pos[0], pos[1],Pars.numPointsW,Pars.numPointsH,Pars.hexSide,Pars.hexDiag);
        while(Field.noMatter[hX][hY]){
            //if outside of tissue, try another angle
            cA=diceRoller.nextInt(360);
            angs[0] = cA;
            angs[1] = (cA+180>360) ? cA+180-360 : (cA+180<0) ? cA+180+360 : cA+180;

            //redefine cell position
            angRad = (float) Math.toRadians(angs[1]);
            pos[0]=(float) (cell.x+2*Pars.rad*(Math.cos(angRad)));
            pos[1]=(float) (cell.y+2*Pars.rad*(Math.sin(angRad)));

            //update hex mesh point to check within tissue
            hX=Functions.getHexCoordinates('x',pos[0], pos[1],Pars.numPointsW,Pars.numPointsH,Pars.hexSide,Pars.hexDiag);
            hY=Functions.getHexCoordinates('y',pos[0], pos[1],Pars.numPointsW,Pars.numPointsH,Pars.hexSide,Pars.hexDiag);
        }

        return pos;
    }

    /*********************
     * Death functions
     *********************/

    private static void killCells(){
        Collections.sort(dList);//sort death list
        if(dList.size()>0){
            shiftTrackIndexes();//shift track indexes to account for dying cells
            for (int q = dList.size()-1; q >=0; q--) {
                Cell cellToDie = (Cell) cells.get(dList.get(q));
                cells.remove(cellToDie);
            }
        }
    }

    private static void shiftTrackIndexes(){
        ArrayList<Integer> indRem = new ArrayList<>();
        for (Integer aCurrentlyTracked : currentlyTracked) {
            Track thisTrack = (Track) tracks.get(aCurrentlyTracked);
            thisTrack.shift = 0;//reset shift amount to zero
        }
        //add shift amounts to tracks to account for removed cells
        for (Integer aDList : dList) {
            for (int q = 0; q < currentlyTracked.size(); q++) {
                Track thisTrack = (Track) tracks.get(currentlyTracked.get(q));
                if (thisTrack.ind > aDList) {//shift all indexes greater than index of removed cell down
                    thisTrack.shift++;
                }
                if (thisTrack.ind == aDList) {
                    indRem.add(q);//add cell to remove list
                }
            }
        }
        //now shift track indexes
        for (Integer aCurrentlyTracked : currentlyTracked) {
            Track thisTrack = (Track) tracks.get(aCurrentlyTracked);
            thisTrack.ind = thisTrack.ind - thisTrack.shift;
        }
        if(indRem.size()>0){//remove tracks in descending order
            for (int q = currentlyTracked.size()-1; q >=0; q--) {
                for(int v=indRem.size()-1;v>=0;v--){
                    if(q==indRem.get(v)){
                        currentlyTracked.remove(q);//remove track
                        indRem.remove(v);//remove index from list
                    }
                }
            }
        }
    }

	/******************
	*
    * 		DATA
    *
    *******************/

    private static void initializeMets(Cell cell, int i){
        int dX = (int) (Math.floor(cell.x/(Data.meshR+0.f)));
        int dY = (int) (Math.floor(cell.y/(Data.meshR+0.f)));
        dX=(dX<0)?0:(dX>=Data.meshRNx)?Data.meshRNx-1:dX;
        dY=(dY<0)?0:(dY>=Data.meshRNy)?Data.meshRNy-1:dY;
        cell.tracked=false;
        if(cell.labeled && !cell.quiescent && !cell.toDie && Data.densCells[dX][dY]<=Pars.rimPercent*Data.gwCap[dX][dY]){
            trackableList.add(i);
        }
    }

    private static void setupTracksCollection(int frameNum){
        if((frameNum == tracksTime) && Pars.tracks){
            System.out.println("tracksOn");
            tracksOn=true;//turn on tracks to record & write
            trackFrame0=frameNum;//record time turned on
        }
        if((frameNum==Pars.collectTimesI[indexI] && Pars.txType==0) || (frameNum==Pars.collectTimesBTx[indexI] && Pars.txType!=0)){
            setTracks();//set Tracks for collection if metrics on;
            Pars.collect=true;//set to collect data
        }

    }

    private static void metricsCollection(int frameNum){
        if(Pars.txType==0 && frameNum==Pars.collectTimesB[indexB]){
            float tempIR=(Data.sampleIR[1]==0) ? 0 : Math.round(100*(Data.sampleIR[0]+0.f)/(Data.sampleIR[1]+0.f))/100f;
            float[] tempS = new float[]{Data.meanR,tempIR};
            Functions.writeFloatVectorVert(Pars.outFile+"data/metricsS.txt", tempS);

            if(indexB < 2){ indexB++;}//go to next collection time point
        }
        if(Pars.collect){
            if((Pars.txType==0 && frameNum==Pars.collectTimesI[indexI]+Pars.dataTime[datTimeInd])||
                    Pars.txType!=0 && frameNum==Pars.collectTimesBTx[indexI]+Pars.dataTime[datTimeInd]){
                //at end of collection: turn off and record metrics
                Pars.collect=false;
                tracksOn=false;
                Data.getM(cells);
                Data.writeM(frameNum);
                //update indexes to set for next collection time
                indexI +=((Pars.txType==0 && indexI <2)||(Pars.txType!=0 && indexI <2))?1:0;
                datTimeInd+=(Pars.txType==0 && datTimeInd<1)?1:0;
            }
            else if((frameNum-Pars.collectTimesI[indexI])%30.0f==0){//while on, get tracks and increase tx time
                Data.getDistTrav(cells);
            }
        }
    }

	private static void collectData(Cell cell){//find tumor size at different angles
        float distMidX = midX-cell.x;
        float distMidY = midY-cell.y;
        float thisDist = (float) (Math.sqrt(distMidX*distMidX+distMidY*distMidY));
        int thisAng = (int) (180.f*Math.atan2(distMidY,distMidX)/Math.PI);
        thisAng = (thisAng<0) ? thisAng+360 : (thisAng>360) ? thisAng-360 : thisAng;
        float tempT = (cell.quiescent) ? Pars.divMax/Pars.divConv : 1.f/(cell.envRespDiv()*Pars.divConv);
        Data.findDistStats(thisDist,thisAng,tempT,cell.envRespSp()/Pars.speedConv,cell.prevDiv/Pars.divConv,cell.prevSp/Pars.speedConv);//collect -> Data.popR
        Data.findIR(cell.x, cell.y, cell.pop);//collect infected and recruited sample
    }

    private static void finalizeData(int frameNum){
        float sumA=0;
        for(int ang=0;ang<12;ang++){
            for(int rad=0;rad<Data.meshRN;rad++){
                float cellsPAR = Data.popR[ang][rad]/Data.areaR[ang][rad];//calculated density
                float cellsPR2 = 0.1f*Pars.maxMesh/(Data.meshR*Data.meshR);//10% of max density
                if(cellsPAR>cellsPR2){//find regions with at least 10% of max density
                    Data.maxdR[ang]=(rad+1)*Data.meshR;//find max radial distance from center over 12 angles
                }
                Data.divR[ang][rad]=Data.divR[ang][rad]/Data.popR[ang][rad];
                Data.speR[ang][rad]=Data.speR[ang][rad]/Data.popR[ang][rad];
                Data.divPR[ang][rad]=Data.divPR[ang][rad]/Data.popR[ang][rad];
                Data.spePR[ang][rad]=Data.spePR[ang][rad]/Data.popR[ang][rad];
            }
        }
        for(int ang=0;ang<6;ang++){
            sumA+=Data.maxdR[ang]+Data.maxdR[ang+6];//find diameter at each of 6 angles and sum
        }
        Data.meanR=Math.round(10*sumA/(6.f*Pars.micMToPix))/10000.f;//find micron value of diameter
	    if(frameNum%(Pars.movTime)==0){
            int getTime=frameNum*Pars.frameTime/60;
            int[] tempPops = new int[]{getTime,numInf, numRec, aList.size(), numLab, numMass};
            Functions.writeIntVectorVert(Pars.outFile+"data/popTime.txt",tempPops);
	    	Functions.writeFloat(Pars.outFile+"data/concRad.txt",Data.meanR);
		}
	}

    private static void findMaxKi67(){//find region of maximum activity of Ki67, and calculate %
        int cellNum=0;//total cells in region
        Ki67=0;//total Ki67 in region
        float maxdivR=0;//tracks maximum Ki67 over all regions
        int[] inds = {0,0};//tracks indexes of maximum ki67
        for(int i=0;i<Data.aMX;i++){//loop over regions and find max activity
            for (int j = 0; j < Data.aMY; j++) {
                if(Data.activityMesh[i][j]>maxdivR){
                    maxdivR=1.f*Data.activityMesh[i][j];
                    inds[0]=i;
                    inds[1]=j;
                }
            }
        }
        //redefine histology section based on maximum activity area
        Pars.histX=inds[0]*Pars.histSizeX;
        Pars.histY=inds[1]*Pars.histSizeY;
        for (int i = 0; i <cells.size(); i++) {//find %ki67 within max region
            Cell cell = (Cell) cells.get(i);
            //only count mass cells within the region of interest
            if(cell.mass && cell.x>=Pars.histX && cell.x<=Pars.histX+Pars.histSizeX && cell.y>=Pars.histY && cell.y<=Pars.histY+Pars.histSizeY) {
                Ki67+=cell.ki67;//add Ki67 val in area
                cellNum+=1;//add total cells in area
            }
        }
        Functions.writeFloat(Pars.outFile+"data/Ki.txt",Ki67/cellNum);//write to file
    }

	private static void recordTracks(int i, Cell cell, int fN){//tracks are written to files
		for(int d=0;d<numTracks;d++){
            Track track = (Track) tracks.get(d);
		   	if(track.ind==i){
           		track.addx(cell.x);
           		track.addy(cell.y);
                int[] tempTracks1 = new int[]{fN-trackFrame0,d};
                float[] tempTracks2 = new float[]{cell.x/Pars.micMToPix,cell.y/Pars.micMToPix};
                int tempDiv = (track.div==true)?1:0;
                Functions.writeIntVectorHoriz(Pars.outFile+"tracks/"+cell.pop+"/track"+d+".txt",tempTracks1);
                Functions.writeFloatVectorHoriz(Pars.outFile+"tracks/"+cell.pop+"/track"+d+".txt",tempTracks2);
                Functions.writeInt(Pars.outFile+"tracks/"+cell.pop+"/track"+d+".txt",tempDiv);
                if(tempDiv==1){
                    System.out.println(i);
                }
           	}
        }
	}

	/********************
	*
	* GRAPHICS
	*
	*********************/

    static BufferedImage drawIt(int trait){
		BufferedImage bi = (trait==6)?new BufferedImage(Pars.histScale*Pars.histSizeX,Pars.histScale*Pars.histSizeY, BufferedImage.TYPE_INT_RGB):
                new BufferedImage(Pars.sizeW, Pars.sizeH, BufferedImage.TYPE_INT_RGB);
	  	Graphics2D g1 = bi.createGraphics();

		//draw background
	  	if(trait==6){
	  	    Color col=new Color(255,250,240);
	  	    g1.setColor(col);
	  	    g1.fillRect(0,0,Pars.histScale*Pars.histSizeX,Pars.histScale*Pars.histSizeY);
	  	}
	  	else{
	  	    g1.setColor(Color.black);
	  	    g1.fillRect(0,0,Pars.sizeW,Pars.sizeH);
	  	}


		//draw field
	  	for (int i = 0; i < Pars.numPointsW; i++) { 
    		for(int j = 0; j < Pars.numPointsH; j++) {
    			if(trait==1 || trait==4) Field.drawMatter(g1,i,j);
            	if(trait==2) Field.drawGF(g1,i,j);
    		}
    	}

    	//draw cells
        for (int i = cells.size()-1; i >= 0; i--) {
         	Cell cell = (Cell) cells.get(i);
         	boolean noShow = !(cell.mass || cell.labeled);

         	if(trait==1 || trait ==4){//draw cells with traits
                noShow = !(cell.mass || cell.labeled);
                float newSp = (trait==4) ? cell.prevSp : cell.envRespSp();
                float tempS = 1*Pars.sBins*(newSp-Pars.spMin)/(Pars.spMax-Pars.spMin);
                tempS = (tempS>Pars.sBins) ? Pars.sBins-1 : (tempS<0) ? 0 : (int) tempS;

                float tempT = (trait==4) ? (cell.prevDiv) : (cell.quiescent) ? Pars.divMax : 1.f/cell.envRespDiv();
                float tempD =1*Pars.dBins*(tempT-Pars.divMin)/(Pars.divMax-Pars.divMin);
                tempD =(tempD>Pars.dBins) ? Pars.dBins-1 : (tempD<0) ? 0 : (int) tempD;

         		cell.color=Functions.colorSet2D(tempD,Pars.dBins,0,tempS,Pars.sBins,0,
                    1.0f, 0.5f, 1.1f);
                cell.colorO = Functions.colorSet2D(tempD,Pars.dBins,0,tempS,Pars.sBins,0,
                        1.0f,0.5f,0.5f);
                cell.draw(g1, trait, noShow,Pars.cellScale);
         	}
         	else if(trait==3){//draw cells as infected or recruited
         		cell.color = (cell.pop==0) ? Color.green : Color.red;
                cell.draw(g1, trait, noShow,Pars.cellScale);
         	}
         	if(trait==6){//draw cells according to Ki67+ or not
                if(!noShow && cell.x>=Pars.histX && cell.x<=Pars.histX+Pars.histSizeX && cell.y>=Pars.histY && cell.y<=Pars.histY+Pars.histSizeY){
                    cell.color=Functions.colorGrad(cell.ki67,1,0.5f,153, 204, 255, 102, 51, 0);
                    cell.drawPart(g1, noShow);
                }
            }
        }
	    return bi;
	}

	static BufferedImage drawConc(int scale){//draws the concentration of PDGF
		int thisSX = scale*Data.meshRNx;
		int thisSY = scale*Data.meshRNy;
        BufferedImage bi = new BufferedImage(thisSX, thisSY,BufferedImage.TYPE_INT_RGB);
	  	Graphics2D g1 = bi.createGraphics();
        Functions.drawConc(g1,Data.meshRNx,Data.meshRNy,scale,Data.densCells, Pars.maxMesh);
        return bi;
    }

    /***********************
     *
     *     INITIALIZATION
     *
     ************************/

    static void setPars(int[] input){//set parameter values from the input file
        startNumCells[1] = 100;//number starting recruited labeled
        startNumCells[0] = 100;//number starting infected

        ////////PARAMETERS
        //finds # of cells from density, assuming 80% brain tissue, 1000 converts # to %, maxMesh*meshRNX*meshRNY is total possible cells in domain
        startNumCells[2] = (int) (0.80f*(input[0]/1000.f)*Pars.maxMesh*Data.meshRNx*Data.meshRNy);
        Pars.theEr = input[1];//angle deviation in white matter
        Pars.startCons=input[2];//in ng/mL

        Pars.Dc = (float) ((input[3])*Math.pow(10,-6)*Pars.dCConv);//diffusion coefficient for PDGF
        Pars.pDecay = (input[4]/1000.f)*Pars.consRConv;//decay of PDGF
        Pars.secR=(input[5]/1.f)*Pars.consRConv;//infected secretion rate PDGF
        Pars.consR=(input[6]/100.f)*Pars.secR;//consumption rate

        Pars.boost = input[7]/10.f;//autocrine boost for infected cells
        Pars.halfMaxP=input[8]/1.f;//half max proliferation response to PDGF
        Pars.halfMaxM=input[9]/1.f;//half max migration response to PDGF

        Pars.recP=input[10]/100.f;//adjustment of activation for recruited cells for proliferation
        Pars.recM=input[11]/100.f;//adjustment of activation for recruited cells for migration

        Pars.div = (input[12]);Pars.divMean = Pars.div*Pars.divConv;//division rate (in hours) and converted to frames
        Pars.divEr = input[13];Pars.divError = Pars.divEr*Pars.divConv;//div error and converted
        Pars.sp = input[14];Pars.spMean = Pars.sp*Pars.speedConv;//migration speed (in microns/hr) and converted to pix per frame
        Pars.spEr = input[15];Pars.spError = Pars.spEr*Pars.speedConv;//speed error and converted

        //load data - persistence times (go and stop) and rat brain atlas
        Pars.goArray=Functions.getVector("../source/input/persGo.dat",Pars.sizeGArray);//array in minutes
        Pars.stopArray=Functions.getVector("../source/input/persStop.dat",Pars.sizeSArray);//array in minutes
        Pars.fArray=Functions.getArray("../source/input/whiteMatter.txt", Pars.sizeW, Pars.sizeH);//white matter array
        Pars.f2Array=Functions.getArray("../source/input/noMatter.txt", Pars.sizeW, Pars.sizeH);//not white or gray matter array
    }

    private static void setTracks(){
        currentlyTracked.clear();
        tracks.clear();
        Track.numPops[0]=0;//infected
        Track.numPops[1]=0;//recruited

        Collections.shuffle(trackableList);//randomize the trackable cells and just take the first cells up to Pars.maxTracks
        numTracks=(trackableList.size()>Pars.maxTracks) ? Pars.maxTracks : trackableList.size();
        for(int i=0;i<numTracks;i++){
            currentlyTracked.add(i);
            int index = trackableList.get(i);
            Cell cell =(Cell) cells.get(index);
            cell.x0 = cell.x;//record x at the start of the track
            cell.y0 = cell.y;//record y at the start of the track
            cell.divided=false;
            cell.tracked=true;

            //define new track
            Track track = new Track(index);
            tracks.add(track);
            track.addx(cell.x);//record position
            track.addy(cell.y);
            track.pop=cell.pop;//record cell type, 0=infected, 1=recruited
            Track.numPops[track.pop]+=1;
        }
    }

    static void setTiming(){//set the simulation kill time and the time to collect track info
        killTime=(Pars.txType!=0) ? Pars.txSt+Pars.txDur+60 : Pars.collectTimesB[2]+Pars.dataTime[1]+60;//time to kill simulation
        tracksTime=(Pars.tracks && Pars.txType==0) ? Pars.collectTimesI[1]:(Pars.tracks)? Pars.collectTimesBTx[1]:100*24*60;//time to record tracks
    }

    static void setField(){
        for (int i=0; i<Pars.numPointsW; i++) {
            for (int j=0; j<Pars.numPointsH; j++){
                Field.setField(i,j);//set initial field white/gray matter, PDGF concentration
                Data.setMesh(i,j,Field.whiteMatter[i][j]);//set density mesh based on gray/white
            }
        }
        Data.findGW();//set carrying capacity due to average gray/white matter
    }

    static void setCells(){
        midX=Pars.centerX;
        midY=Pars.centerY;
        //set cell positions - infected
        for (int h = 0; h < startNumCells[0]; h++){
            int[] pts = Functions.random2DNormal(Pars.centerX,Pars.centerY,Pars.disp);
            Cell cell = new Cell(0, pts[0], pts[1],true);
            cells.add(cell);
        }
        //set cell positions - recruited
        for (int h = startNumCells[0]; h<(startNumCells[2]+startNumCells[0]); h++){
            int[] pts = (h<(startNumCells[1]+startNumCells[0])) ?
                    Functions.random2DNormal(Pars.centerX,Pars.centerY,Pars.disp):
                    Functions.random2D(Pars.sizeW,Pars.sizeH,10);
            int indX = Functions.getHexCoordinates('x', pts[0], pts[1], Pars.numPointsW, Pars.numPointsH, Pars.hexSide, Pars.hexDiag);
            int indY = Functions.getHexCoordinates('y', pts[0], pts[1], Pars.numPointsW, Pars.numPointsH, Pars.hexSide, Pars.hexDiag);
            while(Field.noMatter[indX][indY]){
                pts = (h<(startNumCells[1]+startNumCells[0])) ?
                        Functions.random2DNormal(Pars.centerX,Pars.centerY,Pars.disp):
                        Functions.random2D(Pars.sizeW,Pars.sizeH,10);
                indX = Functions.getHexCoordinates('x', pts[0], pts[1], Pars.numPointsW, Pars.numPointsH, Pars.hexSide, Pars.hexDiag);
                indY = Functions.getHexCoordinates('y', pts[0], pts[1], Pars.numPointsW, Pars.numPointsH, Pars.hexSide, Pars.hexDiag);
            }
            boolean bool = h < (startNumCells[1] + startNumCells[0]);//only the starting cells are labeled, the rest start unlabeled and inactivated
            Cell cell = new Cell(1, pts[0], pts[1],bool);
            cells.add(cell);
        }
        //set attributes of cells
        for (int h = cells.size(); --h >= 0; ){
            Cell cell = (Cell) cells.get(h);

            cell.findEnvironment();
            cell.init();
            cell.reset();
            Data.findAngArea();
        }

    }

	
}