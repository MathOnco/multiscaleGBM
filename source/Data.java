import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;

class Data{

	//mesh properties
	private static int micMesh = 100;//mesh size in microns
	static final int meshR = (int) (micMesh*Pars.micMToPix);//mesh size in pixels
	static final int meshRNx = (int) (Pars.sizeW/(0.f+meshR));//number of x mesh points
	static final int meshRNy = (int) (Pars.sizeH/(0.f+meshR));//number of y mesh points
	static int meshRN = (int) (.75f*meshRNx);//mesh size for radial metrics
	static int[][] densCells = new int[meshRNx][meshRNy];//density of all tumor cells
	static int[][][] idCellsMesh = new int[meshRNx][meshRNy][Pars.maxMesh*10];//save id's to easily

	//histology-Ki67 parameters and mesh
	static int aMX = (int) (1.f*Pars.sizeW/Pars.histSizeX);//size of histology section for Ki67
	static int aMY = (int) (1.f*Pars.sizeH/Pars.histSizeY);
	static float[][] activityMesh = new float[aMX][aMY];//to record Ki67 activity
	static float[][] activityDens = new float[aMX][aMY];

	//gray white values
	private static int[][] gwMesh = new int[meshRNx][meshRNy];
	private static int[][] gwNum = new int[meshRNx][meshRNy];
    static int[][] gwCap = new int[meshRNx][meshRNy];

	//measured
    static float meanR;//records mean radius of tumor
    static int[] sampleIR = new int[2];//find infected and recruited cells within sample region
	static float[][] popR = new float[12][meshRN];//population over angles and radii
	static float[][] areaR = new float[12][meshRN];//area of each angle/radii region
	static float[] maxdR = new float[12];//records max radius at each of 12 angles
	static float[][] divR = new float[12][meshRN];//records division rate at different angles and radii
	static float[][] speR = new float[12][meshRN];//records migration rate at different angles and radii
	static float[][] divPR = new float[12][meshRN];//records potential division rate at different angles and radii
	static float[][] spePR = new float[12][meshRN];//records potential migration rate at different angles and radii

	static int[] numDivs = new int[2];//number of divisions for infected, recruited
	static float[] divRate = new float[2];//% divisions per hr for infected, recruited
	static float[] meanSpeed = new float[2];//mean speed for infected, recruited
	static float[] stdSpeed = new float[2];//standard deviation of speed for infected, recruited
	static float[] meanStop = new float[2];//mean stop time for infected, recruited
	static float[][] spList = new float[Pars.MAX_CELLS][2];//list of all speeds
	static float[][] timeS = new float[Pars.MAX_CELLS][2];//list of stop times
	static int[][] divisions = new int[Pars.MAX_CELLS][2];//list of division events

	/******************
	 * setup
	 ******************/

   	static void setMesh(int ix, int iy, boolean env){
   	    //find x and y mesh pts from hex pts
   		int dX = (int) (Math.floor(Field.Fx[ix][iy]/(meshR+0.f)));
		int dY = (int) (Math.floor(Field.Fy[ix][iy]/(meshR+0.f)));
		dX=(dX<0)?0:(dX>=meshRNx)?meshRNx-1:dX;
		dY=(dY<0)?0:(dY>=meshRNy)?meshRNy-1:dY;
   		gwMesh[dX][dY]+=(env)? 2 : 1; //find GW of mesh from hex, white=2, gray=1
   		gwNum[dX][dY]+=1;//calculate hex points within mesh
   	}

   	static void findGW(){
		for(int i=0; i<meshRNx;i++){
			for(int j=0; j<meshRNy;j++){
				if(gwNum[i][j]>0){
					gwMesh[i][j]=Math.round(gwMesh[i][j]/(1.f*gwNum[i][j]));//finds the average environment - gray/white
					//sets the white carrying capacity to be a fraction of the gray
                    gwCap[i][j]=(Data.gwMesh[i][j]==2) ? (int) (Pars.diffGW*Pars.maxMesh) : Pars.maxMesh;
				}
			}
		}
	}

	/******************
	 * INIT and reset
	 ******************/

	static void initialize(){//reset values
   		sampleIR[0]=0;sampleIR[1]=0;

	   	for(int i=0;i<12;i++){
	   		maxdR[i]=0;
	   		for(int j=0;j<meshRN;j++){
	   			popR[i][j]=0;
	   			divR[i][j]=0;
	   			speR[i][j]=0;
				divPR[i][j]=0;
				spePR[i][j]=0;
	   		}
	   	}
	}

	static void resetMetricArrays(){//reset metric arrays
		Arrays.fill(divRate, 0.f);
		Arrays.fill(meanSpeed, 0.f);
		Arrays.fill(stdSpeed, 0.f);
		spList=Functions.clearArray2Float(spList);
		timeS=Functions.clearArray2Float(timeS);
		divisions=Functions.clearArray2Int(divisions);
		Arrays.fill(numDivs, 0);
	}

	static void clearArrays(int frameNum){
		densCells=Functions.clearArray2Int(densCells);
		activityMesh=Functions.clearArray2Float(activityMesh);
		idCellsMesh=Functions.clearArray3Int(idCellsMesh);

		if(frameNum%(Pars.movTime)==0){
			Data.initialize();
		}

	}

	/**************
	* CALCULATIONS
	**************/
	static void populateMesh(float x, float y, int id, int pop, float act){
		//find the corresponding mesh indexes
	    int dX = (int) (Math.floor(x/(meshR+0.f)));
		int dY = (int) (Math.floor(y/(meshR+0.f)));
		dX=(dX<0)?0:(dX>=meshRNx)?meshRNx-1:dX;
		dY=(dY<0)?0:(dY>=meshRNy)?meshRNy-1:dY;

		idCellsMesh[dX][dY][densCells[dX][dY]]=id;//records the id's of cells at each mesh pt
		densCells[dX][dY] += 1;//records total density of cells at each mesh pt

		int dX2 = (int) (Math.floor(x/(Pars.histSizeX+0.f)));//find indexes for Ki67 histology record
		int dY2 = (int) (Math.floor(y/(Pars.histSizeY+0.f)));
		dX2=(dX2<0)?0:(dX2>=meshRNx)?meshRNx-1:dX2;
		dY2=(dY2<0)?0:(dY2>=meshRNy)?meshRNy-1:dY2;
		activityMesh[dX2][dY2]+=act;//sum up activity
		activityDens[dX2][dY2]+=1;//sum up density

        //Add mesh point to quiescent list if not already there
		int q=meshRNx*dY+dX;
		if(!World.qList.contains(q)){
			World.qList.add(q);
		}
	}

	static void findAngArea(){
		for(int i=0;i<12;i++){//12 angles
			for(int j=0;j<meshRN;j++){//radial distances
				areaR[i][j]=(float)(Math.PI*(meshR*(j+1)*meshR*(j+1)-(meshR*j)*(meshR*j))/12);//find areas of each mesh point
			}
		}
	}

	//this finds the proliferation rate and migration speeds within radial distances of the center of the tumor mass
	static void findDistStats(float distR, int ang, float div, float spe, float divP, float speP){
		int dR = (int) (Math.ceil(distR/(meshR+0.f)));//find radial lattice distance
		dR=(dR>=meshRN) ? meshRN-1 : dR;
		int thisA = (int) (Math.floor(ang/30.f));//find angular lattice pt
		popR[thisA][dR] += 1;//populate density
		//sum up trait values
		divR[thisA][dR]+=div;
		speR[thisA][dR]+=spe;
		divPR[thisA][dR]+=divP;
		spePR[thisA][dR]+=speP;
	}

	static void findIR(float x, float y, int pop){
		int dX = (int) (Math.floor(x/(meshR+0.f)));
		int dY = (int) (Math.floor(y/(meshR+0.f)));
		dX=(dX<0)?0:(dX>=meshRNx)?meshRNx-1:dX;
		dY=(dY<0)?0:(dY>=meshRNy)?meshRNy-1:dY;

		if(dX>=111 && dX<=114 && dY>=30 && dY<=33){//only calculate the I/R ratio within a small section of the tumor core.
			sampleIR[pop]+=1;
		}
	}

	static void getDistTrav(ArrayList cells){
		for (Integer aCurrentlyTracked : World.currentlyTracked) {
			Track thisTrack = (Track) World.tracks.get(aCurrentlyTracked);
			Cell cell = (Cell) cells.get(thisTrack.ind);
			float distTr = (float) (Math.sqrt((cell.x - cell.x0) * (cell.x - cell.x0) + (cell.y - cell.y0) * (cell.y - cell.y0)) / Pars.micMToPix);
			if (distTr > 5) { //moving
				thisTrack.adddelr(distTr);//add distance traveled
				thisTrack.addtimeM();//add to moving time
				cell.x0 = cell.x;//reset new position
				cell.y0 = cell.y;
			} else { //staying still
				thisTrack.addtimeS();//add to stopped time
			}
		}
	}

	static void getM(ArrayList cells){
		//loop over tracked cells
		for (Integer aCurrentlyTracked : World.currentlyTracked) {
			Track thisTrack = (Track) World.tracks.get(aCurrentlyTracked);
			spList[aCurrentlyTracked][thisTrack.pop] = (thisTrack.timeM>0) ? thisTrack.delR /thisTrack.timeM : 0;//record mean speed
			timeS[aCurrentlyTracked][thisTrack.pop]=thisTrack.timeS/(Pars.dataTime[World.datTimeInd]/60.f);//record amount of time stopped to total time
			Cell cell =(Cell) cells.get(thisTrack.ind);
			if(cell.divided){
				numDivs[cell.pop] +=1;//add up divisions 0=I or 1=R
				divisions[aCurrentlyTracked][thisTrack.pop]+=1;//record divisions
			}
		}
		//finalize #s to correct format and write to file
		getSpeeds();
		getDivs(Pars.dataTime[World.datTimeInd]);
	}

	private static void getSpeeds(){
		for(int j=0;j<2;j++){
			//set sums to zero
			float sumSp = 0.f;//speeds
			float sumSt = 0.f;//stops
			int tempi = 0;//track index
			for(int i=0;i<World.numTracks;i++){
				sumSp+=(spList[i][j]>0)?spList[i][j]:0;//sum up tracked cells that are moving
				tempi+=(spList[i][j]>0)?1:0;//find number of cells that are moving
				sumSt+=timeS[i][j];//sum up stop times
			}
			meanSpeed[j]=(tempi>0)?sumSp/(tempi+0.f):0;//divide by number moving
			meanStop[j]=sumSt/(Track.numPops[j]+0.f);//divide by total number of cells
			float stdSpI=0;
			for(int i=0;i<World.numTracks;i++){
				stdSpI+=(spList[i][j]>0)?(spList[i][j]-meanSpeed[j])*(spList[i][j]-meanSpeed[j]):0;
			}
			stdSpeed[j] = (tempi>0) ? (float) (Math.sqrt(stdSpI/(tempi+0.f))) : 0.f;
		}
	}

	private static void getDivs(int dTime){
		divRate[0] = (Track.numPops[0]>0) ? 60*100.f*(numDivs[0])/((dTime+.0f)*(Track.numPops[0]+0.f)) : 0;
		divRate[1] = (Track.numPops[1]>0) ? 60*100.f*(numDivs[1])/((dTime+.0f)*(Track.numPops[1]+0.f)) : 0;
	}

	/********************
	 *
	 * OUTPUT
	 *
	 *********************/

	static void writeAllData(int index, ArrayList cells){
		writeDataTraits(index,cells.size(), cells);
		Functions.writeIntMatrix(Pars.outFile+"/data/AllDens"+index+".txt",Data.densCells);

		Functions.writeFloatMatrix(Pars.outFile+"/data/divR"+index+".txt",Data.divR);
		Functions.writeFloatMatrix(Pars.outFile+"/data/speR"+index+".txt",Data.speR);
		Functions.writeFloatMatrix(Pars.outFile+"/data/divPR"+index+".txt",Data.divPR);
		Functions.writeFloatMatrix(Pars.outFile+"/data/spePR"+index+".txt",Data.spePR);
	}

	static void writeM(int frameNum){//finalize metrics and output to file
		float[] tempMets = new float[]{divRate[0],divRate[1],meanSpeed[0],stdSpeed[0],meanSpeed[1],stdSpeed[1],meanStop[0],meanStop[1]};
		int tempFT = (frameNum-Pars.dataTime[World.datTimeInd])/(24*60);
		Functions.writeFloatVectorVert(Pars.outFile+"data/metrics"+tempFT+".txt",tempMets);
	}

	static void writeDataTraits(int ind, int numCells, ArrayList cells){//write single cell data
		try{
			BufferedWriter fout1 = new BufferedWriter(new FileWriter(Pars.outFile+"/data/divspro"+ind+".txt",true));
			BufferedWriter fout2 = new BufferedWriter(new FileWriter(Pars.outFile+"/data/speespro"+ind+".txt",true));
            BufferedWriter fout3 = new BufferedWriter(new FileWriter(Pars.outFile+"/data/divsppro"+ind+".txt",true));
            BufferedWriter fout4 = new BufferedWriter(new FileWriter(Pars.outFile+"/data/speesppro"+ind+".txt",true));

			for(int i=0;i<Pars.MAX_CELLS;i++){
				if(i<numCells){
					Cell cell = (Cell) cells.get(i); 
					if(cell.labeled){
						if(cell.quiescent){
						}
						else if(cell.vDiv>0){
							fout1.write("		"+1.f/(cell.vDiv*Pars.divConv));
                            fout2.write("	"+cell.speed/Pars.speedConv);
						}
						else{
                            fout1.write("		"+Pars.dMax);
                            fout2.write("	"+cell.speed/Pars.speedConv);
						}
						if(!cell.quiescent){
                            fout3.write("		"+cell.prevDiv/(Pars.divConv));
                            fout4.write("	"+cell.prevSp/Pars.speedConv);
						}
					}
				}
			}
            fout1.write("\n");
            fout2.write("\n");
            fout3.write("\n");
            fout4.write("\n");
            fout1.close();
            fout2.close();
            fout3.close();
            fout4.close();
		}
		catch(IOException e){
			System.out.println("There was a problem"+e);
		}
	}


}