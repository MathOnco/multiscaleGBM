import java.awt.Color;
import java.awt.Graphics;

class Field {
    private static final float fSizeX = (float) ((Pars.sizeW+.0)/Pars.numPointsW);//width of hex lattice
    private static final float fSizeY = (float) ((Pars.sizeH+.0)/Pars.numPointsH);//length of hex lattice
    static float[][] Fx = new float[Pars.numPointsW][Pars.numPointsH]; //stores x value for hex positions
    static float[][] Fy = new float[Pars.numPointsW][Pars.numPointsH];  //stores y value for hex positions
    static boolean[][] whiteMatter = new boolean[Pars.numPointsW][Pars.numPointsH]; //stores whether lattice point is white matter
    static boolean[][] noMatter = new boolean[Pars.numPointsW][Pars.numPointsH]; //stores whether lattice point is empty space
    static float[][] conc = new float[Pars.numPointsW][Pars.numPointsH]; //stores local PDGF concentration

	/**************
	* SET FIELD
	*************/
	static void setField(int i, int j){
        //set matter from input
        int i2 = (int)(i* fSizeX); //convert hex index to array values
        int j2 = (int)(j* fSizeY);
        whiteMatter[i][j] = Pars.fArray[i2][j2] == 1; //defines white matter boolean from imported array
        noMatter[i][j] = Pars.f2Array[i2][j2] == 1; //defines empty matter boolean from imported array

        //set initial PDGF = 0 in tissue
        Fx[i][j]=Functions.getXYFromHexInd('x',i,j,Pars.hexSide,Pars.hexDiag);
        Fy[i][j]=Functions.getXYFromHexInd('y',i,j,Pars.hexSide,Pars.hexDiag);
        conc[i][j] = 0.f;

	    //set initial PDGF in wound area
        float bSize =2.0f*(5000.f/Pars.sizeMicm);//radius of PDGF bolus
	    float distEll=(i2-Pars.centerX)*(i2-Pars.centerX)/(Pars.sizeW/2.f)+(j2-Pars.centerY)*(j2-Pars.centerY)/(1.f*Pars.sizeH); //find distance of hex from center point
	    if(distEll<=bSize){
            conc[i][j] = conc[i][j]+Pars.startCons;
   		} 
	}

	/*****************
	* UPDATE FIELD
	*****************/

	static void update(Cell cell){//secrete & consume PDGF
        int i = Functions.getHexCoordinates('x', cell.x, cell.y, Pars.numPointsW, Pars.numPointsH, Pars.hexSide, Pars.hexDiag);
        int j = Functions.getHexCoordinates('y', cell.x, cell.y, Pars.numPointsW, Pars.numPointsH, Pars.hexSide, Pars.hexDiag);

        conc[i][j] += (cell.pop==0) ? Pars.secR : 0;//secrete if infected
        conc[i][j] -= (conc[i][j]>=Pars.consR*Pars.frameTime) ? Pars.consR*Pars.frameTime : conc[i][j];//consume at rate if possible, otherwise consume all
	}

	/**************
	* DRAW
	*************/

	static void drawGF(Graphics g, int i, int j) {//draw PDGF according to concentration
        int tempNum = (int) (conc[i][j]*255.f/100.f);
        tempNum = (tempNum > 255) ? 255 : ((tempNum < 0) ? 0 : tempNum);
		Color fillColor = new Color((160)*tempNum/255,(32)*tempNum/255,(240)*tempNum/255);
        Functions.drawHexField(g, Fx[i][j], Fy[i][j], Pars.hexSide, fillColor);
	}

	static void drawMatter(Graphics g, int i, int j) {//draw white/gray/no tissue in background
        int col = (noMatter[i][j])? 0 : (whiteMatter[i][j]) ? 170 : 110;
		Color fillColor = new Color(col,col,col);
        Functions.drawHexField(g, Fx[i][j], Fy[i][j], Pars.hexSide, fillColor);
	}
		     

}
