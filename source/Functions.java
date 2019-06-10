import java.awt.*;
import java.util.Random;
import java.io.*;
import java.util.*;
public class Functions {

    private static Random diceRoller = new Random();

    /***********************
     * Input
     ***********************/

    static int[] getVector(String filename, int fileLength){
        int[] thisVector = new int[fileLength];
        try {
            Scanner input=new Scanner(new File(filename));

            for (int i=0;i<fileLength;i++){
                thisVector[i] = input.nextInt();
            }

        } catch (IOException ignored) {
        }
        return thisVector;
    }

    static int[][] getArray(String filename, int rows, int cols){
        int[][] thisArray = new int[rows][cols];
        try {
            Scanner input;
            input = new Scanner(new File(filename));

            for (int j=0;j<cols;j++){
                for (int i=0;i<rows;i++){
                    thisArray[i][j] = input.nextInt();
                }
            }
        } catch (IOException ignored) {
        }
        return thisArray;
    }

    /******************
     * Checks
     ******************/

    static boolean checkOutOfBounds(int i, int x, int y, int maxX, int maxY){
        boolean temp = false;
        if(x<=0 || y<=0 || x>=maxX || y>=maxY){
            System.out.println("out of bounds "+i+" "+x+" "+y);
            temp=true;
        }
        return temp;
    }

    /***********************
     * Output
     ***********************/

    static void writeInt(String filename, int val){
        try{
            BufferedWriter fout0 = new BufferedWriter(new FileWriter(filename,true));
            fout0.write(val+"\n");
            fout0.close();
        }
        catch(IOException e){
            System.out.println("There was a problem"+e);
        }
    }

    static void writeFloat(String filename, float val){
        try{
            BufferedWriter fout0 = new BufferedWriter(new FileWriter(filename,true));
            fout0.write(val+"\n");
            fout0.close();
        }
        catch(IOException e){
            System.out.println("There was a problem"+e);
        }
    }

    static void writeIntVectorVert(String strF, int[] vector){
        try{
            BufferedWriter fout0 = new BufferedWriter(new FileWriter(strF,true));
            for (int aVector : vector) {
                fout0.write(aVector + " ");
            }
            fout0.write("\n");
            fout0.close();
        }
        catch(IOException e){
            System.out.println("There was a problem"+e);
        }
    }

    static void writeIntVectorHoriz(String strF, int[] vector){
        try{
            BufferedWriter fout0 = new BufferedWriter(new FileWriter(strF,true));
            for (int aVector : vector) {
                fout0.write(aVector + " ");
            }
            fout0.close();
        }
        catch(IOException e){
            System.out.println("There was a problem"+e);
        }
    }

    static void writeFloatVectorVert(String strF, float[] vector){
        try{
            BufferedWriter fout0 = new BufferedWriter(new FileWriter(strF,true));
            for (float aVector : vector) {
                fout0.write(aVector + " ");
            }
            fout0.write("\n");
            fout0.close();
        }
        catch(IOException e){
            System.out.println("There was a problem"+e);
        }
    }

    static void writeFloatVectorHoriz(String strF, float[] vector){
        try{
            BufferedWriter fout0 = new BufferedWriter(new FileWriter(strF,true));
            for (float aVector : vector) {
                fout0.write(aVector + " ");
            }
            fout0.close();
        }
        catch(IOException e){
            System.out.println("There was a problem"+e);
        }
    }

    static void writeIntMatrix(String strF, int[][] matrix){
        try{
            BufferedWriter fout0 = new BufferedWriter(new FileWriter(strF,true));
            for (int row=0;row<matrix.length;row++) {
                for (int col=0;col<matrix[0].length;col++){
                    fout0.write(matrix[row][col] + " ");
                }
                fout0.write("\n");
            }
            fout0.close();
        }
        catch(IOException e){
            System.out.println("There was a problem"+e);
        }
    }

    static void writeFloatMatrix(String strF, float[][] matrix){
        try{
            BufferedWriter fout0 = new BufferedWriter(new FileWriter(strF,true));
            for (int row=0;row<matrix.length;row++) {
                for (int col=0;col<matrix[0].length;col++){
                    fout0.write(matrix[row][col] + " ");
                }
                fout0.write("\n");
            }
            fout0.close();
        }
        catch(IOException e){
            System.out.println("There was a problem"+e);
        }
    }

	/***********************
 	*distributions
	***********************/
	static float boundedGaussian(float mean, float dev, float max, float min){
		float  val;
		float gauss = (float) diceRoller.nextGaussian();
		val = dev*gauss+mean;
		while(val>max || val<min){
			gauss = (float) diceRoller.nextGaussian();
			val = dev*gauss+mean;
		}
		return val;
	}

	static int[] random2DNormal(int centerX, int centerY, int dev){
		int[] rander = {0,0};
		float gauss = (float) diceRoller.nextGaussian();
		rander[0] = (int)(dev*gauss+centerX);
		gauss = (float) diceRoller.nextGaussian();
		rander[1] = (int)(dev*gauss+centerY);

		return rander;
	}

	static int[] random2D(int width, int height, int cushion){
		int[] rander = {0,0};
		rander[0] = cushion+diceRoller.nextInt(width-2*cushion);
		rander[1] = cushion+diceRoller.nextInt(height-2*cushion);

		return rander;
	}

	static float skewedGaussian(float mean, float dev, float max, float min){
		float gauss = (float) diceRoller.nextGaussian();
		int sign = (gauss<0) ? -1 : 1;
		float val = (float)(mean+sign*Math.pow(gauss,2)*dev);
		while(val>max || val<min){
			gauss = (float) diceRoller.nextGaussian();
			sign = (gauss<0) ? -1 : 1;
			val = (float)(mean+sign*Math.pow(gauss,2)*dev);
		}
		return val;
	}

	/***********************
 	*coordinate conversions
	***********************/
	static int getHexCoordinates(char dim, float xPos, float yPos, int width, int height, float hexSide, float hexDiag){
		int ind[]={0,0};
		float tempy = yPos/(hexDiag);
		ind[1] = ((tempy-Math.floor(tempy))*(tempy-Math.floor(tempy))>(tempy-Math.ceil(tempy))*(tempy-Math.ceil(tempy)))
				? (int) (Math.ceil(tempy)) : (int) (Math.floor(tempy));
		ind[1]=(ind[1]>=height) ? height-1 : (ind[1]<0) ? 0 : ind[1];

		float tempx = (ind[1]%2==0) ? (int) ((xPos)/(3*hexSide)) :
			(int)((xPos-1.5f*hexSide)/(3*(hexSide)));
		ind[0] = Math.round(tempx);
		ind[0]=(ind[0]>=width) ? width-1 : (ind[0]<0) ? 0 : ind[0];

		if(dim=='x'){return ind[0];}
		else{return ind[1];}
	}

	static float getXYFromHexInd(char dim, int ind1, int ind2, float hexSide, float hexDiag){
		float Fx = (ind2%2==0) ? 3*ind1*hexSide:(float) (hexSide*(1.5+3*ind1));
		float Fy = hexDiag*ind2;

		if(dim=='x'){return Fx;}
		else{return Fy;}
	}

	/***********************
 	*array operations
	***********************/

    static int[][][] clearArray3Int(int[][][] array){
        for (int[][] innerRow: array)
        {
            for (int[] innerInnerRow: innerRow)
            {
                Arrays.fill(innerInnerRow, 0);
            }
        }
        return array;
    }

    static int[][] clearArray2Int(int[][] array){
        for (int[] row : array)
            Arrays.fill(row, 0);
        return array;
    }

    static float[][] clearArray2Float(float[][] array){
        for (float[] row : array)
            Arrays.fill(row, 0);
        return array;
    }

	/***********************
 	*inheritance
	***********************/
	static float[] inheritDirect(float prev1, float prev2){
		float[] vals = {0,0};
		vals[0] = prev1;
		vals[1] = prev2;

		return vals;
	}

	/***********************
 	*diffusion
	***********************/

    static void fieldDiffusionBounded(int hexW, int hexH, float[][] conc, float hexDiag,
                                   float diffC, int frameTime, float decay, boolean[][] noMatter){
        for(int k=0; k<hexW;k++){
		   	for(int j=0;j<hexH;j++){
                float sumTemp = 0;
                if(Field.Fx[k][j]>0. && Field.Fy[k][j]>0. && Field.Fx[k][j]<Pars.sizeW && Field.Fy[k][j]<Pars.sizeH && !noMatter[k][j]
                        && k>=2 && j>=2 && k<hexW-2 && j<hexH-2){
                    float tempConcCheck= (j%2==0) ? conc[k][j]+conc[k-1][j-1]+conc[k][j-2]+conc[k][j-1]+conc[k-1][j+1]+conc[k][j+2]+conc[k][j+1]
                            : conc[k][j]+conc[k][j-1]+conc[k][j-2]+conc[k+1][j-1]+conc[k][j+1]+conc[k][j+2]+conc[k+1][j+1];
                    if(tempConcCheck>=.00000001){
                        Integer list1AX[][] = {{k-1,k, k},{k-1,k, k}};
                        Integer listAY[][] = {{j-1,j-2,j-1},{j+1,j+2,j+1}};
                        Integer list2AX[][] = {{k, k, k+1},{k, k, k+1}};
                        for(int q=0;q<2;q++){
                            for(int p=0;p<3;p++){
                                if(j%2==0){
                                    if(!noMatter[list1AX[q][p]][listAY[q][p]]){
                                        sumTemp += conc[list1AX[q][p]][listAY[q][p]]-conc[k][j];
                                    }
                                }
                                else{
                                    if(!noMatter[list2AX[q][p]][listAY[q][p]]){
                                        sumTemp += conc[list2AX[q][p]][listAY[q][p]]-conc[k][j];
                                    }
                                }
                            }
                        }
                        float D2field = sumTemp/(4.f*hexDiag*hexDiag);
                        conc[k][j] = (conc[k][j]+diffC*D2field*frameTime);
                        conc[k][j] = (conc[k][j]<decay*frameTime) ? 0 : conc[k][j]-decay*frameTime;
                    }
                }
            }
        }
    }

    /*****************
     * draw
     *****************/

    static void drawHexField(Graphics g, float x, float y, float hexSide, Color col) {
        Polygon p = new Polygon();
        for (int i = 0; i <= 5; i++){
            p.addPoint((int) (x+hexSide*Math.cos((Math.PI/180)*(60*i))),
                    (int) (y-hexSide*Math.sin((Math.PI/180)*(60*i))));
        }
        g.setColor(col);
        g.fillPolygon(p);
    }

    static void drawConc(Graphics g1,int meshNx, int meshNy, int binSize, int[][] densCells, int carryCap){

        for (int i = 0; i < meshNx; i++) {
            for(int j = 0; j < meshNy; j++) {
                float dC = densCells[i][j];
                dC=(dC>=carryCap)?1.f:1.f*dC/carryCap;
                g1.setColor(new Color(dC,dC,dC));
                g1.fillRect(i*binSize,j*binSize,binSize,binSize);
            }
        }
    }

    /*****************
     * colors
     *****************/

    static int[] angles={20,100,180,210,270,320};
    static Color colorSet2D(float current1, float max1, float min1, float current2, float max2, float min2,
                              float adjustP, float adjustM, float adjustDark){

        Color color;
        int angSect;
        float thetaCol;
        float p1;

        float xCol = 2*(current1-(adjustP*max1-min1)/2-min1)/(adjustP*max1-min1);
        float current2a = (current2>adjustM*max2) ? adjustM*max2 : current2;
        float yCol = 2*(current2a-adjustM*max2/2-min2)/(adjustM*max2-min2);

        float colDist = (float) (adjustDark*Math.sqrt(xCol*xCol+yCol*yCol));
        colDist=(colDist>1) ? 1.f : colDist;
        thetaCol=(float) ((180.f/Math.PI)*Math.atan2(yCol,xCol));
        thetaCol=(thetaCol>360) ? thetaCol-360.f : (thetaCol<0) ? thetaCol+360.f: thetaCol;

        int[] mm = angles;
        int[] pp = {360-mm[5]+mm[0],mm[1]-mm[0],mm[2]-mm[1],mm[3]-mm[2],mm[4]-mm[3],mm[5]-mm[4]};
        angSect = (colDist<0.0) ? 0 : (thetaCol<=mm[0] && thetaCol>mm[5]) ?  1 :
                (thetaCol>mm[0] && thetaCol<=mm[1]) ? 2 :
                        (thetaCol>mm[1] && thetaCol<=mm[2]) ? 3 :
                                (thetaCol>mm[2] && thetaCol<=mm[3]) ? 4 :
                                        (thetaCol>mm[3] && thetaCol<=mm[4]) ? 5 :
                                                (thetaCol>mm[4] && thetaCol<mm[5]) ? 6 : 1;
        switch (angSect){
            case 0: {
                color = new Color(1-colDist/.5f,1-colDist/.5f,1-colDist/.5f);
                break;}
            case 1: {//cyan->green
                p1 = (thetaCol<=mm[0]) ? (360-mm[5]+thetaCol)/(pp[0]) : (thetaCol-mm[5])/(pp[0]);
                color = new Color(0,colDist,colDist*(1-p1));
                break;}
            case 2:{//green->yellow
                p1 = (thetaCol-mm[0])/pp[1];
                color = new Color(colDist*(p1),colDist,0);
                break;
            }
            case 3:{//yellow->red
                p1 = (thetaCol-mm[1])/pp[2];
                color = new Color(colDist,colDist*(1-p1),0);
                break;
            }
            case 4:{//red->magenta
                p1 = (thetaCol-mm[2])/pp[3];
                color = new Color(colDist,0,colDist*(p1));
                break;
            }
            case 5:{//magenta->blue
                p1 = (thetaCol-mm[3])/pp[4];
                color = new Color(colDist*(1-p1),0,colDist);
                break;
            }
            case 6:{//blue->cyan
                p1 = (thetaCol-mm[4])/pp[5];
                color = new Color(0,colDist*(p1),colDist);
                break;
            }
            default:{
                color = new Color(0,0,0);
            }
        }
        return color;
    }

    static Color colorGrad(float val, float max, float min, int r1, int b1, int g1, int r2, int b2, int g2){
        if (val <= min) {
            Color col=new Color(r1,b1,g1);
            return col;
        }
        if (val >= max) {
            Color col=new Color(r2,b2,g2);
            return col;
        }
        val = (val - min) / (max - min);
        int rnew = (int) (r1 + val*(r2-r1));
        int bnew = (int) (b1 + val*(b2-b1));
        int gnew = (int) (g1 + val*(g2-g1));

        Color col=new Color(rnew,bnew,gnew);
        return col;
    }


}
