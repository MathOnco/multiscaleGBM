import java.awt.Color;
import java.awt.Graphics;
/**
 * The dividing Cell.
 */
class Cell {
    //cell states
    boolean mass;//cell is part of tumor mass if 1) infected, or 2) recruited and divided once
    boolean toDie;//cell is marked to die
    boolean vState;//cell is in moving state
    boolean divided;//cell has divided within recording period
    final boolean labeled;//cell is labeled (can be tracked)
    boolean tracked;//cell is currently tracked during recording period
    float ki67;
  
    //cell attributes
    float x, y, x0, y0;//current and previous positions
    float vDiv, xDiv;//proliferation rate and position in cell cycle
    float speed, prevSp;//current speed, previous (maximum, inheritable) speed
    final int pop;//population (0=infected, 1=recruited)
    int prevDiv;//previous (maximum, inheritable) intermitotic time
    int walkTime;//time for persistent walk (or time to stay stopped)
    int angleInDegree;//angle of movement in degrees

    //environment attributes
    boolean activated;//cell is activated if 1) infected, or 2) recruited and in PDGF>Pars.concCutoff
    boolean quiescent;//cell does not have enough space to divide
    float conc;//keeps local concentration
    boolean whiteMatter;//keeps local matter status

    //other
    Color color;//inner cell color
    Color colorO;//outer cell color
  
    /*****************
    * CONSTRUCTOR: 
    ****************/
    Cell(int pop, float x, float y, boolean lab) {
        this.x = x;
        this.y = y;
        this.x0 = x;
        this.y0 = y;

        this.pop = pop;
        this.labeled = lab;
        this.activated= this.pop == 0;
        this.mass= this.pop == 0;
        this.quiescent=false;
        this.divided = false;
        this.tracked = false;
        this.toDie=false;
    }

    /*****************
    * SET CELL
    *****************/

    void init(){
        this.prevDiv = (int) (Functions.skewedGaussian(Pars.divMean, Pars.divError, Pars.divMax, Pars.divMin));//set intermitotic time
        this.xDiv = (World.diceRoller.nextInt(60)+0.f)/60f;//start cell cycle position randomly
        this.vDiv = 1.f/(this.prevDiv+0.f);//set proliferation rate from intermitotic time

        this.vState = true;//set moving=true, not moving=false
        this.speed = Functions.skewedGaussian(Pars.spMean, Pars.spError, Pars.spMax, Pars.spMin);//set speed
        this.prevSp = this.speed;

        this.angleInDegree = World.diceRoller.nextInt(360);//start movment in random direction
        this.walkTime = envRespWalk();//set persistence time
    }

    /*****************
     * DIVISION
     **************/

    void divNewParams(float div, float spe) {
        this.xDiv = 0;//start new cell at top of cell cycle
        this.vState = true;//start new cell moving
        this.mass=true;//start cell as part of tumor mass

        //set daughter traits = parent traits
        float[] tempT=Functions.inheritDirect(div,spe);
        this.prevDiv = (int) (tempT[0]);
        this.prevSp = tempT[1];

        //reset speed, div in accordance to environment, and persistence time
        this.speed = envRespSp();
        this.vDiv = envRespDiv();
        this.walkTime = envRespWalk();
    }

    /*******************************
    * Environmental interactions
    *******************************/

  float envRespDiv(){//find division rate based on PDGF and type
      float fC=this.conc;
      float cRespP = (pop==0) ? (fC+Pars.boost)/(fC+Pars.boost+Pars.halfMaxP) : (fC)/(fC+Pars.recP*Pars.halfMaxP);
      return 1.f*cRespP/(this.prevDiv);
  }

  float envRespSp(){//find speed based on PDGF and type
      float fC=this.conc;
      float cRespM = (pop==0) ? (fC+Pars.boost)/(fC+Pars.boost+Pars.halfMaxM) : (fC)/(fC+Pars.recM*Pars.halfMaxM);
      return Treatment.AMval*cRespM*(this.prevSp);//returns calculated speed adjusted with environment and treatment
  }

  private int envRespWalk(){//find persistence from arrays
      float walkWeight = (this.whiteMatter) ? 3.f : 2.f;//more persistence in white matter
      return (this.vState) ?
              (int) ((Pars.goArray[World.diceRoller.nextInt(Pars.sizeGArray)])*walkWeight):
              Pars.stopArray[World.diceRoller.nextInt(Pars.sizeSArray)];
  }

  void findEnvironment(){//find corresponding white/gray matter and PDGF conc from hex field
      int indX = Functions.getHexCoordinates('x', this.x, this.y, Pars.numPointsW, Pars.numPointsH, Pars.hexSide, Pars.hexDiag);
      int indY = Functions.getHexCoordinates('y', this.x, this.y, Pars.numPointsW, Pars.numPointsH, Pars.hexSide, Pars.hexDiag);
      this.whiteMatter = Field.whiteMatter[indX][indY];
      this.conc = Field.conc[indX][indY];
  }
    

  /*****************
  * UPDATING
  *****************/

  void reset(){//reset vState, angle, and persistence
      this.vState = this.quiescent || (World.diceRoller.nextInt(2) == 0);//pick v state randomly 50% chance each way, quiescent = moving
      //pick angle: if in white use distribution, else random
      int tempAng=(this.whiteMatter) ? (int) (this.angleInDegree+Functions.boundedGaussian(0, Pars.theEr, 180, -180))
              : World.diceRoller.nextInt(360);
      this.angleInDegree = (tempAng>360) ? tempAng-360 : (tempAng<0) ? tempAng+360 : tempAng;
      this.speed = envRespSp();
      this.vDiv = envRespDiv();
      this.walkTime = envRespWalk();
  }

  float[] moveIntoTissue(){
      float[] tempXY=new float[2];
      int k=0;
      int j=0;
      int hX, hY;
      boolean out = true;

      while(out && k<10){//find angle and distance that is still within tissue
          this.angleInDegree=World.diceRoller.nextInt(360);
          tempXY[0] = (float) (this.x + 2*Pars.rad*k*Math.cos(this.angleInDegree*Math.PI/180.f));
          tempXY[1] = (float) (this.y + 2*Pars.rad*k*Math.sin(this.angleInDegree*Math.PI/180.f));
          hX = Functions.getHexCoordinates('x', tempXY[0], tempXY[1], Pars.numPointsW, Pars.numPointsH, Pars.hexSide, Pars.hexDiag);
          hY = Functions.getHexCoordinates('y', tempXY[0], tempXY[1], Pars.numPointsW, Pars.numPointsH, Pars.hexSide, Pars.hexDiag);
          out=Field.noMatter[hX][hY];

          j=j+1;//check new angle
          if(j>10){
              k=k+1;//will increase distance from original pt by a cell diameter
              j=0;//reset angle increment
              if(k>=3){System.out.println("over");}//quit and notify
          }
      }
      return tempXY;
  }

    void move(){//move while avoiding overly dense regions
        float[] tempXY = new float[2];//location
        int mX,mY, hX,hY;//mesh indexes
        int mX0,mY0;//mesh indexes
        if(this.vState){
            mX0 = (int) (Math.floor(this.x/(Data.meshR+0.f)));
            mY0 = (int) (Math.floor(this.y/(Data.meshR+0.f)));
            // temporary move
            tempXY[0] = (float) (this.x + this.speed*Math.cos(this.angleInDegree*Math.PI/180.f)*Pars.frameTime);
            tempXY[1] = (float) (this.y - this.speed*Math.sin(this.angleInDegree*Math.PI/180.f)*Pars.frameTime);

            //check density
            mX = (int) (Math.floor(tempXY[0]/(Data.meshR+0.f)));
            mY = (int) (Math.floor(tempXY[1]/(Data.meshR+0.f)));
            //get hex coords
            hX=Functions.getHexCoordinates('x',tempXY[0], tempXY[1],Pars.numPointsW,Pars.numPointsH,Pars.hexSide,Pars.hexDiag);
            hY=Functions.getHexCoordinates('y',tempXY[0], tempXY[1],Pars.numPointsW,Pars.numPointsH,Pars.hexSide,Pars.hexDiag);

            if(mX<0 || mY<0 || mX>=Data.meshRNx || mY>=Data.meshRNy || Field.noMatter[hX][hY]){
                tempXY=moveIntoTissue();//move cell into tissue
            }
            else{
                //if density of new position is greater than current AND density is greater than carrying capaciity, get new moving angle and don't move
                if(Data.densCells[mX][mY]/(Data.gwCap[mX][mY]+0.f)>Data.densCells[mX0][mY0]/(Data.gwCap[mX0][mY0]+0.f)
                        && Data.densCells[mX][mY]/(Data.gwCap[mX][mY]+0.f)>1.f){
                    this.angleInDegree=World.diceRoller.nextInt(360);
                    tempXY[0] = this.x;
                    tempXY[1] = this.y;

                }
            }
            //set new position
            this.x=tempXY[0];
            this.y=tempXY[1];
        }
        this.walkTime-=Pars.frameTime;//count down persistence time
    }

  void divUpdate(){//update the division timer
      this.xDiv += (!this.quiescent)?this.vDiv*Pars.frameTime:0;
  }

  void getKi67(){//get Ki67, cell is + if within the last 10 hours of the cell cycle
      float tempT=(this.xDiv>(1.f*this.prevDiv-0.5f*Pars.divMin)/this.prevDiv)? 1: 0;
      this.ki67 = tempT;
    }
   
  /***************************
  * GRAPHICS
  ******************************/

  void draw(Graphics g, int trait, boolean noShow, int scale) {//this draws the cells
      if(!noShow){
          Color topColor = color;
          g.setColor(topColor);
          g.fillOval((int)(x - Pars.rad), (int)(y - Pars.rad), (int)(scale * Pars.rad),(int)(scale * Pars.rad));
          g.setColor(colorO);
          if(trait==1 || trait==4){
              g.drawOval((int)(x - Pars.rad), (int)(y - Pars.rad), (int)(scale * Pars.rad),(int)(scale * Pars.rad));
          }
      }
  }

    void drawPart(Graphics g, boolean noShow) {//this draws cells for the Ki67 graphics
        if(!noShow){
            Color topColor = color;
            g.setColor(topColor);
            float xConvert = (x-Pars.histX)*(1.f*Pars.histScale);
            float yConvert = (y-Pars.histY)*(1.f*Pars.histScale);
            g.fillOval((int)(xConvert- Pars.rad), (int)(yConvert- Pars.rad), (int)(Pars.histScale * Pars.rad),(int)(Pars.histScale * Pars.rad));
        }
    }

}