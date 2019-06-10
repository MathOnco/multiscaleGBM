import java.util.ArrayList;
class Treatment{
    static float AMval = 1f;//this multiplier to migration speed is adjusted due to anti-migratory Tx

    static void get(int txType, int frameNum, ArrayList cells){
      if(txType==1){//anti-proliferative drug
        antiPro(cells, 60, false);//kill proliferative cells faster than 60h IMT
      }
      else if(txType==2){//anti-migratory drug
          if(frameNum==Pars.txSt) {
              cellSlow();//slow migration of cells
          }
      }
      else if(txType==3){
          if(frameNum==Pars.txSt) {//anti-proliferative and anti-migratory drug combination
              cellSlow();
          }
          antiPro(cells, 60, false);
      }
  }

  /***********************
  * Treatments
  ***********************/
  private static void antiPro(ArrayList cells, int killThres, boolean envR){
      for(int i = cells.size()-1; i >= 0; i--){
          Cell cell = (Cell) cells.get(i);
          float tempT = (!envR) ? cell.prevDiv/(1.f*Pars.divConv) : 1.f/(Pars.divConv*cell.envRespDiv());
          if(tempT<=killThres && !cell.quiescent && cell.mass){//check that cell is below threshold, not quiescent, and part of the mass
              cell.toDie=true;
          }
      }
  }

  private static void cellSlow(){ AMval=.10f;}//multiplier for migration set to 10% of original


}
