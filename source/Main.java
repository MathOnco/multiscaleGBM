import java.io.*;
import javax.imageio.*;
import java.awt.image.*;

public class Main {//test
  private static int[] input = new int[16];

  public static void main(String[] args) {
      for(int i=0;i<=15;i++){
          input[i] = Integer.parseInt(args[i]);
      }
      Pars.tumN=Integer.parseInt(args[16]);
      Pars.outFile = "../gbm"+Pars.tumN+"/";
      Pars.txType = Integer.parseInt(args[17]); // which tx type?
      Pars.movie = Integer.parseInt(args[18]) != 0; //movie on?
      Pars.tracks = Integer.parseInt(args[19]) != 0; //record tracks?
      Pars.phenos = Integer.parseInt(args[20]) != 0; //record phenotypes?

      noDispInit();
      Functions.writeIntVectorVert(Pars.outFile+"data/params.txt",input);

      //set up cell world
      World.setTiming();
      World.setPars(input);
      World.setField();
      World.setCells();

      int frameNum = 0;
      while(World.cells.size()<Pars.MAX_CELLS && !World.killIt){
          World.frameUpdate(frameNum);

          //write graphics files
          if(frameNum%(Pars.movTime)==0 && Pars.movie) {
              int frHour = frameNum/60;

              for(int i=1;i<=6;i++){
                  BufferedImage bi=(i==5) ? World.drawConc(Pars.concScale) : World.drawIt(i);
                  File f = new File(Pars.outFile+"movie/"+i+"/"+frHour+".gif");
                  try {ImageIO.write(bi, "gif", f);}
                  catch (IOException ex) {ex.printStackTrace();}
              }
          }
          frameNum++;
      }
  }

  private static void noDispInit(){
      //print run-specific info
      String strTx=(Pars.txType==0)?"NO TREATMENT":(Pars.txType==1)?"ANTI-PROLIFERATIVE TX":(Pars.txType==2)?"ANTI-MIGRATORY TX":"ANTI-PROLIFERATIVE+ANTI-MIGRATORY TX";
      String strTum=(Pars.tumN==1)?"NODULAR":(Pars.tumN==2)?"INTERMEDIATE":(Pars.tumN==3)?"DIFFUSE":(Pars.tumN==4)?"HETEROGENEOUS":"HOMOGENEOUS";
      System.out.println(".............");
      System.out.println("Initializing simulation with "+strTum+" tumor");
      System.out.println("under "+strTx);
      System.out.println(".............");
      System.out.println("time(days)  #activated cells  diameter(mm)");

      System.setProperty("java.awt.headless", "true");//no head

      //graphics
      new File(Pars.outFile).mkdir();
      new File(Pars.outFile+"data").mkdir();
      if(Pars.movie){
        new File(Pars.outFile+"movie").mkdir();
        new File(Pars.outFile+"movie/1").mkdir();
        new File(Pars.outFile+"movie/2").mkdir();
        new File(Pars.outFile+"movie/3").mkdir();
        new File(Pars.outFile+"movie/4").mkdir();
        new File(Pars.outFile+"movie/5").mkdir();
        new File(Pars.outFile+"movie/6").mkdir();
      }
      if(Pars.tracks){
        new File(Pars.outFile+"tracks").mkdir();
        new File(Pars.outFile+"tracks/0").mkdir();
        new File(Pars.outFile+"tracks/1").mkdir();
      }
  }


}
