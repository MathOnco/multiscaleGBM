import java.util.ArrayList;

class Track{

  ArrayList<Float> x = new ArrayList<>();
  ArrayList<Float> y = new ArrayList<>();
  float delR;
  float timeM;
  float timeS;
  int sampleT=30;//sample time in minutes
  int ind;
  int pop;
  boolean div;
  int shift;
  static int[] numPops=new int[2];

  Track(int i) {
      delR =0;
      timeM=0;
      timeS=0;
      ind=i;
      div=false;
  }

  void addx(float xs){
    this.x.add(xs);
  }
  void addy(float ys){
   	this.y.add(ys);
  }
  void adddelr(float dr) {this.delR =this.delR +dr;}
  void addtimeM() {this.timeM=this.timeM+sampleT/60.f;}//add migration time in hours
  void addtimeS() {this.timeS=this.timeS+sampleT/60.f;}//add stop time in hours

}