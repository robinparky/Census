package scripts.histogram;
////********************************************************************
// Histogram.java 
//
//********************************************************************
//   A simple histogram class. The setData(double f) finds in which bin
// the value falls for nBins between the given minimum and maximum
// values. An integer array keeps track of the number of times the input
// value fell into a particular bin.
//    The DataChart class is used to display the histogram. This class
// in turn uses the BarChart tool which needs the data passed as a 
// string array. This is not optimal but OK for demonstration purposes.
//
//********************************************************************
import java.util.*;	

//********************************************************************
public class Histogram {

  int [] freqArr = null;
  int nBins;
  double xLow,xHigh;
  double delBin;
 
  int overFlows=0,underFlows=0;

  String dataString=null;


  public static void main(String[] args)
  {
    Histogram h = new Histogram(20, -4, 6);

    double[] arr = {-2.0, -3.4, 0.2, -4.5, 2.3,4.6, 3.2, 3.1,8.6,5.7,9.1,9.7, 11, 12};

    for(double d:arr)
    {
	h.setData(d);
    }
/*
    int[] freqArr = h.getFreqArr();
    double[] bins = h.getBins();


    System.out.println("0\t" + h.getUnderFlows());
    for(int i=0;i<freqArr.length;i++)
	System.out.println(bins[i] + "\t" + freqArr[i]);

    System.out.println((bins[bins.length-1] + h.getDelBin()) + "\t" + h.getOverFlows());
    System.out.println(h.getDelBin());
*/
    h.print();
  }

  //----------------------------------------------------------------
  public Histogram (int nBins, double xLow, double xHigh){
  
   this.nBins = nBins;
   this.xLow  = xLow;
   this.xHigh = xHigh;

   freqArr = new int[nBins];
   delBin = (xHigh-xLow)/(double)nBins;

   //chart = new DataChart(dataString, _width, _height);

  }
  
  public void print()
  {
    double[] bins = getBins();
    System.out.println("under flow\t" + underFlows);
    for(int i=0;i<freqArr.length-1;i++)
	System.out.println(bins[i] + "\t" + freqArr[i]);

    System.out.println((bins[bins.length-1]) + "\t" + (freqArr[freqArr.length-1] + overFlows));

  }
  
  public String getResult()
  {
    StringBuffer sb = new StringBuffer();
    double[] bins = getBins();
    sb.append("0\t" + underFlows).append("\n");
    for(int i=0;i<freqArr.length-1;i++)
	sb.append(bins[i] + "\t" + freqArr[i]).append("\n");

    sb.append((bins[bins.length-1]) + "\t" + (freqArr[freqArr.length-1] + overFlows)).append("\n");

    return sb.toString();

  }

  public int[] getFreqArr()
  {
    return freqArr;
  }

  public double[] getBins()
  {
    double[] arr = new double[freqArr.length];

    for(int i=0;i<arr.length;i++)
    {
	arr[i] = delBin*i + xLow;
    }

    return arr;
  }

  public int getOverFlows()
  {
    return overFlows;
  }

  public int getUnderFlows()
  {
    return underFlows;
  }
  
  public double getDelBin()
  {
    return delBin;
  }

  public void loadData(List<Double> l)
  {
      for(Iterator<Double> itr=l.iterator(); itr.hasNext(); )
      {
	setData( itr.next() );
      }
  }

  //----------------------------------------------------------------
  void setData(double data){
   
   if( data < xLow)
     underFlows++;
   else if ( data >= xHigh) 
     overFlows++;
   else{
     int bin = (int)((data-xLow)/delBin);
     if(bin >=0 && bin < nBins) freqArr[bin]++;
   }
  }    

  //----------------------------------------------------------------
  // To display the histogram in a chart, we need to pass the data
  // as a string. 
  public void graphIt(){
    dataString = "";
    for (int i=0; i<nBins; i++){
      dataString += freqArr[i] + " ";
    }
  }
  

 	
}

