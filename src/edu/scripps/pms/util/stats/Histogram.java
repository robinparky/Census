package edu.scripps.pms.util.stats;

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
import java.io.*;	
import org.apache.commons.lang3.ArrayUtils;

//********************************************************************
public class Histogram {

  int [] freqArr = null;
  int nBins;
  double xLow,xHigh;
  double delBin;
 
  int overFlows=0,underFlows=0;

  String dataString=null;


  public static void main(String[] args) throws Exception
  {
	  if(args.length<4) {
		  System.out.println("Usage: histogram file_name bin_size low_bound high_bound");
		  System.exit(0);
	  }
	  int nBins = Integer.parseInt(args[1]);
	  double low = Double.parseDouble(args[2]);
	  double high = Double.parseDouble(args[3]);

	  BufferedReader br = new BufferedReader(new FileReader(args[0]));

	  String eachLine;

	  List<Double> list = new ArrayList<Double>();

	  while( (eachLine = br.readLine()) != null )
	  {
		  
		  list.add(Double.parseDouble(eachLine));
	  }

	  Double[] doArr = list.toArray(new Double[list.size()]);
	  double[] darr = ArrayUtils.toPrimitive(doArr);

	  //Histogram h = new Histogram(20, -4, 6);
	  Histogram h = new Histogram(nBins, low, high);
	  h.setData(darr);
	  int[] freqArr = h.getFreqArr();
	  double[] bins = h.getBins();

	  System.out.println("under flow\t" + h.getUnderFlows());
	  System.out.println("over flow\t" + h.getUnderFlows());
	  System.out.println("bin size\t" + h.getDelBin());
	  for(int i=0;i<freqArr.length;i++)
		  System.out.println(bins[i] + "\t" + freqArr[i]);
//	  h.print();
/*
    Histogram h = new Histogram(20, -4, 6);

    double[] arr = {-2.0, -3.4, 0.2, -4.5, 2.3,4.6, 3.2, 3.1,8.6,5.7,9.1,9.7, 11, 12};

    for(double d:arr)
    {
	h.setData(d);
    }
    h.print();
*/
/*
    int[] freqArr = h.getFreqArr();
    double[] bins = h.getBins();


    System.out.println("0\t" + h.getUnderFlows());
    for(int i=0;i<freqArr.length;i++)
	System.out.println(bins[i] + "\t" + freqArr[i]);

    System.out.println((bins[bins.length-1] + h.getDelBin()) + "\t" + h.getOverFlows());
    System.out.println(h.getDelBin());
*/
  }

  //----------------------------------------------------------------
  public Histogram (int nBins, double xLow, double xHigh){
  
   this.nBins = nBins;
   this.xLow  = xLow;
   this.xHigh = xHigh;

   freqArr = new int[nBins];
   delBin = (xHigh-xLow)/(double)nBins;
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

  public int[] getFreqArrSmooth(int smoothWindow) {
	  int arrSize = freqArr.length;

	  if(arrSize<=smoothWindow)
		  return freqArr;

	  int[] resultArr = new int[arrSize];

	  for(int i=0;i<arrSize;i++) {

		  int start=i-smoothWindow/2;
		  start = (start>=0)?start:0;
		  int end=i+smoothWindow/2;
		  end = (end>=arrSize)?arrSize-1:end;


		  for(int j=start;j<=end;j++)
			  resultArr[i] += freqArr[j];

		  //      System.out.println(i + "\t" + start + "\t" + end + "\t" + (i+smoothWindow));

	  }

	  return resultArr;
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

  public void setData(double[] arr) {
    for(double d:arr) {
	setData(d);
    }
  }
  
  public void setData(long[] arr) {
    for(long l:arr) {
	setData( (double)l );
    }
  }
  

  //----------------------------------------------------------------
  public void setData(double data){
  
   if( data < xLow) { underFlows++; }
   else if ( data >= xHigh) { overFlows++; }
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

