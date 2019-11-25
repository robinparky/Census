package scripts.labelfree;

import java.util.*;
import java.text.*;
import java.io.*;
import java.math.*;

public class NormalizeLabelFree 
{
    public static void main(String[] args) throws IOException
    {
	if(args.length<1) {
	    System.out.println("Usage: java NormalizeLabelFree input_file");
	    System.exit(0);
	}

	//int experimentNum = Integer.parseInt(args[1]);
	//NormalizeLabelFree q = new NormalizeLabelFree(args[0]);
	NormalizeLabelFree q = new NormalizeLabelFree();
	double[] correctArr = q.getNormalizeFactor(args[0]);

	//for(double d : correctArr)
	  //  System.out.println(d);

	normalize(args[0], correctArr);
    }

    public double[] getNormalizeFactor(String fileName) throws IOException {

        BufferedReader br = new BufferedReader(new FileReader(fileName));            

        String eachLine = br.readLine();

	String[] headArr = eachLine.split("\t");
	int experimentNum = headArr.length;

	double[] totalRatio = new double[experimentNum];
	int indexStartSpC=3;
	gnu.trove.TDoubleArrayList[] ratioArr = new gnu.trove.TDoubleArrayList[experimentNum];

	double[] correctArr = new double[experimentNum];

	for(int i=0;i<ratioArr.length;i++)
	    ratioArr[i] = new gnu.trove.TDoubleArrayList();

	while( (eachLine = br.readLine())!=null )
	{
	    if(!eachLine.startsWith("["))
		continue;

	    String[] arr = eachLine.split("\t");

	    boolean skip = false;
	    for(int i=arr.length-experimentNum;i<arr.length;i++)
	    {
		double d = Double.parseDouble(arr[i]);
		if(d<=0)
		    skip = true;
	    }

	    if(skip)
		continue;

	    for(int i=arr.length-experimentNum;i<arr.length;i++)
	    {
		double tmpDouble = Double.parseDouble(arr[i]);
		//System.out.println("==" + tmpDouble + "\t");
		tmpDouble = round(tmpDouble,3);
		//System.out.println("==" + (i-arr.length+experimentNum));
		totalRatio[i-arr.length+experimentNum] += tmpDouble;
		totalRatio[i-arr.length+experimentNum] = round(totalRatio[i-arr.length+experimentNum], 3);

		double logRatio = Math.log(tmpDouble)/Math.log(2);

		ratioArr[i-arr.length+experimentNum].add(logRatio);
	    }
	}

	int count=0;
	for(gnu.trove.TDoubleArrayList t:ratioArr) {

	    edu.scripps.pms.util.stats.Histogram hist = new edu.scripps.pms.util.stats.Histogram(120, -3.0, 3.0);
	    hist.setData( t.toNativeArray() );

	    double[] bins = hist.getBins();
	    int[] freqArr = hist.getFreqArr();

	    double maxBin = 0;
	    double maxValue = 0;
	    //for(int i=0;i<freqArr.length-1;i++) {
	    for(int i=1;i<freqArr.length-2;i++) {
		bins[i] = edu.scripps.pms.util.MathUtil.getScaled(bins[i], 3);

		if(maxValue<freqArr[i]) {
		    maxValue = freqArr[i];
		    maxBin = bins[i];
		}
	    }

	    int underFlow = hist.getUnderFlows()+ freqArr[0];

	    correctArr[count] = maxBin;

	    count++;

/*
	    System.out.println("median: " + maxBin + "\t" + maxValue);
	    System.out.println(bins[0] + "\t" + underFlow);
	    for(int i=1;i<freqArr.length-1;i++) {
		System.out.println(bins[i] + "\t" + freqArr[i]);
	    }
	    System.out.println(bins[bins.length-1] + "\t" + (freqArr[freqArr.length-1] + hist.getOverFlows()));
*/
	}

	return correctArr;
    }



    public NormalizeLabelFree() {
    }

    public static void normalize(String fileName, double[] correctArr) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(fileName));            

        String eachLine = br.readLine();

	String[] headArr = eachLine.split("\t");
	int experimentNum = headArr.length;

	int indexStartSpC=3;

	while( (eachLine = br.readLine())!=null )
	{
	    System.out.print(eachLine.trim());
	    String[] arr = eachLine.split("\t");

	    if(eachLine.startsWith("P")) {
		int indexCount = 0;
		for(int i=arr.length-correctArr.length-1;i<arr.length-1;i++) {
		    Double d = Double.parseDouble(arr[i]);
		    double logRatio = Math.log(d)/Math.log(2);
		    logRatio -= correctArr[indexCount];

		    double correctR = Math.pow(2,logRatio);

		    System.out.print("\t" + correctR);

		    indexCount++;
		}

	    } else {

		int indexCount = 0;
		for(int i=arr.length-correctArr.length;i<arr.length;i++) {
		    Double d = Double.parseDouble(arr[i]);
		    double logRatio = Math.log(d)/Math.log(2);
		    logRatio -= correctArr[indexCount];

		    double correctR = Math.pow(2,logRatio);

		    System.out.print("\t" + correctR);

		    indexCount++;
		}
	    }

	    System.out.println("");
	}
    }

    // positive value only.
    public double round(double value, int decimalPlace)
    {
	double power_of_ten = 1;
	// floating point arithmetic can be very tricky.
	// that's why I introduce a "fudge factor"
	double fudge_factor = 0.05;
	while (decimalPlace-- > 0) {
	    power_of_ten *= 10.0d;
	    fudge_factor /= 10.0d;
	}
	return Math.round((value + fudge_factor)* power_of_ten)  / power_of_ten;
    }

    class Protein {
	private String locus;
	private String description;
	private String pepCount;
	List<Double> ratioList = new ArrayList<Double>();
	List<Double> normRatioList = new ArrayList<Double>();
	List<Integer> spCList = new ArrayList<Integer>();
	List<Double> spCRatioList = new ArrayList<Double>();

	public void addSpC(int i) {
	    spCList.add(i);
	}

	public void addSpCRatio(double d) {
	    spCRatioList.add(d);
	}

	public void printSpC() {
	    for(Iterator<Integer> itr=spCList.iterator(); itr.hasNext(); ) {
		Integer spC = itr.next();
		System.out.print(spC);
		System.out.print("\t");
	    }
	}

	public void printSpCRatio() {
	    for(Iterator<Double> itr=spCRatioList.iterator(); itr.hasNext(); ) {
		Double spCRatio = itr.next();
		System.out.print(spCRatio);
		System.out.print("\t");
	    }
	}
	public void addRatio(double r) {
	    ratioList.add(r);
	}

	public void normalizeRatio(double[] normFactor) {
	    for(int i=0;i<ratioList.size();i++) {
		normRatioList.add( round(ratioList.get(i)*normFactor[i], 3));
	    }
	}

	public List<Double> getNormRatioList() {
	    return normRatioList;
	}

	public List<Double> getRatioList() {
	    return ratioList;
	}

	public void setLocus(String locus) {
	    this.locus = locus;
	}

	public String getLocus() {
	    return this.locus;
	}

	public void setPepCount(String pepCount) {
	    this.pepCount = pepCount;
	}

	public String getPepCount() {
	    return this.pepCount;
	}

	public void setDescription(String description) {
	    this.description = description;
	}

	public String getDescription() {
	    return this.description;
	}
	

    }
}
