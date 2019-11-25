package scripts.singleton;

import java.util.*;
import java.io.*;


public class HistogramGenerator
{
    public static Set<String> getList(String fileName, double comScoreThresh) throws IOException
    {
	BufferedReader br = new BufferedReader(new FileReader(fileName));
        String eachLine;

	double lowerBound=0.1;
	double upperBound=10.0;
	double linearRegThreshold=0.5;

	//List<Double> list = new ArrayList<Double>();
	//List<Double> tmplist = new ArrayList<Double>();
	Set<String> set = new HashSet<String>();
	Set<String> tmpSet = new HashSet<String>();

	boolean isProtein=true;

	StringBuffer sb = new StringBuffer();

	int totalPepCount = 0;
	int totalProCount = 0;

	while( (eachLine = br.readLine())!=null )
	{
	    if(eachLine.startsWith("P"))
	    {
		isProtein=true;

		if(tmpSet.size()>=3)
		{
		    set.addAll(tmpSet);

		    totalPepCount += tmpSet.size();
		    totalProCount++;

		    //System.out.print(sb.toString());

		}
		sb = new StringBuffer();

		tmpSet.clear();

	    }

	    if(!eachLine.startsWith("S"))
		continue;

	    isProtein=false;

	    String[] arr = eachLine.split("\t");

	    double comScore = Double.parseDouble(arr[arr.length-2]);
	    double area1 = Double.parseDouble(arr[arr.length-5]);
	    double area2 = Double.parseDouble(arr[arr.length-4]);
	    double detValue = Double.parseDouble(arr[5]);
	    String areaRatio = arr[arr.length-3];
	    String seq = arr[2];
	    String desc = arr[arr.length-1];

	    if(detValue>=linearRegThreshold)
		continue;

	    //if(area1<=0 || area2<=0)
	//	continue;

	    double ratio = area1/area2;

	    if(area2<=0 || (ratio>=upperBound || ratio<=lowerBound))
	    {
		if(comScore<comScoreThresh)
		    continue;

//		System.out.println(comScore + "&" + eachLine);
//		System.out.println(detValue + "\t" +  linearRegThreshold);
		//System.out.println(area1 + "\t" + area2 + "\t" + ratio + "\t" +  comScore);
		tmpSet.add(seq + desc);
		sb.append(eachLine).append("\n");
		//System.out.println(eachLine);
	    }
	}

	    System.out.println(comScoreThresh + "\t" + totalProCount + "\t" + totalPepCount + "\t" + set.size());
	return set;
    }

    public static void main(String[] args) throws IOException
    {
	double thresh = 0;

	//Set s = getList(args[0], 0.95);
	//System.out.println(s.size());

	//if(true)
	//return;

	for(double d=0;d<=1.0;d+=0.01)
	{
	    d = round(d, 2);
	    Set<String> set= getList(args[0], d);
	 //   System.out.println(d + "\t" + set.size());
	}

    }

    public static double round(double val, int places) {
	long factor = (long)Math.pow(10,places);

	// Shift the decimal the correct number of places
	// to the right.
	val = val * factor;

	// Round to the nearest integer.
	long tmp = Math.round(val);

	// Shift the decimal the correct number of places
	// back to the left.
	return (double)tmp / factor;
    }
}
