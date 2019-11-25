package scripts.xmhan;

import java.io.*;
import java.util.*;
import edu.scripps.pms.census.model.*;
import edu.scripps.pms.census.util.*;
import edu.scripps.pms.census.util.io.FastaReader;
import edu.scripps.pms.util.seq.Fasta;

import gnu.trove.TDoubleArrayList;
import edu.scripps.pms.util.stats.*;

/**
 *
 * @author rpark
 * @version $Id:$
 */
public class TMTReportParser {
    
    public static void main(String args[]) throws Exception
    {
//	String fileName = "/home/rpark/rpark_on_data/project/xmhan/TMT_MEF_6plex-PQD2-013009/census-out.txt";
	String fileName = args[0];
	double normValue = 0.251155649;  //log scale
	System.out.println("threshold\tmean\tstdev\tmean-2*sigma\tmean+2*sigma\tmean-1\tmean+1");


	    analyze(fileName, 0, normValue, 500000);
	    analyze(fileName, 500000, normValue, 1500000);
	    analyze(fileName, 1500000, normValue, 2500000);
	    analyze(fileName, 2500000, normValue, 3500000);

	int num=0;
	for(int i=0;i<40;i++) {
	    analyze(fileName, num, normValue, num+1000000);
		num+=1000000;

if(true)  return;
	}


	double threshold = 0;
	for(int i=1;i<150;i++) {
	    boolean isDone = analyze(fileName, threshold, normValue,-1);
//	    if(true)
//		break;

	    threshold += 500000;

	}
    }

    public static boolean analyze(String fileName, double threshold, double normValue, double maxValue) throws Exception {
        BufferedReader br = new BufferedReader(new FileReader(fileName));

        String eachLine;
        while( (eachLine = br.readLine()).startsWith("H\t") );

	TDoubleArrayList list = new TDoubleArrayList();
        
        while( null != (eachLine = br.readLine()) )
	{
	    if(!eachLine.startsWith("S"))
		continue;

	    String[] arr = eachLine.split("\t");
//	    double m126 = Double.parseDouble(arr[3]);
//	    double m131 = Double.parseDouble(arr[8]);
	    double m1 = Double.parseDouble(arr[3]);
	    double m2 = Double.parseDouble(arr[5]);
	    //double m131 = Double.parseDouble(arr[4]);

	    if(0==m1 || 0==m2)
		continue;

	    double sum = m1 + m2;


	    if(threshold > sum)
		continue;

if(maxValue>0 && sum>maxValue) continue;

//	    System.out.println(m126 + " " + m131);

//		continue;
	    double ratio = m1/m2;
	    double logRatio = Math.log(ratio)/Math.log(2);
	    logRatio += normValue;
	    list.add(logRatio);
	    System.out.println(logRatio);
	}

	    //System.out.println(list.size() + "\t");
	double[] darr = list.toNativeArray();
	StatCalc calc = new StatCalc(darr);
	double stdev = calc.getStandardDeviation();
	double mean = calc.getMean();

	System.out.print(threshold);
	System.out.print("\t");
	System.out.print(mean);
	System.out.print("\t");
	System.out.print(stdev);
	System.out.print("\t");
	System.out.print((mean-2*stdev));
	System.out.print("\t");
	System.out.print((mean+2*stdev));
	System.out.print("\t");
	System.out.print((mean-1));
	System.out.print("\t");
	System.out.println((mean+1));

	
//	System.out.println( (mean-2*stdev) + "\t" + (mean-1) + "\t" + ((mean-2*stdev)>(mean-1)) );
//	System.out.println( (mean+2*stdev) + "\t" + (mean+1) + "\t" + ((mean+2*stdev)<(mean+1)) );
	if( (mean-2*stdev)>(mean-1) && (mean+2*stdev)<(mean+1) )
	{
	    //for(double d:darr)
	//	System.out.println(d);
	    return true;
	}
	else 
	    return false;
    }
}
