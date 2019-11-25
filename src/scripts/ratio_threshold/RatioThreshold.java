package scripts.ratio_threshold;

import java.io.*;
import java.util.*;
import edu.scripps.pms.census.model.*;
import edu.scripps.pms.util.stats.*;

public class RatioThreshold
{
    public static void main(String args[]) throws Exception
    {
	//args1 = /home/rpark/rpark_on_data/project/bingwlu/silac_expression_pvalue/20070713/light/census-out.txt
	//args2 = 100

	if(args.length<2)
	{
	    System.out.println("Usage: java RatioThreshold census-out.txt batchSize");
	    System.exit(0);
	}

	int batchSize = Integer.parseInt(args[1]);
	String fileName = args[0];
	runPeptide(fileName, batchSize);

//	while( !(eachLine = br.readLine()).startsWith("H\tPLINE\tLOC") );
//	System.out.println(eachLine);
//P       IPI00218621.1   0.45            0.45    1       3       1.4132804E7     Tax_Id=9606 Isoform J of SON protein

    }

    public static void runPeptide(String fileName, int batchSize) throws Exception {

	BufferedReader br = new BufferedReader(new FileReader(fileName));
	String eachLine;
	List<ChroPeptide> list = new ArrayList<ChroPeptide>();
	while( (eachLine = br.readLine()) != null )
	{
	    if(!eachLine.startsWith("S\t")) 
		continue;

	    String[] arr = eachLine.split("\t");

	    ChroPeptide pep = new ChroPeptide();
	    pep.setRatio(Double.parseDouble(arr[3]));
	    pep.setFileName(arr[12]);
	    pep.setSequence(arr[2]);
	    pep.setScanNum(Integer.parseInt(arr[13]));
	    pep.setChargeState(Integer.parseInt(arr[14]));

	    //pro.setAverageRatio(Double.parseDouble(arr[2]));
	    pep.setTotalIntensity(Double.parseDouble(arr[8]) + Double.parseDouble(arr[9]));

	    list.add(pep);
        }

	Collections.sort(list, new ChroPeptideByIntensity());

	List<ChroPeptide> batchList = new ArrayList<ChroPeptide>();

	int count=0;
	for(Iterator<ChroPeptide> itr=list.iterator(); itr.hasNext(); ) {
	    ChroPeptide pep = itr.next();

	    count++;

	    batchList.add(pep);

	    if(count%100 == 0) {
		int remain = list.size()-count;
		if(remain>batchSize) {
		    runStat(batchList);
		    batchList.clear();
		}
	    }
	}

	runStat(batchList);

	Hashtable<String, ChroPeptide> ht = new Hashtable<String, ChroPeptide>();

	for(Iterator<ChroPeptide> itr=list.iterator(); itr.hasNext(); ) {

	    ChroPeptide pep = itr.next();
	    String key = pep.getSequence() + pep.getScanNum() + pep.getChargeState();
	    ht.put(key, pep);
	}

	printPeptide(fileName, ht);
    }

    public static void printPeptide(String fileName, Hashtable<String, ChroPeptide> ht) throws Exception {

	BufferedReader br = new BufferedReader(new FileReader(fileName));
	String eachLine;

	while( (eachLine = br.readLine()) != null )
	{
	    if(!eachLine.startsWith("S\t")) {
		System.out.println(eachLine);
		continue;
	    }

	    String[] arr = eachLine.split("\t");
	    String key = arr[2] + arr[13] + arr[14];

	    ChroPeptide pep = ht.get(key);
	    System.out.println(eachLine + "\t" + pep.getPvalue());
	}
    }

    public static void runStat(List<ChroPeptide> list) throws Exception {
	int size = list.size();
	double[] arr = new double[size];
	double sum = 0;

	for(int i=0;i<arr.length;i++) {
	    ChroPeptide pep = list.get(i);
	    arr[i] = Math.log(pep.getRatio())/Math.log(2);
	    sum += arr[i];
	}

	double mean = sum / size;

	double stdev = StatCalc.getStandardDeviation(arr);

	for(int i=0;i<arr.length;i++) {
	    ChroPeptide pep = list.get(i);

	    double z = (pep.getRatio() - mean) / stdev;
	    double p = StatCalc.zScore2PValue(z);
	    pep.setPvalue(p);

	}
    }

    public static void runProtein(String fileName, int batchSize) throws Exception {

	BufferedReader br = new BufferedReader(new FileReader(fileName));
	String eachLine;
	List<ChroProtein> list = new ArrayList<ChroProtein>();
	while( (eachLine = br.readLine()) != null )
	{
	    if(!eachLine.startsWith("P\t")) 
		continue;

	    String[] arr = eachLine.split("\t");

	    ChroProtein pro = new ChroProtein();
	    pro.setLocus(arr[1]);
	    pro.setAverageRatio(Double.parseDouble(arr[2]));
	    pro.setIntensity(Double.parseDouble(arr[arr.length-2]));

	    list.add(pro);
        }

	Collections.sort(list, new ChroProteinByIntensity());

	for(Iterator<ChroProtein> itr=list.iterator(); itr.hasNext(); ) {
	    ChroProtein pro = itr.next();
//	    System.out.println(pro.getIntensity());
	}
    }
}

