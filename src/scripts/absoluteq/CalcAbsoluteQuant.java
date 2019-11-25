
/*
 * CalcAbsoluteQuant.java
 *
 * Created on February 23, 2006, 3:03 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

//package scripts;

import java.io.*;
import java.util.*;
import edu.scripps.pms.census.model.*;
import edu.scripps.pms.census.util.*;
import edu.scripps.pms.census.util.io.FastaReader;
import edu.scripps.pms.util.seq.Fasta;

/**
 *
 * @author rpark
 * @version $Id: CalcAbsoluteQuant.java,v 1.1 2008/10/30 21:15:15 rpark Exp $
 */
public class CalcAbsoluteQuant {

    private static String refProtein;
    private static double refProteinAmount;

    public void readParam(String paramFile) throws Exception {
        //BufferedReader br = new BufferedReader(new FileReader("census.params"));
        BufferedReader br = new BufferedReader(new FileReader(paramFile));

        String eachLine;
        while( (eachLine = br.readLine()) != null ) {
	    String[] arr = eachLine.split("\t");
		    
	    if(eachLine.startsWith("ref_protein\t"))
		refProtein = arr[1];
		refProteinAmount = Double.parseDouble(arr[2]);
	}
    }
    
    public static void main(String args[]) throws Exception
    {
//	String refProtein = args[0];
//
	if(args.length<1)
	{
	    System.out.println("Usave: java CalcAbsoluteQuant quant_compare_file census_param");
	    System.exit(0);
	}

	String filename = args[0];

	CalcAbsoluteQuant quant = new CalcAbsoluteQuant();

	quant.readParam(args[1]);
	int[] indexArr = quant.getIntensityIndexList(filename);
	double[] refProteinIntArr = quant.refProteinValue(filename, indexArr);

	quant.calcAbsQuant(refProteinIntArr, filename, indexArr);

    }

    public int[] getIntensityIndexList(String filename) throws Exception {
        BufferedReader br = new BufferedReader(new FileReader(filename));

        String eachLine;
        while( !(eachLine = br.readLine()).startsWith("SLINE\t") );

	gnu.trove.TIntArrayList intIndexList = new gnu.trove.TIntArrayList();
	String[] tmpArr = eachLine.split("\t");
	for(int i=0;i<tmpArr.length;i++) {
	    if(tmpArr[i].startsWith("INTENSITY"))
		intIndexList.add(i);
	}

	return intIndexList.toNativeArray();
    }




    public double[] refProteinValue(String filename, int[] indexArr) throws Exception {
        BufferedReader br = new BufferedReader(new FileReader(filename));

        String eachLine;
        
        while( (eachLine = br.readLine()) != null  )
	{
	    if(eachLine.startsWith("P\t" + refProtein))
		break;
	}


	double[] intSum = new double[indexArr.length];
	int[] intCount = new int[indexArr.length];
        while( (eachLine = br.readLine()) != null && !eachLine.startsWith("P\t") ) {
	    String[] arr = eachLine.split("\t");

	    for(int i=0;i<indexArr.length;i++) {
		//System.out.print(arr[indexArr[i]] + "\t");

		if("NA".equals(arr[indexArr[i]]))
			continue;

		double d = Double.parseDouble(arr[indexArr[i]]);
		intSum[i] += d;
		intCount[i]++;

		
	    }

	}

//	for(int i:intCount)
//	    System.out.println(i);
//	for(double i:intSum)
//	    System.out.println(i);

	double[] refProteinIntArr = new double[intSum.length];
	for(int i=0;i<intSum.length;i++) 
	    refProteinIntArr[i] = intSum[i]/intCount[i];

	return refProteinIntArr;

    }

    public double[] calcAbsQuant(double[] refProteinIntArr, String filename, int[] indexArr) throws Exception {
        BufferedReader br = new BufferedReader(new FileReader(filename));

        String eachLine;
	List<String> proteinList=new ArrayList<String>();
	List<String> peptideList=new ArrayList<String>();

	while( (eachLine = br.readLine())!=null )
	{
	    if(eachLine.startsWith("P\t")) { //protein line
		proteinList.add(eachLine);

		while ( (eachLine = br.readLine()) != null && eachLine.startsWith("P\t"))
		    proteinList.add(eachLine);

		peptideList.add(eachLine);
                
		while ( (eachLine = br.readLine()) != null && eachLine.startsWith("S\t"))
		    peptideList.add(eachLine);

                absCalc(indexArr, peptideList, proteinList, refProteinIntArr);
		proteinList.clear();
		peptideList.clear();

		if(null != eachLine && eachLine.startsWith("P\t"))
		    proteinList.add(eachLine);

	    } else if(eachLine.startsWith("S\t")) {
		peptideList.add(eachLine);
		while ( (eachLine = br.readLine()) != null && eachLine.startsWith("S\t")) {
		    peptideList.add(eachLine);
		}

                absCalc(indexArr, peptideList, proteinList, refProteinIntArr);

		proteinList.clear();
		peptideList.clear();

		if(null != eachLine && eachLine.startsWith("P\t"))
		    proteinList.add(eachLine);
	    }
	}

	return null;
    }

    private double[] absCalc(int[] indexArr, List<String> peptideList, List<String> proteinList, double[] refProteinIntArr) throws Exception {

	double[] intSum = new double[indexArr.length];
	int[] intCount = new int[indexArr.length];

	for(Iterator<String> itr=peptideList.iterator(); itr.hasNext(); ) {

	    String[] arr = itr.next().split("\t");

	    for(int i=0;i<indexArr.length;i++) {
		if("NA".equals(arr[indexArr[i]]))
			continue;

		double d = Double.parseDouble(arr[indexArr[i]]);
		intSum[i] += d;
		intCount[i]++;
	    }
	}

	double[] avgIntArr = new double[intSum.length];

	StringBuffer sb = new StringBuffer();
	for(int i=0;i<avgIntArr.length;i++) {
	    avgIntArr[i] = intSum[i]/intCount[i];



	    double absAmt = avgIntArr[i] * refProteinAmount / refProteinIntArr[i] ;
	    sb.append(absAmt).append("\t");
//	    System.out.println("==" + avgIntArr[i] + "\t" + absAmt);
	}

	for(Iterator<String> itr=proteinList.iterator(); itr.hasNext(); ) {

	    System.out.print(itr.next());
	    System.out.print(sb.toString());
	}

	System.out.println("");
	return null;
    }
}
