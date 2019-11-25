package scripts.labelfree;

import java.util.*;
import java.text.*;
import java.io.*;
import java.math.*;

import edu.scripps.pms.util.stats.*;
import edu.scripps.pms.census.model.ChroPeptide;
import edu.scripps.pms.census.tools.Formatter;
import gnu.trove.TDoubleArrayList;
import org.jfree.ui.LengthLimitingDocument;

public class NormalizeLabelFreeIntensityBase
{
    private static List<Sample> sampleList = new ArrayList<Sample>();

    public static void main(String[] args) throws IOException
    {

	/*if(args.length<1) {
	    System.out.println("Usage: java NormalizeLabelFreeIntensityBase input_file num_of_min_id_per_replicates");
	    System.exit(0);
	}*/

	//int experimentNum = Integer.parseInt(args[1]);
	//NormalizeLabelFree q = new NormalizeLabelFree(args[0]);
	

	
	NormalizeLabelFreeIntensityBase q = new NormalizeLabelFreeIntensityBase();
	//double[] correctArr = q.getNormalizeFactor(args[0]);
	//int numIdPerReplicates = Integer.parseInt(args[1]);
	//args[0] = "/data/2/rpark/ip2_data/crestani/Iron_deprivation_Cgattii/labelfree_quant/census_labelfree_out_1085.txttmp";
        //args[0] = "/data/2/rpark/ip2_data/pankows/CFTR/labelfree_quant/census_labelfree_out_1123.txttmp";
        //args[0] = "/data/2/rpark/ip2_data/xmhan/FAd_labelfree/labelfree_quant/census_labelfree_out_1596_tempstat.txt";
        //getNormalizeWithLocus(args[0], "out.txt", "IPI00220030.1", 0);

        //args[0] = "/data/2/rpark/ip2_data/xmhan/FAd_labelfree/labelfree_quant/census_labelfree_out_1855_stat.txt";
        //getNormalizeWithLocusSpecCount(args[0], "out.txt", "IPI00220030", 0);


	//String outputFilenameTmp = "/home/rpark/pms/Census/build/classes/census_labelfree_out.txttmp";
        String outputFilenameTmp = "/home/rpark/test.txttmp";
        System.out.println("normalize...");
        double[] correctArr = q.getNormalizeFactor(outputFilenameTmp); //"census_label_free_result_with_id_temp.txt");
        q.normalize("/home/rpark/out.txt", outputFilenameTmp, correctArr, 0);



	//double[] correctArr = q.getNormalizeFactor(args[0], "IPI00302383.3");
	//int numIdPerReplicates = Integer.parseInt("0");

	//q.normalize("outputfilename", args[0], correctArr, numIdPerReplicates);
    }

    public int[] getIntensityIndexArr(int experimentNum, String sline) {
	    
	String[] arr = sline.split("\t");
	int[] indexArr = new int[experimentNum];
	int tmpIndex=0;

	for(int i=0;i<arr.length;i++) {
	    
	    if("INTENSITY".equals(arr[i])) {
		indexArr[tmpIndex++] = i;
	    }
	}

	return indexArr;
    }

    public int[] getIndexArrByPrefix(int experimentNum, String sline, String prefix) {
	    
	String[] arr = sline.split("\t");
	int[] indexArr = new int[experimentNum];
	int tmpIndex=0;
        
	for(int i=0;i<arr.length;i++) {
	    
	    if(arr[i].startsWith(prefix)) {
                
		indexArr[tmpIndex++] = i;
	    }
	}

	return indexArr;
    }

    public class Sample {

        private int replicateCount;
        
        public Sample(int replicateCount) {
            this.replicateCount = replicateCount;

        }

        public int getReplicateCount() {
            return replicateCount;
        }

        public void setReplicateCount(int replicateCount) {
            this.replicateCount = replicateCount;
        }
    }

    public double[] getNormalizeFactor(String fileName) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(fileName));

	sampleList = new ArrayList<Sample>();

        String eachLine; // = br.readLine();
        int experimentNum=0;
	while( (eachLine = br.readLine())!=null )
	{
            if(!eachLine.startsWith("H\tGROUP_SAMPLE"))
                break;

            String[] headArr = eachLine.split("\t");
            int tmpNum = headArr.length-3;
            experimentNum += tmpNum;


	    Sample s = new Sample(tmpNum);
            sampleList.add(s);
        }

	while( (eachLine = br.readLine())!=null && !eachLine.startsWith("S") );

	//int[] intensityIndexArr = getIntensityIndexArr(experimentNum, eachLine);
	int[] intensityIndexArr = getIndexArrByPrefix(experimentNum, eachLine, "NORM_INTENSITY");

	double[] totalRatio = new double[experimentNum];
	//int intenseIndex=6;
	//int increment=7;
	gnu.trove.TDoubleArrayList[] intRatioArr = new gnu.trove.TDoubleArrayList[experimentNum];

	for(int i=0;i<intRatioArr.length;i++)
	    intRatioArr[i] = new gnu.trove.TDoubleArrayList();

	while( (eachLine = br.readLine())!=null )
	{

	    if(!eachLine.startsWith("S\t")) {
		continue;
	    }

	    String[] arr = eachLine.split("\t");

	    int expCount = 0;


	    //double referInt = Double.parseDouble( arr[intensityIndexArr[0]] );
	    //if(0 == referInt)
		//`continue;

	    for(int i=0;i<intensityIndexArr.length;i++) {

		String intenseCol = arr[intensityIndexArr[i]];
                if(intenseCol.equals("NA")){
                    continue;
                }
		double d = Double.parseDouble(intenseCol);
		if(0 == d || experimentNum==d) {
		    //remove extreme values
		//    System.out.println(d);
		    continue;
		}


		intRatioArr[i].add(d);
	    }

	}


	//double[] tmpCorrectArr = new double[experimentNum-1];

	double[] correctArr = new double[experimentNum];
	int count=0;
	double totalCorrect=0;
	for(gnu.trove.TDoubleArrayList t:intRatioArr) {

	    org.apache.commons.math.stat.descriptive.DescriptiveStatistics desc = new org.apache.commons.math.stat.descriptive.DescriptiveStatistics();
	    for(double dd:t.toNativeArray())
		desc.addValue(dd);

	    correctArr[count++] = desc.getMean();
	    totalCorrect += desc.getMean();

	}

	double averageCorrect = totalCorrect/correctArr.length;
	//for(int i=0;i<correctArr.length;i++)
	//    System.out.println("==" + correctArr[i]);

	//System.out.println(averageCorrect);

	for(int i=0;i<correctArr.length;i++)
	    correctArr[i] = averageCorrect-correctArr[i];

	//for(int i=0;i<correctArr.length;i++)
	//    System.out.println("==" + correctArr[i]);

	return correctArr;
        
    }

    public static void getNormalizeWithLocusSpecCount(String inputFilename, String outputFilename, String locus, int numIdPerReplicates) throws IOException {


        String tmpOutfilename = inputFilename + "2";
	BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(tmpOutfilename));
	PrintStream p = new PrintStream(out);

        BufferedReader br = new BufferedReader(new FileReader(inputFilename));
        String eachLine; // = br.readLine();
	while( (eachLine = br.readLine())!=null ) {
		if(eachLine.startsWith("PLINE\tLOCUS")) break;
	}

	gnu.trove.TIntArrayList spcIndexList = new gnu.trove.TIntArrayList();
	gnu.trove.TDoubleArrayList spcIndexCorrect = new gnu.trove.TDoubleArrayList();
	String[] tmpArr = eachLine.split("\t");
	for(int i=0;i<tmpArr.length;i++) {
	    if(tmpArr[i].startsWith("NORM_SPEC_COUNT"))
		spcIndexList.add(i);
	}

	int[] indexList = spcIndexList.toNativeArray();

	String searchkey = "P\t" + locus;
	while( (eachLine = br.readLine())!=null )
	{
		if(eachLine.startsWith(searchkey)) break;
	}

	tmpArr = eachLine.split("\t");
	for(int ind:indexList) {
		String[] narr = tmpArr[ind].split(",");
		for(String s:narr)
			spcIndexCorrect.add( Double.parseDouble(s) );
	}


	
//System.out.println("===========" + spcIndexCorrect);
	br.close();


        br = new BufferedReader(new FileReader(inputFilename));
	while( (eachLine = br.readLine())!=null ) {
		if(!eachLine.startsWith("P\t")) {
			p.println(eachLine);

			continue;
		}

		String[] tarr = eachLine.split("\t");

		for(int i=0;i<indexList[0];i++) {
			p.print(tarr[i]);
			p.print("\t");
		}

		int tindex=0;
		for(int i:indexList) {
			String[] narr = tarr[i].split(",");
			for(String s:narr) {
				double corrd = Double.parseDouble(s)/spcIndexCorrect.get(tindex);


				if(0==spcIndexCorrect.get(tindex)) {
					p.print(0);
					p.print(",");
				}
				else {
					p.print(corrd);
					p.print(",");
				}

//				System.out.println( s + "\t" + spcIndexCorrect.get(tindex) + "\t" + corrd);
				tindex++;

			}

			p.print("\t");

		}
	
		for(int i=indexList[indexList.length-1]+1;i<tarr.length;i++) {
			p.print(tarr[i]);
			p.print("\t");
		}

		

		p.println("");

	}

	br.close();
	p.close();
/*
//	List sList = new ArrayList<Sample>();

        int expSize=0;
	int sampleSize=0;
	while( (eachLine = br.readLine())!=null )
	{
            p.println(eachLine);

	    if(eachLine.startsWith("H\tNorm"))
		continue;

            
            if(!eachLine.startsWith("H\tGROUP_SAMPLE") && !eachLine.startsWith("H\tsample_group"))
                break;

            String[] headArr = eachLine.split("\t");
            int tmpNum = headArr.length-3;
            expSize += tmpNum;

	
//	    Sample s = new Sample(tmpNum);
//	    Sample s1 = new Sample(2);
            //sList.add(s);
	    sampleSize++;
        }
        
       // System.out.println("--" + eachLine + inputFilename);

        if(!eachLine.startsWith("PLINE\t")) {
            while( (eachLine = br.readLine()) != null ) {

                if(eachLine.startsWith("PLINE"))
                    break;

                p.println(eachLine);
            }
        }

        //System.out.println("--" + eachLine + inputFilename);
        p.print(eachLine);

        for(int i=0;i<sampleSize;i++) {
            p.print("SAMPLE_GROUP_INTENSITY_AVG_");
            p.print( (i+1) );
            p.print("\t");
        }

        p.println("");

	//sline
	eachLine = br.readLine();
	p.println(eachLine);

	gnu.trove.TIntArrayList intIndexList = new gnu.trove.TIntArrayList();
        //System.out.println(eachLine);
	String[] tmpArr = eachLine.split("\t");
	for(int i=0;i<tmpArr.length;i++) {
	    if(tmpArr[i].startsWith("INTENSITY"))
		intIndexList.add(i);
	}

	int[] indexList = intIndexList.toNativeArray();

	List<String> proteinList=new ArrayList<String>();
	List<String> peptideList=new ArrayList<String>();

	double[] pepNormRatioSumArr = new double[expSize];

//	boolean prevProtein=false;

//	while ( (eachLine = br.readLine()) != null && eachLine.startsWith("P\t"));

        String searchKey = "P\t" + locus;
        boolean foundKey = false;
        double[] intArr = null;
	while( (eachLine = br.readLine())!=null )
	{

	    if(eachLine.startsWith("P\t")) { //protein line
		
		proteinList.add(eachLine);
                if(eachLine.startsWith(searchKey))
                    foundKey = true;

		while ( (eachLine = br.readLine()) != null && eachLine.startsWith("P\t")) {

                    if(eachLine.startsWith(searchKey))
                        foundKey = true;

		    proteinList.add(eachLine);
		}

		peptideList.add(eachLine);

		while ( (eachLine = br.readLine()) != null && eachLine.startsWith("S\t")) {
		    peptideList.add(eachLine);
		}


                if(foundKey) {
                  //  System.out.println("aaaa" + eachLine);
                  //  System.out.println("aaaa" + proteinList);

                    intArr = getCorrectFactorByLocusSpecCount(indexList, peptideList, proteinList, expSize, pepNormRatioSumArr, numIdPerReplicates);

                    break;
                }
		//proteinCalc(indexList, peptideList, proteinList, correctArr, p, expSize, pepNormRatioSumArr, numIdPerReplicates);
		proteinList.clear();
		peptideList.clear();

		if(null != eachLine && eachLine.startsWith("P\t")) {

                    if(eachLine.startsWith(searchKey))
                        foundKey = true;

		    proteinList.add(eachLine); // + "<<"  );
                }

	    } else if(eachLine.startsWith("S\t")) {
		peptideList.add(eachLine);
		while ( (eachLine = br.readLine()) != null && eachLine.startsWith("S\t")) {
		    peptideList.add(eachLine);
		}


                if(foundKey) {
                    //System.out.println("aaaa" + eachLine);
                    //System.out.println("aaaa" + proteinList);
                    intArr = getCorrectFactorByLocusSpecCount(indexList, peptideList, proteinList, expSize, pepNormRatioSumArr, numIdPerReplicates);
                    break;
                }
		//proteinCalc(indexList, peptideList, proteinList, correctArr, p, expSize, pepNormRatioSumArr, numIdPerReplicates);

		proteinList.clear();
		peptideList.clear();

		if(null != eachLine && eachLine.startsWith("P\t")) {
                    if(eachLine.startsWith(searchKey))
                        foundKey = true;

                    proteinList.add(eachLine);

                }
		    
	    }
	}

        
	if(null != p)
	    p.close();


        normalize(outputFilename, inputFilename, intArr, numIdPerReplicates, sampleSize);

*/
        //intArr : correctArr

        //normalizeCorrectArr(intArr);

    }


    public static void getNormalizeWithLocus(String inputFilename, String outputFilename, String locus, int numIdPerReplicates) throws IOException {
        String tmpOutfilename = inputFilename + "2";
	BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(tmpOutfilename));
	PrintStream p = new PrintStream(out);

        BufferedReader br = new BufferedReader(new FileReader(inputFilename));
//	List sList = new ArrayList<Sample>();

        String eachLine; // = br.readLine();
        int expSize=0;
	int sampleSize=0;
	while( (eachLine = br.readLine())!=null )
	{
            p.println(eachLine);

	    if(eachLine.startsWith("H\tNorm"))
		continue;

            
            if(!eachLine.startsWith("H\tGROUP_SAMPLE") && !eachLine.startsWith("H\tsample_group"))
                break;

            String[] headArr = eachLine.split("\t");
            int tmpNum = headArr.length-3;
            expSize += tmpNum;

	
//	    Sample s = new Sample(tmpNum);
//	    Sample s1 = new Sample(2);
            //sList.add(s);
	    sampleSize++;
        }
        
       // System.out.println("--" + eachLine + inputFilename);

        if(!eachLine.startsWith("PLINE\t")) {
            while( (eachLine = br.readLine()) != null ) {

                if(eachLine.startsWith("PLINE"))
                    break;

                p.println(eachLine);
            }
        }

        //System.out.println("--" + eachLine + inputFilename);
        p.print(eachLine);

        for(int i=0;i<sampleSize;i++) {
            p.print("SAMPLE_GROUP_INTENSITY_AVG_");
            p.print( (i+1) );
            p.print("\t");
        }

        /*
        for(int i=0;i<expSize;i++) {
            p.print("NORM_INTENSITY_CORRECT_");
            p.print( (i+1) );
            p.print("\t");
	}

        for(int i=0;i<sampleList.size();i++) {
            p.print("SAMPLE_GROUP_INTENSITY_CORRECT_AVG_");
            p.print( (i+1) );
            p.print("\t");
        }*/

        p.println("");

	//sline
	eachLine = br.readLine();
	p.println(eachLine);

	gnu.trove.TIntArrayList intIndexList = new gnu.trove.TIntArrayList();
        //System.out.println(eachLine);
	String[] tmpArr = eachLine.split("\t");
	for(int i=0;i<tmpArr.length;i++) {
	    if(tmpArr[i].startsWith("INTENSITY"))
		intIndexList.add(i);
	}

	int[] indexList = intIndexList.toNativeArray();

	List<String> proteinList=new ArrayList<String>();
	List<String> peptideList=new ArrayList<String>();

	double[] pepNormRatioSumArr = new double[expSize];

//	boolean prevProtein=false;

//	while ( (eachLine = br.readLine()) != null && eachLine.startsWith("P\t"));

        String searchKey = "P\t" + locus;
        boolean foundKey = false;
        double[] intArr = null;
	while( (eachLine = br.readLine())!=null )
	{

	    if(eachLine.startsWith("P\t")) { //protein line
		
		proteinList.add(eachLine);
                if(eachLine.startsWith(searchKey))
                    foundKey = true;

		while ( (eachLine = br.readLine()) != null && eachLine.startsWith("P\t")) {

                    if(eachLine.startsWith(searchKey))
                        foundKey = true;

		    proteinList.add(eachLine);
		}

		peptideList.add(eachLine);

		while ( (eachLine = br.readLine()) != null && eachLine.startsWith("S\t")) {
		    peptideList.add(eachLine);
		}


                if(foundKey) {
                  //  System.out.println("aaaa" + eachLine);
                  //  System.out.println("aaaa" + proteinList);

                    intArr = getCorrectFactorByLocus(indexList, peptideList, proteinList, expSize, pepNormRatioSumArr, numIdPerReplicates);

                    break;
                }
		//proteinCalc(indexList, peptideList, proteinList, correctArr, p, expSize, pepNormRatioSumArr, numIdPerReplicates);
		proteinList.clear();
		peptideList.clear();

		if(null != eachLine && eachLine.startsWith("P\t")) {

                    if(eachLine.startsWith(searchKey))
                        foundKey = true;

		    proteinList.add(eachLine); // + "<<"  );
                }

	    } else if(eachLine.startsWith("S\t")) {
		peptideList.add(eachLine);
		while ( (eachLine = br.readLine()) != null && eachLine.startsWith("S\t")) {
		    peptideList.add(eachLine);
		}


                if(foundKey) {
                    //System.out.println("aaaa" + eachLine);
                    //System.out.println("aaaa" + proteinList);
                    intArr = getCorrectFactorByLocus(indexList, peptideList, proteinList, expSize, pepNormRatioSumArr, numIdPerReplicates);
                    break;
                }
		//proteinCalc(indexList, peptideList, proteinList, correctArr, p, expSize, pepNormRatioSumArr, numIdPerReplicates);

		proteinList.clear();
		peptideList.clear();

		if(null != eachLine && eachLine.startsWith("P\t")) {
                    if(eachLine.startsWith(searchKey))
                        foundKey = true;

                    proteinList.add(eachLine);

                }
		    
	    }
	}

        
	if(null != p)
	    p.close();

//for(double d:intArr)
//System.out.println("==" + d);

        normalize(outputFilename, inputFilename, intArr, numIdPerReplicates, sampleSize);

        //intArr : correctArr

        //normalizeCorrectArr(intArr);

    }

    /*
    private static double[] normalizeCorrectArr(double[] inputArr) {

        double[] tempArr = inputArr;
        if(inputArr.length<=3)
            return inputArr;
        
        //double[] filteredArr = inputArr;
        int currentSize = tempArr.length;

        //System.out.println(inputArr.length);

        //for(int i=0;i<inputArr.length;i++)
        //System.out.println("---" + i+ "\t" +inputArr[i]);
        int count=0;
        while(true)
        {
                //dArr = edu.scripps.pms.stats.GrubbsTest.filter(tmpArr, 0.1);
            tempArr = edu.scripps.pms.stats.GrubbsTest.filterExcludingNegative(tempArr, 0.1);

            if(currentSize == tempArr.length)
                break;

            currentSize = tempArr.length;

            count++;

            if(count>10)
                break;
        }


//System.out.println(inputArr.length);


        count=0;
        double intSum=0;
        for(double d:tempArr) {
            if(d<=0)
                continue;
            
            intSum += d;
            count++;
        }
        
        double average = intSum/count;

        System.out.println(intSum);
        System.out.println(average);

        double[] correctFactorArr = new double[tempArr.length];
        for(int i=0;i<tempArr.length;i++) {
            if(tempArr[i]<=0)
                correctFactorArr[i] = 0;

            correctFactorArr[i] = inputArr[i]/average;
        }

        for(int i=0;i<inputArr.length;i++)
            System.out.println("==" + correctFactorArr[i]+ "\t" + tempArr[i] + "\t" + inputArr[i]);
        //for(double d:correctFactorArr)
            
        
        return null;
    }*/


    //we are not using this method for now
    public double[] getNormalizeFactor(String fileName, String locus) throws IOException {
        

        BufferedReader br = new BufferedReader(new FileReader(fileName));            

	sampleList = new ArrayList<Sample>();

        String eachLine; // = br.readLine();
        int experimentNum=0;
	while( (eachLine = br.readLine())!=null )
	{
            if(!eachLine.startsWith("H\tGROUP_SAMPLE"))
                break;

            String[] headArr = eachLine.split("\t");
            int tmpNum = headArr.length-3;
            experimentNum += tmpNum;


	    Sample s = new Sample(tmpNum);
            sampleList.add(s);
        }	
    
        System.out.println(eachLine);

        
        
	int[] intensityIndexArr = getIndexArrByPrefix(experimentNum, eachLine, "NORM_INTENSITY");

        String searchKey = "P\t" + locus;
	while( (eachLine = br.readLine())!=null && !eachLine.startsWith(searchKey) );

System.out.println("==" + eachLine);
System.out.println("==" + experimentNum);
        String[] arr = eachLine.split("\t");
        System.out.println(arr.length);
        for(int i:intensityIndexArr) {
            System.out.println(i + "\t" + arr[i]);
        }


        System.out.println(eachLine);
	double[] normIntArr = new double[experimentNum];
	String[] targetProteinArr = eachLine.split("\t");
	double intSum = 0;
	for(int i:intensityIndexArr) {
		double normInt = Double.parseDouble(targetProteinArr[i]);
		intSum += normInt;

	}
	double avgInt = intSum/experimentNum;
		System.out.println(avgInt + " " + intSum);


	double[] correctArr = new double[experimentNum];




	return correctArr;

    }

    public NormalizeLabelFreeIntensityBase() {
    }


    public static void normalize(String outputFilename, String fileName, double[] correctArr, int numIdPerReplicates) throws IOException {
	normalize(outputFilename, fileName, correctArr, numIdPerReplicates, sampleList.size());
    }

    public static void normalize(String outputFilename, String fileName, double[] correctArr, int numIdPerReplicates, int sampleSize) throws IOException {

	int expSize = correctArr.length;
	BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(outputFilename));
	PrintStream p = new PrintStream(out);
	p.print("H\tNormalization factor\t");
	for(double d:correctArr) {
	    p.print(d);
	    p.print("\t");
	}
	p.println("");

        BufferedReader br = new BufferedReader(new FileReader(fileName));
//System.out.println("==========" + fileName);

        String eachLine = "";

	while( (eachLine = br.readLine()) != null ) {
	
	    if(eachLine.startsWith("PLINE")) 
		break;

	    p.println(eachLine);
	}

        p.print(eachLine);
        
        for(int i=0;i<sampleSize;i++) {
            p.print("SAMPLE_GROUP_INTENSITY_AVG_");
            p.print( (i+1) );
            p.print("\t");
        }

        for(int i=0;i<expSize;i++) {
            p.print("NORM_INTENSITY_CORRECT_");
            p.print( (i+1) );
            p.print("\t");
	}

        //for(int i=0;i<sampleList.size();i++) {
        for(int i=0;i<sampleSize;i++) {
            p.print("SAMPLE_GROUP_INTENSITY_CORRECT_AVG_");
            p.print( (i+1) );
            p.print("\t");
        }

        p.println("");

	//sline
	eachLine = br.readLine();
	p.println(eachLine);


	gnu.trove.TIntArrayList intIndexList = new gnu.trove.TIntArrayList();
	String[] tmpArr = eachLine.split("\t");
	for(int i=0;i<tmpArr.length;i++) {
	    if(tmpArr[i].startsWith("INTENSITY"))
		intIndexList.add(i);
	}

	int[] indexList = intIndexList.toNativeArray();

	List<String> proteinList=new ArrayList<String>();
	List<String> peptideList=new ArrayList<String>();

	double[] pepNormRatioSumArr = new double[expSize];

//	boolean prevProtein=false;

//	while ( (eachLine = br.readLine()) != null && eachLine.startsWith("P\t"));

	while( (eachLine = br.readLine())!=null )
	{

	    if(eachLine.startsWith("P\t")) { //protein line
		//if(!prevProtein) {
//System.out.println("==========" + eachLine);

		proteinList.add(eachLine);
		while ( (eachLine = br.readLine()) != null && eachLine.startsWith("P\t")) {
		    proteinList.add(eachLine);
		}

		peptideList.add(eachLine);

                
		while ( (eachLine = br.readLine()) != null && eachLine.startsWith("S\t")) {
		    peptideList.add(eachLine);
		}

             //   System.out.println("==========" + eachLine);
                

		//List<String> newPepList = pepCalc(peptideList, pepNormRatioSumArr, correctArr, expSize);
		//proteinCalc(indexList, newPepList, proteinList, correctArr, p, expSize, pepNormRatioSumArr);

            //    if(null == peptideList)
                    proteinCalc(indexList, peptideList, proteinList, correctArr, p, expSize, pepNormRatioSumArr, numIdPerReplicates);
		proteinList.clear();
		peptideList.clear();

		if(null != eachLine && eachLine.startsWith("P\t"))
		    proteinList.add(eachLine);

	    } else if(eachLine.startsWith("S\t")) {
		peptideList.add(eachLine);
		while ( (eachLine = br.readLine()) != null && eachLine.startsWith("S\t")) {
		    peptideList.add(eachLine);
		}

		//List<String> newPepList = pepCalc(peptideList, pepNormRatioSumArr, correctArr, expSize);
		//proteinCalc(indexList, newPepList, proteinList, correctArr, p, expSize, pepNormRatioSumArr);
                proteinCalc(indexList, peptideList, proteinList, correctArr, p, expSize, pepNormRatioSumArr, numIdPerReplicates);

		proteinList.clear();
		peptideList.clear();

		if(null != eachLine && eachLine.startsWith("P\t"))
		    proteinList.add(eachLine);
	    }
	}

	if(null != p)
	    p.close();
    }


    private static double[] getCorrectFactorByLocus (
	int[] indexList,
	List<String> peptideList,
	List<String> proteinList,	
	int expSize,
	double[] pepNormRatioSumArr,
        int numIdPerReplicates
	) {
	
	double[] avgIntensity = new double[indexList.length];

	int pepCount = 0;

	for(Iterator<String> pepItr=peptideList.iterator(); pepItr.hasNext(); ) {
	    String pepLine = pepItr.next();

	    //sb.append(pepLine).append("\n");

	    String[] peplineArr = pepLine.split("\t");
            boolean addToProtein = true;

            int tmpCount=peplineArr.length-expSize;
            for(Iterator<Sample> itr=sampleList.iterator(); itr.hasNext(); ) {
                Sample s = itr.next();

                int repCount = s.getReplicateCount();

		int numValidId=0;
                for(int i=0;i<repCount;i++) {
                    String eachv = peplineArr[tmpCount++];

                    if(!"0".equals(eachv)) {
			numValidId++;
                      //  break;
                    }
                }

		if(numValidId<numIdPerReplicates) {
                    addToProtein = false;
                    break;
		}
            }

            if(!addToProtein)
                continue;

	    for(int i=0;i<expSize;i++) {
		String eachv = peplineArr[peplineArr.length-expSize+i];

		if(!"NA".equals(eachv)) {

		    double d = Double.parseDouble(eachv);
		    avgIntensity[i] += d;
		}
	    }

            pepCount++;
	}

	//StringBuffer proteinRatioSb = new StringBuffer();
        TDoubleArrayList intArr = new TDoubleArrayList();




	if(pepCount<=0) {

	    for(int i=0;i<avgIntensity.length;i++) {
                intArr.add(0.0);
	    }

            for(Iterator<Sample> itr=sampleList.iterator(); itr.hasNext(); ) {
                Sample s = itr.next();

                intArr.add(0.0);
          
            }
	}
	else {
	    for(int i=0;i<avgIntensity.length;i++) {

		avgIntensity[i] /= pepCount;

		intArr.add(avgIntensity[i]);

	    }
        }


        return intArr.toNativeArray();
    }


    private static void proteinCalc(
	int[] indexList, 
	List<String> peptideList, 
	List<String> proteinList, 
	double[] correctArr, 
	PrintStream p, 
	int expSize, 
	double[] pepNormRatioSumArr,
	int numIdPerReplicates
	) {

	StringBuffer sb = new StringBuffer();
	double[] avgIntensity = new double[indexList.length];

	int pepCount = 0;
        

	for(Iterator<String> pepItr=peptideList.iterator(); pepItr.hasNext(); ) {
	    String pepLine = pepItr.next();
	    sb.append(pepLine).append("\n");

	    String[] peplineArr = pepLine.split("\t");

            boolean addToProtein = true;

            int tmpCount=peplineArr.length-expSize;


            for(Iterator<Sample> itr=sampleList.iterator(); itr.hasNext(); ) {
                Sample s = itr.next();

                int repCount = s.getReplicateCount();
                boolean sampleValidation=false;

		int numValidId=0;
                for(int i=0;i<repCount;i++) {
                    String eachv = peplineArr[tmpCount++];

	    /*
                    if(!"0".equals(eachv)) {
                        sampleValidation = true;
                      //  break;
                    }
		    */
                    if(!"0".equals(eachv)) {
			numValidId++;
                      //  break;
                    }
                }


		if(numValidId<numIdPerReplicates) {
                    addToProtein = false;
                    break;
		}

                /*if(!sampleValidation) {
                    addToProtein = false;
                    break;
                }*/

            }
             
            if(!addToProtein)
                continue;

	    for(int i=0;i<expSize;i++) {                
		String eachv = peplineArr[peplineArr.length-expSize+i];                              

		if(!"NA".equals(eachv) && !"NaN".equals(eachv)) {

		    double d = Double.parseDouble(eachv);
		    avgIntensity[i] += d;                    
		}
	    }

            pepCount++;
	}

	StringBuffer proteinRatioSb = new StringBuffer();
	StringBuffer normalizedSb = new StringBuffer();

	//System.out.println("===" + proteinList);
	if(pepCount<=0) {
	    for(int i=0;i<avgIntensity.length;i++) {
		proteinRatioSb.append("\tNA");
                normalizedSb.append("NA\t");
	    }

            proteinRatioSb.append("\t");

            for(Iterator<Sample> itr=sampleList.iterator(); itr.hasNext(); ) {
                Sample s = itr.next();

                proteinRatioSb.append("NA\t");
                normalizedSb.append("NA\t");
            }
	}
	else {
	    for(int i=0;i<avgIntensity.length;i++) {

		avgIntensity[i] /= pepCount;
                
		proteinRatioSb.append("\t").append( Formatter.formatFourDecimal(avgIntensity[i]) );


		if(avgIntensity[i]<=0 || avgIntensity[i]>=correctArr.length)
			normalizedSb.append( Formatter.formatFourDecimal((avgIntensity[i])) ).append("\t");
		else {
			//double tmpd = avgIntensity[i]/correctArr[i];
			double tmpd = avgIntensity[i]+correctArr[i];
			//if(tmpd<=0 || tmpd>=correctArr.length || avgIntensity[i]==0)
			//if(tmpd<=0 || avgIntensity[i]==0)
			//normalizedSb.append( "==>>" + tmpd + " " + correctArr[i] + " " + avgIntensity[i] + "--++" + Formatter.formatFourDecimal(tempd) ).append("\t");
			normalizedSb.append( Formatter.formatFourDecimal(tmpd) ).append("\t");
		}
	    }

            proteinRatioSb.append("\t");

            int tmpCount=0;
            //for(Iterator<Sample> itr=sampleList.iterator(); itr.hasNext(); ) {
	    for(int i=0;i<sampleList.size();i++) {
                //Sample s = itr.next();
                Sample s = sampleList.get(i); 

                int repCount = s.getReplicateCount();

                double tmpSum1=0;
                double tmpSum2=0;
                for(int j=0;j<repCount;j++) {

                    tmpSum1 += avgIntensity[tmpCount];

		    double tmpd = avgIntensity[tmpCount]/correctArr[tmpCount];
		    if(tmpd<=0 || tmpd>=correctArr.length || avgIntensity[tmpCount]==0) {
			    tmpSum2 += avgIntensity[tmpCount];
		    } else {
			    tmpSum2 += tmpd;
		    }
                    tmpCount++;

                }

                tmpSum1 /= repCount;
		tmpSum2 /= repCount;

                proteinRatioSb.append( Formatter.formatFourDecimal(tmpSum1) ).append("\t");
		normalizedSb.append( Formatter.formatFourDecimal(tmpSum2) ).append("\t");
            }
        }

	/*
	StringBuffer normProSb = new StringBuffer();
	for(double d:pepNormRatioSumArr) {
	    normProSb.append("\t").append(d/peptideList.size());
	}
	*/


	for(Iterator<String> proItr=proteinList.iterator(); proItr.hasNext(); ) {
	    p.print(proItr.next());
	    p.print(proteinRatioSb.toString());
	    p.println(normalizedSb.toString());
	//System.out.println("===" + proteinRatioSb.toString());
	//System.out.println("---" + normalizedSb.toString());
	}

	p.print(sb.toString());
    }

    private static List<String> pepCalc(List<String> peptideList, double[] pepNormRatioSumArr, double[] correctArr, int expNum) {

	List<String> newPepList = new ArrayList<String>();

	for(Iterator<String> itr=peptideList.iterator(); itr.hasNext(); ) {

	    String eachLine = itr.next();
            
            
	    String[] arr = eachLine.split("\t");

	    int intIndex=6;
	    int count=0;
	    double totalIntSum = 0;

	   // for(int i=intIndex;i<=arr.length-expNum;i+=intIndex) {

	    double[] intArr = new double[expNum];

	    for(int i=0;i<expNum;i++) {

		int tmpIndex = intIndex + i*intIndex;

		if("NA".equals(arr[tmpIndex])) {
		    eachLine += "NA\t";
		    count++;
		    continue;
		}

		Double d = Double.parseDouble(arr[tmpIndex]);
		double tmpValue = Math.log(d)/Math.log(2);

		tmpValue += correctArr[count];

		double tmpPow = Math.pow(2,tmpValue);
		eachLine += tmpPow; 
		eachLine += "\t";
		intArr[i] = tmpPow;

		pepNormRatioSumArr[count] += tmpPow;

		totalIntSum += tmpPow; 

		count++;
	    }

	    double averageIntensity = totalIntSum/expNum;

	    for(int i=0;i<expNum;i++) {
		double pepRatio = intArr[i]/averageIntensity;
		eachLine += pepRatio;
		eachLine += "\t";
	    }


	    newPepList.add(eachLine);
	}

	return newPepList;
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
}
