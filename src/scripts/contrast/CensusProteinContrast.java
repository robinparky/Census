/*
* Copyright (c) 2008 Proteomics Integrated Solutions. All rights reserved.  
*/


import java.util.*;
import java.io.*;

public class CensusProteinContrast
{

    private static Hashtable<String, Set> descHt = new Hashtable<String, Set>();

    public static Hashtable<String, String> readFile(String fileName, boolean includeRedundant) throws Exception
    {
	Hashtable<String, String> proHt = new Hashtable<String, String>();

        BufferedReader br = new BufferedReader(new FileReader(fileName));

        String eachLine;

	Hashtable<String, Set> ht = new Hashtable<String, Set>();

	String prevLine="";

	while( (eachLine = br.readLine()) != null)
	{       
	    String curLine = eachLine;

	    if(!includeRedundant && prevLine.startsWith("P\t")) {
		prevLine = curLine;
		continue;
	    }

	    if( eachLine.startsWith("P\t") )
	    {
		String[] arr = eachLine.split("\t");                

		//proHt.put(arr[1] + "\t" + arr[2], eachLine);
		proHt.put(arr[1], eachLine);
	    }            

	    prevLine = curLine;
	}

	return proHt;
    }

    public static void main(String[] args) throws IOException, Exception
    {       
/*	if(args.length<=0)
	{
	    System.out.println("it needs census-out.txt files as paramters");
	    System.exit(0);
	}
*/

	CensusProteinContrast contrast = new CensusProteinContrast();

	String[] input = new String[2];
	//input[0] = "/home/rpark/ip2_tomcat/webapps/ip2/ip2_data/jhprieto/Fansidar/fansidar_4th_rep1_2009_01_15_18_36/quant/2009_02_15_01_107/census-out.txt";
	//input[1] = "/home/rpark/ip2_tomcat/webapps/ip2/ip2_data/jhprieto/Fansidar/fansidar_4th_rep1_2009_01_15_18_36/quant/2009_02_15_17_124/census-out.txt";
	input[0] = "/home/rpark/Dropbox/TSRI/collaborators/AnChi/HDLD_reverse.txt";
	input[1] = "/home/rpark/Dropbox/TSRI/collaborators/AnChi/HDLL_Forward.txt";

	contrast.run(input, true);
    }

    public void run(String[] inputFileArr, boolean includeRedundant) throws Exception {
	run(null, inputFileArr, null, "", includeRedundant);
    }
    
    public void run(String[] nameArr, String[] inputFileArr, String resultPath, String suffix, boolean includeRedundant) throws Exception {
	if(null == resultPath)
	    resultPath = ".";
	
	if(null == nameArr)
	    nameArr = inputFileArr;

	String outputFile = resultPath + File.separator + "quant_compare" + suffix + ".txt";
        PrintStream out = new PrintStream(new FileOutputStream(outputFile));
        
	List<Hashtable<String, String>> proHtList = new ArrayList<Hashtable<String, String>>();
//	List<Hashtable<String, Set>> setList = new ArrayList<Hashtable<String, Set>>();
    
	StringBuffer sbHeader = new StringBuffer();
	sbHeader.append("PLINE").append("\t");
	sbHeader.append("LOCUS").append("\t");

	for(int i=0;i<nameArr.length;i++)
	{
	    out.print(inputFileArr[i] + "\t");

	    int count=i+1;
	    sbHeader.append("AVERAGE_RATIO_").append(count).append("\t");
	    sbHeader.append("STANDARD_DEVIATION_").append(count).append("\t");
	    sbHeader.append("WEIGHTED_AVERAGE_").append(count).append("\t");
	    sbHeader.append("PEPTIDE_NUM_").append(count).append("\t");
	    sbHeader.append("SPEC_COUNT_").append(count).append("\t");
	    sbHeader.append("AREA_RATIO_").append(count).append("\t");
	}

	out.println("");

	sbHeader.append("DESCRIPTION").append("\n");
	out.println(sbHeader.toString());


	for(int i=0;i<inputFileArr.length;i++)
	{
	    Hashtable<String, String> proHt = readFile(inputFileArr[i], includeRedundant);
	    proHtList.add(proHt);
	}

	Set<String> analyzedSet = new HashSet<String>();
	
	for(int i=0;i<proHtList.size();i++)
	{
	    Hashtable<String, String> proHt = proHtList.get(i);
	    Set s = proHt.keySet();

	    int countPro = 0;
	    int colSize=0;
	    int maxColSize=0;
	    for(Iterator<String> itr=s.iterator(); itr.hasNext(); )
	    {
		String key = itr.next();
		String str = proHt.get(key);

		if(analyzedSet.contains(key))
		    continue;

		String[] arr = str.split("\t");
		colSize = arr.length;

		if(colSize>maxColSize)
		    maxColSize = colSize;
	    }

	    for(Iterator<String> itr=s.iterator(); itr.hasNext(); )
	    {
		String key = itr.next();
		String str = proHt.get(key);

		if(analyzedSet.contains(key))
		    continue;

		analyzedSet.add(key);

		String[] arr = str.split("\t");

		out.print("P\t"); out.print(arr[1]); out.print("\t");

		for(int tmpC=0;tmpC<i;tmpC++)
		{
		    //for(String eachStr:arr)
		    for(int tmpC2=2; tmpC2<maxColSize-1; tmpC2++)
			out.print("X\t");
		}

		for(int tmpC=2; tmpC<arr.length-1; tmpC++) {
		    out.print(arr[tmpC]); out.print("\t");
		}
		//out.print(str + "\t");



		for(int tmpC=i+1;tmpC<proHtList.size();tmpC++)
		{
		    Hashtable<String, String> tmpPepHt = proHtList.get(tmpC);
		    String tValue = tmpPepHt.get(key);

		//    tmpPepHt.remove(key);
//		    Hashtable<String, Set> tht = setList.get(tmpC);
//		    Set tset = tht.get(key);

		    //out.print(((tset==null)?"X\tX\tX\tX\tX":tmpPepHt.get(key)) + "\t");

		    if(tValue == null)
		    {
			//for(String eachStr:arr)
			for(int tmpC2=2; tmpC2<maxColSize-1; tmpC2++)
			    out.print("X\t");
		    }
		    else
		    {
			String[] tmpArr = tValue.split("\t");
			for(int j=2; j<tmpArr.length-1; j++)  {
			    out.print(tmpArr[j]); out.print("\t");
			}
//			out.print(tValue + "\t");
		    }
		}

		out.print(arr[arr.length-1]); out.println("");
	    }
	}

        out.close();


    }

    class MergeProteinModel {

	private Set<String> proteinNames = new HashSet<String>();
	private Hashtable<String, String> proteinDescHt = new Hashtable<String, String>();

	public MergeProteinModel() {
	}

	public Hashtable<String, String> getProteinDescHt()
	{
	    return proteinDescHt;
	}

	public String getProteinDesc(String proteinName)
	{
	    return proteinDescHt.get(proteinName);
	}

	public Set getProteinNames() {
	    return proteinNames;
	}

    }
}

