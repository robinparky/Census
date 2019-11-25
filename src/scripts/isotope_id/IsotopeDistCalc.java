package scripts.isotope_id;

import java.util.*;
import java.io.*;

import edu.scripps.pms.census.io.*;
import edu.scripps.pms.census.hash.*;
import edu.scripps.pms.census.*;

import java.io.BufferedReader;
import java.util.*;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.FileReader;
import java.io.RandomAccessFile;

import edu.scripps.pms.census.util.dtaselect.Protein;
import edu.scripps.pms.census.util.dtaselect.Peptide;

import edu.scripps.pms.census.model.*;
import edu.scripps.pms.census.ElementComposition;

import org.jdom.*;
import org.jdom.output.*;
import org.jdom.input.*;
import edu.scripps.pms.census.util.*;
import edu.scripps.pms.census.conf.*;

import gnu.trove.*;


public class IsotopeDistCalc 
{

    private final static double PROTON_MASS = 1.00728;

    public static void main(String args[]) throws Exception
    {
	if(args.length<2)
	{
	    System.out.println("Usage: java IsotopeDistCalc path DTASelect-filter.txt");
	    System.exit(0);
	}
	//String path = "/home/rpark/pms/Census/data";
	String path = args[0];
	Configuration conf = Configuration.getInstance();
	conf.setSpectrumFormat(0);
	ChroGenerator.setConfiguration(conf);
	Hashtable<String, IndexedFile> ht = ChroGenerator.createIndexedFiles(path, CensusConstants.MS1_FILE);

	//Convert DTASelect to pepXML

        BufferedReader br = new BufferedReader(new FileReader(args[1]));            
	String eachLine;
	while( !(eachLine = br.readLine()).startsWith("Unique") );

	while( (eachLine = br.readLine())!=null )
	{
	    String[] arr = eachLine.split("\t");

	    if(arr.length<11)
		continue;

	    String peptide = arr[arr.length-1];


	    String seq = peptide.substring(2, peptide.length()-2);
	    String csStr = arr[1];
	    int cs = Integer.parseInt(String.valueOf(csStr.charAt(csStr.length()-1)));

	    try {
		String scan = arr[1].substring(arr[1].indexOf(".")+1);
		scan = scan.substring(0, scan.indexOf("."));

		String fileName = arr[1].substring(0, arr[1].indexOf("."));
		fileName = path + File.separator + fileName;
		fileName += ".ms1";

		IndexedFile iFile = ht.get(fileName);

		RandomAccessFile rfile = iFile.getFile();
		TIntLongHashMap indexMap = iFile.getMsIndex();

		int[] keys = iFile.getKeys();
		int keyIndex = Arrays.binarySearch(keys, Integer.parseInt(scan));

		if(keyIndex<0) //Cannot find index
		    keyIndex=-(++keyIndex); //Math.abs(++keyIndex);

		if(keyIndex>=keys.length)
		    keyIndex--;

		double[][] calcArr = getIsotopeFile(seq);
		double[] samIsoArr = calcArr[0];

		for(int i=0;i<samIsoArr.length;i++) {
		    samIsoArr[i] = (samIsoArr[i]+cs*PROTON_MASS)/cs;
		}

		
		double[] relabun = calcArr[1];


		double[][] massIntArr  = CalcUtil.readSpectrumOnly(keyIndex, iFile);
		double[] massArr = massIntArr[0];
		double[] intArr = massIntArr[1];
//	for(double d : samIsoArr)
//	    System.out.println( "==" + d);



		double[] specIntArr = CalcUtil.intensityArr(massArr, intArr, samIsoArr, 5.0);
		double[] smallRelabun = new double[specIntArr.length];
/*
		for(int i=0;i<massArr.length;i++)
		    System.out.println(massArr[i] + "\t" + intArr[i]);

		for(int i=0;i<samIsoArr.length;i++)
		    System.out.println(samIsoArr[i]);


*/
		for(int i=0;i<specIntArr.length;i++)
		    smallRelabun[i] = relabun[i];
		LinearRegressionDouble reg = new LinearRegressionDouble(smallRelabun, specIntArr, 0, specIntArr.length, 0);
		double slope = reg.getSlope();
		double intercept = reg.getIntercept();
		double corr = reg.getCorr();

		System.out.println("S\t" + peptide + "\t" + scan + "\t" + slope + "\t" + intercept + "\t" + corr);
		for(int i=0;i<smallRelabun.length;i++)
		    System.out.println("I\t" + samIsoArr[i] + "\t" + specIntArr[i] + "\t" + smallRelabun[i]);


//	for(double d : specIntArr)
//	    System.out.println( d);
		    

	    } catch (Exception e) {
		System.out.println(eachLine);
		System.out.println(seq + " " + e.toString());
		e.printStackTrace();
	    }
	}

    }

    public static double[][] getIsotopeFile(String seq) throws Exception {

        SAXBuilder sb = new SAXBuilder();
	Document doc = sb.build("/data/1/rpark/deploy/census_config_repository/census_config_silac.xml");
	Element root = doc.getRootElement();

	IsotopeReader isoReader = new IsotopeReader(root);
	char[] ch = seq.toCharArray();


	ElementComposition element = new ElementComposition(ch, 0, ch.length, isoReader.getIsotope());
	element.calculate();

	IsotopeDist sampleDist = new IsotopeDist(element.getElementSampleArr(), element.getModShift(), true);

	double[] samIsoArr = sampleDist.getHighMassList();
	double[] relabun = sampleDist.getRelabun();

	double[][] calcArr = new double[2][];
	calcArr[0] = samIsoArr;
	calcArr[1] = relabun;

	return calcArr;

/*
	for(double d : samIsoArr)
	    System.out.println( d);

	System.out.println( "===" );
	for(double d : relabun)
	    System.out.println( d);
	System.out.println( "===" );
*/

    }
}
