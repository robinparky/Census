package scripts.emily;

import org.jdom.*;
import org.jdom.output.XMLOutputter;
import org.jdom.input.*;
import java.math.BigInteger;
import java.io.*;
import java.util.*;
import java.text.DecimalFormat;
import java.text.NumberFormat;

//find all nothing from filtered_chro.xml
public class FindAllNothing
{
    private static DecimalFormat format = new DecimalFormat("0.00");
    private static DecimalFormat format2 = new DecimalFormat("0.0000");
    public static void main(String args[]) throws Exception
    {

	if(args.length<1)
	{
	    System.out.println("Usage: FindAllNothing filtered_chor.xml");
	    ///data/2/gcantin/lee/wei_yi/SILAC_phos_expers_112006/L-EGF_H_112206/light_parc/all_nothing_candidates_03_20_chro.xml

	    System.exit(0);
	}

	System.out.println("H\tAll and Nothing peptides analysis result");
	System.out.println("H\tSung Kyu (Robin) Park rpark@scripps.edu");
	System.out.println("H\tlocus\tseq_ct\tspec_ct\tseq_cov\tlength\tmolwt\tpi\tval\tdesc");
	System.out.println("H\tseq\tcharge\tlightSum\theavySum\tratio\tR\tR^2");
	SAXBuilder builder = new SAXBuilder();

	Document doc = builder.build(new File(args[0]));
//	String skipWord = args[2];

	//Document doc = builder.build(new File("/data/2/gcantin/lee/wei_yi/SILAC_phos_expers_112006/L-EGF_H_112206/light_parc/all_nothing_candidates_03_20_chro.xml"));
	Element rootEle = doc.getRootElement();

	StringBuffer proSb = new StringBuffer();

	for(Iterator<Element> itr=rootEle.getChildren("protein").iterator(); itr.hasNext(); )
	{
	    Element proEle = itr.next();

//	    if((null != skipWord && proEle.getAttributeValue("desc").contains(skipWord)) || proEle.getAttributeValue("locus").startsWith("Reve"))
//		continue;

	    proSb.append("P").append("\t");
	    proSb.append(proEle.getAttributeValue("locus")).append("\t");
	    proSb.append(proEle.getAttributeValue("seq_ct")).append("\t");
	    proSb.append(proEle.getAttributeValue("spec_ct")).append("\t");
	    proSb.append(proEle.getAttributeValue("seq_cov")).append("\t");
	    proSb.append(proEle.getAttributeValue("length")).append("\t");
	    proSb.append(proEle.getAttributeValue("molwt")).append("\t");
	    proSb.append(proEle.getAttributeValue("pi")).append("\t");
	    proSb.append(proEle.getAttributeValue("val")).append("\t");
	    proSb.append(proEle.getAttributeValue("desc")).append("\n");

	    StringBuffer pepSb = new StringBuffer();
	    for(Iterator<Element> pepitr=proEle.getChildren("peptide").iterator(); pepitr.hasNext(); )
	    {
		Element pepEle = pepitr.next();
		String[] chroArr = pepEle.getChildText("chro").split(";");

		String[] peakRangeArr = chroArr[0].split(" ");
		int peakStart = Integer.parseInt(peakRangeArr[1]);
		int peakEnd = Integer.parseInt(peakRangeArr[2]);

		ArrayList<Long> lightArr = new ArrayList<Long>();
		ArrayList<Long> heavyArr = new ArrayList<Long>();

		long lightSum = 0;	
		long heavySum = 0;	
                double massDiff=0;
                int peakNoneZeroCount=0;

		for(int i=1;i<chroArr.length;i++)
		{
		    String[] tmpArr = chroArr[i].split(" ");

		    int scanNum = Integer.parseInt(tmpArr[0]);

		    if(scanNum<peakStart || scanNum>peakEnd)
			continue;

		    long lightValue = Long.parseLong(tmpArr[1]);
		    long heavyValue = Long.parseLong(tmpArr[2]);

		    lightSum += lightValue;
		    heavySum += heavyValue; 

		    lightArr.add(lightValue);
		    heavyArr.add(heavyValue);

		    double d = Double.parseDouble(tmpArr[7]);
		    if(d>=0) //if tolerance is less than zero (-1.0), no peaks are found
		    {
			peakNoneZeroCount++;
			massDiff += d;
		    }

		    d = Double.parseDouble(tmpArr[8]);
		    if(d>=0) //if tolerance is less than zero (-1.0), no peaks are found
		    {
			peakNoneZeroCount++;
			massDiff += d;
		    }
		}

		if(peakNoneZeroCount<=0)
		    massDiff = -1;
		else
		    massDiff = massDiff/peakNoneZeroCount;

		double average = 0;
		double stdev = 0;
		double max = 0;
		double[] arr = null;
		double[] normArr = null;

		if(lightSum > heavySum)
		{
		    average = getCOG(lightArr);

		    arr = new double[lightArr.size()];
		    for(int i=0;i<lightArr.size();i++)
		    {
			arr[i] = lightArr.get(i);
		    }
		}
		else
		{
		    average = getCOG(heavyArr);

		    arr = new double[heavyArr.size()];
		    for(int i=0;i<heavyArr.size();i++)
		    {
			arr[i] = heavyArr.get(i);
		    }
		}
//		arr = normalizeOnArea(arr, normArr);
/*
		String seq = pepEle.getAttributeValue("seq");
		String ch = pepEle.getAttributeValue("charge");
		String sc = pepEle.getAttributeValue("scan");
		System.out.println("===" + sc);
		if(seq.equals("K.VPNACLFTINK.E") && ch.equals("2") && sc.equals("02381"))
		{
		    for(int ii=0;ii<arr.length;ii++)
			System.out.println(arr[ii] + "\t"); // + normArr[ii] + "\t&");
			//System.out.println(arr[ii] + "\t" + normArr[ii] + "\t&");

		    System.out.println(lightSum + "\t" + heavySum);
		    for(int ii=0;ii<lightArr.size();ii++)
			System.out.println(lightArr.get(ii) + "\t" + heavyArr.get(ii) + "\t*");

		}
		System.out.println("===");
*/

		stdev = getStandardDeviationFromDist(arr, average);
//		normArr = getNormalDistFrequency(getXArr(arr), average, stdev);

		/*************************  Calculate Pearson's chi-square test *******************/
		//chi test does not work well on all-nothing case
    //public static double[] getNormalDistFrequency(double[] arr, double mean, double sigma)
//		normArr = changeFrequency(normArr, getMax(arr));		

//		getChiScore(arr, normArr);

		/*************************  End of Calculate Pearson's chi-square test ************/

		arr = getNormalize(arr);


		double[] xArr = getXArr(arr);
		normArr = getNormalDistFrequency(xArr, average, stdev);

		LinearRegression2 reg = new LinearRegression2(arr, normArr, 0, arr.length, 0);
		double corr = reg.getCorr();
		double threshold = Double.parseDouble(args[1]);
		if(threshold < (corr*corr) && corr>0)
		{
//		pepSb.append(threshold + " " + (corr*corr));
		    if(heavySum==0 || (lightSum/heavySum>5.0 || lightSum/heavySum<0.2))
		    {
			pepSb.append("S").append("\t").append(pepEle.getAttributeValue("seq") + "\t" + pepEle.getAttributeValue("charge"));
			pepSb.append("\t" + lightSum);
			pepSb.append("\t" + heavySum);
			pepSb.append("\t" + ((heavySum>0)?format.format((double)lightSum/heavySum):"NA"));
			pepSb.append("\t" + format2.format(massDiff));
			pepSb.append("\t" + format.format(reg.getCorr()));
			pepSb.append("\t" + ((reg.getCorr()>0)?(format.format(reg.getCorr()*reg.getCorr())):"0.00")).append("\n");
		    }
		}
	    }

	    if(proEle.getChildren("peptide").size()>0)
	    {
		System.out.print(proSb.toString());
		System.out.print(pepSb.toString());
		proSb = new StringBuffer();
	    }
	}
    }

    public static double[] normalizeOnArea(double[] arr, double[] normArr)
    {
	double normArea = 0;
	double arrArea = 0;

	for(double d:normArr)
	    normArea += d;

	for(double d:arr)
	    arrArea += d;

	if(normArea>arrArea)
	{
	    double factor = 0.0001;
	    double[] tempArr = arr;

	    while( Math.abs(normArea-arrArea) > 0.01 )
	    {
		tempArr = changeFrequency(arr, (1+factor));
		arrArea = 0;

		for(double d:tempArr)
		    arrArea += d;

		double diff = normArea - arrArea;
		factor = factor+0.0001;
//	System.out.println(normArea + " --- " + arrArea + " " +  factor);
		
	    }

	    arr = tempArr;
	}
	else
	{
	    double factor = 0.9999;
	    double[] tempArr = arr;

	    while( Math.abs(normArea-arrArea) > 0.01 )
	    {
		tempArr = changeFrequency(arr, (factor));
		arrArea = 0;

		for(double d:tempArr)
		    arrArea += d;
//	System.out.println(normArea + " - " + arrArea + " " +  factor);

		double diff = normArea - arrArea;
		factor = factor-0.0001;
		
	    }

	    arr = tempArr;
	}

	//System.out.println(normArea + " - " + arrArea);


	return arr;

    }

    public static double[] changeFrequency(double[] arr, double factor)
    {
	double[] normArr = new double[arr.length];
	for(int i=0;i<arr.length;i++)
	{
	    normArr[i] = arr[i]*factor;
	}

	return normArr;
    }
    
    public static double[] getChiScore(double[] arr, double[] expectArr)
    {
	double chiSquare=0;
	for(int i=0;i<arr.length;i++)
	{
	    chiSquare = chiSquare + (arr[i]-expectArr[i]) * (arr[i]-expectArr[i]) / expectArr[i];
	    System.out.println(arr[i] + "-" + expectArr[i]);

	}

	    System.out.println(chiSquare);
	    System.out.println(chiSquare);
	System.exit(0);

	return null;
    }
    
    //get arr of x axis starting from 0
    public static double[] getXArr(double[] arr)
    {
	double[] xArr = new double[arr.length];
	for(int j=0;j<xArr.length;j++)
	{
	    xArr[j] = j;
	}

	return xArr;
    }
   
    public static double[] getNormalDistFrequency(double[] arr, double mean, double sigma)
    {
	double[] normArr = new double[arr.length];

	for(int i=0;i<arr.length;i++)
	{
	    normArr[i] = Math.exp(-Math.pow(arr[i] - mean, 2)/(2 * Math.pow(sigma,2))) / (sigma * Math.sqrt(2 * Math.PI));
	}
	
	return getNormalize(normArr);
    }
    
    public static double[] getNormalize(double[] arr) {

	double max = getMax(arr);

	double[] normArr = new double[arr.length];

	for(int i=0;i<arr.length;i++)
	{
	    normArr[i] = arr[i]/max;
	}

	return normArr;
    }
    
    public static double getMax(double[] arr) {
	double max=0;

	for(double value: arr)
	    if(max<value)
		max = value;
	
	return max;
    }

    //get standard deviation from frequency distribution
    public static double getStandardDeviationFromDist(double[] values, double mean) {

	int num = 0;
	double sumDiffSqr = 0;
	
	for(int i=0;i<values.length;i++) { //double d : values) {

	    double diff = i - mean;
	    sumDiffSqr += (diff*diff*values[i]);
	    num += values[i];
	}
	return Math.sqrt(sumDiffSqr/(num-1));

    }
	
    //get center of gravity
    public static double getCOG(ArrayList list)
    {
	int i=0;
	long sumMass=0;
	long sumArea=0;
	for(Iterator<Long> itr=list.iterator(); itr.hasNext(); )
	{
	    Long l = itr.next();
	    sumMass += l.longValue();
	    sumArea += i*l.longValue();

	    i++;
	}

	return (double)sumArea/sumMass;
    }
}
