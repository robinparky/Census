/*
 * AllNoneUtil.java
 *
 * Created on November 1, 2007, 9:55 AM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package edu.scripps.pms.census.util;

import java.util.*;
import edu.scripps.pms.census.model.*;
import org.apache.commons.lang3.ArrayUtils;

/**
 *
 * @author rpark
 */
public class AllNoneUtil {

    /** Creates a new instance of AllNoneUtil */
    public AllNoneUtil() {
    }

    //get center of gravity
    public static double getCOG(ArrayList list) {
      int i = 0;
      long sumMass = 0;
      long sumArea = 0;
      for (Iterator<Long> itr = list.iterator(); itr.hasNext(); ) {
        Long l = itr.next();
        sumMass += l.longValue();
        sumArea += i * l.longValue();

        i++;
      }
      return (double) sumArea / sumMass;
    }
      //get center of gravity
  public static double getCOG(double[] arr)
  {

    double sumMass=0;
    double sumArea=0;
    for(int i=0;i<arr.length;i++)
    {
      double d = arr[i];
      sumMass += d;
      sumArea += i*d;

      i++;
    }

	  return sumArea/sumMass;
    }

    //get center of gravity
    public static double getCOG(long[] list)
    {
	long sumMass=0;
	long sumArea=0;
	for(int i=0;i<list.length;i++)
	{
	    long l = list[i];
	    sumMass += l;
	    sumArea += i*l;

	    i++;
	}

	return (double)sumArea/sumMass;
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

    //************** For all None analysis for Census *********************//
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

    public static void getANScore(ChroPeptide peptide) {
	List<ChroData> l = peptide.getDataList();

	int startRange = Integer.parseInt(peptide.getStartRange());
	int endRange = Integer.parseInt(peptide.getEndRange());
	getANScore(peptide, l, startRange, endRange);
    }

    public static void getANScore(ChroPeptide peptide, List<ChroData> dataList, int peakStart, int peakEnd) {

	ArrayList<Long> lightArr = new ArrayList<Long>();
	ArrayList<Long> heavyArr = new ArrayList<Long>();

	long lightSum = 0;
	long heavySum = 0;
	double massDiffLight=0;
	double massDiffHeavy=0;
	int peakNoneZeroCountLight=0;
	int peakNoneZeroCountHeavy=0;

	long maxInt = 0;
	int totalIsoSam=0;
	int totalIsoRef=0;
	int isoFoundSam=0;
	int isoFoundRef=0;

	for(Iterator<ChroData> dataItr=dataList.iterator(); dataItr.hasNext(); )
	{
	    ChroData eachData = dataItr.next();

	    totalIsoSam += eachData.getTotalIsoPeaksSam();
	    totalIsoRef += eachData.getTotalIsoPeaksRef();
	    isoFoundSam += eachData.getFoundIsoPeaksSam();
	    isoFoundRef += eachData.getFoundIsoPeaksRef();

	    //String[] tmpArr = chroArr[i].split(" ");
	    int scanNum = eachData.getScanNum();

	    if(scanNum<peakStart || scanNum>peakEnd)
		continue;

	    long lightValue = eachData.getSampleIntensity();
	    long heavyValue = eachData.getRefIntensity();

	    long[] intenArr = eachData.getIntensityArr();

	    if(maxInt<lightValue)
		maxInt = lightValue;

	    if(maxInt<heavyValue)
		maxInt = heavyValue;

	    lightSum += lightValue;
	    heavySum += heavyValue;

	    lightArr.add(lightValue);
	    heavyArr.add(heavyValue);

	    //System.out.println(lightValue + "\t" + heavyValue);
	    double d = eachData.getMassToleranceLight();
	    if(d>=0) //if tolerance is less than zero (-1.0), no peaks are found
	    {
		peakNoneZeroCountLight++;
		massDiffLight += d;
	    }

	    d = eachData.getMassToleranceHeavy();
	    if(d>=0) //if tolerance is less than zero (-1.0), no peaks are found
	    {
		peakNoneZeroCountHeavy++;
		massDiffHeavy += d;
	    }
	}

	if(peakNoneZeroCountLight<=0)
	    massDiffLight = -1;
	else
	    massDiffLight = massDiffLight/peakNoneZeroCountLight;

	if(peakNoneZeroCountHeavy<=0)
	    massDiffHeavy = -1;
	else
	    massDiffHeavy = massDiffHeavy/peakNoneZeroCountHeavy;


	double average = 0;
	double stdev = 0;
	double max = 0;
	double[] arr = null;
	double[] normArr = null;

	if(lightSum > heavySum)
	{
	    average = AllNoneUtil.getCOG(lightArr);

	    arr = new double[lightArr.size()];
	    for(int i=0;i<lightArr.size();i++)
	    {
		arr[i] = lightArr.get(i);
	    }
	}
	else
	{
	    average = AllNoneUtil.getCOG(heavyArr);

	    arr = new double[heavyArr.size()];
	    for(int i=0;i<heavyArr.size();i++)
	    {
		arr[i] = heavyArr.get(i);
	    }
	}

	stdev = AllNoneUtil.getStandardDeviationFromDist(arr, average);

	arr = AllNoneUtil.getNormalize(arr);

	int pepDataPoints = arr.length;
	for(int i=0;i<arr.length;i++)
	{
	    if(arr[i]<=0)
		pepDataPoints--;
	}

	double[] xArr = AllNoneUtil.getXArr(arr);
	normArr = AllNoneUtil.getNormalDistFrequency(xArr, average, stdev);


	LinearRegressionDouble regAN = new LinearRegressionDouble(arr, normArr, 0, arr.length, 0);
        //for(int i=0;i<arr.length;i++) {
        //    System.out.println(arr[i] + "\t" + normArr[i]);
       // }

	double corrAN = regAN.getCorr();

	double massTolerance = -1;

	if(massDiffLight>=0 && massDiffHeavy>=0)
	    massTolerance = (massDiffLight+massDiffHeavy)/(peakNoneZeroCountLight+peakNoneZeroCountHeavy);
	else if(massDiffLight<=0 && massDiffHeavy>=0)
	    massTolerance = massDiffHeavy/peakNoneZeroCountHeavy;
	else if(massDiffLight>=0 && massDiffHeavy<=0)
	    massTolerance = massDiffLight/peakNoneZeroCountLight;

	//System.out.println("===>>" + proAN.getLocus() + " " + pepAN.getSequence() + " " + corrAN + " " + pepDataPoints + " " + massTolerance);
	peptide.setCorrToNorm(corrAN);
	peptide.setSpectraDataPoints(pepDataPoints);
	peptide.setMassTolerance(massTolerance);

/*
	if(peptide.getSequence().equals("R.KTGQAAGFSYTDANK.N"))
	{
for(int i=0;i<arr.length;i++)
    System.out.println(arr[i]);
	    System.out.println( "\t" + "R.KTGQAAGFSYTDANK.N" + corrAN + "\t" + massTolerance + "\t" + massDiffLight + " " + massDiffHeavy + " " + peakNoneZeroCountLight + " " + peakNoneZeroCountHeavy);

	}

*/
    if(corrAN<0 || massTolerance<0)
        peptide.setAnCompositeScore(0.0);
    else
        peptide.setAnCompositeScore(corrAN*1.0 - massTolerance*0.026564);
  }


  public static void getANScoreTriple(ChroPeptide peptide) {

    List<ChroData> dataList = peptide.getDataList();

    int peakStart = Integer.parseInt(peptide.getStartRange());
    int peakEnd = Integer.parseInt(peptide.getEndRange());
    //getANScore(peptide, l, startRange, endRange);


    double maxInt = 0;
    double totalInt = 0;
    int count=0;

   // System.out.println("aaaaaaaaaaa");
    List<Double> intensityList = new ArrayList<>();
    for(Iterator<ChroData> dataItr=dataList.iterator(); dataItr.hasNext(); )
    {
      ChroData eachData = dataItr.next();


      //String[] tmpArr = chroArr[i].split(" ");
      int scanNum = eachData.getScanNum();

      if(scanNum<peakStart || scanNum>peakEnd)
        continue;

      //long lightValue = eachData.getSampleIntensity();
      //long heavyValue = eachData.getRefIntensity();
      double[] intArr = eachData.getdIntensityArr();
      double sum = intArr[0] + intArr[1] + intArr[2];

     //System.out.println("==============--\t" + intArr[0] + "\t" + intArr[1] + "\t" + intArr[2]);

      totalInt += sum/1000.0;

      intensityList.add(sum);
     // System.out.println(sum);
      //System.out.println(intArr[0] + "\t" + intArr[1] + "\t" + intArr[2]);

      if(maxInt<sum)
        maxInt = sum;

      count++;
    }
    /*

    if(peakNoneZeroCountLight<=0)
      massDiffLight = -1;
    else
      massDiffLight = massDiffLight/peakNoneZeroCountLight;

    if(peakNoneZeroCountHeavy<=0)
      massDiffHeavy = -1;
    else
      massDiffHeavy = massDiffHeavy/peakNoneZeroCountHeavy;
    */


    //double averageInt = totalInt/count*1000;
    Double[] intArr = intensityList.toArray(new Double[intensityList.size()]);
    double[] arr = ArrayUtils.toPrimitive(intArr);
    double averageInt = AllNoneUtil.getCOG(arr);;
    double[] normArr = null;

    double stdev = AllNoneUtil.getStandardDeviationFromDist(arr, averageInt);
    arr = AllNoneUtil.getNormalize(arr);


    int pepDataPoints = arr.length;
    for(int i=0;i<arr.length;i++)
    {
      if(arr[i]<=0)
        pepDataPoints--;
    }

    double[] xArr = AllNoneUtil.getXArr(arr);
    normArr = AllNoneUtil.getNormalDistFrequency(xArr, averageInt, stdev);
    //normArr = AllNoneUtil.getNormalDistFrequency(arr, averageInt, stdev);

    LinearRegressionDouble regAN = new LinearRegressionDouble(arr, normArr, 0, arr.length, 0);
    //for(int i=0;i<arr.length;i++) {
    //    System.out.println(arr[i] + "\t" + normArr[i]);
    // }

    double corrAN = regAN.getCorr();

    //System.out.println("===>>" + proAN.getLocus() + " " + pepAN.getSequence() + " " + corrAN + " " + pepDataPoints + " " + massTolerance);
    peptide.setCorrToNorm(corrAN);
    peptide.setSpectraDataPoints(pepDataPoints);
    peptide.setProfileScore(corrAN);

  }

  public static void main(String[] args) {
    int size = 60;
    double[] arr = new double[size];

    for(int i=0;i<size;i++)
      arr[i] = i;

    double mean = 1216714.2379310345;
    double sigma = 1216703.0856637058;

    double[] resArr = getNormalDistFrequency(arr, mean, sigma);
    System.out.println(resArr);

  }

}
