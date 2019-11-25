package edu.scripps.pms.stats;

import edu.scripps.pms.census.model.ChroPeptide;
import edu.scripps.pms.census.model.MergeProteinModel;
import edu.scripps.pms.census.tmtFilter.Peptide;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import java.util.*;
import java.io.*;
import javax.print.attribute.HashAttributeSet;
/**
 *
 * @author  Robin Park
 * @version $Id: GrubbsTest.java,v 1.20 2014/08/06 05:29:53 rpark Exp $
*/

public class GrubbsTest
{
    public static void main(String args[]) throws Exception
    {

/*
	BufferedReader br = new BufferedReader(new FileReader("/data/1/rpark/out.txt"));

	String eachLine;

	gnu.trove.TDoubleArrayList inputArr = new gnu.trove.TDoubleArrayList();

	while( (eachLine = br.readLine()) != null )
	{
//	    System.out.println(eachLine);
	    inputArr.add(Double.parseDouble(eachLine));
	}

	double[] arr = inputArr.toNativeArray();

	double stdev = STDev.getStdev(arr);
	double mean = STDev.getMean(arr);

	for(int i=0;i<arr.length;i++)
	{
	    double diff = Math.abs(mean-arr[i]);
	    double pValue = TTest.T_p(diff/stdev, arr.length-1);

	    if(pValue<=0.9)
	    System.out.println( arr[i] ); // + "\t" + pValue );
	}
*/

	//System.out.println( dev.T_p(1.96, 7) );

        double[] arr = {
		0.030224608, 0.012895734, 0.018672025,
        };


        	filterAndRemove(arr, 0.05);

        System.exit(0);
	while(true) {
		double[] tmpArr = filterAndRemove(arr, 0.1);
		if(arr.length == tmpArr.length) break;

		arr = tmpArr;
	}

        for(double d:arr)
	{
System.out.println(d);
}
/*
	double[] arr =
	{
	3.359893581,
	13.1229677,
	1.731213731,
	1.0935822,
	2.691398138,
	0.352295687,
	};
*/

/*

//	for(int i=0;i<arr.length;i++)
//	    newArr[i] = Math.log( arr[i] );

	for(int i=0;i<arr.length;i++)
	    System.out.println(arr[i]);
	double[] result = filterAndRemove(arr, 0.2);
//	double[] result = filterExcludingNegative(arr, 0.2);

	System.out.println("========");
	for(int i=0;i<result.length;i++)
	    System.out.println(result[i]);

	*/
    }

    public static List<MergeProteinModel.Peptide> filterMerge(List<MergeProteinModel.Peptide> list, double pValue) {

		double sum =0;

			for(Iterator<MergeProteinModel.Peptide> itr=list.iterator(); itr.hasNext(); )
		{
				MergeProteinModel.Peptide peptide = itr.next();
				sum+=peptide.getRatio();
			}

		double mean = sum/list.size();

		double devSum=0;
		for(Iterator<MergeProteinModel.Peptide> itr=list.iterator(); itr.hasNext(); )
		{
				MergeProteinModel.Peptide peptide = itr.next();
			devSum += (mean-peptide.getRatio())*(mean-peptide.getRatio());
		}

		double stdev = Math.sqrt(devSum/(list.size()-1));

		for(Iterator<MergeProteinModel.Peptide> itr=list.iterator(); itr.hasNext(); )
		{
				MergeProteinModel.Peptide peptide = itr.next();
			double diff = Math.abs(mean-peptide.getRatio());
			double temp = TTest.T_p(diff/stdev, list.size()-1);

				peptide.setProbability(temp);

			if(temp<=pValue)
					peptide.setFilterOut(true);
		}

        return list;
    }

/*
Computing an approximate P value
You can also calculate an approximate P value as follows.
Calculate  t
N is the number of values in the sample, Z is calculated for the suspected outlier as shown above.
Look up the two-tailed P value for the student t distribution with the calculated value of T and N-2 degrees of freedom. Using Excel, the formula is =TDIST(T,DF,2) (the '2' is for a two-tailed P value).
Multiply the P value you obtain in step 2 by N. The result is an approximate P value for the outlier test. This P value is the chance of observing one point so far from the others if the data were all sampled from a Gaussian distribution. If Z is large, this P value will be very accurate. With smaller values of Z, the calculated P value may be too large.
*/

    public static double getGrubbsPvalue(double d, double mean, double stdev, int n) {
            double z = Math.abs(d-mean)/stdev;

            double t = Math.sqrt( (n*(n-2)*z*z)/((n-1)*(n-1)-n*z*z) );
            TTest ttest = new TTest();
            double tdistvalue = ttest.T_p(t, n-2);
//System.out.println("===" + d + " " + z + " " + t + " " + tdistvalue + " " +tdistvalue*n);


        return tdistvalue*n;
    }

    public static List<ChroPeptide> filter(List<ChroPeptide> list, double pValue) {

	double sum =0;

        for(Iterator<ChroPeptide> itr=list.iterator(); itr.hasNext(); )
	{
            ChroPeptide peptide = itr.next();
            sum+=peptide.getSlope();
	}

	double mean = sum/list.size();

	double devSum=0;
	for(Iterator<ChroPeptide> itr=list.iterator(); itr.hasNext(); )
	{
            ChroPeptide peptide = itr.next();
//	System.out.println("======" + peptide.getSlope() + "\t" + mean);
	    devSum += (mean-peptide.getSlope())*(mean-peptide.getSlope());
	}

	double stdev = Math.sqrt(devSum/(list.size()-1));

	for(Iterator<ChroPeptide> itr=list.iterator(); itr.hasNext(); )
	{
            ChroPeptide peptide = itr.next();
	    double diff = Math.abs(mean-peptide.getSlope());
	    double temp = TTest.T_p(diff/stdev, list.size()-1);

            peptide.setProbability(temp);

	    if(temp<=pValue)
                peptide.setFilterOut(true);
	}

	for(int i=list.size()-1;i>=0;i--)
	{
            ChroPeptide peptide = list.get(i);
	    if(peptide.isFilterOut())
	    {
		list.remove(peptide);
	    }
	}

        return list;
    }

    //calculate p value for array
    //sourcearr is 1st arr
    //targetArr is 2nd arr. calculate p value for targetArr
    public static double[] calculatePValueExcludingNegative(double[] sourceArr, double[] targetArr) {

        //System.out.println(sourceArr.length + " " + targetArr.length);

        double[] allArr = new double[sourceArr.length + targetArr.length];

        for(int i=0;i<sourceArr.length;i++)
            allArr[i] = sourceArr[i];

        for(int i=sourceArr.length;i<allArr.length;i++)
            allArr[i] = targetArr[i-sourceArr.length];


	double stdev = STDev.getStdevWithoutNegative(allArr);
	double mean = STDev.getMeanExcludingNegative(allArr);

	double[] result = new double[targetArr.length];

	for(int i=0;i<result.length;i++)
	{
	    if(targetArr[i]<0)
	    {
                result[i] = -1;
                continue;
	    }

	    double diff = Math.abs(mean-targetArr[i]);
	    double pValue = TTest.T_p(diff/stdev, allArr.length-1);

            result[i] = pValue;
	}

        //for(int i=0;i<result.length;i++)
          //  System.out.println(result[i]);

	return result;

    }

    public static double[] filterExcludingNegative(double[] arr, double pValue) {

	double stdev = STDev.getStdevWithoutNegative(arr);
	double mean = STDev.getMeanExcludingNegative(arr);

	double[] result = new double[arr.length];

	for(int i=0;i<arr.length;i++)
	{
	    if(arr[i]<0)
	    {
		    result[i] = -1;
		    continue;
	    }

	    double diff = Math.abs(mean-arr[i]);

	    double temp = TTest.T_p(diff/stdev, arr.length-1);

	    if(temp<=pValue)
		result[i] = -1;
	    else
		result[i] = arr[i];
            //result[i] = temp;
	}

	return result;

    }

    public static double[] filter(double[] arr, double pValue) {

	double stdev = STDev.getStdev(arr);
	double mean = STDev.getMean(arr);

	double[] result = new double[arr.length];

	for(int i=0;i<arr.length;i++)
	{
	    double diff = Math.abs(mean-arr[i]);
	    double temp = TTest.T_p(diff/stdev, arr.length-1);

	    if(temp<=pValue)
		result[i] = -1;
	    else
		result[i] = arr[i];
	}

	return result;
    }

    public static void filterAndRemove(List<ChroPeptide> list, double pValue) {

        DescriptiveStatistics stat = new DescriptiveStatistics();
        for(Iterator<ChroPeptide> itr=list.iterator(); itr.hasNext(); )
        {
            ChroPeptide peptide = itr.next();
            double logR = Math.log(peptide.getSlope())/Math.log(2);
              stat.addValue(logR);
        }

        double stdev = stat.getStandardDeviation();
        double mean = stat.getMean();
        int n = list.size();
	if(stdev==0) return;

        gnu.trove.TDoubleArrayList result = new gnu.trove.TDoubleArrayList();
        for(Iterator<ChroPeptide> itr=list.iterator(); itr.hasNext(); )
        {
                  ChroPeptide peptide = itr.next();
            double logR = Math.log(peptide.getSlope())/Math.log(2);
            double gp = getGrubbsPvalue(logR, mean, stdev, n);

                  peptide.setProbability(gp);

            if(gp<=pValue)
                      peptide.setFilterOut(true);
        }

	for(int i=list.size()-1;i>=0;i--)
	{
            ChroPeptide peptide = list.get(i);
	    if(peptide.isFilterOut())
	    {
		list.remove(peptide);
	    }
	}

    }
//Added by harshil Shah-- for the peptideLevel Analysis........
     public static int filterAndRemovePeptideNoLog(Peptide peptideData, double pValue) {

         int counter =0;
         Set removingIndex = new HashSet<>();
//         DescriptiveStatistics stat = new DescriptiveStatistics();
         int size1 = peptideData.comparedGroup.size();
         int size2 = peptideData.comparedGroup.get(0).size();
         for(int z=0;z<size2;z++)
         {
              DescriptiveStatistics stat = new DescriptiveStatistics();
             for(int k=0;k<size1;k++)
             {
//                 double logR = Math.log(peptideData.comparedGroup.get(k).get(z).getMean()) / Math.log(2);

                    double val = peptideData.comparedGroup.get(k).get(z).getMean();
                    stat.addValue(val);
             }
             double stdev = stat.getStandardDeviation();
            double mean = stat.getMean();
            int n = size1;
            if (stdev == 0) {
                return 0;
            }
            gnu.trove.TDoubleArrayList result = new gnu.trove.TDoubleArrayList();
             for(int k=0;k<size1;k++)
             {
//                double logR = Math.log(peptideData.comparedGroup.get(k).get(z).getMean()) / Math.log(2);
                    double val = peptideData.comparedGroup.get(k).get(z).getMean();


                double gp = getGrubbsPvalue(val, mean, stdev, n);
                if (gp <= pValue) {
//                    removingIndex.add(z);
                    removingIndex.add(k);
                    counter++;
                }
             }

          }

         peptideData.removeEntry(removingIndex);
         return counter;
    }
// This is for protein level analysis.....
     public static Set filterAndRemovePeptideFromProteinNoLog(List list, double pValue) {

        int counter = 0;
        Set removingIndex = new HashSet<>();
//         DescriptiveStatistics stat = new DescriptiveStatistics();


        DescriptiveStatistics stat = new DescriptiveStatistics();
        for (int k = 0; k < list.size(); k++) {
         //   double logR = Math.log((double) list.get(k)) / Math.log(2);

            stat.addValue((double) list.get(k));
        }
        double stdev = stat.getStandardDeviation();
        double mean = stat.getMean();

//            if (stdev == 0) {
//                return 0;
//            }
        gnu.trove.TDoubleArrayList result = new gnu.trove.TDoubleArrayList();
        for (int k = 0; k < list.size(); k++) {

            double gp = getGrubbsPvalue((double) list.get(k), mean, stdev, list.size());
            if (gp <= pValue) {
                removingIndex.add(k);
                counter++;
            }
        }
        return removingIndex;
    }


    public static double[] filterAndRemove(double[] arr, double pValue) {

	List<Double> l = new ArrayList<Double>();

        DescriptiveStatistics stat = new DescriptiveStatistics();
        for(double d:arr)
	{
            stat.addValue(d);
	}

        double stdev = stat.getStandardDeviation();
        double mean = stat.getMean();

	if(stdev==0 || Double.isNaN(stdev)) {
		return arr;
	}

        int n = arr.length;

        gnu.trove.TDoubleArrayList result = new gnu.trove.TDoubleArrayList();
        for(double d:arr)
	{
	    double gp = getGrubbsPvalue(d, mean, stdev, n);

	    if(gp>pValue) l.add(d);
	}

	double[] narr = convertArray(l);

	//System.out.println("------------" + l);
//        for(double d:narr)
//		System.out.print(d + "\t");

	return narr;

    }

    public static double[] convertArray(List<Double> l)
    {
	    double[] ret = new double[l.size()];
	    for (int i=0; i < ret.length; i++)
	    {
		    ret[i] = l.get(i).doubleValue();
	    }
	    return ret;
    }

    //This is TMT outlier filter under testing for Dmcclat yet.
    public static void filterAndRemoveTMT(List<ChroPeptide> list, double pValue) {

        DescriptiveStatistics stat1 = new DescriptiveStatistics();
        DescriptiveStatistics stat2 = new DescriptiveStatistics();
        DescriptiveStatistics stat3 = new DescriptiveStatistics();
        /*double[] arr1 = new double[list.size()];
        double[] arr2 = new double[list.size()];
        double[] arr3 = new double[list.size()];

        int count=0;*/
        for(Iterator<ChroPeptide> itr=list.iterator(); itr.hasNext(); )
	{
            ChroPeptide peptide = itr.next();
            long[] larr = peptide.getTotalIntArr();

            double[] normarr = new double[larr.length];
            for(int i=0;i<normarr.length;i++) {
                //normarr[i] = larr[i] * averageInt / sumIntArr[i];
                normarr[i] = larr[i];// * averageInt / sumIntArr[i];

            //    System.out.println(normarr[i]);

            }

            double avg1 = (normarr[0]+normarr[3])/2;
            double avg2 = (normarr[1]+normarr[4])/2;
            double avg3 = (normarr[2]+normarr[5])/2;

            double v1 = Math.log(avg1/avg2)/Math.log(2);
            double v2 = Math.log(avg1/avg3)/Math.log(2);
            double v3 = Math.log(avg2/avg3)/Math.log(2);

           // System.out.println("===" + v1);
           // System.out.println("===" + v2);
           // System.out.println("===" + v3);

            stat1.addValue(v1);
            stat2.addValue(v2);
            stat3.addValue(v3);
        }

        double stdev1 = stat1.getStandardDeviation();
        double mean1 = stat1.getMean();
        double stdev2 = stat2.getStandardDeviation();
        double mean2 = stat2.getMean();
        double stdev3 = stat3.getStandardDeviation();
        double mean3 = stat3.getMean();

	if(stdev1==0 || stdev2==0 || stdev3==0) return;

        //gnu.trove.TDoubleArrayList result = new gnu.trove.TDoubleArrayList();
	for(Iterator<ChroPeptide> itr=list.iterator(); itr.hasNext(); )
	{
            ChroPeptide peptide = itr.next();

            long[] larr = peptide.getTotalIntArr();
            double avg1 = (larr[0]+larr[3])/2;
            double avg2 = (larr[1]+larr[4])/2;
            double avg3 = (larr[2]+larr[5])/2;

            double v1 = Math.log(avg1/avg2)/Math.log(2); //126 & 129
            double v2 = Math.log(avg1/avg3)/Math.log(2); //127 & 130
            double v3 = Math.log(avg2/avg3)/Math.log(2); //128 & 131


	    double gp = getGrubbsPvalue(v1, mean1, stdev1, list.size());
           // System.out.println("==" + v1/v2 + " " + mean1 + " " + stdev1 + " " + gp);


            if(gp<=pValue) peptide.setFilterOut(true);
            gp = getGrubbsPvalue(v2, mean2, stdev2, list.size());

           // System.out.println("==" + v1/v3 + " " + mean2 + " " + stdev2 + " " + gp);

            if(gp<=pValue) peptide.setFilterOut(true);
            gp = getGrubbsPvalue(v3, mean3, stdev3, list.size());
            //System.out.println("==" + v2/v3 + " " + mean3 + " " + stdev3 + " " + gp);
            if(gp<=pValue) peptide.setFilterOut(true);


	}

	for(int i=list.size()-1;i>=0;i--)
	{
            ChroPeptide peptide = list.get(i);
	    if(peptide.isFilterOut())
	    {
		list.remove(peptide);
	    }
	}

    }

/*
    public static double[] filterAndRemove(double[] arr, double pValue) {

	double stdev = STDev.getStdev(arr);
	double mean = STDev.getMean(arr);

//	double[] result = new double[arr.length];
	gnu.trove.TDoubleArrayList result = new gnu.trove.TDoubleArrayList();

	for(int i=0;i<arr.length;i++)
	{
	    double diff = Math.abs(mean-arr[i]);
//	    System.out.println( diff/stdev);

	    double temp = TTest.T_p(diff/stdev, arr.length-1);

	    if(temp>pValue)
		result.add(arr[i]);
	}

	return result.toNativeArray();
    }
*/
}
