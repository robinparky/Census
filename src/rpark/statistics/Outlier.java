package rpark.statistics;

import edu.scripps.pms.census.model.MergeProteinModel;
import edu.scripps.pms.stats.TTest;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * Created by rpark on 6/27/16.
 */
public class Outlier {
  public static void main(String[] args) {



    double[] d1 = {0.001,	0.00939,	0.07186,	0.001,	0.001,	0.001,	0.00621,	0.001,	0.05701}; //{-6.64385619,	-6.64385619,	-6.64385619,	-2.395928676};

    double[] d2 ={56.08041,	4.56619,	0.11519,	0.18281,	18.68084,	0.0207,	0.07164,	4.45359}; //{7.639449246,	1.695993813,	5.524502537,	4.859969548,	5.601102203,	3.672425342,	2.432959407,
      //3.876762491,	5.743084056,};
    double[] d3 = {64.41268,	64.6906,	58.19795,	50.67256,	50.55953,	0.04439,	0.16213}; //{8.254697964,	5.277984747,	5.64587455};

    calc15N(d1, 0.35);
    calc15N(d2, 0.35);
    calc15N(d3, 0.35);
  }

  public static void calc15N(double[] arr) {
	calc15N(arr, 0.1);
  }


  public static void calc15N(double[] arr, double pvalue) {

    List<Double> l = Outlier.convertToList(arr);
    //List<Double> l = Outlier.convertToList(d1);
    //List<Double> l = Outlier.convertToList(d1);

    while(true) {

      int count1 = l.size();

      List<Double> newL = Outlier.removeOutlierByLog2(l, pvalue);
      if(newL.size()==count1) break;

      l = newL;
    }

    List<Double> enrichmentList = new ArrayList<>();
    for(double d:l)
     // enrichmentList.add( d );
      enrichmentList.add( 100*(1/(1+d)) );

    double mean = CommonStat.getMeanValue(enrichmentList);
    double[] farr = new double[enrichmentList.size()];
    for(int i=0;i<enrichmentList.size();i++)
      farr[i] = enrichmentList.get(i);

    double standardDev = CommonStat.getStandardDeviation(farr);
    double start = mean - standardDev;
    double end = mean + standardDev;
  //  System.out.println("new\t" + l);
    System.out.println("mean\t" + mean);
    System.out.println("stdev\t" + standardDev);
    System.out.println("new enr\t" + enrichmentList);
    System.out.println("error bar\t" + start + "\t" + end);

  }


  public static List<Double> convertToList(double[] arr) {
    List<Double> l = new ArrayList<>();
    for(double d:arr)
      l.add(d);

    return l;
  }


    public static DescriptiveStatistics removeOutlier(double[] arr, double pValue) {
        double sum =0;

        for(double d:arr)
        {

            sum+=d;
        }

        double mean = sum/arr.length;

        double devSum=0;
        for(double d:arr) {
            double diff = mean-d;

            devSum += diff*diff;
        }

        double stdev = Math.sqrt(devSum/(arr.length-1));

        for(double d:arr)
        {
            double diff = Math.abs(mean-d);
            double prob = TTest.T_p(diff / stdev, arr.length - 1);

            System.out.println(d + "\t" + prob);
//            if(temp<=pValue)
  //              peptide.setFilterOut(true);
        }

    //    return list;
        return null;

    }

    public static List<Double> removeOutlierByLog2(List<Double> list, double pValue) {
      List<Double> l = new ArrayList<>();
      for(Double d:list) {
        l.add(Math.log(d)/Math.log(2));
      }

      List<Double> resultList = new ArrayList<>();
      for(Double d:removeOutlier(l, pValue))
        resultList.add( Math.exp(d * Math.log(2)) );

      return resultList;
    }

    public static List<Double> removeOutlier(List<Double> list, double pValue) {
        double sum =0;
        List<Double> result = new ArrayList<>();
        for(double d:list)
        {

            sum+=d;
        }

        double mean = sum/list.size();

        double devSum=0;
        for(double d:list) {
            double diff = mean-d;

            devSum += diff*diff;
        }

        double stdev = Math.sqrt(devSum/(list.size()-1));


        for(double d:list)
        {
            if(!Double.isNaN(d) && !Double.isInfinite(d)) {
                double diff = Math.abs(mean - d);
                double prob = TTest.T_p(diff / stdev, list.size() - 1);

           //     System.out.println(d + "\t" + prob);

                if (Double.compare(prob, pValue) > 0 ) {
                    result.add(d);
                }
            }

        }

        return result;

    }
    public static List<Double> removeOutlierZero(List<Double> list, double pValue) {
        boolean  removedOutlier = true;
        List<Double> result = new ArrayList<>();
        while(removedOutlier)
        {
            double sum =0;
            for(double d:list)
            {

                sum+=d;
            }

            double mean = sum/list.size();

            double devSum=0;
            for(double d:list) {
                double diff = mean-d;

                devSum += diff*diff;
            }

            double stdev = Math.sqrt(devSum/(list.size()-1));


            for(double d:list)
            {
                if(!Double.isNaN(d) && !Double.isInfinite(d) &&(Double.compare(d,0.0)!=0)) {
                    double diff = Math.abs(mean - d);
                    double prob = TTest.T_p(diff / stdev, list.size() - 1);

                    //System.out.println(d + "\t" + prob);

                    if (Double.compare(prob, pValue) > 0 ) {
                        result.add(d);
                    }
                }

            }
            if(list.size() == result.size()) removedOutlier = false;
            else
            {
                list.clear();
                list.addAll(result);
                result.clear();
            }
        }


        return result;

    }
    public static List<double[]> removeOutlierDoubleArr(List<double[]> list, double pValue) {
        double sum =0;
        List<double []> result = new ArrayList<>();
        for(double[] darr:list)
        {

            for(double d:darr)
            {
                if(Double.isInfinite(d) || Double.isNaN(d)) continue;

                sum+=d;
            }

        }

        double mean = sum/list.size();

        double devSum=0;
        for(double[] darr:list) {
            for(double d:darr)
            {
                if(Double.isInfinite(d) || Double.isNaN(d)) continue;

                double diff = mean-d;
                devSum += diff*diff;
            }

        }

        double stdev = Math.sqrt(devSum/(list.size()-1));


        for(double[]  darr:list)
        {
            List<Double> temp = new ArrayList<>();
            for(double d:darr)
            {
                if(Double.isInfinite(d) || Double.isNaN(d)) continue;
                double diff = Math.abs(mean - d);
                double prob = TTest.T_p(diff / stdev, list.size() - 1);

                //System.out.println(d + "\t" + prob);
                if (Double.compare(prob, pValue) >= 0) {
                    temp.add(d);
                }
            }
            double[] add = new double[temp.size()];
            for(int i=0; i<temp.size();i++)
            {
                add[i] = temp.get(i);
            }
            result.add(add);
        }

        return result;

    }

    public static DescriptiveStatistics removeOutlier(DescriptiveStatistics list, double pValue) {

        return removeOutlier(list.getValues(), pValue);


    }

}
