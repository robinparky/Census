/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rpark.statistics;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 *
 * @author rpark
 */
public class CommonStat {


  public static double getStandardDeviation(List<Double> list){
    DescriptiveStatistics stat = new DescriptiveStatistics();
    for(Double d:list)
      stat.addValue(d);

    return stat.getStandardDeviation();
  }

    public static double getStandardDeviation(double[] arr){
        DescriptiveStatistics stat = new DescriptiveStatistics();
        for(double d:arr)
            stat.addValue(d);

        return stat.getStandardDeviation();
    }

    public static double getMedianValue(double[] arr) {
        DescriptiveStatistics stats = new DescriptiveStatistics();

        for(double d:arr)
            stats.addValue(d);

        double median = stats.getPercentile(50);
        return median;
    }


    public static double getMedianValue(List<Double> list) {
        DescriptiveStatistics stats = new DescriptiveStatistics();

        for(Double d : list)
            stats.addValue(d);

        double median = stats.getPercentile(50);
        return median;
    }

    public static double getMeanValue(List<Double> list) {
        DescriptiveStatistics stats = new DescriptiveStatistics();

        for(Double d : list)
            stats.addValue(d);

        return stats.getMean();
    }

    public static double getMeanValue(double[] arr) {
        DescriptiveStatistics stats = new DescriptiveStatistics();

        for(double d:arr)
            stats.addValue(d);

	    return stats.getMean();
    }

    /*
    convert ratio to log
    calculate stdev
    delog
 */
    public static double getRatioStdevValueWithMax(List<Double> list, double min,double max){
        org.apache.commons.math3.stat.descriptive.DescriptiveStatistics stats = new org.apache.commons.math3.stat.descriptive.DescriptiveStatistics();

        for(Double d : list) {
            if(d<min){
                d = min;
            }
            if(d>max){
                d=max;
            }

            stats.addValue(Math.log(d) / Math.log(2));
        }

        double stdevLog = stats.getStandardDeviation();

        return Math.pow(2, stdevLog);
    }


    public static double getRatioStdevValue(List<Double> list, double nonZeroValue) {
        org.apache.commons.math3.stat.descriptive.DescriptiveStatistics stats = new org.apache.commons.math3.stat.descriptive.DescriptiveStatistics();

        for(Double d : list) {
            if(d<=0) d = nonZeroValue;

            stats.addValue(Math.log(d) / Math.log(2));
        }

        double stdevLog = stats.getStandardDeviation();

        return Math.pow(2, stdevLog);
    }

    public static double getRatioMeanValue(List<Double> list, double nonZeroValue) {
        org.apache.commons.math3.stat.descriptive.DescriptiveStatistics stats = new org.apache.commons.math3.stat.descriptive.DescriptiveStatistics();

        for(Double d : list) {
            if(d<=0) d = nonZeroValue;

            stats.addValue( Math.log(d)/Math.log(2) );
        }

        return Math.pow(2, stats.getMean());
    }

    public static double getRatioMedianValue(List<Double> list) {
        DescriptiveStatistics stats = new DescriptiveStatistics();

        for(Double d : list) {
        //    if(d<=0) d = nonZeroValue;

            stats.addValue( Math.log(d)/Math.log(2) );
        }

        double median = stats.getPercentile(50);

        return Math.pow(2, median);

    }


    //relative standard deviation
    public static double getRelativeStandardDeviation(List<Double> list) {

      double mean = getMeanValue(list);
      double stdeve = getStandardDeviation(list);

      return stdeve/mean*100;
    }

    //coefficient of variance
    public static double getCoefficientOfVariance(List<Double> list) {

      double mean = getMeanValue(list);
      double stdeve = getStandardDeviation(list);

      return Math.abs(stdeve/mean*100);
      //   stdev/mean
    }



  public static void main(String[] args) {
        List<Double> l = new ArrayList<>();
        l.add(300.0);
        l.add(20.0);
        l.add(50.0);
        l.add(5.0);

       // CommonStat.getRatioStdevValue(l, 0.001);
       double value = CommonStat.getRatioStdevValueWithMax(l, 0.033, 30);
        System.out.println("");
    }

}
