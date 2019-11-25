/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.scripps.pms.stats;

/**
 *
 * @author rpark
 */
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

public class StatUtil {
    
    public static double getMedian(double[] arr) {
        DescriptiveStatistics statIntensityValues = new DescriptiveStatistics(arr);         
        double median = statIntensityValues.getPercentile(50);
        return median;        
    }
    
}
