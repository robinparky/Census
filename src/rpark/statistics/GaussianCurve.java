/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rpark.statistics;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.fitting.GaussianCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoints;
import org.apache.commons.math3.special.Erf;
import org.apache.commons.math3.analysis.function.Gaussian;

/**
 *
 * @author rpark
 */
public class GaussianCurve {
    
    //Math.sqrt(2)
    public static final double SQRT_2= 1.41421356237;
    public static void main(String[] args) {
        
        double sigma = 0.22925237900504442;
        double x = 38.488688946275126;
        double y = 594399.7092468902;
        
        /*
        NormalDistribution norm = new NormalDistribution(x, sigma);
        */
        
        Gaussian g = new Gaussian(y, x, sigma);
        double start = x-4*sigma;
        double end = x+4*sigma;

        for(double i=start;i<=end;i+=0.1) {
            //System.out.println( i + "\t" + GaussianCurve.cdf( (double)i, x, sigma ) );
            System.out.println( i + "\t" + g.value(i));
            //System.out.println( i + "\t" + norm.cumulativeProbability(i) );
        }

        
    }
    
    //draw gaussian graph
    public static double cdf(double x, double u, double sigma) {
        double cfd = Erf.erf( (x-u)/(sigma*SQRT_2) ) +1;
        cfd = cfd/2;
        //double cfd = Erf.
        
        return cfd;
    }
    
}
