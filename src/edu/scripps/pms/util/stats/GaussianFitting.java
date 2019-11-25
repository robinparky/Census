/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.scripps.pms.util.stats;

/**
 *
 * @author rpark
 */

import edu.scripps.pms.census.labelFree.PeakRange;
import edu.scripps.pms.census.model.LabelingResultModel;
import edu.scripps.pms.util.spectrum.Range;
import java.util.Arrays;
import org.apache.commons.math3.fitting.*;
import org.apache.commons.math3.fitting.WeightedObservedPoints;
//import flanagan.analysis.CurveSmooth;
import java.util.ArrayList;
import java.util.List;
import java.util.Hashtable;

public class GaussianFitting {
    /*
    public static double[] smooth(double[] yArr, int width) {
        CurveSmooth csm = new CurveSmooth(yArr);
        
        return csm.movingAverage(width);           
    }
    
    public static void smooth(double[] xArr, double[] yArr, int width) {
        CurveSmooth csm = new CurveSmooth(xArr, yArr);
        csm.movingAverage(width);        
    }
    */
    
    public static void main(String[] args) {
        
        
        
        double[] yArr = { 
            531026, 
 984167, 
 1887233, 
 2687152, 
 3461228, 
 3580526, 
 3439750, 
 2877648, 
 2175960, 
 1447024, 
 717104, 
 620014, 
 620015, 
 620016, 
 620017, 
 620018, 
 620019, 
3580526, 
3439750, 
2877648, 
2175960, 
        };
        
        double[] xArr = new double[yArr.length];
        for(int i=1;i<=xArr.length;i++) {
            xArr[i-1] = i;
        }
        
        /*
        double[] dArr = GaussianFitting.smooth(yArr, 9);
        
        for(double d: dArr)
            System.out.println(d);

        Range range = GaussianFitting.getGaussianPeakRange(xArr, yArr);
        
        System.out.println(range.getLowBound() + "\t" + range.getHighBound());
        
        if(true) return;
        
        WeightedObservedPoints obs = new WeightedObservedPoints();
        
        obs.add(1, 531026);
        obs.add(2, 984167);
        obs.add(3, 1887233);
        obs.add(4, 2687152);
        obs.add(5, 3461228);
        obs.add(6, 3580526);
        obs.add(7, 3439750);
        obs.add(8, 2877648);
        obs.add(9, 2175960);
        obs.add(10, 1447024);
        obs.add(11, 717104);
        obs.add(12, 620014);
        obs.add(13, 620015);
        obs.add(14, 620016);
        obs.add(15, 620017);
        obs.add(16, 620018);
        obs.add(17, 620019);
        obs.add(18,3580526);
        obs.add(19,3439750);
        obs.add(20,2877648);
        obs.add(21,2175960);
        */
        
/*
        obs.add(4.0254623,  531026.0);
        obs.add(4.03128248, 984167.0);
        obs.add(4.03839603, 1887233.0);
        obs.add(4.04421621, 2687152.0);
        obs.add(4.05132976, 3461228.0);
        obs.add(4.05326982, 3580526.0);
        obs.add(4.05779662, 3439750.0);
        obs.add(4.0636168,  2877648.0);
        obs.add(4.06943698, 2175960.0);
        obs.add(4.07525716, 1447024.0);
        obs.add(4.08237071, 717104.0);
        obs.add(4.08366408, 620014.0);
        obs.add(5.08366408,620015);
        obs.add(6.08366408,620016);
        obs.add(7.08366408,620017);
        obs.add(8.08366408,620018);
        obs.add(9.08366408,620019);
        obs.add(10.08366408,3580526);
        obs.add(11.08366408,3439750);
        obs.add(12.08366408,2877648);
        obs.add(13.08366408,2175960);
        */
        //double[] parameters = GaussianCurveFitter.create().fit(obs.toList());
        //System.out.println("done");
        
    }
    
    //return parameters: y, x, and sigma
    //peak ranges from -3 x sigma to 3 x sigma
    public static double[] getGaussianCurveFit(double[] xArr, double[] yArr) {
        
        WeightedObservedPoints obs = new WeightedObservedPoints();
        
        for(int i=0;i<xArr.length;i++) {
            obs.add(xArr[i], yArr[i]);
            
            //System.out.println(xArr[i] + "\t" + yArr[i]);
        }
      //  System.out.println("x-" + Arrays.toString(xArr));
      //  System.out.println("y-"+Arrays.toString(yArr));

      //  long startTime = System.currentTimeMillis(); //fetch starting time
        double[] parameters = null;
        
        //GaussianCurveFitter cf = GaussianCurveFitter.create().withMaxIterations(1000);
          
        try {
            parameters = GaussianCurveFitter.create().withMaxIterations(1000).fit(obs.toList());
        } catch(Exception ex) {
            
        }

        return parameters;        
    }
    
    
    //return parameters: y, x, and sigma
    //peak ranges from -3 x sigma to 3 x sigma
    public static Range getGaussianPeakRange(double[] xArr, double[] yArr) {
        
        double[] params = getGaussianCurveFit(xArr, yArr);
        
        double start = -2 * params[2] + params[1];
        double end = 2 * params[2] + params[1];
        
        int startIndex ;
        int endIndex;
        int i =0;
       
        while( i < xArr.length && xArr[i] < start )
            i++;
        if(i>1)
            start = xArr[i-1];
        else
            start = xArr[0];
        
        while(i<xArr.length && xArr[i] < end )
            i++;
        if(i<xArr.length-1)
            end = xArr[i+1];
        else
            end = xArr[xArr.length-1];
        
        return new Range(start, end);
    }
    
    //return parameters: y, x, and sigma
    //peak ranges from -3 x sigma to 3 x sigma
    public static Range getGaussianPeakRangeIndex(double[][] resultArr, int startIndex, int endIndex, int includeIndex) {
        
       // List<Double> xList = new ArrayList<Double>();
       // List<Double> yList = new ArrayList<Double>();
          
        double[] xArr = new double[endIndex-startIndex+1];
        double[] yArr = new double[xArr.length];

        int count=0;
        for(int i=startIndex;i<=endIndex;i++) {
            xArr[count] = resultArr[i][0]; 
            yArr[count] = resultArr[i][1] + resultArr[i][2];
            
           // System.out.println(xArr[count] + "\t" + yArr[count]);
            
            count++;
        }
        
        double[] params = getGaussianCurveFit(xArr, yArr);
        
        if(params==null) return new Range(startIndex, endIndex);
        
        double start = -2 * params[2] + params[1];
        double end = 2 * params[2] + params[1];
        
        int i =0;
       
        while( i < xArr.length && xArr[i] < start )
            i++;
        if(i>1)
            start = xArr[i-1];
        else
            start = xArr[0];
        
        while(i<xArr.length && xArr[i] < end )
            i++;
        if(i<xArr.length-1)
            end = xArr[i+1];
        else
            end = xArr[xArr.length-1];
        
        for(int j=startIndex;j<=endIndex;j++) {
            if(resultArr[j][0]>=start) {
                if((j-1)>=0) start = resultArr[j-1][0];
                else start = resultArr[j][0];
                
                break;
            }
        }
        
        for(int j=startIndex;j<=endIndex;j++) {
            if(resultArr[j][0]>=end && (j+1)>=0) {
                if( (j+1)<resultArr.length ) end = resultArr[j+1][0];
                else end = resultArr[j][0];
                
                break;
            }            
        }
        
        //include identified peak
        if(start>includeIndex) start = includeIndex;
        if(end<includeIndex) end = includeIndex;
        
        return new Range(start, end); 
    }

    
    //return parameters: y, x, and sigma
    //peak ranges from -3 x sigma to 3 x sigma
    public static Range getGaussianPeakRangeIndex(LabelingResultModel[] resultArr, int startIndex, int endIndex, int includeScan) {
        
        List<Double> xList = new ArrayList<Double>();
        List<Double> yList = new ArrayList<Double>();
                
        Hashtable<Double, Integer> scanIndexHt = new Hashtable<Double, Integer>();
        for(int i=startIndex;i<=endIndex;i++) {
            
        //    System.out.println("====" + (double)resultArr[i].getScanNum());
            if(null == resultArr[i]) continue;
            
            scanIndexHt.put((double)resultArr[i].getScanNum(), i);
            xList.add((double)resultArr[i].getScanNum());
            yList.add(resultArr[i].getTotalIntensity());            
        }
        
        double[] xArr = new double[xList.size()];
        double[] yArr = new double[xArr.length];
        
        for(int i=0;i<xList.size();i++) {
            xArr[i] = xList.get(i);
            yArr[i] = yList.get(i);
        }
        
        double[] params = getGaussianCurveFit(xArr, yArr);
        /////////////////HOLD ON THIS................
        if(params==null) return new Range(startIndex, endIndex);
        
        double start = -3 * params[2] + params[1];
        double end = 3 * params[2] + params[1];
//        System.out.println("start...........");
        
        int i =0;
       
        while( i < xArr.length && xArr[i] < start )
            i++;
        if(i>1)
            start = xArr[i-1];
        else
            start = xArr[0];
        
        while(i<xArr.length && xArr[i] < end )
            i++;
        if(i<xArr.length-1)
            end = xArr[i+1];
        else
            end = xArr[xArr.length-1];


        //if( resultArr[(int)start].getScanNum() >includeScan) start = includeScan;
        if( start>includeScan && start>resultArr[startIndex].getScanNum() ) {

            //start = includeScan;
            start = (includeScan>resultArr[startIndex].getScanNum())?includeScan:resultArr[startIndex].getScanNum();
        }
        //if(resultArr[(int)end].getScanNum()<includeScan) end = includeScan;
        if(end<includeScan && end<resultArr[endIndex].getScanNum()) {

            end = (includeScan<resultArr[endIndex].getScanNum())?includeScan:resultArr[endIndex].getScanNum();
        }

        //System.out.println("==========" + includeScan);
        //System.out.println("==========" + scanIndexHt);
//System.out.println("==========" + start + " " + end + " " + scanIndexHt.get(start) + "\t" +  scanIndexHt.get(end));
        return new Range(scanIndexHt.get(start), scanIndexHt.get(end));

    }
    
    
    
}
