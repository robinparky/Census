/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rpark.statistics;

/**
 *
 * @author rpark
 * @author rohan
 */

//import edu.scripps.pms.census.labelFree.PeakRange;
//import edu.scripps.pms.census.model.LabelingResultModel;

import edu.scripps.pms.census.util.dtaselect.Peptide;
import java.util.Arrays;
import org.apache.commons.math3.fitting.*;
import org.apache.commons.math3.fitting.WeightedObservedPoints;
//import flanagan.analysis.CurveSmooth;
import java.util.ArrayList;
import java.util.List;
import java.util.Hashtable;
import org.apache.commons.math3.analysis.function.Gaussian;
import rpark.statistics.model.GaussianPeakModel;

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
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                773376,
                1836191,
                738394,
                1727735,
                1243792,
                1537636,
                1131183,
                1658066,
                1588888,
                961123,
                1551891,
                1419243,
                2309482,
                1292776,
                1212351,
                1238360,
                1105376,
                1046221,
                3204661,
                330720,
                1520253,
                754539,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,


        };

        double[] xArr = new double[yArr.length];
        for(int i=1;i<xArr.length;i++) {
            xArr[i-1] = i;

      //      System.out.println(i + "\t" + yArr[i]);
        }

        GaussianPeakModel gm = GaussianFitting.getGaussianPeakRangeIndex(xArr, yArr, 0, yArr.length-1);

        int start = gm.getPeakStartIndex();
        int end = gm.getPeakEndIndex();


      //  System.out.println("-=-----" + gm + " " + start + " " + end + " " + yArr[start] + " " + yArr[end]);
      //System.out.println(gm.getPeakArea() + " " + gm.getSigma());


/*
      int[] scanArr = peakModel.getScanArr();
      double[] retArr = peakModel.getRetArr();
      double[] peakArr = peakModel.getPeakArr();

      double[] gxArr = peakModel.getGaussianXArr();
      double[] gyArr = peakModel.getGaussianYArr();

      expPep.setChroData(sb.toString());
      expPep.setGaussianPeakString(gPeakSb.toString());
      expPep.setPeakSigma(peakModel.getSigma());
      expPep.setPeakx(peakModel.getX());
      expPep.setPeaky(peakModel.getY());
      expPep.setIsoArr(peakModel.getIsoArr());



        /*
        double[] dArr = GaussianFitting..smooth(yArr, 9);

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
    public static GaussianPeakModel getGaussianPeakRange(double[] xArr, double[] yArr) {

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

        return new GaussianPeakModel(start, end);
    }

    //return parameters: y, x, and sigma
    //peak ranges from -3 x sigma to 3 x sigma
    public static GaussianPeakModel getGaussianPeakRangeIndex(double[][] resultArr, int startIndex, int endIndex, int includeIndex) {

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

        if(params==null) return new GaussianPeakModel(startIndex, endIndex);

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

        return new GaussianPeakModel(start, end);
    }

        //return parameters: y, x, and sigma
    //peak ranges from -3 x sigma to 3 x sigma
    public static double[] getGaussianCurveFitRange(double[] xArr, double[] yArr, int startIndex, int endIndex) {

        WeightedObservedPoints obs = new WeightedObservedPoints();

        for(int i=startIndex+1;i<=endIndex;i++) {
          if(xArr[i]<=0) continue;

        //  System.out.println(startIndex + " " + endIndex + " " + i + "\t==\t" + xArr[i] + "\t" +  yArr[i]);

          if(yArr[i]<0) yArr[i]=0;
            obs.add(xArr[i], yArr[i]);

            //System.out.println(xArr[i] + "\t" + yArr[i]);
        }

        double[] parameters = null;

        try {
            parameters = GaussianCurveFitter.create().withMaxIterations(1000).fit(obs.toList());
        } catch(Exception ex) {

        }
        return parameters;
    }

    public static GaussianPeakModel getGaussianPeakRangeIndex(long[] xArr, long[] yArr, int startIndex, int endIndex) {

      double[] dxArr = new double[xArr.length];
      for(int i=0;i<xArr.length;i++) {
        dxArr[i] = xArr[i];
      }

      double[] dyArr = new double[yArr.length];
      for(int i=0;i<yArr.length;i++) {
        dyArr[i] = yArr[i];
      }

      return getGaussianPeakRangeIndex(dxArr, dyArr, 0, dyArr.length-1);

    }

    //return parameters: y, x, and sigma
    //peak ranges from -3 x sigma to 3 x sigma
    public static GaussianPeakModel getGaussianPeakRangeIndex(double[] xArr, double[] yArr, int startIndex, int endIndex) {

        double[] params = getGaussianCurveFitRange(xArr, yArr, startIndex, endIndex);


        /////////////////HOLD ON THIS................
        if(params==null) {
            GaussianPeakModel model = new GaussianPeakModel();
            model.setPeakStartIndex(startIndex);;
            model.setPeakEndIndex(endIndex);
            return model;
        }

        double start = -3 * params[2] + params[1];
        double end = 3 * params[2] + params[1];


        int peakStartIndex=-1;
        int peakEndIndex=-1;


//        System.out.println("start...........");

        int i =0;

        while( i < xArr.length && xArr[i] < start )
            i++;
        if(i>1) {
            start = xArr[i-1];
            peakStartIndex=i-1;
        } else {
            start = xArr[0];
            peakStartIndex=0;
        }

        while(i<xArr.length && xArr[i] < end )
            i++;
        if(i<xArr.length-1) {
            end = xArr[i+1];
            peakEndIndex=i+1;
        } else {
            end = xArr[xArr.length-1];
            peakEndIndex=xArr.length-1;
        }

        GaussianPeakModel peakModel = new GaussianPeakModel(start, end);
        peakModel.setX(params[1]);
        peakModel.setY(params[0]);
        peakModel.setSigma(params[2]);
        peakModel.setPeakStartIndex(peakStartIndex);
        peakModel.setPeakEndIndex(peakEndIndex);


        //draw gaussian peaks
        Gaussian g = new Gaussian(peakModel.getY(), peakModel.getX(), peakModel.getSigma());
        double gstart = -4 * params[2] + params[1];
        double gend = 4 * params[2] + params[1];

        List<Double> gxList = new ArrayList<>();
        List<Double> gyList = new ArrayList<>();

        for(double d:xArr) {
            if(d<gstart || d>gend)
                continue;

            gxList.add(d);
            gyList.add(g.value(d));

        }

        double[] gxArr = new double[gxList.size()];
        double[] gyArr = new double[gxList.size()];
        for(int j=0;j<gxList.size();j++) {
            gxArr[j] =gxList.get(j).doubleValue();
            gyArr[j] =gyList.get(j).doubleValue();

        }

        peakModel.setGaussianXArr(gxArr);
        peakModel.setGaussianYArr(gyArr);

        return peakModel;

    }
    public static GaussianPeakModel getGaussianPeakRangeIndex2(Peptide peptide,double[] xArr, double[] yArr,
                                                               int startIndex, int endIndex,double[] xIonArr,
                                                               double[] yIonArr, int startIonIndex, int endIonIndex) {

      //  if(peptide.getScanNum().equals("140847") && peptide.getFileName().equals("20141101_HeLa_1ug_BEH60_140min_35ms_IC_DE5_5e3_2")){
      //      System.out.println("");
      // }


        double[] params = getGaussianCurveFitRange(xArr, yArr, startIndex, endIndex);
        double [] paramsIon= getGaussianCurveFitRange(xArr, yIonArr, startIonIndex, endIonIndex);


        /////////////////HOLD ON THIS................
        if(params==null) {

            GaussianPeakModel model = new GaussianPeakModel();
            model.setPeakStartIndex(startIndex);;
            model.setPeakEndIndex(endIndex);
            return model;
        }

        double start = -3 * params[2] + params[1];
        double end = 3 * params[2] + params[1];


        int peakStartIndex=-1;
        int peakEndIndex=-1;


//        System.out.println("start...........");

        int i =0;

        while( i < xArr.length && xArr[i] < start )
            i++;
        if(i>1) {
            start = xArr[i-1];
            peakStartIndex=i-1;
        } else {
            start = xArr[0];
            peakStartIndex=0;
        }

        while(i<xArr.length && xArr[i] < end )
            i++;
        if(i<xArr.length-1) {
            end = xArr[i+1];
            peakEndIndex=i+1;
        } else {
            end = xArr[xArr.length-1];
            peakEndIndex=xArr.length-1;
        }

        GaussianPeakModel peakModel = new GaussianPeakModel(start, end);
        peakModel.setX(params[1]);
        peakModel.setY(params[0]);
        peakModel.setSigma(params[2]);

      /*
System.out.println("===" + params[1]);
System.out.println("===" + params[0]);
System.out.println("===" + params[2]);
System.out.println("===" + start);
System.out.println("===" + end);
System.out.println("peak area.===========" + peakModel.getPeakArea());; */


        peakModel.setPeakStartIndex(peakStartIndex);
        peakModel.setPeakEndIndex(peakEndIndex);
        if (paramsIon == null){
          //  for(int j =0;j<yArr.length;j++){
          //  System.out.println(""+yArr[j]);
          //  }

            peakModel.setxIonInjectionCorrection(-1);
            peakModel.setyIonInjectionCorrection(-1);
            peakModel.setSigmaIonInjectCorrection(-1);

        }else{
            peakModel.setxIonInjectionCorrection(paramsIon[1]);
        peakModel.setyIonInjectionCorrection(paramsIon[0]);
        peakModel.setSigmaIonInjectCorrection(paramsIon[2]);
        }


        //draw gaussian peaks
        Gaussian g = new Gaussian(peakModel.getY(), peakModel.getX(), peakModel.getSigma());
        double gstart = -4 * params[2] + params[1];
        double gend = 4 * params[2] + params[1];

        List<Double> gxList = new ArrayList<>();
        List<Double> gyList = new ArrayList<>();

        for(double d:xArr) {
            if(d<gstart || d>gend)
                continue;

            gxList.add(d);
            gyList.add(g.value(d));

        }

        double[] gxArr = new double[gxList.size()];
        double[] gyArr = new double[gxList.size()];
        for(int j=0;j<gxList.size();j++) {
            gxArr[j] =gxList.get(j).doubleValue();
            gyArr[j] =gyList.get(j).doubleValue();

        }

        peakModel.setGaussianXArr(gxArr);
        peakModel.setGaussianYArr(gyArr);


    //  System.out.println("remove this..........................");
    //  peakModel.getPeakArea();

        return peakModel;

    }


}
