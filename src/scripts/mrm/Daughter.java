
/*
* Copyright (c) 2008 Integrated Proteomics Applications.  All rights reserved.  
*/

package scripts.mrm;

/**
 *
 * @author Sung Kyu, Robin, Park
 * @email robinparky@yahoo.com
 * Created on Feb 1, 2010 
 * $Revision:$
 * $Date:$
 */

import java.util.*;

public class Daughter {

    private double mass;
    private List<Double> intensityList = new ArrayList<Double>();
    private List<Integer> scanList = new ArrayList<Integer>();
    private List<Double> rtList = new ArrayList<Double>();
    private int startIndex;
    private int endIndex;
    
    public Daughter(double mass) {

        this.mass = mass;
    }

    public void addIntensity(double d, int i, double rt) {
        this.intensityList.add(d);
        this.scanList.add(i);
        this.rtList.add(rt);
    }

    public void findPeakArea(int peakIndex) {
        startIndex=peakIndex-5;
        while(true) {
            double area1=0;
            double area2=0;

            for(int i=0;i<3;i++) {
                area1 += intensityList.get(startIndex-i);
            }

            for(int i=1;i<4;i++) {
                
                area2 += intensityList.get(startIndex-i);
            }


            if(area2<=0 || area2/area1>=0.95)
                break;

            startIndex--;
        }

        

        endIndex=peakIndex+5;
        while(true) {
            double area1=0;
            double area2=0;

            for(int i=0;i<3;i++) {
                area1 += intensityList.get(endIndex+i);
            }

            for(int i=1;i<4;i++) {
                area2 += intensityList.get(endIndex+i);
            }

//System.out.println( (startIndex) + " " + area1 + " " + area2 + " " + (area1/area2));

            if(area2<=0 || area2/area1>=0.95)
                break;

            endIndex++;
        }


    }
   
    public double getIntensity(int index) {
	return this.intensityList.get(index);
    }

    public List<Double> getIntensityList() {
        return intensityList;
    }

    public void setIntensityList(List<Double> intensityList) {
        this.intensityList = intensityList;
    }

    public int getScan(int index) {
	return this.scanList.get(index);	
    }

    public List<Integer> getScanList() {
        return scanList;
    }

    public void setScanList(List<Integer> scanList) {
        this.scanList = scanList;
    }


    public double getMass() {
        return mass;
    }

    public void setMass(double mass) {
        this.mass = mass;
    }

    public double getRt(int index) {
	return this.rtList.get(index);	
    }
    
    public List<Double> getRtList() {
        return rtList;
    }

    public void setRtList(List<Double> rtList) {
        this.rtList = rtList;
    }

    public int getEndIndex() {
        return endIndex;
    }

    public void setEndIndex(int endIndex) {
        this.endIndex = endIndex;
    }

    public int getStartIndex() {
        return startIndex;
    }

    public void setStartIndex(int startIndex) {
        this.startIndex = startIndex;
    }

    public int getStartScan() {
        return this.getScan(this.startIndex);
    }

    public int getEndScan() {
        return this.getScan(this.endIndex);
    }

}
