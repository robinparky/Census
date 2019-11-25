/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package edu.scripps.pms.census.tools;

/**
 *
 * @author Harshil
 */
public class IsotopeModel {
    
    double[] intensityArr;
    double retentionStartTime;
    double retentionEndTime;
    double[] isoArr; //mass list
    private double retentionTime;
    double ionInjectionTime;
    private int scanNumber=0;
   
    public IsotopeModel( double retentionStartTime, double retentionEndTime, double[] isoArr,double actualRetentionTime, double ionInjectionTime,double[] intensityArr,int scanNumber) {
        this.intensityArr = intensityArr;
        this.retentionStartTime = retentionStartTime;
        this.retentionEndTime = retentionEndTime;
        this.isoArr = isoArr;
        this.retentionTime = actualRetentionTime;
        this.ionInjectionTime = ionInjectionTime;
        this.scanNumber = scanNumber;
    }

    /**
     * 
     * @param isoArr
     * @param actualRetentionTime
     * @param ionInjectionTime
     * @param intensityArr
     * @param scanNumber 
     */
    public IsotopeModel(double[] isoArr,double actualRetentionTime, double ionInjectionTime,double[] intensityArr,int scanNumber) {
        this.intensityArr = intensityArr;
        this.retentionStartTime = retentionStartTime;
        this.retentionEndTime = retentionEndTime;
        this.isoArr = isoArr;
        this.retentionTime = actualRetentionTime;
        this.ionInjectionTime = ionInjectionTime;
        this.scanNumber = scanNumber;
    }
    
    public IsotopeModel() {
    }
    
    
    
    public void setIntensityArr(double[] intensityArr) {
        this.intensityArr = intensityArr;
    }

    public void setRetentionStartTime(double retentionStartTime) {
        this.retentionStartTime = retentionStartTime;
    }

    public void setRetentionEndTime(double retentionEndTime) {
        this.retentionEndTime = retentionEndTime;
    }

    public void setIsoArr(double[] isoArr) {
        this.isoArr = isoArr;
    }

    public void setIonInjectionTime(double ionInjectionTime) {
        this.ionInjectionTime = ionInjectionTime;
    }

    public double[] getIntensityArr() {
        return intensityArr;
    }

    public double getRetentionStartTime() {
        return retentionStartTime;
    }

    public double getRetentionEndTime() {
        return retentionEndTime;
    }

    public double[] getIsoArr() {
        return isoArr;
    }

    public double getIonInjectionTime() {
        return ionInjectionTime;
    }

    public int getScanNumber() {
        return scanNumber;
    }
    
            
    public String getIsoString()
    {
        StringBuffer sb = new StringBuffer();
        sb.append(scanNumber).append(":").append(retentionTime).append(":").append(ionInjectionTime).append(" ");
        for(double value : intensityArr)
        {
            sb.append(value+" ");
        }
        
        
        return sb.toString();
    }
    
    public void setScanNumber(int scanNumber) {
        this.scanNumber = scanNumber;
    }

    public double getRetentionTime() {
        return retentionTime;
    }

    public void setRetentionTime(double retentionTime) {
        this.retentionTime = retentionTime;
    }
    public double getIntensitySum()
    {
        double total =0;
        for(double x : intensityArr)
            total+=x;
        return total;
    }
}
