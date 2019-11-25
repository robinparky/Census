/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.scripps.pms.census.model;

/**
 *
 * @author rpark
 */
public class IsotopeModel {
    double sumIntensity;
    int foundIsoNum; 
    int totalPeakFound;
    double toleranceSum;
    double[] isoArr;
    double[] isoFoundArr;
    
    public double getSumIntensity() {
        return sumIntensity;
    }

    public void setSumIntensity(double sumIntensity) {
        this.sumIntensity = sumIntensity;
    }    

    public int getFoundIsoNum() {
        return foundIsoNum;
    }

    public void setFoundIsoNum(int foundIsoNum) {
        this.foundIsoNum = foundIsoNum;
    }

    public int getTotalPeakFound() {
        return totalPeakFound;
    }

    public void setTotalPeakFound(int totalPeakFound) {
        this.totalPeakFound = totalPeakFound;
    }

    public double getToleranceSum() {
        return toleranceSum;
    }

    public void setToleranceSum(double toleranceSum) {
        this.toleranceSum = toleranceSum;
    }

    public double[] getIsoArr() {
        return isoArr;
    }

    public void setIsoArr(double[] isoArr) {
        this.isoArr = isoArr;
    }

    public double[] getIsoFoundArr() {
        return isoFoundArr;
    }

    public void setIsoFoundArr(double[] isoFoundArr) {
        this.isoFoundArr = isoFoundArr;
    }
    
    
    
}
