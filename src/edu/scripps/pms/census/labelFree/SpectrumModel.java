/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.scripps.pms.census.labelFree;

import gnu.trove.TDoubleArrayList;
import gnu.trove.TIntArrayList;

/**
 *
 * @author harshil
 */
public class SpectrumModel {
    private double[] mass;
    private double[] intensity;
    private int[] csArray = null;
    private boolean hasCsArray = false;
    private double retentionTime;
    private double ionInjectionTime;
    private int scanNumber;
    private long precursorPeakIntensity;

    private TDoubleArrayList massArrayList;
    private TIntArrayList intArrayList;

    
    public SpectrumModel() {
        
    }

    public SpectrumModel(double[] mass, double[] intensity, double retentionTime, double ionInjectionTime, int scanNumber) {
        this.mass = mass;
        this.intensity = intensity;
        this.retentionTime = retentionTime;
        this.ionInjectionTime = ionInjectionTime;
        this.scanNumber = scanNumber;
    }
    
    /** Creates a new instance of IrisDataModel */
    
    
    public double[] getMass() {
        return mass;
    }

    public void setMass(double[] mass) {
        this.mass = mass;
    }

    public double[] getIntensity() {
        return intensity;
    }

    public void setIntensity(double[] intensity) {
        this.intensity = intensity;
    }

    public double getRetentionTime() {
        return retentionTime;
    }

    public void setRetentionTime(double retentionTime) {
        this.retentionTime = retentionTime;
    }

    public double getIonInjectionTime() {
        return ionInjectionTime;
    }

    public void setIonInjectionTime(double ionInjectionTime) {
        this.ionInjectionTime = ionInjectionTime;
    }

    public int getScanNumber() {
        return scanNumber;
    }

    public void setScanNumber(int scanNumber) {
        this.scanNumber = scanNumber;
    }

    public long getPrecursorPeakIntensity() {
        return precursorPeakIntensity;
    }

    public void setPrecursorPeakIntensity(long precursorPeakIntensity) {
        this.precursorPeakIntensity = precursorPeakIntensity;
    }


    public TIntArrayList getIntArrayList() {
        return intArrayList;
    }

    public void setIntArrayList(TIntArrayList intArrayList) {
        this.intArrayList = intArrayList;
    }

    public TDoubleArrayList getMassArrayList() {
        return massArrayList;
    }

    public void setMassArrayList(TDoubleArrayList massArrayList) {
        this.massArrayList = massArrayList;
    }

    public int[] getCsArray() {
        return csArray;
    }

    public void setCsArray(int[] csArray) {
        this.csArray = csArray;
    }

    public boolean isHasCsArray() {
        return hasCsArray;
    }

    public void setHasCsArray(boolean hasCsArray) {
        this.hasCsArray = hasCsArray;
    }
}
