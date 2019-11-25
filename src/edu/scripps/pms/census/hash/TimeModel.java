/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.scripps.pms.census.hash;

/**
 *
 * @author harshil
 */
public class TimeModel {
    
    private double retentionTime ;
    private double ionInjectionTime;
    private int scanNumber;

    public TimeModel(double retentionTime, double ionInjectionTime, int scanNumber) {
        this.retentionTime = retentionTime;
        this.ionInjectionTime = ionInjectionTime;
        this.scanNumber = scanNumber;
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
    
}
