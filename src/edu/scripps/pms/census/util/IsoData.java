/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package edu.scripps.pms.census.util;

/**
 *
 * @author Harshil
 */
public class IsoData {

    private int scanNumber;
    private double[] lightData;
    private double[] heavyData;

//    public IsoData( int scanNumber, String[] lightData, String[] heavyData) {
//        this.lightData = new double[lightData.length];
//        this.heavyData = new double[heavyData.length];
//        for(int i =0;i<lightData.length;i++)
//        {
//            this.lightData[i]= Double.parseDouble(lightData[i]);
//            this.heavyData[i]= Double.parseDouble(heavyData[i]);
//        }
//        this.scanNumber = scanNumber;
//    }
    public IsoData(int scanNumber, double[] lightData, double[] heavyData) {
        this.lightData = lightData;
        this.heavyData = heavyData;
        this.scanNumber = scanNumber;
    }

    IsoData() {
    }
    
    
    public static void getIsoObject(double[] arr)
    {
        
    }
    

    public double[] getLightData() {
        return lightData;
    }

    public void setLightData(double[] data1) {
        this.lightData = data1;
    }

    public double[] getHeavyData() {
        return heavyData;
    }

    public void setHeavyData(double[] data2) {
        this.heavyData = data2;
    }

    public int getScanNumber() {
        return scanNumber;
    }

    public void setScanNumber(int scanNumber) {
        this.scanNumber = scanNumber;
    }
    
}
