/*
 * ChroiTRAQLabelData.java
 *
 * Created on August 30, 2006, 11:31 AM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package edu.scripps.pms.census.model;

/**
 *
 * @author rpark
 */
public class ChroiTRAQLabelData extends ChroData {
    
    private int scanNum;
    private long[] intensityArr;
    
    /** Creates a new instance of ChroiTRAQLabelData*/
    public ChroiTRAQLabelData() {
        
    }

    public void setFullScanData(String[] arr)
    {
        scanNum = (int)Double.parseDouble(arr[0]);        
        intensityArr = new long[arr.length-1];        
        
        for(int i=0;i<intensityArr.length;i++)
            intensityArr[i] = Long.parseLong(arr[i+1]);
    }

    public int getScanNum() {
        return scanNum;
    }

    public void setScanNum(int scanNum) {
        this.scanNum = scanNum;
    }

    public long[] getIntensityArr() {
        return intensityArr;
    }

    public void setIntensityArr(long[] intensityArr) {
        this.intensityArr = intensityArr;
    }


    
}
