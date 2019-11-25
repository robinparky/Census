/*
 * ChroNonLabelData.java
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
public class ChroNonLabelData extends ChroData {
    
    private int[] scanNumArr;
    private long[] intensityArr;
    
    /** Creates a new instance of ChroNonLabelData */
    public ChroNonLabelData() {
        
    }

    public void setFullScanData(String[] arr)
    {
        scanNumArr = new int[arr.length/2];
        intensityArr = new long[scanNumArr.length];
        
        for(int i=0;i<scanNumArr.length;i++)
        {
            scanNumArr[i] = Integer.parseInt(arr[i]);
            intensityArr[i] = Long.parseLong(arr[i+intensityArr.length]);
        }        

    }



    //public ChroNonLabelData(int[] scanArr, long[] intensityArr) {        
    //}

    public int[] getScanNumArr() {
        return scanNumArr;
    }

    public void setScanNumArr(int[] scanNumArr) {
        this.scanNumArr = scanNumArr;
    }

    public long[] getIntensityArr() {
        return intensityArr;
    }

    public void setIntensityArr(long[] intensityArr) {
        this.intensityArr = intensityArr;
    }

    
}
