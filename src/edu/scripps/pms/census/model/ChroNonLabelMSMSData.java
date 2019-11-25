/*
 * ChroNonLabelMSMSData.java
 *
 * Created on August 30, 2006, 11:31 AM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package edu.scripps.pms.census.model;

import java.util.*;

/**
 *
 * @author rpark
 */
public class ChroNonLabelMSMSData extends ChroData {
    
    private int[] scanArr;  //multiple scan numbers at the same retention time
    private long[] totalIntArr; //multiple intensities at the same retention time
        
    private ArrayList bList = new ArrayList();
    private ArrayList yList = new ArrayList();

    public void addBList(String[] arr)
    {
        bList.add(arr);
    }

    public void addYList(String[] arr)
    {
        yList.add(arr);
    }

    public ArrayList getBList() {
        return bList;
    }

    public ArrayList getYList() {
        return yList;
    }
        
    
    /** Creates a new instance of ChroNonLabelData */
    public ChroNonLabelMSMSData() {
        
    }

    
    public void setTandemData(String[] eachScanArr)
    {
        String[] scanTmpArr = eachScanArr[0].split(" ");
        String[] intTmpArr = eachScanArr[1].split(" ");
        this.scanArr = new int[scanTmpArr.length]; //scan num arr
        this.totalIntArr = new long[scanTmpArr.length]; //frag ion intensity sum arr
        
        for(int i=0;i<getScanArr().length;i++)
        {
            this.scanArr[i] = Integer.parseInt(scanTmpArr[i]);
            this.totalIntArr[i] = Long.parseLong(intTmpArr[i]);
            
            String[] bArr = eachScanArr[i+2].split(" ");
            String[] yArr = eachScanArr[i+3].split(" ");
            
            bList.add(bArr);
            yList.add(yArr);
        }
    }

    public int[] getScanArr() {
        return scanArr;
    }

    public void setScanArr(int[] scanArr) {
        this.scanArr = scanArr;
    }

    public long[] getTotalIntArr() {
        return totalIntArr;
    }

    public void setTotalIntArr(long[] totalIntArr) {
        this.totalIntArr = totalIntArr;
    }

    
}
