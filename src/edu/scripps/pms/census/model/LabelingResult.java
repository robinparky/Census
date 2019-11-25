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
public class LabelingResult {
    private int peakStart;
    private int peakEnd;
    private int leftIndex;
    private int rightIndex;
    LabelingResultModel[] resultArr;

    public int getPeakStart() {
        return peakStart;
    }

    public void setPeakStart(int peakStart) {
        this.peakStart = peakStart;
    }

    public int getPeakEnd() {
        return peakEnd;
    }

    public void setPeakEnd(int peakEnd) {
        this.peakEnd = peakEnd;
    }

    public int getLeftIndex() {
        return leftIndex;
    }

    public void setLeftIndex(int leftIndex) {
        this.leftIndex = leftIndex;
    }

    public int getRightIndex() {
        return rightIndex;
    }

    public void setRightIndex(int rightIndex) {
        this.rightIndex = rightIndex;
    }

    public LabelingResultModel[] getResultArr() {
        return resultArr;
    }

    public void setResultArr(LabelingResultModel[] resultArr) {
        this.resultArr = resultArr;
    }
    
    
    
}
