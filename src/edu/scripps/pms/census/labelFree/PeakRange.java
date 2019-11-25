/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.scripps.pms.census.labelFree;

/**
 *
 * @author harshil
 */
public class PeakRange {
    private int start=0;
    private int end =0;

    public PeakRange(int start,int end) {
        this.start = start;
        this.end = end;
    }

    public int getEnd() {
        return end;
    }

    public void setEnd(int end) {
        this.end = end;
    }

    public int getStart() {
        return start;
    }

    public void setStart(int start) {
        this.start = start;
    }

    @Override
    public String toString() {
        
       return  String.valueOf(start)+ ":" + String.valueOf(end);
    }
    public boolean isInRange(int value) {
        return value >= start && value <= end;
    }
    
    
}
