/*
 * BestFragIon.java
 *
 * Created on October 3, 2005, 6:55 PM
 */

package edu.scripps.pms.census.model;

/**
 *
 * @author rpark
 */
public class BestFragIon {

    private int index;
    private long[] sArr;
    private long[] rArr;
    
    /** Creates a new instance of BestFragIon */
    public BestFragIon(int index, long[] sArr, long[] rArr) {
        this.setIndex(index);
        this.setSArr(sArr);
        this.setRArr(rArr);
    }

    public int getIndex() {
        return index;
    }

    public void setIndex(int index) {
        this.index = index;
    }

    public long[] getSArr() {
        return sArr;
    }

    public void setSArr(long[] sArr) {
        this.sArr = sArr;
    }

    public long[] getRArr() {
        return rArr;
    }

    public void setRArr(long[] rArr) {
        this.rArr = rArr;
    }
    
}
