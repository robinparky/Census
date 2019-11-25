/*
 * FragIon.java
 *
 * Created on October 3, 2005, 6:55 PM
 */

package edu.scripps.pms.census.model;

/**
 *
 * @author rpark
 */
public class FragIon implements Comparable {

    private int index;
    private double regScore;
    private boolean isBion;

    private long[] sArr;
    private long[] rArr;
    
    /** Creates a new instance of FragIon */
    public FragIon(int index, long[] sArr, long[] rArr, double regScore, boolean isBion) {
        this.setIndex(index);
        this.setSArr(sArr);
        this.setRArr(rArr);
        this.setRegScore(regScore);
        this.setIsBion(isBion);
    }

    public int getIndex() {
        return index;
    }

    public void setIsBion(boolean isBion) {
        this.isBion = isBion;
    }

    public boolean isBion() {
        return this.isBion;
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

    public void setRegScore(double regScore) {
        this.regScore = regScore;
    }

    public double getRegScore() {
        return this.regScore;
    }
    
    public int compareTo(Object o) throws ClassCastException {

	FragIon ion = (FragIon)o;

        double temp = this.regScore - ion.getRegScore();

        if(temp<0)
            return -1;
        else if(temp>0)
            return 1;
        else 
            return 0;
    }
}
