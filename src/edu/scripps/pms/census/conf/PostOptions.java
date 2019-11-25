/*
 * PostOptions.java
 *
 * Created on September 27, 2005, 11:04 AM
 */

package edu.scripps.pms.census.conf;

/**
 *
 * @author rpark
 */
public class PostOptions {
    
    
    private static PostOptions options=null;
    
    private boolean filterFragmentIons=true;
    private boolean displayFragmentIons=false;
    private boolean detFactorCheck=true;
    private float detFactorValue=0.0f;    
    //private double bestFragIonTolerance=0.5;
    private double bestFragIonTolerance=0.7;
    
    private PostOptions() {
    }

    public static PostOptions getInstance()
    {
        if (options == null)
            options = new PostOptions();

        return options;
    }    
    
    public boolean isFilterFragmentIons() {
        return filterFragmentIons;
    }

    public void setFilterFragmentIons(boolean filterFragmentIons) {
        this.filterFragmentIons = filterFragmentIons;
    }

    public boolean isDisplayFragmentIons() {
        return displayFragmentIons;
    }

    public void setDisplayFragmentIons(boolean displayFragmentIons) {
        this.displayFragmentIons = displayFragmentIons;
    }

    public float getDetFactorValue() {
        return detFactorValue;
    }

    public void setDetFactorValue(float detFactorValue) {
        this.detFactorValue = detFactorValue;
    }

    public boolean isDetFactorCheck() {
        return detFactorCheck;
    }

    public void setDetFactorCheck(boolean detFactorCheck) {
        this.detFactorCheck = detFactorCheck;
    }

    public double getBestFragIonTolerance() {
        return bestFragIonTolerance;
    }

    public void setBestFragIonTolerance(double bestFragIonTolerance) {
        this.bestFragIonTolerance = bestFragIonTolerance;
    }
    
}
