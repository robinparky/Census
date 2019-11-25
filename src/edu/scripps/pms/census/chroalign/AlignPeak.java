/*
 * AlignPeak.java
 *
 * Created on June 28, 2006, 10:09 AM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package edu.scripps.pms.census.chroalign;

/**
 *
 * @author jvenable
 */
public class AlignPeak implements Comparable<AlignPeak> {
    
    private double mz;
    private double intensity;
    private boolean refSpectrum;

    /** Creates a new instance of AlignPeak */
    public AlignPeak(double mz, double intensity, boolean refSpectrum) {
        this.setMz(mz);
        this.setIntensity(intensity);
        this.setRefSpectrum(refSpectrum);
    }
    
    public int compareTo(AlignPeak peak)
    {
        double temp = this.mz - peak.getMz();
        
        if(temp>0)
            return 1;
        else if(temp<0)
            return -1;
        else return 0;        
    }

    public double getMz() {
        return mz;
    }

    public void setMz(double mz) {
        this.mz = mz;
    }
    
    public double getIntensity() {
        return intensity;
    }

    public void setIntensity(double intensity) {
        this.intensity = intensity;
    }
    
    public boolean getRefSpectrum() {
        return refSpectrum;
    }

    public void setRefSpectrum(boolean refSpectrum) {
        this.refSpectrum = refSpectrum;
    }
}
