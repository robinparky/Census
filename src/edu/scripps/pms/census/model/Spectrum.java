/*
 * Spectrum.java
 *
 * Created on March 25, 2005, 3:58 PM
 */

package edu.scripps.pms.census.model;

import gnu.trove.TDoubleIntHashMap;
import gnu.trove.TDoubleArrayList;

import java.util.Arrays;

/**
 *
 * @author rpark
 */
public class Spectrum extends TDoubleIntHashMap 
{
    
    private double mw;
    
    /** Creates a new instance of Spectrum */
    public Spectrum() {
        super();
    }
  
    

/*
    private TDoubleArrayList massList = new TDoubleArrayList();
    public int binarySearch(double mw)
    {
        return massList.binarySearch(mw);        
    }    
    
    public int getIntensity(int mwKey)
    {
        return get( massList.get(mwKey) );        
    }

    public void add(String mw, String intensity)
    {
        this.mw = Double.parseDouble(mw);
        super.put(this.mw, Integer.parseInt(intensity));
        massList.add(this.mw);
    }

    public  TDoubleArrayList getMassIntensityList()
    {
	return massList;

    }

    */
    private double[] massList;

    public int binarySearch(double mw)
    {
        return Arrays.binarySearch(massList, mw);        
    }    
    
    public int getIntensity(int mwKey)
    {
        return super.get( massList[mwKey] );        
    }

    public void add(double mw, int intensity)
    {
        super.put(mw, intensity);
    }

    public double[] getMassIntensityList()
    {
	return massList;

    }
    public void setMassList(double[] massList)
    {
        this.massList = massList;
    }


}
