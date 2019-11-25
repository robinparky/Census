/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.scripps.pms.census.util;

/**
 *
 * @author rpark
 */

import java.util.*;

public class Peak 
{
    
    double mass;
    double intensity;
    double isoPeak;
    boolean valid;
    int peakMatch=0;
    Set<Integer> indexArr = new HashSet<Integer>();
    double totalIsoIntensity;

    public Peak(double mass, double intensity) 
    {
        this.mass = mass;
        this.intensity = intensity;
       
    }
    public double getMass()
    {
        return mass;
    }
    public double getIntensity()
    {
        return intensity;
    }

   public void setPeakMatch(int p) {
	peakMatch = p;
   }

   public int getPeakMatch() {
	   return peakMatch;
   }

   public void addIndex(int i) {
	indexArr.add(i);
   }

  public Set<Integer> getIndexSet() {
	return indexArr;

  }

  public void setTotalIsoIntensity(double d) {
	  totalIsoIntensity = d;
  }

public double getTotalIsoIntensity() {
	return totalIsoIntensity;
}
}
