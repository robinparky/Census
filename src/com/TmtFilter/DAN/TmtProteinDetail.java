/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package com.TmtFilter.DAN;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

/**
 *
 * @author Harshil
 */
public class TmtProteinDetail 
{
    String protein = null;
    int peptideCount =0;
    private DescriptiveStatistics average[] = new DescriptiveStatistics[3];
    private DescriptiveStatistics stdDev[] = new  DescriptiveStatistics[3];
    String proteinDesc = null;

    public TmtProteinDetail() {
        for(int i =0;i<3;i++)
        {
            average[i]= new DescriptiveStatistics();
            stdDev[i]= new DescriptiveStatistics();
        }
    }
    
    public DescriptiveStatistics[] getAverage() {
        return average;
    }

    public void setAverage(DescriptiveStatistics[] average) {
        this.average = average;
    }

    public DescriptiveStatistics[] getStdDev() {
        return stdDev;
    }

    public void setStdDev(DescriptiveStatistics[] stdDev) {
        this.stdDev = stdDev;
    }
    
   

    public String getProteinDesc() {
        return proteinDesc;
    }

    public void setProteinDesc(String proteinDesc) {
        this.proteinDesc = proteinDesc;
    }
     
    public String getProtein() {
        return protein;
    }

    public void setProtein(String protein) {
        this.protein = protein;
    }

    public int getPeptideCount() {
        return peptideCount;
    }

    public void setPeptideCount(int peptideCount) {
        this.peptideCount = peptideCount;
    }
     public void addPeptideCount() {
        this.peptideCount++;
    }



    
    
}
