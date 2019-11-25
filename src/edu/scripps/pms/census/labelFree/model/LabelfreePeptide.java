/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.scripps.pms.census.labelFree.model;

import edu.scripps.pms.census.conf.Configuration;
import edu.scripps.pms.census.hash.IndexedFile;
import edu.scripps.pms.census.model.ChroPeptide;
import edu.scripps.pms.census.util.IsotopeDist;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.LinkedHashMap;
import java.util.List;

/**
 *
 * @author rpark
 */
public class LabelfreePeptide {
    private List<ChroPeptide> peptideList = new ArrayList<ChroPeptide>();
    private double startRetTime=10000000;
    private double endRetTime=-1;
    private String sequence;
    private int chargeState;
    private List<Double> normIntensityList =new ArrayList<>();
    private double [] spectra;
    private IsotopeDist isotopeDist;
    private long [] chromeSpectra =null;
    private double [] isoArrCS;
   // private boolean findRetTime=false;
    
    /*
    public void addPeptide(ChroPeptide peptide) {
        this.peptideList.add(peptide);
        
        findRetRange();
    }*/

    public IsotopeDist getIsotopeDist() {
        return isotopeDist;
    }

    public void setIsotopeDist(IsotopeDist isotopeDist) {
        this.isotopeDist = isotopeDist;
    }

    public double[] getSpectra() {
        return spectra;
    }

    public void setSpectra(double[] spectra) {
        this.spectra = spectra;
    }

    public void addChroPeptide(ChroPeptide cp) {
        this.peptideList.add(cp);
    }

    public List<Double> getNormIntensityList() {
        return normIntensityList;
    }

    public void setNormIntensityList(List<Double> normIntensityList) {
        this.normIntensityList = normIntensityList;
    }
    
    public void addNormIntensityList(Double d){
        this.normIntensityList.add(d);
    }
    
    public void findRetRange() {
        Configuration conf = Configuration.getInstance();
        //Hashtable<String, IndexedFile> indexHt = conf.getIndexHt();
        
        //System.out.println("---");
        for(ChroPeptide pep:this.peptideList) {
         //   IndexedFile iFile = indexHt.get(pep.getFileName() + ".ms1");
            //pep.getRetentionTime()
            if(pep.getScanNum()<=0) continue;
            double low = pep.getRetentionTime()-conf.getRetentionTimeWindow();
            double high = pep.getRetentionTime()+conf.getRetentionTimeWindow();
            
            if(low<startRetTime)
                this.startRetTime=low;
            if(high>endRetTime)
                this.endRetTime=high;
            
        //    System.out.println(low + " " + high + "  + ret ==" + pep.getScanNum() + " " + pep.getRetentionTime());
            
        }
     //   System.out.println(this.startRetTime + " " + this.endRetTime);
        
    }
    
    public List<ChroPeptide> getPeptideList() {
        return peptideList;
    }

    public void setPeptideList(List<ChroPeptide> peptideList) {
        
        this.peptideList = peptideList;
        
    }

    public double getStartRetTime() {
        
        return startRetTime;
    }

    public void setStartRetTime(double startRetTime) {
        this.startRetTime = startRetTime;
    }

    public double getEndRetTime() {
        
        return endRetTime;
    }

    public void setEndRetTime(double endRetTime) {
        this.endRetTime = endRetTime;
    }

    public String getSequence() {
        return sequence;
    }

    public void setSequence(String sequence) {
        this.sequence = sequence;
    }

    public int getChargeState() {
        return chargeState;
    }

    public void setChargeState(int chargeState) {
        this.chargeState = chargeState;
    }

    public long[] getChromeSpectra() {
        return chromeSpectra;
    }

    public void setChromeSpectra(long[] chromeSpectra) {
        this.chromeSpectra = chromeSpectra;
    }

    public double[] getIsoArrCS() {
        return isoArrCS;
    }

    public void setIsoArrCS(double[] isoArrCS) {
        this.isoArrCS = isoArrCS;
    }
}
