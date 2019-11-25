/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package edu.scripps.ptm;

import edu.scripps.pms.util.seq.Fasta;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 *
 * @author Harshil
 */
public class PeptidePhosphoModel {
    
    
    private String unique=null;
    private String fileName ;
    private double xCorr;
    private double deltaCN;
    private double conf;
    private double mh;
    private double calcMh;
    private double ppm;
    private double totalIntensity;
    private int spr;
    private double probScore;
    private double pi;
    private double ionProportion;
    private int redundancy;
    private String sequence;
    private List<String> proteinAccession = new ArrayList<>();
    private List<String> proteinDescription = new ArrayList<>();
//    private String proteinAccession;
//    private String proteinDescription;
    private String modSequence;
    private List<Double> localizationScore = new ArrayList<>();
    private double debunkerScore;

    /**
     * @return the unique
     */
    public String getUnique() {
        return unique;
    }

    /**
     * @param unique the unique to set
     */
    public void setUnique(String unique) {
        this.unique = unique;
    }

    /**
     * @return the fileName
     */
    public String getFileName() {
        return fileName;
    }

    /**
     * @param fileName the fileName to set
     */
    public void setFileName(String fileName) {
        this.fileName = fileName;
    }

    /**
     * @return the xCorr
     */
    public double getxCorr() {
        return xCorr;
    }

    /**
     * @param xCorr the xCorr to set
     */
    public void setxCorr(double xCorr) {
        this.xCorr = xCorr;
    }

    /**
     * @return the deltaCN
     */
    public double getDeltaCN() {
        return deltaCN;
    }

    /**
     * @param deltaCN the deltaCN to set
     */
    public void setDeltaCN(double deltaCN) {
        this.deltaCN = deltaCN;
    }

    /**
     * @return the conf
     */
    public double getConf() {
        return conf;
    }

    /**
     * @param conf the conf to set
     */
    public void setConf(double conf) {
        this.conf = conf;
    }

    /**
     * @return the mh
     */
    public double getMh() {
        return mh;
    }

    /**
     * @param mh the mh to set
     */
    public void setMh(double mh) {
        this.mh = mh;
    }

    /**
     * @return the calcMh
     */
    public double getCalcMh() {
        return calcMh;
    }

    /**
     * @param calcMh the calcMh to set
     */
    public void setCalcMh(double calcMh) {
        this.calcMh = calcMh;
    }

    /**
     * @return the ppm
     */
    public double getPpm() {
        return ppm;
    }

    /**
     * @param ppm the ppm to set
     */
    public void setPpm(double ppm) {
        this.ppm = ppm;
    }

    /**
     * @return the totalIntensity
     */
    public double getTotalIntensity() {
        return totalIntensity;
    }

    /**
     * @param totalIntensity the totalIntensity to set
     */
    public void setTotalIntensity(double totalIntensity) {
        this.totalIntensity = totalIntensity;
    }

    /**
     * @return the spr
     */
    public int getSpr() {
        return spr;
    }

    /**
     * @param spr the spr to set
     */
    public void setSpr(int spr) {
        this.spr = spr;
    }

    /**
     * @return the probScore
     */
    public double getProbScore() {
        return probScore;
    }

    /**
     * @param probScore the probScore to set
     */
    public void setProbScore(double probScore) {
        this.probScore = probScore;
    }

    /**
     * @return the pi
     */
    public double getPi() {
        return pi;
    }

    /**
     * @param pi the pi to set
     */
    public void setPi(double pi) {
        this.pi = pi;
    }

    /**
     * @return the ionProportion
     */
    public double getIonProportion() {
        return ionProportion;
    }

    /**
     * @param ionProportion the ionProportion to set
     */
    public void setIonProportion(double ionProportion) {
        this.ionProportion = ionProportion;
    }

    /**
     * @return the redundancy
     */
    public int getRedundancy() {
        return redundancy;
    }

    /**
     * @param redundancy the redundancy to set
     */
    public void setRedundancy(int redundancy) {
        this.redundancy = redundancy;
    }

    /**
     * @return the sequence
     */
    public String getSequence() {
        return sequence;
    }

    /**
     * @param sequence the sequence to set
     */
    public void setSequence(String sequence) {
        this.sequence = sequence;
    }

//    /**
//     * @return the proteinAccession
//     */
//    public String getProteinAccession() {
//        return proteinAccession;
//    }
//
//    /**
//     * @param proteinAccession the proteinAccession to set
//     */
//    public void setProteinAccession(String proteinAccession) {
//        this.proteinAccession = proteinAccession.split("[\\[\\] ]")[1];
//        System.out.println(this.proteinAccession);
//    }
//
//    /**
//     * @return the proteinDescription
//     */
//    public String getProteinDescription() {
//        return proteinDescription;
//    }
//
//    /**
//     * @param proteinDescription the proteinDescription to set
//     */
//    public void setProteinDescription(String proteinDescription) {
//        this.proteinDescription = proteinDescription;
//    }

    /**
     * @return the modSequence
     */
    public String getModSequence() {
        return modSequence;
    }

    /**
     * @param modSequence the modSequence to set
     */
    public void setModSequence(String modSequence) {
        this.modSequence = modSequence;
    }

    

    /**
     * @return the debunkerScore
     */
    public double getDebunkerScore() {
        return debunkerScore;
    }

    /**
     * @param debunkerScore the debunkerScore to set
     */
    public void setDebunkerScore(double debunkerScore) {
        this.debunkerScore = debunkerScore;
    }

    public List<Double> getLocalizationScore() {
        return localizationScore;
    }
    
    public Double getLocalizationScoreAvg() {
        Double sum=0.0;
        for (Double value : localizationScore)
        {
            sum += value;
        }
        
        return sum/localizationScore.size()*1.0;
    }

    public void setLocalizationScore(List<Double> localizationScore) {
        this.localizationScore = localizationScore;
    }
    
    
    
      public void addLocalizationScore(Double localizationScore) {
        this.localizationScore.add(localizationScore);
    }     

    public List<String> getProteinAccession() {
        return proteinAccession;
    }

    public void setProteinAccession(List<String> proteinAccession) {
        this.proteinAccession = proteinAccession;
    }
    public void addProteinAccession(String proteinAccession) {
        this.proteinAccession.add(Fasta.getAccession(proteinAccession.split("[\\[\\] ]")[1]));
                }

    public List<String> getProteinDescription() {
        return proteinDescription;
    }

    public void setProteinDescription(List<String> proteinDescription) {
        this.proteinDescription = proteinDescription;
    }
    public void addProteinDescription(String proteinDescription) {
        this.proteinDescription.add(proteinDescription);
    }
    
          
    
    
}
