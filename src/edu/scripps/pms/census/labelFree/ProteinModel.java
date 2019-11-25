/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package edu.scripps.pms.census.labelFree;

import edu.scripps.pms.census.conf.Configuration;
import edu.scripps.pms.census.labelFree.model.LabelfreePeptide;
import edu.scripps.pms.census.model.ChroPeptide;
import edu.scripps.pms.census.model.ChroProtein;

import java.util.*;

import rpark.statistics.CommonStat;
import rpark.statistics.TTestUtil;

/**
 *
 * @author harshil
 * Robin
 */
public class ProteinModel {
  //  private ChroProtein chroProtein = null;
    private HashMap<String, List<ChroPeptide>> peptideMap = new LinkedHashMap<>();// Key is sequence+chargeState
  //  private List<ChroPeptide> chroPeptideList=null;
    private List<LabelfreePeptide> peptideList = new ArrayList<LabelfreePeptide>();

    private List<ChroProtein> redundnatProteinList = new ArrayList<ChroProtein>();
    private double pvalue=-1;
    private double qvalue=-1;
  private double iitPvalue = -1;
    private double[] peptideMedianLogRatioArr = null;
    private List[] logRatioArr;
    private double multipleRatiopValue;
    private List<Long> bestCorrelationIntensityList=new ArrayList<>();
    private List<Long> bestCorrelationEachpeptideIntensityList=new ArrayList<>();
    private List<Long> bestCorrelationIntensityListIIT=new ArrayList<>();
    private List<Long> bestCorrelationEachpeptideIntensityListIIT=new ArrayList<>();
    private List<Double> normIntensityList =new ArrayList<>();
    private double totalPeptideIntensity;
    private double totalIntensityDividedByProteinLengthLog10;
    private int proteinLength;
    private double[] intensitySumArr;
    private List<Double> intensityEachSampleSumList =new ArrayList<>();

    private List<Double> ionInjectionTimeIntensitySumList = new ArrayList<>();
    private List<Double> ionInjectionTimeIntensityAvgGroupList = new ArrayList<>();

    private List<List<Double>> allIITPeptideList = new ArrayList<>();

    public List<Double> getNormIntensityList() {
        return normIntensityList;
    }

    public void setNormIntensityList(List<Double> normIntensityList) {
        this.normIntensityList = normIntensityList;
    }
    public void addNormIntensityList(Double d){
        this.normIntensityList.add(d);
    }

    public double getMultipleRatiopValue() {
        return multipleRatiopValue;
    }

    public void setMultipleRatiopValue(double multipleRatiopValue) {
        this.multipleRatiopValue = multipleRatiopValue;
    }

    public void addRedundantProtein(ChroProtein pmodel) {
        this.redundnatProteinList.add(pmodel);
    }


    public void addPeptideList(List<ChroPeptide> chroPeptideList)
    {
//        this.chroPeptideList = chroPeptideList;
        LabelfreePeptide lfpep = new LabelfreePeptide();

        lfpep.setPeptideList(chroPeptideList);
        lfpep.findRetRange();

        this.peptideList.add(lfpep);

        String key = null ;
        for(ChroPeptide peptide : chroPeptideList)
        {
            if(peptide.getFileName() != null)
            {
                key = peptide.getSequence()+"_"+peptide.getChargeState();
                lfpep.setSequence(peptide.getSequence());
                lfpep.setChargeState(peptide.getChargeState());
                break;
            }
        }
        peptideMap.put(key, chroPeptideList);

    }
    /**
     * @return the chroProtein
     */

    public ChroProtein getChroProtein() {
        System.err.println("Error: This function is deprecated. Don't use this anymore");

        return null; //return chroProtein;
    }

    /**
     * @param chroProtein the chroProtein to set
     *//*
    public void setChroProtein(ChroProtein chroProtein) {
        this.chroProtein = chroProtein;
    }*/

    /**
     * @return the peptideMap
     */
    public HashMap<String, List<ChroPeptide>> getPeptideMap() {
        return peptideMap;
    }

    /**
     * @param peptideMap the peptideMap to set
     */
    public void setPeptideMap(HashMap<String, List<ChroPeptide>> peptideMap) {
        this.peptideMap = peptideMap;
    }

    /**
     *
     * @return List of peptide as sequence_chargeState
     */
    public List<String> getPeptideKey()
    {
        return new ArrayList<>(peptideMap.keySet());
    }

    public List<LabelfreePeptide> getPeptideList() {
        return peptideList;
    }

    public void setPeptideList(List<LabelfreePeptide> peptideList) {
        this.peptideList = peptideList;
    }

    public List<ChroProtein> getRedundnatProteinList() {
        return redundnatProteinList;
    }

    public void setRedundnatProteinList(List<ChroProtein> redundnatProteinList) {
        this.redundnatProteinList = redundnatProteinList;
    }

    public double getPvalue() {
        return pvalue;
    }

    public void setPvalue(double pvalue) {
        this.pvalue = pvalue;
    }

    public double getQvalue() {
        return qvalue;
    }

    public void setQvalue(double qvalue) {
        this.qvalue = qvalue;
    }

    public double[] getPeptideMedianLogRatioArr() {
        return peptideMedianLogRatioArr;
    }

    public void setPeptideMedianLogRatioArr(double[] peptideMedianLogRatioArr) {
        this.peptideMedianLogRatioArr = peptideMedianLogRatioArr;
    }


    public List[] getLogRatioArr() {
        return logRatioArr;
    }

    public void setLogRatioArr(List<Double>[] logRatioArr) {

        this.peptideMedianLogRatioArr = new double[logRatioArr.length];

        for(int i=0;i<logRatioArr.length;i++) {
            double[] dArr = new double[logRatioArr[i].size()];

            int count=0;
            for(Double d:logRatioArr[i]) {
                dArr[count] = d;
                count++;

            }

            this.peptideMedianLogRatioArr[i] = CommonStat.getMedianValue(dArr);

        //    System.out.println("aa");
        }

        try {
            this.pvalue = TTestUtil.oneSampleTTestBasedOnLog(this.peptideMedianLogRatioArr);
            //this.pvalue = TTestUtil.oneSampleTTestBasedOnLog(this.peptideMedianLogRatioArr);
        } catch (Exception e) {
            e.printStackTrace();
        }


       // this.logRatioArr = logRatioArr;
    }

    public List<Long> getBestCorrelationIntensityList() {
        return bestCorrelationIntensityList;
    }

    public void setBestCorrelationIntensityList(List<Long> bestCorrelationIntensityList) {
        this.bestCorrelationIntensityList = bestCorrelationIntensityList;
    }

    public double getTotalPeptideIntensity() {
        return totalPeptideIntensity;
    }

    public void setTotalPeptideIntensity(double totalPeptideIntensity) {
        this.totalPeptideIntensity = totalPeptideIntensity;
    }


    public List<Long> getBestCorrelationEachpeptideIntensityList() {
        return bestCorrelationEachpeptideIntensityList;
    }

    public void setBestCorrelationEachpeptideIntensityList(List<Long> bestCorrelationEachpeptideIntensityList) {
        this.bestCorrelationEachpeptideIntensityList = bestCorrelationEachpeptideIntensityList;
    }

    public List<Long> getBestCorrelationEachpeptideIntensityListIIT() {
        return bestCorrelationEachpeptideIntensityListIIT;
    }

    public void setBestCorrelationEachpeptideIntensityListIIT(List<Long> bestCorrelationEachpeptideIntensityListIIT) {
        this.bestCorrelationEachpeptideIntensityListIIT = bestCorrelationEachpeptideIntensityListIIT;
    }

    public List<Long> getBestCorrelationIntensityListIIT() {
        return bestCorrelationIntensityListIIT;
    }

    public void setBestCorrelationIntensityListIIT(List<Long> bestCorrelationIntensityListIIT) {
        this.bestCorrelationIntensityListIIT = bestCorrelationIntensityListIIT;
    }
    public void addBestCorrelationIntensityListIIT(long num){
        this.bestCorrelationIntensityListIIT.add(num);
    }
    public void addBestCorrelationEachpeptideIntensityListIIT(long num){
        this.bestCorrelationEachpeptideIntensityListIIT.add(num);
    }
    public void addBestCorrelationEachpeptideIntensityList(long num){
        this.bestCorrelationEachpeptideIntensityList.add(num);
    }
    public void addBestCorrelationIntensityList(long num){
        this.bestCorrelationIntensityList.add(num);
    }

    public double getTotalIntensityDividedByProteinLengthLog10() {
        return totalIntensityDividedByProteinLengthLog10;
    }

    public void setTotalIntensityDividedByProteinLengthLog10(double totalIntensityDividedByProteinLengthLog10) {
        this.totalIntensityDividedByProteinLengthLog10 = totalIntensityDividedByProteinLengthLog10;
    }

    public int getProteinLength() {
        return proteinLength;
    }

    public void setProteinLength(int proteinLength) {
        this.proteinLength = proteinLength;
    }

    public double[] getIntensitySumArr() {
        return intensitySumArr;
    }

    public void setIntensitySumArr(double[] intensitySumArr) {
        this.intensitySumArr = intensitySumArr;
    }

    public void setIonInjectionTimeIntensitySumList(List<Double> ionInjectionTimeIntensitySumList)
    {
        this.ionInjectionTimeIntensitySumList = ionInjectionTimeIntensitySumList;
    }

    public List<Double> getIonInjectionTimeIntensitySumList()
    {
        return ionInjectionTimeIntensitySumList;
    }

  public List<Double> getIntensityEachSampleSumList() {
    return intensityEachSampleSumList;
  }

  public void setIntensityEachSampleSumList(List<Double> intensityEachSampleSumList) {
    this.intensityEachSampleSumList = intensityEachSampleSumList;
  }

  public void addIntensityEachSampleSumList(Double intensity){
    this.intensityEachSampleSumList.add(intensity);
  }

  public double getIitPvalue() {
    return iitPvalue;
  }

  public void setIitPvalue(double iitPvalue) {
    this.iitPvalue = iitPvalue;
  }

  public List<Double> getIonInjectionTimeIntensityAvgGroupList() {
    return ionInjectionTimeIntensityAvgGroupList;
  }

  public void setIonInjectionTimeIntensityAvgGroupList(List<Double> ionInjectionTimeIntensityAvgGroupList) {
    this.ionInjectionTimeIntensityAvgGroupList = ionInjectionTimeIntensityAvgGroupList;
  }

  public List<List<Double>> getAllIITPeptideList() {
    return allIITPeptideList;
  }

  public void setAllIITPeptideList(List<List<Double>> allIITPeptideList) {
    this.allIITPeptideList = allIITPeptideList;
  }

  public void addEachIITPeptideToList(int sampleIndex, double intensity) {

    if(this.allIITPeptideList.size()<=sampleIndex) {

      this.allIITPeptideList.add(sampleIndex, new ArrayList<Double>());
    }


    List<Double> iitPeptideList = this.allIITPeptideList.get(sampleIndex);
    iitPeptideList.add(intensity);
  }

  public List<Double> getEachIITPeptideList(int sampleIndex) {
    return allIITPeptideList.get(sampleIndex);
  }



}
