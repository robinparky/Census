package edu.scripps.pms.census;

import java.util.*;
import edu.scripps.pms.census.model.*;
import edu.scripps.pms.census.conf.Configuration;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import javax.swing.*;

public class ReportParam {
    private ArrayList<ChroProtein> proteinList;
    private Configuration conf;
    private boolean noFilter;
    private boolean isGui;
    private JProgressBar aJProgressBar;
    private double detValue;
    private boolean filterFragmentIons;
    private double correctFactorValue;
    private boolean discardUnlabeledPeptide;
    private boolean discardReverseProtein = true;

    private boolean removeNegative;
    private boolean detSelect = true;
    private boolean isUniquePeptide = false;
    private boolean pValueSelect;
    private double pValue;
    private double allNoneLowerBound;
    private double allNoneUpperBound;
    private double allNoneCompositeScore;
    private int allNoneMinPeptide;
    private boolean smoothingPeaks=false;
    private boolean iterateOutlier=true;
    private boolean useProfileScore=false;
    private double profileScore=0.0;
    private boolean discardAN;
    private int maxSpectrumShift = 0;
    private int intensityThresholdApproach=0;  //0 don't use, 1 statistical approach (z score), 2 hard threshold
    private double intensityThreshold=0;

    private String n15Source="peptide";
    private double n15ApeThreshold=0.0;
    private double n15ProfileThreshold=0.5;
  //  private double n15RatioThreshold=0.0;
    private double n15OutlierThreshold=0.1;
    private int minimumPeptidePerProtein=1;

    private double purityThreshold = -1;
    private double signalToNoiseThreshold = Double.MIN_VALUE;
    private int minimumNumberOfUniquePeptides = -1;

    public double getPurityThreshold() {
        return purityThreshold;
    }

    public void setPurityThreshold(double purityThreshold) {
        this.purityThreshold = purityThreshold;
    }

    public double getSignalToNoiseThreshold() {
        return signalToNoiseThreshold;
    }

    public void setSignalToNoiseThreshold(double signalToNoiseThreshold) {
        this.signalToNoiseThreshold = signalToNoiseThreshold;
    }

    public boolean isUseProfileScore() {
	return useProfileScore;
    }

    public void setUseProfileScore(boolean useProfileScore) {
	this.useProfileScore = useProfileScore;
    }

    public boolean isDiscardAN() {
	return discardAN;
    }

    public void setDiscardAN(boolean discardAN) {
	this.discardAN = discardAN;
    }

    public boolean isNoFilter() {
	return noFilter;
    }

    public void setNoFilter(boolean noFilter) {
	this.noFilter = noFilter;
    }

    public boolean isIsGui() {
	return isGui;
    }

    public void setIsGui(boolean isGui) {
	this.isGui = isGui;
    }

    public boolean isFilterFragmentIons() {
	return filterFragmentIons;
    }

    public void setFilterFragmentIons(boolean filterFragmentIons) {
	this.filterFragmentIons = filterFragmentIons;
    }

    public double getCorrectFactorValue() {
	return correctFactorValue;
    }

    public void setCorrectFactorValue(double correctFactorValue) {
	this.correctFactorValue = correctFactorValue;
    }

    public boolean isDiscardUnlabeledPeptide() {
	return discardUnlabeledPeptide;
    }

    public void setDiscardUnlabeledPeptide(boolean discardUnlabeledPeptide) {
	this.discardUnlabeledPeptide = discardUnlabeledPeptide;
    }

    public ArrayList<ChroProtein> getProteinList() {
	return proteinList;
    }

    public void setProteinList(ArrayList<ChroProtein> proteinList) {
	this.proteinList = proteinList;
    }

    public Configuration getConf() {
	return conf;
    }

    public void setConf(Configuration conf) {
	this.conf = conf;
    }

    public JProgressBar getAJProgressBar() {
	return aJProgressBar;
    }

    public void setAJProgressBar(JProgressBar aJProgressBar) {
	this.aJProgressBar = aJProgressBar;
    }

    public boolean isRemoveNegative() {
	return removeNegative;
    }

    public void setRemoveNegative(boolean removeNegative) {
	this.removeNegative = removeNegative;
    }

    public boolean isDetSelect() {
	return detSelect;
    }

    public void setDetSelect(boolean detSelect) {
	this.detSelect = detSelect;
    }

    public boolean isIsUniquePeptide() {
	return isUniquePeptide;
    }

    public void setIsUniquePeptide(boolean isUniquePeptide) {
	this.isUniquePeptide = isUniquePeptide;
    }

    public boolean isPValueSelect() {
	return pValueSelect;
    }

    public void setPValueSelect(boolean pValueSelect) {
	this.pValueSelect = pValueSelect;
    }

    public double getPValue() {
	return pValue;
    }

    public void setPValue(double pValue) {
	this.pValue = pValue;
    }

    public double getDetValue() {
	return detValue;
    }

    public void setDetValue(double detValue) {
	this.detValue = detValue;
    }

    public double getAllNoneLowerBound() {
	return allNoneLowerBound;
    }

    public void setAllNoneLowerBound(double allNoneLowerBound) {
	this.allNoneLowerBound = allNoneLowerBound;
    }

    public double getAllNoneUpperBound() {
	return allNoneUpperBound;
    }

    public void setAllNoneUpperBound(double allNoneUpperBound) {
	this.allNoneUpperBound = allNoneUpperBound;
    }

    public double getAllNoneCompositeScore() {
	return allNoneCompositeScore;
    }

    public void setAllNoneCompositeScore(double allNoneCompositeScore) {
	this.allNoneCompositeScore = allNoneCompositeScore;
    }

    public int getAllNoneMinPeptide() {
	return allNoneMinPeptide;
    }

    public void setAllNoneMinPeptide(int allNoneMinPeptide) {
	this.allNoneMinPeptide = allNoneMinPeptide;
    }

    public boolean isDiscardReverseProtein() {
        return discardReverseProtein;
    }

    public void setDiscardReverseProtein(boolean discardReverseProtein) {
        this.discardReverseProtein = discardReverseProtein;
    }

    public boolean isSmoothingPeaks() {
        return smoothingPeaks;
    }

    public void setSmoothingPeaks(boolean smoothingPeaks) {
        this.smoothingPeaks = smoothingPeaks;
    }

    public boolean isIterateOutlier() {
        return iterateOutlier;
    }

    public void setIterateOutlier(boolean iterateOutlier) {
        this.iterateOutlier = iterateOutlier;
    }

    public double getProfileScore() {
	return profileScore;
    }

    public void setProfileScore(double profileScore) {
	if(profileScore>0) setUseProfileScore(true);

	this.profileScore = profileScore;
    }

    public int getMaxSpectrumShift() {
        return maxSpectrumShift;
    }

    public void setMaxSpectrumShift(int maxSpectrumShift) {
        this.maxSpectrumShift = maxSpectrumShift;
    }

    public int getIntensityThresholdApproach() {
        return intensityThresholdApproach;
    }

    public void setIntensityThresholdApproach(int intensityThresholdApproach) {
        this.intensityThresholdApproach = intensityThresholdApproach;
    }

    public double getIntensityThreshold() {
        return intensityThreshold;
    }

    public void setIntensityThreshold(double intensityThreshold) {
        this.intensityThreshold = intensityThreshold;
    }

    public void addN15Params(String paramPath) throws Exception {
        String n15Param = paramPath.replace("census.param", "census_n15.param");

        File f = new File(n15Param);

        String eachLine;

        if(f.exists()) {
            BufferedReader fbr = new BufferedReader(new FileReader(n15Param));
            Hashtable<String, String> paramHt = new Hashtable<String, String>();
            while((eachLine=fbr.readLine()) != null)  {
		if("".equals(eachLine.trim())) continue;

                String[] arr = eachLine.split("\t");
                paramHt.put(arr[0].trim(), arr[1].trim());
            }

            String value = paramHt.get("source");
            if(null != value) this.n15Source = value;

            value = paramHt.get("ape_threshold");
            if(null != value) this.n15ApeThreshold= Double.parseDouble(value);

            value = paramHt.get("outlier_threshold");
            if(null != value) this.n15OutlierThreshold=Double.parseDouble(value);

            value = paramHt.get("profile_threshold");
            if(null != value) this.n15ProfileThreshold=Double.parseDouble(value);

          //  value = paramHt.get("ratio_threshold");
          //  if(null != value) this.n15RatioThreshold=Double.parseDouble(value);
        }
    }

    public JProgressBar getaJProgressBar() {
        return aJProgressBar;
    }

    public void setaJProgressBar(JProgressBar aJProgressBar) {
        this.aJProgressBar = aJProgressBar;
    }

    public boolean ispValueSelect() {
        return pValueSelect;
    }

    public void setpValueSelect(boolean pValueSelect) {
        this.pValueSelect = pValueSelect;
    }

    public double getpValue() {
        return pValue;
    }

    public void setpValue(double pValue) {
        this.pValue = pValue;
    }

    public String getN15Source() {
        return n15Source;
    }

    public void setN15Source(String n15Source) {
        this.n15Source = n15Source;
    }

    public double getN15ApeThreshold() {
        return n15ApeThreshold;
    }

    public void setN15ApeThreshold(double n15ApeThreshold) {
        this.n15ApeThreshold = n15ApeThreshold;
    }

    public double getN15ProfileThreshold() {
        return n15ProfileThreshold;
    }

    public void setN15ProfileThreshold(double n15ProfileThreshold) {
        this.n15ProfileThreshold = n15ProfileThreshold;
    }

    public double getN15OutlierThreshold() {
        return n15OutlierThreshold;
    }

    public void setN15OutlierThreshold(double n15OutlierThreshold) {
        this.n15OutlierThreshold = n15OutlierThreshold;
    }

    public int getMinimumPeptidePerProtein() {
      return minimumPeptidePerProtein;
    }

    public void setMinimumPeptidePerProtein(int minimumPeptidePerProtein) {
      this.minimumPeptidePerProtein = minimumPeptidePerProtein;
    }

    public int getMinimumNumberOfUniquePeptides() {
        return minimumNumberOfUniquePeptides;
    }

    public void setMinimumNumberOfUniquePeptides(int minimumNumberOfUniquePeptides) {
        this.minimumNumberOfUniquePeptides = minimumNumberOfUniquePeptides;
    }


}
