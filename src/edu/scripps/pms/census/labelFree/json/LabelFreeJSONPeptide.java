package edu.scripps.pms.census.labelFree.json;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class LabelFreeJSONPeptide {

	private String auc = "";

	private String rowId = "";
	
	private String count = "";

	private String unique = "";
	
	private String MHplus  = "";
	
	private String PROFILE_SCORE  = "";
	
	private String calcMHplus  = "";
	
	private String mass_list  = "";
	
	private String intensity  = "";
	
	private String spScore  = "";
	
	private String deltaCN  = "";
	
	private String end_scan  = "";
	
	private String ionInjectionTime  = "";
	
	private String heavyAvgMass  = "";
	
	private String file  = "";
	
	private String spRank  = "";
	
	private String deltaMass  = "";
	
	private String xcorr  = "";
	
	private String redundancy  = "";
	
	private String start_scan  = "";
	
	private String charge  = "";
	
	private String lightAvgMass  = "";
	
	private String chro_iso  = "";
	
	private String peaks ="";
	
	private String seq  = "";
	
	private String scan  = "";
	
	private String retentionTime  = "";
	
	private String missed  = "";
	
	private String spC  = "";
	
	private String startRt = "";
	
	private String endRt = "";

	private String sampleName = "";
	
	private String[] maxPeakValue;
	private String groupName = "";

	private double[] isotopeArr;
	private int rank;

	public String getStartRt() {
		return startRt;
	}

	public void setStartRt(String startRt) {
		this.startRt = startRt;
	}

	public String getEndRt() {
		return endRt;
	}

	public void setEndRt(String endRt) {
		this.endRt = endRt;
	}

	public String getUnique() {
		return unique;
	}

	public void setUnique(String unique) {
		this.unique = unique;
	}

	public String getMHplus() {
		return MHplus;
	}

	public void setMHplus(String mHplus) {
		MHplus = mHplus;
	}

	public String getPROFILE_SCORE() {
		return PROFILE_SCORE;
	}

	public void setPROFILE_SCORE(String pROFILE_SCORE) {
		PROFILE_SCORE = pROFILE_SCORE;
	}

	public String getCalcMHplus() {
		return calcMHplus;
	}

	public void setCalcMHplus(String calcMHplus) {
		this.calcMHplus = calcMHplus;
	}

	public String getMass_list() {
		return mass_list;
	}

	public void setMass_list(String mass_list) {
		this.mass_list = mass_list;
	}

	public String getIntensity() {
		return intensity;
	}

	public void setIntensity(String intensity) {
		this.intensity = intensity;
	}

	public String getSpScore() {
		return spScore;
	}

	public void setSpScore(String spScore) {
		this.spScore = spScore;
	}

	public String getDeltaCN() {
		return deltaCN;
	}

	public void setDeltaCN(String deltaCN) {
		this.deltaCN = deltaCN;
	}

	public String getEnd_scan() {
		return end_scan;
	}

	public void setEnd_scan(String end_scan) {
		this.end_scan = end_scan;
	}

	public String getIonInjectionTime() {
		return ionInjectionTime;
	}

	public void setIonInjectionTime(String ionInjectionTime) {
		this.ionInjectionTime = ionInjectionTime;
	}

	public String getHeavyAvgMass() {
		return heavyAvgMass;
	}

	public void setHeavyAvgMass(String heavyAvgMass) {
		this.heavyAvgMass = heavyAvgMass;
	}

	public String getFile() {
		return file != null ? file : "";
	}

	public void setFile(String file) {
		this.file = file;
	}

	public String getSpRank() {
		return spRank;
	}

	public void setSpRank(String spRank) {
		this.spRank = spRank;
	}

	public String getDeltaMass() {
		return deltaMass;
	}

	public void setDeltaMass(String deltaMass) {
		this.deltaMass = deltaMass;
	}

	public String getXcorr() {
		return xcorr;
	}

	public void setXcorr(String xcorr) {
		this.xcorr = xcorr;
	}

	public String getRedundancy() {
		return redundancy;
	}

	public void setRedundancy(String redundancy) {
		this.redundancy = redundancy;
	}

	public String getStart_scan() {
		return start_scan;
	}

	public void setStart_scan(String start_scan) {
		this.start_scan = start_scan;
	}

	public String getCharge() {
		return charge != null ? charge : "";
	}

	public void setCharge(String charge) {
		this.charge = charge;
	}

	public String getLightAvgMass() {
		return lightAvgMass;
	}

	public void setLightAvgMass(String lightAvgMass) {
		this.lightAvgMass = lightAvgMass;
	}

	public String getChro_iso() {
		return chro_iso != null ? chro_iso : "";
	}

	public void setChro_iso(String chro_iso) {
		this.chro_iso = chro_iso;
	}

	public String getSeq() {
		return seq != null ? seq : "";
	}

	public void setSeq(String seq) {
		this.seq = seq;
	}

	public String getScan() {
		return scan != null ? scan : "";
	}

	public void setScan(String scan) {
		this.scan = scan;
	}

	public String getRetentionTime() {
		return retentionTime;
	}

	public void setRetentionTime(String retentionTime) {
		this.retentionTime = retentionTime;
	}

	public String getMissed() {
		return missed;
	}

	public void setMissed(String missed) {
		this.missed = missed;
	}

	public String getSpC() {
		return spC;
	}

	public void setSpC(String spC) {
		this.spC = spC;
	}

	public String getCount() {
		return count;
	}

	public void setCount(String count) {
		this.count = count;
	}

	public String getRowId() {
		return rowId;
	}

	public void setRowId(String rowId) {
		this.rowId = rowId;
	}

	public String getPeaks() {
		return peaks;
	}

	public void setPeaks(String peaks) {
		this.peaks = peaks;
	}

	public String[] getMaxPeakValue() {
		return maxPeakValue;
	}

	public void setMaxPeakValue(String[] maxPeakValue) {
		this.maxPeakValue = maxPeakValue;
	}

	public double[] getIsotopeArr() {
		return isotopeArr;
	}

	public void setIsotopeArr(double[] isotopeArr) {
		this.isotopeArr = Arrays.copyOf(isotopeArr,isotopeArr.length);
	}

	public String getAuc() {
		return auc;
	}

	public void setAuc(String auc) {
		this.auc = auc;
	}


	public String getSampleName() {
		return sampleName;
	}

	public void setSampleName(String sampleName) {
		this.sampleName = sampleName;
	}

	public String getGroupName() {
		return groupName;
	}

	public void setGroupName(String groupName) {
		this.groupName = groupName;
	}

	public int getRank() {
		return rank;
	}

	public void setRank(int rank) {
		this.rank = rank;
	}
}
