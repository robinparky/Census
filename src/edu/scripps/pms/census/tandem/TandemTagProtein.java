package edu.scripps.pms.census.tandem;

import java.util.ArrayList;
import java.util.List;

public class TandemTagProtein {

	private String locus = null;

	private String specCount = null;

	private String pepNum = null;

	private String average = null;

	private String description = null;
	
	private List<TandemTagPeptide> peptideList = new ArrayList<>();

	public List<TandemTagPeptide> getPeptideList() {
		return peptideList;
	}

	public void setPeptideList(List<TandemTagPeptide> peptideList) {
		this.peptideList = peptideList;
	}

	public String getLocus() {
		return locus;
	}

	public void setLocus(String locus) {
		this.locus = locus;
	}

	public String getSpecCount() {
		return specCount;
	}

	public void setSpecCount(String specCount) {
		this.specCount = specCount;
	}

	public String getPepNum() {
		return pepNum;
	}

	public void setPepNum(String pepNum) {
		this.pepNum = pepNum;
	}

	public String getAverage() {
		return average;
	}

	public void setAverage(String average) {
		this.average = average;
	}

	public String getDescription() {
		return description;
	}

	public void setDescription(String description) {
		this.description = description;
	}

}
