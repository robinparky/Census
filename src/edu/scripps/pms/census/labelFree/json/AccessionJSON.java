package edu.scripps.pms.census.labelFree.json;

import java.util.ArrayList;
import java.util.List;

public class AccessionJSON {
	
	private LabelFreeJSONProtein protein;
	
	private List<List<LabelFreeJSONPeptide>> peptideList = new ArrayList<>();

	public LabelFreeJSONProtein getProtein() {
		return protein;
	}

	public void setProtein(LabelFreeJSONProtein protein) {
		this.protein = protein;
	}

	public List<List<LabelFreeJSONPeptide>> getPeptideList() {
		return peptideList;
	}

	public void setPeptideList(List<List<LabelFreeJSONPeptide>> peptideList) {
		this.peptideList = peptideList;
	}

}
