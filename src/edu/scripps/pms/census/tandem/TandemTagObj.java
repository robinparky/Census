package edu.scripps.pms.census.tandem;

import java.util.List;

public class TandemTagObj {
	
	private List<TandemTagProtein> proteins = null;
	
	private List<TandemTagPeptide> peptides = null;

	public List<TandemTagProtein> getProteins() {
		return proteins;
	}

	public void setProteins(List<TandemTagProtein> proteins) {
		this.proteins = proteins;
	}

	public List<TandemTagPeptide> getPeptides() {
		return peptides;
	}

	public void setPeptides(List<TandemTagPeptide> peptides) {
		this.peptides = peptides;
	}
	
	

}
