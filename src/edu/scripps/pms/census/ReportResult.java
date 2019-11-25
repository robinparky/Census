package edu.scripps.pms.census;

public class ReportResult {
	private StringBuffer result = new StringBuffer();
	private StringBuffer singletonResult = new StringBuffer();

	int totalCount=0;
	int quantifiedCount=0;
	int redunProteinCount=0;
	int uniqueProteinCount=0;
	int proteinGroupCount=0;
	int quantifiedCountWithSingleton=0;


	public StringBuffer getResult() {
	    return result;
	}

	public void setResult(StringBuffer result) {
	    this.result = result;
	}

	public StringBuffer getSingletonResult() {
	    return singletonResult;
	}

	public void setSingletonResult(StringBuffer singletonResult) {
	    this.singletonResult = singletonResult;
	}

	public int getTotalCount() {
	    return totalCount;
	}

	public void setTotalCount(int totalCount) {
	    this.totalCount = totalCount;
	}

	public int getQuantifiedCount() {
	    return quantifiedCount;
	}

	public void setQuantifiedCount(int quantifiedCount) {
	    this.quantifiedCount = quantifiedCount;
	}

	public int getRedunProteinCount() {
	    return redunProteinCount;
	}

	public void setRedunProteinCount(int redunProteinCount) {
	    this.redunProteinCount = redunProteinCount;
	}

	public int getUniqueProteinCount() {
	    return uniqueProteinCount;
	}

	public void setUniqueProteinCount(int uniqueProteinCount) {
	    this.uniqueProteinCount = uniqueProteinCount;
	}

	public int getProteinGroupCount() {
	    return proteinGroupCount;
	}

	public void setProteinGroupCount(int proteinGroupCount) {
	    this.proteinGroupCount = proteinGroupCount;
	}


	public int getQuantifiedCountWithSingleton() {
		return quantifiedCountWithSingleton;
	}

	public void setQuantifiedCountWithSingleton(int quantifiedCountWithSingleton) {
		this.quantifiedCountWithSingleton = quantifiedCountWithSingleton;
	}
}


