package edu.scripps.pms.census.tandem;

import edu.scripps.pms.census.util.XYPoint;

import java.util.ArrayList;
import java.util.List;


public class TandemTagPeptide {

	private boolean unique;
	
	private String sequence = null;
	
	private String spc = null;
	
	private String scanNum = null;
	
	private String cstate = null;
	
	private String filename = null;
	
	private String startMass = "";
	
	private String endMass = "";
	
	private List<MZValues> mzValues = new ArrayList<>();
	
	private List<XYPoint> xy = new ArrayList<>();

	public boolean isUnique() {
		return unique;
	}

	public void setUnique(boolean unique) {
		this.unique = unique;
	}

	public String getSequence() {
		return sequence;
	}

	public void setSequence(String sequence) {
		this.sequence = sequence;
	}

	public String getSpc() {
		return spc;
	}

	public void setSpc(String spc) {
		this.spc = spc;
	}

	public String getScanNum() {
		return scanNum;
	}

	public void setScanNum(String scanNum) {
		this.scanNum = scanNum;
	}


	public String getFilename() {
		return filename;
	}

	public void setFilename(String filename) {
		this.filename = filename;
	}

	public List<MZValues> getMzValues() {
		return mzValues;
	}

	public void setMzValues(List<MZValues> mzValues) {
		this.mzValues = mzValues;
		try {
			this.startMass = mzValues.get(0).getMzHeader();
			this.endMass = mzValues.get(mzValues.size() -1).getMzHeader();
			
			this.startMass = this.startMass.substring(4, this.startMass.length()-4);
			this.endMass  =  this.endMass.substring(4, this.endMass.length()-4);
			
		} catch (Exception e) {
			System.out.println("Error st mass & end Mass ::"+ e.getMessage());
		}
		
	}

	public String getCstate() {
		return cstate;
	}

	public void setCstate(String cstate) {
		this.cstate = cstate;
	}

	public List<XYPoint> getXy() {
		return xy;
	}

	public void setXy(List<XYPoint> xy) {
		this.xy = xy;
	}

	public String getStartMass() {
		return startMass;
	}

	public void setStartMass(String startMass) {
		this.startMass = startMass;
	}

	public String getEndMass() {
		return endMass;
	}

	public void setEndMass(String endMass) {
		this.endMass = endMass;
	}
		
}
