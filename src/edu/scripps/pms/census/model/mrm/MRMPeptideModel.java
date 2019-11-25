/*
 * MRMPeptideModel.java
 *
 * Created on September 10, 2007, 1:42 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package edu.scripps.pms.census.model.mrm;

import java.util.*;
import edu.scripps.pms.census.util.dtaselect.Peptide;
/**
 *
 * @author rpark
 */
public class MRMPeptideModel extends Peptide {

    private String sequence = null;
    private double parentMass;    
    private String name = null;
    private String desc = null;
    private List<Daughter> daughters = new ArrayList<Daughter>();
    private boolean labeled;
    private String fileName;
    
    private double rtTolerance; //retention time tolerance
    private double snTolerance; //scan number tolerance
    
    private double[][] bionArr;
    private double[][] yionArr;
    
    private int keyIndex;
    //private double[][] bionHeavyArr;
    //private double[][] yionHeavyArr;    
    
    private int startScan; //scan # with biggest peak. used for peak finding.
    private double startRt; //RT # with biggest peak. used for peak finding.
    private String rt;
    
    public static class Daughter
    {
        private double rt;
        private double mass;

        public Daughter(double mass, double rt)
        {
            this.mass = mass;
            this.rt = rt;
        }
        
        public double getRt() {
            return rt;
        }

        public void setRt(double rt) {
            this.rt = rt;
        }

        public double getMass() {
            return mass;
        }

        public void setMass(double mass) {
            this.mass = mass;
        }        
    }
    
    /** Creates a new instance of MRMPeptideModel */
    public MRMPeptideModel() {
	super();
    }

    public void setRt(String rt) {
	this.rt = rt;
    }

    public String getRt() {
	return this.rt;
    }
    public void addDaughter(double mass, double rt)
    {
        this.daughters.add(new Daughter(mass, rt));
    }
    
    public String getSequence() {
        return sequence;
    }

    public void setSequence(String sequence) {
        this.sequence = sequence;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String getDesc() {
        return desc;
    }

    public void setDesc(String desc) {
        this.desc = desc;
    }

    public List<Daughter> getDaughters() {
        return daughters;
    }

    public void setDaughters(List<Daughter> daughters) {
        this.daughters = daughters;
    }

    public double getParentMass() {
        return parentMass;
    }

    public void setParentMass(double parentMass) {
        this.parentMass = parentMass;
    }

    public boolean isLabeled() {
        return labeled;
    }

    public void setLabeled(boolean labeled) {
        this.labeled = labeled;
    }

    public double getRtTolerance() {
        return rtTolerance;
    }

    public void setRtTolerance(double rtTolerance) {
        this.rtTolerance = rtTolerance;
    }

    public double getSnTolerance() {
        return snTolerance;
    }

    public void setSnTolerance(double snTolerance) {
        this.snTolerance = snTolerance;
    }

    public String getFileName() {
        return fileName;
    }

    public void setFileName(String fileName) {
        this.fileName = fileName;
    }

    public double getStartRt() {
        return startRt;
    }

    public void setStartRt(double startRt) {
        this.startRt = startRt;
    }    

    public double[][] getBionArr() {
        return bionArr;
    }

    public void setBionArr(double[][] bionArr) {
        this.bionArr = bionArr;
    }

    public double[][] getYionArr() {
        return yionArr;
    }

    public void setYionArr(double[][] yionArr) {
        this.yionArr = yionArr;
    }

    public int getKeyIndex() {
        return keyIndex;
    }

    public void setKeyIndex(int keyIndex) {
        this.keyIndex = keyIndex;
    }
}
