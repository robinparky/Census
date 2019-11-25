/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package edu.scripps.pms.census.labelFree;

import edu.scripps.pms.census.tools.IsotopeModel;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 *
 * @author Harshil
 */
public class MissedPeptide {
    private String sequence;
    private double startRange;// its a retention time
    private double endRange; // it is a retention time
    private int experimentIndex;// starts with 0..
    private int cs;
    private List<Double> originalRetentionTime = new ArrayList<>();
    private List<IsotopeModel> isoTopeList = new ArrayList<IsotopeModel>();
    private String fileName = null;
    private int basePeakScanNumber = -1;
    private double totalIntensitySum = 0;
    
    public MissedPeptide(String sequence, double startRange, double endRange, int experimentIndex, double originalRetentionTime) {
        this.sequence = sequence.split("_")[0];
        this.cs = Integer.parseInt(sequence.split("_")[1]);
        this.startRange = startRange;
        this.endRange = endRange;
        this.experimentIndex = experimentIndex;
        this.originalRetentionTime.add(originalRetentionTime);
    }
    
    public void merge(MissedPeptide newPeptide)
    {
        this.setStartRange(Math.max(this.getStartRange(), newPeptide.getStartRange()));
        this.setEndRange(Math.min(this.getEndRange(), newPeptide.getEndRange()));
        getOriginalRetentionTime().addAll(newPeptide.getOriginalRetentionTime());
        
    }
    
    public String getSequenceOnly()
    {
        return this.getSequence().split("[.]")[1];
    }
    
    public double getRetentionEndTime(double tolerance)
    {
        double endtime = Collections.max(getOriginalRetentionTime());
        return endtime + tolerance;
    }
            
    public double getRetentionStartTime(double tolerance)
    {
        double startTime=Collections.min(getOriginalRetentionTime());
        startTime -= tolerance;
        if(startTime<0)
            startTime=0;
        return startTime;
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

    /**
     * @return the start Retention time
     */
    public double getStartRange() {
        
        return startRange;
    }

    /**
     * @param startRange the startRange to set
     */
    public void setStartRange(double startRange) {
        this.startRange = startRange;
    }

    /**
     * @return the endRange : the end retention time
     */
    public double getEndRange() {
        return endRange;
    }

    /**
     * @param endRange the end retention time
     */
    public void setEndRange(double endRange) {
        this.endRange = endRange;
    }

    /**
     * @return the experimentIndex
     */
    public int getExperimentIndex() {
        return experimentIndex;
    }

    /**
     * @param experimentIndex the experimentIndex to set
     */
    public void setExperimentIndex(int experimentIndex) {
        this.experimentIndex = experimentIndex;
    }

    /**
     * @return the cs
     */
    public int getCs() {
        return cs;
    }

    /**
     * @param cs the cs to set
     */
    public void setCs(int cs) {
        this.cs = cs;
    }

    /**
     * @return the originalRetentionTime
     */
    public List<Double> getOriginalRetentionTime() {
        return originalRetentionTime;
    }

    /**
     * @param originalRetentionTime the originalRetentionTime to set
     */
    public void setOriginalRetentionTime(List<Double> originalRetentionTime) {
        this.originalRetentionTime = originalRetentionTime;
    }
    
       /**
     * @return the isoTopeList
     */
    public List<IsotopeModel> getIsoTopeList() {
        return isoTopeList;
    }

    /**
     * @param isoTopeList the isoTopeList to set
     */
    public void setIsoTopeList(List<IsotopeModel> isoTopeList) {
        this.isoTopeList = isoTopeList;
       
    }
    
    /**
     * generates the BasePeak ScanNumber form the list of isotopes....
     * Sets the basePeakScanNumber value.....
     * 
     */
    public void generateBasePeakScan()
    {
        double max_total = 0;
        int maxScanNumber = 0;
        for (IsotopeModel isoModel : isoTopeList) 
        {

            if (max_total < isoModel.getIntensitySum()) {
                max_total = isoModel.getIntensitySum();
                maxScanNumber = isoModel.getScanNumber();
            }
        }
        this.basePeakScanNumber = maxScanNumber;
    }

    public String getFileName() {
        return fileName;
    }

    public void setFileName(String fileName) {
        this.fileName = fileName;
    }

    public int getBasePeakScanNumber() {
        return basePeakScanNumber;
    }

    public void setBasePeakScanNumber(int basePeakScanNumber) {
        this.basePeakScanNumber = basePeakScanNumber;
    }

    public double getTotalIntensitySum() {
        return totalIntensitySum;
    }

    public void setTotalIntensitySum(double totalIntensitySum) {
        this.totalIntensitySum = totalIntensitySum;
    }

}
