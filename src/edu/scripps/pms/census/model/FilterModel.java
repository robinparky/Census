/*
 * FilterModel.java
 *
 * Created on July 11, 2006, 10:36 AM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package edu.scripps.pms.census.model;

/**
 *
 * @author rpark
 */
public class FilterModel {

    private boolean detSelect;
    private double detValue;
    private boolean pValueSelect;
    private double pValue;
    private boolean filterFragmentIons;
    private boolean removeNegative;
    private boolean uniquePeptide;
    
    
    /** Creates a new instance of FilterModel */
    public FilterModel() {
    }

    public boolean isDetSelect() {
        return detSelect;
    }

    public void setDetSelect(boolean detSelect) {
        this.detSelect = detSelect;
    }

    public double getDetValue() {
        return detValue;
    }

    public void setDetValue(double detValue) {
        this.detValue = detValue;
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

    public boolean isFilterFragmentIons() {
        return filterFragmentIons;
    }

    public void setFilterFragmentIons(boolean filterFragmentIons) {
        this.filterFragmentIons = filterFragmentIons;
    }

    public boolean isRemoveNegative() {
        return removeNegative;
    }

    public void setRemoveNegative(boolean removeNegative) {
        this.removeNegative = removeNegative;
    }

    public boolean isUniquePeptide() {
        return uniquePeptide;
    }

    public void setUniquePeptide(boolean uniquePeptide) {
        this.uniquePeptide = uniquePeptide;
    }
    
}
