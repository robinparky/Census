/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.scripps.pms.census.model;

import edu.scripps.pms.census.model.ChroData;
/**
 *
 * @author rpark
 */
public class CompositeScoreModel implements Comparable {

    private double intRatioLog;
    private boolean outlier=false;
    
    
    
    public CompositeScoreModel(double intRatioLog) {
        this.intRatioLog = intRatioLog;
    //    this.cdata = cdata;
    }
    
    @Override
    public int compareTo(Object o) {
        
        CompositeScoreModel c = (CompositeScoreModel)o;
        final int BEFORE =-1;
        final int EQUALS =0;
        final int AFTER =1;
        
        if(this.intRatioLog < c.getIntRatioLog()) return BEFORE;
        else if(this.intRatioLog ==c.getIntRatioLog()) return EQUALS;
        else if(this.intRatioLog > c.getIntRatioLog()) return AFTER;
        else
            return 0;
    }

    public double getIntRatioLog() {
        return intRatioLog;
    }

    public void setIntRatioLog(double intRatioLog) {
        this.intRatioLog = intRatioLog;
    }

    public boolean isOutlier() {
        return outlier;
    }

    public void setOutlier(boolean outlier) {
        this.outlier = outlier;
    }
}
