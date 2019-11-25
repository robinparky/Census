/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.scripps.pms.census.model;

/**
 *
 * @author rpark
 */
import java.util.*;

public class LabelingResultModel {
    private List<IsotopeModel> resultList = new ArrayList<IsotopeModel>();
    private int scanNum;
    private double retTime;

    public List<IsotopeModel> getResultList() {
        return resultList;
    }

    public void setResultList(List<IsotopeModel> resultList) {
        this.resultList = resultList;
    }
    

    public void addResult(IsotopeModel model) {
        this.resultList.add(model);
    }
    
    public double getTotalIntensity() {
        double totalIntensity=0;
        for(Iterator<IsotopeModel> itr=resultList.iterator(); itr.hasNext(); ) {
            IsotopeModel each = itr.next();
            totalIntensity += each.getSumIntensity();            
        }
        
        return totalIntensity;
    }
        

    public String getContent() {
        StringBuffer sb = new StringBuffer();
        sb.append(scanNum).append(" ").append(this.retTime).append(" ");
        
        for(Iterator<IsotopeModel> itr=resultList.iterator(); itr.hasNext(); ) {
            IsotopeModel each = itr.next();
            sb.append(each.getSumIntensity()).append(" ");  
        }
        
        return sb.toString();
    }
    
    public int getScanNum() {
        return scanNum;
    }

    public void setScanNum(int scanNum) {
        this.scanNum = scanNum;
    }

    public double getRetTime() {
        return retTime;
    }

    public void setRetTime(double retTime) {
        this.retTime = retTime;
    }
    
    
    
}
