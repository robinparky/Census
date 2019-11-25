/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.scripps.pms.census.model;

/**
 *
 * @author rpark
 */
import java.util.*;

public class ReportIon {
    private double mass;
    private List<Double> cPlusList = new ArrayList<Double>();
    private List<Double> cMinusList = new ArrayList<Double>();

    public ReportIon() { }
        
    public ReportIon(double mass, List cPlusList, List<Double> cMinusList) {
        this.mass = mass;
        this.cPlusList = cPlusList;
        this.cMinusList = cMinusList;
    }
    
    public double getMass() {
        return mass;
    }

    public void setMass(double mass) {
        this.mass = mass;
    }

    public List<Double> getcPlusList() {
        return cPlusList;
    }

    public void setcPlusList(String str) {
        String[] arr = str.split(",");
        for(String each:arr)            
            this.cPlusList.add(new Double(each));
    }

    public void setcPlusList(List<Double> cPlusList) {
        this.cPlusList = cPlusList;
    }

    public List<Double> getcMinusList() {
        return cMinusList;
    }

    public void setcMinusList(String str) {
        String[] arr = str.split(",");
        for(String each:arr)            
            this.cMinusList.add(new Double(each));
    }

    public void setcMinusList(List<Double> cMinusList) {
        this.cMinusList = cMinusList;
    }
    
    
}
