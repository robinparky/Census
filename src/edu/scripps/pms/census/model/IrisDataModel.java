/*
 * IrisDataModel.java
 *
 * Created on December 22, 2005, 11:24 AM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package edu.scripps.pms.census.model;

/**
 *
 * @author rpark
 */
public class IrisDataModel {
    
    private double[] mass;
    private double[] intensity;
            
    /** Creates a new instance of IrisDataModel */
    public IrisDataModel() {
    }

    public IrisDataModel(double[] mass, double[] intensity) {
        this.mass = mass;
        this.intensity = intensity;
    }

    public double[] getMass() {
        return mass;
    }

    public void setMass(double[] mass) {
        this.mass = mass;
    }

    public double[] getIntensity() {
        return intensity;
    }

    public void setIntensity(double[] intensity) {
        this.intensity = intensity;
    }
    
}
