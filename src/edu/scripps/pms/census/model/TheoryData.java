/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package edu.scripps.pms.census.model;

/**
 *
 * @author Harshil
 */
public class TheoryData {

    private double[] lightMass;
    private double[] heavyMass;
    private double[] lightIntensity;
    private double[] heavyIntensity;
    private double[] normHeavyIntensity;
    private double[] normLightIntensity;

    public double[] getNormHeavyIntensity() {
        return normHeavyIntensity;
    }

    public void setNormHeavyIntensity(double[] normHeavyIntensity) {
        this.normHeavyIntensity = normHeavyIntensity;
    }

    public double[] getNormLightIntensity() {
        return normLightIntensity;
    }

    public void setNormLightIntensity(double[] normLightIntensity) {
        this.normLightIntensity = normLightIntensity;
    }
    
    
    
    public void generateNormalizedValues()
    {
        normHeavyIntensity = new double[heavyIntensity.length];
        normLightIntensity = new double[lightIntensity.length];
        double total =0;
        for(double val : lightIntensity)
        {
            total+=val;
        }
        for(int i=0;i<lightIntensity.length;i++)
        {
            normLightIntensity[i]=total/lightIntensity[i];
        }
        total = 0;
        for(double val : heavyIntensity)
        {
            total+=val;
        }
        for(int i=0;i<heavyIntensity.length;i++)
        {
            normHeavyIntensity[i]=total/heavyIntensity[i];
        }
    }
    public double[] getLightMass() {
        return lightMass;
    }

    public void setLightMass(double[] lightMass) {
        this.lightMass = lightMass;
    }

    public double[] getHeavyMass() {
        return heavyMass;
    }

    public void setHeavyMass(double[] heavyMass) {
        this.heavyMass = heavyMass;
    }

    public double[] getLightIntensity() {
        return lightIntensity;
    }

    public void setLightIntensity(double[] lightIntensity) {
        this.lightIntensity = lightIntensity;
    }

    public double[] getHeavyIntensity() {
        return heavyIntensity;
    }

    public void setHeavyIntensity(double[] heavyIntensity) {
        this.heavyIntensity = heavyIntensity;
    }
    
    
}
