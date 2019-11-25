/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.scripps.pms.census.model;

import java.util.ArrayList;
import java.util.Iterator;

import java.util.List;

/**
 *
 * @author Rohan Rampuria <rampuria@scripps.edu>
 */
public class ChroPeptideModel {
    
    private String sequence;
    private int cs;
    
    public void getcPeptide(ChroPeptide currentPeptide){
        
             if(currentPeptide.getSequence() != null){
                setSequence(currentPeptide.getSequence());
            }
            if(currentPeptide.getChargeState() != 0){
                setCstate(currentPeptide.getChargeState());  
            }
            
        
    }
    
    public String getSequence() {
        return sequence;
    }

    public void setSequence(String sequence) {
        this.sequence = sequence;
    }

    public int getCstate() {
        return cs;
    }

    public void setCstate(int cs) {
        this.cs = cs;
    }
    
    
    public String getSequenceOnly()
    {
        return this.getSequence().split("[.]")[1];
    }
    
    
}
