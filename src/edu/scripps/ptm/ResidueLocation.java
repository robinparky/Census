/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package edu.scripps.ptm;

/**
 *
 * @author rpark
 */
public class ResidueLocation {
    
    private char residue;
    private int relativeLocation;
    private final double lineLength=400;

    public char getResidue() {
        return residue;
    }

    public void setResidue(char residue) {
        this.residue = residue;
    }

    public int getRelativeLocation() {
        return relativeLocation;
    }

    public void setRelativeLocation(int relativeLocation) {
        this.relativeLocation = relativeLocation;
    }
    
    
    
    
    
}
