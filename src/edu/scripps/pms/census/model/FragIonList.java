/*
 * FragIonList.java
 *
 * Created on October 15, 2005, 3:36 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package edu.scripps.pms.census.model;

import java.util.Vector;

/**
 *
 * @author rpark
 */
public class FragIonList extends Vector<FragIon> {
    
    /** Creates a new instance of FragIonList */
    private int bestIndex;
    
    public FragIonList() {
        super();
    }

    public int getBestIndex() {
        return bestIndex;
    }

    public void setBestIndex(int bestIndex) {
        this.bestIndex = bestIndex;
    }
    
}
