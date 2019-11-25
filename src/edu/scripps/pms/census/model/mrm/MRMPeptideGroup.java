/*
 * MRMPeptideGroup.java
 *
 * Created on September 10, 2007, 2:16 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package edu.scripps.pms.census.model.mrm;

import java.util.*;

/**
 *
 * @author rpark
 */
public class MRMPeptideGroup {

    List<MRMPeptideModel> peptideList = new ArrayList<MRMPeptideModel>();
    
    /** Creates a new instance of MRMPeptideGroup */
    public MRMPeptideGroup() {
    }
     
    public void addPeptide(MRMPeptideModel peptide)
    {
        this.peptideList.add(peptide);
    }

    public List<MRMPeptideModel> getPeptideList() {
        return peptideList;
    }

    public void setPeptideList(List<MRMPeptideModel> peptideList) {
        this.peptideList = peptideList;
    }
    
}
