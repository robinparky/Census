
/*
* Copyright (c) 2008 The Scripps Research Institute, Yates Lab.  All rights reserved.  
*/

package edu.scripps.pms.census.model;

/**
 *
 * @author Sung Kyu, Robin, Park
 * @email rpark@scripps.edu
 * Created on May 23, 2008 
 * $Revision: 1.1 $
 * $Date: 2008/09/09 22:31:05 $
 */

import java.util.*;

public class ChroLabelfreeMergeProteinModel {

    private List<ChroProtein> proteinList;
    private List<PeptideArray> peptideList;
    
    public ChroLabelfreeMergeProteinModel() {
        
    }

    public List<ChroProtein> getProteinList() {
        return proteinList;
    }

    public void setProteinList(List<ChroProtein> proteinList) {
        this.proteinList = proteinList;
    }

    public List<PeptideArray> getPeptideList() {
        return peptideList;
    }

    public void setPeptideList(List<PeptideArray> peptideList) {
        this.peptideList = peptideList;
    }
    
    class PeptideArray {
        
        private List<ChroPeptide> peptideArr;

        public List<ChroPeptide> getPeptideArr() {
            return peptideArr;
        }

        public void setPeptideArr(List<ChroPeptide> peptideArr) {
            this.peptideArr = peptideArr;
        }
        
        
    }
}
