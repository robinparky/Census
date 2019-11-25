/*
 * PostCalculation.java
 *
 * Created on June 9, 2005, 4:01 PM
 */

package edu.scripps.pms.census.util;

import edu.scripps.pms.census.model.ChroData;
import edu.scripps.pms.census.model.ChroPeptide;

import java.util.*;

/**
 *
 * @author rpark
 */
public class PostCalculation {
    
    private ChroPeptide peptide;
    private int scanNum;
    private double[] samInt;
    private double[] refInt;
    private int size;
    
    /** Creates a new instance of PostCalculation */
    public PostCalculation(ChroPeptide peptide, double[] samInt, double[] refInt) {
        this.peptide = peptide;
        //this.scanNum = Integer.parseInt(peptide.getScanNum());
        
        this.samInt = samInt;
        this.refInt = refInt;
        
//        init();
    }
    /*
    private void init()
    {
        int i=0;
        
        for(Iterator<ChroData> itr=peptide.getData().iterator(); itr.hasNext(); )
        {
            ChroData data = itr.next();
        }
    }
    */
    
}
