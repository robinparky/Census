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

public class TMTOutlierModel {
    //outlier model for tmt to discard spectra for each peptides - useful for modification analysis
    
    private String sequence;
    private List<ChroPeptide> pepList = new ArrayList<ChroPeptide>();
    private double[] referencePattern = null;
    private boolean isCalculated=false;
    
    public TMTOutlierModel(int refSize) {
        referencePattern = new double[refSize];
    }
            
    public void addPeptide(ChroPeptide pep) {
        this.pepList.add(pep);
    }
    
    public double[] getReferencePattern(double intensityThreshold) {
        
        if(isCalculated) return referencePattern;
        this.isCalculated= true;
        
        for(Iterator<ChroPeptide> itr=pepList.iterator(); itr.hasNext(); ) {
            ChroPeptide each = itr.next();
            
            long[] larr = each.getTotalIntArr();
            double avgPepInt = 0;

            long tmpSum =0;
            for(long l:larr)
                tmpSum +=l;

            //double tmpAvg = tmpSum/larr.length;

            if(tmpSum<=0 || tmpSum<=intensityThreshold) continue;

            for(int i=0;i<larr.length;i++) {                            
                referencePattern[i] += larr[i]/(double)tmpSum;

                //System.out.println("----" + i + " " + referencePattern[i] + " " + larr[i] + " " + tmpSum);
            }          
            
            
            for(int i=0;i<referencePattern.length;i++) {
                referencePattern[i] = referencePattern[i]/pepList.size();
            }

         
        }
        
        //for(double d:referencePattern)
        //    System.out.println(d);
        
         //  System.out.println("----------");
        return referencePattern;
    }
}
