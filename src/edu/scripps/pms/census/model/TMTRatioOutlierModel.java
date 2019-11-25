/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.scripps.pms.census.model;

import edu.scripps.pms.stats.GrubbsTest;
import java.util.ArrayList;
import java.util.List;
import java.util.Iterator;
import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;

/**
 *
 * @author rpark
 */
public class TMTRatioOutlierModel {
    
    private String sequence;
    private List<ChroPeptide> pepList = new ArrayList<ChroPeptide>();
    private double averageInt;
    private long[] sumIntArr;

    private double stdev1;
    private double mean1;
    private double stdev2;
    private double mean2;
    private double stdev3;
    private double mean3;

    private boolean calcuateMean=false;
    
    private DescriptiveStatistics stat1 = new DescriptiveStatistics();
    private DescriptiveStatistics stat2 = new DescriptiveStatistics();
    private DescriptiveStatistics stat3 = new DescriptiveStatistics();
    
    //public double getMean() {
   // }
    //first approach with Dan.  Average of two ions
    public boolean isOutlier(ChroPeptide peptide, double pValue) {
        if(!calcuateMean) {                        
            for(Iterator<ChroPeptide> itr=pepList.iterator(); itr.hasNext(); ) {
                ChroPeptide each = itr.next();                
                long[] larr = each.getTotalIntArr();
                double[] normarr = new double[larr.length];
                for(int i=0;i<normarr.length;i++) {
                    normarr[i] = larr[i] * averageInt / sumIntArr[i];
               //     System.out.println("norm " + normarr[i]);
                }

                double avg1 = (normarr[0]+normarr[3])/2;
                double avg2 = (normarr[1]+normarr[4])/2;
                double avg3 = (normarr[2]+normarr[5])/2;

                double v1 = Math.log(avg1/avg2)/Math.log(2);
                double v2 = Math.log(avg1/avg3)/Math.log(2);
                double v3 = Math.log(avg2/avg3)/Math.log(2);
                       
                stat1.addValue(v1);
                stat2.addValue(v2);
                stat3.addValue(v3);
                
        //        System.out.println(v1 + " @ " + v2 +" @ " + v3);
            
            }
            
            stdev1 = stat1.getStandardDeviation();
            mean1 = stat1.getMean();
            stdev2 = stat2.getStandardDeviation();
            mean2 = stat2.getMean();
            stdev3 = stat3.getStandardDeviation();
            mean3 = stat3.getMean();            
            
            this.calcuateMean = true;            
        }        
        
        
	if(stdev1==0 || stdev2==0 || stdev3==0) return false;    
        
        
        long[] larr = peptide.getTotalIntArr();
        double[] normarr = new double[larr.length];
        for(int i=0;i<normarr.length;i++) {
            normarr[i] = larr[i] * averageInt / sumIntArr[i];
            //System.out.println("No " + normarr[i] + " " + larr[i] + " " +  averageInt + " " + sumIntArr[i]);
        }
        
        double avg1 = (normarr[0]+normarr[3])/2;
        double avg2 = (normarr[1]+normarr[4])/2;
        double avg3 = (normarr[2]+normarr[5])/2;

        
        double v1 = Math.log(avg1/avg2)/Math.log(2); //126 & 129
        double v2 = Math.log(avg1/avg3)/Math.log(2); //127 & 130
        double v3 = Math.log(avg2/avg3)/Math.log(2); //128 & 131
        //System.out.println("norm\t" + normarr[0] + " " + normarr[1] +" " + normarr[2]);                
        //System.out.println(v1 + " " + v2 +" " + v3);
        //System.out.println(v1 + " " + v2 +" " + v3);
        double gp = GrubbsTest.getGrubbsPvalue(v1, mean1, stdev1, this.pepList.size());                         
        //System.out.println("V1\t" + v1 + "\t" + gp + "\t" +  pValue + "\t" + mean1 + "\t" + stdev1);
        if(gp<=pValue) return true;
        
        gp = GrubbsTest.getGrubbsPvalue(v2, mean2, stdev2, this.pepList.size());             
        //System.out.println("V2\t" + v2 + "\t" + gp + "\t" +  pValue + "\t" + mean2 + "\t" + stdev2);
        if(gp<=pValue) return true;
        
        gp = GrubbsTest.getGrubbsPvalue(v3, mean3, stdev3, this.pepList.size()); 
        //System.out.println("V3\t" + v3 + "\t" + gp + "\t" +  pValue + "\t" + mean3 + "\t" + stdev3);
        if(gp<=pValue) return true;

            
	
            
        return false;
    }
    
    
    //second approach 126/129, 127/130, 128/131
    public boolean isOutlier2(ChroPeptide peptide, double pValue) {
        if(!calcuateMean) {                        
            for(Iterator<ChroPeptide> itr=pepList.iterator(); itr.hasNext(); ) {
                ChroPeptide each = itr.next();                
                long[] larr = each.getTotalIntArr();
                double[] normarr = new double[larr.length];
                for(int i=0;i<normarr.length;i++) {
                    normarr[i] = larr[i] * averageInt / sumIntArr[i];
               //     System.out.println("norm " + normarr[i]);
                }

                double ratio1 = normarr[0]/normarr[3];
                double ratio2 = normarr[1]/normarr[4];
                double ratio3 = normarr[2]/normarr[5];

                double v1 = Math.log(ratio1)/Math.log(2);
                double v2 = Math.log(ratio2)/Math.log(2);
                double v3 = Math.log(ratio3)/Math.log(2);
                       
                stat1.addValue(v1);
                stat2.addValue(v2);
                stat3.addValue(v3);
                
        //        System.out.println(v1 + " @ " + v2 +" @ " + v3);
            
            }
            
            stdev1 = stat1.getStandardDeviation();
            mean1 = stat1.getMean();
            stdev2 = stat2.getStandardDeviation();
            mean2 = stat2.getMean();
            stdev3 = stat3.getStandardDeviation();
            mean3 = stat3.getMean();            
            
            this.calcuateMean = true;            
        }        
        
        
	if(stdev1==0 || stdev2==0 || stdev3==0) return false;    
        
        
        long[] larr = peptide.getTotalIntArr();
        double[] normarr = new double[larr.length];
        for(int i=0;i<normarr.length;i++) {
            normarr[i] = larr[i] * averageInt / sumIntArr[i];
            //System.out.println("No " + normarr[i] + " " + larr[i] + " " +  averageInt + " " + sumIntArr[i]);
        }
        
        
        double ratio1 = normarr[0]/normarr[3];
        double ratio2 = normarr[1]/normarr[4];
        double ratio3 = normarr[2]/normarr[5];

        double v1 = Math.log(ratio1)/Math.log(2);
        double v2 = Math.log(ratio2)/Math.log(2);
        double v3 = Math.log(ratio3)/Math.log(2);

        //System.out.println("norm\t" + normarr[0] + " " + normarr[1] +" " + normarr[2]);                
        //System.out.println(v1 + " " + v2 +" " + v3);
        //System.out.println(v1 + " " + v2 +" " + v3);
        double gp = GrubbsTest.getGrubbsPvalue(v1, mean1, stdev1, this.pepList.size());                         
        //System.out.println("V1\t" + v1 + "\t" + gp + "\t" +  pValue + "\t" + mean1 + "\t" + stdev1);
        if(gp<=pValue) return true;
        
        gp = GrubbsTest.getGrubbsPvalue(v2, mean2, stdev2, this.pepList.size());             
        //System.out.println("V2\t" + v2 + "\t" + gp + "\t" +  pValue + "\t" + mean2 + "\t" + stdev2);
        if(gp<=pValue) return true;
        
        gp = GrubbsTest.getGrubbsPvalue(v3, mean3, stdev3, this.pepList.size()); 
        //System.out.println("V3\t" + v3 + "\t" + gp + "\t" +  pValue + "\t" + mean3 + "\t" + stdev3);
        if(gp<=pValue) return true;

            
        return false;
    }
    
    
    public TMTRatioOutlierModel(double averageInt, long[] sumIntArr) {
        this.averageInt = averageInt;
        this.sumIntArr = sumIntArr;
    }
    
    public String getSequence() {
        return sequence;
    }

    public void setSequence(String sequence) {
        this.sequence = sequence;
    }

    public List<ChroPeptide> getPepList() {
        return pepList;
    }

    public void setPepList(List<ChroPeptide> pepList) {
        this.pepList = pepList;
    }

          
    public void addPeptide(ChroPeptide pep) {
        this.pepList.add(pep);
    }

    public double getAverageInt() {
        return averageInt;
    }

    public void setAverageInt(double averageInt) {
        this.averageInt = averageInt;
    }

    public long[] getSumIntArr() {
        return sumIntArr;
    }

    public void setSumIntArr(long[] sumIntArr) {
        this.sumIntArr = sumIntArr;
    }


    public boolean isCalcuateMean() {
        return calcuateMean;
    }

    public void setCalcuateMean(boolean calcuateMean) {
        this.calcuateMean = calcuateMean;
    }

    
    
}
