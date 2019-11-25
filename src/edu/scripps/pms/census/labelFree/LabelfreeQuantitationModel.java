/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package edu.scripps.pms.census.labelFree;

/*
* Copyright (c) 2008 Integrated Proteomics Applications.  All rights reserved.  
*/


//import com.ipa.ip2.model.LabelfreeSampleAnova;


/**
 *
 * @author Sung Kyu, Robin, Park
 * @email robinparky@yahoo.com
 * Created on Oct 19, 2009 
 * $Revision: 1.6 $
 * $Date: 2012/05/11 18:35:44 $
 */

import java.util.*;

import org.apache.commons.math.stat.inference.*;
import org.apache.commons.math.stat.descriptive.*;
import org.apache.commons.math.distribution.*;
//import com.ipa.ip2.util.Formatter;


public class LabelfreeQuantitationModel {
    private String accession;
    private String description;
//    private List<LabelfreeSampleAnova> sampleList = new ArrayList<LabelfreeSampleAnova>();

    private double anovaPvalue=-1;
    private double anovaPvalueLog10=-1;
    private double anovaFvalue=-1;

    private double specAnovaPvalue=-1;
    private double specAnovaFvalue=-1;
//    private boolean dataUpdate = false;

    private ArrayList<Double> ratioPvalueList = new ArrayList<Double>();
    private ArrayList<String> specCountList = new ArrayList<String>();
    private ArrayList<String> normSpecCountList = new ArrayList<String>();
    private ArrayList<String> normIntensityList = new ArrayList<String>();
    private ArrayList<String> normIntensityCorrectList = new ArrayList<String>();
    private ArrayList<Double> logRatioChangeList = new ArrayList<Double>();//about norm intersity log value,
    																									//log2 ration change 1/n=log2(norm intensity avg1/norm intensity avgn)
    private ArrayList<Double> normSpecCountlogRatioChangeList = new ArrayList<Double>();//about spec Count log value,
    private ArrayList<Double> normSpecCountRatioChangeList = new ArrayList<Double>();//about spec Count log value,
    																									//specCount log2 ration change 1/n=log2(spec count avg1/spec count avgn)
    private ArrayList<Double> correctedIntensityAverageList = new ArrayList<Double>();
 
private ArrayList<Double> normSpecCountAvgList= new ArrayList<Double>();
       
   private ArrayList<Double> specCountAvgList= new ArrayList<Double>();   
    private Double gscore=-1.0;
    private Double gscorePvalue=-1.0;
    private Double gscorePvalueLog10=-1.0;
    private Double specCountTscore=-1.0;
    private Double specCountTscoreLog10=-1.0;

    //Added by Harshil Shah
         
    private String sortCode = null;
    private String sortSCFillered = null;
/*
    public void updateStatValues() throws Exception {
        List intenClasses = new ArrayList();
        List specClasses = new ArrayList();

        
        for(Iterator<LabelfreeSampleAnova> itr=sampleList.iterator(); itr.hasNext(); ) {
            LabelfreeSampleAnova s = itr.next();

            if(s.getNormIntensityArr().length<2)
                return;


            intenClasses.add(s.getNormIntensityArr());
            specClasses.add(s.getNormSpecCountArr());



        }

        //System.out.println(TestUtils.oneWayAnovaPValue(intenClasses));
        //System.out.println(TestUtils.oneWayAnovaFValue(intenClasses));
        //System.out.println(TestUtils.oneWayAnovaPValue(specClasses));
        //System.out.println(TestUtils.oneWayAnovaFValue(specClasses));

        
        double tmpScore = TestUtils.oneWayAnovaPValue(intenClasses);
        if(!Double.isNaN(tmpScore) && !Double.isInfinite(tmpScore))
            this.anovaPvalue = Formatter.round(tmpScore , 5);

        tmpScore = TestUtils.oneWayAnovaFValue(intenClasses);
        if(!Double.isNaN(tmpScore) && !Double.isInfinite(tmpScore))
            this.anovaFvalue = Formatter.round(tmpScore , 5);

        tmpScore = TestUtils.oneWayAnovaPValue(specClasses);
        if(!Double.isNaN(tmpScore) && !Double.isInfinite(tmpScore))
            this.specAnovaPvalue = Formatter.round(tmpScore , 5);

        tmpScore = TestUtils.oneWayAnovaFValue(specClasses);
        if(!Double.isNaN(tmpScore) && !Double.isInfinite(tmpScore))
            this.specAnovaFvalue = Formatter.round(tmpScore , 5);

        dataUpdate = true;
    }*/


    public double getSpecAnovaFvalue() throws Exception {

        return specAnovaFvalue;
    }

    public void setSpecAnovaFvalue(double specAnovaFvalue) {
        this.specAnovaFvalue = specAnovaFvalue;
    }

    public double getSpecAnovaPvalue() throws Exception {

        return specAnovaPvalue;
    }

    public void setSpecAnovaPvalue(double specAnovaPvalue) {
        this.specAnovaPvalue = specAnovaPvalue;
    }
   public double getAnovaPvalueLog10() throws Exception {

        return anovaPvalueLog10;
    }

    public void setAnovaPvalueLog10(double anovaPvalueLog10) {
        this.anovaPvalueLog10 = anovaPvalueLog10;
    }

    public double getAnovaFvalue() throws Exception  {

        return anovaFvalue;
    }

    public void setAnovaFvalue(double anovaFvalue) {
        this.anovaFvalue = anovaFvalue;
    }

    public double getAnovaPvalue() throws Exception {

        return anovaPvalue;
    }

    public void setAnovaPvalue(double anovaPvalue) {
        this.anovaPvalue = anovaPvalue;
    }
    

//    public List<LabelfreeSampleAnova> getSampleList() {
//        return sampleList;
//    }
//
//
//
//    public void setSampleList(List<LabelfreeSampleAnova> sampleList) {
//        this.sampleList = sampleList;
//    }

    /*
    public void addSample(String name,
            int[] specCountArr,
            double[] normSpecCountArr,
            double[] normIntensityArr) {

        LabelfreeSampleAnova s = new LabelfreeSampleAnova();
        s.setName(name);
        s.setSpecCountArr(specCountArr);
        s.setNormSpecCountArr(normSpecCountArr);
        s.setNormIntensityArr(normIntensityArr);

        this.sampleList.add(s);

    }
*/



    public String getAccession() {
        return accession;
    }

    public void setAccession(String accession) {
        this.accession = accession;
    }

    public String getDescription() {
        return description;
    }

    public void setDescription(String description) {
        this.description = description;
    }

    /*
    public List<Integer> getSpecCountList() {
        return specCountList;
    }

    public void setSpecCountList(List<Integer> specCountList) {
        this.specCountList = specCountList;
    }

    public List<Double> getNormIntensityList() {
        return normIntensityList;
    }

    public void setNormIntensityList(List<Double> normIntensityList) {
        this.normIntensityList = normIntensityList;
    }

    public List<Double> getNormSpecCountList() {
        return normSpecCountList;
    }

    public void setNormSpecCountList(List<Double> normSpecCountList) {
        this.normSpecCountList = normSpecCountList;
    }

    public void addSpecCount(int specCount) {
        this.specCountList.add(specCount);
    }

    public void addNormSpecCount(double normSpecCount) {
        this.normSpecCountList.add(normSpecCount);
    }

    public void addNormIntensity(double normIntensity) {
        this.normIntensityList.add(normIntensity);
    }
    */
    //private List<Integer> specCountList = new ArrayList<Integer>();
    //private List<Double> normSpecCountList = new ArrayList<Double>();
    //private List<Double> normIntensityList = new ArrayList<Double>();


    /*
    class RedundantProtein {
        //SCOUNT_12   PEP_COUNT       NORM_INTENSITY_1        NORM_INTENSITY_2        CALCUATED_PEPTIDE_NUM

        //PLINE   ACCESSION       DESCRIPTION     SCOUNT_1        SCOUNT_2        NORM_SCOUNT_1   NORM_SCOUNT_2   PEP_COUNT       NORM_INTENSITY_1        NORM_INTENSITY_2        CALCUATED_PEPTIDE_NUM


    }*/


    public void addSpecCount(String specArr) {
        this.specCountList.add(specArr);
    }

    public void addRatioPvalueList(double d) {
        this.ratioPvalueList.add(d);
    }

    public ArrayList<Double> getRatioPvalueList() {
        return ratioPvalueList;
    }

    public void setRatioPvalueList(ArrayList<Double> ratioPvalueList) {
        this.ratioPvalueList = ratioPvalueList;
    }

    public ArrayList<String> getSpecCountList() {
        return specCountList;
    }

    public void setSpecCountList(ArrayList<String> specCountList) {
        this.specCountList = specCountList;
    }

    public void addNormSpecCountList(String normSpecCountArr) {
        this.normSpecCountList.add(normSpecCountArr);
    }
    public void editNormSpecCountList(int index,String normSpecCount) {
        this.normSpecCountList.set(index, normSpecCount);
    }


    public ArrayList<String> getNormSpecCountList() {
        return normSpecCountList;
    }

    public void setNormSpecCountList(ArrayList<String> normSpecCountList) {
        this.normSpecCountList = normSpecCountList;
    }

    public void addNormIntensityList(String normIntensityArr) {
        this.normIntensityList.add(normIntensityArr);
    }

    public void addNormIntensityCorrectList(String normIntensityCorrectArr) {
        this.normIntensityCorrectList.add(normIntensityCorrectArr);
    }


    public void addLogRatioChangeList(double d) {
        this.logRatioChangeList.add(d);
    }
    public ArrayList<Double> getLogRatioChangeList() {
        return logRatioChangeList;
    }

    public void setLogRatioChangeList(ArrayList<Double> logRatioChangeList) {
        this.logRatioChangeList = logRatioChangeList;
    }
    
     public void addNormSpecCountlogRatioChangeList(double d) {
        this.normSpecCountlogRatioChangeList.add(d);
    }
    
     public ArrayList<Double> getNormSpecCountlogRatioChangeList() {
        return normSpecCountlogRatioChangeList;
    }

    public void setNormSpecCountlogRatioChangeList(ArrayList<Double> normSpecCountlogRatioChangeList) {
        this.normSpecCountlogRatioChangeList = normSpecCountlogRatioChangeList;
    }

     public void addNormSpecCountRatioChangeList(double d) {
        this.normSpecCountRatioChangeList.add(d);
    }
    
     public ArrayList<Double> getNormSpecCountRatioChangeList() {
        return normSpecCountRatioChangeList;
    }

    public void setNormSpecCountRatioChangeList(ArrayList<Double> normSpecCountRatioChangeList) {
        this.normSpecCountRatioChangeList = normSpecCountRatioChangeList;
    }



    public ArrayList<String> getNormIntensityCorrectList() {
        return normIntensityCorrectList;
    }

    public void setNormIntensityCorrectList(ArrayList<String> normIntensityCorrectList) {
        this.normIntensityCorrectList = normIntensityCorrectList;
    }

    public ArrayList<String> getNormIntensityList() {
        return normIntensityList;
    }

    public void setNormIntensityList(ArrayList<String> normIntensityList) {
        this.normIntensityList = normIntensityList;
    }

    public void addCorrectedIntensityAverageList(double d) {
        this.correctedIntensityAverageList.add(d);
    }

    public ArrayList<Double> getCorrectedIntensityAverageList() {
        return correctedIntensityAverageList;
    }

    public void setCorrectedIntensityAverageList(ArrayList<Double> correctedIntensityAverageList) {
        this.correctedIntensityAverageList = correctedIntensityAverageList;
    }

  public ArrayList<Double> getNormSpecCountAvgList() {
        return this.normSpecCountAvgList;
    }

    public void setNormSpecCountAvgList(ArrayList<Double> normSpecCountAvgList) {
        this.normSpecCountAvgList =normSpecCountAvgList;
    }
     public void addNormSpecCountAvgList(double d) {
        this.normSpecCountAvgList.add(d);
    }
     public void editNormSpecCountAvgList(int index,double value) {
        this.normSpecCountAvgList.set(index,value);
    }
   
   
    public ArrayList<Double> getSpecCountAvgList() {
        return specCountAvgList;
    }

    public void setSpecCountAvgList(ArrayList<Double> specCountAvgList) {
        this.specCountAvgList = specCountAvgList;
    }
     public void addSpecCountAvgList(double d) {
        this.specCountAvgList.add(d);
    }
 
		 public Double getGscore() {
			return gscore;
		}
		public void setGscore(Double gscore) {
			this.gscore = gscore;
		}
		public Double getGscorePvalue() {
			return gscorePvalue;
		}
		public void setGscorePvalue(Double gscorePvalue) {
			this.gscorePvalue = gscorePvalue;
		}
		public Double getGscorePvalueLog10(){
			gscorePvalueLog10 = -Math.log(gscorePvalue)/Math.log(10);
			return gscorePvalueLog10;
		}
		public void setGscorePvalueLog10(Double gpl){
			gscorePvalueLog10 = gpl;
		}
		public Double getSpecCountTscore() {
			return specCountTscore;
		}
		public void setSpecCountTscore(Double specCountTscore) {
			this.specCountTscore = specCountTscore;
		}
		public Double getSpecCountTscoreLog10() {
			specCountTscoreLog10 = -Math.log(specCountTscore)/Math.log(10);
			return specCountTscoreLog10;
		}
		public void setSpecCountTscoreLog10(Double sctl) {
			this.specCountTscoreLog10 = sctl;
		}

    public String getSortCode() {
        return sortCode;
    }

    public void setSortCode(String sortCode) {
        this.sortCode = sortCode;
    }

    public String getSortSCFillered() {
        return sortSCFillered;
    }

    public void setSortSCFillered(String sortSCFillered) {
        this.sortSCFillered = sortSCFillered;
    }
		
    
}
