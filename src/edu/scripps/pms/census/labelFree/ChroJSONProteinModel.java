/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package edu.scripps.pms.census.labelFree;

import edu.scripps.pms.census.model.ChroPeptide;
import edu.scripps.pms.census.model.ChroProtein;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;

/**
 *
 * @author Harshil
 * @version $Id:
 */
public class ChroJSONProteinModel {

    private ChroProtein protein = null;
    private List<List<ChroPeptide>> peptideGroup= new ArrayList<>();//each will have alist of 4(Depending on experiments) peptides
    private List<ChroProtein> redundantProtein = new ArrayList<>();
    private List<List<Double>> peptideRatioIntensity = new ArrayList<>();
    
    private double medianIntensity = -1;

    
    protected static final String PEPTIDE_LIST = "peptideList";
    protected static final String REDUNDANT_PROTEIN = "redundantProtein";
    protected static final String PROTEIN = "protein";
    public ChroJSONProteinModel() {
        protein = new ChroProtein();
        peptideGroup = new ArrayList<>();
    }
    
    
    public void readJsonObject(JSONObject proteinObj)
    {
        if(proteinObj.containsKey(PROTEIN))
            this.protein = ChroProtein.readJsonObject((JSONObject)proteinObj.get(PROTEIN));
        if(proteinObj.containsKey(PEPTIDE_LIST))
        {
            for(int i=0;i<((JSONArray) proteinObj.get(PEPTIDE_LIST)).size();i++)
            {
                JSONArray arraqy = (JSONArray) ((JSONArray) proteinObj.get(PEPTIDE_LIST)).get(i);
                List<ChroPeptide> pepList = new ArrayList<>();
                {
                    for(int j=0;j<arraqy.size();j++)
                    {
                        ChroPeptide peptide = ChroPeptide.readJsonObject((JSONObject) arraqy.get(j));
                        pepList.add(peptide);
                    }
                }
                this.peptideGroup.add(pepList);
            }            
        }
        if(proteinObj.containsKey(REDUNDANT_PROTEIN))
        {
            JSONArray array = (JSONArray) proteinObj.get(REDUNDANT_PROTEIN);
            for(int i =0;i<array.size();i++)
            {
                this.redundantProtein.add(ChroProtein.readJsonObject((JSONObject) array.get(i)));
            }
        }
            
      
//        
//        if(proteinObj.containsKey("val"))
//            protein.setValidation((String) proteinObj.get("val"));
//        if(proteinObj.containsKey("molwt"))
//            protein.setMolWt((String) proteinObj.get("molwt"));
//        if(proteinObj.containsKey("desc"))
//            protein.setDescription((String) proteinObj.get("desc"));
//        if(proteinObj.containsKey("length"))
//            protein.setLength((String) proteinObj.get("length"));
//        if(proteinObj.containsKey("seq_ct"))
//            protein.setSeqCount((String) proteinObj.get("seq_ct"));
//        if(proteinObj.containsKey("accession"))
//            protein.setLocus((String) proteinObj.get("accession"));
//        
//        if(proteinObj.containsKey("peptides"))
//        {
//            JSONArray pepList =  (JSONArray) proteinObj.get("peptides");
//
//            for(Iterator expetimentItr = pepList.iterator();expetimentItr.hasNext();)
//            {
//                JSONArray experiementPeptide = (JSONArray) expetimentItr.next();
//                List<ChroPeptide> peptideList = new ArrayList<>();
//                for(Iterator peptideItr = experiementPeptide.iterator();peptideItr.hasNext();)
//                {
//                    ChroPeptide peptide = new ChroPeptide();
//                    peptide.readJsonObject((JSONObject) peptideItr.next());
//                    peptideList.add(peptide);
//                }
//                this.peptideGroup.add(peptideList);
//            }
//        }
    }

    public ChroProtein getProtein() {
        return protein;
    }

    public void setProtein(ChroProtein protein) {
        this.protein = protein;
    }

    public List<List<ChroPeptide>> getPeptideGroup() {
        return peptideGroup;
    }

    public void setPeptideGroup(List<List<ChroPeptide>> peptideGroup) {
        this.peptideGroup = peptideGroup;
    }
    
    public void addPeptideGroup(List<ChroPeptide> peptideGroup) {
        this.peptideGroup.add(peptideGroup);
    }
    
    public JSONObject getJSONObject()
    {
        JSONObject jsonObj = new JSONObject();
        jsonObj.put(PROTEIN, this.protein.getJSONObj());
        JSONArray array = new JSONArray();
        for(List<ChroPeptide> peptideList : this.peptideGroup)
        {
            JSONArray exprimentArray = new JSONArray();
            for(ChroPeptide peptide : peptideList)
            {
                exprimentArray.add(peptide.getJSONObj());
            }
            
            array.add(exprimentArray);
        }
        jsonObj.put(PEPTIDE_LIST, array);
        JSONArray proteinArray = new JSONArray();
        for(ChroProtein redundantProtein :this.redundantProtein)
        {
            proteinArray.add(redundantProtein.getJSONObj());
        }
        jsonObj.put(REDUNDANT_PROTEIN, proteinArray);
        return jsonObj;
    }


    public List<ChroProtein> getRedundantProtein() {
        return redundantProtein;
    }

    public void setRedundantProtein(List<ChroProtein> redundantProtein) {
        this.redundantProtein = redundantProtein;
    }
    public void addRedundantProtein(ChroProtein redundantProtein) {
        this.redundantProtein.add(redundantProtein);
    }

    /**
     * Peptide level does 1-1 group ratio of avg intensity.
     * 1/0 -- 2/0 -- 3/0....etc
     */
    public List<List<Double>> getPeptideRatioIntensity() {
        return peptideRatioIntensity;
    }

    public void setPeptideRatioIntensity(List<List<Double>> peptideRatioIntensity) {
        this.peptideRatioIntensity = peptideRatioIntensity;
    }
    
    /**
     * Peptide level does 1-1 group ratio of avg intensity.
     * 1/0 -- 2/0 -- 3/0....etc
     */
    public void addPeptideRatioIntensity(List<Double> peptideRatioList) {
        this.peptideRatioIntensity.add(peptideRatioList);
    }

    public double getMedianIntensity() {
        return medianIntensity;
    }

    public void setMedianIntensity(double medianIntensity) {
        this.medianIntensity = medianIntensity;
    }
}
