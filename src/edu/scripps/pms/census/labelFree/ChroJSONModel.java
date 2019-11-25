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
public class ChroJSONModel {

    private ChroProtein protein = null;
    private List<List<ChroPeptide>> peptideGroup= null;//each will have alist of 4 peptides

    public ChroJSONModel() {
        protein = new ChroProtein();
        peptideGroup = new ArrayList<>();
    }
    
    
    public void readJsonObject(JSONObject proteinObj)
    {
        
        if(proteinObj.containsKey("val"))
            protein.setValidation((String) proteinObj.get("val"));
        if(proteinObj.containsKey("molwt"))
            protein.setMolWt((String) proteinObj.get("molwt"));
        if(proteinObj.containsKey("desc"))
            protein.setDescription((String) proteinObj.get("desc"));
        if(proteinObj.containsKey("length"))
            protein.setLength((String) proteinObj.get("length"));
        if(proteinObj.containsKey("seq_ct"))
            protein.setSeqCount((String) proteinObj.get("seq_ct"));
        if(proteinObj.containsKey("accession"))
            protein.setLocus((String) proteinObj.get("accession"));
        
        if(proteinObj.containsKey("peptides"))
        {
            JSONArray pepList =  (JSONArray) proteinObj.get("peptides");

            for(Iterator expetimentItr = pepList.iterator();expetimentItr.hasNext();)
            {
                JSONArray experiementPeptide = (JSONArray) expetimentItr.next();
                List<ChroPeptide> peptideList = new ArrayList<>();
                for(Iterator peptideItr = experiementPeptide.iterator();peptideItr.hasNext();)
                {
                    ChroPeptide peptide = new ChroPeptide();
                    peptide.readJsonObject((JSONObject) peptideItr.next());
                    peptideList.add(peptide);
                }
                this.peptideGroup.add(peptideList);
            }
        }
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
}
