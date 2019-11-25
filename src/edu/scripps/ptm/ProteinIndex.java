/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package edu.scripps.ptm;

import java.util.*;

/**
 *
 * @author Harshil
 */
public class ProteinIndex {
    private final int proteinScaledLength=400;
    private String accession;
    private int[] unModifiedArr;
    private int[] modifiedArr;
    private Set<Integer> changedIndex = new LinkedHashSet<Integer>();
    private String proteinSequence;
    private String modSiteString;//format is S:27,T:32,
    private int proteinLength;
    
    private HashMap<Integer, List<Double>> indexToLocalization = new HashMap<>();
    private List<ResidueLocation> residueList = null;
    private List<String> ptmPeptide = new ArrayList<>();
    
    public ProteinIndex(String accession,String proteinSequence) {
      
        this.accession = accession;
        this.proteinSequence = proteinSequence;
        proteinLength = proteinSequence.length();
        unModifiedArr = new int[proteinSequence.length()];
        modifiedArr = new int[proteinSequence.length()];
    }
    
    
    public void setUnModifiedArr(int index,int value) {
        unModifiedArr[index] = value;
    }

    public void setModifiedArr(int index,int value) {
        modifiedArr[index] = value;
    }
    
    public List<ResidueLocation> getResidueList() {
        
        if(null != residueList)
            return residueList;
        
        residueList = new ArrayList<>();
        List<Integer> l = new ArrayList(changedIndex);
        for(Iterator<Integer> itr=l.iterator(); itr.hasNext(); ) {
            int i = itr.next();
            ResidueLocation r = new ResidueLocation();
            r.setResidue(proteinSequence.charAt(i));
            int relativeLocation = (int)(proteinScaledLength*i/proteinLength);
            r.setRelativeLocation(relativeLocation);
            residueList.add(r);            
        }
        
        
        return residueList;
    }

    public void setResidueList(List<ResidueLocation> residueList) {
        this.residueList = residueList;
    }

    
    public int getProteinLength() {
        return this.proteinSequence.length();        
    }

    public void setProteinLength(int proteinLength) {
        this.proteinLength = proteinLength;
    }
    
    

    public String getModSiteString() {
        
        //build string like S34 (34/43), T32 (53,32)
//        return "dummy site index";
        return modSiteString;
    }

    public void setModSiteString(String modSiteString) {
        this.modSiteString = modSiteString;
    }
    
    
    
    private void updateModSiteString()
    {
       String oldData[] = this.modSiteString.split(";");
       StringBuffer sb = new StringBuffer();
       Iterator<Integer> changedIndexItr = changedIndex.iterator();
       for(int i =0; i<oldData.length;i++)
       {
           String eachString = oldData[i];
           sb.append(eachString.substring(0, eachString.length()-1));
           List<Double> localizationScore = indexToLocalization.get(changedIndexItr.next());
           double averageLocalizationScore =0.0;
           if(localizationScore != null)
           {
                for(Double value : localizationScore)
                     averageLocalizationScore+=value;   
                averageLocalizationScore/=localizationScore.size();
           }
//           else
//                System.err.println("Not found..................");

           sb.append(",").append(averageLocalizationScore);
           sb.append(");");
       }
       this.modSiteString=sb.toString();
    }
    
    

    public void addChangedIndex(int value)
    {
        getChangedIndex().add(value);
        addPtmPeptide(value);
        
    }
    
    
    public void addModifiedIndex(int index, int value)
    {
        modifiedArr[index] += value;
    }
    
    public void addUnModifiedIndex(int index, int value)
    {
        unModifiedArr[index] += value;
    }
    public void setUnModifiedIndex(int index, int value)
    {
        unModifiedArr[index]=value;
    }
    public String getAccession() {
        return accession;
    }

    public void setAccession(String accession) {
        this.accession = accession;
    }

    /**
     * @return the unModifiedArr
     */
    public int[] getUnModifiedArr() {
        return unModifiedArr;
    }

    /**
     * @return the modifiedArr
     */
    public int[] getModifiedArr() {
        return modifiedArr;
    }

    /**
     * @return the changedIndex
     */
    public Set<Integer> getChangedIndex() {
        return changedIndex;
    }

    /**
     * @return the proteinSequence
     */
    public String getProteinSequence() {
        return proteinSequence;
    }

    /**
     * @return the ptmPeptide
     */
    public List<String> getPtmPeptide() {
        return ptmPeptide;
    }

    /**
     * @param ptmPeptide the ptmPeptide to set
     */
    public void setPtmPeptide(List<String> ptmPeptide) {
        this.ptmPeptide = ptmPeptide;
    }

   
    public void addPtmPeptide(int index) {
        String s = getProteinSequence();
        int start =(index-7 <0) ? 0:index-7;
        int end =(index+7 > s.length()) ? s.length()-1:index+7;
        
        StringBuffer sb = new StringBuffer();
        sb.append(s.subSequence(start,index-1)).append(".");
        sb.append(s.charAt(index)).append(".");
        sb.append(s.subSequence(index+1,end));
        this.ptmPeptide.add(sb.toString());
    }

    /**
     * @return the indexToLocalization
     */
    public HashMap<Integer, List<Double>> getIndexToLocalization() {
        return indexToLocalization;
    }

    /**
     * @param indexToLocalization the indexToLocalization to set
     */
    public void setIndexToLocalization(HashMap<Integer, List<Double>> indexToLocalization) {
        this.indexToLocalization = indexToLocalization;
        updateModSiteString();

    }
    
    
}
