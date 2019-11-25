/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package edu.scripps.ptm;

import edu.scripps.pms.util.seq.Fasta;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * This reader will read the  DTASelect-filter.txt.phospho file.....
 * @author Harshil
 */
public class DTASelectPhosphoReader {
    
    
    private int uniqueIndex=-1;
    private int fileNameIndex=-1 ;
    private int xCorrIndex=-1;
    private int deltaCNIndex=-1;
    private int confIndex=-1;
    private int mhIndex=-1;
    private int calcMhIndex=-1;
    private int ppmIndex=-1;
    private int totalIntensityIndex=-1;
    private int sprIndex=-1;
    private int probScoreIndex=-1;
    private int piIndex=-1;
    private int ionProportionIndex=-1;
    private int redundancyIndex=-1;
    private int sequenceIndex=-1;
    private int proteinAccessionIndex=-1;
    private int proteinDescriptionIndex=-1;
    private int modSequenceIndex=-1;
    private int localizationScoreIndex=-1;
    private int debunkerScoreIndex=-1;
    
    
    private BufferedReader br = null;
    private String headerLine = null;

    public DTASelectPhosphoReader(String fileName) {
        try {
            br = new BufferedReader(new FileReader(fileName));
            init();
        } catch (FileNotFoundException ex) {
            Logger.getLogger(DTASelectPhosphoReader.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    
    private void init(){
        try {
            this.headerLine = br.readLine();
        } catch (IOException ex) {
            Logger.getLogger(DTASelectPhosphoReader.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        
        String headerWords[] = this.headerLine.split("\t");
        for(int i =0;i<headerWords.length;i++)
        {
            String word = headerWords[i];
            if(word.equalsIgnoreCase("Unique") && uniqueIndex == -1)
                uniqueIndex= i;
            else if(fileNameIndex == -1 && word.equalsIgnoreCase("FileName")  )
                fileNameIndex= i;
            else if(xCorrIndex == -1 && word.equalsIgnoreCase("XCorr")  )
                xCorrIndex= i;
            else if(deltaCNIndex == -1 && word.equalsIgnoreCase("DeltCN")  )
                deltaCNIndex= i;
            else if(confIndex == -1 && word.equalsIgnoreCase("Conf%")  )
                confIndex= i;
            else if(mhIndex == -1 && word.equalsIgnoreCase("M+H+")  )
                mhIndex= i;
            else if(calcMhIndex == -1 && word.equalsIgnoreCase("CalcM+H+")  )
                calcMhIndex= i;
            else if(ppmIndex == -1 && word.equalsIgnoreCase("PPM")  )
                ppmIndex= i;
            else if(totalIntensityIndex == -1 && word.equalsIgnoreCase("TotalIntensity")  )
                totalIntensityIndex= i;
            else if(sprIndex == -1 && word.equalsIgnoreCase("SpR")  )
                sprIndex= i;
            else if(probScoreIndex == -1 && word.equalsIgnoreCase("Prob Score")  )
                probScoreIndex= i;
            else if(piIndex == -1 && word.equalsIgnoreCase("pI")  )
                piIndex= i;
            else if(ionProportionIndex == -1 && word.equalsIgnoreCase("IonProportion")  )
                ionProportionIndex= i;
            else if(redundancyIndex == -1 && word.equalsIgnoreCase("Redundancy")  )
                redundancyIndex= i;
            else if(sequenceIndex == -1 && word.equalsIgnoreCase("Sequence")  )
                sequenceIndex= i;
            else if(proteinAccessionIndex == -1 && word.equalsIgnoreCase("Proteins")  )
                proteinAccessionIndex= i;
            else if(proteinDescriptionIndex == -1 && word.equalsIgnoreCase("Protein Descriptions")  )
                proteinDescriptionIndex= i;
            else if(modSequenceIndex == -1 && word.equalsIgnoreCase("Mod Sequence")  )
                modSequenceIndex= i;
            else if(localizationScoreIndex == -1 && word.equalsIgnoreCase("Localization Score")  )
                localizationScoreIndex= i;
            else if(debunkerScoreIndex == -1 && word.equalsIgnoreCase("Debunker Score")  )
                debunkerScoreIndex= i;
								
                
        }
        
    }
    
    public Iterator<PeptidePhosphoModel> getPeptides()
    {
        return new Iterator<PeptidePhosphoModel>() {

            PeptidePhosphoModel peptide = null;
            @Override
            public boolean hasNext() {
                return getProtein();
            }

            @Override
            public PeptidePhosphoModel next() {
                return peptide;
            }

            @Override
            public void remove() {
                throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
            }
            
            private boolean getProtein()
            {
                String line=null;
                try {
                    line = br.readLine();
                } catch (IOException ex) {
                    Logger.getLogger(DTASelectPhosphoReader.class.getName()).log(Level.SEVERE, null, ex);
                }
                if(line == null)
                    return false;
                
                peptide = new PeptidePhosphoModel();
                String words[] = line.split("\t");
                if(uniqueIndex!=-1)
                    peptide.setUnique(words[uniqueIndex]);
                if(fileNameIndex!=-1)
                    peptide.setFileName(words[fileNameIndex]);
                if(xCorrIndex!=-1)
                    peptide.setxCorr(Double.parseDouble(words[xCorrIndex]));
                if(deltaCNIndex!=-1)
                    peptide.setDeltaCN(Double.parseDouble(words[deltaCNIndex]));
                if(confIndex!=-1)
                    peptide.setConf(Double.parseDouble(words[confIndex]));
                if(mhIndex!=-1)
                    peptide.setMh(Double.parseDouble(words[mhIndex]));
                if(calcMhIndex!=-1)
                    peptide.setCalcMh(Double.parseDouble(words[calcMhIndex]));
                if(ppmIndex!=-1)
                    peptide.setPpm(Double.parseDouble(words[ppmIndex]));
                if(totalIntensityIndex!=-1)
                    peptide.setTotalIntensity(Double.parseDouble(words[totalIntensityIndex]));
                if(sprIndex!=-1)
                    peptide.setSpr(Integer.parseInt(words[sprIndex]));
                if(probScoreIndex!=-1)
                    peptide.setProbScore(Double.parseDouble(words[probScoreIndex]));
                if(piIndex!=-1)
                    peptide.setPi(Double.parseDouble(words[piIndex]));
                if(ionProportionIndex!=-1)
                    peptide.setIonProportion(Double.parseDouble(words[ionProportionIndex]));
                if(redundancyIndex!=-1)
                    peptide.setRedundancy(Integer.parseInt(words[redundancyIndex]));
                if(sequenceIndex!=-1)
                    peptide.setSequence(words[sequenceIndex]);
                if(proteinAccessionIndex!=-1)
                {
                    for(String value : words[proteinAccessionIndex].split(","))
                        peptide.addProteinAccession(value);
                }
                if(proteinDescriptionIndex!=-1)
                {
                    for(String value : words[proteinDescriptionIndex].split(","))
                        peptide.addProteinDescription(value);
                }
                if(modSequenceIndex!=-1)
                    peptide.setModSequence(words[modSequenceIndex]);
                if(localizationScoreIndex!=-1)
                {
                    words[localizationScoreIndex] = words[localizationScoreIndex].replaceAll("\\[", "");
                    words[localizationScoreIndex] = words[localizationScoreIndex].replaceAll("\\]", "");
                    
                    for(String word : words[localizationScoreIndex].split(","))
                    {
                        peptide.addLocalizationScore(Double.parseDouble(word));
                    }
                    
                    
                }
                if(debunkerScoreIndex!=-1)
                    peptide.setDebunkerScore(Double.parseDouble(words[debunkerScoreIndex]));
                
                return true;
            }
            
            
        };
    }
    
    
    public HashMap<String,List<PeptidePhosphoModel>> getGroupedByProtein()
    {
        HashMap<String,List<PeptidePhosphoModel>> proteinToPeptList = new HashMap<>();
        
        for(Iterator<PeptidePhosphoModel> itr = getPeptides();itr.hasNext();)
        {
            PeptidePhosphoModel peptide = itr.next();
            List<PeptidePhosphoModel> tempList = new ArrayList<>();
            for(String proteinAccession :peptide.getProteinAccession() )
            {
                if(proteinToPeptList.containsKey(proteinAccession))
                {
                    tempList = proteinToPeptList.get(proteinAccession);
                }
                    tempList.add(peptide);
                    proteinToPeptList.put(proteinAccession,tempList);
            }
        }
        return proteinToPeptList;
    }
    
    
//    public void proteinToSiteLocalizationScore()
//    {
//        HashMap<String,List<PeptidePhosphoModel>> proteinToPeptide = getGroupedByProtein();
//        HashMap<String,HashMap<Integer, List<Double>>> proteinToSiteTable = new HashMap<>();
//        for(String ProteinAccession : proteinToPeptide.keySet())
//        {
//            for(PeptidePhosphoModel peptideSequence :proteinToPeptide.get(ProteinAccession))
//            {
//                int index = getIndex(peptideSequence.getSequence(), ProteinAccession);
//                List<Double> localizationList = new ArrayList<>();
//                HashMap<Integer, List<Double>> siteTable = new HashMap<>();
//                if(proteinToSiteTable.containsKey(ProteinAccession))
//                {
//                   siteTable = proteinToSiteTable.get(ProteinAccession);
//                   
//                   if(siteTable.containsKey(index)) 
//                   {
//                       localizationList = siteTable.get(index);
//                   }
//                }
//                localizationList.add(peptideSequence.getLocalizationScoreAvg());
//                siteTable.put(index, localizationList);
//                proteinToSiteTable.put(ProteinAccession, siteTable);
//            }
//            
//        }
//        
//        System.out.println("asds");
//        
//    }
    
    public void proteinToSiteLocalizationScore(List<ProteinIndex> proteinIndexList)
    {// Calculates the Localization score for the given PRotein.......
        HashMap<String,List<PeptidePhosphoModel>> proteinToPeptide = getGroupedByProtein();
        HashMap<String,HashMap<Integer, List<Double>>> proteinToSiteTable = new HashMap<>();
        
        for(ProteinIndex proteinIndex: proteinIndexList)
        {
//            if(! proteinToPeptide.containsKey(proteinIndex.getAccession()))
//            {
////                System.out.println("--------"+proteinIndex.getAccession());
////                System.out.println("");
//            }
//            System.out.println("--------"+proteinIndex.getAccession());
            List<PeptidePhosphoModel> pepmodelList = proteinToPeptide.get(proteinIndex.getAccession());
            if(null == pepmodelList)
                pepmodelList = proteinToPeptide.get(Fasta.getAccession(proteinIndex.getAccession()));
                
            for(PeptidePhosphoModel peptideSequence :pepmodelList)
            {
                int index = getIndex(peptideSequence.getSequence(), proteinIndex.getProteinSequence());
                List<Double> localizationList = new ArrayList<>();
                HashMap<Integer, List<Double>> siteTable = new HashMap<>();
                if(proteinToSiteTable.containsKey(proteinIndex.getAccession()))
                {
                   siteTable = proteinToSiteTable.get(proteinIndex.getAccession());
                   
                   if(siteTable.containsKey(index)) 
                   {
                       localizationList = siteTable.get(index);
                   }
                }
                localizationList.add(peptideSequence.getLocalizationScoreAvg());
                siteTable.put(index, localizationList);
                proteinToSiteTable.put(proteinIndex.getAccession(), siteTable);
            }
            proteinIndex.setIndexToLocalization(proteinToSiteTable.get(proteinIndex.getAccession()));
        }
        
        
    }
    
    public static void main(String args[])
    {
        //DTASelectPhosphoReader phosphoReader = new DTASelectPhosphoReader("/data/2/rpark/ip2_data/benstein/Mammalian_LKB1_AMPK_Interactome/20131122_EVh_SandP_LKB1l_SandP_10ug_each_n1_2013_11_25_18_20659/search/projects2013_11_28_14_53260/phospho/DTASelect-filter.txt.phospho");
        DTASelectPhosphoReader phosphoReader = new DTASelectPhosphoReader("/data/2/rpark/ip2_data//benstein/STIM2_14_3_3_Target/20141025_PT_endog_S2_IP_fl_fl_wt_phen_gelband_4hr_single_phase_2014_10_26_10_28252/search/projects2014_10_26_10_72084//phospho/DTASelect-filter.txt.phospho");
        
        
        Iterator<PeptidePhosphoModel> itr = phosphoReader.getPeptides();
        for(;itr.hasNext();)
        {
            System.out.println(itr.next().getSequence());
        }
//        phosphoReader.proteinToSiteLocalizationScore();
    }
    
    
    public int getIndex(String peptideSequence,String proteinSequence )
    {
        
        String[] splitter =  peptideSequence.split("[().]");
        String searchingString="";
         for(int i=1;i<splitter.length-1;i+=3)
             searchingString += splitter[i];
        int startingIndex = proteinSequence.indexOf(searchingString);
        searchingString ="";
        splitter = peptideSequence.split("[.]");
        for(int i=1;i<splitter.length-1;i++)
             searchingString += splitter[i];
        for(int i=0;i< searchingString.length();i++)
        {
            if(i<searchingString.length()-1  && searchingString.charAt(i+1) == '(')
            {
                //Its a modifide string 
                return startingIndex;
                
//                proteinIndex.addModPeptideList(peptideSequence.getSequence());
            }
            startingIndex++;
        }
        
        return ++startingIndex;
    
    }
}
