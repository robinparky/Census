package edu.scripps.pms.census.labelFree;

import edu.scripps.pms.census.model.ChroPeptide;
import edu.scripps.pms.census.model.ChroProtein;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;


/**
 * txttmp reader.
 * Reads the file return  the model with Protein object and list of peptide, where the key is 
 * seq_chargeState.
 * @author harshil
 */
public class LabelFreeReader 
{
    
    private BufferedReader br = null;
    private String txttmpFile;
            
    private final int eachExperimentLength = 19;
    private int totalExperiments = 0;

    private List<Integer> sequenceIndexList = new ArrayList<Integer>();
    private List<Integer> fileNameIndexList = new ArrayList<Integer>();
    private List<Integer> scanIndexList = new ArrayList<Integer>();
    private List<Integer> csIndexList = new ArrayList<Integer>();
    private List<Integer> intensityIndexList = new ArrayList<Integer>();
    private List<Integer> profileScoreIndexList = new ArrayList<Integer>();
    private List<Integer> mhPlusIndexList = new ArrayList<Integer>();
    private List<Integer> calcMHPlusIndexList = new ArrayList<Integer>();
    private List<Integer> totalIntensityIndexList = new ArrayList<Integer>();
    private List<Integer> xCorIndexList = new ArrayList<Integer>();
    private List<Integer> dcnIndexList = new ArrayList<Integer>();
    private List<Integer> dMassIndexList = new ArrayList<Integer>();
    private List<Integer> sprankIndexList = new ArrayList<Integer>();
    private List<Integer> spScoreIndexList = new ArrayList<Integer>();
    private List<Integer> redundancyIndexList = new ArrayList<Integer>();
    private List<Integer> startIndexList = new ArrayList<Integer>();
    private List<Integer> endIndexList = new ArrayList<Integer>();
    private List<Integer> retentionIndexList = new ArrayList<Integer>();
    private List<Integer> ionInjectionIndexList = new ArrayList<Integer>();
    private List<Integer> scountIndexList = new ArrayList<Integer>();
    private List<Integer> normIntensityIndexList = new ArrayList<Integer>();
    private int accessionIndex = -1;
    private int descriptionIndex  =-1;
    private int pepCountIndex  =-1;
    
    
    public LabelFreeReader(String txttmpFile) {
        this.txttmpFile = txttmpFile;
        try {
            br = new BufferedReader(new FileReader(txttmpFile));
        } catch (FileNotFoundException ex) {
            Logger.getLogger(LabelFreeReader.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    public static void main(String args[])
    {
        LabelFreeReader reader = new LabelFreeReader("/home/harshil/census_labelfree_out_7841.txttmp");
        System.out.println(reader.readWholeFile().size());
        
    }
    
    
    /**
     * Reads whole file.
     * @return List of Protein and peptide in a Model (ProteinModel)
     */
    public List<ProteinModel> readWholeFile() {

        List<ProteinModel> proteinList = new ArrayList<>();
        try
        {
             String currentLine = null;
             ChroProtein chroProtein = null;
             ProteinModel proteinModel = null;
             
             while( (currentLine = br.readLine()) !=null)
             {
                 List<ChroPeptide> chroPepList = new ArrayList<>();
                
                 if(currentLine.startsWith("PLINE"))
                     parseProteinHeader(currentLine);
                 else if(currentLine.startsWith("SLINE"))
                     parsePeptideHeader(currentLine);
                 else if(currentLine.startsWith("P\t"))
                 {
                     if(proteinModel!= null)
                         proteinList.add(proteinModel);
                    proteinModel = new ProteinModel();
                    chroProtein = parseProteinLine(currentLine);
                    proteinModel.addRedundantProtein(chroProtein);
                 }
                 
                 else if(currentLine.startsWith("S\t"))
                 {
                     
                     chroPepList = parsePeptideLine(currentLine);
                     proteinModel.addPeptideList(chroPepList);
                 }
             }
             proteinList.add(proteinModel);
            
        } catch (Exception ex) {
            Logger.getLogger(ChroJSONGenerator.class.getName()).log(Level.SEVERE, null, ex);
        }
        finally{
            try {
                br.close();
            } catch (IOException ex) {
                Logger.getLogger(ChroJSONGenerator.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        return proteinList;
    }
    
    
    
    
    
    public ChroProtein parseProteinLine(String currentLine)
    {
        ChroProtein protein = new ChroProtein();
        String words[] = currentLine.split("\t");
        if(accessionIndex != -1)
            protein.setLocus(words[accessionIndex]);
        if(descriptionIndex != -1)
            protein.setDescription(words[descriptionIndex]);
        if(pepCountIndex != -1)
            protein.setPepCount(Integer.parseInt(words[pepCountIndex]));
        return protein;
    }
    
    public void parseProteinHeader(String currentLine)
    {
      //  ACCESSION	DESCRIPTION	SCOUNT_1	SCOUNT_2	SCOUNT_3	SCOUNT_4	PEP_COUNT	NORM_INTENSITY_1	NORM_INTENSITY_2	NORM_INTENSITY_3	NORM_INTENSITY_4	
        String words[] = currentLine.split("\t");
        
        for (int i = 0; i < words.length; i++) 
        {
            if (words[i].equalsIgnoreCase("ACCESSION")) 
                accessionIndex= i;
            else if (words[i].equalsIgnoreCase("DESCRIPTION")) 
                descriptionIndex= i;
            else if (words[i].contains("SCOUNT_")) 
                scountIndexList.add(i);
            else if (words[i].equalsIgnoreCase("PEP_COUNT")) 
                pepCountIndex = i;
            
        }
    }
    
    public List<ChroPeptide> parsePeptideLine(String currentLine) {
        String words[] = currentLine.split("\t");
        int wordCounter = 0;
        
        List<ChroPeptide> cPeptideList = new ArrayList<>();
        for (int i = 0; i < totalExperiments; i++) {
            
            ChroPeptide currentPeptide = new ChroPeptide();
            if (sequenceIndexList.size() > 0) {
                if (words[sequenceIndexList.get(i)].equalsIgnoreCase("NA")) {
                    cPeptideList.add(currentPeptide);
                    continue;
                } else {
                    currentPeptide.setSequence(words[sequenceIndexList.get(i)]);
                }
            }
            if (fileNameIndexList.size() > 0) {
                currentPeptide.setFileName(words[fileNameIndexList.get(i)]);
//                if ( !words[fileNameIndexList.get(i)].equalsIgnoreCase("na")) {
//                    msFileList.set(i,words[fileNameIndexList.get(i)]);
//                }
            }
            if (scanIndexList.size() > 0) {
                currentPeptide.setScanNum(Integer.parseInt(words[scanIndexList.get(i)]));
            }
            if (csIndexList.size() > 0) {
                currentPeptide.setChargeState(words[csIndexList.get(i)]);
            }
            if (intensityIndexList.size() > 0) {
                currentPeptide.setAverageIntensity(Double.parseDouble(words[intensityIndexList.get(i)]));
            }
            if (profileScoreIndexList.size() > 0) {
                currentPeptide.setAnCompositeScore(Double.parseDouble(words[profileScoreIndexList.get(i)]));
            }
            if (mhPlusIndexList.size() > 0) {
                currentPeptide.setMhPlus(words[mhPlusIndexList.get(i)]);
            }
            if (calcMHPlusIndexList.size() > 0) {
                currentPeptide.setCalcMHplus(words[calcMHPlusIndexList.get(i)]);
            }
            if (totalIntensityIndexList.size() > 0) {
                currentPeptide.setTotalIntensity(Double.parseDouble(words[totalIntensityIndexList.get(i)]));
            }
            if (xCorIndexList.size() > 0) {
                currentPeptide.setXCorr(words[xCorIndexList.get(i)]);
            }
            if (dcnIndexList.size() > 0) {
                currentPeptide.setDeltCN(words[dcnIndexList.get(i)]);
            }
            if (dMassIndexList.size() > 0) {
                currentPeptide.setDeltMass(words[dMassIndexList.get(i)]);
            }
            if (sprankIndexList.size() > 0) {
                currentPeptide.setSpRank(words[sprankIndexList.get(i)]);
            }
            if (spScoreIndexList.size() > 0) {
                currentPeptide.setSpScore(words[spScoreIndexList.get(i)]);
            }
            if (redundancyIndexList.size() > 0) {
                currentPeptide.setRedundancy(words[redundancyIndexList.get(i)]);
            }
            if (startIndexList.size() > 0) {
                currentPeptide.setStartRange(words[startIndexList.get(i)]);
            }
            if (endIndexList.size() > 0) {
                currentPeptide.setEndRange(words[endIndexList.get(i)]);
            }
            if (retentionIndexList.size() > 0) {
                currentPeptide.setRetentionTime(Double.parseDouble(words[retentionIndexList.get(i)]));
            }
            if (ionInjectionIndexList.size() > 0) {
                currentPeptide.setIonInjectionTime(Double.parseDouble(words[ionInjectionIndexList.get(i)]));
            }

            cPeptideList.add(currentPeptide);
        }

        return cPeptideList;
    }

    public void parsePeptideHeader(String currentLine) {
        String words[] = currentLine.split("\t");
        
        for (int i = 0; i < words.length; i++) {
//         
            
            if (words[i].equalsIgnoreCase("SEQUENCE")) {
                sequenceIndexList.add(i);
                totalExperiments++;
            } else if (words[i].equalsIgnoreCase("FILENAME")) {
                fileNameIndexList.add(i);
            } else if (words[i].equalsIgnoreCase("SCAN")) {
                scanIndexList.add(i);
            } else if (words[i].equalsIgnoreCase("CSTATE")) {
                csIndexList.add(i);
            } else if (words[i].equalsIgnoreCase("INTENSITY")) {
                intensityIndexList.add(i);
            } else if (words[i].equalsIgnoreCase("PROFILE_SCORE")) {
                profileScoreIndexList.add(i);
            } else if (words[i].equalsIgnoreCase("MHPLUS")) {
                mhPlusIndexList.add(i);
            } else if (words[i].equalsIgnoreCase("CALCMHPLUS")) {
                calcMHPlusIndexList.add(i);
            } else if (words[i].equalsIgnoreCase("TOTALINTENSITY")) {
                totalIntensityIndexList.add(i);
            } else if (words[i].equalsIgnoreCase("XCORR")) {
                xCorIndexList.add(i);
            } else if (words[i].equalsIgnoreCase("DCN")) {
                dcnIndexList.add(i);
            } else if (words[i].equalsIgnoreCase("DMASS")) {
                dMassIndexList.add(i);
            } else if (words[i].equalsIgnoreCase("SPRANK")) {
                sprankIndexList.add(i);
            } else if (words[i].equalsIgnoreCase("SPSCORE")) {
                spScoreIndexList.add(i);
            } else if (words[i].equalsIgnoreCase("REDUNDANCY")) {
                redundancyIndexList.add(i);
            } else if (words[i].equalsIgnoreCase("STARTRANGE")) {
                startIndexList.add(i);
            } else if (words[i].equalsIgnoreCase("ENDRANGE")) {
                endIndexList.add(i);
            } else if (words[i].equalsIgnoreCase("RETENTIONTIME")) {
                retentionIndexList.add(i);
            } else if (words[i].equalsIgnoreCase("IONINJECTIONTIME")) {
                ionInjectionIndexList.add(i);
            }
        }
    }

    /**
     * If we found on peptide_ChargeState in one protein but its not there in 
     * another protein then this method will add that missing peptide into that peptide too.
     * @return 
     */
    private void addMissingPeptide(List<ProteinModel> proteinList){
        
        for(ProteinModel proteinModel : proteinList)
        {
            for(List<ChroPeptide> currentPeptideGroup : proteinModel.getPeptideMap().values())
            {
//                for(ChroPeptide peptide : )
            }
        }
        
    }
}