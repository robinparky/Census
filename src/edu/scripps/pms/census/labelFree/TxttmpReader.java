package edu.scripps.pms.census.labelFree;

import com.mongodb.util.Hash;
import edu.scripps.pms.census.labelFree.json.LabelFreeJSONPeptide;
import edu.scripps.pms.census.labelFree.json.LabelFreeJSONProtein;
import edu.scripps.pms.census.labelFree.model.LabelfreePeptide;
import edu.scripps.pms.census.model.ChroPeptide;
import edu.scripps.pms.census.model.ChroProtein;

import java.io.*;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;


/**
 * txttmp reader.
 * Reads the file return  the model with Protein object and list of peptide, where the key is 
 * seq_chargeState.
 * @author harshil
 */
public class TxttmpReader 
{
    
    private BufferedReader br = null;
    private String txttmpFile;
            
    protected final int eachExperimentLength = 19;
    protected int totalExperiments = 0;

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
    private List<Integer> corrInjectionIntensityIndexList = new ArrayList<Integer>();
    private List<Integer> scountIndexList = new ArrayList<Integer>();
    private List<Integer> normIntensityIndexList = new ArrayList<Integer>();
    private int accessionIndex = -1;
    private int descriptionIndex  =-1;
    private int pepCountIndex  =-1;
    private int[] AvgNormIntensityIndex;
    
    //PLINE   ACCESSION       DESCRIPTION     AvgNormIntensity_1      AvgNormIntensity_2     
    
    public TxttmpReader(String txttmpFile) {
        this.txttmpFile = txttmpFile;
        try {
            br = new BufferedReader(new FileReader(txttmpFile));
        } catch (FileNotFoundException ex) {
            Logger.getLogger(TxttmpReader.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public TxttmpReader() {
    }
    
    
    public static void main(String args[])
    {
        //TxttmpReader reader = new TxttmpReader("/data/2/rpark/ip2_data//rpark/fusion_labelfree/labelfree_quant/labelfree_10809.txt");
        //TxttmpReader reader = new TxttmpReader("/data/2/rpark/ip2_data//rpark/fusion_labelfree/labelfree_quant/labelfree_10986/census_labelfree_result_10986.txt");
        
//        reader.regenerateFile();
      //  System.out.println(reader.readWholeFile().size());
        TxttmpReader txtReader = new TxttmpReader("/data/2/rpark/ip2_data/hyukim/nomogram/labelfree_quant/labelfree_18910/small_test/small.txttmp");
        List<ProteinModel> proteinList = txtReader.readWholeFile();
        System.out.println("");
        
    }
    
    
    
    /**
     * Reads whole file.
     * This read redundnat prtoeins on separate entry with separate peptide list.  
     * Robin will retire this function
     * @return List of Protein and peptide in a Model (ProteinModel)
     */
    /*
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
                    proteinModel.setChroProtein(chroProtein);
                 }
                 
                 else if(currentLine.startsWith("S\t"))
                 {
                     //System.out.println("==" + curre);
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
    */
    public List<ProteinModel> readWholeFile() {

        List<ProteinModel> proteinList = new ArrayList<>();
        try
        {
            String eachLine = null;
            ChroProtein chroProtein = null;
            ProteinModel proteinModel = null;
            List<ChroPeptide> chroPepList;

          
            boolean isNewProtein=true;

                 
            while( (eachLine = br.readLine()) != null)
            {
                if(eachLine.startsWith("H\t")) continue;
                
                if(eachLine.startsWith("PLINE"))
                     parseProteinHeader(eachLine);
                 else if(eachLine.startsWith("SLINE"))
                     parsePeptideHeader(eachLine);
                 else if(eachLine.startsWith("P\t"))                
                {
                  //  arr = eachLine.split("\t");

                    if(!isNewProtein)
                    {
                   //     proteinModel = new ProteinModel();
                        chroProtein = parseProteinLine(eachLine);
                        proteinModel.addRedundantProtein(chroProtein);
                    }
                    else
                    {
                        proteinModel = new ProteinModel();                        
                        chroProtein = parseProteinLine(eachLine);
                        proteinModel.addRedundantProtein(chroProtein);

                        proteinList.add(proteinModel);
                        isNewProtein = false;
                    }

                } else {
                    isNewProtein = true;

                    chroPepList = chroPepList = parsePeptideLine(eachLine);
                    proteinModel.addPeptideList(chroPepList);                     
                   
                }
            }


            int numOfExperiments = proteinList.get(0).getPeptideList().get(0).getPeptideList().size();
            List<Hashtable<String, ChroPeptide>> foundPeptideHtList = new ArrayList<>();
            for(int i=0;i<numOfExperiments;i++)
                foundPeptideHtList.add(new Hashtable<>());


            for (Iterator<ProteinModel> itr = proteinList.iterator(); itr.hasNext();) {
                ProteinModel p = itr.next();


                for (LabelfreePeptide each : p.getPeptideList()) {

                    int count=0;
                    for (ChroPeptide expPep : each.getPeptideList()) {
                        if(null == expPep.getFileName()) {
                            count++;
                            continue;
                        }

                        Hashtable<String, ChroPeptide> tmpHt = foundPeptideHtList.get(count);
                        String key = expPep.getSequence() + expPep.getChargeState();
                        if(null == tmpHt.get(key))
                            tmpHt.put(key, expPep);

                        //System.out.println(expPep.getFileName() + " " + count + " " + expPep.getSequence());

                        count++;
                    }
                }
            }


            for (Iterator<ProteinModel> itr = proteinList.iterator(); itr.hasNext();) {
                ProteinModel p = itr.next();


                for (LabelfreePeptide each : p.getPeptideList()) {

                    String pepKey = null;
                    for (ChroPeptide expPep : each.getPeptideList()) {
                        if (null != expPep.getSequence()) {
                            pepKey = expPep.getSequence() + expPep.getChargeState();
                            break;
                        }
                    }

                  //  System.out.println(pepKey);


                    int count=0;
                    for (ChroPeptide expPep : each.getPeptideList()) {
                        if(null != expPep.getFileName()) {
                            count++;
                            continue;
                        }

                        Hashtable<String, ChroPeptide> tmpHt = foundPeptideHtList.get(count);
                        ChroPeptide chroPeptide = tmpHt.get(pepKey);

                        if(null != chroPeptide)
                            each.getPeptideList().set(count, chroPeptide);

                        count++;
                    }
                }
            }



            //proteinList.get(0).getPeptideList().get(0).
            
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
        if(scountIndexList.size()!= 0)
        {
            for (int index : scountIndexList)
            {
                if(words[index].equalsIgnoreCase("NA"))
                    protein.addSpecCountList(-1);
                else
                    protein.addSpecCountList(Integer.parseInt(words[index]));
            }
        }
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
        List<ChroPeptide> cPeptideList = new ArrayList<>();
        String words[] = currentLine.split("\t");
       // int wordCounter = 0;
        
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
            if (corrInjectionIntensityIndexList.size() > 0) {
                currentPeptide.setCorrIonInjectionIntensity(Double.parseDouble(words[corrInjectionIntensityIndexList.get(i)]));
            }

            cPeptideList.add(currentPeptide);
        }
        

        return cPeptideList;
    }

    public void parsePeptideHeader(String currentLine) {
        String words[] = currentLine.split("\t");
        
        for (int i = 0; i < words.length; i++) {
//         
            
            if (words[i].equalsIgnoreCase(ChroJSONReader.SEQUENCE)) {
                sequenceIndexList.add(i);
                totalExperiments++;
            } else if (words[i].equalsIgnoreCase(ChroJSONReader.FILENAME)) {
                fileNameIndexList.add(i);
            } else if (words[i].equalsIgnoreCase(ChroJSONReader.SCAN)) {
                scanIndexList.add(i);
            } else if (words[i].equalsIgnoreCase(ChroJSONReader.CSTATE)) {
                csIndexList.add(i);
            } else if (words[i].equalsIgnoreCase(ChroJSONReader.INTENSITY)) {
                intensityIndexList.add(i);
            } else if (words[i].equalsIgnoreCase(ChroJSONReader.PROFILE_SCORE)) {
                profileScoreIndexList.add(i);
            } else if (words[i].equalsIgnoreCase(ChroJSONReader.MHPLUS)) {
                mhPlusIndexList.add(i);
            } else if (words[i].equalsIgnoreCase(ChroJSONReader.CALCMHPLUS)) {
                calcMHPlusIndexList.add(i);
            } else if (words[i].equalsIgnoreCase(ChroJSONReader.TOTALINTENSITY)) {
                totalIntensityIndexList.add(i);
            } else if (words[i].equalsIgnoreCase(ChroJSONReader.XCORR)) {
                xCorIndexList.add(i);
            } else if (words[i].equalsIgnoreCase(ChroJSONReader.DCN)) {
                dcnIndexList.add(i);
            } else if (words[i].equalsIgnoreCase(ChroJSONReader.DMASS)) {
                dMassIndexList.add(i);
            } else if (words[i].equalsIgnoreCase(ChroJSONReader.SPRANK)) {
                sprankIndexList.add(i);
            } else if (words[i].equalsIgnoreCase(ChroJSONReader.SPSCORE)) {
                spScoreIndexList.add(i);
            } else if (words[i].equalsIgnoreCase(ChroJSONReader.REDUNDANCY)) {
                redundancyIndexList.add(i);
            } else if (words[i].equalsIgnoreCase(ChroJSONReader.STARTRANGE)) {
                startIndexList.add(i);
            } else if (words[i].equalsIgnoreCase(ChroJSONReader.ENDRANGE)) {
                endIndexList.add(i);
            } else if (words[i].equalsIgnoreCase(ChroJSONReader.RETENTIONTIME)) {
                retentionIndexList.add(i);
            } else if (words[i].equalsIgnoreCase(ChroJSONReader.IONINJECTIONTIME)) {
                ionInjectionIndexList.add(i);
	    }else if (words[i].equalsIgnoreCase("CORRIONINJECTION_INTENISTY") || words[i].equalsIgnoreCase("CORRIONINJECTION_INTENSITY")) {
                corrInjectionIntensityIndexList.add(i);
            }
        }
    }

    /**
     * If we found on peptide_ChargeState in one protein but its not there in 
     * another protein then this method will add that missing peptide into that peptide too.
     * @return 
     */
    private void addMissingPeptide(List<ProteinModel> proteinList){
        
        HashMap<String, List<ChroPeptide>> seqToPep = getPeptideMap(proteinList);
        //add Missing information from seqToPep map.....
        for(ProteinModel proteinModel : proteinList)
        {
            for(List<ChroPeptide> currentPeptideGroup : proteinModel.getPeptideMap().values())
            {
                for(int i=0;i<currentPeptideGroup.size();i++)
                {
                    ChroPeptide peptide = currentPeptideGroup.get(i);
                    String key = getkey(currentPeptideGroup);
                     
                    if(peptide.getSequence()== null)
                    {
                        if(seqToPep.containsKey(key))
                        {
                            peptide = seqToPep.get(key).get(i);
                            currentPeptideGroup.set(i, peptide);
                            //DEBUG purpose
//                            if(peptide.getSequence()!=null)
//                                System.out.println(peptide.getSequence() +" is added to the new File");
                                
                        }
                    }
                }
            }
        }
        
    }
    /**
     * this will return all the peptide Map from the txttmp file.
     */
    public HashMap<String, List<ChroPeptide>> getPeptideMap(List<ProteinModel> proteinList)
    {
        HashMap<String, List<ChroPeptide>> seqToPep = new HashMap<>();
        //Generates the seqToPep hashMap
        for(ProteinModel proteinModel : proteinList)
        {
            for(List<ChroPeptide> currentPeptideGroup : proteinModel.getPeptideMap().values())
            {
                for(int i=0;i<currentPeptideGroup.size();i++)
                {
                    ChroPeptide peptide = currentPeptideGroup.get(i);
                    String key = getkey(currentPeptideGroup);
                    if(!seqToPep.containsKey(key))
                    {
                        seqToPep.put(key, currentPeptideGroup);
                    }
                        
                     if(peptide.getSequence() != null)
                     {
                         seqToPep.get(key).set(i, peptide);
                     }
//                    if(peptide.getSequence()== null)
//                    {
//                        String key = getkey(currentPeptideGroup);
//                        if(seqToPep.containsKey(key))
//                    }
                }
            }
        }
        return seqToPep;
    }
    
    /**
     * 
     * @param pepList
     * @return key will always be sequence_chargeState
     */
    private String getkey(List<ChroPeptide> pepList)
    {
        for(ChroPeptide peptide : pepList)
        {
            if(peptide.getSequence()!= null)
                return  peptide.getSequence()+"_"+peptide.getChargeState();
        }
        return null;
    }
    
    /**
     * If we found on peptide_ChargeState in one protein but its not there in 
     * another protein then this method will add that missing peptide into that peptide too.
     * 
     * Robin terminated this method.  label free needs to be cleaned up with new algorithm
     */
    public void regenerateFile()
    {
        /*
        List<ProteinModel> proteinList = readWholeFile();
        addMissingPeptide(proteinList);
        StringBuffer output = new StringBuffer();
        output.append(getHeader());
        for(ProteinModel proteinModel : proteinList)
        {
           output.append(toProteinSting(proteinModel.getChroProtein()));
           for(List<ChroPeptide> pepList : proteinModel.getPeptideMap().values())
           {
               StringBuffer temp = new StringBuffer();
               temp.append("S");
               for(int i =0;i<pepList.size(); i++)
               {
                   ChroPeptide peptide = pepList.get(i);
                   temp.append(toPeptideString(peptide, i));
               }
               output.append(temp).append("\n");
           }
        }
        BufferedWriter bw = null;
        try {
            bw = new BufferedWriter(new FileWriter(txttmpFile));
            bw.write(output.toString());
            System.out.println("Txttmp file is regenerated " + txttmpFile);
        } catch (IOException ex) {
            Logger.getLogger(TxttmpReader.class.getName()).log(Level.SEVERE, null, ex);
        } finally{
            try {
                bw.close();
            } catch (IOException ex) {
                Logger.getLogger(TxttmpReader.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        */
        
    }
    
    private String toPeptideString(ChroPeptide chroPeptide,int index)
  {
      StringBuffer sb = new StringBuffer();
      sb.append("\t[").append(index+1).append("]");
      if(chroPeptide.getSequence() == null)
      {
          for(int i =0;i<eachExperimentLength;i++)
              sb.append("\tNA");
          return sb.toString();
      }
      sb.append("\t").append(chroPeptide.getSequence());
      sb.append("\t").append(chroPeptide.getFileName());
      sb.append("\t").append(chroPeptide.getScanNum());
      sb.append("\t").append(chroPeptide.getChargeState());
      
      sb.append("\t").append(chroPeptide.getAverageIntensity());
      sb.append("\t").append(chroPeptide.getAnCompositeScore());
      
      
      sb.append("\t").append(chroPeptide.getMhPlus());
      sb.append("\t").append(chroPeptide.getCalcMHplus());
      sb.append("\t").append(chroPeptide.getTotalIntensity());
      sb.append("\t").append(chroPeptide.getXCorr());
      sb.append("\t").append(chroPeptide.getDeltCN());
      sb.append("\t").append(chroPeptide.getDeltMass());
      
      sb.append("\t").append(chroPeptide.getSpRank());
      sb.append("\t").append(chroPeptide.getSpScore());
      sb.append("\t").append(chroPeptide.getRedundancy());
      sb.append("\t").append(chroPeptide.getStartRange());
      sb.append("\t").append(chroPeptide.getEndRange());
      sb.append("\t").append(chroPeptide.getRetentionTime());
      sb.append("\t").append(chroPeptide.getIonInjectionTime());
      
      return sb.toString();
  }
    
    private String toProteinSting(ChroProtein chroProtein)
    {
    //    PLINE   ACCESSION       DESCRIPTION     SCOUNT_1        SCOUNT_2        SCOUNT_3        SCOUNT_4        PEP_COUNT 
        StringBuffer sb = new StringBuffer();
        sb.append("P").append("\t");
        sb.append(chroProtein.getLocus()).append("\t");
        sb.append(chroProtein.getDescription()).append("\t");
        for(int value : chroProtein.getSpecCountList())
            sb.append(value).append("\t");
        sb.append(chroProtein.getPepCount()).append("\n");
        
        
        
        return sb.toString();

    }
    private  String getHeader()
    {
        BufferedReader br = null;
        StringBuffer header = new StringBuffer();
        String line = null;
        try {
            br = new BufferedReader(new FileReader(txttmpFile));
            while((line = br.readLine()) != null )
            {
                if(line.startsWith("P\t"))
                    break;
                header.append(line).append("\n");
            }
        } catch (FileNotFoundException ex) {
            Logger.getLogger(TxttmpReader.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(TxttmpReader.class.getName()).log(Level.SEVERE, null, ex);
        }
        finally{
            try {
                br.close();
            } catch (IOException ex) {
                Logger.getLogger(TxttmpReader.class.getName()).log(Level.SEVERE, null, ex);
            }
            
        }
        return header.toString();
    }
    
    
   //getter and setter
    public BufferedReader getBr() {
        return br;
    }

    public void setBr(BufferedReader br) {
        this.br = br;
    }

    public String getTxttmpFile() {
        return txttmpFile;
    }

    public void setTxttmpFile(String txttmpFile) {
        this.txttmpFile = txttmpFile;
    }

    public int getTotalExperiments() {
        return totalExperiments;
    }

    public void setTotalExperiments(int totalExperiments) {
        this.totalExperiments = totalExperiments;
    }

    public List<Integer> getSequenceIndexList() {
        return sequenceIndexList;
    }

    public void setSequenceIndexList(List<Integer> sequenceIndexList) {
        this.sequenceIndexList = sequenceIndexList;
    }

    public List<Integer> getFileNameIndexList() {
        return fileNameIndexList;
    }

    public void setFileNameIndexList(List<Integer> fileNameIndexList) {
        this.fileNameIndexList = fileNameIndexList;
    }

    public List<Integer> getScanIndexList() {
        return scanIndexList;
    }

    public void setScanIndexList(List<Integer> scanIndexList) {
        this.scanIndexList = scanIndexList;
    }

    public List<Integer> getCsIndexList() {
        return csIndexList;
    }

    public void setCsIndexList(List<Integer> csIndexList) {
        this.csIndexList = csIndexList;
    }

    public List<Integer> getIntensityIndexList() {
        return intensityIndexList;
    }

    public void setIntensityIndexList(List<Integer> intensityIndexList) {
        this.intensityIndexList = intensityIndexList;
    }

    public List<Integer> getProfileScoreIndexList() {
        return profileScoreIndexList;
    }

    public void setProfileScoreIndexList(List<Integer> profileScoreIndexList) {
        this.profileScoreIndexList = profileScoreIndexList;
    }

    public List<Integer> getMhPlusIndexList() {
        return mhPlusIndexList;
    }

    public void setMhPlusIndexList(List<Integer> mhPlusIndexList) {
        this.mhPlusIndexList = mhPlusIndexList;
    }

    public List<Integer> getCalcMHPlusIndexList() {
        return calcMHPlusIndexList;
    }

    public void setCalcMHPlusIndexList(List<Integer> calcMHPlusIndexList) {
        this.calcMHPlusIndexList = calcMHPlusIndexList;
    }

    public List<Integer> getTotalIntensityIndexList() {
        return totalIntensityIndexList;
    }

    public void setTotalIntensityIndexList(List<Integer> totalIntensityIndexList) {
        this.totalIntensityIndexList = totalIntensityIndexList;
    }

    public List<Integer> getxCorIndexList() {
        return xCorIndexList;
    }

    public void setxCorIndexList(List<Integer> xCorIndexList) {
        this.xCorIndexList = xCorIndexList;
    }

    public List<Integer> getDcnIndexList() {
        return dcnIndexList;
    }

    public void setDcnIndexList(List<Integer> dcnIndexList) {
        this.dcnIndexList = dcnIndexList;
    }

    public List<Integer> getdMassIndexList() {
        return dMassIndexList;
    }

    public void setdMassIndexList(List<Integer> dMassIndexList) {
        this.dMassIndexList = dMassIndexList;
    }

    public List<Integer> getSprankIndexList() {
        return sprankIndexList;
    }

    public void setSprankIndexList(List<Integer> sprankIndexList) {
        this.sprankIndexList = sprankIndexList;
    }

    public List<Integer> getSpScoreIndexList() {
        return spScoreIndexList;
    }

    public void setSpScoreIndexList(List<Integer> spScoreIndexList) {
        this.spScoreIndexList = spScoreIndexList;
    }

    public List<Integer> getRedundancyIndexList() {
        return redundancyIndexList;
    }

    public void setRedundancyIndexList(List<Integer> redundancyIndexList) {
        this.redundancyIndexList = redundancyIndexList;
    }

    public List<Integer> getStartIndexList() {
        return startIndexList;
    }

    public void setStartIndexList(List<Integer> startIndexList) {
        this.startIndexList = startIndexList;
    }

    public List<Integer> getEndIndexList() {
        return endIndexList;
    }

    public void setEndIndexList(List<Integer> endIndexList) {
        this.endIndexList = endIndexList;
    }

    public List<Integer> getRetentionIndexList() {
        return retentionIndexList;
    }

    public void setRetentionIndexList(List<Integer> retentionIndexList) {
        this.retentionIndexList = retentionIndexList;
    }

    public List<Integer> getIonInjectionIndexList() {
        return ionInjectionIndexList;
    }

    public void setIonInjectionIndexList(List<Integer> ionInjectionIndexList) {
        this.ionInjectionIndexList = ionInjectionIndexList;
    }

    public List<Integer> getScountIndexList() {
        return scountIndexList;
    }

    public void setScountIndexList(List<Integer> scountIndexList) {
        this.scountIndexList = scountIndexList;
    }

    public List<Integer> getNormIntensityIndexList() {
        return normIntensityIndexList;
    }

    public void setNormIntensityIndexList(List<Integer> normIntensityIndexList) {
        this.normIntensityIndexList = normIntensityIndexList;
    }

    public int getAccessionIndex() {
        return accessionIndex;
    }

    public void setAccessionIndex(int accessionIndex) {
        this.accessionIndex = accessionIndex;
    }

    public int getDescriptionIndex() {
        return descriptionIndex;
    }

    public void setDescriptionIndex(int descriptionIndex) {
        this.descriptionIndex = descriptionIndex;
    }

    public int getPepCountIndex() {
        return pepCountIndex;
    }

    public void setPepCountIndex(int pepCountIndex) {
        this.pepCountIndex = pepCountIndex;
    }

    public int[] getAvgNormIntensityIndex() {
        return AvgNormIntensityIndex;
    }

    public void setAvgNormIntensityIndex(int[] AvgNormIntensityIndex) {
        this.AvgNormIntensityIndex = AvgNormIntensityIndex;
    }
}
