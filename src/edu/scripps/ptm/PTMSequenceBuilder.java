package edu.scripps.ptm;

import edu.scripps.pms.census.util.dtaselect.Peptide;
import edu.scripps.pms.census.util.dtaselect.Protein;
import edu.scripps.pms.census.util.io.DTASelectFilterReader;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Harshil
 */
public class PTMSequenceBuilder 
{
    public List<ProteinIndex> getProteins(String filePath,String inputFileName,String outputFileName)
    {
        List<ProteinIndex> proteinIndexList = new ArrayList();
         try {
            DTASelectFilterReader reader = new DTASelectFilterReader(filePath+File.separator+inputFileName);
            Iterator<Protein> itr =reader.getProteins();
            Protein protein=null;
//            HashMap<String,String> proteinSequenceToAccession = new HashMap();
            while(itr.hasNext())
            {                   
                protein=itr.next();
                String accession = protein.getAccession();
                FastaProteinReader fr = new FastaProteinReader();
                String proteinSequence =fr.getProteinSeq(reader.getDbFilePathAndName(), accession);
//                if(null == proteinSequence)
//                proteinSequenceToAccession.put(proteinSequence,accession);
//                  String proteinSequence =fr.getProteinSeq(filePath + File.separator+"UniProt_Human_sprot_11-08-2010_reversed.fasta", accession);
                List<Peptide> peptideList = protein.getPeptideList();
                Iterator<Peptide> peptideIterator = peptideList.iterator();
                ProteinIndex proteinIndex = new ProteinIndex(protein.getAccession(),proteinSequence);

                while(peptideIterator.hasNext())
                {
                    Peptide peptide = peptideIterator.next();
                    getIndex(peptide, proteinSequence,proteinIndex);
                }
//                proteinIndex.generateModSiteString();
                if(proteinIndex.getChangedIndex().size()!=0)
                    proteinIndexList.add(proteinIndex);
            } 
//            DTASelectPhosphoReader phosphoReader = new DTASelectPhosphoReader(filePath + File.separator+"DTASelect-filter.txt.phospho");        
//            phosphoReader.proteinToSiteLocalizationScore(proteinIndexList);            
            
            
        } catch (IOException ex) {
            Logger.getLogger(PTMSequenceBuilder.class.getName()).log(Level.SEVERE, null, ex);
        }

         write(filePath,outputFileName,proteinIndexList);
//         read(filePath,outputFileName);
//        return proteinIndexList;
         return read(filePath,outputFileName);
    }
    
    public static void main(String args[])
    {
        PTMSequenceBuilder ptmBuilder = new PTMSequenceBuilder();
       //List<ProteinIndex> plist= ptmBuilder.read("C:\\Users\\Harshil\\Desktop\\", "PTM-filter.txt");
//       List<ProteinIndex> plist= ptmBuilder.read("C:\\Users\\Harshil\\Desktop\\", "PTM-filter.txt");
//       List<ProteinIndex> plist= ptmBuilder.read("C:\\Users\\Harshil\\Desktop\\", "PTM-filter.txt");
//        List<ProteinIndex> plist= ptmBuilder.read("/data/2/rpark/ip2_data/benstein/Mammalian_LKB1_AMPK_Interactome/20131122_EVh_SandP_LKB1l_SandP_10ug_each_n1_2013_11_25_18_20659/search/projects2013_11_28_14_53260/phospho", "ptm_results.txt");        
        List<ProteinIndex> plist= ptmBuilder.read("/home/harshil/", "ptm_results.txt");        
      
      ///data/2/rpark/ip2_data//benstein/Mammalian_LKB1_AMPK_Interactome/20140824_flag_IP_S1FL_wt_FM_293T_2014_08_25_10_26833/search/projects2014_08_25_11_68511/
        
      for(ProteinIndex pr : plist)
      {
          System.out.println(pr.getModSiteString());
          System.out.println(pr.getPtmPeptide());
          
//          System.out.println(pr.getAccession()+"-->" +pr.getModSiteString().charAt(0)+ pr.getChangedIndex() + "(" + pr.getModifiedArr()[209] + "/" + pr.getUnModifiedArr()[209] ) ;
      }
        System.out.println("done");
    }
    
    public List<ProteinIndex> read(String filePath,String inputFileName){
        BufferedReader br = null;
        List<ProteinIndex> proteinIndexdList = new ArrayList<>();
        try {
            File f = new File(filePath + File.separator + inputFileName);
            br = new BufferedReader(new FileReader(f));
            if(f.exists())
            {
                String currentLine = br.readLine();
                while(currentLine != null)
                {
                    if(currentLine.startsWith("P"))
                    {
                        String words[] = currentLine.split("\t");
                        ProteinIndex currentProteinIndex = new ProteinIndex(words[1], words[4]);
                        
                        StringBuffer modStringBuffer = new StringBuffer();
                        
                        // set the changed Index....
                        String[] commaWords = words[2].split(",");//S:24,T:14
                        for(String s:commaWords)
                        {
                            currentProteinIndex.addChangedIndex(Integer.parseInt(s.split(":")[1])-1);
                        }
                        //Set modified and unmodified array
                        commaWords = words[3].split(",");
                        List<Integer> changedIndex = new ArrayList<>(currentProteinIndex.getChangedIndex());
                        for(int i=0;i<commaWords.length;i++)
                        {
                            String[] s = commaWords[i].split(":");
                            currentProteinIndex.setModifiedArr(changedIndex.get(i),Integer.parseInt(s[0]));
                            currentProteinIndex.setUnModifiedArr(changedIndex.get(i),Integer.parseInt(s[1]));
                            modStringBuffer.append(words[2].split(",")[i]+"("+Integer.parseInt(s[0])+"/"+Integer.parseInt(s[1])+");");
                            
                        }
                        currentProteinIndex.setModSiteString(modStringBuffer.toString());
                        
//                        commaWords = words[5].split("[,\\[\\] ]");
//                        for(String value : commaWords)
//                        {
//                            if(value.length()>0)
//                                currentProteinIndex.addModPeptideList(value);
//                        }
//                        System.out.println(modStringBuffer);
                        proteinIndexdList.add(currentProteinIndex);
                    }
                    currentLine = br.readLine();
                }
            }
        } 
        catch (FileNotFoundException ex) 
        {
            Logger.getLogger(PTMSequenceBuilder.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(PTMSequenceBuilder.class.getName()).log(Level.SEVERE, null, ex);
        } 
        finally 
        {
            try {
                br.close();
            } 
            catch (IOException ex) {
                Logger.getLogger(PTMSequenceBuilder.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        DTASelectPhosphoReader phosphoReader = new DTASelectPhosphoReader(filePath +File.separator+"DTASelect-filter.txt.phospho");        
        phosphoReader.proteinToSiteLocalizationScore(proteinIndexdList);            
            
       return proteinIndexdList;

    }
    public void write(String filePath,String outputFileName,List<ProteinIndex> proteinIndexList)
    {
        
        BufferedWriter bw= null;
        try {
            StringBuffer sb= new StringBuffer();
            File f = new File(filePath + File.separator + outputFileName);
            bw = new BufferedWriter(new FileWriter(f));
            sb.append("H\tProteinAccession\tPTM:index\tseq_count\tproteinSequence\tModifiedSequence\n");
            for(ProteinIndex currentProteinIndex : proteinIndexList)
            {
                List<Integer> s=new ArrayList<>(currentProteinIndex.getChangedIndex());
                Collections.sort(s);
                sb.append("P\t" + currentProteinIndex.getAccession()+"\t");
                for(int index : s)
                {
                    sb.append(currentProteinIndex.getProteinSequence().charAt(index) +":" + (index+1) + ",");
                }
                sb.append("\t");
                for(int index : s)
                {
                    sb.append(currentProteinIndex.getModifiedArr()[index] + ":" + currentProteinIndex.getUnModifiedArr()[index] +",");
                }
                sb.append("\t");
                sb.append(currentProteinIndex.getProteinSequence() );
//                sb.append("\t").append(currentProteinIndex.getModPeptideList());
                sb.append("\n");
                
            }   
            bw.write(sb.toString());
            sb = new StringBuffer();
        } 
        catch (FileNotFoundException ex) {
            Logger.getLogger(PTMSequenceBuilder.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(PTMSequenceBuilder.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            try {
                bw.close();
            } catch (IOException ex) {
                Logger.getLogger(PTMSequenceBuilder.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }
    
    
    public void getIndex(Peptide peptideSequence,String proteinSequence ,ProteinIndex proteinIndex)
    {
        String[] splitter =  peptideSequence.getSequence().split("[().]");
        String searchingString="";
         for(int i=1;i<splitter.length-1;i+=3)
             searchingString += splitter[i];
        int startingIndex = proteinSequence.indexOf(searchingString);
        
        searchingString ="";
        splitter = peptideSequence.getSequence().split("[.]");
        for(int i=1;i<splitter.length-1;i++)
             searchingString += splitter[i];
        for(int i=0;i< searchingString.length();i++)
        {
            if(i<searchingString.length()-1  && searchingString.charAt(i+1) == '(')
            {
                //Its a modifide string 
                proteinIndex.addModifiedIndex(startingIndex,peptideSequence.getSpectralCount());
//                proteinIndex.setUnModifiedIndex(startingIndex,0);
                proteinIndex.addChangedIndex(startingIndex);
                while(true)
                {
                    i++;
                    if(searchingString.charAt(i) == ')')
                        break;
                }
//                proteinIndex.addModPeptideList(peptideSequence.getSequence());
            }
            else
            {
                //add the getseq to the arrray index......
                proteinIndex.addUnModifiedIndex(startingIndex,peptideSequence.getSpectralCount());
            }
            startingIndex++;
        }
    }


}
