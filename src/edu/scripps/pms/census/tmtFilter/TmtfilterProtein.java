/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.scripps.pms.census.tmtFilter;

import edu.scripps.pms.census.compare.ProteinCompare;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

/**
 *
 * @author Harshil
 */
public class TmtfilterProtein extends AbstractTmt
{
    Map<String, Protein> proteinMap = new HashMap<>();
    private int locusIndex,specCountIndex,pepCountIndex,descriptionIndex,startExperimentIndex,endExperimentIndex;
    private List<Integer> normMedianIndexList = new ArrayList<Integer>();
    private List<Integer> normAverageIndexList = new ArrayList<Integer>();
    
//    private int   count_intensity , position_spc ,position_scn ,position_cstate , position_file ;
//    private String headerSline[];
    public static void main(String args[]) throws Exception {
         Map<String, List> groupNameMap = new LinkedHashMap<>();
       List groupList= new ArrayList();
/*
        groupList.add("126.1277");
        groupList.add("127.1248");
       groupNameMap.put("one",groupList);
       
       groupList = new ArrayList();
       groupList.add("128.1344");
       groupList.add("129.1315");
        groupNameMap.put("two", groupList);
       
       groupList = new ArrayList();
       groupList.add("130.1411");
       groupList.add("131.1382");
        groupNameMap.put("three", groupList);
       
*/
       groupList.add("126.1277");
        groupList.add("127.1248");
        groupNameMap.put("control", groupList);

        groupList = new ArrayList();
        groupList.add("129.1315");
        groupList.add("130.1411");
        groupNameMap.put("test", groupList);
        
        groupList = new ArrayList();
        groupList.add("128.1344");
        groupList.add("131.1382");
        groupNameMap.put("test2", groupList);
        
        
        TmtfilterProtein tmtProtein = new TmtfilterProtein();
        tmtProtein.runTmtProteinAnalysis(groupNameMap,"C:\\Users\\Harshil\\Documents\\NetBeansProjects\\TmtFilter\\input_data","C:\\Users\\Harshil\\Documents\\NetBeansProjects\\TmtFilter\\input_data","cesus-out-proteinAnalysis.txt",true);
//        tmtProtein.runTmtProteinAnalysis(groupNameMap,args[0], "TMT_protein_comare.txt"); 
        
    }
    
    protected void printToFile(String inputPath,String directory, String fileName)
    {
        NumberFormat formatter = new DecimalFormat("##.###");
        StringBuffer output = new StringBuffer();
        int compareGroupSize =0;
//        int compareGroupSize = peptideData.get(peptideNameList.get(0)).getRatioAvg().size();
        output.append("H\t" + new Date() + "\n");
        
        output.append("H\tTMT Protein Level\n");
        output.append("H\tPath\t"+inputPath +"\n");

        List groupNameList = new ArrayList(groupKeyMap.keySet());
        
        for(int i=0;i<groupNameList.size();i++)
        {   
            output.append("H\tSample\t"+ groupNameList.get(i));
            for(int j=0;j<groupKeyMap.get(groupNameList.get(i)).size();j++)
                 output.append("\t" +groupKeyMap.get(groupNameList.get(i)).get(j));
            output.append("\n");
        }
        for(int i=0;i<groupNameList.size();i++)
        {
            for(int j=i+1;j<groupNameList.size();j++)
            {
                int counter = i+j;
                compareGroupSize++;
                output.append("H\tCompare" + counter +"\t" + groupNameList.get(i) + "/" + groupNameList.get(j)+"\n");
            }
        }
        
        output.append("H\tLocus\tSPEC_COUNT\tPEP_NUM");
        int totalGroup = ((Protein)((Collection<Protein>)proteinMap.values()).toArray()[0]).getpValueList().size();
        for(int i =0 ;i<totalGroup;i++)
        {
            output.append("\tpvalue_"+i);
        }
        
        for(int i =0 ;i<totalGroup;i++)
        {
            output.append("\tpvalueRatio_"+i);
        }
        
        for(int i =0 ;i<totalGroup;i++)
        {
            output.append("\tpvalueRatioData_"+i);
        }
        
         for(int i =0 ;i<totalGroup;i++)
        {
            output.append("\tBHCorrection_"+i);
        }
         
         for(int i =0 ;i<totalGroup;i++)
        {
            output.append("\tBHCorrectionRatio_"+i);
        }
         
        
        for(int i =0;i<groupNameList.size();i++)
        {
            output.append("\t"+groupNameList.get(i)+"_intensity_sum");
        }
        
       
//   
        
        for(int i =0;i<compareGroupSize;i++)
        {
            output.append("\tCom_" + (i+1) + "_ratio");
            output.append("\tCom_" + (i+1) + "_avg");
            output.append("\tCom_" + (i+1) + "_stdDev");
            output.append("\tCom_" + (i+1) + "_RSDValue");
        }
        
        
        output.append("\tOutlier");
        output.append("\tDESCRIPTION");
        output.append("\n");
        
        
        Iterator proteinItr = proteinMap.keySet().iterator();
        while(proteinItr.hasNext())
        {
            Protein protein = proteinMap.get(proteinItr.next());
            if(protein.isNil)
                continue;
            output.append("P\t"+protein.getLocus());
            output.append("\t"+ protein.getSpecCount());
            output.append("\t"+ protein.getPepNumber());
//            output.append("\t"+ formatter.format(protein.getpValue()));

            for(double value : protein.getpValueList())
            {
                output.append("\t"+ formatter.format(value));
            }
            
            for(double value : protein.getpValueRatioList())
            {
                output.append("\t"+ formatter.format(value));
            }
            
            for(List<Double> currentGroup : (List<List<Double>> )protein.getpValueRatioData())
            {
                output.append("\t");
                for(Double value : currentGroup)
                 {
                    output.append(formatter.format(value)+",");
                }
            }
            
            for(double value : protein.getBHCorrectionList())
            {
                output.append("\t"+ formatter.format(value));
            }
            
            for(double value : protein.getBHCorrectionRatioList())
            {
                output.append("\t"+ formatter.format(value));
            }
            
            for(List<Double> currentIntensityGroup :protein.getIntensityAvg())
            {
                output.append("\t");
                for(Double value : currentIntensityGroup)
                {
                    output.append(formatter.format(value)+",");
                }
            }
            
            
            
            for(int i =0;i<compareGroupSize;i++)
             {
                output.append("\t"+formatter.format(protein.getIntensitySumRatio().get(i)));
                output.append("\t"+formatter.format(protein.getRatioAvg().get(i)));
                output.append("\t"+formatter.format(protein.getStdDev().get(i)));
                output.append("\t"+formatter.format(protein.getRsdValue().get(i)));
             }
        
             output.append("\t"+ protein.getOutlier());
            output.append("\t"+ protein.getDescription());
            output.append("\n");
        }
        
         BufferedWriter bufferedWriter = null;
        try{
            if(directory == "")
                bufferedWriter = new BufferedWriter(new FileWriter(fileName));
            else
                bufferedWriter = new BufferedWriter(new FileWriter(directory+ File.separator+fileName));
            bufferedWriter.write(output.toString());
            bufferedWriter.close();
            System.out.println("See the output File at " + directory + File.separator +fileName );
        }
        catch(Exception e)
        {
            Logger.getLogger(TmtfilterPeptide.class.getName()).log(Level.SEVERE, null, e);
        }
    }
//    public  void runTmtProteinAnalysis(Map groupMap,String path, String outputFileName)
//    {
//        TmtfilterPeptide tmtfilter = new TmtfilterPeptide();
//        Map<String,Peptide> allPeptideMap = tmtfilter.getPeptideMap(groupMap, path);
//        Peptide pep = allPeptideMap.get("K.ASAVYQALQK.S");
//        parseInputFile(path,allPeptideMap);
//         groupNameToKeys(groupMap);
////        generateGroupInfo();
//        printToFile(path, outputFileName);
//    }
  
// only for the proteinSimple update later on    
    public  void runTmtProteinAnalysis(Map groupMap,String inputPath,String outputPath, String outputFileName, boolean normalizedValue)  {
        isNormalized = normalizedValue;
        TmtfilterPeptide tmtfilter = new TmtfilterPeptide();
        Map<String,Peptide> allPeptideMap = tmtfilter.getPeptideMap(groupMap, inputPath);
//        Peptide pep = allPeptideMap.get("K.ASAVYQALQK.S");
        parseInputFile(inputPath,allPeptideMap, groupMap);
        groupNameToKeys(groupMap);
//        generateGroupInfo();
        printToFile(inputPath,outputPath, outputFileName);
    }
    
    public void setHeaderLine(String header,Map groupMap)
    {
        if(header.contains("PLINE"))
        {
            int flag = 0;
            String words[] = header.split("\t");
            header_line = words;
            for (int i = 0; i < words.length; i++) 
            {
                if (words[i].equals("LOCUS")) {
                    locusIndex = i - 1;
                } else if (words[i].equals("SPEC_COUNT")) {
                    specCountIndex = i - 1;
                } else if (words[i].equals("PEP_NUM")) {
                    pepCountIndex = i - 1;
                } else if (words[i].equals("DESCRIPTION")) {
                    descriptionIndex = i - 1;
                } else if (flag == 0 && words[i].contains("m/z")) {
                    startExperimentIndex = i - 1;
                    flag = 1;
                } else if (flag == 1 && !words[i].contains("m/z")) {
                    endExperimentIndex = i - 1;
                    flag = -1;
                } else if (words[i].startsWith("norm_average")) {
                    normAverageIndexList.add(i-1);
                } else if (words[i].startsWith("norm_median")) {
                    normMedianIndexList.add(i-1);
                }                              
            }
        }
        else if(header.contains("SLINE"))
        {
            headerSline = header.split("\t");
        }
//       {
//           
//           headerSline = header.split("[\t]");
//        if (headerSline[1].equals("SLINE"))
//        {
//            for (int temp = 4; temp < headerSline.length; temp++) 
//            {
//                String splitter[] = new String[10];
//                splitter = headerSline[temp].split("[_()]");
////                if (splitter[0].equals("norm") && !splitter[1].equals("ratio"))
//                    normintensity_name.add(headerSline[temp]);
//                if (headerSline[temp].contains("_") && !headerSline[temp].contains("norm")) 
//                    intensityNames.add(headerSline[temp]);
//                if (headerSline[temp].equals("SpC")) 
//                {
//                    count_intensity = temp - 4; // count total number of intensity.
//                    position_spc = temp - 1;
//                    position_scn = position_spc + 1;
//                    position_cstate = position_scn + 1;
//                    position_file = position_cstate + 1;
//                    break;
//                }
////                ar.add(headerSline[temp]);
//            }
//        }
//        
//
//       }
    }
    public void parseInputFile(String path,Map<String, Peptide> allPeptideMap, Map groupMap) {
        BufferedReader br = null;
         Map<String,Peptide> peptideMap = new LinkedHashMap();
//        ProteinDetails proteinDetails = new ProteinDetails();
        try {
            File f = new File(path + File.separator + "census-out.txt");
//            File f = new File(path + File.separator + "test.txt");
            br = new BufferedReader(new FileReader(f));
            Protein proteinDetail =null;
            String currentLine = br.readLine();
            
            while (currentLine != null) {

                if (currentLine.contains("PLINE") ||currentLine.contains("SLINE") ) {
                    setHeaderLine(currentLine, groupMap);
//                    groupNameToKeys(groupKeyMap);
                }
                if (currentLine.startsWith("P")) {
                    String words[] = currentLine.split("\t");
                     proteinDetail = new Protein();
                    proteinDetail.setLocus(words[locusIndex]);
                    proteinDetail.setSpecCount(words[specCountIndex]);
                    proteinDetail.setPepNumber(words[pepCountIndex]);
                    proteinDetail.setDescription(words[descriptionIndex]);
                    
                    for(Iterator<Integer> itr=this.normAverageIndexList.iterator(); itr.hasNext(); ) {
                        itr.next();
                    }
                    
                    for(Iterator<Integer> itr=this.normAverageIndexList.iterator(); itr.hasNext(); ) {
                        itr.next();
                    }
//                    proteinDetail.setExperimentGroup(words);
                     peptideMap = new LinkedHashMap();

                }
                if(currentLine.startsWith("S"))
                {
                   
                    while( currentLine != null && currentLine.startsWith("S") )
                    {
                        
                        
                      //  System.out.println(currentLine);
                        
                       String[] words = currentLine.split("\t");
                       String currentPeptideName = words[2];
                       if(allPeptideMap.containsKey(currentPeptideName))
                            peptideMap.put(currentPeptideName, allPeptideMap.get(currentPeptideName) );
                       
                       
                               
                       currentLine=br.readLine();
                    }
                    
                    groupKeyMap.get(0);
                    
                    
                    proteinDetail.setPeptideMap(peptideMap);
                    proteinDetail.filterValues(0.05);
                  //   proteinDetail.generateRatio();
                     
//                     proteinDetail.generatePValue();
                    proteinMap.put(proteinDetail.getLocus(), proteinDetail);
                    continue;
                }
//                if(currentLine != null && currentLine.startsWith("S"))
                     currentLine = br.readLine();
            }
        } 
        catch (FileNotFoundException ex) {
            Logger.getLogger(ProteinCompare.class.getName()).log(Level.SEVERE, null, ex);
        } 
        catch (IOException ex) {
            Logger.getLogger(ProteinCompare.class.getName()).log(Level.SEVERE, null, ex);
        } 
        finally 
        {
            try {
                br.close();
            } catch (IOException ex) {
                Logger.getLogger(ProteinCompare.class.getName()).log(Level.SEVERE, null, ex);
                }
        }
        
        
        
        
    }

  /*  public String setPeptide(String line , Map<String, Peptide> peptideData)
    {
        String sequence = null;
        boolean checkzero = false;
        Peptide pep_details = null;
        String[] words_peptide = line.split("[\t]");
        for (int temp = 0; temp < count_intensity; temp++) {
            if (Float.parseFloat(words_peptide[temp + 4 - 1]) == 0) {
                checkzero = true;
                break;
            }
        }
        if (checkzero) {
            return null;
        }
        if (peptideData.containsKey(words_peptide[2])) {
            pep_details = peptideData.get(words_peptide[2]);
            pep_details.setPeptideCounter();
            if (!checkValidEntryForpeptide(pep_details, Integer.parseInt(words_peptide[position_scn]), words_peptide[position_file])) {
                return null;
            }
            for (int temp = 0; temp < count_intensity; temp++) {
                pep_details.setallIntensity(headerSline[temp + 4], Float.parseFloat(words_peptide[temp + 4 - 1]));
            }
            setPeptideDetails(pep_details, words_peptide);
        } else {
            for (int temp = 0; temp < count_intensity; temp++) {
                if (Float.parseFloat(words_peptide[temp + 4 - 1]) == 0) {
                    checkzero = true;
                    break;
                }
            }
            if (checkzero) {
                return null;
            } //                         if(!peptideNameList.contains(words_peptide[2]))
            //                            peptideNameList.add(words_peptide[2]);
            else {
                System.out.println("Error in addinf peptide at the peptide list");
            }
            for (int temp = 0; temp < count_intensity; temp++) {
                pep_details.setallIntensity(headerSline[temp + 4], Float.parseFloat(words_peptide[temp + 4 - 1]));
            }
            pep_details.setPeptideCounter();
            setPeptideDetails(pep_details, words_peptide);
            pep_details.setSequence(words_peptide[2]);
        }
        pep_details.addProteinDetails(words[1]);
        peptideData.put(words_peptide[2], pep_details);
//                         hash.add_data(words_peptide[2], words[1]);
        return sequence;

                 }
    */
//public void setPeptideDetails(Peptide pep_details,String[] words_peptide) {
//        pep_details.setSpc(Integer.parseInt(words_peptide[position_spc]));
//        pep_details.addScannum(Integer.parseInt(words_peptide[position_scn]));
//        pep_details.setCstate(Integer.parseInt(words_peptide[position_cstate]));
//     //   pep_details.setratio(normintensity_name);
//        pep_details.setFilename(words_peptide[position_file]);
//    }
}
