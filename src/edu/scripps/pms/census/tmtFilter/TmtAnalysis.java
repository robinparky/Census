///*
// * To change this template, choose Tools | Templates
// * and open the template in the editor.
// */
//package edu.scripps.pms.census.tmtFilter;
//
//import static edu.scripps.pms.census.tmtFilter.TmtfilterPeptide.intensityNames;
//import java.io.BufferedReader;
//import java.io.BufferedWriter;
//import java.io.File;
//import java.io.FileReader;
//import java.io.FileWriter;
//import java.io.IOException;
//import java.util.ArrayList;
//import java.util.HashMap;
//import java.util.List;
//import java.util.Map;
//import java.util.Set;
//import java.util.logging.Level;
//import java.util.logging.Logger;
//import javax.print.attribute.HashAttributeSet;
//import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
//
///**
// *
// * @author Harshil
// * This script is for brain(other lab) to make a specified TMT version printToFile only... 
// * its only for his work... 
// */
//public class TmtAnalysis 
//{
//    multihash hash = null;
//     Map groupMap =null;
//     Map<String,Peptide> peptide_data = null;
//    ArrayList<String> normintensity_name = null;
//    Map intensityNameToIndex = new HashMap();
//    public static void main(String args[])
//    { 
//           TmtAnalysis tmtAnalysis = new TmtAnalysis();
//          String path = "C:\\Users\\Harshil\\Documents\\NetBeansProjects\\censusFiles\\2014_02_11_11_9198";
//         tmtAnalysis.groupMap= new HashMap();
//         tmtAnalysis.generateGroupList();
//        
//        tmtAnalysis.parseCensusFile(path);
//        tmtAnalysis.generateIntensityMap();
//        System.out.println("Data are stored in the object....");
////        tmtAnalysis.findStdDev();
//        tmtAnalysis.printToFile(path);
//        
//    }
//    public void generateGroupList()
//    {
//        List groupList= new ArrayList();
//        groupList.add("126.1277");
//       groupMap.put("1",groupList);
//       groupList = new ArrayList();
//       groupList.add("127.1248");
//        groupMap.put("two", groupList);
//        groupList = new ArrayList();
//       groupList.add("128.1344");
//        groupMap.put("two1", groupList);
//        groupList = new ArrayList();
//       groupList.add("129.1315");
//        groupMap.put("two2", groupList);
//        groupList = new ArrayList();
//       groupList.add("130.1411");
//        groupMap.put("two3", groupList);
//        groupList = new ArrayList();
//       groupList.add("131.1382");
//        groupMap.put("two4", groupList);                
//    }
//    public void printToFile(String path)
//    {   
//        BufferedWriter bufferedWriter = null;
//        try 
//        {
//             bufferedWriter = new BufferedWriter(new FileWriter(new File(path + File.separator + "census_peptide_analysis.txt")));
//             bufferedWriter.write("H"+"\tPeptide_Name\t");
//             for(String currentGroup : (Set<String>)groupMap.keySet())
//             {
//                    for(int z=0;z<(((List) groupMap.get(currentGroup)).size()) ;z++)
//                    {
//                        bufferedWriter.write("Avg_" +  ((List)groupMap.get(currentGroup)).get(z)+ "\t");
//                        bufferedWriter.write("STD_" + ((List)groupMap.get(currentGroup)).get(z) + "\t");
//                        bufferedWriter.write("Avg_norm_" +  ((List)groupMap.get(currentGroup)).get(z)+ "\t");
//                        bufferedWriter.write("STD_norm_" + ((List)groupMap.get(currentGroup)).get(z) + "\t");
//     
//                    }
//                    
//                      
//                }
//                bufferedWriter.write("SpC"+"\t"+"ScanNum"+"\t"+"Cstate"+"\t"+"Filename"+"\t"+"Protein_Set");
//                bufferedWriter.write("\n");
//                Peptide currentPeptideData = new Peptide();
//                for(String currentPeptideName : peptide_data.keySet())
//                {
//                    currentPeptideData = peptide_data.get(currentPeptideName);
//                   bufferedWriter.write("S\t"+ currentPeptideName);
//                   for(String currentGroup : (Set<String>)groupMap.keySet())
//                     {
//                        for(int z=0;z<(((List) groupMap.get(currentGroup)).size()) ;z++)
//                        {
//                            String currentIntensity = (String) ((List)groupMap.get(currentGroup)).get(z);
//                            bufferedWriter.write("\t"  + currentPeptideData.getstat(currentIntensity).getMean()  );
//                            bufferedWriter.write("\t" + currentPeptideData.getstat(currentIntensity).getStandardDeviation()  );
//                            bufferedWriter.write("\t"  + currentPeptideData.getstat("norm_"+currentIntensity).getMean()  );
//                            bufferedWriter.write("\t" + currentPeptideData.getstat("norm_" +currentIntensity).getStandardDeviation()  );
//                            
//                        }
//                    }
//                bufferedWriter.write("\t"+currentPeptideData.getSpc() +"\t"+ currentPeptideData.getScannum()+"\t"+currentPeptideData.getCstate()+"\t"+currentPeptideData.getFilename());
//                bufferedWriter.write("\t" + hash.getProteinSet(currentPeptideName));
//                bufferedWriter.write("\n");
//                        
//                }
//                
//                 bufferedWriter.close();
//        }
//        catch (IOException ex) {
//            Logger.getLogger(TmtAnalysis.class.getName()).log(Level.SEVERE, null, ex);
//        }
//      
//    }
//    //This method will map the grouped Intensity to its actual position in teh List
//    public void generateIntensityMap()
//    {
//        Map newGroupMap = new HashMap();
////       List keyset = (String[]) intensityNameToIndex.keySet().toArray();
//       List keyset = new ArrayList();
//       keyset.addAll(intensityNameToIndex.keySet());
//        for(String currentGroup : (Set<String>)groupMap.keySet())
//             {
//                 List ar= new ArrayList();
//                    for(int z=0;z<(((List) groupMap.get(currentGroup)).size()) ;z++)
//                    {
//                        for(int j=0;j<keyset.size();j++)
//                        {
////                            if(intensityNameToIndex.get(keyset.get(j)  ))
//                            
//                            
//                            
//                            if(((String)(keyset.get(j))).contains((String)((List) groupMap.get(currentGroup)).get(z)) && !(((String)(keyset.get(j))).contains("norm")))
//                            {
//                                ar.add(((String)keyset.get(j)));
//                                break;
//                            }
//                        }
////                        bufferedWriter.write("Avg_" +  ((List)groupNameMap.get(currentGroup)).get(z)+ "\t");
////                        bufferedWriter.write("STD_" + ((List)groupNameMap.get(currentGroup)).get(z) + "\t");
//                    }
//                    newGroupMap.put(currentGroup, ar);
//                }
//        groupMap = newGroupMap;
//        System.out.println("");
//    }
//    
//            public void parseCensusFile(String path)
//    {
//        BufferedReader br = null;
//        String protine,peptide,directory="";
//        String words[]=new String[10000];
//        String words_peptide[]=new String[10000];
//        String header_line[]=new String[150];
//         hash=new multihash();
//        ArrayList<String> testList=new ArrayList<>();
//        ArrayList<String> peptide_list=new ArrayList<>();
//         normintensity_name=new ArrayList<>();
//        peptide_data=new HashMap<>(); // The main map for contains the peptide sequence name and all the datails of each peptide sequence.
//        DescriptiveStatistics stats = new DescriptiveStatistics();
//
//        int count_intensity = 0;
//        boolean checkzero=false;
//        int position_spc = 0,position_scn = 0,position_cstate = 0,position_file = 0;
//	System.out.println("Please wait...");
//        try
//        {
//          File f = new File(path + File.separator + "census-out.txt");
//            br = new BufferedReader(new FileReader(f));
//            String line;
//            line=br.readLine();
//            Peptide pep_details = new Peptide();        
//           while(line != null)
//            {
//                checkzero=false;
//                pep_details=new Peptide();
//                stats = new DescriptiveStatistics();
//                if(line.startsWith("H") || line.startsWith("h")) 
//                {
//                   // System.out.println("Its a header line.");
//                    header_line= line.split("[\t]");
//                        if(header_line[1].equals("SLINE"))
//                        {
//                           for(int temp=4;temp<header_line.length;temp++)
//                           {
//                               String splitter[] =new String[10];
//                               splitter=header_line[temp].split("[_()]");
//                               if(splitter[0].equals("norm") && !splitter[1].equals("ratio"))
//                                   normintensity_name.add(header_line[temp]);
//                               if(header_line[temp].contains("_") && !header_line[temp].contains("norm"))
//                                   intensityNames.add(header_line[temp]);
//                               if(header_line[temp].contains("m/z"))
//                                   intensityNameToIndex.put(header_line[temp], temp);
//                               if(header_line[temp].equals("SpC"))
//                               {    
//                                   count_intensity=temp-4; // count total number of intensity.
//                                   position_spc=temp-1;
//                                   position_scn=position_spc+1;
//                                   position_cstate=position_scn+1;
//                                   position_file=position_cstate+1;  
//                                   break;
//                               }
//                               testList.add(header_line[temp]);
//                            }
//                        }
//                }
//                if(line.startsWith("P") || line.startsWith("p"))
//                    words=line.split("[\t]");
//                if(line.startsWith("S") || line.startsWith("s"))
//                 {
//                     //System.out.print("Its a Peptide line.");
//                     words_peptide=line.split("[\t]");
//                     for(int temp=0;temp<count_intensity;temp++)
//                     {
//                          if(Float.parseFloat(words_peptide[temp+4-1]) == 0)
//                            {
//                                checkzero=true;
//                                break;
//                            }
//                     }
//                     if(checkzero)
//                         {
//                            line=br.readLine();
//                            continue;
//                         }   
//                     if(peptide_data.containsKey(words_peptide[2]))
//                     {
//                         pep_details=peptide_data.get(words_peptide[2]);
//                         
//                         int peptide_count=pep_details.getPeptide_counter(words_peptide[2]);
//                         pep_details.setPeptide_counter(words_peptide[2]);
//                         if(!checkValidEntryForpeptide(pep_details,Integer.parseInt(words_peptide[position_scn]),words_peptide[position_file]))
//                         {
//                             line=br.readLine();
//                             continue;
//                         }
//                         for(int temp=0;temp<count_intensity;temp++)
//                        {
//                   
//                            float z=pep_details.getAverage(header_line[temp+4]) ;
//                            float new_average=(float) ((z*peptide_count + Float.parseFloat(words_peptide[temp+4-1]))/(peptide_count+1.0));
//                            pep_details.setAverage(header_line[temp+4], new_average );
//                            pep_details.setallIntensity(header_line[temp+4], Float.parseFloat(words_peptide[temp+4-1]));
//                            // Adding the data to the standard deviation to find it.   
//                            if(header_line[temp+4].equals("m/z_126.1277_int"))
//                            {
//                                System.out.println(z+"st time of "+ new_average);
//                            }
//                            stats=pep_details.getstat(header_line[temp+4]);
//                            stats.addValue(Float.parseFloat(words_peptide[temp+4-1]));
//                            pep_details.setstat(header_line[temp+4],stats);
//                         }
//                         pep_details.setSpc(Integer.parseInt(words_peptide[position_spc]));
//                         pep_details.addScannum(Integer.parseInt(words_peptide[position_scn]));
//                         pep_details.setCstate(Integer.parseInt(words_peptide[position_cstate]));
//                         pep_details.setratio(normintensity_name);
//                         pep_details.setFilename(words_peptide[position_file]);
//                     }
//                     else
//                     {
//                         for (int temp = 0; temp < count_intensity; temp++) {
//                             if (Float.parseFloat(words_peptide[temp + 4 - 1]) == 0) {
//                                 checkzero = true;
//                                 break;
//                             }
//                         }
//                         if (checkzero)
//                         {
//                            line=br.readLine();
//                            continue;
//                         }   
//                         if(!peptide_list.contains(words_peptide[2]))
//                            peptide_list.add(words_peptide[2]);
//                         else
//                           System.out.println("Error in addinf peptide at the peptide list");
//                         for(int temp=0;temp<count_intensity;temp++)
//                        {
//                            stats = new DescriptiveStatistics();
//                            pep_details.setAverage(header_line[temp+4],Float.parseFloat(words_peptide[temp+4-1] ));
//                            stats.addValue(Float.parseFloat(words_peptide[temp+4-1]));
//                            pep_details.setstat(header_line[temp+4],stats);
//                            pep_details.setallIntensity(header_line[temp+4], Float.parseFloat(words_peptide[temp+4-1]));
//                            if(header_line[temp+4].equals("m/z_126.1277_int"))
//                            {
//                                System.out.println("1st time of "+ words_peptide[temp+4-1]);
//                            }
//                         }
//                         pep_details.setPeptide_counter(words_peptide[2]);
//                         pep_details.setSpc(Integer.parseInt(words_peptide[position_spc]));
//                         pep_details.addScannum(Integer.parseInt(words_peptide[position_scn]));
//                         pep_details.setCstate(Integer.parseInt(words_peptide[position_cstate]));
//                         pep_details.setratio(normintensity_name);
//                         pep_details.setFilename(words_peptide[position_file]);
//                     }
//                         peptide_data.put(words_peptide[2],pep_details);
//                         hash.add_data(words_peptide[2], words[1]);// PepSeq----IPI name
//                 }
//                line=br.readLine();
//            }
//        }   
//        catch (IOException ex)
//        {
//            Logger.getLogger(TmtAnalysis.class.getName()).log(Level.SEVERE, null, ex);
//        }
//    }
//      public boolean checkValidEntryForpeptide(Peptide pep_detail,int newScanNumber,String newFileName)
//      {
//          ArrayList scanNumber = pep_detail.getScannum();
//          Set fileName =  pep_detail.getFilename();
//          if(scanNumber.contains(newScanNumber) && fileName.contains(newFileName))
//              return false;
//          else
//              return true;
//      }
//}
