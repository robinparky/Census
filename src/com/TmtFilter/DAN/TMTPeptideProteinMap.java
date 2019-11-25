/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package com.TmtFilter.DAN;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

/**
 *String 
 * @author Harshil
 * This is script is for Dan who need to convert TMT peptides to proteins
 */
public class TMTPeptideProteinMap {

    public Map proteinMap = new HashMap();
    int counter = 0;

    public static void main(String args[]) {
        System.out.println("Generating the file..");
        TMTPeptideProteinMap pf = new TMTPeptideProteinMap();
 // parameters should be the filePath and inputFileName........       
        
        if(args.length<=0) {
            System.out.println("Usage: TMTPeptideProteinMap path census-out.txt");
            System.exit(0);
        }
        
        pf.parseInputFile(args[0],args[1]);
        pf.setProteinDescription(args[0]);
        pf.writeTotheFile(args[0],"ProteinAnalysis.txt");
        System.out.println("Look the file in " + args[0]+File.separator+"ProteinAnalysis.txt");
//        System.out.println(((TmtProteinDetail) pf.proteinMap.get("IPI00129268.1")).getProteinDesc());
//        System.out.println(((TmtProteinDetail) pf.proteinMap.get(" IPI00654755.3")));

    }
    /*    
     //if U wish to use Census File---------    
     public void setProteinDescription (String path)
     {
     BufferedReader br= null;
     try {      

     br = new BufferedReader(new FileReader(new File(path + File.separator + "census-out.txt")));
     String currentLine = br.readLine();
     while(currentLine!= null)
     {
     if(currentLine.startsWith("P"))
     {
     String words[] = currentLine.split("\t");
     String proteinDesc = words[words.length -1];
     String proteinName = words[1];
     if(proteinMap.containsKey(proteinName))
     {
     ((TmtProteinDetail)proteinMap.get(proteinName)).setProteinDesc(proteinDesc);
     }
     }
     currentLine= br.readLine();
     }
     } catch (FileNotFoundException ex) {
     Logger.getLogger(TMTPeptideProteinMap.class.getName()).log(Level.SEVERE, null, ex);
     } catch (IOException ex) {
     Logger.getLogger(TMTPeptideProteinMap.class.getName()).log(Level.SEVERE, null, ex);
     }
            
            
     }
     */

    public void setProteinDescription(String path) {
// it should receive the DTASelect-Fileter.txt file......         
        BufferedReader br = null;
        try {

            br = new BufferedReader(new FileReader(new File(path + File.separator + "DTASelect-filter.txt")));
            String currentLine = br.readLine();
            while (currentLine != null) {
                if (currentLine.startsWith("IPI")) {
                    String words[] = currentLine.split("\t");
                    String proteinDesc = words[words.length - 1];
                    String proteinName = words[0].split(":|\\|")[1];
                    if (proteinName.equals("IPI00654755.3")) {
                        System.out.println("");
                    }
                    if (proteinMap.containsKey(proteinName)) {
                        TmtProteinDetail pd = (TmtProteinDetail) proteinMap.get(proteinName);
                        pd.setProteinDesc(proteinDesc);
                        proteinMap.put(proteinName, pd);

//                        ((TmtProteinDetail)proteinMap.get(proteinName)).setProteinDesc(proteinDesc);
                    }

                }
                currentLine = br.readLine();
            }
        } catch (FileNotFoundException ex) {
            Logger.getLogger(TMTPeptideProteinMap.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(TMTPeptideProteinMap.class.getName()).log(Level.SEVERE, null, ex);
        }


    }

    public void parseInputFile(String path,String fileName) {
// It should receive the census-out.txt file.......        
        BufferedReader br = null;
        String proteinNames[];
        int indexOfprotein = 0;
        try {

            br = new BufferedReader(new FileReader(new File(path + File.separator + fileName)));
            String currentLine = br.readLine();

            while (currentLine != null) {
                String words[] = currentLine.split("\t");
                if (currentLine.startsWith("H")) 
                {
                        indexOfprotein = words.length -1;
//                    indexOfprotein = 42;
                }
                if (currentLine.startsWith("S")) {
                    proteinNames = words[indexOfprotein].split("\\[|,|\\]| ");
                    TmtProteinDetail tmtProteinDetail = new TmtProteinDetail();

                    for (String currentProteinName : proteinNames) {
                        if (currentProteinName.equals(" IPI00654755.3")) {
                            System.out.println("");
                        }
                        if (currentProteinName.contains(" ")) {
                            System.out.println("");
                        }
                        if (currentProteinName.equals("")) {
                            continue;
                        }
                     
                        DescriptiveStatistics avg[] = new DescriptiveStatistics[3];
                        DescriptiveStatistics stdDev[] = new DescriptiveStatistics[3];
                        if (proteinMap.containsKey(currentProteinName)) {
                            tmtProteinDetail = (TmtProteinDetail) proteinMap.get(currentProteinName);
//                           tmtProteinDetail.addPeptideCount();
                            avg = tmtProteinDetail.getAverage();
                            counter++;

                            stdDev = tmtProteinDetail.getStdDev();
//                           for(int i=1;i<=3;i++)
//                           {
//                                avg.addValue(Float.parseFloat(words[indexOfprotein+i]));
//                                stdDev.addValue(Float.parseFloat(words[indexOfprotein+i+3]));
//                           }
//                           tmtProteinDetail.setAverage(avg);
//                           tmtProteinDetail.setStdDev(stdDev);
                        } else {
                            tmtProteinDetail = new TmtProteinDetail();
                            for(int i =0;i<3;i++)
                            {
                                avg[i]= new DescriptiveStatistics();
                                stdDev[i]= new DescriptiveStatistics();
                            }
                        }
                        tmtProteinDetail.addPeptideCount();
                        tmtProteinDetail.setProtein(currentProteinName);
                        for (int i = 1; i <= 3; i++) {
                            avg[i - 1].addValue(Float.parseFloat(words[indexOfprotein + i]));
                            stdDev[i - 1].addValue(Float.parseFloat(words[indexOfprotein + i + 3]));
                        }
                        
                        tmtProteinDetail.setAverage(avg);
                        tmtProteinDetail.setStdDev(stdDev);
//                        System.out.println(currentProteinName);
                        proteinMap.put(currentProteinName, tmtProteinDetail);
                    }

                }
                currentLine = br.readLine();

            }
        } catch (FileNotFoundException ex) {
            Logger.getLogger(TMTPeptideProteinMap.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(TMTPeptideProteinMap.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            try {
                br.close();
            } catch (IOException ex) {
                Logger.getLogger(TMTPeptideProteinMap.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }

    public void writeTotheFile(String path,String fileName) {
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(new File(path + File.separator + fileName)));
            StringBuffer sb = new StringBuffer();
            sb.append("H\tProtein_Name\tPeptideCounts\t");
            for (int i = 0; i < 3; i++) {
                sb.append("Avg_" + (i + 1) + "\t");
            }
            for (int i = 0; i < 3; i++) {
                sb.append("StdDev_" + (i + 1) + "\t");
            }
            sb.append("description");
            sb.append("\n");
            for (String currentProteinId : (Set<String>) proteinMap.keySet()) {
                sb.append("P\t" + currentProteinId);
                TmtProteinDetail details = (TmtProteinDetail) proteinMap.get(currentProteinId);
                sb.append("\t" + details.getPeptideCount());
                for (int i = 0; i < 3; i++) {
                    sb.append("\t" + details.getAverage()[i].getMean());
                }
                for (int i = 0; i < 3; i++) {
                    sb.append("\t" + details.getAverage()[i].getStandardDeviation());
                }
                sb.append("\t" + details.getProteinDesc());
                sb.append("\n");

            }
            bw.write(sb.toString());
            bw.close();
            
        } catch (IOException ex) {
            Logger.getLogger(TMTPeptideProteinMap.class.getName()).log(Level.SEVERE, null, ex);
        }

    }
}
