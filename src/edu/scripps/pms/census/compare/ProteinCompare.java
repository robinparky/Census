/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.scripps.pms.census.compare;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.math.stat.inference.TestUtils;

/**
 *
 * @author Harshil
 */
public class ProteinCompare {

    public static ProteinDetails proteinDetails = new ProteinDetails();
    public static ProteinDetails deletedProteinDetails;
    public static List experimentGroupNames = new ArrayList();
    public static List experimentGroupAttribute = new ArrayList();;

    public static void main(String args[]) {
/*      
*/
        /*
        pc = new ProteinCompare();
        
        pc.runTmtCompare(map,"C:\\Users\\Harshil\\Documents\\NetBeansProjects\\TmtFilter\\input_data", "census-anova.txt");
        */
        Map map= new HashMap();
        List ls= new ArrayList();
        ls.add("126.1283 mz");
        ls.add("127.1316 mz");
//        ls.add("128.135 mz");
        map.put("one", ls);
        List ls1= new ArrayList();
//        ls1.add("129.1383 mz");
        ls1.add("130.1417 mz");
        ls1.add("131.1387 mz");
        map.put("two", ls1);
        
        ProteinCompare.runTmtCompare(map,"C:\\Users\\Harshil\\Documents\\NetBeansProjects\\2012", "census-anova.txt");
    }
    public static void runTmtCompare(Map groupList, String path, String outputFileName)
    {
        System.out.println("Please Wait......");
        experimentGroupAttribute.addAll((Collection)groupList.keySet());
        Iterator itr = experimentGroupAttribute.iterator();
        while(itr.hasNext())
        {
            List ls = (List) groupList.get(itr.next());
            double[] d = new double[ls.size()];
            for (int i = 0; i < ls.size(); i++) {

                d[i] = (double) Double.parseDouble(((String) ls.get(i)).split(" ")[0]);
            }
            experimentGroupNames.add(d);
        }
        parseInputFile(path);
        proteinDetails.generateAnovaValues();
        solveNaNValues();
        proteinDetails.setVariancePvalueBH();
       // addDeletedEntries();
        generateOutputFile(path,outputFileName);
    }
    public static void solveNaNValues()
    {
        deletedProteinDetails = new ProteinDetails();
        proteinDetails.removeNaNEntries(deletedProteinDetails);

    }
  
    public static void parseInputFile(String path) {
        BufferedReader br = null;
//        ProteinDetails proteinDetails = new ProteinDetails();

        try {
            File f = new File(path + File.separator + "census-out.txt");
//            File f = new File(path + File.separator + "test.txt");
            br = new BufferedReader(new FileReader(f));
            String currentLine = br.readLine();
            while (currentLine != null) {

                if (currentLine.contains("PLINE")) {
                    proteinDetails.setHeaderLine(currentLine);
                    proteinDetails.parseHeaderLine(currentLine);
                    proteinDetails.setExperimentGroupNameIndexs(experimentGroupNames);
                    proteinDetails.setExperimentGroupAttribute(experimentGroupAttribute);
                }
                if (currentLine.startsWith("P")) {
                    String words[] = currentLine.split("\t");
                    proteinDetails.setLocus(words);
                    proteinDetails.setSpecCount(words);
                    proteinDetails.setPepCount(words);
                    proteinDetails.setDescriptions(words);
                    proteinDetails.setExperimentGroup(words);

                }
                currentLine = br.readLine();
            }
            proteinDetails.setStdDev();
//            proteinDetails.setVariancePvalueBH();
        } catch (FileNotFoundException ex) {
            Logger.getLogger(ProteinCompare.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(ProteinCompare.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            try {
                br.close();
            } catch (IOException ex) {
                Logger.getLogger(ProteinCompare.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }

    public static void generateOutputFile(String path , String fileName) {
        BufferedWriter bw = null;
        try {
            File f = new File(path + File.separator + fileName);
            bw = new BufferedWriter(new FileWriter(f));
            StringBuffer sb = new StringBuffer();
            int totalNumberOfGroups = experimentGroupNames.size();
            List groupNames = new ArrayList();
            // sb.add
            for (int i = 0; i < proteinDetails.getExperimentGroupAttribute().size(); i++) {
                sb.append("H\tGroup\t");
                sb.append(proteinDetails.getExperimentGroupAttribute().get(i).toString() + "\t");
                StringBuffer sb1 = new StringBuffer();
                for (int j = 0; j < ((double[]) experimentGroupNames.get(i)).length; j++) {
                    sb1.append(((double[]) experimentGroupNames.get(i))[j]);
                    sb1.append(",");
                }
                groupNames.add(sb1.toString());
                sb.append(sb1 + "\n");
            }
            sb.append("H\tLocus");
            for (int j = 0; j < totalNumberOfGroups; j++) {
                sb.append("\tnorm_intensity_" + proteinDetails.getExperimentGroupAttribute().get(j).toString());
            }
            sb.append("\tpValue\tfValue\tpValue_BHCorrections");
            for (int j = 0; j < totalNumberOfGroups; j++) {
                sb.append("\tStdDev_" + proteinDetails.getExperimentGroupAttribute().get(j).toString());
            }
            for (int j = 0; j < totalNumberOfGroups; j++) {
                sb.append("\tAverage_" + proteinDetails.getExperimentGroupAttribute().get(j).toString());
            }
            for (int j = 0; j < totalNumberOfGroups; j++) {
                sb.append("\tRsdValue_" + proteinDetails.getExperimentGroupAttribute().get(j).toString());
            }
            sb.append("\tDescription");
            sb.append("\n");
            for (int i = 0; i < proteinDetails.getLocus().size(); i++) //            for(int i=0;i<50;i++)
            {

                sb.append("P\t" + proteinDetails.getLocus().get(i) + "\t");
                for (int j = 0; j < totalNumberOfGroups; j++) {
                    for (int k = 0; k < ((double[]) ((List) proteinDetails.getExperimentGroup().get(i)).get(j)).length; k++) {
                        sb.append(((double[]) ((List) proteinDetails.getExperimentGroup().get(i)).get(j))[k]);
                        sb.append(",");
                    }
                    sb.append("\t");
                }

                sb.append(proteinDetails.getpValues().get(i));
                sb.append("\t" + proteinDetails.getfValues().get(i));
                sb.append("\t" + proteinDetails.getVariancePvalueBH().get(i));
                for (int j = 0; j < totalNumberOfGroups; j++) {
                    sb.append("\t" + ((ArrayList) proteinDetails.getStdDev().get(i)).get(j));
                }
                for (int j = 0; j < totalNumberOfGroups; j++) {
                    sb.append("\t" + ((ArrayList) proteinDetails.getAverage().get(i)).get(j));
                }
                for (int j = 0; j < totalNumberOfGroups; j++) {
                    sb.append("\t" + ((ArrayList) proteinDetails.getRsdValue().get(i)).get(j));
                }
                sb.append("\t" + proteinDetails.getDescriptions().get(i));
                sb.append("\n");
            }
            bw.write(sb.toString());
//            System.out.println(sb);
            System.out.println("Find the File at ->" + path + File.separator + "census-anova.txt");
        } catch (IOException ex) {
            Logger.getLogger(ProteinCompare.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            try {
                bw.close();
            } catch (IOException ex) {
                Logger.getLogger(ProteinCompare.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

    }
}
