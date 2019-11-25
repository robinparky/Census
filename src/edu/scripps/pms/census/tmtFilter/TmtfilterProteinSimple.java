/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.scripps.pms.census.tmtFilter;

//import static edu.scripps.pms.util.stats.Fisher.combinedProbabilityTest;
import edu.scripps.pms.stats.TTest;
import java.sql.Array;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
//import org.apache.commons.math.MathException;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.inference.TestUtils;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;

/**
 *Census2 : Isobaric Labelling Data Analysis in an Automatic way. Research Paper
 * for more details......
 * @author Harshil
 */
public class TmtfilterProteinSimple extends TmtfilterProtein {

    public static void main(String args[]) throws Exception {
        Map<String, List> groupNameMap = new LinkedHashMap<>();

        List groupList = new ArrayList();
        groupList.add("126.127725");
        groupList.add("127.12476");
        groupList.add("128.134433");
        groupNameMap.put("control1", groupList);
        
        
        
        

//        groupList.add("130.141141");
        
        groupList = new ArrayList();
        groupList.add("129.131468");
        groupList.add("130.141141");
        groupList.add("131.138176");
//        groupList.add("131.138176"); 

        groupNameMap.put("bell", groupList);

        groupList = new ArrayList();
        
        

//        groupList = new ArrayList();
//        groupList.add("130.141141");
//        groupList.add("131.138176"); 
//        groupNameMap.put("dell", groupList);

                
        java.util.Set<String> s = groupNameMap.keySet();
        for(Iterator<String> itr=s.iterator(); itr.hasNext(); ) {
            String str = itr.next();
            System.out.println(groupNameMap.get(str));
        }

        TmtfilterProteinSimple tmtProteinSimple = new TmtfilterProteinSimple();
//        TmtfilterProteinSi tmtProtein = new TmtfilterProtein();

//        tmtProteinSimple.runTmtProteinAnalysis(groupNameMap, "E:\\", "E:\\","cesus-out-proteinAnalysisSimple.txt",false);
//        tmtProteinSimple.runTmtProteinAnalysis(groupNameMap, "E:\\tmt\\newData", "E:\\tmt\\newData","cesus-out-proteinAnalysisSimple.txt",false);        
//tmtProteinSimple.runTmtProteinAnalysis(groupNameMap, "/data/2/rpark/ip2_data/dmcclat/IPS_Pranav/expt1_ExoTMT_2015_07_10_09_33407/quant/2015_07_13_15_11243", "/home/rpark","cesus-out-proteinAnalysisSimple.txt",false);                
tmtProteinSimple.runTmtProteinAnalysis(groupNameMap, "/home/rpark", "/home/rpark","cesus-out-proteinAnalysisSimple.txt",false);                
        

//        tmtProteinSimple.runTmtProteinAnalysis(groupNameMap,"/home/rpark/test_data/meha_tmt", "new_compare.txt");
    }

    @Override
    public void runTmtProteinAnalysis(Map groupMap, String inputpath,String outputPath, String outputFileName, boolean normalizedValue)  {
        isNormalized = normalizedValue;
        TmtfilterPeptideSimple tmtfilter = new TmtfilterPeptideSimple();
        Map<String, Peptide> allPeptideMap = tmtfilter.getPeptideMap(groupMap, inputpath);
        
//        Peptide pep = allPeptideMap.get("K.ASAVYQALQK.S");
        parseInputFile(inputpath, allPeptideMap, allPeptideMap);
        groupNameToKeys(groupMap);
//        generateGroupInfo();
        generatePValues();
        generateBHCorrection();
        generatePValuesRatio();
        generateBHCorrectionRatio();
        printToFile(inputpath,outputPath, outputFileName);
    }

    private void generateBHCorrection()
    {
        int totalGroup = ((Protein)((Collection<Protein>)proteinMap.values()).toArray()[0]).getpValueList().size();
        for(int i=0;i<totalGroup ; i++)
        {
           Iterator proteinItr = proteinMap.keySet().iterator();
           List<Double> adjestedPValue = new ArrayList<>();
           while (proteinItr.hasNext()) {
            Protein currentProtein = proteinMap.get(proteinItr.next());
            adjestedPValue.add(currentProtein.getpValueList().get(i));
            }
            adjestedPValue = TmtfilterPeptideSimple.performBHCorrection(adjestedPValue);
            proteinItr = proteinMap.keySet().iterator();
            for(double value : adjestedPValue)
            {
                Protein currentProtein = proteinMap.get(proteinItr.next());
                currentProtein.addBHCorrectionList(value);
            }
        }
        
        System.out.println("BHCorrection Done.....");
    }
    
    private void generateBHCorrectionRatio()
    {
        int totalGroup = ((Protein)((Collection<Protein>)proteinMap.values()).toArray()[0]).getpValueRatioList().size();
        for(int i=0;i<totalGroup ; i++)
        {
           Iterator proteinItr = proteinMap.keySet().iterator();
           List<Double> adjestedPValue = new ArrayList<>();
           while (proteinItr.hasNext()) {
            Protein currentProtein = proteinMap.get(proteinItr.next());
            adjestedPValue.add(currentProtein.getpValueRatioList().get(i));
            }
            adjestedPValue = TmtfilterPeptideSimple.performBHCorrection(adjestedPValue);
            proteinItr = proteinMap.keySet().iterator();
            for(double value : adjestedPValue)
            {
                Protein currentProtein = proteinMap.get(proteinItr.next());
                currentProtein.addBHCorrectionRatioList(value);
            }
            
            
        }
        
        System.out.println("BHCorrection Done.....");
    }
    
    public void generatePValuesRatio()
    {
        Iterator proteinItr = proteinMap.keySet().iterator();
        while (proteinItr.hasNext()) 
         {
             Protein currentProtein = (Protein) proteinMap.get(proteinItr.next());
             Iterator<Peptide> peptideItr = currentProtein.getPeptideMap().values().iterator();
             Peptide currentPeptide = peptideItr.next();
             int dataSize = ((List)currentPeptide.ratioListForPvalue.get(0)).size();
             int totalGroupSize =  currentPeptide.ratioListForPvalue.size();
             for (int j = 0; j < totalGroupSize; j++)
             {
                double data[] = new double[dataSize];

                 peptideItr = currentProtein.getPeptideMap().values().iterator();
                 while (peptideItr.hasNext()) 
                 {

                     currentPeptide = peptideItr.next();

                     List<List<DescriptiveStatistics>> pValueRatioData = currentPeptide.ratioListForPvalue;
                     List<DescriptiveStatistics> currentGroup = pValueRatioData.get(j);
                     for (int i = 0; i < currentGroup.size(); i++)
                     {
                         data[i] += currentGroup.get(i).getMean();
                     }
                 }
                currentProtein.generatepValueRatioList(data, currentProtein.getPeptideMap().size());// to find average divide by totalpeptpide
                

             }
             
             
             
         }
             
             System.out.println("pvalue ration done..");
   /*      while (proteinItr.hasNext()) 
         {
             Protein currentProtein = (Protein) proteinMap.get(proteinItr.next());
             Iterator<Peptide> peptideItr = currentProtein.getPeptideMap().values().iterator();
            double data[]= new double[peptideItr.next().getRatioPvalueData().size()];
            peptideItr = currentProtein.getPeptideMap().values().iterator();
             while (peptideItr.hasNext())
             {
                    Peptide currentPeptide = peptideItr.next();
                    List<Double> pValueRatioData = currentPeptide.getRatioPvalueData();
                    for(int i =0;i<pValueRatioData.size();i++)
                    {
                        data[i]+=pValueRatioData.get(i);
                    }
             }
             currentProtein.generatepValueRatioList(data,currentProtein.getPeptideMap().size());// to find average divide by totalpeptpide
             
             
         }
             
             System.out.println("pvalue ration done..");
             */
             
//             while(peptideItr.hasNext())
//             {
//                 Peptide currentPeptide = peptideItr.next();
//                 currentPeptide.ratioListForPvalue.clear();
//                 List<List<Double>> intensityList = currentPeptide.getIntensityAverage();
//                 List<DescriptiveStatistics> statList = new ArrayList<>();
//                 for(int i =0 ;i<intensityList.size();i++)
//                 {
//                     List<Double> currentGroup = intensityList.get(i);
//                     for(int j= i +1; j<intensityList.size();j++)
//                     {
//                         List<Double> nextGroup = intensityList.get(i);
//                         DescriptiveStatistics stat = new DescriptiveStatistics();
//                         for(int z =0; z<currentGroup.size();z++)
//                         {
//                            double value = Math.log(currentGroup.get(z)/nextGroup.get(z))/Math.log(2);
//                            stat.addValue(value);
//                         }
//                         statList.add(stat);
//                     }
//                     currentPeptide.ratioListForPvalue.add(statList);
//                 }
//                 System.out.println("");
//             }
//             peptideItr = currentProtein.getPeptideMap().values().iterator();
//             ----------------------------------------------------
//             List<DescriptiveStatistics> proteinStat = new ArrayList<>();
//             int length = peptideItr.next().ratioListForPvalue.size();
//             peptideItr = currentProtein.getPeptideMap().values().iterator();
//            double[] ttestData  = new double[length];
//
//             for(int i=0;i<length;i++)
//             {
//                 peptideItr = currentProtein.getPeptideMap().values().iterator();
//                DescriptiveStatistics currentPeptideValue = new DescriptiveStatistics();
//
//                 while(peptideItr.hasNext())
//                 {
//                    Peptide currentPeptide = peptideItr.next();
//                    for(DescriptiveStatistics ls : (List<DescriptiveStatistics>)currentPeptide.ratioListForPvalue.get(i))
//                    {
//                        currentPeptideValue.addValue(ls.getMean());
//                    }
//                 }
//                 ttestData[i]=currentPeptideValue.getMean();
//             }
//              make approp[rita changes///..
//             
//             
             
             
             
    }
    
    
    
      public void generatePValues()
    {
        Iterator proteinItr = proteinMap.keySet().iterator();
         while (proteinItr.hasNext()) {
             Protein currentProtein = (Protein) proteinMap.get(proteinItr.next());
             currentProtein.generateProteinIntAvg();
//            Map<String, Peptide> peptideData = currentProtein.getPeptideMap();
//            Iterator peptideItr = peptideData.keySet().iterator();
//            List<Double> pValueList = new ArrayList<>();

            List<List<Double>> intensityAverage = currentProtein.getIntensityAvg();
            int totalGroup = intensityAverage.size();
            int minSizedGroupLength = intensityAverage.get(0).size();
            for (List currentGroup : intensityAverage) {
                if (currentGroup.size() < minSizedGroupLength) {
                    minSizedGroupLength = currentGroup.size();
                }
            }

           
                double pValueData[] = new double[minSizedGroupLength];
                for (int j = 0; j < totalGroup; j++) 
                {
                   List<Double> currentIntensity = intensityAverage.get(j);
                    for (int k = j + 1; k < totalGroup; k++) 
                    {
                        List<Double> nextIntensity = intensityAverage.get(k);

                        for (int l = 0; l < minSizedGroupLength; l++)
                            pValueData[l] = Math.log(currentIntensity.get(l) / nextIntensity.get(l)) / Math.log(2);;
                        try {
                            double d = TestUtils.tTest(0, pValueData);
                            if(Double.isNaN(d) )
                                d=1.0;
                            currentProtein.addPValueList(d);
                        } catch (IllegalArgumentException ex) {
                            Logger.getLogger(TmtfilterPeptideSimple.class.getName()).log(Level.SEVERE, null, ex);
                        }
                    }
                }
            
         }
    }
      
/*    
    public void generatePValues()
    {
        Iterator proteinItr = proteinMap.keySet().iterator();
         while (proteinItr.hasNext()) {
             Protein currentProtein = (Protein) proteinMap.get(proteinItr.next());
             currentProtein.generateProteinIntAvg();
            Map<String, Peptide> peptideData = currentProtein.getPeptideMap();
            Iterator peptideItr = peptideData.keySet().iterator();
            List<Double> pValueList = new ArrayList<>();

            while (peptideItr.hasNext()) {
                Peptide currentPeptide = peptideData.get(peptideItr.next());
                double value = currentPeptide.getpValue();
                
                pValueList.add(value);
            }
            if(pValueList.size() >0)
            {
//                System.out.println("----"+ currentProtein.getLocus() + " :"+ pValueList.toString());
                double pValue = edu.scripps.pms.util.stats.Fisher.combinedProbabilityTest(pValueList);
                currentProtein.setpValue(pValue);
            }
            
            
            
         }
    }
    */

}
