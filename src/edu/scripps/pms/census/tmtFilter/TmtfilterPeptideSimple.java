/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.scripps.pms.census.tmtFilter;

import java.math.BigDecimal;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.inference.TestUtils;
import rpark.statistics.AnovaUtil;

import static edu.scripps.pms.util.isotope.RatioIntensityUtils.checkRatio;

/**
 * Census2 : Isobaric Labelling Data Analysis in an Automatic way. Research
 * Paper for more details......
 *
 * @author Harshil
 */
public class TmtfilterPeptideSimple extends TmtfilterPeptide {

    public static void main(String args[]) throws Exception {
        TmtfilterPeptideSimple tmtfilterSimple = new TmtfilterPeptideSimple();
        Map<String, List> groupNameMap = new LinkedHashMap<>();
        List groupList = new ArrayList();
        groupList.add("128.128116");
        groupList.add("129.131471");
        groupList.add("130.134825");
        groupList.add("131.13818");

        groupNameMap.put("CTR", groupList);

        groupList = new ArrayList();
        groupList.add("127.131081");
        groupList.add("128.134436");
        groupList.add("129.13779");
        groupList.add("130.141145");

        groupNameMap.put("CDKL5", groupList);

//        groupList = new ArrayList();
//        groupList.add("130.141141");
//        groupList.add("131.138176"); 
//        groupNameMap.put("dell", groupList);
        
//          tmtfilterSimple.runTmtPeptideAnalysis(groupNameMap,"/home/rpark/test_data/meha_tmt", "new_compare_peptide.txt");
        String path = "E:\\";
//          String path = "E:\\";
//         String path = "/data/2/rpark/ip2_data//yrc/EvoProteomics/20140607_HR4_6_multinotchTMT_6hrUPLC_2014_06_09_16_25485/quant/2014_06_10_15_9811/";

//           String path = "/data/2/rpark/ip2_data//yrc/EvoProteomics/20140607_HR4_6_multinotchTMT_6hrUPLC_2014_06_09_16_25485/quant/2014_06_10_15_9811";
//        String path = "/home/rpark";
        //tmtfilterSimple.runTmtPeptideAnalysis(groupNameMap, path, path, "census-out-peptideAnalysisSimple.txt", false);
//        tmtfilterSimple.runTmtPeptideAnalysis(groupNameMap, "/data/2/rpark/ip2_data//xmhan/FAd_TMT/HGX_TMT6_20140725_Fusion_1_2014_07_29_09_26292/quant/2014_07_30_11_10034/", "/data/2/rpark/ip2_data//xmhan/FAd_TMT/HGX_TMT6_20140725_Fusion_1_2014_07_29_09_26292/quant/2014_07_30_11_10034/", "census-out_compare.txt" , false);
     /*   tmtfilterSimple.runTmtPeptideAnalysis(groupNameMap,"C:\\Users\\Harshil\\Desktop",
                "C:\\Users\\Harshil\\Desktop",
                "census-out-peptideAnalysisSimple.txt",true);*/
        tmtfilterSimple.runTmtPeptideAnalysis(groupNameMap,"/home/yateslab/project_data/census/1804tmtDebug",
                "/home/yateslab/project_data/census/1804tmtDebug",
                "census-out-peptideAnalysisSimple.txt",false);
    }

    @Override
    public void runTmtPeptideAnalysis(Map groupMap, String inputPath, String outputPath, String outPutFileName, boolean normalizeValue) throws Exception {
        isNormalized = normalizeValue;
        parseCensusFile(inputPath);
        groupNameToKeys(groupMap);
        generateGroupInfo();
        filterValues(0.05);
        generateRatio();
        generatePValues();
        generateBHCorrection();
        generatePValuesRatio();
        generateBHCorrectionRatio();

        printToFile(inputPath, outputPath, outPutFileName);
    }

    @Override
    public Map getPeptideMap(Map groupMap, String inputPath)  {
        parseCensusFile(inputPath);
        groupNameToKeys(groupMap);
        generateGroupInfo();
        filterValues(0.05);
        generateRatio();
        generatePValues();
        generateBHCorrection();
        generatePValuesRatio();
        generateBHCorrectionRatio();        
        return peptideData;
    }

    private void generateBHCorrectionRatio()
    {
     
        int totalGroup = ((Peptide)((Collection<Peptide>)peptideData.values()).toArray()[0]).getpValueRatioList().size();
        for(int i=0;i<totalGroup ; i++)
        {
           Iterator peptideItr = peptideData.keySet().iterator();
           List<Double> adjustedPValue = new ArrayList<>();
           List<Double> adjustedNormPvalue = new ArrayList<>();
           while (peptideItr.hasNext()) {
            Peptide currentPeptide = peptideData.get(peptideItr.next());
            adjustedPValue.add(currentPeptide.getpValueRatioList().get(i));
            adjustedNormPvalue.add(currentPeptide.getpValueNormRatioList().get(i));
            }
            adjustedPValue = performBHCorrection(adjustedPValue);
            adjustedNormPvalue = performBHCorrection(adjustedNormPvalue);
            peptideItr = peptideData.keySet().iterator();
            for(double value : adjustedPValue)
            {
                Peptide currentPeptide = peptideData.get(peptideItr.next());
                currentPeptide.addBHCorrectionRatio(value);
            }
            peptideItr = peptideData.keySet().iterator();
            for(double value : adjustedNormPvalue)
            {
                Peptide currentPeptide = peptideData.get(peptideItr.next());
                currentPeptide.addBHCorrectionNormRatio(value);
            }
        }
        
        System.out.println("BHCorrection-Ratio Done.....");
    }
    private void generatePValuesRatio()  
   {
        Iterator peptideItr = peptideData.keySet().iterator();
        List groupNameList = new ArrayList(groupKeyMap.keySet());
        while (peptideItr.hasNext()) {
            String currentPeptideName = (String) peptideItr.next();
            Peptide currentPeptide = peptideData.get(currentPeptideName);
            Map<String, List> allIntensity = currentPeptide.getallIntensity();
            Map<String, List<Double>> normAllIntensity = currentPeptide.getNormAllIntensity();
           // if(currentPeptide.getSequence().equals("K.TAVNALWGK.V"))
         //       System.out.println("");
            List statList = new ArrayList();
            for (int i = 0; i < groupNameList.size() - 1; i++) 
            {
                List tempStat = new ArrayList();
                List<DescriptiveStatistics> normTempStat = new ArrayList<>();
                List currrentGroupKey = groupKeyMap.get(groupNameList.get(i));
                for (int j = i + 1; j < groupNameList.size(); j++) 
                {
                    List nextGroupKey = groupKeyMap.get(groupNameList.get(j));
                    tempStat = new ArrayList();
                    normTempStat = new ArrayList<>();
                    for (int z = 0; z < currrentGroupKey.size(); z++) 
                    {
                        DescriptiveStatistics stat = new DescriptiveStatistics();
                        DescriptiveStatistics normStat = new DescriptiveStatistics();
                        for (int t1 = 0; t1 < allIntensity.get(normintensity_name.get(0)).size(); t1++) 
                        {
                            double value1 = (double) ((List) allIntensity.get(currrentGroupKey.get(z))).get(t1);
                            double value2 = (double) ((List) allIntensity.get(nextGroupKey.get(z))).get(t1);
                            double ratio = value1/value2;

                            double normValue1 = normAllIntensity.get(currrentGroupKey.get(z)).get(t1);
                            double normValue2 = normAllIntensity.get(nextGroupKey.get(z)).get(t1);
                            double normRatio = normValue1/normValue2;

                            ratio = checkRatio(ratio);

                            if(!Double.isNaN(ratio))
                            {
                                double logRatio = Math.log(ratio) / Math.log(2);
                                stat.addValue(logRatio);
                            }


                            normRatio = checkRatio(normRatio);

                            if(!Double.isNaN(normRatio))
                            {
                                double normLogRatio = Math.log(normRatio)/ Math.log(2);
                                normStat.addValue(normLogRatio);
                            }



//                                    double logRatio = value1/value2;
//                                  stat.addValue(((float)((List)allIntensity.get(currrentGroupKey.get(t1))).get(z)) / ((float)((List)allIntensity.get(nextGroupKey.get(t2))).get(z)) );
                            //System.out.print(value1 +" " + value2 +" "+logRatio+" ");
                        }
                        tempStat.add(stat);
                        normTempStat.add(normStat);
                    }
                    currentPeptide.ratioListForPvalue.add(tempStat);
                    currentPeptide.normRatioListForPValue.add(normTempStat);
                }
            }
            for (List<DescriptiveStatistics> currentGroup : (List<List<DescriptiveStatistics>>) currentPeptide.ratioListForPvalue)
            {
                double tTestData[] = new double[currentGroup.size()];
                DescriptiveStatistics tempStat = new DescriptiveStatistics();
                for (int i = 0; i < currentGroup.size(); i++)
                {
                    DescriptiveStatistics stat = (DescriptiveStatistics) currentGroup.get(i);
                    tTestData[i] = stat.getSum();
                    tempStat.addValue(stat.getMean());
//                    currentPeptide.getRatioPvalueData().add(stat.getMean());
                }

                double x= tempStat.getStandardDeviation()/tempStat.getMean()*100;
               // System.out.println(">>ratio calc "+x+" "+tempStat.getStandardDeviation());

                currentPeptide.addRsdPaired(tempStat.getStandardDeviation()/tempStat.getMean()*100);
                currentPeptide.addstdevPaired(tempStat.getStandardDeviation());

                try {
                    Double d = TestUtils.tTest(0, tTestData);
                    if(Double.isNaN(d))
                        d=1.0;
                    currentPeptide.getpValueRatioList().add(d);
                } catch (Exception ex) {
                    Logger.getLogger(TmtfilterPeptideSimple.class.getName()).log(Level.SEVERE, null, ex);
                }
            }

            for (List<DescriptiveStatistics> normCurrentGroup :  currentPeptide.normRatioListForPValue)
            {
                double tTestData[] = new double[normCurrentGroup.size()];
                DescriptiveStatistics normTempStat = new DescriptiveStatistics();
                for (int i = 0; i < normCurrentGroup.size(); i++)
                {
                    DescriptiveStatistics stat = (DescriptiveStatistics) normCurrentGroup.get(i);
                    tTestData[i] = stat.getSum();
                    normTempStat.addValue(stat.getMean());
//                    currentPeptide.getRatioPvalueData().add(stat.getMean());
                }

                double x= normTempStat.getStandardDeviation()/normTempStat.getMean()*100;
               // System.out.println(">>norm ratio calc "+x+" "+normTempStat.getStandardDeviation());

                currentPeptide.addNormRsdPaired(normTempStat.getStandardDeviation()/normTempStat.getMean()*100);
                currentPeptide.addNormstdevPaired(normTempStat.getStandardDeviation());

                try {
                    Double d = TestUtils.tTest(0, tTestData);
                    if(Double.isNaN(d))
                        d=1.0;
                    currentPeptide.getpValueNormRatioList().add(d);
                } catch (Exception ex) {
                    Logger.getLogger(TmtfilterPeptideSimple.class.getName()).log(Level.SEVERE, null, ex);
                }
            }

//            System.out.println("");
        }
    }
    
    public void generatePValues()  {
        Iterator peptideItr = peptideData.keySet().iterator();
        List groupNameList = new ArrayList(groupKeyMap.keySet());

        while (peptideItr.hasNext()) {
            Peptide currentPeptide = peptideData.get(peptideItr.next());
            List<Double> pValueList = new ArrayList<>();
        /*     if(currentPeptide.getSequence().equals("K.TAVNALWGK.V"))
                System.out.println("");*/
            List<List<Double>> intensitySum = currentPeptide.getIntensityAverage();
            List<List<Double>> normIntensitySum = currentPeptide.getNormIntensitySum();
            int totalGroup = intensitySum.size();
            int minSizedGroupLength = intensitySum.get(0).size();
            for (List currentGroup : intensitySum) {
                if (currentGroup.size() < minSizedGroupLength) {
                    minSizedGroupLength = currentGroup.size();
                }
            }
            double pValueData[] = new double[minSizedGroupLength];
            for (int j = 0; j < totalGroup; j++) {
                List<Double> currentIntensity = intensitySum.get(j);
                List<Double> currentNormIntensity = normIntensitySum.get(j);
                double [] ref  = convertToArray(currentIntensity);
                double [] normRef = convertToArray(currentNormIntensity);
                for (int k = j + 1; k < totalGroup; k++) {
                    List<Double> nextIntensity = intensitySum.get(k);
                    List<Double> normIntensity = normIntensitySum.get(k);

                    double [] next  = convertToArray(nextIntensity);
                    double [] normNext = convertToArray(normIntensity);

                    List<double []> classes = new ArrayList<>();
                    List<double []> normClasses = new ArrayList<>();

                    normClasses.add(normRef);
                    normClasses.add(normNext);

                    classes.add(ref);
                    classes.add(next);
                    double pvalue = Double.NaN;
                    double normPValue = Double.NaN;
                    try {
                        pvalue = AnovaUtil.calculateAnovaPvalue(classes);
                        normPValue = AnovaUtil.calculateAnovaPvalue(normClasses);
                    }
                    catch(Exception exc)
                    {

                    }


                    for (int l = 0; l < minSizedGroupLength; l++) {
                        pValueData[l] = Math.log(currentIntensity.get(l) / nextIntensity.get(l)) / Math.log(2);
                    };
                    try {

                        if(Double.isNaN(pvalue))
                        {
                            currentPeptide.addPValueList(1.0);

                        }
                        else
                        {
                            currentPeptide.addPValueList(pvalue);
                        }

                        if(Double.isNaN(normPValue))
                        {
                            currentPeptide.addNormPValueList(1.0);

                        }
                        else
                        {
                            currentPeptide.addNormPValueList(normPValue);
                        }



                        /*
                        if(Double.isNaN(TestUtils.tTest(0, pValueData)))
                            currentPeptide.addPValueList(1.0);
                        else
                            currentPeptide.addPValueList(TestUtils.tTest(0, pValueData));
*/
                    } 
                    catch (IllegalArgumentException ex) {
                        Logger.getLogger(TmtfilterPeptideSimple.class.getName()).log(Level.SEVERE, null, ex);
                    } 
                   
                }
            }
        }

    }
    public static double [] convertToArray(List<Double> list)
    {
        double[] result = new double[list.size()];
        for(int i=0; i<result.length; i++)
        {
            result[i]=list.get(i);
        }
        return result;
    }
    
    private void generateBHCorrection()
    {
        int totalGroup = ((Peptide)((Collection<Peptide>)peptideData.values()).toArray()[0]).getpValueList().size();
        for(int i=0;i<totalGroup ; i++)
        {
           Iterator peptideItr = peptideData.keySet().iterator();
           List<Double> adjustedPValue = new ArrayList<>();
           List<Double> adjustedNormPValue = new ArrayList<>();
           while (peptideItr.hasNext()) {
            Peptide currentPeptide = peptideData.get(peptideItr.next());
            adjustedPValue.add(currentPeptide.getpValueList().get(i));
            adjustedNormPValue.add(currentPeptide.getNormPValueList().get(i));
            }
            adjustedPValue = performBHCorrection(adjustedPValue);
           adjustedNormPValue =performBHCorrection(adjustedNormPValue);
            peptideItr = peptideData.keySet().iterator();
            for(double value : adjustedPValue)
            {
                Peptide currentPeptide = peptideData.get(peptideItr.next());
                currentPeptide.addBHCorrection(value);
            }
            peptideItr = peptideData.keySet().iterator();
            for(double value : adjustedNormPValue)
            {
                Peptide currentPeptide = peptideData.get(peptideItr.next());
                currentPeptide.addBHCorrectionNorm(value);
            }
        }
        
        System.out.println("BHCorrection Done.....");
    }
    
    protected static List<Double> performBHCorrection(List pValues) {
        int len = pValues.size();
        double[] orderedPValues = new double[len];
        double[] adjustedpValues = new double[len];
        List<Double> BhCorrection = new ArrayList<>();
        int[] indexOfValues = new int[len];
        final int RESULT_SCALE = 10;

        for (int i = 0; i < len; i++) {
            orderedPValues[i] = (double) pValues.get(i);
        }
        // sort the values
        java.util.Arrays.sort(orderedPValues);
        for (int i = 0; i < len; i++) {
            indexOfValues[i] = getIndexOf(orderedPValues, (double) pValues.get(i));
        }
        // calculate the post hoc adjustment
        BigDecimal min = new BigDecimal("" + 1);
        BigDecimal mkprk;
//        System.out.println("");
        for (int i = len; i > 0; i--) {
            mkprk = (new BigDecimal("" + len).multiply(new BigDecimal(orderedPValues[i - 1]))).divide(new BigDecimal("" + i), RESULT_SCALE, BigDecimal.ROUND_HALF_UP);
            if (mkprk.compareTo(min) < 0) {
                min = mkprk;
            }
            adjustedpValues[i - 1] = min.doubleValue();
        }
        // adjust the sequence
        len = pValues.size();
        int j = 0;
        for (int i = 0; i < len; i++) {
            try {
                double tmp = (double) pValues.get(i);

                if (tmp > -100000.0) {
                    NumberFormat formatter = new DecimalFormat("##.#####");
                    String apvString = formatter.format(adjustedpValues[indexOfValues[j]]);
                    //                      proteinList.get(i).setPostHocP(apvString);
                    BhCorrection.add(adjustedpValues[indexOfValues[j]]);
                    j++;
                } else {
                    //                      proteinList.get(i).setPostHocP("-1");
                    BhCorrection.add(-1.0);
                }
            } catch (Exception e) {
                BhCorrection.add(-1.0);
            }
        }
        return BhCorrection;
    }

    private static int getIndexOf(double[] orderedPvalues, double value) {
        int index = -1;
        for (int i = 0; i < orderedPvalues.length; i++) {
            if (orderedPvalues[i] == value) {
                return i;
            }
        }
        return index;
    }

    /*
     public void generatePValues() {
     Iterator peptideItr = peptideData.keySet().iterator();
     List groupNameList = new ArrayList(groupKeyMap.keySet());

     while (peptideItr.hasNext()) {
     Peptide currentPeptide = peptideData.get(peptideItr.next());
     List<Double> pValueList = new ArrayList<>();
     //            System.out.println(currentPeptide.getSequence());
     Map<String, List> allIntensity = currentPeptide.getallIntensity();
     for (int i = 0; i < allIntensity.get(normintensity_name.get(0)).size(); i++) {
     List tTestList = new ArrayList();
     for (int j = 0; j < groupNameList.size(); j++) {
     List currrentGroupKey = groupKeyMap.get(groupNameList.get(j));
     double[] group = new double[currrentGroupKey.size()];
     for (int t1 = 0; t1 < currrentGroupKey.size(); t1++) {
     group[t1] = (double) allIntensity.get(currrentGroupKey.get(t1)).get(i);

     }
     tTestList.add(group);
                    
                    

     }
     double pValue = 1.0;
     try {
     pValue=TestUtils.oneWayAnovaPValue(tTestList);
     } catch (IllegalArgumentException ex) {
     Logger.getLogger(TmtfilterPeptideSimple.class.getName()).log(Level.SEVERE, null, ex);
     } catch (MathException ex) {
     Logger.getLogger(TmtfilterPeptideSimple.class.getName()).log(Level.SEVERE, null, ex);
     }
                
     pValueList.add(pValue);
                

     if (pValueList.size() <= 0) {
     continue;
     }

     pValue = combinedProbabilityTest(pValueList);
     //Logger.getLogger(TmtfilterPeptideSimple.class.getName()).log(Level.SEVERE, null, ex);

     if (pValue >= 1.0) {
     continue;
     }

     if (pValue <= 0.0) {
     pValue = 0.0000001;
     }

     pValueList.add(pValue);
     }

     double pValue = 1.0d;

     if (pValueList.size() > 0) {
     pValue = combinedProbabilityTest(pValueList);
     }
     //else
     //    pValue = pValueList.get(0);
     if (pValue <= 0.0) {
     pValue = 0.0000001;
     }

     currentPeptide.setpValue(pValue);
     currentPeptide.setpValueList(pValueList);
     //            System.out.println("-->Pvalue-"+pValue+"-->pValueList-"+pValueList.toString());                    
     }
     }
     */
}
