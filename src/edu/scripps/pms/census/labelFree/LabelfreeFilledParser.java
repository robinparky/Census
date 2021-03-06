/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.scripps.pms.census.labelFree;

import edu.scripps.pms.census.conf.Configuration;
import edu.scripps.pms.census.io.FastaFileProteinLengthCalc;
import edu.scripps.pms.census.model.ChroPeptide;
import edu.scripps.pms.census.model.ChroProtein;
import edu.scripps.pms.census.model.SampleModel;
import edu.scripps.pms.census.util.io.FastaReader;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import rpark.statistics.AnovaUtil;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.math.BigDecimal;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Created by rpark on 4/15/16.
 * modified by rohan rampuria , Titus Jung
 */
public class LabelfreeFilledParser {

    private BufferedReader br = null;
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
    private List<Integer> corrInjectionIntensityIndexList = new ArrayList<>();
    private List<Integer> scountIndexList = new ArrayList<Integer>();
    private List<Integer> normIntensityIndexList = new ArrayList<Integer>();
    private List<Integer> pepLogratioIndexList = new ArrayList<>();
    private List<Integer> bestintpepIndexList = new ArrayList<>();
    private List<Integer> bestintIndexList = new ArrayList<>();
    private List<Integer> bestintpepIITIndexList = new ArrayList<>();
    private List<Integer> bestintIITIndexList = new ArrayList<>();
    private List<Integer> pepmedianlogratioIndexList = new ArrayList<>();
    private List<Integer> groupIntensitySumIndexList = new ArrayList<>();
    private List<Integer> bestPepIonInjectionIndexList = new ArrayList<>();
    private List<Integer> bestPepGroupIonInjectionIndexList = new ArrayList<>();
    private List<Integer> totolPeptideIntensityIndex = new ArrayList<>();
    private int totolPeptideIntensityDividedByProteinIndex = -1;
    private int proteinLengthIndex = -1;

    private List<Integer> normPeptideIntensityIndexList = new ArrayList<>();

    private int accessionIndex = -1;
    private int descriptionIndex  =-1;
    private int pepCountIndex  =-1;
    private int[] AvgNormIntensityIndex;
    private final String filledFile;

    public static void main(String[] args) throws Exception{

        //String fname = "/data/2/rpark/ip2_data/xudong/TMTMS3/labelfree_quant/labelfree_12334/census_labelfree_out_12334_filled.txt.small";


        //String fname = "/data/2/rpark/ip2_data/xudong/TMTMS3/labelfree_quant/labelfree_12334/census_labelfree_out_12334_filled.txt";
  /*     String fname = "/home/rampuria/data/labelfree_12575/census_labelfree_out_12575_filled.txt";
       //String fname = "/home/rampuria/data/labelfree_11912/census_labelfree_out_11912_filled.txt";
       // String confFile = "/home/rampuria/data/labelfree_11912/census_config_labelfree_11912.xml";
        String confFile = "/home/rampuria/data/labelfree_12575/census_config_labelfree_12575.xml";
        List<Integer> sample1 = new ArrayList<>();
        sample1.add(0);
        sample1.add(1);
        sample1.add(2);
        List<Integer> sample2 = new ArrayList<>();
        sample2.add(3);
        sample2.add(4);
        sample2.add(5);

        LabelfreeFilledParser l = new LabelfreeFilledParser(fname);
        List<ProteinModel> list = l.readWholeFile();
        LabelfreeMissingPeptideBuilder.proteinCompareBasedOnRatios(list,sample1,sample2);
        LabelfreeMissingPeptideBuilder.generateLabelfreeOutputFile(list,"/home/rampuria/data/labelfree_12575/newstat.txt",confFile);

*/
       // String fname = "/data/2/rpark/ip2_data/yrc/Dan_Salomon_Saby/labelfree_quant/labelfree_13259/census_labelfree_out_13259_filled.txt";
        //String fname = "/data/2/rpark/ip2_data/rpark/an_chi/labelfree_quant/labelfree_13378/census_labelfree_out_13378_filled.txt";
       // String confFile = "/data/2/rpark/ip2_data/rpark/an_chi/labelfree_quant/labelfree_13378/census_config_labelfree_13378.xml";
       // String confFile = "/data/2/rpark/ip2_data/yrc/Dan_Salomon_Saby/labelfree_quant/labelfree_13259/census_config_labelfree_13259.xml";

        String confFile = args[0];
       String fname = args[1];
       String fastaFile = args[2];

        //String confFile = "/data/2/rpark/ip2_data/drocha/Foc_Brazil/labelfree_quant/labelfree_13521/test/census_config_labelfree_13521.xml";
        //String fname = "/data/2/rpark/ip2_data/drocha/Foc_Brazil/labelfree_quant/labelfree_13521/test/census_labelfree_out_13521_filled.txt";

        /*
        Configuration conf = Configuration.getInstance();
        if (!conf.isReadConfigFile()) {
            conf.setLabelfree(true);
            conf.readXMLParam(confFile);
        }
        List<org.jdom.Element> samGroupEleList = conf.getRootConfEle().getChildren("sample");
        List<Integer> sample1IndexList = new ArrayList<Integer>();
        List<Integer> sample2IndexList = new ArrayList<Integer>();
        int count = 0;

        List<org.jdom.Element> sampleEleList = samGroupEleList.get(0).getChildren("each_sample");

        for (org.jdom.Element eachSam : sampleEleList) {
            sample1IndexList.add(count);
            count++;
        }

        sampleEleList = samGroupEleList.get(1).getChildren("each_sample");

        for (org.jdom.Element eachSam : sampleEleList) {
            sample2IndexList.add(count);
            count++;
        }
        */
        //sample2.add(3);

        LabelfreeFilledParser l = new LabelfreeFilledParser(fname);
        List<ProteinModel> list = l.readWholeFile();
        processProteinList(list,confFile,fastaFile);
        LabelfreeMissingPeptideBuilderSplit.generateLabelfreeOutputFile(list, fname.substring(0, fname.length() - 11) + "_stat.txt",confFile);

    }
    public LabelfreeFilledParser(String fileName) {
        this.filledFile = fileName;
        try {
            br = new BufferedReader(new FileReader(filledFile));
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }
    public List<ProteinModel> readWholeFile() {

        List<ProteinModel> proteinList = new ArrayList<>();
        List<Integer> sampleCount = new ArrayList<>();
        try
        {
            String eachLine = null;
            ChroProtein chroProtein = null;
            ProteinModel proteinModel = null;
            List<ChroPeptide> chroPepList;


            boolean isNewProtein=true;


            while( (eachLine = br.readLine()) != null)
            {
                if(eachLine.startsWith("H\t")){
                    String [] words = eachLine.split("\t");
                    int j=0;
                    if(sampleCount.size() != 0) continue;
                    for(int i =3;i<words.length;i++){
                        sampleCount.add(j);
                        j++;
                    }
                    System.out.println("");
                }

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
                        String [] words = eachLine.split("\t");
                        if (bestintpepIndexList.size() != 0) {
                            for (int index : bestintpepIndexList) {
                                if (words[index].equalsIgnoreCase("NA") || words[index].equalsIgnoreCase("NaN")) {
                                    proteinModel.addBestCorrelationEachpeptideIntensityList(-1);
                                } else {
                                    long l = new BigDecimal(words[index]).longValue();
                                    proteinModel.addBestCorrelationEachpeptideIntensityList(l);
                                }
                            }
                        }
                        if (bestintIndexList.size() != 0) {
                            for (int index : bestintIndexList) {
                                if (words[index].equalsIgnoreCase("NA") || words[index].equalsIgnoreCase("NaN")) {
                                    proteinModel.addBestCorrelationIntensityList(-1);
                                } else {
                                    proteinModel.addBestCorrelationIntensityList(Long.parseLong(words[index]));
                                }
                            }
                        }
                        if (bestintpepIITIndexList.size() != 0) {
                            for (int index : bestintpepIITIndexList) {
                                    if (words[index].equalsIgnoreCase("NA") || words[index].equalsIgnoreCase("NaN")) {
                                        proteinModel.addBestCorrelationEachpeptideIntensityListIIT(-1);
                                    } else {
                                        long lng = new BigDecimal(words[index]).longValue();
                                        proteinModel.addBestCorrelationEachpeptideIntensityListIIT(lng);
                                    }

                            }
                        }

                        if (bestintIITIndexList.size() != 0) {
                            for (int index : bestintIITIndexList) {
                                if(index <words.length)
                                {
                                    if (words[index].equalsIgnoreCase("NA") || words[index].equalsIgnoreCase("NaN")) {
                                        proteinModel.addBestCorrelationIntensityListIIT(-1);
                                    } else {
                                        long in = new BigDecimal(words[index]).longValue();
                                        proteinModel.addBestCorrelationIntensityListIIT(in );
                                    }
                                }

                            }
                        }

                        if (normIntensityIndexList.size() != 0) {
                            for (int index : normIntensityIndexList) {
                                if (words[index].equalsIgnoreCase("NA") || words[index].equalsIgnoreCase("NaN")) {
                                    proteinModel.addNormIntensityList(-1.0);
                                } else {
                                    proteinModel.addNormIntensityList(Double.parseDouble(words[index]));
                                }
                            }
                        }
                        if (pepmedianlogratioIndexList.size() != 0) {
                            double [] peplogArr = new double[pepmedianlogratioIndexList.size()];
                            int k=0;
                            for (int index : pepmedianlogratioIndexList) {
                                if (words[index].equalsIgnoreCase("NA") || words[index].equalsIgnoreCase("NaN")) {
                                    peplogArr[k]=Double.NaN;
                                    k++;
                                } else {
                                    if(words[index].equals("Infinity")){
                                        peplogArr[k]=Double.POSITIVE_INFINITY;
                                        k++;
                                    }else{
                                        peplogArr[k]= Double.parseDouble(words[index]);
                                        k++;
                                    }

                                }
                            }
                            proteinModel.setPeptideMedianLogRatioArr(peplogArr);
                        }
                        if(groupIntensitySumIndexList.size() != 0)
                        {
                            double [] intensitySumArr = new double[groupIntensitySumIndexList.size()];
                            int k =0;
                            for (int index : groupIntensitySumIndexList) {
                                if (words[index].equalsIgnoreCase("NA") || words[index].equalsIgnoreCase("NaN")) {

                                    intensitySumArr[k++] =-1;
                                } else {
                                    long entry = new BigDecimal(words[index]).longValue();
                                    intensitySumArr[k++] = entry;
                                }
                            }
                            proteinModel.setIntensitySumArr(intensitySumArr);
                        }
                        if (bestPepIonInjectionIndexList.size() != 0) {

                            for (int index : bestPepIonInjectionIndexList) {
                                if (words[index].equalsIgnoreCase("NA") || words[index].equalsIgnoreCase("NaN")) {
                                    proteinModel.addBestCorrelationEachpeptideIntensityListIIT(-1);
                                } else {
                                    long entry = new BigDecimal(words[index]).longValue();
                                    proteinModel.addBestCorrelationEachpeptideIntensityListIIT(entry);
                                }
                            }
                        }
                        if (bestPepGroupIonInjectionIndexList.size() != 0) {
                            //long [] pepGroupIonInjectionArr = new long[bestPepGroupIonInjectionIndexList.size()];
                            //int k =0;
                            for (int index : bestPepGroupIonInjectionIndexList) {
                                if (words[index].equalsIgnoreCase("NA") || words[index].equalsIgnoreCase("NaN")) {
                                    proteinModel.addBestCorrelationIntensityListIIT(-1);
                                    //pepGroupIonInjectionArr[k++] =-1;
                                } else {
                                    long entry = new BigDecimal(words[index]).longValue();
                                    proteinModel.addBestCorrelationIntensityListIIT(entry);

                                    // pepGroupIonInjectionArr[k++] = entry;
                                }
                            }
                            //proteinModel.setIntensitySumArr(pepGroupIonInjectionArr);
                        }


                        proteinList.add(proteinModel);
                        isNewProtein = false;
                    }

                } else if(eachLine.startsWith("S\t")) {
                    isNewProtein = true;

                    chroPepList = parsePeptideLine(eachLine);
                    proteinModel.addPeptideList(chroPepList);

                }

            }



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
        if(pepLogratioIndexList.size()!= 0)
        {
            for (int index : pepLogratioIndexList)
            {
                if(words[index].equalsIgnoreCase("NA") || words[index].equalsIgnoreCase("NaN"))
                    protein.addPepRatioList(-1);
                else
                    protein.addPepRatioList(Double.parseDouble(words[index]));
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
            else if (words[i].contains("PEP_LOG_RATIO_"))
                pepLogratioIndexList.add(i);
            else if (words[i].startsWith("INTENSITY_"))
                bestintpepIndexList.add(i);
            else if (words[i].startsWith("AVG_GROUP_INTENSITY_VALUE_"))
                bestintIndexList.add(i);
            else if (words[i].startsWith("AVG_ION_INJECTION_TIME_INTENSITY_"))
                bestintpepIITIndexList.add(i);
            else if (words[i].startsWith("AVG_GROUP_ION_INJECTION_TIME_INTENSITY_VALUE_"))
                bestintIITIndexList.add(i);
            else if (words[i].startsWith("NORM_INTENSITY_"))
                normIntensityIndexList.add(i);
            else if (words[i].startsWith("MEDIAN_LOG_RATIO_"))
                pepmedianlogratioIndexList.add(i);
            else if (words[i].startsWith("GROUP_INTENSITY_SUM_"))
                groupIntensitySumIndexList.add(i);
            else if (words[i].startsWith("BEST_PEP_ION_INJECTION_TIME_INTENSITY_"))
                bestPepIonInjectionIndexList.add(i);
            else if (words[i].startsWith("BEST_PEP_GROUP_ION_INJECTION_TIME_INTENSITY_VALUE_"))
                bestPepGroupIonInjectionIndexList.add(i);
        }

    //    System.out.println("aaa");
    }


    public List<ChroPeptide> parsePeptideLine(String currentLine) {
        List<ChroPeptide> cPeptideList = new ArrayList<>();
        String words[] = currentLine.split("\t");
       // int wordCounter = 0;

        for (int i = 0; i < totalExperiments; i++) {

            ChroPeptide currentPeptide = new ChroPeptide();
            if (sequenceIndexList.size() > 0) {

                    currentPeptide.setSequence(words[sequenceIndexList.get(i)]);

            }
            if (intensityIndexList.size() > 0) {
                currentPeptide.setAverageIntensity(Double.parseDouble(words[intensityIndexList.get(i)]));
            }
            if (fileNameIndexList.size() > 0) {
                 if (words[fileNameIndexList.get(i)].equalsIgnoreCase("NA")) {
                    cPeptideList.add(currentPeptide);
                   // continue;
                 }
                 else{
                currentPeptide.setFileName(words[fileNameIndexList.get(i)]);
                 }
//                if ( !words[fileNameIndexList.get(i)].equalsIgnoreCase("na")) {
//                    msFileList.set(i,words[fileNameIndexList.get(i)]);
//                }
            }


            if (normPeptideIntensityIndexList.size() > 0) {
                if (words[normPeptideIntensityIndexList.get(i)].equalsIgnoreCase("NA")) {
                    currentPeptide.setNormIntensity(0.0);
                }
                else{
                    currentPeptide.setNormIntensity( Double.parseDouble(words[normPeptideIntensityIndexList.get(i)]) );
                }
            }




            if (scanIndexList.size() > 0) {
                if(!"NA".equals(words[scanIndexList.get(i)]) && "NaN".equals(words[scanIndexList.get(i)]))
                    currentPeptide.setScanNum(Integer.parseInt(words[scanIndexList.get(i)]));
            }
            if (csIndexList.size() > 0) {
                if(!"NA".equals(words[csIndexList.get(i)]) && "NaN".equals(words[csIndexList.get(i)]))
                    currentPeptide.setChargeState(words[csIndexList.get(i)]);
            }

            if (profileScoreIndexList.size() > 0) {
                if(!"NA".equals(words[profileScoreIndexList.get(i)]) && "NaN".equals(words[profileScoreIndexList.get(i)]))
                    currentPeptide.setAnCompositeScore(Double.parseDouble(words[profileScoreIndexList.get(i)]));
            }
            if (mhPlusIndexList.size() > 0) {
                if(!"NA".equals(words[mhPlusIndexList.get(i)]) && "NaN".equals(words[mhPlusIndexList.get(i)]))
                    currentPeptide.setMhPlus(words[mhPlusIndexList.get(i)]);
            }
            if (calcMHPlusIndexList.size() > 0) {
                if(!"NA".equals(words[calcMHPlusIndexList.get(i)]) && "NaN".equals(words[calcMHPlusIndexList.get(i)]))
                    currentPeptide.setCalcMHplus(words[calcMHPlusIndexList.get(i)]);
            }
            if (totalIntensityIndexList.size() > 0) {
                if(!"NA".equals(words[totalIntensityIndexList.get(i)]) && "NaN".equals(words[totalIntensityIndexList.get(i)]))
                    currentPeptide.setTotalIntensity(Double.parseDouble(words[totalIntensityIndexList.get(i)]));
            }
            if (xCorIndexList.size() > 0) {
                if(!"NA".equals(words[xCorIndexList.get(i)]) && "NaN".equals(words[xCorIndexList.get(i)]))
                    currentPeptide.setXCorr(words[xCorIndexList.get(i)]);
            }
            if (dcnIndexList.size() > 0) {
                if(!"NA".equals(words[dcnIndexList.get(i)]) && "NaN".equals(words[dcnIndexList.get(i)]))
                    currentPeptide.setDeltCN(words[dcnIndexList.get(i)]);
            }
            if (dMassIndexList.size() > 0) {
                if(!"NA".equals(words[dMassIndexList.get(i)]) && "NaN".equals(words[dMassIndexList.get(i)]))
                    currentPeptide.setDeltMass(words[dMassIndexList.get(i)]);
            }
            if (sprankIndexList.size() > 0) {
                if(!"NA".equals(words[sprankIndexList.get(i)]) && "NaN".equals(words[sprankIndexList.get(i)]))
                    currentPeptide.setSpRank(words[sprankIndexList.get(i)]);
            }
            if (spScoreIndexList.size() > 0) {
                if(!"NA".equals(words[spScoreIndexList.get(i)]) && "NaN".equals(words[spScoreIndexList.get(i)]))
                    currentPeptide.setSpScore(words[spScoreIndexList.get(i)]);
            }
            if (redundancyIndexList.size() > 0) {
                if(!"NA".equals(words[redundancyIndexList.get(i)]) && "NaN".equals(words[redundancyIndexList.get(i)]))
                    currentPeptide.setRedundancy(words[redundancyIndexList.get(i)]);
            }
            if (startIndexList.size() > 0) {
                if(!"NA".equals(words[startIndexList.get(i)]) && "NaN".equals(words[startIndexList.get(i)]))
                    currentPeptide.setStartRange(words[startIndexList.get(i)]);
            }
            if (endIndexList.size() > 0) {
                if(!"NA".equals(words[endIndexList.get(i)]) && "NaN".equals(words[endIndexList.get(i)]))
                    currentPeptide.setEndRange(words[endIndexList.get(i)]);
            }
            if (retentionIndexList.size() > 0) {
                if(!"NA".equals(words[retentionIndexList.get(i)]) && "NaN".equals(words[retentionIndexList.get(i)]))
                    currentPeptide.setRetentionTime(Double.parseDouble(words[retentionIndexList.get(i)]));
            }
            if (ionInjectionIndexList.size() > 0) {
                if(!"NA".equals(words[ionInjectionIndexList.get(i)]) && "NaN".equals(words[ionInjectionIndexList.get(i)]))
                    currentPeptide.setIonInjectionTime(Double.parseDouble(words[ionInjectionIndexList.get(i)]));
            }
            if (corrInjectionIntensityIndexList.size() > 0) {
                if(!"NA".equals(words[corrInjectionIntensityIndexList.get(i)]) && "NaN".equals(words[corrInjectionIntensityIndexList.get(i)]))
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
	        }else if (words[i].equalsIgnoreCase("CORRIONINJECTION_INTENISTY") || words[i].equalsIgnoreCase("CORRIONINJECTION_INTENSITY")) {
                corrInjectionIntensityIndexList.add(i);
            }else if (words[i].startsWith("NORM_INTENSITY"))
                normPeptideIntensityIndexList.add(i);




        }
    }

    public static void processProteinList(List<ProteinModel> proteinList, String configFile, String fastaFile) throws Exception
    {


        Iterator f = proteinList.iterator();
        List<Double> pvalueList = new ArrayList<>();
        Configuration conf = Configuration.getInstance();

        if (!conf.isReadConfigFile()) {
            conf.setLabelfree(true);
            conf.softReadXMLParam(configFile);
        }

        List<org.jdom.Element> samGroupEleList = conf.getRootConfEle().getChildren("sample");
        ArrayList sampleListt = new ArrayList();
        int count = samGroupEleList.size();
        int ctr = 0;
        for (int i = 0; i < count; i++) {
            List<org.jdom.Element> sampleEleList = samGroupEleList.get(i).getChildren("each_sample");
            List<Integer> sampleIndexList = new ArrayList<Integer>();
            for (org.jdom.Element eachSam : sampleEleList) {
                sampleIndexList.add(ctr);
                ctr++;
            }
            sampleListt.add(sampleIndexList);
        }
        while (f.hasNext()) {
            ProteinModel protein = (ProteinModel) f.next();
            int d = 0;
            int size1 = protein.getPeptideList().size();
            List classes = new ArrayList<>();
            List<List<Long>> intList = new ArrayList<>();
            List<List<Long>> totalList = new ArrayList<>();
            List<List<Long>> intIndividualPeptideList = new ArrayList<>();

            List<List<Long>> intListIIT = new ArrayList<>();
            List<List<Long>> intIndividualPeptideListIIT = new ArrayList<>();
             ctr = 0;

            List<Long> pepAvgIntList = null;
            List<Long> pepIntList = null;
            List<Long> pepAvgIntListIIT = null;
            List<Long> pepIntListIIT = null;
            double[] intensitySumArr = new double[samGroupEleList.size()];
            List<Double> ionInjectionTimeIntensitySumList =new ArrayList<>();
            long totalPeptideIntensity=0;

            for (ctr = 0; ctr < size1; ctr++) {

                int size2 = protein.getPeptideList().get(ctr).getPeptideList().size();
                if(ctr==0) ionInjectionTimeIntensitySumList.addAll(Collections.nCopies(size2,new Double(0)));
                //List<Double> intensity1 = new ArrayList<>();
                // List<Double> intensity2 = new ArrayList<>();
                List<Double> intensity = new ArrayList<>();
                List<Double> iitIntensity = new ArrayList<>();

                Iterator it = sampleListt.iterator();

                pepAvgIntList = new ArrayList<>();
                pepIntList = new ArrayList<>();

                pepAvgIntListIIT = new ArrayList<>();
                pepIntListIIT = new ArrayList<>();

                boolean hasReplicates=true;

                int sampleCount=0;

                while (it.hasNext()) {
                    List<Integer> sampleIndexList = (List<Integer>) it.next();
                    for (int u = 0; u < size2; u++) {
                        if (sampleIndexList.contains(u)) {
                            intensity.add(protein.getPeptideList().get(ctr).getPeptideList().get(u).getAverageIntensity());
                            iitIntensity.add(protein.getPeptideList().get(ctr).getPeptideList().get(u).getCorrIonInjectionIntensity());

                        }
                    }
                    double[] intensityy = new double[intensity.size()];
                    double[] iitIntensityy = new double[intensity.size()];

                    //Iterator<Double> iterator = intensity.iterator();
                    int t = 0;
                    long intSum = 0;
                    long intSumIIT=0;
                    //while (iterator.hasNext()) {
                    for(int i=0;i<intensity.size();i++) {

                        intensityy[t] = intensity.get(i);
                        //System.out.println("-------------------" + intensityy[t]);
                        intensitySumArr[sampleCount] += intensityy[t];

                        intSum += (long)intensityy[t];
                        pepIntList.add((long)intensityy[t]);


                        iitIntensityy[t] = iitIntensity.get(i);
                        intSumIIT += (long)iitIntensityy[t];
                        pepIntListIIT.add((long)iitIntensityy[t]);
                        t++;
                    }

                    //intensity sum     System.out.println("sum=============" + (long)intSum);

                    long averageIntensity = (long)intSum/t;
                    long averageIntensityIIT = (long)intSumIIT/t;

                    pepAvgIntList.add(averageIntensity);
                    pepAvgIntListIIT.add(averageIntensityIIT);
                    totalPeptideIntensity += intSum;
                    //System.out.println("---->" + intSum);
                    if(intensityy.length<=1)
                        hasReplicates=false;

                    classes.add(intensityy);
                    intensity.clear();
                    iitIntensity.clear();

                    sampleCount++;

                }

                //    System.out.println("===\t" + protein.getPeptideList().get(ctr).getPeptideList().get(0).getSequence() + "\t" + totalPeptideIntensity);

                intList.add(pepAvgIntList);
                intIndividualPeptideList.add(pepIntList);

                intListIIT.add(pepAvgIntListIIT);
                intIndividualPeptideListIIT.add(pepIntListIIT);
                for(int i =0; i<pepIntListIIT.size(); i++)
                {
                    ionInjectionTimeIntensitySumList.set(i,ionInjectionTimeIntensitySumList.get(i)+pepIntListIIT.get(i));
                }

                double pvalue = 1.0;
                if(hasReplicates)
                    pvalue = AnovaUtil.calculateAnovaPvalue(classes);

                pvalueList.add(pvalue);
            }

            //calculate most common pattern of changes for finding protein intensity
            int bestIndex=0;
            double currentSum=0;

            protein.setIntensitySumArr(intensitySumArr);

            /*
            for(long inten:intensitySumArr) {
                System.out.print(inten + "\t");
            }
            System.out.println("=============");
            */

            //if only two peptides, then calculate average
            if(intList.size()==2) {
                List<Long> averageList = new ArrayList<>();
                List<Long> averageIndPeptideList = new ArrayList<>();
                List<Long> averageIITList = new ArrayList<>();
                List<Long> averageIndPeptideIITList = new ArrayList<>();

                //   for(int i=0;i<2;i++) {

                for(int j=0;j<intList.get(0).size();j++) {
                    averageList.add(j, (intList.get(0).get(j) + intList.get(1).get(j))/2 );
                }

                for(int j=0;j<intIndividualPeptideList.get(0).size();j++) {
                    averageIndPeptideList.add(j, (intIndividualPeptideList.get(0).get(j) + intIndividualPeptideList.get(1).get(j))/2 );
                }

                for(int j=0;j<intListIIT.get(0).size();j++) {
                    averageIITList.add(j, (intListIIT.get(0).get(j) + intListIIT.get(1).get(j))/2 );
                }

                for(int j=0;j<intIndividualPeptideListIIT.get(0).size();j++) {
                    averageIndPeptideIITList.add(j, (intIndividualPeptideListIIT.get(0).get(j) + intIndividualPeptideListIIT.get(1).get(j))/2 );
                }
                // }

                //List<Long> bestList = intList.get(bestIndex);
                protein.setBestCorrelationIntensityList(averageList);
                /*protein.setBestCorrelationEachpeptideIntensityList(averageIndPeptideList);
                protein.setBestCorrelationIntensityListIIT(averageIITList);
                protein.setBestCorrelationEachpeptideIntensityListIIT(averageIndPeptideIITList);
                protein.setIonInjectionTimeIntensitySumList(ionInjectionTimeIntensitySumList);*/
                //avg_ion_injection_time_intensity calc
                //     List<Long> eachList = intList.get(0);


            } else {
                for (int i = 0; i < intList.size(); i++) {

                    List<Long> currentList = intList.get(i);
                    double corrSum = 0;
                    for (int j = 0; j < intList.size(); j++) {
                        if (i == j) continue;

                        List<Long> eachList = intList.get(j);

                        SimpleRegression regression = new SimpleRegression();
                        for (int k = 0; k < currentList.size(); k++) {
                            regression.addData(currentList.get(k), eachList.get(k));
                        }
                        double corr = regression.getR();

                        corrSum += corr;

                    }

                    if (corrSum > currentSum) {
                        bestIndex = i;
                        currentSum = corrSum;
                    }
                }

                List<Long> bestList = intList.get(bestIndex);
                protein.setBestCorrelationIntensityList(bestList);

                List<Long> bestIndividualPeptideList = intIndividualPeptideList.get(bestIndex);
                //protein.setBestCorrelationEachpeptideIntensityList(bestIndividualPeptideList);
                //each_exp_intensity_sum
                List<Long> bestListIIT = intListIIT.get(bestIndex);
                //protein.setBestCorrelationIntensityListIIT(bestListIIT);

                List<Long> bestIndividualPeptideListIIT = intIndividualPeptideListIIT.get(bestIndex);
                //protein.setBestCorrelationEachpeptideIntensityListIIT(bestIndividualPeptideListIIT);
                protein.setIonInjectionTimeIntensitySumList(ionInjectionTimeIntensitySumList);
                //avg_ion_injection_time_intensity calc

            }

            protein.setTotalPeptideIntensity( totalPeptideIntensity );
            HashMap<String, Integer> proteinLengthMap = getProteinLengthMap(fastaFile);

            int proteinLength = proteinLengthMap.get(protein.getRedundnatProteinList().get(0).getLocus());

            protein.setTotalIntensityDividedByProteinLengthLog10( Math.log10(totalPeptideIntensity/(double)proteinLength) );
            //totol_peptide_intensity_dividedby_length_log10
            protein.setProteinLength(proteinLength);

            intList.clear();
        }
    }

    public static HashMap<String, Integer> getProteinLengthMap(String fastaFile) throws Exception {

        String proteinLengthFile = fastaFile + ".length";
        if(!new File(proteinLengthFile).exists())
            FastaFileProteinLengthCalc.generateProteinLength(fastaFile);

        HashMap<String, Integer> proteinLengthMap = FastaReader.getProteinLengthMap(proteinLengthFile);
        return proteinLengthMap;
    }

}

