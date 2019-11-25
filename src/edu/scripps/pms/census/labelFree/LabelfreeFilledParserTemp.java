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
import edu.scripps.pms.census.util.io.FastaReader;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import rpark.statistics.AnovaUtil;

import java.io.*;
import java.math.BigDecimal;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Created by rpark on 4/15/16.
 * modified by rohan rampuria Titus
 */
public class LabelfreeFilledParserTemp {

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
  //  private List<Integer> intensitySumEachSampleList = new ArrayList<>();
    private List<Integer> peptideIntensitySumIndexList = new ArrayList<>();
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

    private int accessionIndex = -1;
    private int descriptionIndex  =-1;
    private int pepCountIndex  =-1;
    private int[] AvgNormIntensityIndex;
    private final String filledFile;

    private int[] sampleIndex = null;

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
        String confFile = args[0];
        String fname = args[1];
        String fastaFile = args[2];


        LabelfreeFilledParserTemp l = new LabelfreeFilledParserTemp(fname);
        List<ProteinModel> proteinList = l.readWholeFile(confFile);
        String baseFileName = fname.substring(0, fname.length() - 11);

        l.processProteinList(proteinList,confFile,fastaFile, baseFileName + "_IIT_intensity.txt");

        LabelfreeMissingPeptideBuilderSplit.generateLabelfreeOutputFile(proteinList, baseFileName + "_stat.txt",confFile);

    }
    public LabelfreeFilledParserTemp(String fileName) {
        this.filledFile = fileName;
        try {
            br = new BufferedReader(new FileReader(filledFile));
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }



    public List<ProteinModel> readWholeFile(String configFile) throws Exception {


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

      sampleIndex = new int[ctr];
      int indexTmp=0;
      for (int i = 0; i < count; i++) {
        List<org.jdom.Element> sampleEleList = samGroupEleList.get(i).getChildren("each_sample");
        for (org.jdom.Element eachSam : sampleEleList) {
          sampleIndex[indexTmp++] = i;
        }
      }



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
              //System.out.println(eachLine);
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
                        /*if (intensitySumEachSampleList.size() != 0) {
                            for (int index : intensitySumEachSampleList) {
                                if (words[index].equalsIgnoreCase("NA") || words[index].equalsIgnoreCase("NaN")) {
                                    proteinModel.addBestCorrelationEachpeptideIntensityList(-1);
                                } else {
                                    long l = new BigDecimal(words[index]).longValue();
                                    proteinModel.addBestCorrelationEachpeptideIntensityList(l);
                                }
                            }
                        } */
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

                              if(index>= words.length) continue;

                                if (words[index].equalsIgnoreCase("NA") || words[index].equalsIgnoreCase("NaN")) {
                                    proteinModel.addBestCorrelationIntensityListIIT(-1);
                                } else {
                                    long in = new BigDecimal(words[index]).longValue();
                                    proteinModel.addBestCorrelationIntensityListIIT(in );
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


                      if (peptideIntensitySumIndexList.size() != 0) {
                        for (int index : peptideIntensitySumIndexList) {
                          if (words[index].equalsIgnoreCase("NA") || words[index].equalsIgnoreCase("NaN")) {
                            proteinModel.addIntensityEachSampleSumList(0.0);

                          } else {
                            proteinModel.addIntensityEachSampleSumList(Double.parseDouble(words[index]));

                          }
                        }
                        //proteinModel.setIntensitySumArr(pepGroupIonInjectionArr);
                      }


                        proteinList.add(proteinModel);
                        isNewProtein = false;
                    }

                } else if(eachLine.startsWith("S\t")) {
                    isNewProtein = true;

                    chroPepList = parsePeptideLine(eachLine, samGroupEleList.size(), sampleIndex);
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
      protein.setProteinLine(currentLine);
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
           // else if (words[i].startsWith("INTENSITY_"))
           //     intensitySumEachSampleList.add(i);
            //else if (words[i].startsWith("PEPTIDE_INTENSITY_SUM_"))
            else if (words[i].startsWith("INTENSITY_"))
                peptideIntensitySumIndexList.add(i);
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

       // System.out.println("aaa");
    }


    public List<ChroPeptide> parsePeptideLine(String currentLine, int sampleGroupSize, int[] sampleGroupIndex) {
        List<ChroPeptide> cPeptideList = new ArrayList<>();
        String words[] = currentLine.split("\t");
        // int wordCounter = 0;

      double[] iitRatioArr  = new double[sampleGroupSize];
        if (intensityIndexList.size() > 0) {
          for(int i=0;i<intensityIndexList.size();i++) {
          //  System.out.println(intensityIndexList.get(i) );
            String iitIntensity = words[corrInjectionIntensityIndexList.get(i)];
            String intensity = words[intensityIndexList.get(i)];
           // System.out.println(intensity + "\t" +  iitIntensity ); //Double.parseDouble();

            if("NA".equals(iitIntensity) || "NA".equals(intensity)) continue;

            double iitCorrection = Double.parseDouble(intensity)/Double.parseDouble(iitIntensity);
            iitRatioArr[sampleGroupIndex[i]] = iitCorrection;

           // System.out.println(intensity + "\t" +  iitIntensity + "\t" + iitRatioArr[sampleGroupIndex[i]] + "\t" + i); //Double.parseDouble();
          }

          /*
          for(int i=0;i<intensityIndexList.size();i++) {
            String iitIntensity = words[corrInjectionIntensityIndexList.get(i)];
            String intensity = words[intensityIndexList.get(i)];

            double correction = iitRatioArr[sampleGroupIndex[i]];
            System.out.println("corr11" + correction);
            if(correction==0 || !"NA".equals(iitIntensity) ||"NA".equals(intensity)) continue;
            System.out.println("corr22" + correction);



          } */

        }

        for (int i = 0; i < totalExperiments; i++) {

            ChroPeptide currentPeptide = new ChroPeptide();
            if (sequenceIndexList.size() > 0) {

                currentPeptide.setSequence(words[sequenceIndexList.get(i)]);

            }
            if (intensityIndexList.size() > 0) {
                currentPeptide.setAverageIntensity(Double.parseDouble(words[intensityIndexList.get(i)]));
            }

            if (corrInjectionIntensityIndexList.size() > 0) {
              String corrIITIntensity = words[corrInjectionIntensityIndexList.get(i)];
              //System.out.println("-----------\t" + corrIITIntensity);
              if (!"NA".equals(corrIITIntensity)) {
                currentPeptide.setCorrIonInjectionIntensity(Double.parseDouble(corrIITIntensity));
              } else {
                double correction = iitRatioArr[sampleGroupIndex[i]];
                double intensity = Double.parseDouble(words[intensityIndexList.get(i)]);
                double correctedIntensity=0;

                if(correction>0)
                  correctedIntensity= intensity/correction;

                currentPeptide.setCorrIonInjectionIntensity(correctedIntensity);

              }
            }
            if (fileNameIndexList.size() > 0) {
                if (words[fileNameIndexList.get(i)].equalsIgnoreCase("NA")) {
                    cPeptideList.add(currentPeptide);
                    continue;
                }
                else{
                    currentPeptide.setFileName(words[fileNameIndexList.get(i)]);
                }
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
	    }else if (words[i].equalsIgnoreCase("CORRIONINJECTION_INTENISTY") || words[i].equalsIgnoreCase("CORRIONINJECTION_INTENSITY")) {
                corrInjectionIntensityIndexList.add(i);
            }

        }
    }

    public void processProteinList(List<ProteinModel> proteinList, String configFile, String fastaFile, String iitProteinFile) throws Exception
    {


      BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(iitProteinFile));
      PrintStream p = new PrintStream(out);



      HashMap<String, Integer> proteinLengthMap = getProteinLengthMap(fastaFile);


        Iterator f = proteinList.iterator();
        List<Double> pvalueList = new ArrayList<>();
        //List<Double> pvalueListIIT = new ArrayList<>();

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
            //List classesIIT = new ArrayList<>();
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


          StringBuffer sb = new StringBuffer();
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
              boolean hasReplicatesIIT=true;

                int sampleCount=0;

              boolean zeroValueIIT=false;
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

                        if(iitIntensityy[t]==0) zeroValueIIT=true;

                        pepIntListIIT.add((long)iitIntensityy[t]);
                        t++;
                    }

                    /*
                    if(!zeroValueIIT) {
                      for(int i=0;i<intensity.size();i++)
                        protein.addEachIITPeptideToList(sampleCount, iitIntensity.get(i));
                    } */

                    //intensity sum     System.out.println("sum=============" + (long)intSum);

                    long averageIntensity = (long)intSum/t;
                    long averageIntensityIIT = (long)intSumIIT/t;

                    pepAvgIntList.add(averageIntensity);
                    pepAvgIntListIIT.add(averageIntensityIIT);
                    totalPeptideIntensity += intSum;
                    //System.out.println("---->" + intSum);
                    if(intensityy.length<=1)
                        hasReplicates=false;

             //       if(iitIntensityy.length<=1)
               //       hasReplicatesIIT=false;

                    classes.add(intensityy);
                   // classesIIT.add(iitIntensityy);

                    intensity.clear();
                    iitIntensity.clear();

                    sampleCount++;

                }

                //    System.out.println("===\t" + protein.getPeptideList().get(ctr).getPeptideList().get(0).getSequence() + "\t" + totalPeptideIntensity);

                intList.add(pepAvgIntList);
                intIndividualPeptideList.add(pepIntList);

                intListIIT.add(pepAvgIntListIIT);
                intIndividualPeptideListIIT.add(pepIntListIIT);


                if(!zeroValueIIT) {
                  sb.append("S\t");

                  for (int i = 0; i < pepIntListIIT.size(); i++) {
                    sb.append(pepIntListIIT.get(i)).append("\t");

                    ionInjectionTimeIntensitySumList.set(i, ionInjectionTimeIntensitySumList.get(i) + pepIntListIIT.get(i));

                  }

                  int tcount=0;
                  for(int i=0;i<sampleListt.size();i++) {
                    List l = (List)sampleListt.get(i);

                    for(int j=0;j<l.size();j++) {

                      protein.addEachIITPeptideToList(i, pepIntListIIT.get(tcount++));

                    }
                  }

                  sb.append("\n");
                }


              List<Double> ionInjectionTimeIntensityAvgGroupList = new ArrayList<>();
              for(Object obj:samGroupEleList) {
                ionInjectionTimeIntensityAvgGroupList.add(0.0);
              }

              for(int i=0;i<ionInjectionTimeIntensitySumList.size();i++) {
                double tmpInt = ionInjectionTimeIntensityAvgGroupList.get(sampleIndex[i]);
                tmpInt += ionInjectionTimeIntensitySumList.get(i);

                ionInjectionTimeIntensityAvgGroupList.set(sampleIndex[i], tmpInt);
              }

              protein.setIonInjectionTimeIntensityAvgGroupList(ionInjectionTimeIntensityAvgGroupList);



              double pvalue = 1.0;
              //  double pvalueIIT=1.0;
                if(hasReplicates)
                    pvalue = AnovaUtil.calculateAnovaPvalue(classes);

               //  if(hasReplicatesIIT)
               //   pvalueIIT = AnovaUtil.calculateAnovaPvalue(classesIIT);

                pvalueList.add(pvalue);
               // pvalueListIIT.add(pvalueIIT);
            }

            //calculate most common pattern of changes for finding protein intensity
            int bestIndex=0;
            double currentSum=0;

            protein.setIntensitySumArr(intensitySumArr);
          protein.setIonInjectionTimeIntensitySumList(ionInjectionTimeIntensitySumList);


          if(sb.length()>0) {

            List<ChroProtein> plist = protein.getRedundnatProteinList();
            for(ChroProtein cp:plist) {

              p.print("P\t" + cp.getLocus());
              p.print("\t");
            }
            for(int i=0;i<ionInjectionTimeIntensitySumList.size();i++) {
              p.print(ionInjectionTimeIntensitySumList.get(i));
              p.print("\t");
            }

            for(Double avgInt:protein.getIonInjectionTimeIntensityAvgGroupList()) {
              p.print(avgInt);
              p.print("\t");
            }

            p.print("\n");

            p.print(sb.toString());

          }


//              int count = samGroupEleList.size();

          //int[] sampleIndex = new int[ctr];
          List<Double> arrList1 = new ArrayList<>();
          List<Double> arrList2 = new ArrayList<>();
          for (int i = 0; i < ionInjectionTimeIntensitySumList.size(); i++) {
            if(sampleIndex[i] == 0) arrList1.add(ionInjectionTimeIntensitySumList.get(i));
            else if(sampleIndex[i] == 1) arrList2.add(ionInjectionTimeIntensitySumList.get(i));
          }


          double[] arr1 = new double[arrList1.size()];
          double[] arr2 = new double[arrList2.size()];
          for(int i=0;i<arrList1.size();i++)
            arr1[i] = arrList1.get(i);
          for(int i=0;i<arrList2.size();i++)
            arr2[i] = arrList2.get(i);

          List iitClasses = new ArrayList();
          iitClasses.add(0,arr1);
          iitClasses.add(1,arr2);

          double iitPvalue = AnovaUtil.calculateAnovaPvalue(iitClasses);


          protein.setIitPvalue(iitPvalue);

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
                protein.setBestCorrelationEachpeptideIntensityList(averageIndPeptideList);
                protein.setBestCorrelationIntensityListIIT(averageIITList);
                protein.setBestCorrelationEachpeptideIntensityListIIT(averageIndPeptideIITList);
                protein.setIonInjectionTimeIntensitySumList(ionInjectionTimeIntensitySumList);



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
                protein.setBestCorrelationEachpeptideIntensityList(bestIndividualPeptideList);
                //each_exp_intensity_sum
                List<Long> bestListIIT = intListIIT.get(bestIndex);
                protein.setBestCorrelationIntensityListIIT(bestListIIT);

                List<Long> bestIndividualPeptideListIIT = intIndividualPeptideListIIT.get(bestIndex);
                protein.setBestCorrelationEachpeptideIntensityListIIT(bestIndividualPeptideListIIT);
                //protein.setIonInjectionTimeIntensitySumList(ionInjectionTimeIntensitySumList);
                //avg_ion_injection_time_intensity calc

            }

            protein.setTotalPeptideIntensity( totalPeptideIntensity );

          //  int proteinLength = proteinLengthMap.get(protein.getRedundnatProteinList().get(0).getLocus());



          Integer lengthObj = proteinLengthMap.get(protein.getRedundnatProteinList().get(0).getLocus());
          int proteinLength = 1;
          if(null != lengthObj)
            proteinLength = lengthObj.intValue();



            protein.setTotalIntensityDividedByProteinLengthLog10( Math.log10(totalPeptideIntensity/(double)proteinLength) );
            //totol_peptide_intensity_dividedby_length_log10
            protein.setProteinLength(proteinLength);

            intList.clear();
        }

        p.close();
        out.close();

    }

    public static HashMap<String, Integer> getProteinLengthMap(String fastaFile) throws Exception {

        String proteinLengthFile = fastaFile + ".length";
        if(!new File(proteinLengthFile).exists())
            FastaFileProteinLengthCalc.generateProteinLength(fastaFile);

        HashMap<String, Integer> proteinLengthMap = FastaReader.getProteinLengthMap(proteinLengthFile);
        return proteinLengthMap;
    }

}

