/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.scripps.pms.census.labelFree;

import static edu.scripps.pms.census.labelFree.LabelFreeParser.main;
import static edu.scripps.pms.census.labelFree.LabelFreeParser.printUsage;
import static rpark.statistics.GaussianFitting.getGaussianCurveFitRange;

import java.io.*;
import java.util.*;

import edu.scripps.pms.census.CensusConstants;
import edu.scripps.pms.census.ChroGenerator;
import edu.scripps.pms.census.ElementComposition;
import edu.scripps.pms.census.conf.Configuration;
import edu.scripps.pms.census.hash.IndexedFile;
import edu.scripps.pms.census.io.ChroXmlReader;
import edu.scripps.pms.census.io.FastaFileProteinLengthCalc;
import edu.scripps.pms.census.io.IsotopeReader;
import edu.scripps.pms.census.labelFree.json.ChroJSONCreator;
import edu.scripps.pms.census.labelFree.json.LabelFreeJSONPeptide;
import edu.scripps.pms.census.labelFree.json.LabelFreeJSONProtein;
import edu.scripps.pms.census.labelFree.model.LabelfreePeptide;
import edu.scripps.pms.census.labelFree.util.LabelfreeChroUtil;
import edu.scripps.pms.census.model.ChroPeptide;
import edu.scripps.pms.census.model.ChroProtein;
import edu.scripps.pms.census.model.SampleModel;
import edu.scripps.pms.census.util.CalcUtilGeneric;
import edu.scripps.pms.census.util.IsotopeDist;
import edu.scripps.pms.census.util.io.FastaReader;
import edu.scripps.pms.census.util.io.FileUtil;
import edu.scripps.pms.util.sqlite.spectra.SpectraDB;
import gnu.trove.TDoubleArrayList;
import gnu.trove.TDoubleIntHashMap;

import gnu.trove.TObjectIntHashMap;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.inference.TestUtils;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.paukov.combinatorics.Factory;
import org.paukov.combinatorics.Generator;
import org.paukov.combinatorics.ICombinatoricsVector;
import rpark.statistics.*;
import rpark.statistics.model.GaussianPeakModel;
import scripts.MSSplitFolderCreation;
import edu.scripps.pms.census.hash.IndexUtil;

/**
 *
 * @author rpark
 * @author rohan
 * @author titus
 */
public class LabelfreeMissingPeptideBuilderSplit {

    public static int SMOOTH_WINDOW_SIZE = 5;
    public static int SMOOTH_WINDOW_SIZE_7 = 5;

    public static void main(String[] args) throws Exception {

        if (args.length < 3) {
            printUsage();
            return;
        }

        System.out.println("running label free analysis..");

        /***********************************************
         * to re-run without building missling proteins, but just read text output file, refer to
         *
         *  edu.scripps.pms.census.labelFree.LabelfreeFilledParser.main();
         *
         ***********************************************/

        String configFile = args[0];
        String tmpFile = args[1]; //http://192.168.2.9/ip2/viewLabelfree.html?pid=127&projectName=JonB_Cox
        String jsonFile = args[2];


        //   String configFile = "/data/2/rpark/ip2_data/rpark/fusion_labelfree/labelfree_quant/labelfree_10985/temp/census_config_labelfree_10985.xml";
          //       String tmpFile = "/data/2/rpark/ip2_data/rpark/fusion_labelfree/labelfree_quant/labelfree_10985/temp/census_labelfree_out_10985.txttmp";
            //  String jsonFile = "/data/2/rpark/ip2_data/rpark/fusion_labelfree/labelfree_quant/labelfree_10985/temp/census_labelfree_out_10985.json";



        //file /data/2/rpark/ip2_data/rpark/fusion_labelfree/labelfree_quant/labelfree_10985/temp/census_labelfree_out_10985.txttmp
        //String path="/data/2/rpark/ip2_data/rpark/fusion_labelfree/labelfree_quant/labelfree_10985/temp/";
//           String configFile = "/data/2/rpark/ip2_data/rpark/fusion_labelfree/labelfree_quant/labelfree_10985/temp/census_config_labelfree_10985.xml";
  //         String tmpFile = "/data/2/rpark/ip2_data/rpark/fusion_labelfree/labelfree_quant/labelfree_10985/temp/census_labelfree_out_10985.txttmp";
    //      String jsonFile = "/data/2/rpark/ip2_data/rpark/fusion_labelfree/labelfree_quant/labelfree_10985/temp/census_labelfree_out_10985.json";

        // System.out.println(configFile);
        //  System.out.println(tmpFile);
        //  System.out.println(jsonFile);
        // LabelfreeMissingPeptideBuilderSplit.runLabelfree(args[0], args[1], args[2]);
        Configuration conf = Configuration.getInstance();
      //  conf.setLabelfreeCheckChargeState(true);

        if (!conf.isReadConfigFile()) {
            conf.setLabelfree(true);
            conf.readXMLParam(configFile);
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

        LabelfreeMissingPeptideBuilderSplit.runLabelfree(configFile, tmpFile, jsonFile, sample1IndexList, sample2IndexList);
      String mainFileName = tmpFile.substring(0, tmpFile.length() - 7);
      String filledFile = mainFileName + "_filled.txt";
      //fill IIT intensities for unidentified peptides




        LabelfreeFilledParserTemp l = new LabelfreeFilledParserTemp(filledFile);
        List<ProteinModel> proteinList = l.readWholeFile(configFile);

      String baseFileName = filledFile.substring(0, filledFile.length() - 11);

      String fastaFile = getFastaFile(conf);
      l.processProteinList(proteinList,configFile,fastaFile, baseFileName + "_IIT_intensity.txt");


        LabelfreeMissingPeptideBuilderSplit.generateLabelfreeOutputFile(proteinList, mainFileName + "_stat.txt", configFile);

        /*
         LabelFreeParser labelFreeParser = new LabelFreeParser();

         if(!labelFreeParser.checkInputFileExist(tmpFile,configFile)){
         System.out.println("Input File Not Found ");
         return;
         }

         labelFreeParser.setTxtTmpFile(tmpFile);
         */
    }


    public static String getFastaFile(Configuration conf) throws Exception {

      List<SampleModel> sampleList = conf.getSampleList();
      String path = sampleList.get(0).getPathList().get(0);
      BufferedReader br = new BufferedReader(new FileReader(path + File.separator + "sequest.params"));
      String eachLine;
      while( null != (eachLine = br.readLine()) ) {
        if(eachLine.startsWith("database_name")) break;
      }

      String[] arr = eachLine.split("=");
      String fastaFile = arr[1].trim();

      br.close();

      return fastaFile;
    }


    public static HashMap<String, Integer> getProteinLengthMap(Configuration conf) throws Exception {

        List<SampleModel> sampleList = conf.getSampleList();
        String path = sampleList.get(0).getPathList().get(0);
        BufferedReader br = new BufferedReader(new FileReader(path + File.separator + "sequest.params"));
        String eachLine;
        while( null != (eachLine = br.readLine()) ) {
            if(eachLine.startsWith("database_name")) break;
        }

    //    System.out.println("===" + path);
        String[] arr = eachLine.split("=");
        String fastaFile = arr[1].trim();

        String proteinLengthFile = fastaFile + ".length";
        if(!new File(proteinLengthFile).exists())
            FastaFileProteinLengthCalc.generateProteinLength(fastaFile);

        HashMap<String, Integer> proteinLengthMap = FastaReader.getProteinLengthMap(proteinLengthFile);

        br.close();
        return proteinLengthMap;
    }

    public static TObjectIntHashMap getProteinLengthTMap(Configuration conf) throws Exception {

        List<SampleModel> sampleList = conf.getSampleList();
        String path = sampleList.get(0).getPathList().get(0);
        BufferedReader br = new BufferedReader(new FileReader(path + File.separator + "sequest.params"));
        String eachLine;
        while( null != (eachLine = br.readLine()) ) {
            if(eachLine.startsWith("database_name")) break;
        }

        //    System.out.println("===" + path);
        String[] arr = eachLine.split("=");
        String fastaFile = arr[1].trim();

        String proteinLengthFile = fastaFile + ".length";
        if(!new File(proteinLengthFile).exists())
            FastaFileProteinLengthCalc.generateProteinLength(fastaFile);

        TObjectIntHashMap proteinLengthMap = FastaReader.getProteinLengthTMap(proteinLengthFile);

        br.close();
        return proteinLengthMap;
    }


    public static void redoNormalization(List<ProteinModel> proteinList) {
      //List<Double> medianList = new ArrayList<>();
      int length = proteinList.get(0).getIntensityEachSampleSumList().size();
      double [] medianArray = new double[length];
      //double[] intensityArray = new double[proteinList.size()];
        double [] totalIntensitySumArr  = new double[length];

//fix norm intensity
      for (int i = 0; i < length; i++)
      {
          TDoubleArrayList tIntensityList = new TDoubleArrayList();

          //List<Double> intensityList = new ArrayList<>();
        for(int j=0; j<proteinList.size(); j++)
          {
              if(proteinList.get(j).getIntensityEachSampleSumList().get(i)>0){
                  tIntensityList.add(proteinList.get(j).getIntensityEachSampleSumList().get(i));
                  totalIntensitySumArr[i]+=proteinList.get(j).getIntensityEachSampleSumList().get(i);
              }
          }
          double [] intensityArray = tIntensityList.toNativeArray();
        medianArray[i]=CommonStat.getMedianValue(intensityArray);

  //      System.out.println("====\t" + medianArray[i]);

     //   for(double d:intensityArray)
     //     System.out.println(d);

      }
      double average = CommonStat.getMeanValue(totalIntensitySumArr);

    //  System.out.println("avg====\t" + average);
      for (int i = 0; i < length; i++)
      {
        for(int j=0; j<proteinList.size(); j++)
        {
            ProteinModel model = proteinList.get(j);

            double intensity= model.getIntensityEachSampleSumList().get(i);
            //double normIntensity= intensity/medianArray[i]*average;
            double normIntensity= intensity/totalIntensitySumArr[i]*average;

            model.getNormIntensityList().set(i,normIntensity);
        }
      }


    }


    public static void generateLabelfreeOutputFile(List<ProteinModel> proteinList, String filename, String confFile) throws Exception {


        Configuration conf = Configuration.getInstance();
        if (!conf.isReadConfigFile()) {
            conf.setLabelfree(true);
            conf.readXMLParam(confFile);
        }

        redoNormalization(proteinList);

        //HashMap<String,List<Double>> intensityMap = new HashMap<>();
        int sCountSize =0;
        HashMap<String,List<Double>> intensityMap = new HashMap<>();
        HashMap<String,List<Long>> intensityCorrMap = new HashMap<>();
        HashMap<String,List<Double>> iitIntensityMap = new HashMap<>();
        HashMap<String,List<List<Double>>> iitAllPepIntensityMap = new HashMap<>();

        HashMap<String,List<Double>> normIntensityMap = new HashMap<>();
        for (int i = 0; i < proteinList.size(); i++) {
            ProteinModel p = proteinList.get(i);
            List<List<Double>> avgIntensity = new ArrayList<>();
            List<List<Double>> avgCorrIonIntensity = new ArrayList<>();
            List<LabelfreePeptide> pepList = p.getPeptideList();

            for (int j = 0; j < pepList.size(); j++) {
                LabelfreePeptide pep = pepList.get(j);
                List<ChroPeptide> expPep = pep.getPeptideList();

                /*
                boolean validPeptide =true;
                for(ChroPeptide cPep:expPep) {

                    cPep.getPeakArea();

                    String fileName = cPep.getFileName();
                    if(fileName == null) {
                        validPeptide = false;
                        break;
                    }

                }

                if(!validPeptide) continue;
                */


                List<Double> l = new ArrayList<>();
                List<Double> l2 = new ArrayList<>();
                List<Double> norms = new ArrayList<>();
                for (int k = 0; k < expPep.size(); k++) {

                    l.add(expPep.get(k).getAverageIntensity());
                    l2.add(expPep.get(k).getCorrIonInjectionIntensity());
                }
                avgIntensity.add(l);
                avgCorrIonIntensity.add(l2);
            }
            List<Double> avgListIntensity = new ArrayList<>();
            List<Double> avgListCorrIntensity = new ArrayList<>();
            if(avgIntensity.size() == 0 || avgCorrIonIntensity.size() == 0){
                continue;
            }

            for (int h=0;h<avgIntensity.get(0).size();h++) {
                double sum2 =0;
                double sum1 =0;
                double avg =0;
                double avg1 =0;
                int count =0;

                for(int g=0;g<avgIntensity.size();g++) {
                    sum2 = sum2+avgIntensity.get(g).get(h);
                    sum1 = sum1+avgCorrIonIntensity.get(g).get(h);
                    count++;
                }
                avg= sum2/count;
                avg1 = sum1/count;
                avgListIntensity.add(avg);
                avgListCorrIntensity.add(avg1);
            }

            for(int a=0;a<p.getRedundnatProteinList().size();a++){
                //    intensityMap.put(p.getRedundnatProteinList().get(a).getLocus(), avgListIntensity);
                // intensityMap.put(p.getRedundnatProteinList().get(a).getLocus(), p.getBestCorrelationIntensityList());
                // intensityCorrMap.put(p.getRedundnatProteinList().get(a).getLocus(), p.getBestCorrelationIntensityListIIT());
                intensityMap.put(p.getRedundnatProteinList().get(a).getLocus(), p.getIntensityEachSampleSumList());
                intensityCorrMap.put(p.getRedundnatProteinList().get(a).getLocus(), p.getBestCorrelationEachpeptideIntensityListIIT());
                iitIntensityMap.put(p.getRedundnatProteinList().get(a).getLocus(), p.getIonInjectionTimeIntensitySumList());
                iitAllPepIntensityMap.put(p.getRedundnatProteinList().get(a).getLocus(), p.getAllIITPeptideList());
                normIntensityMap.put(p.getRedundnatProteinList().get(a).getLocus(),p.getNormIntensityList());
                sCountSize = p.getRedundnatProteinList().get(0).getSpecCountList().size();
            }
        }

        BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(filename));
        PrintStream p = new PrintStream(out);


        List<org.jdom.Element> samGroupEleList = conf.getRootConfEle().getChildren("sample");
        int expSize = 0;
        p.print("H\tCREATED_DATE\t");
        p.println(new Date());


        List<Integer> indexList = new ArrayList<>();

        for (Iterator<org.jdom.Element> samgItr = samGroupEleList.iterator(); samgItr.hasNext();) {
            int count =0;
            org.jdom.Element groupEle = samgItr.next();
            p.print("H\tGROUP_SAMPLE\t");
            p.print(groupEle.getAttributeValue("group"));
            p.print("\t");

            List<org.jdom.Element> sampleEleList = groupEle.getChildren("each_sample");

            for (Iterator<org.jdom.Element> samItr = sampleEleList.iterator(); samItr.hasNext();) {
                org.jdom.Element eachSample = samItr.next();
                p.print(eachSample.getAttributeValue("name"));
                p.print("\t");
                count++;
                expSize++;
            }

            p.println();
            indexList.add(count);
        }
        int smallestIndex = Integer.MAX_VALUE;
        for(int i: indexList)
        {
            smallestIndex = smallestIndex>i?i:smallestIndex;
        }
        /*
        calculating anova p-value
        */
        List<Double> pvalues = new ArrayList<>();
        List<Double> pvalues2 = new ArrayList<>();
        List<Double> pvaluesIIT = new ArrayList<>();
        List<Double> pvaluesIITAllPep = new ArrayList<>();

        SortedMap<String,Double> pNameMap = new TreeMap<>();
        Map<String,String> key2NameMap = new HashMap<>();
        HashMap<String,Double> pvalueMap = new HashMap<>();
        HashMap<String,Double> pvalueMap2 = new HashMap<>();
        HashMap<String,Double> pvalueMapIIT = new HashMap<>();
      HashMap<String,Double> pvalueMapIITAllPep = new HashMap<>();

        List<Double> normPvalues = new ArrayList<>();
        SortedMap<String,Double> normPvalueMap = new TreeMap<>();
        SortedMap<String,Double> normQvalueMap = new TreeMap<>();

        //SortedMap<String,List<Double>> nameRatioMap = new TreeMap<>();
        SortedMap<String,List<Double>> nameIntensityRatioPMap = new TreeMap<>();
        SortedMap<String,List<Double>>nameNormRatioPMap = new TreeMap<>();
        Map<String,List<Double>> nameIntensityRatioQMap = new HashMap<>();
        Map<String,List<Double>> nameNormRatioQMap = new HashMap<>();

        SortedMap<String,Double> logIntensityPMap = new TreeMap<>();
        SortedMap<String,Double>logNormPMap = new TreeMap<>();
        Map<String,Double> logIntensityQMap = new HashMap<>();
        Map<String,Double> logNormQMap = new HashMap<>();


        //avg_norm_intensity


        SortedMap<String,List<Double>> comparedLogIntensityPMap = new TreeMap<>();
        SortedMap<String, List<Double>> comparedLogNormPMap = new TreeMap<>();
        Map<String,List<Double>>  comparedLogIntensityQMap = new HashMap<>();
        Map<String,List<Double>> comparedLogNormQMap = new HashMap<>();

        Map<String,List<Double>> avgNormIntensityMap = new HashMap<>();
        List<String> nameList = new ArrayList<>();
        List<List<Double>> intensityRatioPValueMatrix = new ArrayList<>();
        List<List<Double>> normIntensityRatioPValueMatrix = new ArrayList<>();

        List<List<Double>> logIntensityRatioPValueMatrix = new ArrayList<>();
        List<List<Double>> logNormIntensityRatioPValueMatrix = new ArrayList<>();


        for(int i=0; i<indexList.size()-1;i++)
        {
            intensityRatioPValueMatrix.add(new ArrayList<Double>());
            normIntensityRatioPValueMatrix.add(new ArrayList<Double>());
            logIntensityRatioPValueMatrix.add(new ArrayList<Double>());
            logNormIntensityRatioPValueMatrix.add(new ArrayList<Double>());
        }


        for(String key:intensityMap.keySet()){
            List classes = new ArrayList();
            List classes2 = new ArrayList();
            List iitClasses = new ArrayList();
            List iitClassesAllPep = new ArrayList();

            List normClasses =new ArrayList();

            List logClasses =new ArrayList();
            List logNormClasses =new ArrayList();

            List<Double> avgNormIntensityList = new ArrayList<>();

            List<Double> dlist = new ArrayList<>();
            int counter=0;
            boolean hasReplicates=true;
            for(int i=0;i<indexList.size();i++){

                double [] intensity =new double[indexList.get(i)];
                double [] normintensity = new double[indexList.get(i)];
                double [] corrIntensity =new double[indexList.get(i)];
                double [] iitIntensity =new double[indexList.get(i)];


                double [] logIntensity = new double[indexList.get(i)];
                double [] logNormIntensity = new double[indexList.get(i)];

                double normSum =0;
                for(int j=0;j<indexList.get(i);j++){

                 // System.out.println("-------------\t" + intensity.length + "\t" + intensityMap.size() + "\t" + key  + "\t" + counter);
                 // System.out.println("-------------\t" + intensityMap);
                    intensity[j] = intensityMap.get(key).get(counter);
                    corrIntensity[j] = intensityCorrMap.get(key).get(counter);
                    iitIntensity[j] = iitIntensityMap.get(key).get(counter);

                   //System.out.println("-------------\t" + iitIntensity[i]);
                    normintensity[j] = normIntensityMap.get(key).get(counter);
                    if(normintensity[j]<=0) normintensity[j] = 1000;
                    normSum+=normintensity[j];
                    logIntensity[j] = Math.log(intensity[j])/Math.log(2);
                    logNormIntensity[j] = Math.log(normintensity[j])/Math.log(2);



                    counter++;
                    dlist.add(intensity[j]);
                }
                avgNormIntensityList.add(normSum/indexList.get(i));

                classes.add(intensity);
                classes2.add(corrIntensity);
                iitClasses.add(iitIntensity);
                normClasses.add(normintensity);

              List allpepList = iitAllPepIntensityMap.get(key);

              if(allpepList.size()<=0) {
                //if no values, then assign zeros from iitIntensity
                iitClassesAllPep.add(iitIntensity);
              } else {

                List<Double> iitAllPepList = iitAllPepIntensityMap.get(key).get(i);
                double[] iitArrAllPep = new double[iitAllPepList.size()];
                //double[] iitArrAllPep = new double[iitPepList.size()];
                for(int iit=0;iit<iitAllPepList.size();iit++)
                  iitArrAllPep[iit] = iitAllPepList.get(iit);

                iitClassesAllPep.add(iitArrAllPep);
              }



                logClasses.add(logIntensity);
                logNormClasses.add(logNormIntensity);

                if(intensity.length<=1 || corrIntensity.length<=1 || normintensity.length<=1)
                    hasReplicates=false;
            }
            avgNormIntensityMap.put(key,avgNormIntensityList);

            List<double[]> intensityRatioList = createLogRatios(classes);
            List<Double> intensityRatioPValueList = createPValuesOneSampleTTest(intensityRatioList);

            List<double[]> normRatioList = createLogRatios(normClasses);
            List<Double> normRatioPValueList = createPValuesOneSampleTTest(normRatioList);

            double pvalue = 1.0;
            double ratioPValue = 1.0;
            double pvalueIITAllPep= 1.0;
            if(hasReplicates) {
                pvalue = AnovaUtil.calculateAnovaPvalue(classes);
              pvalueIITAllPep = AnovaUtil.calculateAnovaPvalue(iitClassesAllPep);
                //ratioPValue = AnovaUtil.calculateAnovaPvalue(intensityRatioList);


            }
            if(Double.isInfinite(pvalue) || Double.isNaN(pvalue)) {
                pvalue = 1;
                ratioPValue = 1;
            }

          if(Double.isInfinite(pvalueIITAllPep) || Double.isNaN(pvalueIITAllPep)) {
            pvalueIITAllPep = 1;
          }


            double pvalue2 = 1.0;
            double pvalueIIT = 1.0;
            if(hasReplicates) {
              pvalue2 = AnovaUtil.calculateAnovaPvalue(classes2);
              pvalueIIT = AnovaUtil.calculateAnovaPvalue(iitClasses);
            }

            if(Double.isInfinite(pvalue2) || Double.isNaN(pvalue2))
                pvalue2 = 1;

          if(Double.isInfinite(pvalueIIT) || Double.isNaN(pvalueIIT))
            pvalueIIT = 1;

            double normPvalue = 1.0;
            double logNormPvalue = 1.0;
            double logPvalue = 1.0;
            if(hasReplicates) {
                normPvalue = AnovaUtil.calculateAnovaPvalue(normClasses);
                //ratioNormPValue = AnovaUtil.calculateAnovaPvalue(normRatioList);
                logNormPvalue = AnovaUtil.calculateAnovaPvalue(logNormClasses);
                logPvalue = AnovaUtil.calculateAnovaPvalue(logClasses);

            }

            if(Double.isInfinite(normPvalue) || Double.isNaN(normPvalue)) {
                normPvalue = 1;
            }
            if(Double.isInfinite(logNormPvalue)|| Double.isNaN(logNormPvalue))
            {
                logNormPvalue = 1;
            }
            if(Double.isInfinite(logPvalue)|| Double.isNaN(logPvalue))
            {
                logPvalue = 1;
            }

            pvalues.add(pvalue);
            pvalues2.add(pvalue2);
            pvaluesIIT.add(pvalueIIT);
            pvaluesIITAllPep.add(pvalueIITAllPep);


            pvalueMap.put(key, pvalue);
            pvalueMap2.put(key, pvalue2);
            pvalueMapIIT.put(key, pvalueIIT);
            pvalueMapIITAllPep.put(key, pvalueIITAllPep);

            String proteinName = createProteinString(dlist);
            List<Double> comparedLogPvalueList = createPValuesCompareWithFirst(logClasses);
            List<Double> comparedLogNormPvalueList = createPValuesCompareWithFirst(logNormClasses);

            if(!pNameMap.containsKey(proteinName))
            {
                nameList.add(proteinName);
                for(int i=0; i<intensityRatioPValueList.size();i++)
                {
                    intensityRatioPValueMatrix.get(i).add(intensityRatioPValueList.get(i));
                    normIntensityRatioPValueMatrix.get(i).add(normRatioPValueList.get(i));
                    logIntensityRatioPValueMatrix.get(i).add(comparedLogPvalueList.get(i));
                    logNormIntensityRatioPValueMatrix.get(i).add(comparedLogNormPvalueList.get(i));
                }
            }


            pNameMap.put(proteinName,pvalue);
            normPvalueMap.put(proteinName,normPvalue);

            nameIntensityRatioPMap.put(proteinName,intensityRatioPValueList);
            nameNormRatioPMap.put(proteinName,normRatioPValueList);
            key2NameMap.put(key,proteinName);

            logIntensityPMap.put(proteinName,logPvalue);
            logNormPMap.put(proteinName,logNormPvalue);





            comparedLogIntensityPMap.put(proteinName,comparedLogPvalueList);
            comparedLogNormPMap.put(proteinName, comparedLogNormPvalueList);

        }

        /*
        BH correction q-value
        */
        List<Double> bhcorr = BHCorrection.runBHCorrectionRemoveRedundant(pvalues);
        List<Double> bhcorr2 = BHCorrection.runBHCorrectionRemoveRedundant(pvalues2);
        List<Double> bhcorrIIT = BHCorrection.runBHCorrectionRemoveRedundant(pvaluesIIT);
        List<Double> bhcorrIITAllPep = BHCorrection.runBHCorrectionRemoveRedundant(pvaluesIITAllPep);
        // List<Double> bhcorr3 = BHCorrection.runBhCorrection(ratiobasedpval);
        HashMap<String,Double> bhcorrMap = new HashMap<>();
        HashMap<String,Double> bhcorrMap2 = new HashMap<>();
        HashMap<String,Double> bhcorrMapIIT = new HashMap<>();
        HashMap<String,Double> bhcorrMapIITAllPep = new HashMap<>();

        int l=0;
        for(String key : intensityMap.keySet()){
            bhcorrMap.put(key,bhcorr.get(l));
            bhcorrMap2.put(key, bhcorr2.get(l));
            bhcorrMapIIT.put(key, bhcorrIIT.get(l));
            bhcorrMapIITAllPep.put(key, bhcorrIITAllPep.get(l));

            l++;
        }

        List<Double> pval = new ArrayList<>();
        pval.addAll(pNameMap.values());
        List<Double> qval = BHCorrection.runBhCorrection(pval);
        int ii =0;
        for(String key: pNameMap.keySet())
        {
            pNameMap.put(key,qval.get(ii++));
        }

        normPvalues.addAll(normPvalueMap.values());
        List<Double> normqval = BHCorrection.runBhCorrection(normPvalues);
        ii=0;
        for(String key: normPvalueMap.keySet())
        {
            normQvalueMap.put(key,normqval.get(ii++));
        }

                /*
        Calculating Intensities Ratios QValue
         */

/*
        for(String key:nameIntensityRatioPMap.keySet())
        {
            List<Double> ratioIntensityPValues = new ArrayList<>();
            ratioIntensityPValues.addAll(nameIntensityRatioPMap.get(key));
            List<Double> ratioIntensityQValues = BHCorrection.runBhCorrection(ratioIntensityPValues);
            nameIntensityRatioQMap.put(key,ratioIntensityQValues);
        }
            */
/*
        List<List<Double>> pValueMatrix = new ArrayList<>();
        for(int i=0; i<indexList.size()-1; i++)
        {
            List<Double> temp = new ArrayList<>();
            for(String key:nameIntensityRatioPMap.keySet())
            {
                temp.add( nameIntensityRatioPMap.get(key).get(i));
            }
            pValueMatrix.add(temp);
        }
        List<List<Double>> qValueMatrix = new ArrayList<>();

        for(List<Double> dList: pValueMatrix)
        {
            qValueMatrix.add(BHCorrection.runBhCorrection(dList));
        }

        int i=0;
        for(String key:nameIntensityRatioPMap.keySet())
        {
           List<Double> temp = new ArrayList<>();

        }*/


/*
        List<List<Double>> intensityRatioQValueMatrix = new ArrayList<>();
        for(List<Double> column:intensityRatioPValueMatrix)
        {
            intensityRatioQValueMatrix.add(BHCorrection.runBhCorrection(column));
        }*/

        List<List<Double>> intensityRatioQValueMatrix = new ArrayList<>();
        List<List<Double>> normIntensityRatioQValueMatrix = new ArrayList<>();
        List<List<Double>> logIntensityRatioQValueMatrix = new ArrayList<>();
        List<List<Double>> logNormIntensityRatioQValueMatrix = new ArrayList<>();


        for(int i=0; i<intensityRatioPValueMatrix.size(); i++)
        {
            intensityRatioQValueMatrix.add(BHCorrection.runBHCorrectionRemoveRedundant(intensityRatioPValueMatrix.get(i)));
            normIntensityRatioQValueMatrix.add(BHCorrection.runBHCorrectionRemoveRedundant(normIntensityRatioPValueMatrix.get(i)));
            logIntensityRatioQValueMatrix.add(BHCorrection.runBHCorrectionRemoveRedundant(logIntensityRatioPValueMatrix.get(i)));
            logNormIntensityRatioQValueMatrix.add(BHCorrection.runBHCorrectionRemoveRedundant(logNormIntensityRatioPValueMatrix.get(i)));
        }




        ii =0;
        for(String name: nameList)
        {

            //for(List<Double> dlist: intensityRatioQValueMatrix)
            List<Double> intTemp = new ArrayList<>();
            List<Double> normTemp = new ArrayList<>();
            List<Double> logTemp = new ArrayList<>();
            List<Double> logNormTemp = new ArrayList<>();
            for(int i=0; i<intensityRatioQValueMatrix.size(); i++)
            {
                intTemp.add(intensityRatioQValueMatrix.get(i).get(ii));
                normTemp.add(normIntensityRatioQValueMatrix.get(i).get(ii));
                logTemp.add(logIntensityRatioQValueMatrix.get(i).get(ii));
                logNormTemp.add(logNormIntensityRatioQValueMatrix.get(i).get(ii));
            }
            nameIntensityRatioQMap.put(name,intTemp);
            nameNormRatioQMap.put(name,normTemp);
            comparedLogIntensityQMap.put(name,logTemp);
            comparedLogNormQMap.put(name,logNormTemp);
            ii++;
        }




/*
        Calculating Norm Intensities Ratios QValue

        for(String key:nameNormRatioPMap.keySet())
        {
            List<Double> normRatioIntensityPValues = new ArrayList<>();
            normRatioIntensityPValues.addAll(nameNormRatioPMap.get(key));
            List<Double> ratioIntensityQValues = BHCorrection.runBhCorrection(normRatioIntensityPValues);
            nameNormRatioQMap.put(key,ratioIntensityQValues);
        }

        /*
        Calculating Log Q Values

        for(String key:comparedLogIntensityPMap.keySet())
        {
            List<Double> comparedLogIntensityPvalues = comparedLogIntensityPMap.get(key);
            List<Double> comparedLogIntensityQvalues = BHCorrection.runBhCorrection(comparedLogIntensityPvalues);
            comparedLogIntensityQMap.put(key,comparedLogIntensityQvalues);
        }
        /*
        Calculating Log Norm Q Values

        for(String key:comparedLogNormPMap.keySet())
        {
            List<Double> comparedLogNormIntensityPvalues = comparedLogNormPMap.get(key);
            List<Double> comparedLogNormIntensityQvalues = BHCorrection.runBhCorrection(comparedLogNormIntensityPvalues);
            comparedLogNormQMap.put(key,comparedLogNormIntensityQvalues);
        }


        /*
        calculating avg of avg_intensity group wise
        */

        HashMap<String,List<Double>> avgofAverageIntensity = new HashMap<>();
        HashMap<String,List<Double>> avgofCorrAverageIntensity = new HashMap<>();
        for(String key:intensityMap.keySet()){
            List<Double> avgList1 = new ArrayList<>();
            List<Double> avgList2 = new ArrayList<>();
            int counter=0;
            for(int i=0;i<indexList.size();i++){
                double [] intensity =new double[indexList.get(i)];
                double [] corrIntensity =new double[indexList.get(i)];
                for(int j=0;j<indexList.get(i);j++){
                    intensity[j] = intensityMap.get(key).get(counter);
                    corrIntensity[j] = intensityCorrMap.get(key).get(counter);
                    counter++;
                }
                double sum1=0;
                double sum2=0;
                double avg1=0;
                double avg2=0;
                for(int k=0;k<intensity.length;k++){
                    sum1 = sum1 + intensity[k];
                    sum2 = sum2 + corrIntensity[k];
                }
                avg1 = sum1/intensity.length;
                avg2 = sum2/corrIntensity.length;
                avgList1.add(avg1);
                avgList2.add(avg2);

            }
            avgofAverageIntensity.put(key, avgList1);
            avgofCorrAverageIntensity.put(key, avgList2);
        }


        p.print("H\tNOTE\tratio p value will be calculated only when number of replicates are same between samples.\n");
        p.print("PLINE\tACCESSION\tINTENSITY_P-VALUE\tIIT_INTENSITY_P-VALUE\tIIT_ALL_PEP_INTENSITY_P-VALUE\tP-VALUE_RATIO_BASED\tINTENSITY_Q-VALUE"
                + "\tIIT_INTENSITY_Q-VALUE\tIIT_ALL_PEP_INTENSITY_Q-VALUE\tINTENSITY_Q-VALUE2\t");
        for (int i = 0; i < smallestIndex; i++) {

            p.print("MEDIAN_LOG_RATIO_" + (i + 1));
            p.print("\t");
        }
        for (int i = 0; i < sCountSize; i++) {
            p.print("SCOUNT_" + (i + 1));
            p.print("\t");
        }
        for (int i = 0; i < expSize; i++) {
            p.print("EACH_EXP_INTENSITY_SUM_" + (i + 1));
            p.print("\t");
        }
        for(int i =0;i<indexList.size();i++){
            p.print("AVG_GROUP_INTENSITY_VALUE_"+(i+1));
            p.print("\t");
        }
        for(int i =1;i<indexList.size();i++){
            p.print("AVG_GROUP_INTENSITY_RATIO_"+(i+1)+"_1");
            p.print("\t");
        }
        for(int i =1;i<indexList.size();i++){
            p.print("LOG_AVG_GROUP_INTENSITY_RATIO_"+(i+1)+"_1");
            p.print("\t");
        }

        for (int i = 0; i < expSize; i++) {
            p.print("AVG_ION_INJECTION_TIME_INTENSITY_" + (i + 1));
            p.print("\t");
        }
        for(int i =0;i<indexList.size();i++){
            p.print("AVG_GROUP_ION_INJECTION_TIME_INTENSITY_VALUE_"+(i+1));
            p.print("\t");
        }
        for(int i =1;i<indexList.size();i++){
            p.print("AVG_GROUP_IIT_INTENSITY_RATIO_"+(i+1)+"_1");
            p.print("\t");
        }
        for(int i =1;i<indexList.size();i++){
            p.print("LOG_AVG_GROUP_IIT_INTENSITY_RATIO_"+(i+1)+"_1");
            p.print("\t");
        }

        for (int i = 0; i < expSize; i++) {
            p.print("NORM_INTENSITY_" + (i + 1));
            p.print("\t");
        }
        for (int i = 0; i < expSize; i++) {
            p.print("LOG2_NORM_INTENSITY_" + (i + 1));
            p.print("\t");
        }

        for (int i = 0; i < indexList.size(); i++) {
            p.print("AVG_GROUP_NORM_INTENSITY_VALUE_" + (i + 1));
            p.print("\t");
        }
        for (int i =1; i < indexList.size(); i++) {
            p.print("AVG_GROUP_NORM_INTENSITY_RATIO_" + (i +1)+"_1");
            p.print("\t");
        }
        for (int i = 1; i < indexList.size(); i++) {
            p.print("LOG_AVG_GROUP_NORM_INTENSITY_RATIO_" + (i +1)+"_1");
            p.print("\t");
        }

        for(int i=0; i< indexList.size()-1; i++)
        {
            p.print("INTENSITY_RATIO_P-VALUE_"+(i+2)+"_1");
            p.print("\t");
        }

        for(int i=0; i< indexList.size()-1; i++)
        {
            p.print("INTENSITY_RATIO_Q-VALUE_"+(i+2)+"_1");
            p.print("\t");
        }

        for(int i=0; i< indexList.size()-1; i++)
        {
            p.print("NORM_INTENSITY_RATIO_P-VALUE_"+(i+2)+"_1");
            p.print("\t");
        }


        for(int i=0; i< indexList.size()-1; i++)
        {
            p.print("NORM_INTENSITY_RATIO_Q-VALUE_"+(i+2)+"_1");
            p.print("\t");
        }

        /*
            Log intensity p values compared with first
         */
        for(int i=0; i< indexList.size()-1; i++)
        {
            p.print("LOG_INTENSITY_P-VALUE_"+(i+2)+"_1");
            p.print("\t");
        }
        for(int i=0; i< indexList.size()-1; i++)
        {
            p.print("LOG_INTENSITY_Q-VALUE_"+(i+2)+"_1");
            p.print("\t");
        }

        for(int i=0; i< indexList.size()-1; i++)
        {
            p.print("LOG_NORM_INTENSITY_P-VALUE_"+(i+2)+"_1");
            p.print("\t");
        }

        for(int i=0; i< indexList.size()-1; i++)
        {
            p.print("LOG_NORM_INTENSITY_Q-VALUE_"+(i+2)+"_1");
            p.print("\t");
        }

        p.print("NORM_P-VALUE");
        p.print("\t");
        p.print("NORM_Q-VALUE");
        p.print("\t");
        p.print("TOTOL_PEPTIDE_INTENSITY\t");
        p.print("TOTAL_PEPTIDE_INTENSITY_DIVIDEDBY_PROTEIN_LENGTH_LOG10\t");
        p.print("PROTEIN_LENGTH\t");

        //p.print("PEP_COUNT\t");

        //for(int i=0;i<expSize;i++) {
        p.println("DESCRIPTION");

        /*
         p.print("SLINE\t");
         for(int i=0;i<expSize;i++) {
         p.print("EXP_" + (i+1));
         p.print("\t");
         p.print("SEQUENCE\tFILENAME\tSCAN\tCSTATE\tINTENSITY\tPROFILE_SCORE\tMHPLUS\tCALCMHPLUS\tTOTALINTENSITY\tXCORR\tDCN\tDMASS\tSPRANK\tSPSCORE\tREDUNDANCY\tSTARTRANGE\tENDRANGE\tRETENTIONTIME\tIONINJECTIONTIME\t");

	}*/
        for (Iterator<ProteinModel> itr = proteinList.iterator(); itr.hasNext();) {
            ProteinModel protein = itr.next();

            for (ChroProtein cprotein : protein.getRedundnatProteinList()) {
              p.print("P\t");
              p.print(cprotein.getLocus());
              p.print("\t");
              p.print(pvalueMap.get(cprotein.getLocus()));
              p.print("\t");
              //p.print(pvalueMap2.get(cprotein.getLocus()));
              p.print(pvalueMapIIT.get(cprotein.getLocus()));
              p.print("\t");
              p.print(pvalueMapIITAllPep.get(cprotein.getLocus()));
              p.print("\t");


              // p.print(protein.getMultipleRatiopValue());
              p.print("1.0\t");
              p.print(bhcorrMap.get(cprotein.getLocus()));
              p.print("\t");
              //p.print(bhcorrMap2.get(cprotein.getLocus()));
              p.print(bhcorrMapIIT.get(cprotein.getLocus()));
              p.print("\t");
              p.print(bhcorrMapIITAllPep.get(cprotein.getLocus()));
              p.print("\t");




              p.print(pNameMap.get(key2NameMap.get(cprotein.getLocus())));
              p.print("\t");

              for (double d : protein.getPeptideMedianLogRatioArr()) {
                p.print(d);
                p.print("\t");
              }
              for (Integer spc : cprotein.getSpecCountList()) {
                if (spc < 0) {
                  p.print(0);
                } else {
                  p.print(spc);
                }
                p.print("\t");
              }

                /*
                for (int i = 0; i < expSize; i++) {
                    List<Double> list = intensityMap.get(cprotein.getLocus());
                    if(list == null || list.isEmpty()){

                        p.print("NA\t");
                    }
                    else{
                       p.print(list.get(i).intValue());
                        p.print("\t");
                    }
                }*/
              //each_exp_intensity_sum

              //for(Long each:protein.getBestCorrelationEachpeptideIntensityList()) {
              for (Double each : protein.getIntensityEachSampleSumList()) {
                p.print(each.longValue());
                p.print("\t");
              }
              //AVG_GROUP_INTENSITY_VALUE_
              double[] intensityList = protein.getIntensitySumArr();
              List<Double> avgIntensityList = new ArrayList<>();
              for (int i = 0; i < indexList.size(); i++) {
                double avg = intensityList[i] / (double) indexList.get(i);
                avgIntensityList.add(avg);
                p.print(avg);
                p.print("\t");
              }



/*
              for(long each:protein.getIntensitySumArr()) {
                p.print(each +"\t");
              }*/
              /*
                for(Long each:protein.getIntensitySumArr())
                    List<Long> list = protein.getBestCorrelationIntensityList();
                    p.print(list.get(i) + "\t");
                }
                */
              //AVG_GROUP_INTENSITY_RATIO_2_1
              for (int i = 1; i < avgIntensityList.size(); i++) {

                //List<Long> list = protein.getBestCorrelationIntensityList();
                double ratio = (double) avgIntensityList.get(i) / avgIntensityList.get(0);
                if (Double.isInfinite(ratio)) ratio = 1000;
                p.print(ratio);
                p.print("\t");

              }
              //LOG_AVG_GROUP_INTENSITY_RATIO_2_1
              for (int i = 1; i < avgIntensityList.size(); i++) {

                //List<Long> list = protein.getBestCorrelationIntensityList();
                double ratio = (double) avgIntensityList.get(i) / avgIntensityList.get(0);
                if (Double.isInfinite(ratio)) ratio = 1000;

                String logRatio = "NA";
                try {
                  double lr = Math.log(ratio) / Math.log(2);
                  if (Double.isInfinite(lr)) lr = 10E-3;

                  p.print(lr);
                } catch (Exception e) {
                  p.print(logRatio);
                }

                p.print("\t");

              }
              //avg_ion_injection_time_intensity


              //  for (double each:protein.getIonInjectionTimeIntensitySumList()) {
              if (protein.getIonInjectionTimeIntensitySumList().size() > 0) {
                for (double each : protein.getIonInjectionTimeIntensitySumList()) {


                  p.print(each);
                  p.print("\t");
                }
              } else {
                for(int i=0;i<expSize;i++) {
                  p.print(0);
                  p.print("\t");
                }

              }
  /*
                for(long each:protein.getBestCorrelationEachpeptideIntensityListIIT()) {
                    p.print(each);
                    p.print("\t");

                }*/
                //AVG_GROUP_ION_INJECTION_TIME_INTENSITY_VALUE_
              if (protein.getIonInjectionTimeIntensityAvgGroupList().size() > 0) {
                for (double each : protein.getIonInjectionTimeIntensityAvgGroupList()) {
                  p.print(each);
                  p.print("\t");
                }
              } else {
                for(int i=0;i<expSize;i++) {
                  p.print(0);
                  p.print("\t");
                }
              }


              if(protein.getIonInjectionTimeIntensityAvgGroupList().size()>0) {
                //AVG_GROUP_IIT_INTENSITY_RATIO_2_1
                for (int i = 1; i < protein.getIonInjectionTimeIntensityAvgGroupList().size(); i++) {
                  List<Double> list = protein.getIonInjectionTimeIntensityAvgGroupList();

                  if (list.get(0) <= 0) p.print(100.0);
                  else if (list.get(i) <= 0) p.print(-100.0);
                  else {
                    double ratio = (double) list.get(i) / list.get(0);
                    p.print(ratio);
                  }


                  //System.out.println(">>----" + ratio);
                  p.print("\t");

                }
              } else {
                for(int i=0;i<expSize;i++) {
                  p.print(0);
                  p.print("\t");
                }
              }

              /*
              if(protein.getBestCorrelationIntensityListIIT().size()>0) {
                //AVG_GROUP_IIT_INTENSITY_RATIO_2_1
                for (int i = 1; i < protein.getBestCorrelationIntensityListIIT().size(); i++) {
                  List<Long> list = protein.getBestCorrelationIntensityListIIT();

                  if (list.get(0) <= 0) p.print(100.0);
                  else if (list.get(i) <= 0) p.print(-100.0);
                  else {
                    double ratio = (double) list.get(i) / list.get(0);
                    p.print(ratio);
                  }


                  //System.out.println(">>----" + ratio);
                  p.print("==========\t");

                }
              } else {
                for(int i=0;i<expSize;i++) {
                  p.print(0);
                  p.print("--------\t");
                }
              }
*/
                //LOG_AVG_GROUP_IIT_INTENSITY_RATIO_2_1
                for (int i = 1; i < indexList.size(); i++) {
                    List<Double> list = protein.getIonInjectionTimeIntensityAvgGroupList();

                    String logRatio = "NA";
                    try {

                        double ratio = (double) list.get(i) / list.get(0);//double ratio = list.get(i)/list.get(0);
                        if (list.get(0) <= 0) p.print(100.0);
                        else if (list.get(i) <= 0) p.print(-100.0);
                        else {
                            double lr = Math.log(ratio) / Math.log(2);
                            if(Double.isInfinite(lr)) lr = 10E-3;
                            //System.out.println("---->>" + list.get(i) + "\t" + list.get(0));
                            //System.out.println("----" + ratio + "\t" + lr);
                            p.print(lr);
                        }
                    } catch (Exception e) {
                        p.print(logRatio);
                    }

                    p.print("\t");
                }
                //norm_intensity
                for(int i=0;i<protein.getNormIntensityList().size();i++){
                    p.print(protein.getNormIntensityList().get(i).longValue()+"\t");
                }
                //norm_intensity
                for(int i=0;i<protein.getNormIntensityList().size();i++){
                    double intensity = protein.getNormIntensityList().get(i);
                    double logInt = Math.log(intensity)/Math.log(2);
                    if(Double.isInfinite(logInt) || Double.isNaN(logInt)) logInt = 10;
                    p.print(logInt+"\t");
                }
                //avg_norm_intensity_value
                for(double d:avgNormIntensityMap.get(cprotein.getLocus()))
                {
                    p.print(d+"\t");
                }

                //avg_norm_intesity_ratio
                for(int i=1; i<avgNormIntensityMap.get(cprotein.getLocus()).size();i++)
                {
                    double ratio = avgNormIntensityMap.get(cprotein.getLocus()).get(i)/avgNormIntensityMap.get(cprotein.getLocus()).get(0);
                    if(Double.isInfinite(ratio)) ratio=1000;
                    p.print(ratio+"\t");
                }
                //log_avg_norm_intesity_ratio
                for(int i=1; i<avgNormIntensityMap.get(cprotein.getLocus()).size();i++) {
                    double ratio = avgNormIntensityMap.get(cprotein.getLocus()).get(i) / avgNormIntensityMap.get(cprotein.getLocus()).get(0);
                    if (Double.isInfinite(ratio)) ratio = 1000;
                    String logRatio = "NA";
                    try {
                        double lr = Math.log(ratio) / Math.log(2);
                        if(Double.isInfinite(lr)) lr = 10E-3;
                        p.print(lr);
                    } catch (Exception e) {
                        p.print(logRatio);
                    }

                    p.print("\t");
                }

                //ratio intensity pvalue
                for(double d:nameIntensityRatioPMap.get(key2NameMap.get(cprotein.getLocus())) )
                {
                    p.print(d+"\t");
                }
                for(double d:nameIntensityRatioQMap.get(key2NameMap.get(cprotein.getLocus())) )
                {
                    p.print(d+"\t");
                }
                //NORM_INTENSITY_RATIO_P-VALUE

                for(double d:nameNormRatioPMap.get(key2NameMap.get(cprotein.getLocus())) )
                {
                    p.print(d+"\t");
                }


                for(double d:nameNormRatioQMap.get(key2NameMap.get(cprotein.getLocus())) )
                {
                    p.print(d+"\t");
                }


                for(double d:comparedLogIntensityPMap.get(key2NameMap.get(cprotein.getLocus())) )
                {
                    p.print(d+"\t");
                }
                for(double d:comparedLogIntensityQMap.get(key2NameMap.get(cprotein.getLocus())) )
                {
                    p.print(d+"\t");
                }

                for(double d:comparedLogNormPMap.get(key2NameMap.get(cprotein.getLocus())) )
                {
                    p.print(d+"\t");
                }
                for(double d:comparedLogNormQMap.get(key2NameMap.get(cprotein.getLocus())) )
                {
                    p.print(d+"\t");
                }


                p.print(normPvalueMap.get(key2NameMap.get(cprotein.getLocus())));
                p.print("\t");
                p.print(normQvalueMap.get(key2NameMap.get(cprotein.getLocus())));
                p.print("\t");
                p.print(protein.getTotalPeptideIntensity());
                p.print("\t");
                //totol_peptide_intensity_dividedby_protein_length_log10
                p.print(protein.getTotalIntensityDividedByProteinLengthLog10());
                p.print("\t");
                p.print(protein.getProteinLength());
                p.print("\t");

                p.println(cprotein.getDescription());

            }
        }

        p.close();
        System.out.println("check file " + filename);

    }

    //public static List<ProteinModel> runLabelfree(
    public static void runLabelfree(
            String configFile, String tmpFile, String jsonFile,
            List<Integer> sample1IndexList,
            List<Integer> sample2IndexList) throws Exception {

        Configuration conf = Configuration.getInstance();

        if (!conf.isReadConfigFile()) {
            conf.setLabelfree(true);
            conf.readXMLParam(configFile);
        }


        TObjectIntHashMap proteinLengthMap = LabelfreeMissingPeptideBuilderSplit.getProteinLengthTMap(conf);
        System.out.println(">>>> finished reading protein length file");
        List<org.jdom.Element> samGroupEleList1 = conf.getRootConfEle().getChildren("sample");
        int expSize = 0;

        List<Integer> indexList = new ArrayList<>();
        for (Iterator<org.jdom.Element> samgItr = samGroupEleList1.iterator(); samgItr.hasNext();) {
            int count =0;
            org.jdom.Element groupEle = samgItr.next();


            List<org.jdom.Element> sampleEleList = groupEle.getChildren("each_sample");

            for (Iterator<org.jdom.Element> samItr = sampleEleList.iterator(); samItr.hasNext();) {
                org.jdom.Element eachSample = samItr.next();
                count++;
                expSize++;
            }
            indexList.add(count);
        }


        IsotopeReader isoReader = new IsotopeReader(conf.getRootConfEle());
        //IdentificationReader idReader = BaseIdentificationReader.getIdentificationInst(isoReader);

        conf.setSpectrumFormat(Configuration.MS_FILE_FORMAT);

        Hashtable<String, ChroPeptide> peptideChroHt = new Hashtable<>();

        List<SampleModel> sampleList = conf.getSampleList();
        // sampleList.get(0).getPathList();
        Hashtable<String, IndexedFile> origMs1FileHt = new Hashtable<>();

        Map<String, HashMap<Integer, Integer>> ms2ToMs1Map = new HashMap<>();

        int smallestIndex = sample1IndexList.size()>sample2IndexList.size() ? sample2IndexList.size() : sample1IndexList.size();

        for (SampleModel eachSample : sampleList) {
            for (String path : eachSample.getPathList()) {

                if(!path.endsWith("/"))
                    path += "/";

                String spectraPath = path + "../../spectra/";
                String splitSpectraPath = path + "../../spectra/split/";

             //   splitSpectraMap.putAll( msp.splitMS1Files(spectraPath, CensusConstants.LABELFREE_MS1_SPLIT_SCAN_NUM, true, conf.isLabelfreeCheckChargeState()) );
             //   splitMs1FileHt.putAll( ChroGenerator.createIndexedFiles(splitSpectraPath, CensusConstants.MS1_FILE) );
                ms2ToMs1Map.putAll( IndexUtil.buildMS2toMS1ScanMapFiles(spectraPath) );

                Hashtable<String, IndexedFile> ht = ChroGenerator.createIndexedFiles(spectraPath, "ms1", true,true);
                origMs1FileHt.putAll(ht);

                //Hashtable<String, IndexedFile> ht = ChroGenerator.createIndexedFiles(path, "ms1", false);
                //indexHt.putAll(ht);

                ChroXmlReader cr = new ChroXmlReader(path + "census_chro_temp.xml");
                ArrayList<ChroProtein> list = cr.getProteinList();
                for (ChroProtein pro : list) {
                    List<ChroPeptide> pepList = pro.getPeptideList();
                    for (ChroPeptide pep : pepList) {
                        //String chroData = pep.getChroData();
                        String key = pep.getFileName() + pep.getSequence() + pep.getChargeState() + pep.getScanNum();

                        peptideChroHt.put(key, pep);

                    }
                }
            }
        }

        //      if(true) System.exit(0);
        conf.setIndexHt(origMs1FileHt);
        TxttmpReader txtReader = new TxttmpReader(tmpFile);
        List<ProteinModel> proteinList = txtReader.readWholeFile();
        int totalProtein = proteinList.size();
        Hashtable<String, ChroPeptide> analyzedPeptideHt = new Hashtable<String, ChroPeptide>();
        int proteinCount = 0;

        ChroJSONCreator chroJSONCreator = new ChroJSONCreator();
        String jsonFilePath = new File(jsonFile).getParent() + File.separator + "JSON_OBJ";
        FileUtil.makeDir(jsonFilePath);
        List<LabelFreeJSONProtein> jsonProteins = new ArrayList<>();
        LabelFreeJSONProtein jsonProtein = null;
        StringBuffer sb = new StringBuffer();
        StringBuffer gPeakSb = new StringBuffer();
        for (Iterator<ProteinModel> itr = proteinList.iterator(); itr.hasNext();) {
            ProteinModel p = itr.next();
            proteinCount++;

            // List<String> keys = p.getPeptideKey();
            //  HashMap<String, List<ChroPeptide>> peptideMap = p.getPeptideMap();

            List[] ratioArr = new List[smallestIndex];
            for (int i = 0; i < ratioArr.length; i++) {
                ratioArr[i] = new ArrayList<Double>();
            }

            jsonProtein = new LabelFreeJSONProtein();
            jsonProtein.setAccession(p.getRedundnatProteinList().get(0).getLocus());
            jsonProtein.setDesc(p.getRedundnatProteinList().get(0).getDescription());
            jsonProteins.add(jsonProtein);

            List<List<LabelFreeJSONPeptide>> jsonAllPeptideList = new ArrayList<>();
            String jsonPeptideListFileName = jsonFilePath + File.separator + p.getRedundnatProteinList().get(0).getLocus() + ".JSON";
            int sampleCount = 0;
            int peptideCount = 0;

            for (LabelfreePeptide each : p.getPeptideList()) {

                double startRt = each.getStartRetTime();
                double endRt = each.getEndRetTime();

                //  System.out.println(each.getPeptideList());
                int count = 0;
                for (ChroPeptide expPep : each.getPeptideList()) {


                 //   if("expPep.getSequence().equals(  )expPep.getChargeState() == 2)

             //       if(!"R.FKDLGEENFK.A".equals(expPep.getSequence()) || (expPep.getChargeState() != 2))
               //         continue;

             //       System.out.println("===\t" + expPep.getSequence() + " " +  expPep.getChargeState() + " " + expPep.getFileName());

                    //not identified missing peptide
                    if (expPep.getFileName() == null && conf.getLabelfreeFillPeptide().equals("true")) { //not identified peptide

                        String key = each.getSequence() + each.getChargeState() + count;
                        ChroPeptide tmpPep = analyzedPeptideHt.get(key);
                        if (null != tmpPep) {
                            expPep.setPeakArea(tmpPep.getPeakArea());
                            expPep.setChroData(tmpPep.getChroData());

                            expPep.setGaussianPeakString(tmpPep.getGaussianPeakString());
                            expPep.setPeakSigma(tmpPep.getPeakSigma());
                            expPep.setPeakx(tmpPep.getPeakx());
                            expPep.setPeaky(tmpPep.getPeaky());

                            // System.out.println(key);
                            // System.out.println(tmpPep.getChroData());
                            count++;

                            continue;
                        }

//                        String fname = .getLabelfreeFilename();
                        List<String> fnameList = sampleList.get(count).getLabelfreeFilenameList();
                        List<String> pathList = sampleList.get(count).getPathList();

                        double peakArea = 0;
                        for(int i=0;i<fnameList.size();i++) {
                            String eachFile = fnameList.get(i);
                            String eachPath = pathList.get(i);
                            String eachKey;

                            if (eachPath.endsWith("/")) {
                                eachKey = eachPath + "../../spectra/" + eachFile;
                            } else {
                                eachKey = eachPath + "/../../spectra/" + eachFile;
                            }

                            IndexedFile origIFile = origMs1FileHt.get(eachKey);
                            TDoubleIntHashMap retentonToScanMap = origIFile.getRetentonToScanMap();

                            int startScan = retentonToScanMap.get(startRt);
                            int endScan = retentonToScanMap.get(endRt);

                            double[] retKeys = origIFile.getRtArr();
                            int startIndex = BinarySearch.binarySearch(retKeys, startRt);

                            if (startScan <= 0) {

                                double rtTime = retKeys[startIndex];

                                startScan = retentonToScanMap.get(rtTime);
                            }
                            int endIndex = BinarySearch.binarySearch(retKeys, endRt);
                            if (endScan <= 0) {

                                double rtTime = retKeys[endIndex];
                                endScan = retentonToScanMap.get(rtTime);
                            }

                            //GaussianPeakModel peakModel = isotopeCalc(startScan, endScan, startIndex, endIndex, isoReader, each.getSequence(), each.getChargeState(), iFile);
                            GaussianPeakModel peakModel = isotopeCalc( startIndex, endIndex, isoReader, each.getSequence(),
                                    each.getChargeState(), origIFile, false);

                            if (null != peakModel) {
                                peakArea += peakModel.getPeakArea();
                               // expPep.setPeakArea(peakModel.getPeakArea());

                                int[] scanArr = peakModel.getScanArr();
                                double[] retArr = peakModel.getRetArr();
                                double[] peakArr = peakModel.getPeakArr();

                                sb.delete(0, sb.length());
                                sb.append("P 0 0e2;");
                                for (int j = 0; j < scanArr.length; j++) {
                                    sb.append(scanArr[j]).append(" ").append(retArr[j]).append(" ").append(peakArr[j]).append(";");
                                }

                                double[] gxArr = peakModel.getGaussianXArr();
                                double[] gyArr = peakModel.getGaussianYArr();
                                gPeakSb.delete(0, gPeakSb.length());


//                                System.out.println(java.util.Arrays.toString(gyArr));
                                if (null != gxArr) {
                                    for (int j = 0; j < gxArr.length; j++) {
                                        gPeakSb.append(gxArr[j]).append(" ").append(gyArr[j]).append(";");
                                    }
                                }

                                expPep.setChroData(sb.toString());
                                expPep.setGaussianPeakString(gPeakSb.toString());
                                expPep.setPeakSigma(peakModel.getSigma());
                                expPep.setPeakx(peakModel.getX());
                                expPep.setPeaky(peakModel.getY());
                            }

                            if(peakArea>0) //sum of peak area if there are more than one spectral files
                            expPep.setPeakArea(peakArea);

                            ChroPeptide tmpChroPeptide = analyzedPeptideHt.get(key);

                            if( !peakModel.isHasPeak() && null!=tmpChroPeptide)
                                continue;
                        }


                        analyzedPeptideHt.put(key, expPep);

                    } else if (expPep.getFileName() == null && conf.getLabelfreeFillPeptide().equals("fill_zero_only")) {


                        String key = each.getSequence() + each.getChargeState() + count;
                        ChroPeptide tmpPep = analyzedPeptideHt.get(key);
                        if (null != tmpPep) {
                            expPep.setPeakArea(tmpPep.getPeakArea());
                            expPep.setChroData(tmpPep.getChroData());

                            expPep.setGaussianPeakString(tmpPep.getGaussianPeakString());
                            expPep.setPeakSigma(tmpPep.getPeakSigma());
                            expPep.setPeakx(tmpPep.getPeakx());
                            expPep.setPeaky(tmpPep.getPeaky());

                            // System.out.println(key);
                            // System.out.println(tmpPep.getChroData());
                            count++;

                            continue;
                        }

//                        String fname = .getLabelfreeFilename();
                        List<String> fnameList = sampleList.get(count).getLabelfreeFilenameList();
                        List<String> pathList = sampleList.get(count).getPathList();

                        double peakArea = 0;
                        for(int i=0;i<fnameList.size();i++) {
                            String eachFile = fnameList.get(i);
                            String eachPath = pathList.get(i);
                            String eachKey;

                            if (eachPath.endsWith("/")) {
                                eachKey = eachPath + "../../spectra/" + eachFile;
                            } else {
                                eachKey = eachPath + "/../../spectra/" + eachFile;
                            }

                            IndexedFile origIFile = origMs1FileHt.get(eachKey);
                            TDoubleIntHashMap retentonToScanMap = origIFile.getRetentonToScanMap();

                            int startScan = retentonToScanMap.get(startRt);
                            int endScan = retentonToScanMap.get(endRt);

                            double[] retKeys = origIFile.getRtArr();
                            int startIndex = BinarySearch.binarySearch(retKeys, startRt);

                            if (startScan <= 0) {

                                double rtTime = retKeys[startIndex];

                                startScan = retentonToScanMap.get(rtTime);
                            }
                            int endIndex = BinarySearch.binarySearch(retKeys, endRt);
                            if (endScan <= 0) {

                                double rtTime = retKeys[endIndex];
                                endScan = retentonToScanMap.get(rtTime);
                            }

                            //GaussianPeakModel peakModel = isotopeCalc(startScan, endScan, startIndex, endIndex, isoReader, each.getSequence(), each.getChargeState(), iFile);
                            GaussianPeakModel peakModel = isotopeCalc( startIndex, endIndex, isoReader, each.getSequence(),
                                    each.getChargeState(), origIFile, true
                            );

                            if (null != peakModel) {
                                peakArea += peakModel.getPeakArea();
                                // expPep.setPeakArea(peakModel.getPeakArea());

                                int[] scanArr = peakModel.getScanArr();
                                double[] retArr = peakModel.getRetArr();
                                double[] peakArr = peakModel.getPeakArr();

                                sb.delete(0, sb.length());
                                sb.append("P 0 0e2;");
                                for (int j = 0; j < scanArr.length; j++) {
                                    sb.append(scanArr[j]).append(" ").append(retArr[j]).append(" ").append(peakArr[j]).append(";");
                                }

                                double[] gxArr = peakModel.getGaussianXArr();
                                double[] gyArr = peakModel.getGaussianYArr();
                                gPeakSb.delete(0, gPeakSb.length());


//                                System.out.println(java.util.Arrays.toString(gyArr));
                                if (null != gxArr) {
                                    for (int j = 0; j < gxArr.length; j++) {
                                        gPeakSb.append(gxArr[j]).append(" ").append(gyArr[j]).append(";");
                                    }
                                }

                                expPep.setChroData(sb.toString());
                                expPep.setGaussianPeakString(gPeakSb.toString());
                                expPep.setPeakSigma(peakModel.getSigma());
                                expPep.setPeakx(peakModel.getX());
                                expPep.setPeaky(peakModel.getY());
                            }

                            if(peakArea>0) //sum of peak area if there are more than one spectral files
                                expPep.setPeakArea(peakArea);

                            ChroPeptide tmpChroPeptide = analyzedPeptideHt.get(key);

                            if( !peakModel.isHasPeak() && null!=tmpChroPeptide)
                                continue;
                        }


                        analyzedPeptideHt.put(key, expPep);


                    } else {

                        expPep.getFileName();
                        String tkey = expPep.getFileName() + expPep.getSequence() + expPep.getChargeState() + expPep.getScanNum();
                        ChroPeptide cPep = peptideChroHt.get(tkey);

                        if (cPep != null) {
                            expPep.setChroData(cPep.getChroData());
                            expPep.setPeakSigma(cPep.getPeakSigma());
                            expPep.setPeakx(cPep.getPeakx());
                            expPep.setPeaky(cPep.getPeaky());
                            expPep.setGaussianPeakString(cPep.getGaussianPeakString());
                        }

                    }
                    count++;

                    //  System.out.print("---" + expPep.getPeakArea()+ "\t");
                }

            //    System.out.println("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa");

                List<ChroPeptide> pepList = each.getPeptideList();
                //  int rcount=0;

                for (int i =0; i<smallestIndex; i++) {
                    double peakArea1 = pepList.get(sample1IndexList.get(i)).getPeakArea();
                    double peakArea2 = pepList.get(sample2IndexList.get(i)).getPeakArea();
                    double ratio;
                    if(Double.compare(peakArea2,peakArea1)==0 && Double.compare(peakArea1,0.0)==0) ratio =0;
                    else if(Double.compare(peakArea2,0.0)==0) ratio = 10;
                    else if(Double.compare(peakArea1,0.0)==0) ratio = -10;
                     else ratio = Math.log(peakArea1 / peakArea2) / Math.log(2);
                    //  double ratio = peakArea1/peakArea2;
                    //   System.out.print(peakArea1  + "\t" + peakArea2 + "\t");

                    ratioArr[i].add(ratio);
                    //   System.out.println(ratio + "\t");
                }
                /*
                for (int index : sample1IndexList) {
                    double peakArea1 = pepList.get(index).getPeakArea();
                    double peakArea2 = pepList.get(index + sample1IndexList.size()).getPeakArea();
                    if(Double.compare(peakArea2,0.0)==0) continue;

                    double ratio = Math.log(peakArea1 / peakArea2) / Math.log(2);
                    //  double ratio = peakArea1/peakArea2;
                    //   System.out.print(peakArea1  + "\t" + peakArea2 + "\t");

                    ratioArr[index].add(ratio);
                    //   System.out.println(ratio + "\t");
                }
*/
                /*
                 * Json - Part
                 */
                sampleCount++;
                List<LabelFreeJSONPeptide> jsonPeptides = getEachSamplePeptides(pepList, sampleCount, peptideCount, each.getChargeState(), each.getSequence(), startRt, endRt);
                peptideCount = peptideCount + jsonPeptides.size();
                jsonAllPeptideList.add(jsonPeptides);

                //         System.out.println("");
            }

            p.setLogRatioArr(ratioArr);

            /*
             * Each Json peptide list wiring - filename is locus
             */
            chroJSONCreator.createJsonForPeptideList(jsonAllPeptideList, jsonPeptideListFileName);

            System.out.println(proteinCount + "\t" + totalProtein + "\t" + ((double) proteinCount) / totalProtein * 100 + "% complete");

        }



        /*
         * Each Json ProteinIndex file  wiring
         */
        chroJSONCreator.createJsonForProteinIndex(jsonProteins, jsonFile);

        //      System.out.println("finding missing peptides are done");
        // adding ANOVA pValue:
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

        ctr = 0;
        Iterator f = proteinList.iterator();
        List<Double> pvalueList = new ArrayList<>();

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

            List<Long> pepAvgIntList = null;
            List<Long> pepIntList = null;
            List<Long> pepAvgIntListIIT = null;
            List<Long> pepIntListIIT = null;
            double[] intensitySumArr = new double[samGroupEleList.size()];
            List<Double> ionInjectionTimeIntensitySumList =new ArrayList<>();
          List<Double> intensitySumList =new ArrayList<>(); //intensity sum for each protein e.g. 6 experiments -> 6 intensity sum list
            double totalPeptideIntensity=0;

            for (ctr = 0; ctr < size1; ctr++) {

               // int size2 = protein.getPeptideList().get(ctr).getPeptideList().size();
                int size2 = expSize;
                if(ctr==0) {
                  ionInjectionTimeIntensitySumList.addAll(Collections.nCopies(size2,new Double(0)));
                  intensitySumList.addAll(Collections.nCopies(size2,new Double(0)));
                }
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
                    double intSum = 0;
                    double intSumIIT=0;
                    int oldT =0;
                    //while (iterator.hasNext()) {
                    for(int i=0;i<intensity.size();i++) {

                        intensityy[t] = intensity.get(i);
                        intensitySumArr[sampleCount] += intensityy[t];

                        intSum += (long)intensityy[t];
                        pepIntList.add((long)intensityy[t]);


                        iitIntensityy[t] = iitIntensity.get(i);
                        intSumIIT += (long)iitIntensityy[t];
                        pepIntListIIT.add((long)iitIntensityy[t]);
                        t++;
                    }


                    long averageIntensity = (long)intSum/t;
                    long averageIntensityIIT = (long)intSumIIT/t;

                    pepAvgIntList.add(averageIntensity);
                    pepAvgIntListIIT.add(averageIntensityIIT);
                    totalPeptideIntensity += intSum;
                    //System.out.println(totalPeptideIntensity + "======---->" + intSum);
                    if(intensityy.length<=1)
                        hasReplicates=false;
                /*    if(classes.size()>0 && oldT !=t)
                        hasReplicates = false;*/
                    oldT = t;

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
                    intensitySumList.set(i, intensitySumList.get(i) + pepIntList.get(i));
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
            protein.setIntensityEachSampleSumList(intensitySumList);


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
                protein.setIonInjectionTimeIntensitySumList(ionInjectionTimeIntensitySumList);
                //avg_ion_injection_time_intensity calc

            }

            protein.setTotalPeptideIntensity( totalPeptideIntensity );

            Integer lengthObj = proteinLengthMap.get(protein.getRedundnatProteinList().get(0).getLocus());
            int proteinLength = 1;
            if(null != lengthObj)
              proteinLength = lengthObj.intValue();


            protein.setTotalIntensityDividedByProteinLengthLog10( Math.log10(totalPeptideIntensity/(double)proteinLength) );
            //totol_peptide_intensity_dividedby_length_log10
            protein.setProteinLength(proteinLength);

            intList.clear();
        }
        //adding norm intensity
        double [] sumIntensityPeptide = null;
        double [] sumIntensityProtein = null;
        boolean check1 = true;
        boolean check2 = true;

        for(int i=0;i<proteinList.size();i++){

            ProteinModel protein  = proteinList.get(i);
            if(check2){
               // sumIntensityProtein = new double[protein.getBestCorrelationEachpeptideIntensityList().size()];
                sumIntensityProtein = new double[protein.getIntensityEachSampleSumList().size()];
                check2=false;
            }
            for(int n=0;n<protein.getIntensityEachSampleSumList().size();n++){
                double sum = sumIntensityProtein[n];
                sum = sum+((double)protein.getIntensityEachSampleSumList().get(n));
                sumIntensityProtein[n] = sum;
            }
            List<LabelfreePeptide> pepList =protein.getPeptideList();
            for(int j=0;j<pepList.size();j++){
                LabelfreePeptide pep = pepList.get(j);
                List<ChroPeptide> expList = pep.getPeptideList();
                for(int k=0;k<expList.size();k++){
                    if(check1){
                        sumIntensityPeptide=new double[expList.size()];
                        check1=false;
                    }

                    ChroPeptide exp = expList.get(k);
                    double sum = sumIntensityPeptide[k];
                    sum = sum+(exp.getAverageIntensity());
                    sumIntensityPeptide[k]=sum;
                }
            }
        }
        double totalsumproteinIntensity =0.0;
        for(int i=0;i<sumIntensityProtein.length;i++){
            totalsumproteinIntensity = totalsumproteinIntensity+sumIntensityProtein[i];
        }
        double totalsumpeptideIntensity =0.0;
        for(int i=0;i<sumIntensityProtein.length;i++){
            totalsumpeptideIntensity = totalsumpeptideIntensity+sumIntensityPeptide[i];
        }
        double avgproteinIntensity = totalsumproteinIntensity/sumIntensityProtein.length;
        double avgpeptideIntensity = totalsumpeptideIntensity/sumIntensityPeptide.length;


        //adding to the model list
        for(int i=0;i<proteinList.size();i++){

            ProteinModel protein  = proteinList.get(i);

            for(int n=0;n<protein.getBestCorrelationEachpeptideIntensityList().size();n++){
                //protein.addNormIntensityList(protein.getBestCorrelationEachpeptideIntensityList().get(n)*(avgpeptideIntensity/sumIntensityProtein[n]));
                protein.addNormIntensityList(protein.getIntensityEachSampleSumList().get(n)*avgpeptideIntensity/sumIntensityProtein[n]);
            }
            List<LabelfreePeptide> pepList =protein.getPeptideList();
            for(int j=0;j<pepList.size();j++){
                LabelfreePeptide pep = pepList.get(j);
                List<ChroPeptide> expList = pep.getPeptideList();
                for(int k=0;k<expList.size();k++){
                    ChroPeptide exp = expList.get(k);
                    pep.addNormIntensityList(exp.getAverageIntensity()*(avgpeptideIntensity/sumIntensityPeptide[k]));
                    //norm_intensity add
                }
            }
        }




        String filledFile = tmpFile.substring(0, tmpFile.lastIndexOf(".")) + "_filled.txt";

        BufferedReader br = new BufferedReader(new FileReader(tmpFile));

        String eachLine;
        int i = -1;
        int k = 0;
        int l = 0;
        int r = 0;
        //print output file
        BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(filledFile));
        PrintStream ps = new PrintStream(out);
        int normIndex =-1;
        while (null != (eachLine = br.readLine()) && !eachLine.startsWith("SLINE\t")) {

            if(eachLine.startsWith("PLINE\t")){

              ps.println("H\tPLINE COLUMN DESCRIPTION\tACCESSION\tprotein accession");
              ps.println("H\tPLINE COLUMN DESCRIPTION\tDESCRIPTION\tprotein description");
              ps.println("H\tPLINE COLUMN DESCRIPTION\tSCOUNT_X\tspectral count");
              ps.println("H\tPLINE COLUMN DESCRIPTION\tPEP_COUNT\tpeptide count");
              ps.println("H\tPLINE COLUMN DESCRIPTION\tNORM_INTENSITY_X\tintensity normalized by total intensity for each experiment");
              ps.println("H\tPLINE COLUMN DESCRIPTION\tINTENSITY_X\tprotein intensity derived from peptide intensity sum for each experiment");
              ps.println("H\tPLINE COLUMN DESCRIPTION\tGROUP_INTENSITY_SUM_X\tintensity sum of replicates for protein");

                String [] words = eachLine.split("\t");
                for(int rr=0;rr<words.length;rr++){
                    if(words[rr].startsWith("NORM")){
                        normIndex=rr;
                        break;
                    }
                }

                ps.print(eachLine);

                for (int a = 0; a < smallestIndex; a++) {

                    ps.print("MEDIAN_LOG_RATIO_" + (a + 1));
                    ps.print("\t");

                }
                //for (int a = 0; a < expSize; a++) {
                //    ps.print("BEST_PEP_INTENSITY_"+(a+1)+"\t");
                //}
                /*
               //put header for total protein intensity here



                 */

                for (int a = 0; a < expSize; a++) {
                    ps.print("INTENSITY_"+(a+1)+"\t");

                }
                for(int a=0;a<indexList.size();a++){
                    ps.print("GROUP_INTENSITY_SUM_"+(a+1)+"\t");

                }
                for(int a=0;a<expSize;a++){
                    ps.print("BEST_PEP_ION_INJECTION_TIME_INTENSITY_"+(a+1)+"\t");

                }
                for(int a=0;a<indexList.size();a++){
                    ps.print("BEST_PEP_GROUP_ION_INJECTION_TIME_INTENSITY_VALUE_"+(a+1)+"\t");

                }
                //modified lines/
                /*
                for(int a=0;a<indexList.size();a++){
                    ps.print("AVG_GROUP_INTENSITY_VALUE_"+(a+1)+"\t");
                }
                for(int a=0;a<indexList.size();a++){
                    ps.print("AVG_ION_INJECTION_TIME_INTENSITY_"+(a+1)+"\t");
                }
                for(int a=0;a<indexList.size();a++){
                    ps.print("AVG_GROUP_ION_INJECTION_TIME_INTENSITY_VALUE_"+(a+1)+"\t");
                }
                //*/
                ps.println();
            }
            else{
                ps.println(eachLine);
            }

        }
        int normIndexS =-1;
        String [] words3 = eachLine.split("\t");
        for(int rr=0;rr<words3.length;rr++){
            if(words3[rr].startsWith("NORM")) {
                normIndexS=rr;
                break;
            }
        }
        ps.println(eachLine + "PVALUE");
        while ((eachLine = br.readLine()) != null) {

            if (eachLine.startsWith("P\t")) {
                i++;

                String [] word = eachLine.split("\t");
                for(int rr=0;rr<normIndex;rr++) {
                    ps.print(word[rr]+"\t");
                }

                for(int rr=0;rr<proteinList.get(i).getNormIntensityList().size();rr++){
                    ps.print(proteinList.get(i).getNormIntensityList().get(rr).longValue()+"\t");

                }
                for (double d : proteinList.get(i).getPeptideMedianLogRatioArr()) {
                    ps.print(d);
                    ps.print("\t");

                }

                //for(int rr=0;rr<proteinList.get(i).getBestCorrelationEachpeptideIntensityList().size();rr++){
                for(double each:proteinList.get(i).getIntensityEachSampleSumList()) {
                  ps.print(each + "\t");
                   // ps.print(proteinList.get(i).getBestCorrelationEachpeptideIntensityList().get(rr).longValue()+"=====\t");
                }

                //peptide intensity sum for protein of each experiment
                /*
                for(int rr=0;rr<proteinList.get(i).getBestCorrelationEachpeptideIntensityList().size();rr++) {
                    List<LabelfreePeptide> lfPepList = proteinList.get(i).getPeptideList();
                    long totalProIntensity = 0;
                    for(Iterator<LabelfreePeptide> lfPepItr=lfPepList.iterator(); lfPepItr.hasNext(); ) {
                        LabelfreePeptide lpep = lfPepItr.next();
                        ChroPeptide cpep = lpep.getPeptideList().get(rr);
                        totalProIntensity += (long)cpep.getAverageIntensity();


                    }
                    ps.print(totalProIntensity +"\t");
                }

*/
                for(double intensity:proteinList.get(i).getIntensitySumArr()) {
                    ps.print((long)intensity +"\t");

                }

                for(int rr=0;rr<proteinList.get(i).getBestCorrelationEachpeptideIntensityListIIT().size();rr++){
                    long lll = proteinList.get(i).getBestCorrelationEachpeptideIntensityListIIT().get(rr).longValue();

                    ps.print(proteinList.get(i).getBestCorrelationEachpeptideIntensityListIIT().get(rr).longValue()+"\t");

                }
                for(int rr=0;rr<proteinList.get(i).getBestCorrelationIntensityListIIT().size();rr++){
                    long lll = proteinList.get(i).getBestCorrelationIntensityListIIT().get(rr).longValue();
                    ps.print(proteinList.get(i).getBestCorrelationIntensityListIIT().get(rr).longValue()+"\t");

                }
                /*
                for(int rr=0;rr<proteinList.get(i).getBestCorrelationIntensityList().size();rr++){
                    ps.print(proteinList.get(i).getBestCorrelationIntensityList().get(rr).longValue()+"\t");
                }*/
                ps.println();

                k = 0;
                l = 0;
                continue;
            }


            String[] words = null;
            words = eachLine.split("\t");

           // System.out.println("=================" + eachLine);


            l = 0;
            for (int j = 0; j < normIndexS; j++) {
                if (words[j].equalsIgnoreCase("NA")) {
                    if (txtReader.getSequenceIndexList().contains(j)) {
                        ps.print(proteinList.get(i).getPeptideList().get(k).getSequence() + "\t");
                    } else if (txtReader.getCsIndexList().contains(j)) {
                        ps.print(proteinList.get(i).getPeptideList().get(k).getChargeState() + "\t");

                    } else if (txtReader.getIntensityIndexList().contains(j)) {
                        ps.print( (double) proteinList.get(i).getPeptideList().get(k).getPeptideList().get(l).getAverageIntensity() + "\t");

                      //  System.out.println("====+++++========\t" +  (double)proteinList.get(i).getPeptideList().get(k).getPeptideList().get(l).getAverageIntensity());
                     //   System.out.println("====+++++========\t" +  proteinList.get(i).getPeptideList().get(k).getPeptideList().get(l).getAverageIntensity());
                        l++;
                    } else {
                        ps.print(words[j] + "\t");
                     //   System.out.println("====&&&&========\t" +  words[j]);

                    }
                } else {
                    if (txtReader.getIntensityIndexList().contains(j)) {
                        l++;
                    }
                    ps.print(words[j] + "\t");

                  //  System.out.println("====****========\t" +  words[j]);
                }

            }
            for(int rr=0;rr<proteinList.get(i).getPeptideList().get(k).getNormIntensityList().size();rr++){
                ps.print(proteinList.get(i).getPeptideList().get(k).getNormIntensityList().get(rr)+"\t");
            }
            //norm intensity
            ps.print(pvalueList.get(r));

            r++;
            ps.print("\n");
            k++;
        }

        br.close();
        ps.close();
        for(Map.Entry<String,IndexedFile> entry: origMs1FileHt.entrySet()){
            IndexedFile indexedFile = entry.getValue();
            SpectraDB db = indexedFile.getSpectraDB(false);
            if(db!=null)
            {
                db.close();;
            }
        }
     //   return proteinList;
    }

    private static List<LabelFreeJSONPeptide> getEachSamplePeptides(List<ChroPeptide> pepList, int count, int pepCount,
                                                                    int chargeState, String sequence, double startRt, double endRt) {

        List<LabelFreeJSONPeptide> eachSamplePeptideList = new ArrayList<>();

        LabelFreeJSONPeptide jsonPeptide = null;

        for (ChroPeptide chroPeptide : pepList) {
            jsonPeptide = new LabelFreeJSONPeptide();
            jsonPeptide.setCount(String.valueOf(count));
            jsonPeptide.setUnique(String.valueOf(chroPeptide.isUnique()));
            jsonPeptide.setChro_iso(formatChroData(chroPeptide.getChroData() == null ? "" : chroPeptide.getChroData()));
            jsonPeptide.setSeq(sequence);
            jsonPeptide.setScan(String.valueOf(chroPeptide.getScanNum()));
            jsonPeptide.setFile(chroPeptide.getFileName() == null ? "" : chroPeptide.getFileName());
            jsonPeptide.setStart_scan(chroPeptide.getStartRange());
            jsonPeptide.setEnd_scan(chroPeptide.getEndRange());
            jsonPeptide.setCharge(String.valueOf(chargeState));
            jsonPeptide.setRowId(String.valueOf(pepCount));
            jsonPeptide.setStartRt(String.valueOf(startRt));
            jsonPeptide.setEndRt(String.valueOf(endRt));
            jsonPeptide.setPeaks(chroPeptide.getGaussianPeakString());
            jsonPeptide.setMaxPeakValue(getHighestPaeakValue(chroPeptide.getGaussianPeakString()));
            eachSamplePeptideList.add(jsonPeptide);
            pepCount++;

            // System.out.println(pepCount);
        }

        return eachSamplePeptideList;
    }

    private static String formatChroData(String chroData) {

        if (chroData != null && !chroData.equals("")) {
            int firstOcurrentCol = chroData.indexOf(";");
            String prefixStr = chroData.substring(0, firstOcurrentCol);
            String prefStrArr[] = prefixStr.split(" ");
//    		chroData = prefStrArr[prefStrArr.length-1] + chroData.substring(firstOcurrentCol,chroData.length() );
            chroData = chroData.substring(firstOcurrentCol + 1, chroData.length());
        }

        return chroData;
    }


    private static String[] getHighestPaeakValue(String paeakData){

        String[] returnValue = new String[2];


        try {

            if(paeakData != null){

                String[] valueArr = paeakData.split(";");

                List<Double> valueList = new ArrayList<>();
                Map<String, Double> peakMap = new HashMap<>();

                for (String xy : valueArr) {
                    String[] xyValue = xy.split(" ");
                    Double intensity = Double.valueOf(xyValue[1]);
                    peakMap.put(xyValue[0], intensity);
                    valueList.add(intensity);
                }

                Double maxVal = Collections.max(valueList);

                for (Map.Entry<String, Double> entry : peakMap.entrySet()) {
                    String key = entry.getKey();
                    Double mapValue = entry.getValue();
                    if(mapValue == maxVal){
                        //	System.out.println(key + " " + mapValue + " " + maxVal);
                        returnValue[0] = key;
                        returnValue[1] = String.valueOf(mapValue);
                    }
                }
            }

        } catch (Exception e) {

//don't generate error message when peak is not found
//			System.out.println("Error in finding highest Peak Value :: " + e.getMessage());
        }

        //System.out.println(" Peak Value :: "+returnValue);

        return returnValue;

    }


    //public static GaussianPeakModel isotopeCalc(int startScan, int endScan, int startIndex, int endIndex, IsotopeReader isoReader, String sequence, int chargeState, IndexedFile iFile) throws Exception {
    public static GaussianPeakModel isotopeCalc(int startScan, int endScan, int startIndex, int endIndex,
                                                IsotopeReader isoReader, String sequence,
                                                int chargeState,
                                                IndexedFile origIFile,
                                                Map<String, String> splitSpectraMap,
                                                Map<String, IndexedFile> splitMs1FileHt

    ) throws Exception {

        return isotopeCalc(startScan, endScan, startIndex, endIndex, isoReader, sequence,
            chargeState, origIFile, splitSpectraMap, splitMs1FileHt, false);

    }
    //public static GaussianPeakModel isotopeCalc(int startScan, int endScan, int startIndex, int endIndex, IsotopeReader isoReader, String sequence, int chargeState, IndexedFile iFile) throws Exception {
    public static GaussianPeakModel isotopeCalc(int startScan, int endScan, int startIndex, int endIndex,
                                                IsotopeReader isoReader, String sequence,
                                                int chargeState,
                                                IndexedFile origIFile,
                                                Map<String, String> splitSpectraMap,
                                                Map<String, IndexedFile> splitMs1FileHt,
                                                boolean fillZeroOnly

    ) throws Exception {

        String origFile = origIFile.getFileName();

        if (sequence.contains(".")) {
            sequence = sequence.substring(2, sequence.length() - 2);
        }

        /**
         * *********************************
         * 1 calculate isotope distribution *********************************
         */
        char[] ch = sequence.toCharArray();

        ElementComposition element = new ElementComposition(ch, 0, ch.length, isoReader.getIsotope());
        element.calculate();

        if (!element.isQuantifiable()) {
            System.out.print("\nError : ");
            System.out.println(sequence + " is not quantifiable.");
            return null;
        }

        Configuration conf = Configuration.getInstance();

        IsotopeDist sampleDist = new IsotopeDist(element.getElementSampleArr(), element.getModShift(), true);

//                        TIntDoubleHashMap retentionTimeMap = iFile.getRetentionTimeMap();
        //                       TIntDoubleHashMap ionInjectionMap = iFile.getIonInjectionMap();
        //  ionInjectionTime = ionInjectionMap.get(tmpScanNumber);

        double[] isoArr = sampleDist.getHighMassList();
        double[] isoIntArr = sampleDist.getRelabun(isoArr.length);
        double pepMass = sampleDist.getHighMassList()[0];

        for (int i = 0; i < isoArr.length; i++) {
            isoArr[i] = (isoArr[i] + chargeState * CensusConstants.PROTON_MASS) / chargeState;
        }

        /**
         * *********************************
         * 2 Re-construct chromatogram *********************************
         */
        //long currentPos = iFile.getPositionByIndex(startIndex);
        //long nextPos = -1;
        long[] chromPeakArr = new long[endIndex - startIndex + 1];
        double[] retArr = new double[chromPeakArr.length];
        int[] scanArr = new int[chromPeakArr.length];
        int count = 0;

        boolean hasPeak=false;
        if(!fillZeroOnly) {
            for (int i = startIndex; i <= endIndex; i++) {

                int eachScan = origIFile.getKeys()[i];
                String fileKey = eachScan + "\t" + origFile;
                String spltiMs1File = splitSpectraMap.get(fileKey);
                IndexedFile splitIFile = splitMs1FileHt.get(spltiMs1File);
                int splitCurIndex = splitIFile.getIndexByScan(eachScan);


                //    System.out.println(currentPos + " " + nextPos + " " + (nextPos-currentPos));
                //diff is for reading byte range

                long nextPos = 0;
                long currentPos = splitIFile.getPositionByIndex(splitCurIndex);


                if ((splitCurIndex + 1) >= splitIFile.getScanPositionMap().size())
                    nextPos = splitIFile.getFile().length();
                else
                    nextPos = splitIFile.getPositionByIndex(splitCurIndex + 1);

                int diff = (int) (nextPos - currentPos);

                SpectrumModel spec = CalcUtilGeneric.labelFreeSpectrumReader(isoArr, currentPos, diff,
                        conf.getMassTolerance(), splitIFile, chargeState, conf, pepMass);/**/
                chromPeakArr[count] = spec.getPrecursorPeakIntensity();
                if (!hasPeak && chromPeakArr[count] > 0)
                    hasPeak = true;

                retArr[count] = spec.getRetentionTime();
                scanArr[count] = spec.getScanNumber();

                count++;
                //  System.out.println(i + "-\t" + resultArr[0]);
                currentPos = nextPos;

            }
        }

        /**
         * *********************************
         * 3 Smooth chromatogram *********************************
         */
        if(!hasPeak) {
            GaussianPeakModel gModel = new GaussianPeakModel();
            gModel.setScanArr(scanArr);
            gModel.setRetArr(retArr);
            gModel.setHasPeak(false);

            int currentScan = origIFile.getKeys()[startIndex];
            String fileKey = currentScan + "\t" + origFile;
            String spltiMs1File = splitSpectraMap.get(fileKey);
            IndexedFile splitIFile = splitMs1FileHt.get(spltiMs1File);
            int splitCurIndex = splitIFile.getIndexByScan(currentScan);

            long startPos = splitIFile.getPositionByScan(currentScan);

            long backGroundNoise = 0;
            long currentPos = splitIFile.getPositionByIndex(splitCurIndex);
            long nextPos = 0;

            if ((splitCurIndex + 1) >= splitIFile.getScanPositionMap().size())
                nextPos = splitIFile.getFile().length();
            else
                nextPos = splitIFile.getPositionByIndex(splitCurIndex + 1);

            int diff = (int) (nextPos - currentPos);
            backGroundNoise = CalcUtilGeneric.getBackGroundNoise(currentPos, diff, splitIFile);


            gModel.setPeakArea(backGroundNoise);
          //  gModel.setPeakArea(0);

            double[] noPeakArr = new double[chromPeakArr.length];
            gModel.setPeakArr(noPeakArr);

            //   System.out.println("peak===\t" + range.getPeakArea());
            //     System.out.println("done..");
            //chromPeakArr, retArr, peakStartIndex, peakEndIndex
            return gModel;
        }

        double[] smoothChromArr = Smooth.smoothAsDouble(chromPeakArr, LabelfreeMissingPeptideBuilderSplit.SMOOTH_WINDOW_SIZE);

        double basePeak = 0;
        int basePeakIndex = 0;

        for (int i = 0; i < smoothChromArr.length; i++) {

            //     System.out.println("==\t" + retArr[i] + "\t" + smoothChromArr[i] + "\t" + chromPeakArr[i]);
            if (basePeak < smoothChromArr[i]) {
                basePeak = smoothChromArr[i];
                basePeakIndex = i;
            }

            // System.out.println(smoothChromArr[i]);
        }

        /**
         * *********************************
         * 4. Find simple/rough peak range (1/3 of base peak) for Gaussian input
         * *********************************
         */


        int[] indexResult = LabelfreeChroUtil.getPeakRange(basePeakIndex, basePeak, smoothChromArr);
        int peakStartIndex = indexResult[0];
        int peakEndIndex = indexResult[1];

        // GaussianPeakModel range = test(retArr, smoothChromArr, peakStartIndex, peakEndIndex);
        GaussianPeakModel gModel = GaussianFitting.getGaussianPeakRangeIndex(retArr, smoothChromArr, peakStartIndex, peakEndIndex);
        gModel.setScanArr(scanArr);
        gModel.setRetArr(retArr);
        gModel.setPeakArr(smoothChromArr);


        double maxPeakIntensity = 0;
        for(double d:smoothChromArr) {
            if(d>maxPeakIntensity)
                maxPeakIntensity = d;
        }

        gModel.setMaxIntensity(maxPeakIntensity);


        //   System.out.println("peak===\t" + gModel.getPeakArea());
       //      System.out.println("done..");
        //chromPeakArr, retArr, peakStartIndex, peakEndIndex
        return gModel;

        //    chroText =
        /*
         //chro.setText( chroText );
         String[] tmpStrArr = chroText.substring(0, chroText.indexOf(";")).split(" ");
         peptideEle.setAttribute("start_scan", tmpStrArr[1]);
         peptideEle.setAttribute("end_scan", tmpStrArr[2]);


         peptideEle.addContent(chro);

         return peptideEle;





         */
    }

    public static GaussianPeakModel isotopeCalc(int startIndex, int endIndex,
                                                IsotopeReader isoReader, String sequence,
                                                int chargeState,
                                                IndexedFile origIFile,
                                                boolean fillZeroOnly

    ) throws Exception {


        if (sequence.contains(".")) {
            sequence = sequence.substring(2, sequence.length() - 2);
        }

        /**
         * *********************************
         * 1 calculate isotope distribution *********************************
         */
        char[] ch = sequence.toCharArray();

        ElementComposition element = new ElementComposition(ch, 0, ch.length, isoReader.getIsotope());
        element.calculate();

        if (!element.isQuantifiable()) {
            System.out.print("\nError : ");
            System.out.println(sequence + " is not quantifiable.");
            return null;
        }

        Configuration conf = Configuration.getInstance();

        IsotopeDist sampleDist = new IsotopeDist(element.getElementSampleArr(), element.getModShift(), true);

        double[] isoArr = sampleDist.getHighMassList();
        double[] isoIntArr = sampleDist.getRelabun(isoArr.length);
        double pepMass = sampleDist.getHighMassList()[0];

        for (int i = 0; i < isoArr.length; i++) {
            isoArr[i] = (isoArr[i] + chargeState * CensusConstants.PROTON_MASS) / chargeState;
        }

        /**
         * *********************************
         * 2 Re-construct chromatogram *********************************
         */

        long[] chromPeakArr = new long[endIndex - startIndex + 1];
        double[] retArr = new double[chromPeakArr.length];
        int[] scanArr = new int[chromPeakArr.length];
        int count = 0;
        SpectraDB spectraDB = origIFile.getSpectraDB();

        boolean hasPeak=false;
        if(!fillZeroOnly) {
            for (int i = startIndex; i <= endIndex; i++) {

                int eachScan = origIFile.getKeys()[i];


                SpectrumModel spec = CalcUtilGeneric.labelFreeSpectrumReader(isoArr, eachScan,
                        conf.getMassTolerance(), spectraDB, chargeState, conf, pepMass);/**/
                chromPeakArr[count] = spec.getPrecursorPeakIntensity();
                if (!hasPeak && chromPeakArr[count] > 0)
                    hasPeak = true;

                retArr[count] = spec.getRetentionTime();
                scanArr[count] = spec.getScanNumber();

                count++;

            }
        }

        /**
         * *********************************
         * 3 Smooth chromatogram *********************************
         */
        if(!hasPeak) {
            GaussianPeakModel gModel = new GaussianPeakModel();
            gModel.setScanArr(scanArr);
            gModel.setRetArr(retArr);
            gModel.setHasPeak(false);

            int currentScan = origIFile.getKeys()[startIndex];


            double backGroundNoise = CalcUtilGeneric.getBackGroundNoise(currentScan, spectraDB);


            gModel.setPeakArea(backGroundNoise);

            double[] noPeakArr = new double[chromPeakArr.length];
            gModel.setPeakArr(noPeakArr);

            return gModel;
        }

        double[] smoothChromArr = Smooth.smoothAsDouble(chromPeakArr, LabelfreeMissingPeptideBuilderSplit.SMOOTH_WINDOW_SIZE);

        double basePeak = 0;
        int basePeakIndex = 0;

        for (int i = 0; i < smoothChromArr.length; i++) {

            if (basePeak < smoothChromArr[i]) {
                basePeak = smoothChromArr[i];
                basePeakIndex = i;
            }

        }

        /**
         * *********************************
         * 4. Find simple/rough peak range (1/3 of base peak) for Gaussian input
         * *********************************
         */


        int[] indexResult = LabelfreeChroUtil.getPeakRange(basePeakIndex, basePeak, smoothChromArr);
        int peakStartIndex = indexResult[0];
        int peakEndIndex = indexResult[1];

        GaussianPeakModel gModel = GaussianFitting.getGaussianPeakRangeIndex(retArr, smoothChromArr, peakStartIndex, peakEndIndex);
        gModel.setScanArr(scanArr);
        gModel.setRetArr(retArr);
        gModel.setPeakArr(smoothChromArr);


        double maxPeakIntensity = 0;
        for(double d:smoothChromArr) {
            if(d>maxPeakIntensity)
                maxPeakIntensity = d;
        }

        gModel.setMaxIntensity(maxPeakIntensity);

        return gModel;
    }


    //return parameters: y, x, and sigma
    //peak ranges from -3 x sigma to 3 x sigma
    public static GaussianPeakModel test(double[] xArr, double[] yArr, int startIndex, int endIndex) {

        /*
         Hashtable<Double, Integer> scanIndexHt = new Hashtable<Double, Integer>();
         for(int i=startIndex;i<=endIndex;i++) {


         scanIndexHt.put((double)resultArr[i].getScanNum(), i);
         }*/
        double[] params = getGaussianCurveFitRange(xArr, yArr, startIndex, endIndex);
        /////////////////HOLD ON THIS................
        if (params == null) {
            return new GaussianPeakModel(startIndex, endIndex);
        }

        double start = -3 * params[2] + params[1];
        double end = 3 * params[2] + params[1];
//        System.out.println("start...........");

        int i = 0;

        while (i < xArr.length && xArr[i] < start) {
            i++;
        }
        if (i > 1) {
            start = xArr[i - 1];
        } else {
            start = xArr[0];
        }

        while (i < xArr.length && xArr[i] < end) {
            i++;
        }
        if (i < xArr.length - 1) {
            end = xArr[i + 1];
        } else {
            end = xArr[xArr.length - 1];
        }

        GaussianPeakModel peakModel = new GaussianPeakModel(start, end);
        peakModel.setX(params[1]);
        peakModel.setY(params[0]);
        peakModel.setSigma(params[2]);

        return peakModel;

    }

    public static void proteinCompareBasedOnRatios(List<ProteinModel> proteinList, List<Integer> sample1List, List<Integer> sample2List) {
        Iterator i = proteinList.iterator();
        int counter = 0;
        sample2List.add(4);
        List<Integer> sampleList2new = sample2List;
        List<Integer> sampleList1new = sample1List;
        if(sample1List.size() < sample2List.size()){
            sampleList2new=new ArrayList<>();
            for(int h=0;h<sample1List.size();h++){
                sampleList2new.add(sample2List.get(h));
            }
        }else{
            sampleList1new=new ArrayList<>();
            for(int h=0;h<sample2List.size();h++){
                sampleList1new.add(sample1List.get(h));
            }
        }
        int reportIonCount = sampleList1new.size() + sampleList2new.size();
        while (i.hasNext()) {
            ProteinModel protein = (ProteinModel) i.next();

            List<LabelfreePeptide> peptideList = protein.getPeptideList();
            List<HashMap<String, Double>> proPeptideInfo = new ArrayList();
            for (int j = 0; j < peptideList.size(); j++) {
                List<ChroPeptide> expList = peptideList.get(j).getPeptideList();
                List<Double> avgIntensityList = new ArrayList<>();
                List permCount = new ArrayList();
                HashMap<String, Double> myHash = new HashMap<>();
                for (int k = 0; k < expList.size(); k++) {
                    ChroPeptide peptide = expList.get(k);
                    //System.out.println("********************");
                    // System.out.println(""+peptide.getAverageIntensity());
                    avgIntensityList.add(peptide.getAverageIntensity());
                }

                ICombinatoricsVector<Integer> initialVector = Factory.createVector(sampleList2new);
                Generator<Integer> generator = Factory.createPermutationGenerator(initialVector);
                for (ICombinatoricsVector<Integer> perm : generator) {
                    permCount.add(perm);
                }
                counter = permCount.size();
                permCount.clear();
                Iterator<ICombinatoricsVector<Integer>> itr = generator.iterator();
                int r = 0;
                int o = 1;
                while (itr.hasNext()) {
                    ICombinatoricsVector<Integer> permutation = itr.next();
                    List<Integer> myvector = permutation.getVector();
                    Iterator vItr = myvector.iterator();
                    while (vItr.hasNext()) {
                        Object obj = vItr.next();
                        double ratio;
                        ratio = avgIntensityList.get((int) sampleList1new.get(r)) / avgIntensityList.get((int) obj);
                        myHash.put("" + o, ratio);
                        r++;
                        o++;
                    }
                    r = 0;
                }
                proPeptideInfo.add(myHash);
                avgIntensityList.clear();
            }

            List expValueList = new ArrayList();
            List pepMedianList = new ArrayList();
            for (int l = 1; l <= ((reportIonCount / 2) * counter); l++) {

                Iterator h = proPeptideInfo.iterator();
                while (h.hasNext()) {

                    HashMap obj = (HashMap) h.next();
                    expValueList.add(obj.get("" + l));

                }
                DescriptiveStatistics stats = new DescriptiveStatistics();

                // Add the data from the array
                for (int z = 0; z < expValueList.size(); z++) {
                    stats.addValue((Double) expValueList.get(z));
                }

                // Compute some statistics
                double median = stats.getPercentile(50);
                pepMedianList.add(median);
                expValueList.clear();

            }
            List<Double> pValues = new ArrayList();
            double[] arrDoubleRatios = new double[reportIonCount / 2];
            int k = 0;
            for (int b = 0; b < pepMedianList.size(); b += reportIonCount / 2) {

                for (int f = 0; f <= ((reportIonCount / 2) - 1); f++) {

                    arrDoubleRatios[f] = (double) pepMedianList.get(k);
                    k++;
                }

                double tmp2;
                if (arrDoubleRatios.length > 1) {
                    tmp2 = TestUtils.tTest(1.0, arrDoubleRatios);
                    pValues.add(tmp2);
                } else {
                    pValues.add(1d);
                }
            }

            DescriptiveStatistics stats2 = new DescriptiveStatistics();

            // Add the data from the array
            for (Object pValue : pValues) {
                stats2.addValue((Double) pValue);
            }

            // Compute some statistics
            double median2 = stats2.getPercentile(50);
            protein.setMultipleRatiopValue(median2);
        }

        System.out.println("Done");

    }        boolean  first = true;



    public static String createProteinString(List<Double> intensityList)
    {
        StringBuilder sb = new StringBuilder();
        for(double d: intensityList)
        {
            sb.append(Double.toString(d));
        }
        return sb.toString();
    }



    public static List<double[]> createLogRatios(List intensities)
    {
        Iterator itr = intensities.iterator();

        double [] farr = (double []) itr.next();
        int controlSize = farr.length; //1st array
        List<double[]> ratioList = new ArrayList<>();
        while(itr.hasNext())
        {

            double [] arr = (double [] ) itr.next();
            int eachSize = arr.length;

            double [] ratioArr = new double[arr.length];
            if(farr.length == arr.length) {
              for (int i = 0; i < arr.length; i++) {
                ratioArr[i] = Math.log(arr[i] / farr[i]) / Math.log(2);
                if(Double.isInfinite(ratioArr[i])) ratioArr[i] =10;
              }
            }
            else
            {
                for (int i = 0; i < arr.length; i++) {
                    ratioArr[i] = 1.0;
                }
            }

            ratioList.add(ratioArr);
        }
        return ratioList;
    }

    public static List<Double> createPValuesOneSampleTTest(List<double []> ratios) throws  Exception
    {
        List<Double> result = new ArrayList<>();

        for(double [] arr: ratios)
        {

            double pvalue = 1.0;
            if(arr.length>1)
              pvalue = TTestUtil.oneSampleTTest(arr);

            if(Double.isNaN(pvalue) || Double.isInfinite(pvalue))
            {
                pvalue = 1.0;
            }
            result.add(pvalue);
        }
        return  result;
    }

    public static List<Double> createPValuesCompareWithFirst(List intensities) throws Exception
    {
        Iterator itr = intensities.iterator();
        double [] farr = (double []) itr.next();
        List<Double> result = new ArrayList<>();
        while(itr.hasNext())
        {
            double [] arr = (double []) itr.next();
            double pvalue = 1.0;
            if(arr.length>1 && farr.length> 1 && arr.length==farr.length)
                pvalue = TTestUtil.calculateUnpairedTTest(farr,arr);
            if(Double.isInfinite(pvalue)|| Double.isNaN(pvalue))
            {
                pvalue = 1.0;
            }
            result.add(pvalue);
        }
        return result;
    }
}
