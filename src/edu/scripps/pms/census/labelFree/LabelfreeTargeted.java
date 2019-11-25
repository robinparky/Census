/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.scripps.pms.census.labelFree;

import static edu.scripps.pms.census.labelFree.LabelFreeParser.printUsage;
import static rpark.statistics.GaussianFitting.getGaussianCurveFitRange;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.*;

import edu.scripps.pms.census.CensusConstants;
import edu.scripps.pms.census.ChroGenerator;
import edu.scripps.pms.census.ElementComposition;
import edu.scripps.pms.census.conf.Configuration;
import edu.scripps.pms.census.hash.IndexedFile;
import edu.scripps.pms.census.io.ChroXmlReader;
import edu.scripps.pms.census.io.IsotopeReader;
import edu.scripps.pms.census.labelFree.json.ChroJSONCreator;
import edu.scripps.pms.census.labelFree.json.LabelFreeJSONPeptide;
import edu.scripps.pms.census.labelFree.json.LabelFreeJSONProtein;
import edu.scripps.pms.census.labelFree.model.LabelfreePeptide;
import edu.scripps.pms.census.labelFree.util.LabelfreeChroUtil;
import edu.scripps.pms.census.model.ChroPeptide;
import edu.scripps.pms.census.model.ChroProtein;
import edu.scripps.pms.census.model.SampleGroup;
import edu.scripps.pms.census.model.SampleModel;
import edu.scripps.pms.census.util.CalcUtilGeneric;
import edu.scripps.pms.census.util.IsotopeDist;
import edu.scripps.pms.census.util.io.FileUtil;
import edu.scripps.pms.util.PhosphoUtil;
import gnu.trove.TDoubleIntHashMap;
import java.io.Reader;

import org.apache.commons.math.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.inference.TestUtils;
import org.apache.commons.math3.stat.regression.RegressionResults;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.jdom.Element;
import org.jfree.util.ArrayUtils;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.paukov.combinatorics.Factory;
import org.paukov.combinatorics.Generator;
import org.paukov.combinatorics.ICombinatoricsVector;
import rpark.statistics.AnovaUtil;
import rpark.statistics.BHCorrection;
import rpark.statistics.BinarySearch;
import rpark.statistics.GaussianFitting;
import rpark.statistics.Smooth;
import rpark.statistics.model.GaussianPeakModel;
import scripts.MSSplitFolderCreation;
import edu.scripps.pms.census.hash.IndexUtil;

/**
 *
 * @author rpark
 * @author rohan
 */
public class LabelfreeTargeted {

    public static int SMOOTH_WINDOW_SIZE = 5;

    public static void main(String[] args) throws Exception {

     //   if (args.length < 3) {
       //     printUsage();
       //     return;
       // }
    	
    	    System.out.println("LabelfreeTargeted - new version");

        //String configFile = args[0];
        String configFile = args[0];
        String peptideFile = args[1];
        String folderName = args[2];
        //args -c /home/rpark/test_data/targeted_lfree/census_config_labelfree_target.xml
        /*
        String configFile = "/home/rpark/test_data/targeted_lfree/census_config_labelfree_target.xml";
  //      String configFile = "/home/rpark/test_data/targeted_lfree/census_config_labelfree_target_small.xml";
        String peptideFile = "/home/rpark/test_data/targeted_lfree/peptides.txt";
        String jsonFile = "/home/rpark/test_data/targeted_lfree/targets.JSON";
        String folderName = "/home/rpark/test_data/targeted_lfree";
        */

        Configuration conf = Configuration.getInstance();
        // conf.setLabelfreeCheckChargeState(true);

        if (!conf.isReadConfigFile()) {
            conf.setLabelfree(true);
            conf.readXMLParam(configFile);
        }
        List<ProteinModel> proteinList = LabelfreeTargeted.targetedAnalysis(folderName, peptideFile,
                conf);
        /*List<ProteinModel> proteinList = LabelfreeTargeted.targetedAnalysisSumPeakArea(folderName, peptideFile,
                conf);*/

/*
        String configFile = "/home/rampuria/projects/deep/test_20December/census_config_labelfree_target.xml";
        String peptideFile = "/home/rampuria/projects/deep/test_20December/peptides.txt";
        String jsonFile = "/home/rampuria/projects/deep/test_20December/out.JSON";
        String folderName = "/home/rampuria/projects/deep/test_20December";
*/



    }

    public static void generateLabelfreeOutputFile(List<ProteinModel> proteinList, String filename, String confFile) throws Exception {
        Configuration conf = Configuration.getInstance();
        if (!conf.isReadConfigFile()) {
            conf.setLabelfree(true);
            conf.readXMLParam(confFile);
        }
        //HashMap<String,List<Double>> intensityMap = new HashMap<>();
        HashMap<String,List<Long>> intensityMap = new HashMap<>();
        HashMap<String,List<Long>> intensityCorrMap = new HashMap<>();

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

                intensityMap.put(p.getRedundnatProteinList().get(a).getLocus(), p.getBestCorrelationEachpeptideIntensityList());
                intensityCorrMap.put(p.getRedundnatProteinList().get(a).getLocus(), p.getBestCorrelationEachpeptideIntensityListIIT());
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
        /*
        calculating anova p-value
        */
        List<Double> pvalues = new ArrayList<>();
        List<Double> pvalues2 = new ArrayList<>();
        HashMap<String,Double> pvalueMap = new HashMap<>();
        HashMap<String,Double> pvalueMap2 = new HashMap<>();
        for(String key:intensityMap.keySet()){
            List classes = new ArrayList();
            List classes2 = new ArrayList();

            int counter=0;
            for(int i=0;i<indexList.size();i++){

                double [] intensity =new double[indexList.get(i)];
                double [] corrIntensity =new double[indexList.get(i)];
                for(int j=0;j<indexList.get(i);j++){

                    intensity[j] = intensityMap.get(key).get(counter);
                    corrIntensity[j] = intensityCorrMap.get(key).get(counter);
                    counter++;
                }
                //counter =counter+indexList.get(i);
                classes.add(intensity);
                classes2.add(corrIntensity);

            }

            double pvalue = AnovaUtil.calculateAnovaPvalue(classes);
            if(Double.isInfinite(pvalue) || Double.isNaN(pvalue))
                pvalue = 1;

            double pvalue2 = AnovaUtil.calculateAnovaPvalue(classes2);
            if(Double.isInfinite(pvalue2) || Double.isNaN(pvalue2))
                pvalue2 = 1;

            pvalues.add(pvalue);
            pvalues2.add(pvalue2);
            pvalueMap.put(key, pvalue);
            pvalueMap2.put(key, pvalue2);

        }

        /*
        BH correction q-value
        */
        List<Double> bhcorr = BHCorrection.runBhCorrection(pvalues);
        List<Double> bhcorr2 = BHCorrection.runBhCorrection(pvalues2);
       // List<Double> bhcorr3 = BHCorrection.runBhCorrection(ratiobasedpval);
        HashMap<String,Double> bhcorrMap = new HashMap<>();
        HashMap<String,Double> bhcorrMap2 = new HashMap<>();
        int l=0;
        for(String key : intensityMap.keySet()){
            bhcorrMap.put(key,bhcorr.get(l));
            bhcorrMap2.put(key, bhcorr2.get(l));
            l++;
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

        p.print("PLINE\tACCESSION\tINTENSITY_P-VALUE\tIIT_INTENSITY_P-VALUE\tP-VALUE_RATIO_BASED\tINTENSITY_Q-VALUE"
                + "\tIIT_INTENSITY_Q-VALUE\t");
        for (int i = 0; i < samGroupEleList.get(0).getChildren().size(); i++) {

            p.print("MEDIAN_LOG_RATIO_" + (i + 1));
            p.print("\t");
        }
        for (int i = 0; i < expSize; i++) {
            p.print("SCOUNT_" + (i + 1));
            p.print("\t");
        }
        for (int i = 0; i < expSize; i++) {
            p.print("INTENSITY_" + (i + 1));
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

        //p.print("PEP_COUNT\t");

        //for(int i=0;i<expSize;i++) {
        //  p.print("NORM_INTENSITY_" + (i+1));
        //  p.print("\t");
//	}
        p.println("DESCRIPTION\t");

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
                p.print(pvalueMap2.get(cprotein.getLocus()));
                p.print("\t");
               // p.print(protein.getMultipleRatiopValue());
                p.print("1.0\t");
                p.print(bhcorrMap.get(cprotein.getLocus()));
                p.print("\t");
                p.print(bhcorrMap2.get(cprotein.getLocus()));
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

                for(Long each:protein.getBestCorrelationEachpeptideIntensityList()) {
                    p.print(each); p.print("\t");
                }

                for (int i = 0; i < indexList.size(); i++) {
                    List<Long> list = protein.getBestCorrelationIntensityList();
                    p.print(list.get(i) + "\t");

                }
                for (int i = 1; i < indexList.size(); i++) {

                    List<Long> list = protein.getBestCorrelationIntensityList();
                    double ratio = (double)list.get(i)/list.get(0);
                    p.print(ratio);
                    p.print("\t");

                }
                for (int i = 1; i < indexList.size(); i++) {

                    List<Long> list = protein.getBestCorrelationIntensityList();
                    double ratio = (double)list.get(i)/list.get(0);
                    String logRatio = "NA";
                    try {
                        double lr = Math.log(ratio)/Math.log(2);
                        p.print(lr);
                    } catch (Exception e) {
                        p.print(logRatio);
                    }

                    p.print("\t");

                }

                for(long each:protein.getBestCorrelationEachpeptideIntensityListIIT()) {
                    p.print(each);
                    p.print("\t");

                }

                for(long each:protein.getBestCorrelationIntensityListIIT()) {
                    p.print(each);
                    p.print("\t");
                }

                for (int i = 1; i < protein.getBestCorrelationIntensityListIIT().size(); i++) {
                    List<Long> list = protein.getBestCorrelationIntensityListIIT();

                    double ratio = (double)list.get(i)/list.get(0);
                    p.print(ratio);
                    p.print("\t");

                }
                for (int i = 1; i < indexList.size(); i++) {
                    List<Long> list = protein.getBestCorrelationIntensityListIIT();

                    String logRatio = "NA";
                    try {

                        double ratio = list.get(i)/list.get(0);
                        double lr = Math.log(ratio)/Math.log(2);
                        p.print(lr);
                    } catch (Exception e) {
                        p.print(logRatio);
                    }

                    p.print("\t");
                }
                for(int i=0;i<protein.getNormIntensityList().size();i++){
                    p.print(protein.getNormIntensityList().get(i).longValue()+"\t");
                }

                p.println(cprotein.getDescription());

            }

            //List<String> keys = p.getPeptideKey();
            //  HashMap<String, List<ChroPeptide>> peptideMap = p.getPeptideMap();
        }

        p.close();

        System.out.println("check file " + filename);

    }
/*
    public static List<ProteinModel> targetedAnalysisSumPeakArea(
            String folder, String peptideFile,
            Configuration conf) throws Exception {





        List<org.jdom.Element> samGroupEleList = conf.getRootConfEle().getChildren("sample");
        int expSize = 0;

        List<Integer> indexList = new ArrayList<>();
        for (Iterator<org.jdom.Element> samgItr = samGroupEleList.iterator(); samgItr.hasNext();) {
            int count =0;
            org.jdom.Element groupEle = samgItr.next();

            List<Element> sampleEleList = groupEle.getChildren("each_sample");

            for (Iterator<Element> samItr = sampleEleList.iterator(); samItr.hasNext();) {
                org.jdom.Element eachSample = samItr.next();
                count++;
                expSize++;
            }
            indexList.add(count);
        }


        IsotopeReader isoReader = new IsotopeReader(conf.getRootConfEle());
        //IdentificationReader idReader = BaseIdentificationReader.getIdentificationInst(isoReader);

        //   conf.setSpectrumFormat(Configuration.MS_FILE_FORMAT);

        //   Hashtable<String, ChroPeptide> peptideChroHt = new Hashtable<>();

        List<SampleModel> allSampleList = conf.getSampleList();
        // sampleList.get(0).getPathList();
        Hashtable<String, IndexedFile> origMs1FileHt = new Hashtable<>();
        ChroJSONCreator chroJSONCreator = new ChroJSONCreator();
        String jsonFilePath = new File(folder) + File.separator + "JSON_OBJ";
        FileUtil.makeDir(jsonFilePath);
        List<LabelFreeJSONProtein> jsonProteins = new ArrayList<>();
        LabelFreeJSONProtein jsonProtein = null;
        StringBuffer sb = new StringBuffer();
        StringBuffer gPeakSb = new StringBuffer();

        jsonProtein = new LabelFreeJSONProtein();
        jsonProtein.setAccession("Targeted peptides");
        jsonProtein.setDesc("Targeted peptides");
        jsonProteins.add(jsonProtein);


        MSSplitFolderCreation msp = new MSSplitFolderCreation();
        Map<String, String> splitSpectraMap = new HashMap<>();
        Map<String, IndexedFile> splitMs1FileHt = new HashMap<>();
        Map<String, HashMap<Integer, Integer>> ms2ToMs1Map = new HashMap<>();

        List<SampleGroup> sampleGroupList = conf.getSampleGroupList();


        BufferedReader br = new BufferedReader(new FileReader(peptideFile));
        String eachLine;
        List<LabelfreePeptide> pepList = new ArrayList<>();
        List<IsotopeDist> isoList = new ArrayList<>();
        while (null != (eachLine = br.readLine())) {


            for (int cs = conf.getTargetedStartCharge(); cs <= conf.getTargetedEndCharge(); cs++) {
                //System.out.println(eachLine);
                LabelfreePeptide pep = new LabelfreePeptide();
                String arr[] = eachLine.split("\t");
                String sequence = arr[0];
                pep.setSequence(arr[0]);
                pep.setChargeState(cs);
                pepList.add(pep);

                char[] ch = sequence.toCharArray();


                ElementComposition element = new ElementComposition(ch, 0, ch.length, isoReader.getIsotope());
                element.calculate();

                if (!element.isQuantifiable()) {
                    System.out.print("\nError : ");
                    System.out.println(sequence + " is not quantifiable.");
                    return null;
                }
                IsotopeDist sampleDist = new IsotopeDist(element.getElementSampleArr(), element.getModShift(), true);
                isoList.add(sampleDist);

                pep.setIsotopeDist(sampleDist);
            }

        }

        br.close();


        List<List<LabelFreeJSONPeptide>> jsonAllPeptideList = new ArrayList<>();
        String jsonPeptideListFileName = jsonFilePath + File.separator + "TARGETED_PEPTIDE.JSON";

        for (SampleGroup sampleGroup : sampleGroupList) {

            for (SampleModel sampleModel : sampleGroup.getSampleModelList()) {

                for (String path : sampleModel.getPathList()) {

                    if (!path.endsWith("/"))
                        path += "/";

                    String spectraPath = path;
                    String splitSpectraPath = path + "split/";

                    System.out.println("===" + path);

                    splitSpectraMap.putAll(msp.splitMS1Files(spectraPath, CensusConstants.LABELFREE_MS1_SPLIT_SCAN_NUM, false));
                    splitMs1FileHt.putAll(ChroGenerator.createIndexedFiles(splitSpectraPath, CensusConstants.MS1_FILE, true, true));
                    ms2ToMs1Map.putAll(IndexUtil.buildMS2toMS1ScanMapFiles(spectraPath));

                    Hashtable<String, IndexedFile> ht = ChroGenerator.createIndexedFiles(spectraPath, "ms1", true, true);

                    origMs1FileHt.putAll(ht);


                }
            }
        }

        conf.setIndexHt(origMs1FileHt);
        // List<String> keys = p.getPeptideKey();
        //  HashMap<String, List<ChroPeptide>> peptideMap = p.getPeptideMap();
        List[] ratioArr = new List[expSize];
        for (int i = 0; i < ratioArr.length; i++) {
            ratioArr[i] = new ArrayList<>();
        }

        for (SampleGroup sampleGroup : sampleGroupList) {

            for (SampleModel sampleModel : sampleGroup.getSampleModelList()) {

                List<String> fnameList = sampleModel.getLabelfreeFilenameList();
                List<String> pathList = sampleModel.getPathList();


                for (int i = 0; i < fnameList.size(); i++) {
                    String eachFile = fnameList.get(i);
                    String eachPath = pathList.get(i);
                    String eachKey = eachFile;

                    if (eachPath.endsWith("/")) {
                        eachKey = eachPath + eachFile;
                    } else {
                        eachKey = eachPath + File.separator + eachFile;
                    }

                    IndexedFile origIFile = origMs1FileHt.get(eachKey);
                    if(null == origIFile) {
                        System.out.println("null error..");
                    }
                 //   getFullSpectra(isoList,origIFile,splitSpectraMap,splitMs1FileHt,);

                }

            }
        }



        Map<String,List<Double>> retTimeListMap = new HashMap<>();
        Map<String,List<Integer>> scanListMap = new HashMap<>();
        Map<String,double []> fileNameSpectraMap = new HashMap<>();
        Map<String,double []> sequenceSpectraMap = new HashMap<>();

        for (LabelfreePeptide each : pepList) {

            double startRt = each.getStartRetTime();
            double endRt = each.getEndRetTime();

            int sampleCount = 0;
            int peptideCount = 0;
            double [] spectraSum = sequenceSpectraMap.get(each.getSequence());


            for (SampleGroup sampleGroup : sampleGroupList) {

                sampleCount++;
                System.out.println("working on sample group\t" + sampleGroup.getName());
                System.out.println("peptide " + each.getSequence() + "\t charge state " + each.getChargeState());
                //         System.out.println("peptide " + each.getSequence() + "\t charge state " + each.getChargeState());


                for (SampleModel sampleModel : sampleGroup.getSampleModelList()) {

//                              System.out.println("sample name\t" + sampleModel.getSampleName());
                    //             if(true) continue;

                    peptideCount++;

                    //for (String path : sampleModel.getPathList()) {
                    //for (int h = 0; h < expSize; h++) {
                    List<String> fnameList = sampleModel.getLabelfreeFilenameList();
                    List<String> pathList = sampleModel.getPathList();


                    for (int i = 0; i < fnameList.size(); i++) {
                        String eachFile = fnameList.get(i);
                        String eachPath = pathList.get(i);
                        String eachKey = eachFile;

                        if (eachPath.endsWith("/")) {
                            eachKey = eachPath + eachFile;
                        } else {
                            eachKey = eachPath + File.separator + eachFile;
                        }

                        IndexedFile origIFile = origMs1FileHt.get(eachKey);
                        if(null == origIFile) {
                            System.out.println("null error..");
                        }

                        // TDoubleIntHashMap retentonToScanMap = origIFile.getRetentonToScanMap();

                        //    int startScan = retentonToScanMap.get(startRt);
                        //  int endScan = retentonToScanMap.get(endRt);

                        List<Double> retList = new ArrayList<>();
                        List<Integer> scanList = new ArrayList<>();
                        double [] spectra = getFullSpectra(isoReader, each.getSequence(),
                                origIFile,
                                splitSpectraMap,
                                splitMs1FileHt, each.getChargeState(),retList,scanList);
                        scanListMap.put(eachKey,scanList);
                        retTimeListMap.put(eachKey,retList);
                       // fileNameSpectraMap.put(eachKey,spectra);
                        if(spectraSum == null)
                        {
              //              spectraSum = spectra;
                        }
                        else
                        {
              //              int length = spectraSum.length< spectra.length? spectraSum.length : spectra.length;
                            for(int j=0 ;j<length; j++ )
                            {
                                spectraSum[j]+=spectra[j];
                            }
                        }
                    }

                }
            }
            sequenceSpectraMap.put(each.getSequence(),spectraSum);
            List<ChroPeptide> chroPepList = each.getPeptideList();

            List<LabelFreeJSONPeptide> jsonPeptides = getEachSamplePeptides(chroPepList, each.getChargeState(), each.getSequence(), startRt, endRt);
            //  peptideCount = peptideCount + jsonPeptides.size();
            jsonAllPeptideList.add(jsonPeptides);
        }

        for(Map.Entry<String,double[]> entry : sequenceSpectraMap.entrySet())
        {
            double[] spectra = entry.getValue();
            int index = 0;
            double max = Double.MIN_VALUE;
            for(int i=0 ;i<spectra.length; i++)
            {
                if(spectra[i]>max)
                {
                    max = spectra[i];
                    index = i;
                }
            }
            System.out.println("max is "+max+" at "+index);
        }

        //  ProteinModel protein = new ProteinModel();
        // protein.setLogRatioArr(ratioArr);

        chroJSONCreator.createJsonForPeptideList(jsonAllPeptideList, jsonPeptideListFileName);

        //System.out.println("Progress: total protein " + totalProtein + "\t" + proteinCount +  "\t" + ((double)proteinCount)/totalProtein*100 + "% complete" );
        //      if(37== proteinCount)
        //        System.out.println("====37");

        System.out.println("100% complete");



        String jsonFile = folder + File.separator + "target.JSON";
        chroJSONCreator.createJsonForProteinIndex(jsonProteins, jsonFile);



        double[] sumIntensityPeptide = null;
        double[] sumIntensityProtein = null;
        boolean check1 = true;
        boolean check2 = true;

        for (int j = 0; j < pepList.size(); j++) {
            LabelfreePeptide pep = pepList.get(j);
            List<ChroPeptide> expList = pep.getPeptideList();
            for (int k = 0; k < expList.size(); k++) {
                if (check1) {
                    sumIntensityPeptide = new double[expList.size()];
                    check1 = false;
                }

                ChroPeptide exp = expList.get(k);
                double sum = sumIntensityPeptide[k];
                sum = sum + (exp.getAverageIntensity() / 1000);
                sumIntensityPeptide[k] = sum;
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



            for(int j=0;j<pepList.size();j++){
                LabelfreePeptide pep = pepList.get(j);
                List<ChroPeptide> expList = pep.getPeptideList();
                for(int k=0;k<expList.size();k++){
                    ChroPeptide exp = expList.get(k);
                    pep.addNormIntensityList(exp.getAverageIntensity()*(avgpeptideIntensity/sumIntensityPeptide[k]));
                }
            }

        String outFile = jsonFile.substring(0, jsonFile.lastIndexOf(File.separator) + 1) + "targeted_peptides.txt";

        //String chroFile = jsonFile.substring(0, jsonFile.lastIndexOf(File.separator)+1) + "chro.txt";

        BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(outFile));
        PrintStream ps = new PrintStream(out);
        int normIndex = -1;

        List<org.jdom.Element> samGroupEleList1 = conf.getRootConfEle().getChildren("sample");
        for(int i=0;i<samGroupEleList1.size();i++){
            ps.print("H\tGROUP_NAME\t"+samGroupEleList1.get(i).getAttributeValue("group")+"\t");
            ps.print(samGroupEleList1.get(i).getChild("each_sample").getAttributeValue("name"));
            ps.print("\n");
        }

        JSONParser parser=new JSONParser();
        Object obj = parser.parse(new FileReader(jsonPeptideListFileName));


        int even=0;

        ps.print("sequence\t"+"chargeState\t"+"AUC1\t"+"AUC2\t"+"isotope_Array\tpeaks\n");

        for (LabelfreePeptide pep : pepList) {
            ps.print(pep.getSequence() + "\t");
            List<ChroPeptide> chroPeptides = pep.getPeptideList();


            ps.print(pep.getChargeState() + "\t");

            for (ChroPeptide cpep : chroPeptides) {
                ps.print(cpep.getPeakArea() + "\t");

            }

            int j=0;

            JSONObject jsonObject =  (JSONObject) obj;
            JSONArray peptideLst = (JSONArray) jsonObject.get("peptideList");

            JSONArray lst=(JSONArray)(peptideLst.get(even));
            for(int k=0;k<lst.size();k++){
                JSONObject lstObj= (JSONObject) lst.get(k);
                JSONArray isotopeArray=(JSONArray) lstObj.get("isotopeArr");
                if(j%2==0) {
                    StringBuilder isoArrBuilder = new StringBuilder();
                    Iterator<java.lang.Double> iterator = isotopeArray.iterator();
                    while (iterator.hasNext()) {
                        isoArrBuilder.append(iterator.next() + ",");
                    }
                    ps.print(isoArrBuilder.substring(0, isoArrBuilder.length() - 1) + "\t");
                }
                ps.print((String)lstObj.get("peaks")+"\t");
                j++;
            }
            even++;
            ps.print("\n");
        }
        /*
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

            for (int j = 0; j < normIndexS; j++) {
                if (words[j].equalsIgnoreCase("NA")) {
                    if (txtReader.getSequenceIndexList().contains(j)) {
                        ps.print(proteinList.get(i).getPeptideList().get(k).getSequence() + "\t");
                    } else if (txtReader.getCsIndexList().contains(j)) {
                        ps.print(proteinList.get(i).getPeptideList().get(k).getChargeState() + "\t");

                    } else if (txtReader.getIntensityIndexList().contains(j)) {
                        ps.print((int) proteinList.get(i).getPeptideList().get(k).getPeptideList().get(l).getAverageIntensity() + "\t");
                        l++;
                    } else {
                        ps.print(words[j] + "\t");
                    }
                } else {
                    if (txtReader.getIntensityIndexList().contains(j)) {
                        l++;
                    }
                    ps.print(words[j] + "\t");
                }

            }
            for(int rr=0;rr<proteinList.get(i).getPeptideList().get(k).getNormIntensityList().size();rr++){
                ps.print(proteinList.get(i).getPeptideList().get(k).getNormIntensityList().get(rr)+"\t");
            }
            ps.print(pvalueList.get(r));
            r++;
            ps.print("\n");
            k++;
        }




        br.close();
        ps.close();

        // LabelfreeMissingPeptideBuilderSplit.proteinCompareBasedOnRatios(proteinList, sample1IndexList, sample2IndexList);

        //return proteinList;
        return null;



    }*/


    public static List<ProteinModel> targetedAnalysis(
            String folder, String peptideFile,
            Configuration conf) throws Exception {

        /*
        Configuration conf = Configuration.getInstance();


        if (!conf.isReadConfigFile()) {
            conf.setLabelfree(true);
            conf.readXMLParam(configFile);
        }
*/
        List<org.jdom.Element> samGroupEleList = conf.getRootConfEle().getChildren("sample");
        int expSize = 0;

        List<Integer> indexList = new ArrayList<>();
        for (Iterator<org.jdom.Element> samgItr = samGroupEleList.iterator(); samgItr.hasNext();) {
            int count =0;
            org.jdom.Element groupEle = samgItr.next();

            List<Element> sampleEleList = groupEle.getChildren("each_sample");

            for (Iterator<Element> samItr = sampleEleList.iterator(); samItr.hasNext();) {
                org.jdom.Element eachSample = samItr.next();
                count++;
                expSize++;
            }
            indexList.add(count);
        }


        IsotopeReader isoReader = new IsotopeReader(conf.getRootConfEle());
        //IdentificationReader idReader = BaseIdentificationReader.getIdentificationInst(isoReader);

     //   conf.setSpectrumFormat(Configuration.MS_FILE_FORMAT);

     //   Hashtable<String, ChroPeptide> peptideChroHt = new Hashtable<>();

        List<SampleModel> allSampleList = conf.getSampleList();
        // sampleList.get(0).getPathList();
	    Hashtable<String, IndexedFile> origMs1FileHt = new Hashtable<>();
        ChroJSONCreator chroJSONCreator = new ChroJSONCreator();
        String jsonFilePath = new File(folder) + File.separator + "JSON_OBJ";
        FileUtil.makeDir(jsonFilePath);
        List<LabelFreeJSONProtein> jsonProteins = new ArrayList<>();
        LabelFreeJSONProtein jsonProtein = null;
        StringBuffer sb = new StringBuffer();
        StringBuffer gPeakSb = new StringBuffer();

        jsonProtein = new LabelFreeJSONProtein();
        jsonProtein.setAccession("Targeted peptides");
        jsonProtein.setDesc("Targeted peptides");
        jsonProteins.add(jsonProtein);


        MSSplitFolderCreation msp = new MSSplitFolderCreation();
        Map<String, String> splitSpectraMap = new HashMap<>();
        Map<String, IndexedFile> splitMs1FileHt = new HashMap<>();
        Map<String, HashMap<Integer, Integer>> ms2ToMs1Map = new HashMap<>();

        List<SampleGroup> sampleGroupList = conf.getSampleGroupList();


        BufferedReader br = new BufferedReader(new FileReader(peptideFile));
        String eachLine;
        List<LabelfreePeptide> pepList = new ArrayList<>();
        while (null != (eachLine = br.readLine())) {


            for (int cs = conf.getTargetedStartCharge(); cs <= conf.getTargetedEndCharge(); cs++) {
                //System.out.println(eachLine);
                LabelfreePeptide pep = new LabelfreePeptide();
                String arr[] = eachLine.split("\t");
                pep.setSequence(arr[0]);
                pep.setChargeState(cs);
                pepList.add(pep);
            }

        }

        br.close();


        List<List<LabelFreeJSONPeptide>> jsonAllPeptideList = new ArrayList<>();
        String jsonPeptideListFileName = jsonFilePath + File.separator + "TARGETED_PEPTIDE.JSON";

        for (SampleGroup sampleGroup : sampleGroupList) {

            for (SampleModel sampleModel : sampleGroup.getSampleModelList()) {

                for (String path : sampleModel.getPathList()) {

                    if (!path.endsWith("/"))
                        path += "/";

                    String spectraPath = path;
                    String splitSpectraPath = path + "split/";

                    System.out.println("===" + path);

                    splitSpectraMap.putAll(msp.splitMS1Files(spectraPath, CensusConstants.LABELFREE_MS1_SPLIT_SCAN_NUM, false));
                    splitMs1FileHt.putAll(ChroGenerator.createIndexedFiles(splitSpectraPath, CensusConstants.MS1_FILE, true, true));
                    ms2ToMs1Map.putAll(IndexUtil.buildMS2toMS1ScanMapFiles(spectraPath));

                    Hashtable<String, IndexedFile> ht = ChroGenerator.createIndexedFiles(spectraPath, "ms1", true, true);

                    origMs1FileHt.putAll(ht);
                }
            }
        }

                    conf.setIndexHt(origMs1FileHt);
                    // List<String> keys = p.getPeptideKey();
                    //  HashMap<String, List<ChroPeptide>> peptideMap = p.getPeptideMap();
                    List[] ratioArr = new List[expSize];
                    for (int i = 0; i < ratioArr.length; i++) {
                        ratioArr[i] = new ArrayList<>();
                    }


                    for (LabelfreePeptide each : pepList) {

                        double startRt = each.getStartRetTime();
                        double endRt = each.getEndRetTime();

                        int sampleCount = 0;
                        int peptideCount = 0;
                        for (SampleGroup sampleGroup : sampleGroupList) {

                            sampleCount++;
                          System.out.println(">>working on sample group\t" + sampleGroup.getName());
                     //     System.out.println(">>peptide " + each.getSequence() + "\t charge state " + each.getChargeState());
                 //         System.out.println("peptide " + each.getSequence() + "\t charge state " + each.getChargeState());


                            for (SampleModel sampleModel : sampleGroup.getSampleModelList()) {

//                              System.out.println("sample name\t" + sampleModel.getSampleName());
                 //             if(true) continue;

                                peptideCount++;

                                //for (String path : sampleModel.getPathList()) {
                                    //for (int h = 0; h < expSize; h++) {
                                List<String> fnameList = sampleModel.getLabelfreeFilenameList();
                                List<String> pathList = sampleModel.getPathList();


                                for (int i = 0; i < fnameList.size(); i++) {
                                    String eachFile = fnameList.get(i);
                                    String eachPath = pathList.get(i);
                                    String eachKey = eachFile;

                                    if (eachPath.endsWith("/")) {
                                        eachKey = eachPath + eachFile;
                                    } else {
                                        eachKey = eachPath + File.separator + eachFile;
                                    }

                                    IndexedFile origIFile = origMs1FileHt.get(eachKey);
                                    if(null == origIFile) {
                                        System.out.println("null error..");
                                    }

                                   // TDoubleIntHashMap retentonToScanMap = origIFile.getRetentonToScanMap();

                                    //    int startScan = retentonToScanMap.get(startRt);
                                    //  int endScan = retentonToScanMap.get(endRt);




                                    GaussianPeakModel peakModel = isotopePeakFinding(isoReader, each.getSequence(),
                                            origIFile,
                                            splitSpectraMap,
                                            splitMs1FileHt, each.getChargeState()
                                    );


                                    ChroPeptide expPep = new ChroPeptide();
                                    expPep.setPeptideIndex(peptideCount);
                                    expPep.setSampleIndex(sampleCount);


                                    if (null != peakModel) {
                                        expPep.setPeakArea(peakModel.getPeakArea());

                                        int[] scanArr = peakModel.getScanArr();
                                        double[] retArr = peakModel.getRetArr();
                                        double[] peakArr = peakModel.getPeakArr();

                                        sb.delete(0, sb.length());
                                        sb.append("P 0 0;");
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
                                        expPep.setIsoArr(peakModel.getIsoArr());
                                        expPep.setFileName(eachFile);
                                        each.addChroPeptide(expPep);
                                    } else {
                                        System.out.println("peak not found");
                                    }

                                }

                            }
                        }

                        List<ChroPeptide> chroPepList = each.getPeptideList();

                        List<LabelFreeJSONPeptide> jsonPeptides = getEachSamplePeptides(chroPepList, each.getChargeState(), each.getSequence(), startRt, endRt);
                        //  peptideCount = peptideCount + jsonPeptides.size();
                        jsonAllPeptideList.add(jsonPeptides);


               /* //  int rcount=0;
                for (int index : sample1IndexList) {
                    double peakArea1 = chroPepList.get(index).getPeakArea();
                    double peakArea2 = chroPepList.get(index + sample1IndexList.size()).getPeakArea();
                    double ratio = Math.log(peakArea1 / peakArea2) / Math.log(2);
                    //  double ratio = peakArea1/peakArea2;
                    //   System.out.print(peakArea1  + "\t" + peakArea2 + "\t");

                    ratioArr[index].add(ratio);
                    //   System.out.println(ratio + "\t");
                }*/


                       // peptideCount = chroPepList.size();

                    }

                  //  ProteinModel protein = new ProteinModel();
                    // protein.setLogRatioArr(ratioArr);

                    chroJSONCreator.createJsonForPeptideList(jsonAllPeptideList, jsonPeptideListFileName);

                    //System.out.println("Progress: total protein " + totalProtein + "\t" + proteinCount +  "\t" + ((double)proteinCount)/totalProtein*100 + "% complete" );
                    //      if(37== proteinCount)
                    //        System.out.println("====37");

                    System.out.println("100% complete");

                    //       System.out.println("p============" + p.getPvalue());
                    //       System.out.println("p============" + p.getPvalue());


                    String jsonFile = folder + File.separator + "target.JSON";
                    chroJSONCreator.createJsonForProteinIndex(jsonProteins, jsonFile);


        /*
                    // adding ANOVA pValue:
                    //List<org.jdom.Element> samGroupEleList = conf.getRootConfEle().getChildren("sample");
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
                    List<Double> pvalueList = new ArrayList<>();

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

                    long totalPeptideIntensity = 0;

                    for (ctr = 0; ctr < size1; ctr++) {
                        int size2 = protein.getPeptideList().get(ctr).getPeptideList().size();
                        //List<Double> intensity1 = new ArrayList<>();
                        // List<Double> intensity2 = new ArrayList<>();
                        List<Double> intensity = new ArrayList<>();
                        List<Double> iitIntensity = new ArrayList<>();

                        Iterator it = sampleListt.iterator();

                        pepAvgIntList = new ArrayList<>();
                        pepIntList = new ArrayList<>();

                        pepAvgIntListIIT = new ArrayList<>();
                        pepIntListIIT = new ArrayList<>();

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
                            long intSumIIT = 0;
                            //while (iterator.hasNext()) {
                            for (int i = 0; i < intensity.size(); i++) {
                                intensityy[t] = intensity.get(i);
                                intSum += (long) intensityy[t];
                                pepIntList.add((long) intensityy[t]);

                                // intensity sum System.out.println(i + "=============" + (long)intensityy[t]);

                                iitIntensityy[t] = iitIntensity.get(i);
                                intSumIIT += (long) iitIntensityy[t];
                                pepIntListIIT.add((long) iitIntensityy[t]);

                                t++;
                            }

                            //intensity sum     System.out.println("sum=============" + (long)intSum);

                            long averageIntensity = (long) intSum / t;
                            long averageIntensityIIT = (long) intSumIIT / t;

                            pepAvgIntList.add(averageIntensity);
                            pepAvgIntListIIT.add(averageIntensityIIT);
                            totalPeptideIntensity += averageIntensity;

                            classes.add(intensityy);
                            intensity.clear();
                            iitIntensity.clear();

                        }

                        intList.add(pepAvgIntList);
                        intIndividualPeptideList.add(pepIntList);

                        intListIIT.add(pepAvgIntListIIT);
                        intIndividualPeptideListIIT.add(pepIntListIIT);

                        double pvalue = AnovaUtil.calculateAnovaPvalue(classes);
                        pvalueList.add(pvalue);
                    }

                    //calculate most common pattern of changes for finding protein intensity
                    int bestIndex = 0;
                    double currentSum = 0;

         *****************************/

//            protein.setTotalPeptideIntensity(totalPeptideIntensity);


                    //adding norm intensity
                    double[] sumIntensityPeptide = null;
                    double[] sumIntensityProtein = null;
                    boolean check1 = true;
                    boolean check2 = true;

                    for (int j = 0; j < pepList.size(); j++) {
                        LabelfreePeptide pep = pepList.get(j);
                        List<ChroPeptide> expList = pep.getPeptideList();
                        for (int k = 0; k < expList.size(); k++) {
                            if (check1) {
                                sumIntensityPeptide = new double[expList.size()];
                                check1 = false;
                            }

                            ChroPeptide exp = expList.get(k);
                            double sum = sumIntensityPeptide[k];
                            sum = sum + (exp.getAverageIntensity() / 1000);
                            sumIntensityPeptide[k] = sum;
                        }
                    }

            /*
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



            for(int j=0;j<pepList.size();j++){
                LabelfreePeptide pep = pepList.get(j);
                List<ChroPeptide> expList = pep.getPeptideList();
                for(int k=0;k<expList.size();k++){
                    ChroPeptide exp = expList.get(k);
                    pep.addNormIntensityList(exp.getAverageIntensity()*(avgpeptideIntensity/sumIntensityPeptide[k]));
                }
            }
*/
                    String outFile = jsonFile.substring(0, jsonFile.lastIndexOf(File.separator) + 1) + "targeted_peptides.txt";

                    //String chroFile = jsonFile.substring(0, jsonFile.lastIndexOf(File.separator)+1) + "chro.txt";

                    BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(outFile));
                    PrintStream ps = new PrintStream(out);
                    int normIndex = -1;

                    List<org.jdom.Element> samGroupEleList1 = conf.getRootConfEle().getChildren("sample");
                    for(int i=0;i<samGroupEleList1.size();i++){
                            ps.print("H\tGROUP_NAME\t"+samGroupEleList1.get(i).getAttributeValue("group")+"\t");
                            ps.print(samGroupEleList1.get(i).getChild("each_sample").getAttributeValue("name"));
                            ps.print("\n");
                    }

                    JSONParser parser=new JSONParser();
                    Object obj = parser.parse(new FileReader(jsonPeptideListFileName));


                    int even=0;

                    ps.print("sequence\t"+"chargeState\t"+"AUC1\t"+"AUC2\t"+"isotope_Array\tpeaks\n");

                    for (LabelfreePeptide pep : pepList) {
                        ps.print(pep.getSequence() + "\t");
                        List<ChroPeptide> chroPeptides = pep.getPeptideList();


                        ps.print(pep.getChargeState() + "\t");

                        for (ChroPeptide cpep : chroPeptides) {
                            ps.print(cpep.getPeakArea() + "\t");

                        }

                        int j=0;

                        JSONObject jsonObject =  (JSONObject) obj;
                        JSONArray peptideLst = (JSONArray) jsonObject.get("peptideList");

                        JSONArray lst=(JSONArray)(peptideLst.get(even));
                        for(int k=0;k<lst.size();k++){
                            JSONObject lstObj= (JSONObject) lst.get(k);
                            JSONArray isotopeArray=(JSONArray) lstObj.get("isotopeArr");
                            if(j%2==0) {
                                StringBuilder isoArrBuilder = new StringBuilder();
                                Iterator<java.lang.Double> iterator = isotopeArray.iterator();
                                while (iterator.hasNext()) {
                                    isoArrBuilder.append(iterator.next() + ",");
                                }
                                ps.print(isoArrBuilder.substring(0, isoArrBuilder.length() - 1) + "\t");
                            }
                            ps.print((String)lstObj.get("peaks")+"\t");
                            j++;
                        }
                        even++;
                        ps.print("\n");
                    }
        /*
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

            for (int j = 0; j < normIndexS; j++) {
                if (words[j].equalsIgnoreCase("NA")) {
                    if (txtReader.getSequenceIndexList().contains(j)) {
                        ps.print(proteinList.get(i).getPeptideList().get(k).getSequence() + "\t");
                    } else if (txtReader.getCsIndexList().contains(j)) {
                        ps.print(proteinList.get(i).getPeptideList().get(k).getChargeState() + "\t");

                    } else if (txtReader.getIntensityIndexList().contains(j)) {
                        ps.print((int) proteinList.get(i).getPeptideList().get(k).getPeptideList().get(l).getAverageIntensity() + "\t");
                        l++;
                    } else {
                        ps.print(words[j] + "\t");
                    }
                } else {
                    if (txtReader.getIntensityIndexList().contains(j)) {
                        l++;
                    }
                    ps.print(words[j] + "\t");
                }

            }
            for(int rr=0;rr<proteinList.get(i).getPeptideList().get(k).getNormIntensityList().size();rr++){
                ps.print(proteinList.get(i).getPeptideList().get(k).getNormIntensityList().get(rr)+"\t");
            }
            ps.print(pvalueList.get(r));
            r++;
            ps.print("\n");
            k++;
        } */




        br.close();
        ps.close();

       // LabelfreeMissingPeptideBuilderSplit.proteinCompareBasedOnRatios(proteinList, sample1IndexList, sample2IndexList);

        //return proteinList;
        return null;



    }

    public static List<LabelFreeJSONPeptide> getEachSamplePeptides(List<ChroPeptide> pepList,
            int chargeState, String sequence, double startRt, double endRt) {

        List<LabelFreeJSONPeptide> eachSamplePeptideList = new ArrayList<>();

        LabelFreeJSONPeptide jsonPeptide = null;

        for (ChroPeptide chroPeptide : pepList) {

            jsonPeptide = new LabelFreeJSONPeptide();
            jsonPeptide.setAuc(Double.toString(chroPeptide.getPeakArea()));
            jsonPeptide.setCount(String.valueOf(chroPeptide.getSampleIndex()));
            jsonPeptide.setUnique(String.valueOf(chroPeptide.isUnique()));
            jsonPeptide.setChro_iso(getFormatedChroISO(formatChroData(chroPeptide.getChroData() == null ? "" : chroPeptide.getChroData())));
            jsonPeptide.setSeq(sequence);
            jsonPeptide.setScan(String.valueOf(chroPeptide.getScanNum()));
            jsonPeptide.setFile(chroPeptide.getFileName() == null ? "" : chroPeptide.getFileName());
            jsonPeptide.setStart_scan(chroPeptide.getStartRange());
            jsonPeptide.setEnd_scan(chroPeptide.getEndRange());
            jsonPeptide.setCharge(String.valueOf(chargeState));
            jsonPeptide.setRowId(String.valueOf(chroPeptide.getSampleIndex()));
            jsonPeptide.setStartRt(String.valueOf(chroPeptide.getStartRt()));
            jsonPeptide.setEndRt(String.valueOf(chroPeptide.getEndRt()));
            jsonPeptide.setRetentionTime(String.valueOf(chroPeptide.getRetentionTime()));
            jsonPeptide.setPeaks(getPeaksNewValue(chroPeptide.getGaussianPeakString()));

	        jsonPeptide.setMaxPeakValue(getHighestPaeakValue(chroPeptide.getGaussianPeakString()));
            jsonPeptide.setIsotopeArr(chroPeptide.getIsoArr());
            //jsonPeptide.setFilename(chroPeptide.getFileName());
            eachSamplePeptideList.add(jsonPeptide);
          //  pepCount++;

            // System.out.println(pepCount);
        }

        return eachSamplePeptideList;
    }
    
    public static String getFormatedChroISO(String chroIso) {

		StringBuffer chroisoDataNewArr = new StringBuffer("");

		try {
			String[] chroisoDataArr = chroIso.split(";");

			for (String peak : chroisoDataArr) {

				String[] valueArr = peak.split(" ");

				if (Double.valueOf(valueArr[2]) > 0.0) {
					chroisoDataNewArr.append(peak).append(";");
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		return chroisoDataNewArr.toString();
	}

	public static String getPeaksNewValue(String peaksVal) {

		StringBuffer peaksNewArr = new StringBuffer("");

		try {

			String[] peaksArr = peaksVal.split(";");

			for (String peak : peaksArr) {
				String[] valueArr = peak.split(" ");
				if (Double.valueOf(valueArr[1]) > 0.0) {
					peaksNewArr.append(peak).append(";");
				}
			}

		} catch (Exception e) {
			// TODO: handle exception
		}

		return peaksNewArr.toString();
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


    public static GaussianPeakModel isotopePeakFinding(//int startIndex, int endIndex,
                                                IsotopeReader isoReader, String sequence,
                                                IndexedFile origIFile,
                                                Map<String, String> splitSpectraMap,
                                                Map<String, IndexedFile> splitMs1FileHt,
                                                int cs
                                                ) throws Exception {

        String origFile = origIFile.getFileName();

      System.out.println("working on " + origFile + "...");
      //  if (sequence.contains(".")) {
      //      sequence = sequence.substring(2, sequence.length() - 2);
      //  }


        /**
         * *********************************
         * 1 calculate isotope distribution *********************************
         */

        char[] ch = sequence.toCharArray();
        /*
        char[] ch = null;

        double ptmMass = 0;
        if(sequence.contains("(")) {
          ptmMass = PhosphoUtil.getTotalPTMMass(sequence);
          String cleanSeq = PhosphoUtil.cleanSequence(sequence);
          ch = cleanSeq.toCharArray();
        } else {
          ch = sequence.toCharArray();
        }
*/

      ElementComposition element = new ElementComposition(ch, 0, ch.length, isoReader.getIsotope());
        element.calculate();

        if (!element.isQuantifiable()) {
            System.out.print("\nError : ");
            System.out.println(sequence + " is not quantifiable.");
            return null;
        }

        Configuration conf = Configuration.getInstance();

        IsotopeDist sampleDist = new IsotopeDist(element.getElementSampleArr(), element.getModShift(), true);

        //  ionInjectionTime = ionInjectionMap.get(tmpScanNumber);
        double[] isoArr = sampleDist.getHighMassList();
        double[] isoIntArr = sampleDist.getRelabun(isoArr.length);
        double pepMass = sampleDist.getHighMassList()[0];



            double[] isoArrCS = new double[isoArr.length];

            for (int i = 0; i < isoArr.length; i++) {
                isoArrCS[i] = (isoArr[i] + cs * CensusConstants.PROTON_MASS) / cs;

            }



      //  System.out.println("uncomment below");
     //   if(true) return null;

            /*
             * *********************************
             * 2 Re-construct chromatogram *********************************
             */

            int[] scanArr = origIFile.getKeys();

            long[] chromPeakArr = new long[scanArr.length];
            double[] retArr = new double[chromPeakArr.length];
            int count = 0;
            boolean hasPeak=false;

            for(int eachScan:scanArr) {

             //   if(eachScan==1051){

                String fileKey = eachScan + "\t" + origFile;
                String spltiMs1File = splitSpectraMap.get(fileKey);
                IndexedFile splitIFile = splitMs1FileHt.get(spltiMs1File);
                try {
                    splitIFile.getIndexByScan(eachScan);
                }
                catch(Exception e) {
                    e.printStackTrace();
                    System.out.println("error");
                }

                int splitCurIndex = splitIFile.getIndexByScan(eachScan);
                long currentPos = splitIFile.getPositionByIndex(splitCurIndex);

                long nextPos = 0;

                if ((splitCurIndex + 1) >= splitIFile.getScanPositionMap().size())
                    nextPos = splitIFile.getFile().length();
                else
                    nextPos = splitIFile.getPositionByIndex(splitCurIndex + 1);

                int diff = (int) (nextPos - currentPos);

                SpectrumModel spec = CalcUtilGeneric.readLabelfreeFullSpectrum(isoArrCS, currentPos, diff,
                        conf.getMassTolerance(), splitIFile, cs, conf, pepMass);


               // System.out.println("each scan\t" + eachScan + "\t" + spec.getPrecursorPeakIntensity());

                chromPeakArr[count] = spec.getPrecursorPeakIntensity();
                if(!hasPeak && chromPeakArr[count]>0)
                    hasPeak=true;

                retArr[count] = spec.getRetentionTime();
                scanArr[count] = spec.getScanNumber();

                count++;
                currentPos = nextPos;
          //  }

            }


                /**
                 * *********************************
                 * 3 Smooth chromatogram *********************************
                 */


            //    System.out.println(hasPeak);
          //      if(!hasPeak) {

                    /*
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
                    if ((splitCurIndex + 1) < splitIFile.getScanPositionMap().size()) {
                        long currentPos = splitIFile.getPositionByIndex(splitCurIndex);
                        long nextPos = splitIFile.getPositionByIndex(splitCurIndex + 1);


                        int diff = (int) (nextPos - currentPos);


                        backGroundNoise = CalcUtilGeneric.getBackGroundNoise(isoArr, currentPos, diff, splitIFile);

                    }


                    gModel.setPeakArea(backGroundNoise);

                    double[] noPeakArr = new double[chromPeakArr.length];
                    gModel.setPeakArr(noPeakArr);

return gModel;
                    */

             //       return null;
               // }

                double[] smoothChromArr = Smooth.smoothAsDouble(chromPeakArr, LabelfreeMissingPeptideBuilderSplit.SMOOTH_WINDOW_SIZE_7);


                double basePeak = 0;
                int basePeakIndex = 0;

                for (int i = 0; i < smoothChromArr.length; i++) {

                    //     System.out.println("==\t" + retArr[i] + "\t" + smoothChromArr[i] + "\t" + chromPeakArr[i]);
                    if (basePeak < smoothChromArr[i]) {
                        basePeak = smoothChromArr[i];
                        basePeakIndex = i;
                    }

              //       System.out.println(smoothChromArr[i]);
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
                gModel.setChargeState(cs);
                gModel.setIsoArr(isoArrCS);





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

    }

}
