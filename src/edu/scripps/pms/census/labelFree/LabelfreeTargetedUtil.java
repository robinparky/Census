package edu.scripps.pms.census.labelFree;

import JSONWriter.JsonWriter;
import edu.scripps.pms.census.CensusConstants;
import edu.scripps.pms.census.ChroGenerator;
import edu.scripps.pms.census.ElementComposition;
import edu.scripps.pms.census.conf.Configuration;
import edu.scripps.pms.census.hash.IndexUtil;
import edu.scripps.pms.census.hash.IndexedFile;
import edu.scripps.pms.census.io.IsotopeReader;
import edu.scripps.pms.census.labelFree.json.ChroJSONCreator;
import edu.scripps.pms.census.labelFree.json.LabelFreeJSONPeptide;
import edu.scripps.pms.census.labelFree.json.LabelFreeJSONProtein;
import edu.scripps.pms.census.labelFree.model.LabelfreePeptide;
import edu.scripps.pms.census.labelFree.util.LabelfreeChroUtil;
import edu.scripps.pms.census.model.ChroPeptide;
import edu.scripps.pms.census.model.SampleGroup;
import edu.scripps.pms.census.model.SampleModel;
import edu.scripps.pms.census.model.Spectrum;
import edu.scripps.pms.census.util.CalcUtilGeneric;
import edu.scripps.pms.census.util.IsotopeDist;
import edu.scripps.pms.census.util.io.FileUtil;
import edu.scripps.pms.census.util.io.SpectrumReader;
import edu.scripps.pms.util.spectrum.PeakList;
import gnu.trove.TDoubleArrayList;
import org.apache.commons.math.stat.descriptive.moment.StandardDeviation;
import org.jdom.Element;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import rpark.statistics.*;
import rpark.statistics.model.GaussianPeakModel;
import scripts.MSSplitFolderCreation;

import java.io.*;
import java.util.*;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static edu.scripps.pms.census.labelFree.LabelfreeTargeted.getEachSamplePeptides;
import static edu.scripps.pms.census.util.CalcUtilGeneric.intensitySumWithIsotopeModeling;

/**
 * Created by Titus Jung titusj@scripps.edu on 5/25/18.
 */
public class LabelfreeTargetedUtil {

    public static final int NUMBER_ADDITIONAL_PEPTIDES = 4;

    public static final Pattern prefixPattern = Pattern.compile("^[A-Z\\-]+\\.");
    public static final Pattern suffixPattern = Pattern.compile("\\.[A-Z\\-]+$");

    public static int [] getPeptideStartEnd(String peptide)
    {
        Matcher prefixMatcher  = prefixPattern.matcher(peptide);
        int start = 0;
        if(prefixMatcher.find())
        {
            start = prefixMatcher.end();
        }

        Matcher suffixMatcher = suffixPattern.matcher(peptide);
        int end =peptide.length();
        if(suffixMatcher.find())
        {
            end = suffixMatcher.start();
        }
        return new int[]{start,end};
    }


    public static void main(String [] args) throws Exception {
        System.out.println("LabelfreeTargeted - new version");

        //String configFile = args[0];
        String configFile = args[0];
        String peptideFile = args[1];
        String folderName = args[2];

        Configuration conf = Configuration.getInstance();
        // conf.setLabelfreeCheckChargeState(true);

        if (!conf.isReadConfigFile()) {
            conf.setLabelfree(true);
            conf.readXMLParam(configFile);
        }
        targetedAnalysisSumPeakArea(folderName,peptideFile,conf);
    }

    public static List<ProteinModel> targetedAnalysisSumPeakArea(
            String folder, String peptideFile,
            Configuration conf) throws Exception {

        List<org.jdom.Element> samGroupEleList = conf.getRootConfEle().getChildren("sample");
        int expSize = 0;

        List<Integer> indexList = new ArrayList<>();
       // System.out.println(retTimeWindow);
        for (Iterator<Element> samgItr = samGroupEleList.iterator(); samgItr.hasNext(); ) {
            int count = 0;
            org.jdom.Element groupEle = samgItr.next();

            List<Element> sampleEleList = groupEle.getChildren("each_sample");

            for (Iterator<Element> samItr = sampleEleList.iterator(); samItr.hasNext(); ) {
                org.jdom.Element eachSample = samItr.next();
                count++;
                expSize++;
            }
            indexList.add(count);
        }


        IsotopeReader isoReader = new IsotopeReader(conf.getRootConfEle());

        List<SampleModel> allSampleList = conf.getSampleList();
        // sampleList.get(0).getPathList();
        Hashtable<String, IndexedFile> origMs1FileHt = new Hashtable<>();
        ChroJSONCreator chroJSONCreator = new ChroJSONCreator();
        String jsonFilePath = new File(folder) + File.separator + "JSON_OBJ";
        FileUtil.makeDir(jsonFilePath);
        File jsonDir = new File(jsonFilePath);
        List<LabelFreeJSONProtein> jsonProteins = new ArrayList<>();
        LabelFreeJSONProtein jsonProtein = null;
        StringBuffer sb = new StringBuffer();
        StringBuffer gPeakSb = new StringBuffer();

        jsonProtein = new LabelFreeJSONProtein();
        jsonProtein.setAccession("Targeted peptides");
        jsonProtein.setDesc("Targeted peptides");
        jsonProteins.add(jsonProtein);


        File [] toDelete =  jsonDir.listFiles(new FilenameFilter() {
            @Override
            public boolean accept(File dir, String name) {
                return name.startsWith("temp_") && name.endsWith(".JSON");
            }
        });
        for(File f: toDelete)
        {
            f.delete();
        }


        MSSplitFolderCreation msp = new MSSplitFolderCreation();
        Map<String, String> splitSpectraMap = new HashMap<>();
        Map<String, IndexedFile> splitMs1FileHt = new HashMap<>();
        Map<String, HashMap<Integer, Integer>> ms2ToMs1Map = new HashMap<>();

        List<SampleGroup> sampleGroupList = conf.getSampleGroupList();


        BufferedReader br = new BufferedReader(new FileReader(peptideFile));
        String eachLine;

        Map<String, List<LabelfreePeptide>> sequencePeptideMap = new HashMap<>();
        Map<String, LabelfreePeptide> seqKeyPeptideMap = new TreeMap<>();
        List<LabelfreePeptide> pepList = new ArrayList<>();

        int numGroups = sampleGroupList.size();
        while (null != (eachLine = br.readLine())) {

            String arr[] = eachLine.split("\t");
            String sequence = arr[0];
            if(sequence.length()>0)
            {
                List<LabelfreePeptide> peptideList = new ArrayList<>();

                for (int cs = conf.getTargetedStartCharge(); cs <= conf.getTargetedEndCharge(); cs++) {
                    //System.out.println(eachLine);
                    edu.scripps.pms.census.labelFree.model.LabelfreePeptide pep = new edu.scripps.pms.census.labelFree.model.LabelfreePeptide();

                    pep.setSequence(arr[0]);
                    pep.setChargeState(cs);

                    char[] ch = sequence.toCharArray();


                    int [] startEnd = getPeptideStartEnd(sequence);
                    // System.out.println("<<> "+sequence +"\t"+sequence.substring(startEnd[0],startEnd[1]));

                    ElementComposition element = new ElementComposition(ch, startEnd[0], startEnd[1]-startEnd[0], isoReader.getIsotope());
                    element.calculate();

                    if (!element.isQuantifiable()) {
                        System.out.print("\nError : ");
                        System.out.println(sequence + " is not quantifiable.");
                        return null;
                    }
                    IsotopeDist sampleDist = new IsotopeDist(element.getElementSampleArr(), element.getModShift(), true);
                    pep.setIsotopeDist(sampleDist);
                    peptideList.add(pep);
                    pepList.add(pep);
                    String seqkey = sequence+cs;
                    seqKeyPeptideMap.put(seqkey,pep);
                }
                sequencePeptideMap.put(sequence, peptideList);
            }


        }

        br.close();


        List<List<LabelFreeJSONPeptide>> jsonAllPeptideList = new ArrayList<>();
        String jsonPeptideListFileName = jsonFilePath + File.separator + "TARGETED_PEPTIDE_OLD.JSON";

        for (SampleGroup sampleGroup : sampleGroupList) {

            for (SampleModel sampleModel : sampleGroup.getSampleModelList()) {

                for (String path : sampleModel.getPathList()) {
/*
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

*/
                }
            }
        }
        conf.setIndexHt(origMs1FileHt);

        double massTolerance = conf.getMassTolerance();
        int sampleCount = 0;
        int peptideCount = 0;
        Map<String, List<LabelFreeJSONPeptide>> sequenceJsonMap = new TreeMap<>();
        Map<String, BufferedWriter> tempBwList = new HashMap<>();

        Map<String,Integer> sampleFileNameMap = new HashMap<>();
        for (SampleGroup sampleGroup : sampleGroupList) {

            System.out.println("working on sample group\t" + sampleGroup.getName());
            //         System.out.println("peptide " + each.getSequence() + "\t charge state " + each.getChargeState());


            for (SampleModel sampleModel : sampleGroup.getSampleModelList()) {
                List<String> fnameList = sampleModel.getLabelfreeFilenameList();
                List<String> pathList = sampleModel.getPathList();


                for (int i = 0; i < fnameList.size(); i++) {
                    String eachFile = fnameList.get(i);
                    String eachPath = pathList.get(i);
                    String eachKey = eachFile;
                    String sampleFileNameKey =sampleGroup.ID+eachFile;
                    //sampleFileNameMap.put(sampleFileNameKey,sampleGroup.ID);
                    //sampleFileNameList.add(sampleFileNameKey);
                 /*   if (eachPath.endsWith("/")) {
                        eachKey = eachPath + eachFile;
                    } else {
                        eachKey = eachPath + File.separator + eachFile;
                    }

                    IndexedFile origIFile = origMs1FileHt.get(eachKey);
                    if (null == origIFile) {
                        System.out.println("null error..");
                    }

                    int[] scanArr = origIFile.getKeys();


                    int count = 0;
                    boolean hasPeak=false;
                    double[] retArr = new double[scanArr.length];
*/

                    SpectrumReader reader = new SpectrumReader(eachPath+File.separator+eachFile,"ms1");
                    int count =0;
                    int[] scanArr = new int[reader.getNumSpectra()];
                    double[] retArr = new double[reader.getNumSpectra()];
                    for(Iterator<PeakList> itr = reader.getSpectra(); itr.hasNext();  )
                    {
                    /*
                    for(int eachScan:scanArr) {

                        //   if(eachScan==1051){

                        String fileKey = eachScan + "\t" + origIFile.getFileName();
                        String spltiMs1File = splitSpectraMap.get(fileKey);
                        IndexedFile splitIFile = splitMs1FileHt.get(spltiMs1File);
                        try {
                            splitIFile.getIndexByScan(eachScan);
                        } catch(Exception e) {
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
                        SpectrumModel model =CalcUtilGeneric.readLabelfreeFullSpectrum(currentPos,diff,splitIFile);
                        */
                        PeakList peakList = itr.next();
                        SpectrumModel model =CalcUtilGeneric.readLabelfreeFullSpectrum(peakList);

                        scanArr[count] = model.getScanNumber();
                        retArr[count] = model.getRetentionTime();
                        double [] massArr = model.getMass();
                        double [] intArr = model.getIntensity();
                        int [] csArr = model.getCsArray();


                        for(Map.Entry<String,List<LabelfreePeptide>> entry: sequencePeptideMap.entrySet())
                        {

                            String sequence = entry.getKey();

                            List<LabelfreePeptide> isotopeDistsList = entry.getValue();
                            for(LabelfreePeptide lpeptide: isotopeDistsList)
                            {

                                long[] chromPeakArr = lpeptide.getChromeSpectra();
                                if(chromPeakArr == null) chromPeakArr= new long[scanArr.length];

                                IsotopeDist sampleDist = lpeptide.getIsotopeDist();
                                int cs = lpeptide.getChargeState();
                                double [] massArrC = massArr;
                                double [] intArrC = intArr;
                                if(csArr!=null)
                                {
                                    double [] tempMassArr = new double[massArr.length];
                                    double [] tempIntArr = new double[massArr.length];
                                    for(int k=0 ;k <massArr.length; k++)
                                    {
                                        if(conf.isLabelfreeCheckChargeState() &&  csArr[k]!=cs)
                                        {
                                            tempMassArr[k] = 0;
                                            tempIntArr[k]=0;
                                        }
                                        else
                                        {
                                            tempMassArr[k] = massArr[k];
                                            tempIntArr[k]=intArr[k];
                                        }
                                    }
                                    massArrC = tempMassArr;
                                    intArrC = tempIntArr;
                                }
                                double pepMass = sampleDist.getHighMassList()[0];

                                //double[] retArr = new double[chromPeakArr.length];
                               // entry.getValue();
                                double[] isoArr = sampleDist.getHighMassList();

                                double[] isoArrCS = new double[isoArr.length];

                                for (int j = 0; j < isoArr.length; j++) {
                                    isoArrCS[j] = (isoArr[j] + cs * CensusConstants.PROTON_MASS) / cs;

                                }
                                double[] tempArr = intensitySumWithIsotopeModeling(massArrC, intArrC, isoArrCS, massTolerance, pepMass);
                                lpeptide.setIsoArrCS(isoArrCS);
                                long peakIntensity = (long) tempArr[0];
                                chromPeakArr[count] = peakIntensity;
                          //      String seq = entry.getKey();
                         //       List<LabelfreePeptide> list = entry.getValue();
                                lpeptide.setChromeSpectra(chromPeakArr);
                                //  double[] smoothChromArr = Smooth.smoothAsDouble(chromPeakArr, LabelfreeMissingPeptideBuilderSplit.SMOOTH_WINDOW_SIZE_7);
                            }
                        }
                        count++;
                    }

                    Map<String,long []> peptidePeakArrMap = new HashMap<>();
                    for (Map.Entry<String, List<LabelfreePeptide>> entry : sequencePeptideMap.entrySet()) {
                            List<LabelfreePeptide> lPeptideList = entry.getValue();
                            long[] chromPeakArr = peptidePeakArrMap.get(entry.getKey());
                            if(chromPeakArr == null) chromPeakArr= new long[scanArr.length];
                            for(LabelfreePeptide labelfreePeptide :lPeptideList)
                            {
                                long[] peakArr =labelfreePeptide.getChromeSpectra();
                                for(int jj=0; jj<peakArr.length; jj++)
                                {
                                    chromPeakArr[jj] += peakArr[jj];
                                }
                            }
                            peptidePeakArrMap.put(entry.getKey(),chromPeakArr);
                    }

                    for(Map.Entry<String,long[]> entry: peptidePeakArrMap.entrySet())
                    {
                        long [] isoArr = entry.getValue();
                        double[] smoothChromArr = Smooth.smoothAsDouble(isoArr, LabelfreeMissingPeptideBuilderSplit.SMOOTH_WINDOW_SIZE_7);
                        double basePeak = 0;
                        int basePeakIndex = 0;

                        for (int jj = 0; jj < smoothChromArr.length; jj++) {

                            //     System.out.println("==\t" + retArr[i] + "\t" + smoothChromArr[i] + "\t" + chromPeakArr[i]);
                            if (basePeak < smoothChromArr[jj]) {
                                basePeak = smoothChromArr[jj];
                                basePeakIndex = jj;
                            }

                            //       System.out.println(smoothChromArr[i]);
                        }

                        List<LabelfreePeptide> lPeptideList = sequencePeptideMap.get(entry.getKey());

                        for(LabelfreePeptide labelfreePeptide :lPeptideList)
                        {
                            int cs = labelfreePeptide.getChargeState();
                            long[] spectra = labelfreePeptide.getChromeSpectra();
                            smoothChromArr = Smooth.smoothAsDouble(spectra, LabelfreeMissingPeptideBuilderSplit.SMOOTH_WINDOW_SIZE_7);
                            int[] indexResult = LabelfreeChroUtil.getPeakRange(basePeakIndex, basePeak, smoothChromArr);
                            int peakStartIndex = indexResult[0];
                            int peakEndIndex = indexResult[1];
                           // sampleCount++;
                           // peptideCount++;

                            List<ChroPeptide> chroPepList = new ArrayList<>();


                            List<int[]> topNPeaks = findTopNPeakRanges(smoothChromArr,indexResult,NUMBER_ADDITIONAL_PEPTIDES);
                            topNPeaks.add(0,indexResult);

                            for(int [] indices: topNPeaks)
                            {
                                ChroPeptide expPep = generatePeakModel(retArr,smoothChromArr,indices[0],indices[1],scanArr,cs,labelfreePeptide);
                                expPep.setPeptideIndex(peptideCount);

                                expPep.setSampleIndex(sampleCount);
                                expPep.setFileName(eachFile);
                                chroPepList.add(expPep);

                            }
                            // chroPepList = labelfreePeptide.getPeptideList();

                            List<LabelFreeJSONPeptide> jsonPeptides = getEachSamplePeptides(chroPepList, labelfreePeptide.getChargeState(),
                                    labelfreePeptide.getSequence(), labelfreePeptide.getStartRetTime(), labelfreePeptide.getEndRetTime());
                            //  peptideCount = peptideCount + jsonPeptides.size();
                            labelfreePeptide.setChromeSpectra(null);
                            String seqkey = entry.getKey()+cs;
                           // List<LabelFreeJSONPeptide> jsonList = sequenceJsonMap.get(seqkey);

                       //     BufferedWriter bw  = tempBwList.get(seqkey);
                            BufferedWriter bw = new BufferedWriter( new FileWriter(jsonFilePath + File.separator +"temp_"+seqkey+".JSON",true));
                           /* if(bw == null)
                            {
                                bw = new BufferedWriter( new FileWriter(jsonFilePath + File.separator +"temp_"+seqkey+".JSON",true));
                                tempBwList.put(seqkey,bw);
                            }*/
                            File tempFile = new File(jsonFilePath + File.separator +"temp_"+seqkey+".JSON");

                            for(int j =0; j<jsonPeptides.size(); j++)
                            {
                                if(tempFile.length()>0)
                                {
                                    bw.append(",");
                                    bw.newLine();
                                }
                                LabelFreeJSONPeptide lfp = jsonPeptides.get(j);
                                lfp.setSampleName(sampleModel.getSampleName());
                                lfp.setGroupName(sampleGroup.getName());
                                lfp.setRowId(Integer.toString(sampleGroup.ID));
                                lfp.setRank(j);
                                chroJSONCreator.createJsonForPeptideIndex(lfp,bw);
                                bw.flush();
                                //if(j<jsonPeptides.size()-1)bw.append(",");
                                //bw.newLine();
                            }
                            bw.flush();
                            bw.close();


                           /* if(jsonList!=null)
                            {
                                jsonList.addAll(jsonPeptides);
                            }
                            else
                            {
                                jsonList = jsonPeptides;
                            }*/
                            //System.out.println(seqkey);

                          //  sequenceJsonMap.put(seqkey,jsonList);
                            //labelfreePeptide.setPeptideList(new ArrayList<ChroPeptide>());
                            //
                            //
                            // jsonAllPeptideList.add(jsonPeptides);
                        }
                    }
                    reader.closeDataFile();

                }
            }
        }
        /*
        for(Map.Entry<String,List<LabelFreeJSONPeptide>> entry : sequenceJsonMap.entrySet() )
        {
            jsonAllPeptideList.add(entry.getValue());
        }*/

      /*  for(BufferedWriter bw: tempBwList.values())
        {
            bw.close();
        }*/




        BufferedWriter bw = new BufferedWriter(new FileWriter(jsonPeptideListFileName));
        bw.append("{");
        bw.newLine();
        bw.append("\t\"peptideList\":[");
        bw.newLine();
        File [] tempFiles =  jsonDir.listFiles(new FilenameFilter() {
            @Override
            public boolean accept(File dir, String name) {
                return name.startsWith("temp_") && name.endsWith(".JSON");
            }
        });
        Arrays.sort(tempFiles, new Comparator<File>() {
            @Override
            public int compare(File o1, File o2) {
                return o1.getName().compareTo(o2.getName());
            }
        });

        for(int i=0; i<tempFiles.length; i++)
        {
            File f = tempFiles[i];
            bw.append("\t\t[");
            bw.newLine();
            BufferedReader btempR = new BufferedReader(new FileReader(f));
            String line;
            while((line= btempR.readLine())!=null)
            {
                bw.append("\t\t\t").append(line);
                bw.newLine();
            }
            if(i<tempFiles.length-1){
                bw.append("\t\t],");
            }
            else
            {
                bw.append("\t\t]");
            }
            bw.newLine();
            btempR.close();
        }
        bw.append("\t]");
        bw.newLine();
        bw.append("}");
        bw.close();

        toDelete =  jsonDir.listFiles(new FilenameFilter() {
            @Override
            public boolean accept(File dir, String name) {
                return name.startsWith("temp_") && name.endsWith(".JSON");
            }
        });
        for(File f: toDelete)
        {
            f.delete();
        }



        // chroJSONCreator.createJsonForPeptideList(jsonAllPeptideList, jsonPeptideListFileName);

        System.out.println("100% complete");

        String jsonFile = folder + File.separator + "target.JSON";
        chroJSONCreator.createJsonForProteinIndex(jsonProteins, jsonFile);


        double[] sumIntensityPeptide = null;
        double[] sumIntensityProtein = null;
        boolean check1 = true;
        boolean check2 = true;

      /*  for (int j = 0; j < pepList.size(); j++) {
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
        }*/


        String outFile = jsonFile.substring(0, jsonFile.lastIndexOf(File.separator) + 1) + "targeted_peptides.txt";


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

        //ps.print("sequence\t"+"chargeState\t"+"AUC1\t"+"AUC2\t"+"isotope_Array\tpeaks\n");
        ps.print("H\tsequence\tchargeState\t");
        for(int i=0; i<numGroups; i++)
        {
            int printI =i +1;
            ps.print("AUC");
            ps.print(printI);
            ps.print("\t");
        }

        for(int i=0; i<numGroups; i++)
        {
            int printI =i +1;
            ps.print("RT_");
            ps.print(printI);
            ps.print("\t");
        }

        for(int i=0; i<numGroups; i++)
        {
            int printI =i +1;
            ps.print("isotope_Array_");
            ps.print(printI);
            ps.print("\t");
        }

        for(int i=0; i<numGroups; i++)
        {
            int printI =i +1;
            ps.print("peaks_");
            ps.print(printI);
            ps.print("\t");
        }
        ps.print("\n");

      //  ps.print("isotope_Array\tpeaks\n");


        JSONObject jsonObject =  (JSONObject) obj;
        JSONArray peptideLst = (JSONArray) jsonObject.get("peptideList");

        for(int j=0; j<peptideLst.size(); j++)
        {

            JSONArray lst=(JSONArray)(peptideLst.get(j));
            String sequence = "";
            Map<String,List<String>> aucMap = new HashMap<>();
            Map<String,List<String>> fileNameMap = new HashMap<>();
            Map<String,List<String>> rtMap = new HashMap<>();
            Map<String,List<String>> peakMap = new HashMap<>();
            Map<String,List<String>> isoMap = new HashMap<>();

            List<String> aucList = new ArrayList<>();
            List<String> fileNamelist = new ArrayList<>();
            List<String> rtList = new ArrayList<>();
            String cs = "";
            List<String> peakList = new ArrayList<>();
            List<String> sampleFileNameList = new ArrayList<>();
            for(int k=0;k<lst.size();k++){
                JSONObject lstObj= (JSONObject) lst.get(k);
                JSONArray isotopeArray=(JSONArray) lstObj.get("isotopeArr");
                String auc = (String)lstObj.get("auc");

                sequence = (String)lstObj.get("seq");

                String fileName = (String)lstObj.get("file");
                String sampleName = (String)lstObj.get("sampleName");
                String groupName = (String)lstObj.get("rowId");
                String key =groupName;
           //     System.out.println("<><>>"+key);
                sampleFileNameList.add(key);


                String startRt = (String)lstObj.get("startRt");
                String endRt = (String)lstObj.get("endRt");
                 cs = (String)lstObj.get("charge");
                String rtRange = startRt+":"+endRt;



                // if(j%2==0) {
                StringBuilder isoArrBuilder = new StringBuilder();
                Iterator<java.lang.Double> iterator = isotopeArray.iterator();
                while (iterator.hasNext()) {
                    isoArrBuilder.append(iterator.next() + ",");
                }
                String isoArrayString = isoArrBuilder.toString();
                String peaks = (String)lstObj.get("peaks");


                if(peakMap.get(key) == null)
                {
                    List<String> temp = new ArrayList<String>();
                    temp.add(peaks);
                    peakMap.put(key,temp);

                    temp = new ArrayList<String>();
                    temp.add(auc);
                    aucMap.put(key,temp);

                    temp = new ArrayList<String>();
                    temp.add(fileName);
                    fileNameMap.put(key,temp);


                    temp = new ArrayList<String>();
                    temp.add(rtRange);
                    rtMap.put(key,temp);


                    temp = new ArrayList<String>();
                    temp.add(rtRange);
                    rtMap.put(key,temp);

                    temp = new ArrayList<String>();
                    temp.add(isoArrayString);
                    isoMap.put(key,temp);

                }
                else
                {
                    peakMap.get(key).add(peaks);
                    aucMap.get(key).add(auc);
                    fileNameMap.get(key).add(fileName);
                    rtMap.get(key).add(rtRange);
                    isoMap.get(key).add(isoArrayString);
                }

               // ps.print(isoArrBuilder.substring(0, isoArrBuilder.length() - 1) + "\t");
                // }
               // ps.print((String)lstObj.get("peaks")+"\t");
            }
            ps.print("S\t");
            ps.print(sequence);
            ps.print("\t");
            ps.print(cs);
            ps.print("\t");
            StringBuilder aucBuilder = new StringBuilder();
            StringBuilder rtBuilder = new StringBuilder();
            StringBuilder peakBuilder = new StringBuilder();
            StringBuilder isoBuilder = new StringBuilder();

            for (SampleGroup sampleGroup : sampleGroupList)
            {
                int i = sampleGroup.ID;
                //System.out.println("<<>><> "+i);
                String key = sampleFileNameList.get(i);
                List<String> temp = aucMap.get(key);
                for(String s :temp)
                {
                    aucBuilder.append(s);
                    aucBuilder.append(",");
                }
                aucBuilder.append("\t");

                temp = rtMap.get(key);
                for(String s :temp)
                {
                    rtBuilder.append(s);
                    rtBuilder.append(",");
                }
                rtBuilder.append("\t");

                temp = peakMap.get(key);
                for(String s :temp)
                {
                    peakBuilder.append(s);
                    peakBuilder.append(",");
                }
                peakBuilder.append("\t");

                temp = isoMap.get(key);
                for(String s :temp)
                {
                    isoBuilder.append(s);
                    isoBuilder.append(",");
                }
                isoBuilder.append("\t");

            }
            ps.print(aucBuilder.toString());
            ps.print(rtBuilder.toString());
            ps.print(isoBuilder.toString());
            ps.print(peakBuilder.toString());

            ps.println();

        }
        String outJSONFile = jsonFilePath + "/TARGETED_PEPTIDE.JSON";
        createAUCTable(jsonPeptideListFileName,outJSONFile,sampleFileNameMap);
       // ps.close();


        /*for (LabelfreePeptide pep : seqKeyPeptideMap.values()) {
            Map<String,StringBuilder>  aucMap = new HashMap<>();
            ps.print(pep.getSequence() + "\t");
            List<ChroPeptide> chroPeptides = pep.getPeptideList();


            ps.print(pep.getChargeState() + "\t");

            for (ChroPeptide cpep : chroPeptides) {
                 aucBuilder = aucMap.get(cpep.getFileName());
                if(aucBuilder == null) aucBuilder = new StringBuilder();
                ps.print(cpep.getPeakArea() + "\t");
                aucBuilder.append(cpep.getPeakArea()).append(";");

                aucMap.put(cpep.getFileName(),aucBuilder);
            }



            int j=0;

            JSONObject jsonObject =  (JSONObject) obj;
            JSONArray peptideLst = (JSONArray) jsonObject.get("peptideList");

            JSONArray lst=(JSONArray)(peptideLst.get(even));
            for(int k=0;k<lst.size();k++){
                JSONObject lstObj= (JSONObject) lst.get(k);
                JSONArray isotopeArray=(JSONArray) lstObj.get("isotopeArr");
               // if(j%2==0) {
                    StringBuilder isoArrBuilder = new StringBuilder();
                    Iterator<java.lang.Double> iterator = isotopeArray.iterator();
                    while (iterator.hasNext()) {
                        isoArrBuilder.append(iterator.next() + ",");
                    }
                    ps.print(isoArrBuilder.substring(0, isoArrBuilder.length() - 1) + "\t");
               // }
                ps.print((String)lstObj.get("peaks")+"\t");
                j++;
            }
            even++;
            ps.print("\n");
        }*/


        br.close();
        out.close();
        return null;
    }

    public static List<int[]> findTopNPeakRanges(double [] smoothChromArr,
                                                 int[] maxPeakIndices, int numPeaks )
    {
        List<int[]> results = new ArrayList<>();
        Set<Integer> excludePeaks = new HashSet<>();
       // System.out.println(1+" place peak is at "+maxPeakIndices[0]+" with intensity "+maxPeakIndices[1]);

        for(int i=maxPeakIndices[0]; i<=maxPeakIndices[1]; i++)
        {

            excludePeaks.add(i);
        }
        for(int i=0; i<numPeaks; i++)
        {
            double basePeak =0;
            int basePeakIndex =-1;
            for(int j=0; j<smoothChromArr.length; j++)
            {
                if(basePeak<smoothChromArr[j] && !excludePeaks.contains(j))
                {
                    basePeak = smoothChromArr[j];
                    basePeakIndex = j;
                }
            }
            if(basePeakIndex<0) break;
            int place = i+2;
         //   System.out.println(place+" place peak is at "+basePeakIndex+" with intensity "+basePeak);
      //      System.out.println("exclude set size is "+excludePeaks.size());
            int[] indexResult = LabelfreeChroUtil.getPeakRange(basePeakIndex, basePeak, smoothChromArr);
            results.add(indexResult);
            for(int j=indexResult[0]; j<=indexResult[1]; j++)
            {
                excludePeaks.add(j);
            }
        }


        return results;
    }


    public static ChroPeptide generatePeakModel(double [] retArr, double[] smoothChromArr,int peakStartIndex, int peakEndIndex,
                                               int[] scanArr, int cs,LabelfreePeptide labelfreePeptide )
    {
        ChroPeptide expPep = new ChroPeptide();

        GaussianPeakModel gModel = GaussianFitting.getGaussianPeakRangeIndex(retArr, smoothChromArr, peakStartIndex, peakEndIndex);
        gModel.setScanArr(scanArr);
        gModel.setRetArr(retArr);
        gModel.setPeakArr(smoothChromArr);
        gModel.setChargeState(cs);
        gModel.setIsoArr(labelfreePeptide.getIsoArrCS());
        expPep.setStartRt(retArr[peakStartIndex]);
        expPep.setEndRt(retArr[peakEndIndex]);
        StringBuffer gPeakSb = new StringBuffer();
        StringBuilder sb  = new StringBuilder();
        if (null != gModel) {
            expPep.setPeakArea(gModel.getPeakAreaTargeted());
            expPep.setRetentionTime(gModel.getX());

            double[] peakArr = gModel.getPeakArr();


            sb.append("P 0 0;");
            for (int j = 0; j < scanArr.length; j++) {
                sb.append(scanArr[j]).append(" ").append(retArr[j]).append(" ").append(peakArr[j]).append(";");
            }

            double[] gxArr = gModel.getGaussianXArr();
            double[] gyArr = gModel.getGaussianYArr();
            gPeakSb.delete(0, gPeakSb.length());


            //                                System.out.println(java.util.Arrays.toString(gyArr));
            if (null != gxArr) {
                for (int j = 0; j < gxArr.length; j++) {
                    gPeakSb.append(gxArr[j]).append(" ").append(gyArr[j]).append(";");
                }
                if(gxArr.length>0)
                {
                    expPep.setStartRt(gxArr[0]);
                    expPep.setEndRt(gxArr[gxArr.length-1]);
                }
            }

            expPep.setChroData(sb.toString());
            expPep.setGaussianPeakString(gPeakSb.toString());
            expPep.setPeakSigma(gModel.getSigma());
            expPep.setPeakx(gModel.getX());
            expPep.setPeaky(gModel.getY());
            expPep.setIsoArr(gModel.getIsoArr());

           // labelfreePeptide.addChroPeptide(expPep);
        } else {
            System.out.println("peak not found");
            return null;
        }
        return expPep;
    }

    public static void createAUCTable(String jsonPath, String outputPath,
                                       Map<String,Integer> sampleFileNameMap) throws Exception {
        JSONParser parser=new JSONParser();
        List<String> sampleFileNameList = new ArrayList<>();
        Set<Integer> idSet = new TreeSet<>();
        for(String s: sampleFileNameMap.keySet())
        {
            sampleFileNameList.add(s);
            idSet.add(sampleFileNameMap.get(s));
        }
        Object obj = parser.parse(new FileReader(jsonPath));
        JSONObject jsonObject =  (JSONObject) obj;
        JSONArray peptideLst = (JSONArray) jsonObject.get("peptideList");


        Map<String,List<AnalysisPeptide>> seqKeyPepList = new TreeMap<>();

        for(int j=0; j<peptideLst.size(); j++) {
            JSONArray lst = (JSONArray) (peptideLst.get(j));
            for(int i=0; i<lst.size(); i++)
            {
                String sequence = "";
                String cs = "";
                JSONObject lstObj = (JSONObject) lst.get(i);
                long rank = (long) lstObj.get("rank");
                if(rank == 0)
                {
                    String aucStr = (String) lstObj.get("auc");
                    double auc = Double.parseDouble(aucStr);
                    sequence = (String) lstObj.get("seq");
                    String fileName = (String) lstObj.get("file");
                    String sampleName = (String) lstObj.get("sampleName");
                    String rowId = (String) lstObj.get("rowId");
                    int groupName = Integer.parseInt(rowId);

                    String isotopeArray=(String) lstObj.get("chro_iso");
                    String startRt = (String)lstObj.get("startRt");
                    String endRt = (String)lstObj.get("endRt");
                    String peaks = (String)lstObj.get("peaks");
                    cs = (String) lstObj.get("charge");
                    String key = sequence+cs;
                    List<AnalysisPeptide> peptideList = seqKeyPepList.getOrDefault(key,new ArrayList<>());

                    AnalysisPeptide peptide = new AnalysisPeptide(sequence,cs,auc,groupName,sampleName,fileName);
                    peptide.setEndRt(endRt);
                    peptide.setStartRt(startRt);
                    peptide.setPeaks(peaks);
                    peptide.setChroIso(isotopeArray);

                    peptideList.add(peptide);
                    seqKeyPepList.put(key,peptideList);
                }
            }
        }
        List<JSONObject> outputList = new ArrayList<>();
        Map<Integer,List<Double>> groupPValueMap = new HashMap<>();
        List<Double> pvalue3_2 = new ArrayList<>();
        Map<Integer,List<JSONObject>> groupOutputMap = new HashMap<>();
        List<JSONObject> outputList3_2 = new ArrayList();

        for(Map.Entry<String,List<AnalysisPeptide>> entry : seqKeyPepList.entrySet())
        {
            List<AnalysisPeptide> peptideList = entry.getValue();
            String sequence = peptideList.get(0).sequence;
            String cs = peptideList.get(0).cs;
            Map<Integer,TDoubleArrayList> aucGroup = new TreeMap<>();
            Map<Integer,List<AnalysisPeptide>> peptideMap = new HashMap<>();


            for(AnalysisPeptide peptide :peptideList)
            {
                TDoubleArrayList aucList = aucGroup.getOrDefault(peptide.groupName,new TDoubleArrayList());
                aucList.add(peptide.auc);
                List<AnalysisPeptide> peptideList1 = peptideMap.getOrDefault(peptide.groupName,new ArrayList<>());
                peptideList1.add(peptide);
                peptideMap.put(peptide.groupName,peptideList1);
                aucGroup.put(peptide.groupName,aucList);

            }
            double [] firstArr = aucGroup.get(0).toNativeArray();

            JSONObject outputJson = new JSONObject();

            outputJson.put("sequence",sequence);
            outputJson.put("cs",cs);
            outputJson.put("rowId",0);

            JSONArray peakArray = new JSONArray();
            JSONArray chroArray = new JSONArray();
            JSONArray startRTArray = new JSONArray();
            JSONArray endRTArray = new JSONArray();
            JSONArray fileNameArray = new JSONArray();



            for(int i : aucGroup.keySet())
            {
                StringBuilder aucBuilder = new StringBuilder();
                List<AnalysisPeptide> peptides = peptideMap.get(i);
                for(AnalysisPeptide peptide: peptides)
                {
                    peakArray.add(peptide.getPeaks());
                    chroArray.add(peptide.getChroIso());
                    startRTArray.add(peptide.getStartRt());
                    endRTArray.add(peptide.getEndRt());
                    fileNameArray.add(peptide.sampleName);
                }
                double [] nextArr = aucGroup.get(i).toNativeArray();
                for(double d: nextArr)
                {
                    aucBuilder.append(d).append(",");
                }
                String aucHeader = "AUC_"+i;
                outputJson.put(aucHeader,aucBuilder.toString());
                if(i!=0)
                {
                    AnalyzeArrayPair pair = new AnalyzeArrayPair(firstArr,nextArr);
                    double pvalue = pair.getPvalue();
                    double ratio = pair.getAverageRatio();
                    int printI = i+1;

                    String ratioHeader = "AVERAGE_RATIO_"+printI+"_1";
                    String pvalueHeader = "PVALUE_"+printI+"_1";
                    outputJson.put(ratioHeader,ratio);
                    String qvalueHeader = "QVALUE_"+printI+"_1";

                    if(!Double.isNaN(pvalue) && !Double.isInfinite(pvalue))
                    {
                        outputJson.put(pvalueHeader,pvalue);
                        List<Double> pvalueList = groupPValueMap.getOrDefault(i,new ArrayList<>());
                        List<JSONObject> outputTempList = groupOutputMap.getOrDefault(i,new ArrayList<>());
                        outputTempList.add(outputJson);
                        groupOutputMap.put(i,outputTempList);
                        pvalueList.add(pvalue);

                        groupPValueMap.put(i,pvalueList);
                    }
                    else
                    {
                        outputJson.put(qvalueHeader,"NA");
                        outputJson.put(pvalueHeader,"NA");
                    }
                }


            }

            if( idSet.size()==3)
            {
                double [] secondGroup = aucGroup.get(1).toNativeArray();
                double [] thirdGroup = aucGroup.get(2).toNativeArray();
                AnalyzeArrayPair pair = new AnalyzeArrayPair(secondGroup,thirdGroup);
                double pvalue = pair.getPvalue();
                double ratio = pair.getAverageRatio();


                String ratioHeader = "AVERAGE_RATIO_3_2";
                String pvalueHeader = "PVALUE_3_2";
                outputJson.put(ratioHeader,ratio);
                String qvalueHeader = "QVALUE_3_2";

                if(!Double.isNaN(pvalue) && !Double.isInfinite(pvalue))
                {
                    outputJson.put(pvalueHeader,pvalue);
                    List<Double> pvalueList = pvalue3_2;

                    pvalueList.add(pvalue);


                }
                else
                {
                    outputJson.put(qvalueHeader,"NA");
                    outputJson.put(pvalueHeader,"NA");
                }

            }
            outputJson.put("chro_iso_arr",chroArray);
            outputJson.put("peaks_arr",peakArray);
            outputJson.put("startRT_arr",startRTArray);
            outputJson.put("endRT_arr",endRTArray);
            outputJson.put("sampleName_arr",fileNameArray);

            outputList.add(outputJson);
        }

        for(Map.Entry<Integer,List<Double>> entry: groupPValueMap.entrySet())
        {
            List<Double> pvalues = entry.getValue();
            List<JSONObject> outputJsonList = groupOutputMap.get(entry.getKey());
            List<Double> qvalues = BHCorrection.runBhCorrection(pvalues);
            for(int i=0; i<outputJsonList.size(); i++)
            {
                JSONObject o = outputJsonList.get(i);
                int outi = entry.getKey()+1;
                String qvalueHeader = "QVALUE_"+outi+"_1";
                o.put(qvalueHeader,qvalues.get(i));
            }
        }

        if(pvalue3_2.size()>0)
        {
            List<Double> qvalues = BHCorrection.runBhCorrection(pvalue3_2);
            for(int i=0; i<outputList3_2.size(); i++)
            {
                JSONObject o = outputList3_2.get(i);
                int outi = i+1;
                String qvalueHeader = "QVALUE_3_2";
                o.put(qvalueHeader,qvalues.get(i));
            }
        }




        BufferedWriter bw = new BufferedWriter(new FileWriter(outputPath));
        bw.append("{\n");
        bw.append("\"peptideList\":[\n");
        for(int i=0; i<outputList.size(); i++)
        {
            JSONObject o = outputList.get(i);
            String jsonData = o.toJSONString();
            String formatedJSONData = JsonWriter.formatJson(jsonData);
            bw.append(formatedJSONData);
            bw.newLine();
            if(i!=outputList.size()-1)bw.append(",");
        }
        bw.append("]\n");
        bw.append("}");
        bw.close();



    }





    public static String createKey(int groupCode, int sampleCode, String fileName)
    {
        return ""+groupCode+sampleCode+fileName;
    }


    public static double generatePvalue(double [] group1, double group2 [], double dev1, double dev2, double avg1, double avg2) throws Exception {
        if(Double.compare(dev1,0)==0 && Double.compare(dev2,0)==0)
        {
            if((Double.compare(avg1,avg2)!=0))
            {
                return 1E-100;
            }
        }
        else
        {
            List classes = new ArrayList();
            classes.add(group1);
            classes.add(group2);
            double pvalue = AnovaUtil.calculateAnovaPvalue(classes);
            return pvalue;
        }
        return Double.NaN;
    }





    public static class AnalyzeArrayPair
    {
        private double [] group1;
        private double [] group2;
        private double averageRatio;
        private double pvalue;
        private double average1;
        private double average2;
        public AnalyzeArrayPair(double [] group1, double [] group2 ) throws Exception {
            this.group1 = group1;
            this.group2 = group2;
            init();
        }
        private void init() throws Exception {


            double total =0;
            for(double d: group1)
            {
                total+=d;
            }
            average1 = total/group1.length;
            StandardDeviation deviation = new StandardDeviation();
            double dev1 = deviation.evaluate(group1);

            double total2 =0;
            for(double d: group2)
            {
                total2+=d;
            }
            double dev2 = deviation.evaluate(group2);
            average2 = total2/group2.length;
            averageRatio = average2/average1;

            if(Double.compare(average1,0)==0)
                averageRatio =20;

            pvalue = generatePvalue(group1,group2,dev1,dev2,average1,average2);
        }

        public double getAverageRatio() {
            return averageRatio;
        }

        public double getPvalue() {
            return pvalue;
        }

        public double getAverage1() {
            return average1;
        }

        public double getAverage2() {
            return average2;
        }
    }



    public static class AnalysisPeptide
    {
        public final String sequence;
        public final String cs;
        public final double auc;
        public final int groupName;
        public final String sampleName;
        public final String fileName;
        private String peaks;
        private String startRt;
        private String endRt;
        private String chroIso;
        private List<String> scanList = new ArrayList<>();
        private List<String> deltaCNList = new ArrayList<>();
        private List<String> redundacy = new ArrayList<>();

        public AnalysisPeptide(String sequence, String cs, double auc, int groupName, String sampleName, String fileName) {
            this.sequence = sequence;
            this.cs = cs;
            this.auc = auc;
            this.groupName = groupName;
            this.sampleName = sampleName;
            this.fileName = fileName;
        }

        public String getPeaks() {
            return peaks;
        }

        public void setPeaks(String peaks) {
            this.peaks = peaks;
        }

        public String getStartRt() {
            return startRt;
        }

        public void setStartRt(String startRt) {
            this.startRt = startRt;
        }

        public String getEndRt() {
            return endRt;
        }

        public void setEndRt(String endRt) {
            this.endRt = endRt;
        }

        public String getChroIso() {
            return chroIso;
        }

        public void setChroIso(String chroIso) {
            this.chroIso = chroIso;
        }
    }



    public static double [] sumArrFromMSFile()
    {

        return null;
    }
}
