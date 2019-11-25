/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.scripps.pms.census.labelFree;

import edu.scripps.pms.census.conf.Configuration;
import edu.scripps.pms.census.hash.IndexedFile;
import edu.scripps.pms.census.model.ChroPeptide;
import edu.scripps.pms.census.model.ChroProtein;
import edu.scripps.pms.census.model.SampleModel;
import edu.scripps.pms.census.tools.IsotopeModel;
import edu.scripps.pms.census.util.RelExFileFilter;
import edu.scripps.pms.util.spectrum.Range;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.jdom.Document;
import org.jdom.Element;
import org.jdom.JDOMException;
import org.jdom.input.SAXBuilder;

/**
 *
 * @author Harshil
 * @version $Id:
 */
public class ChroJSONReader {

    protected static final String NORMINTENSITY = "NORMINTENSITY";
    protected static final String IONINJECTIONTIME = "IONINJECTIONTIME";
    protected static final String RETENTIONTIME = "RETENTIONTIME";
    protected static final String ENDRANGE = "ENDRANGE";
    protected static final String STARTRANGE = "STARTRANGE";
    protected static final String REDUNDANCY = "REDUNDANCY";
    protected static final String SPSCORE = "SPSCORE";
    protected static final String SPRANK = "SPRANK";
    protected static final String DMASS = "DMASS";
    protected static final String DCN = "DCN";
    protected static final String XCORR = "XCORR";
    protected static final String TOTALINTENSITY = "TOTALINTENSITY";
    protected static final String CALCMHPLUS = "CALCMHPLUS";
    protected static final String MHPLUS = "MHPLUS";
    protected static final String PROFILE_SCORE = "PROFILE_SCORE";
    protected static final String INTENSITY = "INTENSITY";
    protected static final String CSTATE = "CSTATE";
    protected static final String SCAN = "SCAN";
    protected static final String SEQUENCE = "SEQUENCE";
    protected static final String FILENAME = "FILENAME";
    protected static final String MISSED = "Missed";

    private String configFile = null;
    private List<SampleModel> sampleList =new ArrayList<>();
    private HashMap<Integer ,List<Integer>> groupNamToIndex = new HashMap<>();
    public ChroJSONReader(String configFile)
    {
        this.configFile = configFile;
        Configuration conf = Configuration.getInstance();
        if(!conf.isReadConfigFile())
            try {
                conf.readXMLParam(configFile);
        } catch (Exception ex) {
            Logger.getLogger(ChroJSONReader.class.getName()).log(Level.SEVERE, null, ex);
        }
        Element root = conf.getRootConfEle();
        List<Element> groupList = root.getChildren("sample");
        for(Element currentGroup : groupList)
        {
            String name = currentGroup.getAttributeValue("group");
            SampleModel sModel = new SampleModel(name);
            List<Element> fileList = currentGroup.getChildren("each_sample");
            for(Element eachFile : fileList)
            {
                sModel.addPath(eachFile.getChild("ms_files").getChild("file").getValue());
            }
            sampleList.add(sModel);
            
        }
        
        
        int counter =0;
        int groupCounter =0;
        List<Integer> tempList = null;
        for(SampleModel sModel  : sampleList)
        {
            for(String path :sModel.getPathList())
            {
                if(groupNamToIndex.containsKey(groupCounter)){
                    
                    tempList  = groupNamToIndex.get(groupCounter);
                    tempList.add(counter);
                }
                else{
                    tempList = new ArrayList<>();
                    tempList.add(counter);
                }
                groupNamToIndex.put(groupCounter, tempList);
                counter++;
            }
            groupCounter++;
        }
    }
    
    public void write(ChroJSONProteinModel jsonModel, String outputFile,int counter) {
        StringBuffer sb = new StringBuffer();
        if(counter == 0)
        {
            int totalExperiment = jsonModel.getPeptideGroup().get(0).size();
            sb.append(createHeader(totalExperiment));
        }
        //read scan tag from json file,... nad pass ot as key to LabelfreeCalcUtil.. do the rest..args.
        ChroProtein protein = jsonModel.getProtein();
        double total = 0;
        sb.append(toProteinSting(protein,jsonModel.getMedianIntensity()));
        List<ChroProtein> redundantProteinList = jsonModel.getRedundantProtein();
        for(ChroProtein redundantProtein : redundantProteinList)
        {
            sb.append(toProteinSting(redundantProtein,jsonModel.getMedianIntensity()));
        }
        for (List<ChroPeptide> peptideList : jsonModel.getPeptideGroup()) {
            sb.append("S");
            for (int i = 0; i < peptideList.size(); i++) {
                ChroPeptide peptide = peptideList.get(i);
                sb.append("\t").append(toPeptideString(peptide, i));

            }
            sb.append("\n");
        }
        System.out.print("peptide Counter : " + counter +"\r");
        BufferedWriter bw = null;

//        String outputFile = outputPath + File.separator + "output.txt";
        try {
            if(counter==0)
                bw = new BufferedWriter(new FileWriter(outputFile));
            else
                bw= new BufferedWriter(new FileWriter(outputFile, true));
            bw.write(sb.toString());
        } catch (Exception e) {
            Logger.getLogger(ChroJSONReader.class.getName()).log(Level.SEVERE, null, e);
        } finally {
            try {
                bw.close();
                
            } catch (IOException ex) {
                Logger.getLogger(ChroJSONReader.class.getName()).log(Level.SEVERE, null, ex);
            }

        }
    }

    private String createHeader(int totalExperiments) {
        StringBuffer headerBuffer = new StringBuffer();
        for(SampleModel sModel : sampleList)
        {
            for(String filePath : sModel.getPathList())
            {
                headerBuffer.append("H\t").append(sModel.getSampleName()).append("\t");
                headerBuffer.append(filePath).append("\n");
            }
        }
        
        
        headerBuffer.append("PLINE\tACCESSION\tDESCRIPTION");
//        int totalExperiments = 4;
        for (int i = 0; i < totalExperiments; i++) 
          headerBuffer.append("\tAvgNormIntensity_" + (i + 1));

        headerBuffer.append("\n");

        headerBuffer.append("SLINE");
        for (int i = 0; i < totalExperiments; i++) {
            headerBuffer.append("\tEXP_" + (i + 1));
            headerBuffer.append("\t"+SEQUENCE);
            headerBuffer.append("\t"+FILENAME);
            headerBuffer.append("\t"+SCAN);
            headerBuffer.append("\t"+CSTATE);
            
            headerBuffer.append("\t"+INTENSITY);
            headerBuffer.append("\t"+PROFILE_SCORE);
            headerBuffer.append("\t"+MHPLUS);
            headerBuffer.append("\t"+CALCMHPLUS);
            
//            headerBuffer.append("\t"+TOTALINTENSITY);
            headerBuffer.append("\t"+XCORR);
            headerBuffer.append("\t"+DCN);
            headerBuffer.append("\t"+DMASS);
            
            headerBuffer.append("\t"+SPRANK);
            headerBuffer.append("\t"+SPSCORE);
            headerBuffer.append("\t"+REDUNDANCY);
            headerBuffer.append("\t"+STARTRANGE);
            
            headerBuffer.append("\t"+ENDRANGE);
            headerBuffer.append("\t"+RETENTIONTIME);
            headerBuffer.append("\t"+IONINJECTIONTIME);
            headerBuffer.append("\t"+NORMINTENSITY);
            headerBuffer.append("\t"+MISSED);
            
        }
        return headerBuffer.append("\n").toString();
    }
    
    

    private String toPeptideString(ChroPeptide chroPeptide, int index) {
        StringBuffer sb = new StringBuffer();
        sb.append("[").append(index + 1).append("]");
        sb.append("\t").append(chroPeptide.getSequence());
        sb.append("\t").append(chroPeptide.getFileName());
        sb.append("\t").append(chroPeptide.getScanNum());
        sb.append("\t").append(chroPeptide.getChargeState());
        
        sb.append("\t").append(chroPeptide.getAverageIntensity());
        sb.append("\t").append(chroPeptide.getAnCompositeScore());
        sb.append("\t").append(chroPeptide.getMhPlus());
        sb.append("\t").append(chroPeptide.getCalcMHplus());
//        sb.append("\t").append(chroPeptide.getTotalIntensity());
        
        sb.append("\t").append(chroPeptide.getXCorr());
        sb.append("\t").append(chroPeptide.getDeltCN());
        sb.append("\t").append(chroPeptide.getDeltMass());
        sb.append("\t").append(chroPeptide.getSpRank());
        
        sb.append("\t").append(chroPeptide.getSpScore());
        sb.append("\t").append(chroPeptide.getRedundancy());
        sb.append("\t").append(chroPeptide.getStartRange());
        sb.append("\t").append(chroPeptide.getEndRange());
        
        sb.append("\t").append(chroPeptide.getRetentionTime());
        sb.append("\t").append(chroPeptide.getIonInjectionTime());
        sb.append("\t").append(chroPeptide.getIonInjectionTimeNormIntensity());
        sb.append("\t").append(chroPeptide.isMissedPeptide());
               
        return sb.toString();
    }

    private String toProteinSting(ChroProtein chroProtein,double meidanIntensity) {
        //    PLINE   ACCESSION       DESCRIPTION     SCOUNT_1        SCOUNT_2        SCOUNT_3        SCOUNT_4        PEP_COUNT       NORM_INTENSITY_1        NORM_INTENSITY_2        NORM_INTENSITY_3        NORM_INTENSITY_4        SAMPLE_GROUP_INTENSITY_AVG_1    SAMPLE_GROUP_INTENSITY_AVG_2    NORM_INTENSITY_CORRECT_1        NORM_INTENSITY_CORRECT_2        NORM_INTENSITY_CORRECT_3        NORM_INTENSITY_CORRECT_4        SAMPLE_GROUP_INTENSITY_CORRECT_AVG_1    SAMPLE_GROUP_INTENSITY_CORRECT_AVG_2    ^M
        StringBuffer sb = new StringBuffer();
        sb.append("P").append("\t");
        sb.append(chroProtein.getLocus()).append("\t");
        sb.append(chroProtein.getDescription());
//        for(Double avgIntensity : chroProtein.getAvgIntensityList())
//        {
//            sb.append("\t").append(avgIntensity);
//
//        }
        sb.append("\t").append(meidanIntensity);
        
        sb.append("\n");
        
        return sb.toString();

    }

    /* 
     * reads the config file gets the ms file list ,then reads the json file and 
     * generate the txt file with normalized intensity.
     */
    public void parse(String jsonFilePath) {

        
        Hashtable<String, IndexedFile> ht = getmsFiles(configFile);
        File file = null;
        IndexedFile iFile = null;
        File f = new File(jsonFilePath);
        String projectId = configFile.split("census_config_labelfree_")[1].split(".xml")[0];
        HashSet<String> set = peptideQualitySet(f.getParent()+ File.separator +"census_labelfree_out_"+projectId+".txttmp" );
        String outputFile = f.getParent() + File.separator + "census_labelfree_result_"+ projectId+".txt";
        HashMap<String, IndexModel> proteinToLocationMap = JsonIndexer.readJsonIndexFile(jsonFilePath + ".index");
        int counter =-1;
        for (IndexModel indexOfProtein : proteinToLocationMap.values()) {
            counter++;
            //read scan tag from json file,... nad pass ot as key to LabelfreeCalcUtil.. do the rest..args.
            ChroJSONProteinModel jsonModel = JsonIndexer.getProteinObject(jsonFilePath, indexOfProtein);
            
            generateNormalizedIntensity(jsonModel);
            generateRatioCombination(jsonModel);
            generateMeidianIntensity(jsonModel,set);
            write(jsonModel, outputFile,counter);
            
        }
        System.out.println("Check the file at " + outputFile);
    }
    
    private HashSet<String> peptideQualitySet(String path)
    {
        BufferedReader br = null;
        HashSet<String> set = new HashSet<String>();
        try {
            br = new BufferedReader(new FileReader(path));
            //      BufferedReader br = new BufferedReader(new FileReader("a.txt"));
            String eachLine;
            
            while ((eachLine = br.readLine()) != null) {
                if (!eachLine.startsWith("S\t")) {
                    continue;
                }
                
                String[] arr = eachLine.split("\t");
                
              //  System.out.println("===" + eachLine);
              //  System.out.println("=----------" + arr[7] + " " + arr[27] +  " " );
              
		//robin come and fix this code to remove hard value
                int count = (arr.length-1)/20;
              
                //if ("NA".equals(arr[7]) || "NA".equals(arr[27]) || "NA".equals(arr[47])) {
                //    continue;
               // }
               //
               
               
               boolean validPeptide=true;
               
               for(int i=1;i<=count;i++) {
                   if("NA".equals(arr[i*7]) || Double.parseDouble(arr[7])<0.5)
                       validPeptide = false;

	       }
                
               if(!validPeptide) continue;
                /*
                double p1 = Double.parseDouble(arr[7]);
                double p2 = Double.parseDouble(arr[27]);
                double p3 = Double.parseDouble(arr[47]);
                
                if (p1 < 0.7 || p2 < 0.7 || p3 < 0.7) {
                    continue;
                }
                */
                
                set.add(arr[2]);
                //System.out.println(arr[2]);
                
            }
        } catch (FileNotFoundException ex) {
            Logger.getLogger(ChroJSONReader.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(ChroJSONReader.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            try {
                br.close();
            } catch (IOException ex) {
                Logger.getLogger(ChroJSONReader.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        
        return set;


    }
    private void generateMeidianIntensity(ChroJSONProteinModel jsonModel,HashSet<String> set)
    {
        DescriptiveStatistics stat = new DescriptiveStatistics();
        List<List<ChroPeptide>> peptideGroup = jsonModel.getPeptideGroup();
        for(List<ChroPeptide> currentPeptide : peptideGroup)
        {
            if(set.contains(currentPeptide.get(0).getSequence() ))
                stat.addValue(currentPeptide.get(0).getTotalIntensity()/currentPeptide.get(1).getTotalIntensity());

        }
        //if(stat.getValues().length != 0)
            
        jsonModel.setMedianIntensity(stat.getPercentile(50));
        
    }
    /**
     * Peptide level does 1-1 group ratio of avg intensity.
     * 1/0 -- 2/0 -- 3/0....etc
     * @param jsonModel 
     */
    private void generateRatioCombination(ChroJSONProteinModel jsonModel)
    {
        List<List<ChroPeptide>> peptideGroup = jsonModel.getPeptideGroup();
        
        for(List<ChroPeptide> currentPeptide : peptideGroup)
        {
            List<Double> ration = new ArrayList<>();
            for(int i=1;i<currentPeptide.size();i++)
            {
                ration.add(currentPeptide.get(i).getTotalIntensity()/currentPeptide.get(0).getTotalIntensity());
            }
            jsonModel.addPeptideRatioIntensity(ration);
        }
    }
    /**
     * Generates Normalized Intensity and AverageIntensioty for proteinlevel
     * @param jsonModel 
     */
    private void generateNormalizedIntensity(ChroJSONProteinModel jsonModel )
    {
        double total = 0;
        if(jsonModel.getPeptideGroup().size()==0){
            System.out.println("no peptides");
            return;}
        int totalExperimentCount=jsonModel.getPeptideGroup().get(0).size();
        int totalPeptideCount=jsonModel.getPeptideGroup().size();
        
        List<DescriptiveStatistics> statList = new ArrayList<>();
        
            
            for (List<ChroPeptide> peptideList : jsonModel.getPeptideGroup()) {
                for (int j=0;j<peptideList.size();j++) {
                    statList.add(new DescriptiveStatistics());
                }
                for (int j=0;j<peptideList.size();j++) {
//                for(  ChroPeptide peptide : peptideList){
                    ChroPeptide peptide = peptideList.get(j);
                    
                    try {
                        
                        Range range = new Range(Double.parseDouble(peptide.getStartRange()), Double.parseDouble(peptide.getEndRange()));
                        for (int i = 0; i < peptide.getIsoTopeModelList().size(); i++) {
                            double ionInjectionScaleFactor = 0;
                            double retentionScalingFactor = 0;

                            IsotopeModel isoTope = peptide.getIsoTopeModelList().get(i);

                            if (range.isInRange(isoTope.getScanNumber())) {
                                for (double value : isoTope.getIntensityArr()) {
                                    ionInjectionScaleFactor += value;
                                }
                            }
                            if (isoTope.getIonInjectionTime() != 0) {
                                ionInjectionScaleFactor /= isoTope.getIonInjectionTime();
                            }
                            if (i != 0) {
                                retentionScalingFactor = (peptide.getIsoTopeModelList().get(i).getRetentionTime() - peptide.getIsoTopeModelList().get(i - 1).getRetentionTime()) / 2;
                            } else {
                                retentionScalingFactor = peptide.getIsoTopeModelList().get(i).getRetentionTime();
                            }
                            if (i == peptide.getIsoTopeModelList().size() - 1) {
                                retentionScalingFactor += peptide.getIsoTopeModelList().get(i).getRetentionTime();
                            } else {
                                retentionScalingFactor += (peptide.getIsoTopeModelList().get(i + 1).getRetentionTime() - peptide.getIsoTopeModelList().get(i).getRetentionTime()) / 2;
                            }

                            total += ionInjectionScaleFactor * retentionScalingFactor;

                        }
                        peptide.setIonInjectionTimeNormIntensity(total);
                        total=0;
//                            System.out.println("--***--" + peptide.getSequence() + "\t" + total);

//                            we will have the scanNumber range .... 
//                            Now traverse it find the area=>sum of intensity in that range,....
//                            print the values to the output file.txt
                    } catch (Exception ex) {
                        Logger.getLogger(ChroJSONReader.class.getName()).log(Level.SEVERE, null, ex);
                    }
                    
                    statList.get(j).addValue(peptide.getIonInjectionTimeNormIntensity());
                }
                
            }
            ChroProtein protein = jsonModel.getProtein();
            for (DescriptiveStatistics stat : statList){
                protein.addAvgIntensityList(stat.getMean());
            }
    }
    private Hashtable<String, IndexedFile> getmsFiles(String configFile) {
        //this will get the MsFilePath and name by combining the config file path and global msfilesName...
        Hashtable<String, IndexedFile> ht = new Hashtable<>();
        try {
            Document doc = new SAXBuilder().build(new File(configFile));
            Element root = doc.getRootElement();
            List<Element> sampleList = root.getChildren("sample");
            for (Element sample : sampleList) {
                List<Element> eachSampleList = sample.getChildren("each_sample");
                for (Element eachSample : eachSampleList) {
                    List<Element> eachmsFilesList = eachSample.getChildren("ms_files");
                    for (Element eachmsFile : eachmsFilesList) {
                        String path = eachmsFile.getChildText("file").split("\\*.")[0];

                        File f = new File(path);
                        String[] list = f.list(new RelExFileFilter("ms1"));
                        for (String value : list) {
                            File indexFile = new File(path + File.separator + value + ".index");
//                            ht.put(value.split(".ms1")[0], new IndexedFile(indexFile, path+File.separator+value));
                            ht.put(value, new IndexedFile(indexFile, path + File.separator + value));
                        }
                    }
                }
            }

        } catch (JDOMException ex) {
            Logger.getLogger(LabelFreeParser.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(LabelFreeParser.class.getName()).log(Level.SEVERE, null, ex);
        }
        return ht;
    }

    private static void printUsage() {
        System.out.println("Usage: ChroJSONReader configFile jsonFile");;
    }

    public static void main(String args[]) {
        if (args.length != 2) {
            printUsage();
            return;
        }

        String configFileName = args[0];
        String jsonFile = args[1];
        ChroJSONReader cr = new ChroJSONReader(configFileName);
        cr.parse(jsonFile);

    }
}
