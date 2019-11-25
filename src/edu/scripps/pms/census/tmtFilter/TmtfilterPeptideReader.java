package edu.scripps.pms.census.tmtFilter;

import edu.scripps.pms.util.seq.Fasta;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * How to add new column....
 * 1) Create the index variable for that column......
 * 2) add the entry in the parseHeaderLine() method......
 * 3) add final entry in the getPeptide() method....
 * @author Harshil
 * @author 04/24/18 Titus Jung
 */
public class TmtfilterPeptideReader 
{
    private String filePath="";
    private String header = "";
    private int totalComparrison=0;
    private BufferedReader br = null;
    private String currentLine = null;
    private List<String> sampleList = new ArrayList<String>();
    private List<String> compareList = new ArrayList<String>();
    private String originalFile;
    private String runDate;
    
    private int sequenceIndex =-1;
    private int scanNoIndex =-1;
    private List<Integer> pValueIndex =new ArrayList<Integer>();
    private List<Integer> normPValueIndex =new ArrayList<Integer>();

    private List<Integer> pValueRatioIndex =new ArrayList<Integer>();
    private List<Integer> normPValueRatioIndex =new ArrayList<Integer>();

    private List<Integer> pValueRatioDataIndex =  new ArrayList<Integer>();
    private List<Integer> normPValueRatioDataIndex =  new ArrayList<Integer>();

    private List<Integer> BhCorrectionIndex=new ArrayList<Integer>();
    private List<Integer> normBhCorrectionIndex=new ArrayList<Integer>();

    private List<Integer> BhCorrectionRatioIndex=new ArrayList<Integer>();
    private List<Integer> normBhCorrectionRatioIndex=new ArrayList<Integer>();

    private int intensityAvgindex =-1;
    private int normIntensityAvgindex =-1;

    private List<Integer> ratioIndex =new ArrayList<Integer>();
    private List<Integer> normRatioIndex =new ArrayList<Integer>();

    private List<Integer> ratioListIndex=new ArrayList<Integer>();
    private List<Integer> normRatioListIndex=new ArrayList<Integer>();

    private List<Integer> logRatioIndex =new ArrayList<Integer>();
    private List<Integer> normLogRatioIndex =new ArrayList<Integer>();

    private List<Integer> logRatioSTDEVIndex =new ArrayList<Integer>();
    private List<Integer> normLogRatioSTDEVIndex =new ArrayList<Integer>();

    private List<Integer> logRatioRSDIndex =new ArrayList<Integer>();
    private List<Integer> normLogRatioRSDIndex =new ArrayList<Integer>();
    
    private List<Integer> stdevPaired =new ArrayList<Integer>();
    private List<Integer> rsdPairedIndex =new ArrayList<Integer>();

    private List<Integer> normStdevPaired =new ArrayList<Integer>();
    private List<Integer> normRsdPairedIndex =new ArrayList<Integer>();
    
    private int spcIndex =-1;
    private int outlierIndex =-1;
    private int fileNameIndex =-1;
    private int accessionIndex =-1;
    private int proteinDescriptionIndex =-1;
    private int geneNameIndex = -1;
    
    
    
    public TmtfilterPeptideReader(String filepath) {
        try {
            this.filePath = filepath;
            br = new BufferedReader(new FileReader(new File(filepath)));
        } catch (FileNotFoundException ex) {
            Logger.getLogger(TmtfilterPeptideReader.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    public static void main(String argsp[])
    {
//        TmtfilterPeptideReader reader = new TmtfilterPeptideReader("E:\\census-out-peptideAnalysisSimple.txt");
                //TmtfilterPeptideReader reader = new TmtfilterPeptideReader("C:\\Users\\Harshil\\Desktop\\quant_compare2600.txt");
       // TmtfilterPeptideReader reader = new TmtfilterPeptideReader("/home/yateslab/project_data/census/1804tmtDebug/census-out-peptideAnalysisSimple.txt");
        TmtfilterPeptideReader reader = new TmtfilterPeptideReader("/Users/rpark/quant_compare8.txt");
//                TmtfilterPeptideReader reader = new TmtfilterPeptideReader("/data/2/rpark/ip2_data//bcsmith/GSNOR_inhibitors/compare/quant/quant_compare2854.txt");
//        TmtfilterPeptideReader reader = new TmtfilterPeptideReader("/data/2/E:\\quant_compare2518.txtrpark/ip2_data//yrc/EvoProteomics/compare/quant/quant_compare2500.txt");
        List<Peptide> peptideLsit = reader.getPeptides();
        System.out.println("Data read success:....");
    }
    public List<Peptide> init()
    {
        
            readHeader();
            List<Peptide> peptideList = new ArrayList<Peptide>();
            while(currentLine!=null)
            {
                try 
                {
                    peptideList.add(getPeptide(currentLine));
                    currentLine=br.readLine();
                } 
                catch (IOException ex) 
                {
                    Logger.getLogger(TmtfilterPeptideReader.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
            return peptideList;

    }
    
    private void readHeader()
    {
        try {
            currentLine = br.readLine();
            this.runDate = currentLine.split("\t")[1];
            
            while(currentLine.startsWith("H"))
            {
                String words[] = currentLine.split("\t");
                if(words[1].toLowerCase().contains("compare"))
                    totalComparrison++;
                
                if(words[1].startsWith("Sample")) {
                    String sampleName="";
                    for(int i=2;i<words.length;i++)
                        sampleName += words[i]+ " ";
                        
                    this.sampleList.add(sampleName);
                } else if(words[1].startsWith("Compare")) {
                    this.compareList.add(words[2]);
                } else if (words[1].startsWith("Path")) {
                    this.originalFile = words[2];
                }
                
                
                if(currentLine.contains("pep_sequence"))
                {
                    parseHeaderLine();
                }
                
                
                    
                    
                currentLine = br.readLine();
            }
        } catch (IOException ex) {
            Logger.getLogger(TmtfilterPeptideReader.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    private void parseHeaderLine()
    {
        String words[] = currentLine.split("\t");
        
        for (int i = 0; i < words.length; i++)
        {
            if(words[i].equals("pep_sequence"))
            {
                header = currentLine;
                sequenceIndex = i;
            }
            else if(words[i].equals("scan_num"))
                scanNoIndex = i;
            else if(!words[i].toLowerCase().contains("norm"))
            {
                if(words[i].contains("pvalue") && !words[i].toLowerCase().contains("ratio"))
                    pValueIndex.add(i);
                else if(words[i].contains("pvalue") && words[i].toLowerCase().contains("ratio_"))
                    pValueRatioIndex.add(i);
                else if(words[i].contains("pvalue") && words[i].toLowerCase().contains("data"))
                    pValueRatioDataIndex.add(i);
                else if(words[i].contains("BHCorrection") && !words[i].toLowerCase().contains("ratio"))
                    BhCorrectionIndex.add(i);
                else if(words[i].contains("BHCorrection") && words[i].toLowerCase().contains("ratio"))
                    BhCorrectionRatioIndex.add(i);
                else if(words[i].contains("intensity_avg") && intensityAvgindex== -1)
                    intensityAvgindex = i;
                else if(words[i].contains("ratio_avg") && !words[i].contains("log")  && !words[i].contains("list") )
                    ratioIndex.add(i);
                else if(words[i].contains("ratio_list") )
                    ratioListIndex.add(i);
                else if(words[i].contains("log_ratio_avg"))
                    logRatioIndex.add(i);
                else if(words[i].contains("log_ratio_STDEV") && !words[i].contains("paired"))
                    logRatioSTDEVIndex.add(i);
                else if(words[i].contains("log_ratio_RSD") && !words[i].contains("paired"))
                    logRatioRSDIndex.add(i);
                else if(words[i].contains("STDEV_paired")  )
                    stdevPaired.add(i);
                else if(words[i].contains("RSD_paired") )
                    rsdPairedIndex.add(i);
                else if(words[i].equals("SpC"))
                    spcIndex = i;
                else if(words[i].equals("outlier_count"))
                    outlierIndex = i;
                else if(words[i].equals("fileName"))
                    fileNameIndex = i;
                else if(words[i].equals("accession"))
                    accessionIndex = i;
                else if(words[i].equals("proteinDescription"))
                    proteinDescriptionIndex = i;
                else if(words[i].equals("Gene_Name"))
                {
                    geneNameIndex = i;
                }
            }
            else
            {
                if(words[i].contains("pvalue") && !words[i].toLowerCase().contains("ratio"))
                    normPValueIndex.add(i);
                else if(words[i].contains("pvalue") && words[i].toLowerCase().contains("ratio_"))
                    normPValueRatioIndex.add(i);
                else if(words[i].contains("pvalue") && words[i].toLowerCase().contains("data"))
                    normPValueRatioDataIndex.add(i);
                else if(words[i].contains("BHCorrection") && !words[i].toLowerCase().contains("ratio"))
                    normBhCorrectionIndex.add(i);
                else if(words[i].contains("BHCorrection") && words[i].toLowerCase().contains("ratio"))
                    normBhCorrectionRatioIndex.add(i);
                else if(words[i].contains("intensity_avg") && normIntensityAvgindex== -1)
                    normIntensityAvgindex = i;
                else if(words[i].contains("ratio_avg") && !words[i].contains("log")  && !words[i].contains("list") )
                    normRatioIndex.add(i);
                else if(words[i].contains("ratio_list") )
                    normRatioListIndex.add(i);
                else if(words[i].contains("log_ratio_avg"))
                    normLogRatioIndex.add(i);
                else if(words[i].contains("log_ratio_STDEV") && !words[i].contains("paired"))
                    normLogRatioSTDEVIndex.add(i);
                else if(words[i].contains("log_ratio_RSD") && !words[i].contains("paired"))
                    normLogRatioRSDIndex.add(i);
                else if(words[i].contains("STDEV_paired")  )
                    normStdevPaired.add(i);
                else if(words[i].contains("RSD_paired"))
                    normRsdPairedIndex.add(i);
            }

            
        }
    }
  
    private Peptide getPeptide(String currentLine)
    {
        Peptide peptide = new Peptide();
        String words[]=currentLine.split("\t");
        if(sequenceIndex != -1)
                peptide.setSequence(words[sequenceIndex]);
        if(scanNoIndex != -1)
        {
            String sacnList[] = words[scanNoIndex].split("[\\[\\], ]");
            for (String scanNumber : sacnList) {
                if (scanNumber.length() > 0) {
                    peptide.addScannum(Integer.parseInt(scanNumber));
                }
            }
        }
        
        if(pValueIndex.size() != 0)
        {
           for (int i = 0; i < pValueIndex.size(); i++) 
                peptide.addPValueList(Double.parseDouble(words[pValueIndex.get(i)]));
        }
        
        if(pValueRatioIndex.size() != 0)
        {
           for (int i = 0; i < pValueRatioIndex.size(); i++) 
                peptide.addPValueRatio(Double.parseDouble(words[pValueRatioIndex.get(i)]));
        }
        
        if(pValueRatioDataIndex.size() != 0)
        {
            for(int i =0 ;i < pValueRatioDataIndex.size();i++)
            {
                List<Double> data = new ArrayList<>();
                for(String value : words[pValueRatioDataIndex.get(i)].split(","))
                   data.add(Double.parseDouble(value));
                peptide.addRatioPvalueData(data);
            }
        }
        
        
        if(BhCorrectionIndex.size() != 0)
        {
            for (int i = 0; i < BhCorrectionIndex.size(); i++)
                peptide.addBHCorrection(Double.parseDouble(words[BhCorrectionIndex.get(i)]));
        }
        
        if(BhCorrectionRatioIndex.size() != 0)
        {
            for (int i = 0; i < BhCorrectionRatioIndex.size(); i++)
                peptide.addBHCorrectionRatio(Double.parseDouble(words[BhCorrectionRatioIndex.get(i)]));
        }
        
        if(intensityAvgindex != -1)
        {
            List<List<Double>> intensityAvgList = new ArrayList<>();

            for (int i = 0; i < sampleList.size(); i++) {

                List<Double> currentIntensityAvg = new ArrayList<>();
                for (String value : words[intensityAvgindex+i].split(",")) {
                    currentIntensityAvg.add(Double.parseDouble(value));
                }
                intensityAvgList.add(currentIntensityAvg);
            }
            peptide.setIntensityAverage(intensityAvgList);
        }

        if(ratioIndex.size() != 0)
        {
            for(int index : ratioIndex)
            {

//System.out.println("===" + index);
                if(words[index].contains(",")) {
                    String tmpStr = words[index].substring(0,words[index].indexOf(","));

                    peptide.addRatioAvg(Double.parseDouble(tmpStr));

                }
                else {

                    if(!"NA".equals(words[index]))
                        peptide.addRatioAvg(Double.parseDouble(words[index]));

                }

            }
        }
        if(ratioListIndex.size() != 0)
        {
            for(int index : ratioListIndex)
            {
                String ratioData[] = words[index].split(",");
                List ratio = new ArrayList();
                for (String currentRatio : ratioData) {
                    ratio.add(Double.parseDouble(currentRatio));
                }
                peptide.addRatio(ratio);
            }
        }
        if(logRatioIndex.size() != 0)
        {
            for(int index : logRatioIndex)
            {
                if(words[index].contains(",")) {
                    String tmpStr = words[index].substring(0,words[index].indexOf(","));

                    peptide.addLogRatioAvg(Double.parseDouble(tmpStr));
		    } else

                if(!"NA".equals(words[index]))
			        peptide.addLogRatioAvg(Double.parseDouble(words[index]));
            }
        }
        if(logRatioSTDEVIndex.size() != 0)
        {
            for(int index : logRatioSTDEVIndex)
            {
                if(words[index].contains(",")) {
                    String tmpStr = words[index].substring(0, words[index].indexOf(","));

                    peptide.addLogRatioStdDev(Double.parseDouble(tmpStr));

                } else
                    peptide.addLogRatioStdDev(Double.parseDouble(words[index]));
            }
        }
        if(logRatioRSDIndex.size() != 0)
        {
            for(int index : logRatioRSDIndex)
            {
                if(words[index].equalsIgnoreCase("NA"))
                    peptide.addLogRatioRSD(-1.0);
                else {

                    if(words[index].contains(",")) {
                        String tmpStr = words[index].substring(0, words[index].indexOf(","));

                        peptide.addLogRatioRSD(Double.parseDouble(tmpStr));

                    } else
                        peptide.addLogRatioRSD(Double.parseDouble(words[index]));

                }
            }
        }
        
        if(rsdPairedIndex.size() != 0)
        {
            for(int index : rsdPairedIndex)
            {
                          //      System.out.println(words[index]+"\t" +peptide.getSequence());

                if(words[index].equalsIgnoreCase("NA"))
                     peptide.addRsdPaired(-1.0);
                else {


                    if(words[index].contains(",")) {
                        String tmpStr = words[index].substring(0, words[index].indexOf(","));
                        peptide.addRsdPaired(Double.parseDouble(tmpStr));

                    } else
                        peptide.addRsdPaired(Double.parseDouble(words[index]));
                }
            }
        }
        
        if(stdevPaired.size() != 0)
        {
            for(int index : stdevPaired)
            {
                if(words[index].equalsIgnoreCase("NA"))
                    peptide.addstdevPaired(-1.0);
                else {

                    if(words[index].contains(",")) {
                        String tmpStr = words[index].substring(0, words[index].indexOf(","));
                        peptide.addstdevPaired(Double.parseDouble(tmpStr));

                    } else
                        peptide.addstdevPaired(Double.parseDouble(words[index]));

                }
            }
        }
        if(spcIndex != -1) {

            if(words[spcIndex].contains(",")) {
                String tmpStr = words[spcIndex].substring(0, words[spcIndex].indexOf(","));
                peptide.setSpc(Integer.parseInt(tmpStr));

            } else
                ;// peptide.setSpc(Integer.parseInt(words[spcIndex]));
        }

        
        if(outlierIndex != -1)
            peptide.setOutlier(Integer.parseInt(words[outlierIndex]));
        
        if(fileNameIndex != -1)
        {
            String fileList[] = words[fileNameIndex].split("\\[\\],");
            for (String fileName : fileList) {
                if (fileName.length() > 0) {
                    peptide.addFilename(fileName);
                }
            }
        }
        
         if(accessionIndex != -1)
        {
            String proteinLocus[] = words[accessionIndex].split("\\[\\],");
            StringBuilder sb = new StringBuilder();
            for (String locus : proteinLocus) {
                if (locus.length() > 0) {
                    String x = locus.replaceAll("[\\[\\]]","");
                    peptide.addProteinLocus(x);
                    sb.append(x).append(",");
                }
            }
            peptide.setProteinLocusListStr(sb.toString());
        }
        if(geneNameIndex != -1)
        {
            String geneNameArr[] = words[geneNameIndex].split("\\[\\],");
            StringBuilder sb = new StringBuilder();

            for (String gene : geneNameArr) {
                if (gene.length() > 0) {
                    String x = gene.replaceAll("[\\[\\]]","");
                    peptide.addGeneName(x);
                    sb.append(x).append(",");
                }
            }
            peptide.setGeneListStr(sb.toString());
        }

         if(proteinDescriptionIndex != -1)
        {
            String proteinDesc[] = words[proteinDescriptionIndex].split("\\[\\],");
            StringBuilder sb = new StringBuilder();
            for (String proteinDescription : proteinDesc) {
                if (proteinDescription.length() > 0) {
                    String x = proteinDescription.replaceAll("[\\[\\]]","");
                    peptide.addProteinDescription(x);
                    sb.append(x).append(",");
                }
            }
            peptide.setProteinDescriptionListStr(sb.toString());
        }

        if(normPValueIndex.size() != 0)
        {
            for (int i = 0; i < normPValueIndex.size(); i++)
                peptide.addNormPValueList(Double.parseDouble(words[normPValueIndex.get(i)]));
        }

        if(normPValueRatioIndex.size() != 0)
        {
            for (int i = 0; i < normPValueRatioIndex.size(); i++)
                peptide.addPValueNormRatio(Double.parseDouble(words[normPValueRatioIndex.get(i)]));
        }

        if(normPValueRatioDataIndex.size() != 0)
        {
            for(int i =0 ;i < normPValueRatioDataIndex.size();i++)
            {
                List<Double> data = new ArrayList<>();
                for(String value : words[normPValueRatioDataIndex.get(i)].split(","))
                    data.add(Double.parseDouble(value));
                peptide.addRationNormPvalueData(data);
            }
        }


        if(normBhCorrectionIndex.size() != 0)
        {
            for (int i = 0; i < normBhCorrectionIndex.size(); i++)
                peptide.addBHCorrectionNorm(Double.parseDouble(words[normBhCorrectionIndex.get(i)]));
        }

        if(normBhCorrectionRatioIndex.size() != 0)
        {
            for (int i = 0; i < normBhCorrectionRatioIndex.size(); i++)
                peptide.addBHCorrectionNormRatio(Double.parseDouble(words[normBhCorrectionRatioIndex.get(i)]));
        }

        if(normIntensityAvgindex != -1)
        {
            List<List<Double>> intensityAvgList = new ArrayList<>();

            for (int i = 0; i < sampleList.size(); i++) {

                List<Double> currentIntensityAvg = new ArrayList<>();
                for (String value : words[normIntensityAvgindex+i].split(",")) {
                    currentIntensityAvg.add(Double.parseDouble(value));
                }
                intensityAvgList.add(currentIntensityAvg);
            }
            peptide.setNormIntensitySum(intensityAvgList);
        }

        if(normRatioIndex.size() != 0)
        {
            for(int index : normRatioIndex)
            {

//System.out.println("===" + index);
                if(!"NA".equals(words[index]))
                    peptide.addNormRatioAvg(Double.parseDouble(words[index]));
                else
                    peptide.addNormRatioAvg(0);

            }
        }
        if(normRatioListIndex.size() != 0)
        {
            for(int index : normRatioListIndex)
            {
                String ratioData[] = words[index].split(",");
                List ratio = new ArrayList();
                for (String currentRatio : ratioData) {
                    ratio.add(Double.parseDouble(currentRatio));
                }
                peptide.addNormRatio(ratio);
            }
        }
        if(normLogRatioIndex.size() != 0)
        {
            for(int index : normLogRatioIndex)
            {
                if(!"NA".equals(words[index]))
                    peptide.addLogRatioNormAvg(Double.parseDouble(words[index]));
                else
                    peptide.addLogRatioNormAvg(0);
            }
        }
        if(normLogRatioSTDEVIndex.size() != 0)
        {
            for(int index : normLogRatioSTDEVIndex)
            {
                if(!"NA".equals(words[index]))
                    peptide.addNormLogRatioStdDev(Double.parseDouble(words[index]));
                else
                    peptide.addNormLogRatioStdDev(0);
            }
        }
        if(normLogRatioRSDIndex.size() != 0)
        {
            for(int index : normLogRatioRSDIndex)
            {
                if(words[index].equalsIgnoreCase("NA"))
                    peptide.addLogRatioNormRSD(-1.0);
                else
                    peptide.addLogRatioNormRSD(Double.parseDouble(words[index]));
            }
        }

        if(normRsdPairedIndex.size() != 0)
        {
            for(int index : normRsdPairedIndex)
            {
            //    System.out.println(words[index]+"\t" +peptide.getSequence());

                if(words[index].equalsIgnoreCase("NA"))
                    peptide.addNormRsdPaired(-1.0);
                else
                    peptide.addNormRsdPaired(Double.parseDouble(words[index]));
            }
        }

        if(normStdevPaired.size() != 0)
        {
            for(int index : normStdevPaired)
            {
                if(words[index].equalsIgnoreCase("NA"))
                    peptide.addNormstdevPaired(-1.0);
                else
                    peptide.addNormstdevPaired(Double.parseDouble(words[index]));
            }
        }

        
        return peptide;
    }
    public List<Peptide> getPeptides()

    {
        return init();
    }

    public List<String> getSampleList() {
        return sampleList;
    }

    public void setSampleList(List<String> sampleList) {
        this.sampleList = sampleList;
    }

    public List<String> getCompareList() {
        return compareList;
    }

    public void setCompareList(List<String> compareList) {
        this.compareList = compareList;
    }

    public String getOriginalFile() {
        return originalFile;
    }

    public void setOriginalFile(String originalFile) {
        this.originalFile = originalFile;
    }

    public String getRunDate() {
        return runDate;
    }

    public void setRunDate(String runDate) {
        this.runDate = runDate;
    }

    public static String extractGeneName (String descriptiveName) {  //single protein
        String gene="";
        String[] arr = descriptiveName.split("GN=");
        if(arr.length<=1) return gene;

        String[] arr2 = arr[1].split(" ");

        gene = arr2[0];
        return gene;

    }


}
