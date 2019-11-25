/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package edu.scripps.pms.census.tmtFilter;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * How to add new column....
 * 1) Create the index variable for that column......
 * 2) add the entry in the parseHeaderLine() method......
 * 3) add final entry in the getProtein() method....
 * @author Harshil
 */
public class TmtfilterProteinReader 
{
    private String filePath;
    private BufferedReader br= null;
    private String currentLine;
    private int totalComparrison=0;
    private List<String> sampleList = new ArrayList<String>();
    private List<String> compareList = new ArrayList<String>();
    private String originalFile;
    private String runDate;
    
    private int locusIndex =-1;
    private int specCountIndex =-1;
    private int pepNumIndex =-1;
    private List<Integer> pValueIndex =new ArrayList<Integer>();
    private List<Integer> bHCorrectionIndex =new ArrayList<Integer>();
    private List<Integer> pValueRatioIndex =new ArrayList<Integer>();
    private List<Integer> pValueRatioDataIndex =new ArrayList<Integer>();
    private List<Integer> bHCorrectionRatioIndex =new ArrayList<Integer>();
    private int intensityAvgIndex =-1;
    private List<Integer> ratioAvg = new ArrayList<Integer>();
    private List<Integer> ratioSumIndex = new ArrayList<Integer>();
    private List<Integer> ratioStdDev = new ArrayList<Integer>();
    private List<Integer> ratioRSD = new ArrayList<Integer>();
    private int outlierIndex =-1;
    private int descriptionIndex =-1;
    
    
    
    private List<ProteinMedian> medianList = new ArrayList<ProteinMedian>();
    ProteinMedian proteinMedian = null;
    
    private List<ProteinRatio> ratioOverList = new ArrayList<ProteinRatio>();
    ProteinRatio proteinRatio = null;
    
    public TmtfilterProteinReader(String filepath) {
        try {
            this.filePath = filepath;
            br = new BufferedReader(new FileReader(new File(filepath)));

        } catch (FileNotFoundException ex) {
            Logger.getLogger(TmtfilterPeptideReader.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
   
        public static void main(String argsp[])
    {
        //TmtfilterProteinReader reader = new TmtfilterProteinReader("E:\\tmt\\jolene\\cesus-out-proteinAnalysisSimple.txt");
        TmtfilterProteinReader reader = new TmtfilterProteinReader("/data/2/rpark/ip2_data//xudong/TMTMS3/compare/quant/quant_compare3389.txt");
//        TmtfilterProteinReader reader = new TmtfilterProteinReader("/data/2/rpark/ip2_data//yrc/EvoProteomics/compare/quant/quant_compare2487.txt");
        //reader.getProteinList();
        reader.init();
        System.out.println(reader.getSampleList());

//        List<Protein> list = reader.getProteinList();
 //       System.out.println(list.get(0).getLocus());
    }
    
    
     public List<Protein> init()
    {
        
            readHeader();
            List<Protein> proteinList = new ArrayList<Protein>();
            while(currentLine!=null)
            {
                try 
                {
                    proteinList.add(getProtein(currentLine));
                    currentLine=br.readLine();
                } 
                catch (IOException ex) 
                {
                    Logger.getLogger(TmtfilterPeptideReader.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
            return proteinList;
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
                    
                if(currentLine.contains("Locus") || currentLine.contains("Accession"))
                    parseHeaderLine();
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
            if(words[i].equalsIgnoreCase("locus") || words[i].equalsIgnoreCase("Accession") )
                locusIndex = i;
            else if(words[i].toLowerCase().equals("spec_count"))
                specCountIndex = i;
            else if(words[i].toLowerCase().contains("pep_num"))
                pepNumIndex = i;
            else if(words[i].toLowerCase().contains("pvalue") && !words[i].toLowerCase().contains("ratio"))
                 pValueIndex.add(i);
            else if(words[i].toLowerCase().contains("pvalueratio_"))
                 pValueRatioIndex.add(i);
            else if(words[i].toLowerCase().contains("pvalueratiodata_") )
                 pValueRatioDataIndex.add(i);
            else if(words[i].toLowerCase().contains("bhcorrection") && !words[i].toLowerCase().contains("ratio"))
                 bHCorrectionIndex.add(i);
            else if(words[i].toLowerCase().contains("bhcorrectionratio") )
                 bHCorrectionRatioIndex.add(i);
            else if(words[i].contains("intensity_avg") && intensityAvgIndex== -1)
                intensityAvgIndex = i;
            else if(words[i].toLowerCase().contains("com") && words[i].toLowerCase().contains("ratio") )
                ratioSumIndex.add(i);
            else if(words[i].contains("avg")  && !words[i].contains("intensity") )
                ratioAvg.add(i);
            else if(words[i].contains("stdDev") )
                ratioStdDev.add(i);
            else if(words[i].contains("RSDValue"))
                ratioRSD.add(i);
            else if(words[i].contains("Outlier"))
                outlierIndex = i;
            else if(words[i].contains("DESCRIPTION"))
                descriptionIndex = i;
            else if(words[i].contains("median")){
            	proteinMedian = new ProteinMedian();
            	proteinMedian.setMedianHeader(words[i]);
            	proteinMedian.setIndex(i);
            	medianList.add(proteinMedian);
            }
            else if(words[i].startsWith("ratio")){
            	proteinRatio = new ProteinRatio();
            	proteinRatio.setRatioHeader(words[i]);
            	proteinRatio.setIndex(i);
            	ratioOverList.add(proteinRatio);
            }
        }
     }
     
     private Protein getProtein(String currentLine)
     {
         Protein protein = new Protein();
         ProteinMedian median = null;
         ProteinRatio ratio = null;
         String words[]=currentLine.split("\t");
         
          if(locusIndex != -1)
                protein.setLocus(words[locusIndex]);
          if(specCountIndex != -1)
                protein.setSpecCount(words[specCountIndex]);
          if(pepNumIndex != -1)
                protein.setPepNumber(words[pepNumIndex]);
          if(pValueIndex.size()!=0)
          {
              for(int i=0; i<pValueIndex.size();i++)
                protein.addPValueList(Double.parseDouble(words[pValueIndex.get(i)]));
          }
          if(pValueRatioIndex.size()!=0)
          {
              for(int i=0; i<pValueRatioIndex.size();i++)
                protein.addPValueRatioList(Double.parseDouble(words[pValueRatioIndex.get(i)]));
          }
          
          if(pValueRatioDataIndex.size()!=0)
          {
               for(int i =0 ;i < pValueRatioDataIndex.size();i++)
            {
                List<Double> data = new ArrayList<>();
                for(String value : words[pValueRatioDataIndex.get(i)].split(","))
                   data.add(Double.parseDouble(value));
                protein.getpValueRatioData().add(data);
            }
          }
          if(bHCorrectionIndex.size()!=0)
          {
              for(int i=0; i<bHCorrectionIndex.size();i++)
                protein.addBHCorrectionList(Double.parseDouble(words[bHCorrectionIndex.get(i)]));
          }
          
          if(bHCorrectionRatioIndex.size()!=0)
          {
              for(int i=0; i<bHCorrectionRatioIndex.size();i++)
                protein.addBHCorrectionRatioList(Double.parseDouble(words[bHCorrectionRatioIndex.get(i)]));
          }
          
          if(intensityAvgIndex != -1)
          {
              List<List<Double>> intensityAvgList = new ArrayList<>();

            for (int i = 0; i < sampleList.size(); i++) {

                List<Double> currentIntensityAvg = new ArrayList<>();
                for (String value : words[intensityAvgIndex+i].split(",")) {
                    currentIntensityAvg.add(Double.parseDouble(value));
                }
                intensityAvgList.add(currentIntensityAvg);
            }
            protein.setIntensityAvg(intensityAvgList);
          }
          
          if(ratioSumIndex.size()!=0)
          {
              for(int i=0; i<ratioSumIndex.size();i++)
                protein.addIntensitySumRatio(Double.parseDouble(words[ratioSumIndex.get(i)]));
          }
          
          if(ratioAvg.size()!=0)
          {
              for(int i=0; i<ratioAvg.size();i++)
                protein.addRatioAvg(Double.parseDouble(words[ratioAvg.get(i)]));
          }
          if(ratioStdDev.size()!=0)
          {
              for(int i=0; i<ratioStdDev.size();i++)
                protein.addStdDev(Double.parseDouble(words[ratioStdDev.get(i)]));
          }
          if(ratioRSD.size()!=0)
          {
              for(int i=0; i<ratioRSD.size();i++)
                protein.addRsdValue(Double.parseDouble(words[ratioRSD.get(i)]));
          }
          if(outlierIndex != -1)
                protein.setOutlier(Integer.parseInt(words[outlierIndex]));
          if(descriptionIndex != -1)
                protein.setDescription(words[descriptionIndex]);
         
        
			if (medianList.size() != 0) {
				List<ProteinMedian> proteinMedianList = new ArrayList<ProteinMedian>();
				for (int i = 0; i < medianList.size(); i++) {
					median =  new ProteinMedian();
					median.setIndex(medianList.get(i).getIndex());
					median.setMedianHeader(medianList.get(i).getMedianHeader());
//					median.getMedianValues().add(words[proteinMedian.getIndex()]);
					median.setMedianValue(words[medianList.get(i).getIndex()]);
					proteinMedianList.add(median);
				}
				protein.setMedianList(proteinMedianList);
			}
			
			if(ratioOverList.size() != 0){
				List<ProteinRatio> proteinRatioList = new ArrayList<ProteinRatio>();
				
				for (int i = 0; i < ratioOverList.size(); i++) {
					ratio = new ProteinRatio();
					ratio.setIndex(ratioOverList.get(i).getIndex());
					ratio.setRatioHeader(ratioOverList.get(i).getRatioHeader());
					ratio.setRationValue(words[ratioOverList.get(i).getIndex()]);
					proteinRatioList.add(ratio);
				}
				
				protein.setRatioOverList(proteinRatioList);
				
			}
          
         return protein;
                 
     }
     
     public Iterator<Protein> getProteins()
    {
        return init().iterator();
    }
    
     public List<Protein> getProteinList()
    {
        return init();
    }

    public String getFilePath() {
        return filePath;
    }

    public void setFilePath(String filePath) {
        this.filePath = filePath;
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


	public List<ProteinMedian> getMedianList() {
		return medianList;
	}


	public void setMedianList(List<ProteinMedian> medianList) {
		this.medianList = medianList;
	}

     
}
