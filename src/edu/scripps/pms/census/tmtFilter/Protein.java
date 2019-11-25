package edu.scripps.pms.census.tmtFilter;

import edu.scripps.pms.stats.GrubbsTest;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.inference.TestUtils;

/**
 *
 * @author Harshil
 */
public class Protein 
{
    private String locus;
    private int specCount;
    private int pepNumber;
    private String description;
    private Map<String,Peptide> peptideMap = new HashMap<>();
    private List ratioAvg = new ArrayList();
    private List ratioList = new ArrayList();
    private List stdDev = new ArrayList();
    private List rsdValue = new ArrayList();
    private List<DescriptiveStatistics> stat = new ArrayList();
    private int outlier = 0;
//    private double pValue=0;
    private List<Double> pValueList = new ArrayList<Double>();
    private List<Double> pValueRatioList = new ArrayList<Double>();
    private List<List<Double>> pValueRatioData = new ArrayList<>();
    private List<Double> BHCorrectionList = new ArrayList<Double>();
    private List<Double> BHCorrectionRatioList = new ArrayList<Double>();
    private List<Double> intensitySumRatio = new ArrayList<Double>();
    private List<List<Double>> intensityAvg = new ArrayList<List<Double>>();
    private int seqcount;
    private String seqcoverage;
    private int length;
    private int molWT;
    private double pI;

    public int getSeqcount() {
        return seqcount;
    }

    public void setSeqcount(int seqcount) {
        this.seqcount = seqcount;
    }

    public String getSeqcoverage() {
        return seqcoverage;
    }

    public void setSeqcoverage(String seqcoverage) {
        this.seqcoverage = seqcoverage;
    }

    public int getLength() {
        return length;
    }

    public void setLength(int length) {
        this.length = length;
    }

    public int getMolWT() {
        return molWT;
    }

    public void setMolWT(int molWT) {
        this.molWT = molWT;
    }

    public double getpI() {
        return pI;
    }

    public void setpI(double pI) {
        this.pI = pI;
    }
    
    
    /*
     * 
     */
    
    private List<ProteinMedian> medianList = new ArrayList<ProteinMedian>();
    private List<ProteinRatio> ratioOverList = new ArrayList<ProteinRatio>();
    
    
/*
    
    This method generate the peptide intensity averagers that we generated during 
    peptide level analysis
    */    
    public void generateProteinIntAvg()
    {
        double total = 0;
        List<Peptide> peptideList = new ArrayList<>(peptideMap.values());
        if(peptideList.size() <=0)
            return;
        int size =peptideList.get(0).getIntensityAverage().size();
        int peptideCount = peptideList.size();
        int groupSize = peptideList.get(0).getIntensityAverage().size();
        
        List<Double> proteinIntAvg = new ArrayList<Double>();
        
            total =0;
            List<List<Double>> intensityList = new ArrayList<>();
            for(int j =0; j<groupSize;j++)
            {
                List<Double> currentGroup = new ArrayList<>();
                for(int k = 0; k < peptideList.get(0).getIntensityAverage().get(j).size(); k++)
                {
                    total =0;
                    int counter=0;
                    for(Peptide currentPeptide : peptideList)
                    {
                     //   counter++;
                        total += currentPeptide.getIntensityAverage().get(j).get(k);
                    }
                //    currentGroup.add(total/counter);
                    currentGroup.add(total);
                }
                intensityList.add(currentGroup);
            }
            setIntensityAvg(intensityList);
            generateIntensitySumRatio();

    }
    
    public void generateIntensitySumRatio ()
    {
        int limit = intensityAvg.size();
         for(int i=0;i<limit;i++)
         {
             double currentTotal = 0;
             for(Double value : intensityAvg.get(i))
                 currentTotal+= value;
             for(int j=i+1;j<limit;j++)
             {
                 double nextTotal=0;
                 for(Double value : intensityAvg.get(j))
                     nextTotal+= value;
                 double ratio = currentTotal/nextTotal;
                 getIntensitySumRatio().add(ratio);
             }
         }
    }

    
    public List<Double> getpValueList() {
        return pValueList;
    }

    public void addPValueList(double pValue)
    {
        pValueList.add(pValue);
    }
    
    
     
    public void addIntensitySumRatio(double pValue)
    {
        intensitySumRatio.add(pValue);
    }
    public void addPValueRatioList(double pValue)
    {
        pValueRatioList.add(pValue);
    }
    
    public void addRatioAvg(double pValue)
    {
        ratioAvg.add(pValue);
    }
    
//    public void addIntensityAvg(List<Double> value)
//    {
//        intensityAvg.add(value);
//    }

    public List<Double> getBHCorrectionList() {
        return BHCorrectionList;
    }
    
     public void addBHCorrectionList(double value)
    {
        BHCorrectionList.add(value);
    }
     
     public void addBHCorrectionRatioList(double value)
    {
        getBHCorrectionRatioList().add(value);
    }

    
    //    public void  generatePValue()
//    {
//        if (peptideMap.size() > 0) {
//            List<Double> pValueList = new ArrayList<>();
//
//            for (String peptideName : peptideMap.keySet()) 
//            {
//                Peptide currentPeptide = peptideMap.get(peptideName);
//                pValueList.add(currentPeptide.getpValue());
//
//            }
//            if(pValueList.size()>1)
//                pValue = Fisher.combinedProbabilityTest(pValueList);
//        }
//    }

    public void setOutlier(int outlier) {
        this.outlier = outlier;
    }
    public boolean isNil = false;
    
//    public void addAverage(double average)
//    {
//        this.ratioAvg.add(average);
//        this.ratioList.add(average);
//    }
    
    public void addStdDev(double stdDev)
    {
        this.stdDev.add(stdDev);
    }
    
    public void addRsdValue(double rsdValue)
    {
        this.rsdValue.add(rsdValue);
    }
    
    public String getLocus() {
        return locus;
    }

  
    public void setLocus(String locus) {
        this.locus = locus;
    }

   
    public int getSpecCount() {
        return specCount;
    }

  
    public void setSpecCount(String specCount) {
        this.specCount = Integer.parseInt(specCount);
    }

  
    public int getPepNumber() {
        return pepNumber;
    }

   
    public void setPepNumber(String pepNumber) {
        this.pepNumber = Integer.parseInt(pepNumber);
    }

 
    public String getDescription() {
        return description;
    }

   
    public void setDescription(String description) {
        this.description = description;
    }

    public Map<String,Peptide> getPeptideMap() {
        return peptideMap;
    }

    
    public void setPeptideMap(Map<String,Peptide> peptideMap) {
        this.peptideMap = peptideMap;
    }
    
    public void generateRatio()
    {
        if(peptideMap.isEmpty())
            return;
         List keys = new ArrayList(peptideMap.keySet());
         Iterator itr = keys.iterator();
//         if(locus.equals("sp|Q6GMV3|CB079_HUMAN"))
//            System.out.println(locus);
         int size = peptideMap.get(keys.get(0)).getRatioAvg().size();
         DescriptiveStatistics descriptiveStatistics = new DescriptiveStatistics();
        for(int i=0;i<size;i++)
        {
            descriptiveStatistics = new DescriptiveStatistics();
            itr = keys.iterator();
             while(itr.hasNext())
            {
                String currentSequence = (String) itr.next();
                if(peptideMap.get(currentSequence) != null)
                {
                     descriptiveStatistics.addValue((double) peptideMap.get(currentSequence).getRatioAvg().get(i));
                }
            }
             stat.add(descriptiveStatistics);
        }
         for(int i=0;i<stat.size();i++)
         {
             getRatioAvg().add(stat.get(i).getMean());
             getStdDev().add(stat.get(i).getStandardDeviation());
             getRsdValue().add(stat.get(i).getStandardDeviation() * 100 / stat.get(i).getMean());
             
//             getRatioAvg().add(Math.pow(2,stat.get(i).getMean()));
//             getStdDev().add(Math.pow(2,stat.get(i).getStandardDeviation()));
//             getRsdValue().add(Math.pow(2,stat.get(i).getStandardDeviation() * 100 / stat.get(i).getMean()));
         }
    }
    public void filterValues(double pValue)
    {
        if(peptideMap.isEmpty())
        {
            isNil = true;
            return;
        }
        List keys = new ArrayList(peptideMap.keySet());
        List<Double> ratios = new ArrayList<>();
        Set removingIndex = new HashSet();
        Iterator itr = keys.iterator();
        int size= peptideMap.get(keys.get(0)).getRatioAvg().size();

        for(int i=0;i<size;i++)
        {
            ratios = new ArrayList<>();
            while(itr.hasNext())
            {
                String currentSequence = (String) itr.next();
                if(peptideMap.get(currentSequence) != null)
                {
                    double value = Math.log((double) peptideMap.get(currentSequence).getRatioAvg().get(i))/ Math.log(2);
                    ratios.add(value);
//                    ratios.add((Double) peptideMap.get(currentSequence).getRatioAvg().get(i));
                }
//                else
//                    System.out.println("");
            }
            removingIndex.addAll(GrubbsTest.filterAndRemovePeptideFromProteinNoLog(ratios,pValue));
        }
        itr = removingIndex.iterator();
        while(itr.hasNext())
        {
            outlier++;
            peptideMap.remove(keys.get((int)itr.next()));
        }
//        System.out.println("Total Entries removed from at protein level" + removingIndex.size());
    }

 
   

    public List getRatioAvg() {
        return ratioAvg;
    }

  
    public List getStdDev() {
        return stdDev;
    }

   
    public List getRsdValue() {
        return rsdValue;
    }

    public int getOutlier() {
        return outlier;
    }

    public List getRatioList() {
        return ratioList;
    }

    public void setRatioList(List ratioList) {
        this.ratioList = ratioList;
    }

    /**
     * @return the intensityAvg
     */
    public List<List<Double>> getIntensityAvg() {
        return intensityAvg;
    }

    
    /**
     * @param proteinIntensityAvg the intensityAvg to set
     */
    public void setIntensityAvg(List<List<Double>> proteinIntensityAvg) {
        this.intensityAvg = proteinIntensityAvg;
    }

    public void setSpecCount(int specCount) {
        this.specCount = specCount;
    }

    public void setPepNumber(int pepNumber) {
        this.pepNumber = pepNumber;
    }

    public void setRatioAvg(List ratioAvg) {
        this.ratioAvg = ratioAvg;
    }

    public void setStdDev(List stdDev) {
        this.stdDev = stdDev;
    }

    public void setRsdValue(List rsdValue) {
        this.rsdValue = rsdValue;
    }

    public void setStat(List<DescriptiveStatistics> stat) {
        this.stat = stat;
    }

    public void setpValueList(List<Double> pValueList) {
        this.pValueList = pValueList;
    }

    public void setBHCorrectionList(List<Double> BHCorrectionList) {
        this.BHCorrectionList = BHCorrectionList;
    }

    public void setIsNil(boolean isNil) {
        this.isNil = isNil;
    }

    /**
     * @return the pValueRatioList
     */
    public List<Double> getpValueRatioList() {
        return pValueRatioList;
    }

    /**
     * @param pValueRatioList the pValueRatioList to set
     */
    public void setpValueRatioList(List<Double> pValueRatioList) {
        this.pValueRatioList = pValueRatioList;
    }
    
    public void generatepValueRatioList(double[] pValueRatioData,int totalPeptide) {
        List<Double> data = new ArrayList<>();
        for(int i =0;i<pValueRatioData.length;i++)
        {
            pValueRatioData[i] /=totalPeptide;
            data.add(pValueRatioData[i]);
        }
        this.getpValueRatioData().add(data);
        try {
            double d = TestUtils.tTest(0, pValueRatioData);
            if(Double.isNaN(d))
                d=1.0;
            this.pValueRatioList.add(d);
        } catch (IllegalArgumentException ex) {
            Logger.getLogger(Protein.class.getName()).log(Level.SEVERE, null, ex);
        } 
        
    }

    /**
     * @return the BHCorrectionRatioList
     */
    public List<Double> getBHCorrectionRatioList() {
        return BHCorrectionRatioList;
    }

    /**
     * @param BHCorrectionRatioList the BHCorrectionRatioList to set
     */
    public void setBHCorrectionRatioList(List<Double> BHCorrectionRatioList) {
        this.BHCorrectionRatioList = BHCorrectionRatioList;
    }

    /**
     * @return the pValueRatioData
     */
    public List<List<Double>> getpValueRatioData() {
        return pValueRatioData;
    }

    /**
     * @param pValueRatioData the pValueRatioData to set
     */
    public void setpValueRatioData(List<List<Double>> pValueRatioData) {
        this.pValueRatioData = pValueRatioData;
    }

    /**
     * @return the intensitySumRatio
     */
    public List<Double> getIntensitySumRatio() {
        return intensitySumRatio;
    }

    /**
     * @param intensitySumRatio the intensitySumRatio to set
     */
    public void setIntensitySumRatio(List<Double> intensitySumRatio) {
        this.intensitySumRatio = intensitySumRatio;
    }

	public List<ProteinMedian> getMedianList() {
		return medianList;
	}

	public void setMedianList(List<ProteinMedian> medianList) {
		this.medianList = medianList;
	}

	public List<ProteinRatio> getRatioOverList() {
		return ratioOverList;
	}

	public void setRatioOverList(List<ProteinRatio> ratioOverList) {
		this.ratioOverList = ratioOverList;
	}
    
    
    
}
