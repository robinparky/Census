package edu.scripps.pms.census.tmtFilter;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

/**
 *
 * @author Harshil
 * 
 * This method contains the Each row of the Peptide sequence from the "census-out.txt" file.
 */
//public class Peptide extends DescriptiveStatistics
public class Peptide
{
    private Map<String, ArrayList> all_intensity = new HashMap();// This one will store all the intensities.
    private Map<String, List<Double>> allNormIntensity = new HashMap<>();
    private List<List<Double>> intensityAverage = new ArrayList();
    private List<List<Double>> normIntensitySum = new ArrayList();

    private Map<String, DescriptiveStatistics> stats=new HashMap();// The key is name of the test carried on the peptide. 
    private int spc=0;
    private String sequence;
    private ArrayList<Integer> scannum = new ArrayList();
    private ArrayList<Integer> cstate = new ArrayList();
    private Set<String> filename = new HashSet();
    private List ratioAvg = new ArrayList();
    private List ratioNormAvg = new ArrayList();
    private List logRatioStdDev = new ArrayList();
    private List logRatioNormStdDev = new ArrayList();
    private List logRatioRSD = new ArrayList();
    private List logRatioNormRSD = new ArrayList();
    private List logRatioAvg = new ArrayList();
    private List logRatioNormAvg = new ArrayList();

    private List<List> ratioList = new ArrayList();//only used for the TmtFiletrPeptideReader
    private List<List> normRatioList = new ArrayList();//only used for the TmtFiletrPeptideReader

    public List ratioListForPvalue = new ArrayList();// these are the data on which we calculate the pValueRatioList
    public List<List<DescriptiveStatistics>> normRatioListForPValue = new ArrayList();// these are the data on which we calculate the pValueRatioList

    private List<List<Double>> ratioPvalueData = new ArrayList();//only used for the TmtFiletrPeptideReader
    private List<List<Double>> ratioNormPvalueData = new ArrayList();//only used for the TmtFiletrPeptideReader

    public List<List<DescriptiveStatistics>> comparedGroup = new ArrayList();// its list of each comparing group
    public List<List<DescriptiveStatistics>> normComparedGroup = new ArrayList();// its list of each comparing group

    private int peptideCounter = 0 ;
    private int discaredCounter = 0;
    private Map<String, ArrayList> ratio=new HashMap();
    private DescriptiveStatistics stdev = new DescriptiveStatistics();
    private Set<String> proteinLocus = new LinkedHashSet<String>();
    private Set<String> proteinDescription = new LinkedHashSet<String>();
    //Added By Harshil Shah 
    private List<Double> pValueList = new ArrayList<Double>();
    private List<Double> normPValueList = new ArrayList<Double>();

    private double pValue=0;
    private List<Double> pValueRatioList = new ArrayList<Double>();
    private List<Double> pValueNormRatioList = new ArrayList<Double>();

    private List<Double> bHCorrection = new ArrayList();
    private List<Double> bHCorrectionNorm = new ArrayList();

    private List<Double> bHCorrectionRatioList = new ArrayList();
    private List<Double> bHCorrectionNormRatioList = new ArrayList();

    //Just for minimum sized stdev and rsd value:....
    private List<Double> rsdPaired = new ArrayList<Double>();
    private List<Double> normRSDPaired = new ArrayList<>();
    private List<Double> stdevPaired = new ArrayList<Double>();
    private List<Double> normStdevPaired = new ArrayList<>();

    private Set<String> geneNameSet = new LinkedHashSet<>();
    private String geneListStr = "";
    private String proteinLocusListStr = "";
    private String proteinDescriptionListStr="";

    public String getGeneListStr() {
        return geneListStr;
    }

    public void setGeneListStr(String geneListStr) {
        this.geneListStr = geneListStr;
    }

    public String getProteinLocusListStr() {
        return proteinLocusListStr;
    }

    public void setProteinLocusListStr(String proteinLocusListStr) {
        this.proteinLocusListStr = proteinLocusListStr;
    }

    public String getProteinDescriptionListStr() {
        return proteinDescriptionListStr;
    }

    public void setProteinDescriptionListStr(String proteinDescriptionListStr) {
        this.proteinDescriptionListStr = proteinDescriptionListStr;
    }

    public Set<String> getGeneNameSet() {
        return geneNameSet;
    }

    public void setGeneNameSet(Set<String> geneNameSet) {
        this.geneNameSet = geneNameSet;
    }

    public void addGeneName(String geneName)
    {
        geneNameSet.add(geneName);
    }

    public void addRsdPaired(Double value)
    {
        rsdPaired.add(value);
    }
    
    public void addstdevPaired(Double value)
    {
        stdevPaired.add(value);
    }

    public void addNormRsdPaired(Double value)
    {
        normRSDPaired.add(value);
    }

    public void addNormstdevPaired(Double value)
    {
        normStdevPaired.add(value);
    }

    public void addRatioPvalueData(List<Double> value)
    {
        ratioPvalueData.add(value);
    }

    public void addRationNormPvalueData(List<Double> value)
    {
        ratioNormPvalueData.add(value);
    }



    public void addBHCorrection(double value)
    {
        bHCorrection.add(value);
    }

    public void addBHCorrectionNorm(double value)
    {
        bHCorrectionNorm.add(value);
    }
    
    public void addBHCorrectionRatio(double value)
    {
        bHCorrectionRatioList.add(value);
    }

    public void addBHCorrectionNormRatio(double value)
    {
        bHCorrectionNormRatioList.add(value);
    }


    public List<Double> getbHCorrection() {
        return bHCorrection;
    }

    public List<Double> getbHCorrectionNorm() {
        return bHCorrectionNorm;
    }

        
    public void addProteinLocus(String proteinName)
    {
        proteinLocus.add(proteinName);
    }
    
    public void addPValueList (double value)
    {
        pValueList.add(value);
    }


    public void addNormPValueList (double value)
    {
        normPValueList.add(value);
    }

    public void addPValueRatio (double value)
    {
        pValueRatioList.add(value);
    }

    public void addPValueNormRatio (double value)
    {
        pValueNormRatioList.add(value);
    }

    public  List<Double> getpValueNormRatioList()
    {
        return pValueNormRatioList;
    }

    public Map<String, List<Double>> getAllNormIntensity() {
        return allNormIntensity;
    }

    public List<List> getNormRatioList() {
        return normRatioList;
    }

    public List getRatioListForPvalue() {
        return ratioListForPvalue;
    }

    public List<List<DescriptiveStatistics>> getNormRatioListForPValue() {
        return normRatioListForPValue;
    }

    public List<List<Double>> getRatioNormPvalueData() {
        return ratioNormPvalueData;
    }

    public List<List<DescriptiveStatistics>> getNormComparedGroup() {
        return normComparedGroup;
    }

    public List<Double> getNormRSDPaired() {
        return normRSDPaired;
    }

    public void addProteinDescription(String description)
    {
        proteinDescription.add(description);
    }
    
    public void removeEntry(Set index)
    {
//        List key = new ArrayList(all_intensity.keySet());
        List key = new ArrayList<>(all_intensity.keySet());
        List indexes = new ArrayList(index);
        Collections.sort(indexes);
        Collections.reverse(indexes);
        Iterator itr= key.iterator();
        Iterator removinIndex = indexes.iterator();
        while(removinIndex.hasNext())
        {
//            System.out.println(sequence);
         
            itr= all_intensity.keySet().iterator();
            int indexToBERemoved = (int) removinIndex.next();
//            if(indexToBERemoved == 11)
//                System.out.println();
            while(itr.hasNext())
            {
                String removeKey = (String)itr.next();
                all_intensity.get(removeKey).remove(indexToBERemoved);
                allNormIntensity.get(removeKey).remove(indexToBERemoved);
            }

            scannum.remove(indexToBERemoved);
            cstate.remove(indexToBERemoved);
//            descriptiveStatisticses.remove(indexToBERemoved);
            comparedGroup.remove(indexToBERemoved);
            peptideCounter--;
            discaredCounter++;
        }
    }



    
    public void setOutlier(int outlier)
    {
        discaredCounter = outlier;
    }

    public void setAllNormIntensity(String key, double value)
    {
        List list = new ArrayList();
        if(allNormIntensity.containsKey(key))
        {
            list = allNormIntensity.get(key);
        }
        list.add(value);
        allNormIntensity.put(key,list);
    }


    public void setallIntensity(String key, double value)
    {
        ArrayList ar=new ArrayList();
        
        if(all_intensity.containsKey(key))
        {
            ar=all_intensity.get(key);
        }
  
        ar.add(value);
        
        all_intensity.put(key, ar);
    }
    
    void printratio()
    {
        System.out.println(ratio);
    }
    public ArrayList getallIntensity(String key)
    {
        return all_intensity.get(key);
    }
     public Map getallIntensity()
    {
        return all_intensity;
    }
    
    public void addFilename(String fileName)
    {
        this.filename.add(fileName);
    }
    
    public void addRatio(List ratio)
    {
        this.ratioList.add(ratio);
    }

    public void addNormRatio(List ratio)
    {
        this.normRatioList.add(ratio);
    }


    
    public void addRatioAvg(double ratioAvg)
    {
        this.ratioAvg.add(ratioAvg);
    }

    public void addNormRatioAvg(double ratioAvg)
    {
        this.ratioNormAvg.add(ratioAvg);
    }


    public void addLogRatioStdDev(double ratioStdDev)
    {
        
        this.logRatioStdDev.add(ratioStdDev);
    }
    public void addNormLogRatioStdDev(double ratioStdDev)
    {

        this.logRatioNormStdDev.add(ratioStdDev);
    }


    public void addLogRatioAvg(double logRatioAvg)
    {
        
        this.logRatioAvg.add(logRatioAvg);
    }

    public void addLogRatioNormAvg(double logRatioAvg)
    {

        this.logRatioNormAvg.add(logRatioAvg);
    }



    public void addLogRatioRSD(double ratioRSD)
    {
        this.logRatioRSD.add(ratioRSD);
    }

        public void addLogRatioNormRSD(double ratioRSD)
    {
        this.logRatioNormRSD.add(ratioRSD);
    }
    
    public int getSpc() 
    {
        return spc;
    }

    public void setSpc(int spc)
    {
        this.spc += spc;
    }

    public ArrayList<Integer> getScannum() 
    {
        return scannum;
    }


    public void addScannum(int number)
    {
        scannum.add(number);
        
    }

 
    public ArrayList<Integer> getCstate() 
    {
        return cstate;
    }

    public void setCstate(int number)
    {
        cstate.add(number);
    }

    public Set<String> getFilename()
    {
        return filename;
    }


    public void setFilename(String filename)
    {
        this.filename.add(filename);
    }


    public DescriptiveStatistics getStat(String key) 
    {
        return stats.get(key);
    }
    public void setstat(String key,  DescriptiveStatistics stat) 
    {
        this.stats.put(key, stat);
    }

 
    public int getPeptideCounter() {
        return peptideCounter;
    }

    public void setPeptideCounter() {
        this.peptideCounter++;
    }

  
    public String getSequence() {
        return sequence;
    }

    public void setSequence(String sequence) {
        this.sequence = sequence;
    }

   
    public List getComparedGroup() {
        return comparedGroup;
    }
    
//    public void addComparedGroup(List stat){
//        comparedGroup.add(stat);
//    }

    public boolean generateRatio()
    {
        /*
            this will generate the details for the ratioAvg and ratioStdDev.
            for a single peptide Entry....
        */
        DescriptiveStatistics stat = new DescriptiveStatistics();
        DescriptiveStatistics unlogStat = new DescriptiveStatistics();
        DescriptiveStatistics normStat = new DescriptiveStatistics();
        
        int size1 = comparedGroup.size();
        if(size1 == 0)
           return false;
        int size2 = comparedGroup.get(0).size();

        for(int z =0;z<size2 ; z++)
        {
            stat = new DescriptiveStatistics();
            unlogStat = new DescriptiveStatistics();
            normStat = new DescriptiveStatistics();
            for(int j =0;j<size1 ; j++)
            {
                stat.addValue(comparedGroup.get(j).get(z).getMean());
                normStat.addValue(normComparedGroup.get(j).get(z).getMean());
                double value = Math.pow(2, comparedGroup.get(j).get(z).getMean());
                unlogStat.addValue(value);
                
            }
            
            double unlogRatio = Math.pow(2, stat.getMean());
            double unlogNormRatio = Math.pow(2,normStat.getMean());
            ratioAvg.add(unlogRatio);
            ratioNormAvg.add(unlogNormRatio);
//            unlogRatio = Math.pow(2, stat.getStandardDeviation());
//            ratioRSD.add(unlogStat.getStandardDeviation()/unlogStat.getMean()*100.00);
              
            //unlogRatio = Math.pow(2, stat.getStandardDeviation());
            
            // All the log values are being averaged...
/*
so if we have 2 groups of 3 experiment then the stat will have all
9 possible combinations from 2 groups..... each is a ratio,,,,
*/
            logRatioAvg.add(stat.getMean());
            logRatioNormAvg.add(normStat.getMean());
            double x = stat.getStandardDeviation()/stat.getMean()*100.00;
            logRatioRSD.add(stat.getStandardDeviation()/stat.getMean()*100.00);
            logRatioNormRSD.add(normStat.getStandardDeviation()/normStat.getMean()*100.00);
          //  System.out.println(">>>From all values "+unlogRatio+" "+stat.getMean()+" "+x+" "+stat.getStandardDeviation());
        /*   System.out.println(sequence);
            for(double value : stat.getValues())
            {
                System.out.println("\t" +value  );
            }*/
            logRatioStdDev.add(stat.getStandardDeviation());
            logRatioNormStdDev.add(normStat.getStandardDeviation());
        }
        return true;
    }

    public void generateIntensityAverage(List<List<String>> headerName)
    {
        
//        for (String key : headerName) 
//        {
//            double total = 0;
//            for(Double value : (List<Double>) all_intensity.get(key))
//            {
//                total += value;
//            }
//            getIntensityAverage().add(total);
//        }
        for(List<String> currentGroup : headerName)
        {
            List<Double> intensityList =new ArrayList<>();;
            List<Double> normIntensityList = new ArrayList<>();
            for(String currentKey : currentGroup)
            {
                double total =0;
                double normTotal =0;
                int i=0;
                for(; i < all_intensity.get(currentGroup.get(0)).size();i++)
                {
                    total+=(Double)all_intensity.get(currentKey).get(i);
                }
                for(int j=0; j<allNormIntensity.get(currentGroup.get(0)).size(); j++)
                {
                    normTotal+=allNormIntensity.get(currentKey).get(j);
                }
                //intensityList.add(total/i);
                intensityList.add(total);
                normIntensityList.add(normTotal);
            }
            getIntensityAverage().add(intensityList);
            getNormIntensitySum().add(normIntensityList);
        }
        
    }


 

    public Set<String> getProteinLocus() {
        return proteinLocus;
    }
    
    public Set<String> getProteinDescription() {
        return proteinDescription;
    }

    /**
     * @return the discaredCounter
     */
    public int getDiscaredCounter() {
        return discaredCounter;
    }

    /**
     * @return the pValueList
     */
    public List<Double> getpValueList() {
        return pValueList;
    }

    public List<Double> getNormPValueList() {
        return normPValueList;
    }
    /**
     * @param pValueList the pValueList to set
     */
    public void setpValueList(List<Double> pValueList) {
        this.pValueList = pValueList;
    }

    /**
     * @return the pValue
     */
    public double getpValue() {
        return pValue;
    }

    /**
     * @param pValue the pValue to set
     */
    public void setpValue(double pValue) {
        this.pValue = pValue;
    }

    public Map<String, ArrayList> getAll_intensity() {
        return all_intensity;
    }

    public void setAll_intensity(Map<String, ArrayList> all_intensity) {
        this.all_intensity = all_intensity;
    }

    public Map<String, DescriptiveStatistics> getStats() {
        return stats;
    }

    public void setStats(Map<String, DescriptiveStatistics> stats) {
        this.stats = stats;
    }

    public List<List> getRatioList() {
        return ratioList;
    }

    
    public void setRatioList(List<List> ratioList) {
        this.ratioList = ratioList;
    }

    /*
    public List<List<DescriptiveStatistics>> getComparedGroup() {
        return comparedGroup;
    }*/

    public void setComparedGroup(List<List<DescriptiveStatistics>> comparedGroup) {
        this.comparedGroup = comparedGroup;
    }

    public Map<String, ArrayList> getRatio() {
        return ratio;
    }

    public void setRatio(Map<String, ArrayList> ratio) {
        this.ratio = ratio;
    }

    public DescriptiveStatistics getStdev() {
        return stdev;
    }

    public void setStdev(DescriptiveStatistics stdev) {
        this.stdev = stdev;
    }

    /**
     * @return the intensityAverage
     */
    public List<List<Double>> getIntensityAverage() {
        return intensityAverage;
    }

    public List<List<Double>> getNormIntensitySum() {
        return normIntensitySum;
    }

    /**
     * @param intensityAverage the intensityAverage to set
     */
    public void setIntensityAverage(List<List<Double>> intensityAverage) {
        this.intensityAverage = intensityAverage;
    }
    public void setNormIntensitySum(List<List<Double>> intensitySum) {
        this.normIntensitySum = intensitySum;
    }

    public List getLogRatioStdDev() {
        return logRatioStdDev;
    }
    public List getLogRatioNormStdDev() {
        return logRatioNormStdDev;
    }
    public List getLogRatioRSD() {
        return logRatioRSD;
    }
    public List getLogRatioNormRSD() {
        return logRatioNormRSD;
    }
    public List getLogRatioAvg() {
        return logRatioAvg;
    }
    public List getLogRatioNormAvg() {
        return logRatioNormAvg;
    }
    public void setScannum(ArrayList<Integer> scannum) {
        this.scannum = scannum;
    }

    public void setCstate(ArrayList<Integer> cstate) {
        this.cstate = cstate;
    }

    public void setFilename(Set<String> filename) {
        this.filename = filename;
    }

    public void setLogRatioStdDev(List logRatioStdDev) {
        this.logRatioStdDev = logRatioStdDev;
    }

    public void setLogRatioRSD(List logRatioRSD) {
        this.logRatioRSD = logRatioRSD;
    }

    public void setLogRatioAvg(List logRatioAvg) {
        this.logRatioAvg = logRatioAvg;
    }

    public void setPeptideCounter(int peptideCounter) {
        this.peptideCounter = peptideCounter;
    }

    public void setDiscaredCounter(int discaredCounter) {
        this.discaredCounter = discaredCounter;
    }

    public void setProteinLocus(Set<String> proteinLocus) {
        this.proteinLocus = proteinLocus;
    }

    public void setProteinDescription(Set<String> proteinDescription) {
        this.proteinDescription = proteinDescription;
    }

    public void setbHCorrection(List<Double> bHCorrection) {
        this.bHCorrection = bHCorrection;
    }

    public List getRatioAvg() {
        return ratioAvg;
    }
    public List<Double> getRatioNormAvg() {
        return ratioNormAvg;
    }
    public void setRatioAvg(List ratioAvg) {
        this.ratioAvg = ratioAvg;
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

    /**
     * @return the bHCorrectionRatioList
     */
    public List<Double> getbHCorrectionRatioList() {
        return bHCorrectionRatioList;
    }
    public List<Double> getbHCorrectionNormRatioList() {
        return bHCorrectionNormRatioList;
    }
    /**
     * @param bHCorrectionRatioList the bHCorrectionRatioList to set
     */
    public void setbHCorrectionRatioList(List<Double> bHCorrectionRatioList) {
        this.bHCorrectionRatioList = bHCorrectionRatioList;
    }

    /**
     * @return the ratioPvalueData
     */
    public List<List<Double>> getRatioPvalueData() {
        return ratioPvalueData;
    }

    /**
     * @param ratioPvalueData the ratioPvalueData to set
     */
    public void setRatioPvalueData(List<List<Double>> ratioPvalueData) {
        this.ratioPvalueData = ratioPvalueData;
    }

    public List<Double> getRsdPaired() {
        return rsdPaired;
    }
    public List<Double> getNormRsdPaired() {
        return normRSDPaired;
    }
    public void setRsdPaired(List<Double> rsdPaired) {
        this.rsdPaired = rsdPaired;
    }

    public List<Double> getStdevPaired() {
        return stdevPaired;
    }
    public List<Double> getNormStdevPaired() {
        return normStdevPaired;
    }
    public void setStdevPaired(List<Double> stdevPaired) {
        this.stdevPaired = stdevPaired;
    }


    public Map<String,List<Double>> getNormAllIntensity() {
        return allNormIntensity;
    }
}
