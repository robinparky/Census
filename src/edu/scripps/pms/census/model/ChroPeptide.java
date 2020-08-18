/*
 * ChroPeptide.java
 *
 * Created on May 23, 2005, 11:35 AM
 */

package edu.scripps.pms.census.model;

import java.util.*;

import java.text.DecimalFormat;
import edu.scripps.pms.census.conf.Configuration;
import edu.scripps.pms.census.CensusConstants;
import edu.scripps.pms.census.tools.IsotopeModel;
import edu.scripps.pms.census.util.IsoData;
import edu.scripps.pms.census.util.io.LuciphorReader;
import org.json.simple.JSONObject;
import edu.scripps.pms.census.util.LinearRegression;
import rpark.statistics.model.GaussianPeakModel;

/**
 *
 * @author rpark
 * @author rohan
 * @version $Id: ChroPeptide.java,v 1.19 2014/08/06 05:29:53 rpark Exp $
 */
    public class ChroPeptide {

    private boolean unique;
    private String fileName;
    private String xCorr;
    private String deltCN;
    private String deltMass;
    private String mhPlus;
    private String calcMHplus;
    private double totalIntensity;
    private double samIntensity=-1;
    private double refIntensity;
    private double ionInjectionTimeNormIntensity;
    private LuciphorReader.LuciphorPeptide luciphorPeptide  =null;

    private String spRank;
    private String spScore;
    private String specCount;

    private String ionProportion;
    private String redundancy;
    private String sequence;
    private String peptideLine;
    private String tmpStr;
    private int scanNum;
    private int chargeState;
    private String[] strArr;
    private String startRange; //start scan for peak range
    private String endRange; //end scan for peak range

    private int dtaStartRange;
    private int dtaEndRange;

    private double[][] bionSample;
    private double[][] bionRef;
    private double[][] yionSample;
    private double[][] yionRef;

    private String bsText;
    private String ysText;
    private String brText;
    private String yrText;

    private double corr; //regression factor
    private double slope=-1;
    private double corrRev; //regression factor
    private double slopeRev=-1;

    private int chroCenter=0; //position of scan number
    private boolean filterOut=false;

    private double probability;
    private double retentionTime;
    private double ionInjectionTime;

    private List<ChroData> dataList;
    private List<IsoData> isoDataList;
    private List<IsoData> isoOrigDataList;

    private TheoryData theoryData;
    private boolean missedPeptide = false;
    private List<IsotopeModel> isoTopeModelList = new ArrayList<>();

    private double snRatio;
    private double ratio;

    //nonlabeling intensity sum array
    private long[] totalIntArr;
    private double[] massMonitorArr;

    private double outlierPValue;

    private DecimalFormat formatter = new DecimalFormat("0.000");


    private double corrToNorm; //correlation to normal distribution
    private int spectraDataPoints;  //for all-none peptides
    private double massTolerance; //for all-none peptides
    private double anCompositeScore;
    private double fragIonNumFoundRatio; // (# of fragment ion found) / (total fragment ion number)
    private double averageIntensity; //intensity / peak count
    private double zscore;
    private double lightMass;
    private double heavyMass;

    private Hashtable<String, String> scoreHt = new Hashtable<String, String>();
    private Configuration conf;

    private double enrichment;

    private double bestEnrichCorr=-1;
    private double bestEnrichDelCN=-1;
    private double corrOnePlus=-1;
    private double corrOneMinus=-1;

    private double pvalue;
    private double zvalue;

    private ArrayList<Double> multiIntensity = new ArrayList<Double>();

    private double maxIntensity = -1;
    private int startRangeInt = -1;
    private int endRangeInt = -1;

    private long precursorLightIntensity=0;
    private long precursorHeavyIntensity=0;
    private List<LinearRegression> regList = new ArrayList<LinearRegression>();
    private double[] intensitySum;
  private double[] diaFragRegressionArr;
    private List<GaussianPeakModel> gaussianPeakModelList = new ArrayList<>();
    private long gaussianPeakArea;
    private long gaussianPeakAreaIonInjectionCorrection;
    private String chroData;
    private double peakSigma;
    private double peakx;
    private double peaky;
    private double startRt;
    private double endRt;
    private String gaussianPeakString;
    private boolean singleton=false;
    private double corrIonInjectionIntensity;
    private double[] isoArr;
    private int peptideIndex;
    private int sampleIndex;
    private double profileScore=0; //anCompositeScore has been used as profile score in 15N, and SILAC.  we start using this profile score in triple silace and dimethyl.
    private double tmtPurity=-1;
    private double signalNoise = -1;
    private double reportIonSignalNoise = -1;

    private String ptmIndex="NA";
    private String ptmIndexProtein="NA";

    //private List<Double> normIntensityList = new ArrayList<Double>();
    private double normIntensity;

    public double getNormIntensity() {
        return normIntensity;
    }

    public void setNormIntensity(double normIntensity) {
        this.normIntensity = normIntensity;
    }

    /*
        public void addNormIntensity(double intensity) {
            normIntensityList.add(intensity);

        }

        public List<Double> getNormIntensityList() {
            return normIntensityList;
        }

        public void setNormIntensityList(List<Double> normIntensityList) {
            this.normIntensityList = normIntensityList;
        }
    */
    public double getReportIonSignalNoise() {
        return reportIonSignalNoise;
    }

    public void setReportIonSignalNoise(double reportIonSignalNoise) {
        this.reportIonSignalNoise = reportIonSignalNoise;
    }

    public double getSignalNoise() {
        return signalNoise;
    }

    public void setSignalNoise(double signalNoise) {
        this.signalNoise = signalNoise;
    }

    public double getCorrIonInjectionIntensity() {
        return corrIonInjectionIntensity;
    }

    public void setCorrIonInjectionIntensity(double corrIonInjectionIntensity) {
        this.corrIonInjectionIntensity = corrIonInjectionIntensity;
    }



    public long getGaussianPeakAreaIonInjectionCorrection() {
        return gaussianPeakAreaIonInjectionCorrection;
    }

    public void setGaussianPeakAreaIonInjectionCorrection(long gaussianPeakAreaIonInjectionCorrection) {
        this.gaussianPeakAreaIonInjectionCorrection = gaussianPeakAreaIonInjectionCorrection;
    }


    public double getMaxIntensity() {
        if(this.maxIntensity>0)
            return this.maxIntensity;

        /*System.out.println("aaaaaa");
        System.out.println("aaaaaa" +this.startRange);
        System.out.println("aaaaaa" +this.scanNum);
        System.out.println("aaaaaa" +this.endRange);
        System.out.println("aaaaaa" + this.dataList.size());
        System.out.println("aaaaaa" + this.dataList.get(0).getScanNum());
        System.out.println("aaaaaa" + this.dataList.get(dataList.size()-1).getScanNum());*/

        if(startRangeInt<0) startRangeInt = (null==startRange)?0:Integer.parseInt(startRange);
        if(endRangeInt<0) endRangeInt = (null==endRange)?0:Integer.parseInt(endRange);

        long maxIntensity=0;
        for(ChroData cd:this.dataList) {
            if(cd.getScanNum()>=this.startRangeInt && cd.getScanNum()<=this.endRangeInt) {
                long tmpSum = cd.getSampleIntensity() + cd.getRefIntensity();

                if(tmpSum>=maxIntensity)
                    maxIntensity=tmpSum;

              ///  System.out.println(tmpSum + "\t" + maxIntensity);
            }
        }


        /*
        for(long l:totalIntArr) {
            System.out.println(l);
        }*/
        //System.out.println("aaaaaabb\t" + maxIntensity);

        return maxIntensity;
    }

    public ChroPeptide()
    {
        dataList = new ArrayList<ChroData>();
	conf = Configuration.getInstance();
    }

    public ChroPeptide(String peptideLine)
    {
        this.peptideLine = peptideLine;
        parseLine();
        dataList = new ArrayList<ChroData>();
	conf = Configuration.getInstance();
    }

    public  static void main(String args[])
    {
        System.err.println("dsada");
    }

    private void parseLine() throws ArrayIndexOutOfBoundsException
    {
        strArr = this.peptideLine.split("\t");

        //scanNum = strArr[1];
        //scanNum = scanNum.substring(0, scanNum.lastIndexOf("."));
        //scanNum = scanNum.substring(scanNum.lastIndexOf(".")+1);

        this.setUnique(!"".equals(strArr[0]));
        this.setFileName(strArr[1]);
        this.setScanNum(Integer.parseInt(strArr[2]));
        this.setSequence(strArr[3]);
        this.setXCorr(strArr[4]);
        this.setDeltCN(strArr[5]);
        this.setChargeState(strArr[6]);

        if(strArr.length>7)
        {
            this.setDtaStartRange( Integer.parseInt(strArr[7]) );
            this.setDtaEndRange( Integer.parseInt(strArr[8]) );
        }

        //this.setMhPlus(strArr[4]);
        //this.setCalcMHplus(strArr[5]);
        //this.setTotalIntensity(strArr[6]);
        //this.setSpRank(strArr[7]);
        //this.setSpScore(strArr[8]);
        //this.setIonProportion(strArr[9]);
        //this.setRedundancy(strArr[10]);
        //this.setSequence(strArr[11]);
        //this.chargeState = Integer.parseInt(
          //      this.fileName.substring(this.fileName.lastIndexOf(".") + 1 )
           //     );
        //this.chargeState = strArr[1].substring( this.strArr[1].lastIndexOf(".") + 1 );
    }

    public void addData(ChroData data)
    {
        dataList.add(data);

        if(data.getScanNum() == scanNum)
            this.chroCenter = dataList.size()-1;
    }

    public List getDataList()
    {
        return dataList;
    }

    public int getDataSize()
    {
        return dataList.size();
    }

    public ChroData getData(int pos)
    {
        return dataList.get(pos);
    }

    public String[] getPeptideData()
    {
	if(conf.getScoreType()==Configuration.PEPXML)
	{
	    List scoreList = conf.getScoreList();
	    String[] strArr = new String[CensusConstants.BASIC_PEPTIDE_COLUMN_SIZE+scoreList.size()];


            strArr[0] = isUnique()?"*":"";
	    strArr[1] = getFileName();
	    strArr[2] = String.valueOf(getScanNum());
	    strArr[3] = getSequence();
	    strArr[4] = String.valueOf(getChargeState());

	    int index=CensusConstants.BASIC_PEPTIDE_COLUMN_SIZE;
	    for(Iterator<String> itr=scoreList.iterator(); itr.hasNext(); )
	    {
		strArr[index++] = scoreHt.get(itr.next());
	    }

	    return strArr;
	}

        if(this.strArr==null)
        {
            //String[] temp = { isUnique()?"*":"", getFileName(), String.valueOf(getScanNum()), getSequence(), (this.getSlope()>0)?formatter.format(this.getSlope()):"NA", (this.getCorr()>0)?formatter.format(this.getCorr()):"NA", (this.getCorr()>0)?formatter.format(this.getCorr()*this.getCorr()):"NA", getXCorr(), String.valueOf(getDeltCN()), String.valueOf(getChargeState()), String.valueOf(getDtaStartRange()), String.valueOf(getDtaEndRange())};
            String[] temp = { isUnique()?"*":"",
                            getFileName(),
                            String.valueOf(getScanNum()),
                            getSequence(),
                            //(this.getSlope()>0)?formatter.format(this.getSlope()):"NA",
                            //(this.getCorr()>0)?formatter.format(this.getCorr()):"NA",
                            //(this.getCorr()>0)?formatter.format(this.getCorr()*this.getCorr()):"NA",
                            getXCorr(),
                            String.valueOf(getDeltCN()),
                            String.valueOf(getChargeState()) };
                           // String.valueOf(getDtaStartRange()),
                            //String.valueOf(getDtaEndRange())};
            return temp;
        }

        return this.strArr;
    }


    public String[] getLabelFreePeptideData()
    {
        if(this.strArr==null)
        {
            String[] temp = { isUnique()?"*":"", getFileName(), String.valueOf(getScanNum()), getSequence(), getXCorr(), String.valueOf(getDeltCN()), String.valueOf(getChargeState()), String.valueOf(getDtaStartRange()), String.valueOf(getDtaEndRange())};
            return temp;
        }

        return this.strArr;
    }

    public int getChargeState()
    {
        return this.chargeState;
    }

    public void setChargeState(String chargeState)
    {
        this.chargeState = Integer.parseInt(chargeState);
    }

    public void setChargeState(int chargeState)
    {
        this.chargeState = chargeState;
    }

    public String getLoScan()
    {
        tmpStr = this.fileName.substring( this.fileName.indexOf(".") +1 );

        return tmpStr.substring(0, tmpStr.indexOf(".") );
    }

    public String getFileName()
    {
        return fileName;
    }

    public String getMS1FileName()
    {
        return fileName+".ms1";
    }


    public String getXCorr()
    {
        return xCorr;
    }

    public String getDeltCN()
    {
        return deltCN;
    }

    public String getMhPlus()
    {
        return mhPlus;
    }

    public String getSpRank()
    {
        return spRank;
    }

    public String getSpScore()
    {
        return spScore;
    }

    public String getIonProportion()
    {
        return ionProportion;
    }

    public String getRedundancy()
    {
        return redundancy;
    }

    public String getSequence()
    {
        return sequence;
    }

    /**
     * return only the sequence by trimming the '.' from both ends
     * @return
     */
    public String getSequenceOnly()
    {
        return this.getSequence().split("[.]")[1];
    }

    public String getPeptideLine()
    {
        return peptideLine;
    }

    public boolean isUnique()
    {
        return unique;
    }

    public String getCalcMHplus()
    {
        return calcMHplus;
    }

    public double getTotalIntensity()
    {
        return totalIntensity;
    }

    public int getScanNum()
    {
        return scanNum;
    }

    public void setFileName(String fileName)
    {
        this.fileName = fileName;
    }

    public void setXCorr(String xCorr)
    {
        this.xCorr = xCorr;
    }

    public void setDeltCN(String deltCN)
    {
        this.deltCN = deltCN;
    }

    public void setMhPlus(String mhPlus)
    {
        this.mhPlus = mhPlus;
    }

    public void setSpRank(String spRank)
    {
        this.spRank = spRank;
    }

    public void setSpScore(String spScore)
    {
        this.spScore = spScore;
    }

    public void setIonProportion(String ionProportion)
    {
        this.ionProportion = ionProportion;
    }

    public void setRedundancy(String redundancy)
    {
        this.redundancy = redundancy;
    }

    public void setSequence(String sequence)
    {
        this.sequence = sequence;
    }

    public void setPeptideLine(String peptideLine)
    {
        this.peptideLine = peptideLine;
    }

    public void setUnique(boolean unique)
    {
        this.unique = unique;
    }

    public void setCalcMHplus(String calcMHplus)
    {
        this.calcMHplus = calcMHplus;
    }

    public void setTotalIntensity(double totalIntensity)
    {
        this.totalIntensity = totalIntensity;
    }

    public void setScanNum(int scanNum)
    {
        this.scanNum = scanNum;
    }

    public int getChroCenter() {
        return chroCenter;
    }

    public void setChroCenter(int chroCenter) {
        this.chroCenter = chroCenter;
    }

    public String getStartRange() {
        if(null==startRange) return "0";

        return startRange;
    }

    public void setStartRange(String startRange) {
        this.startRange = startRange;
    }

    public String getEndRange() {
        if(null==endRange) return "0";

        return endRange;
    }

    public void setEndRange(String endRange) {
        this.endRange = endRange;
    }

    public int getDtaStartRange() {
        return dtaStartRange;
    }

    public void setDtaStartRange(int dtaStartRange) {
        this.dtaStartRange = dtaStartRange;
    }

    public int getDtaEndRange() {
        return dtaEndRange;
    }

    public void setDtaEndRange(int dtaEndRange) {
        this.dtaEndRange = dtaEndRange;
    }

    public double[][] getBionSample()
    {
        return bionSample;
    }

    public double[][] getBionRef()
    {
        return bionRef;
    }

    public double[][] getYionSample()
    {
        return yionSample;
    }

    public double[][] getYionRef()
    {
        return yionRef;
    }

    public void setBionSample(String str)
    {
       // String[] temp = str.split(",");

    }

    public double getSlope() {
        return slope;
    }

    public void setSlope(double slope) {
        this.slope = slope;
    }

    public double getCorr() {
        return corr;
    }

    public void setCorr(double corr) {
        this.corr = corr;
    }
    /**
     * @return the corrRev
     */
    public double getCorrRev() {
        return corrRev;
    }

    /**
     * @param corrRev the corrRev to set
     */
    public void setCorrRev(double corrRev) {
        this.corrRev = corrRev;
    }

    /**
     * @return the slopeRev
     */
    public double getSlopeRev() {
        return slopeRev;
    }

    /**
     * @param slopeRev the slopeRev to set
     */
    public void setSlopeRev(double slopeRev) {
        this.slopeRev = slopeRev;
    }

    public double getDetValue() {
	if(corr<0)
	    return 0.0;

	return corr*corr;

    }

    public boolean isFilterOut() {
        return filterOut;
    }

    public void setFilterOut(boolean filterOut) {
        this.filterOut = filterOut;
    }

    public double getProbability() {
        return probability;
    }

    public void setProbability(double probability) {
        this.probability = probability;
    }

    public double getSamIntensity() {

        if(samIntensity>0)
            return samIntensity;

        int startScan = Integer.parseInt(this.startRange);
        int endScan = Integer.parseInt(this.endRange);

        this.samIntensity = 0.0;
        for(Iterator<ChroData> itr=dataList.iterator(); itr.hasNext(); ) {
            ChroData cData = itr.next();

            cData.getScanNum();

            if(startScan<=cData.getScanNum() && endScan>=cData.getScanNum()) {
                samIntensity += cData.getSampleIntensity();
            }

        }

        return samIntensity;
    }

    public void setSamIntensity(double samIntensity) {
        this.samIntensity = samIntensity;
    }

    public double getRefIntensity() {
        return refIntensity;
    }

    public void setRefIntensity(double refIntensity) {
        this.refIntensity = refIntensity;
    }

    public double getSnRatio() {
        return snRatio;
    }

    public void setSnRatio(double snRatio) {
        this.snRatio = snRatio;
    }

    public long[] getTotalIntArr() {
        return totalIntArr;
    }

    public void setTotalIntArr(long[] totalIntArr) {
        this.totalIntArr = totalIntArr;
    }

    public double getAreaRatio() {

        if(this.refIntensity==0)
            return -1;

        return this.samIntensity/this.refIntensity;
    }

    public double getRatio() {
        return ratio;
    }

    public void setRatio(double ratio) {
        this.ratio = ratio;
    }

    public double[] getMassMonitorArr() {
        return massMonitorArr;
    }

    public void setMassMonitorArr(double[] massMonitorArr) {
        this.massMonitorArr = massMonitorArr;
    }

    public void setMassMonitorArr(String[] massMonitorArr) {
        this.massMonitorArr = new double[massMonitorArr.length];

        for(int i=0;i<massMonitorArr.length;i++)
        {
            this.massMonitorArr[i] = Double.parseDouble( massMonitorArr[i] );
        }
    }

    public double getOutlierPValue() {
        return outlierPValue;
    }

    public void setOutlierPValue(double outlierPValue) {
        this.outlierPValue = outlierPValue;
    }

    public double getCorrToNorm() {
        return corrToNorm;
    }

    public void setCorrToNorm(double corrToNorm) {
        this.corrToNorm = corrToNorm;
    }

    public int getSpectraDataPoints() {
        return spectraDataPoints;
    }

    public void setSpectraDataPoints(int spectraDataPoints) {
        this.spectraDataPoints = spectraDataPoints;
    }

    public double getMassTolerance() {
        return massTolerance;
    }

    public void setMassTolerance(double massTolerance) {
        this.massTolerance = massTolerance;
    }

    public double getAnCompositeScore() {
        return anCompositeScore;
    }

    public void setAnCompositeScore(double anCompositeScore) {
        this.anCompositeScore = anCompositeScore;
    }

    public double getFragIonNumFoundRatio() {
        return fragIonNumFoundRatio;
    }

    public void setFragIonNumFoundRatio(double fragIonNumFoundRatio) {
        this.fragIonNumFoundRatio = fragIonNumFoundRatio;
    }

    public double getAverageIntensity() {
        return averageIntensity;
    }

    public void setAverageIntensity(double averageIntensity) {
        ///System.out.println(""+averageIntensity);
        this.averageIntensity = averageIntensity;
    }

    public double getZscore() {
        return zscore;
    }

    public void setZscore(double zscore) {
        this.zscore = zscore;
    }

    public double getLightMass() {
        return lightMass;
    }

    public void setLightMass(double lightMass) {
        this.lightMass = lightMass;
    }

    public double getHeavyMass() {
        return heavyMass;
    }

    public void setHeavyMass(double heavyMass) {
        this.heavyMass = heavyMass;
    }

    public String getSpecCount() {
        return specCount;
    }

    public void setSpecCount(String specCount) {
        this.specCount = specCount;
    }

    public Hashtable<String, String> getScoreHt() {
	return scoreHt;
    }

    public void setScoreHt( Hashtable<String, String> scoreHt ) {
	this.scoreHt = scoreHt;
    }

    //public void addScore(String score, double value) {
//	this.scoreHt.put(score, value);
    //}


    public double getTmtPurity() {
        return tmtPurity;
    }

    public void setTmtPurity(double tmtPurity) {
        this.tmtPurity = tmtPurity;
    }

    public void addScore(String score, String value) {
	this.scoreHt.put(score, value);
    }

    public double getEnrichment() {
        return enrichment;
    }

    public void setEnrichment(double enrichment) {
        this.enrichment = enrichment;
    }

    public double getPvalue() {
        return pvalue;
    }

    public void setPvalue(double pvalue) {
        this.pvalue = pvalue;
    }

    public double getZvalue() {
        return zvalue;
    }

    public void setZvalue(double zvalue) {
        this.zvalue = zvalue;
    }

    public void addMultiIntensity(double intensity) {
	this.multiIntensity.add(intensity);
    }

    public ArrayList<Double> getMultiIntensity() {
	return this.multiIntensity;
    }

    public double getMultiIntensitySum() {
	double sum = 0;

	for(Iterator<Double> itr=multiIntensity.iterator(); itr.hasNext(); ) {
	    double d = itr.next();
	    sum += d;
	}

	return sum;
    }

    public double getMultiAveIntensity() {
	double sum = 0;

//	if(sequence.equals("R.FDGALNVDLTEFQTNLVPYPR.I"))

	for(Iterator<Double> itr=multiIntensity.iterator(); itr.hasNext(); ) {
	    double d = itr.next();
	    sum += d;
	}

//	System.out.println("===" + chargeState + "\t" + sequence + "\t" + multiIntensity + "\t" + sum + "\t" + multiIntensity.size());

	return sum/multiIntensity.size();
    }

    public String getDeltMass() {
        return deltMass;
    }

    public void setDeltMass(String deltMass) {
        this.deltMass = deltMass;
    }

    public String getBsText() {
        return bsText;
    }

    public void setBsText(String bsText) {
        this.bsText = bsText;
    }

    public String getYsText() {
        return ysText;
    }

    public void setYsText(String ysText) {
        this.ysText = ysText;
    }

    public String getBrText() {
        return brText;
    }

    public void setBrText(String brText) {
        this.brText = brText;
    }

    public String getYrText() {
        return yrText;
    }

    public void setYrText(String yrText) {
        this.yrText = yrText;
    }

    /**
     * @return the isoDataList
     */
    public List<IsoData> getIsoDataList() {
        return isoDataList;
    }

    /**
     * @param isoDataList the isoDataList to set
     */
    public void setIsoDataList(List<IsoData> isoDataList) {
        this.isoDataList = isoDataList;
    }

    /**
     * @return the theoryData
     */
    public TheoryData getTheoryData() {
        return theoryData;
    }

    /**
     * @param theoryData the theoryData to set
     */
    public void setTheoryData(TheoryData theoryData) {
        this.theoryData = theoryData;
        this.theoryData.generateNormalizedValues();
    }

    public List<IsoData> getIsoOrigDataList() {
        return isoOrigDataList;
    }

    public void setIsoOrigDataList(List<IsoData> isoOrigDataList) {
        this.isoOrigDataList = isoOrigDataList;
    }


    public double getBestEnrichCorr() {
        return bestEnrichCorr;
    }

    public void setBestEnrichCorr(double bestEnrichCorr) {
        this.bestEnrichCorr = bestEnrichCorr;
    }

    public double getBestEnrichDelCN() {
        return bestEnrichDelCN;
    }

    public void setBestEnrichDelCN(double bestEnrichDelCN) {
        this.bestEnrichDelCN = bestEnrichDelCN;
    }

    public double getCorrOnePlus() {
        return corrOnePlus;
    }

    public void setCorrOnePlus(double corrOnePlus) {
        this.corrOnePlus = corrOnePlus;
    }

    public double getCorrOneMinus() {
        return corrOneMinus;
    }

    public void setCorrOneMinus(double corrOneMinus) {
        this.corrOneMinus = corrOneMinus;
    }

    public JSONObject getJSONObj() {
        JSONObject obj = new JSONObject();
        obj.put("seq", sequence);
        obj.put("charge", chargeState);
        obj.put("spC", specCount);
        obj.put("unique", unique);
        obj.put("file", fileName);
        obj.put("scan", scanNum);
        obj.put("xcorr", xCorr);
        obj.put("calcMHplus", calcMHplus);
        obj.put("PROFILE_SCORE",anCompositeScore);
        obj.put("MHplus", mhPlus);
        //obj.put("totalIntensity", totalIntensity);
        obj.put("intensity", averageIntensity);
        obj.put("spRank", spRank);
        obj.put("spScore", spScore);
        obj.put("redundancy", redundancy);
        obj.put("deltaCN", deltCN);
        obj.put("deltaMass", deltMass);

//        obj.put("lightStartMass", );
//        obj.put("heavyStartMass", );
        obj.put("lightAvgMass", lightMass);
        obj.put("deltaMass", deltMass);
        obj.put("heavyAvgMass", heavyMass);
        obj.put("start_scan", startRange);
        obj.put("end_scan", endRange);
        obj.put("retentionTime", retentionTime);
        obj.put("ionInjectionTime", ionInjectionTime);
//        StringBuffer sb = new StringBuffer();
//        for (ChroData chroData : dataList) {
//            sb.append(chroData.getJSONString());
//        }
//        obj.put("chro", sb.toString());
        obj.put("missed", missedPeptide);
        if (isoTopeModelList == null) {
            obj.put("chro_iso", 0);
        }
        else{
//        if (missedPeptide){
            StringBuffer sb1 = new StringBuffer();
            StringBuffer massList = new StringBuffer();
//            double max_total = 0;
//            int max_scanNUmber = 0;
            for (IsotopeModel isoModel : isoTopeModelList) {

                sb1.append(isoModel.getIsoString() + ";");
//                if(max_total<isoModel.getIntensitySum())
//                {
//                    max_total=isoModel.getIntensitySum();
//                    max_scanNUmber = isoModel.getScanNumber();
//                }
            }
//            this.scanNum = max_scanNUmber;
            obj.put("scan", scanNum);
            if (isoTopeModelList.size() > 0) {
                obj.put("mass-List", Arrays.toString(isoTopeModelList.get(0).getIsoArr()));
            } else {
                obj.put("mass-List", "0");
            }
            obj.put("chro_iso", sb1.toString());
        }

        return obj;
    }


    public boolean isMissedPeptide() {
        return missedPeptide;
    }

    public void setMissedPeptide(boolean missedPeptide) {
        this.missedPeptide = missedPeptide;
    }

    public List<IsotopeModel> getIsoTopeModelList() {
        return isoTopeModelList;
    }

    /**
     * @param isoTopeModelList the isoTopeModelList to set
     */
    public void setIsoTopeModelList(List<IsotopeModel> isoTopeModelList) {
        this.isoTopeModelList = isoTopeModelList;
    }



    public static ChroPeptide readJsonObject(JSONObject peptideObj)
    {
        ChroPeptide peptide = new ChroPeptide();
        if(peptideObj.containsKey("unique"))
            peptide.unique = ((boolean) peptideObj.get("unique"));
        if(peptideObj.containsKey("MHplus"))
            peptide.setMhPlus((String) peptideObj.get("MHplus"));
        if(peptideObj.containsKey("calcMHplus"))
            peptide.setCalcMHplus((String) peptideObj.get("calcMHplus"));
         if(peptideObj.containsKey("intensity"))
            peptide.setAverageIntensity((double) peptideObj.get("intensity"));
        if(peptideObj.containsKey("PROFILE_SCORE"))
            peptide.setAnCompositeScore((double) peptideObj.get("PROFILE_SCORE"));
        if(peptideObj.containsKey("end_scan"))
            peptide.setEndRange(String.valueOf(peptideObj.get("end_scan")));
        if(peptideObj.containsKey("deltaCN"))
            peptide.setDeltCN((String) peptideObj.get("deltaCN"));
        if(peptideObj.containsKey("spScore"))
            peptide.setSpScore((String) peptideObj.get("spScore"));
        if(peptideObj.containsKey("deltaMass"))
            peptide.setDeltMass((String) peptideObj.get("deltaMass"));
        if(peptideObj.containsKey("redundancy"))
            peptide.setRedundancy((String) peptideObj.get("redundancy"));
        if(peptideObj.containsKey("xcorr"))
            peptide.setXCorr((String) peptideObj.get("xcorr"));
        if(peptideObj.containsKey("start_scan"))
            peptide.setStartRange(String.valueOf(peptideObj.get("start_scan")));
        if(peptideObj.containsKey("lightAvgMass"))
            peptide.setLightMass(  (double) peptideObj.get("lightAvgMass"));
        if(peptideObj.containsKey("charge"))
            peptide.setChargeState(String.valueOf(peptideObj.get("charge")));
        if(peptideObj.containsKey("heavyAvgMass"))
            peptide.setHeavyMass((double) peptideObj.get("heavyAvgMass"));

        if(peptideObj.containsKey("file"))
            peptide.setFileName((String) peptideObj.get("file"));
        if(peptideObj.containsKey("heavyAvgMass"))
            peptide.setHeavyMass((double) peptideObj.get("heavyAvgMass"));
        if(peptideObj.containsKey("seq"))
            peptide.setSequence((String) peptideObj.get("seq"));
        if(peptideObj.containsKey("scan"))
            peptide.setScanNum( Integer.parseInt(String.valueOf(peptideObj.get("scan"))));
//        if(peptideObj.containsKey("totalIntensity"))
//            peptide.setTotalIntensity((double) peptideObj.get("totalIntensity"));
        if(peptideObj.containsKey("spRank"))
            peptide.setSpRank((String) peptideObj.get("spRank"));
        if(peptideObj.containsKey("missed"))
            peptide.setMissedPeptide((boolean) peptideObj.get("missed"));
        if(peptideObj.containsKey("spC"))
            peptide.setSpecCount((String) peptideObj.get("spC"));
        if(peptideObj.containsKey("retentionTime"))
            peptide.setRetentionTime((Double) peptideObj.get("retentionTime"));
        if(peptideObj.containsKey("ionInjectionTime"))
            peptide.setIonInjectionTime((Double) peptideObj.get("ionInjectionTime"));

        if(peptideObj.containsKey("scan"))
        {
            peptide.setScanNum( Integer.parseInt(String.valueOf(peptideObj.get("scan"))));
        }
//        if(peptideObj.containsKey("chro"))
//            this.setMissedPeptide((boolean) peptideObj.get("spC"));
        IsotopeModel isoModel = new IsotopeModel();
//        if(this.missedPeptide)
//        {
            if(peptideObj.containsKey("chro_iso"))
            {
                String isoTag = ((String) peptideObj.get("chro_iso"));
                String massTag[] =((String) peptideObj.get("mass-List")).split("[\\[\\],]");
                double[] mass = new double[massTag.length-1];
                for(int i =1 ;i<massTag.length;i++)
                {
                    String value = massTag[i];
                    if(value.length()>0)
                        mass[i-1] = Double.parseDouble(value);

                }
//                List<Double> massList = new ArrayList<>();
//                for (String value : massTag)
//                {
//                    if(value.length()>0)
//                        massList.add(Double.parseDouble(value));
//                }


                for(String eachScan : isoTag.split(";"))
                {
                    isoModel = new IsotopeModel();
                    String words[] = eachScan.split(" ");
                    if(words[0].length() ==0)
                        break;
                    String subString[] = words[0].split(":");
                    isoModel.setScanNumber(Integer.parseInt(subString[0]));
                    isoModel.setIonInjectionTime(Double.parseDouble(subString[2]));
                    isoModel.setRetentionTime(Double.parseDouble(subString[1]));
                    double[] values = new double[words.length-1];
                    for(int i=1;i<words.length;i++)
                    {
                        values[i-1] = Double.parseDouble(words[i]);
                    }
                    isoModel.setIntensityArr(values);
                    isoModel.setIsoArr(mass);
                    peptide.isoTopeModelList.add(isoModel);
                }
            }
//        }
            return peptide;
    }

    public void setDataList(List<ChroData> dataList) {
        this.dataList = dataList;
    }


    public double getRetentionTime() {
        return retentionTime;
    }

    public void setRetentionTime(double retentionTime) {
        this.retentionTime = retentionTime;
    }

    public double getIonInjectionTime() {
        return ionInjectionTime;
    }

    public void setIonInjectionTime(double ionInjectionTime) {
        this.ionInjectionTime = ionInjectionTime;
    }

    public double getIonInjectionTimeNormIntensity() {
        return ionInjectionTimeNormIntensity;
    }

    public void setIonInjectionTimeNormIntensity(double ionInjectionTimeNormIntensity) {
        this.ionInjectionTimeNormIntensity = ionInjectionTimeNormIntensity;
    }

    public List<LinearRegression> getRegList() {
        return regList;
    }

    public void setRegList(List<LinearRegression> regList) {
        this.regList = regList;
    }

    public double[] getIntensitySum() {
        return intensitySum;
    }

    public void setIntensitySum(double[] intensitySum) {
        this.intensitySum = intensitySum;
    }

    public void setPrecursorLightIntensity(long precursorLightIntensity) {
        this.precursorLightIntensity = precursorLightIntensity;
    }

    public long getPrecursorHeavyIntensity() {
        return precursorHeavyIntensity;
    }

    public void setPrecursorHeavyIntensity(long precursorHeavyIntensity) {
        this.precursorHeavyIntensity = precursorHeavyIntensity;
    }

    public long getGaussianPeakArea() {
        return gaussianPeakArea;
    }

    public void setGaussianPeakArea(long gaussianPeakArea) {
        this.gaussianPeakArea = gaussianPeakArea;
    }

    public double getPeakArea() {
        return this.averageIntensity;
    }

    public void setPeakArea(double averageIntensity) {
        this.averageIntensity = averageIntensity;
    }

    public String getChroData() {
        return chroData;
    }

    public void setChroData(String chroData) {
        this.chroData = chroData;
    }

    public double getPeakSigma() {
        return peakSigma;
    }

    public void setPeakSigma(double peakSigma) {
        this.peakSigma = peakSigma;
    }

    public double getPeakx() {
        return peakx;
    }

    public void setPeakx(double peakx) {
        this.peakx = peakx;
    }

    public double getPeaky() {
        return peaky;
    }

    public void setPeaky(double peaky) {
        this.peaky = peaky;
    }

    public double getStartRt() {
        return startRt;
    }

    public void setStartRt(double startRt) {
        this.startRt = startRt;
    }

    public double getEndRt() {
        return endRt;
    }

    public void setEndRt(double endRt) {
        this.endRt = endRt;
    }

    public String getGaussianPeakString() {
        return gaussianPeakString;
    }

    public void setGaussianPeakString(String gaussianPeakString) {
        this.gaussianPeakString = gaussianPeakString;
    }

    public boolean isSingleton() {
        return singleton;
    }

    public void setSingleton(boolean singleton) {
        this.singleton = singleton;
    }

    public double[] getIsoArr() {
        return isoArr;
    }

    public void setIsoArr(double[] isoArr) {
        this.isoArr = isoArr;
    }

    public int getPeptideIndex() {
        return peptideIndex;
    }

    public void setPeptideIndex(int peptideIndex) {
        this.peptideIndex = peptideIndex;
    }

    public int getSampleIndex() {
        return sampleIndex;
    }

    public void setSampleIndex(int sampleIndex) {
        this.sampleIndex = sampleIndex;
    }

  public double getProfileScore() {
    return profileScore;
  }

  public void setProfileScore(double profileScore) {
    this.profileScore = profileScore;
  }

  public List<GaussianPeakModel> getGaussianPeakModelList() {
    return gaussianPeakModelList;
  }

  public void setGaussianPeakModelList(List<GaussianPeakModel> gaussianPeakModelList) {
    this.gaussianPeakModelList = gaussianPeakModelList;
  }

  public void addGaussianPeakModel(GaussianPeakModel gaussianPeakModel) {
    this.gaussianPeakModelList.add(gaussianPeakModel);
  }


  public double[] getDiaFragRegressionArr() {
    return diaFragRegressionArr;
  }

  public void setDiaFragRegressionArr(double[] diaFragRegressionArr) {
    this.diaFragRegressionArr = diaFragRegressionArr;
  }

    public String getPtmIndex() {
        return ptmIndex;
    }

    public void setPtmIndex(String ptmIndex) {
        this.ptmIndex = ptmIndex;
    }

    public String getPtmIndexProtein() {
        return ptmIndexProtein;
    }

    public void setPtmIndexProtein(String ptmIndexProtein) {
        this.ptmIndexProtein = ptmIndexProtein;
    }

    public LuciphorReader.LuciphorPeptide getLuciphorPeptide() {
        return luciphorPeptide;
    }

    public void setLuciphorPeptide(LuciphorReader.LuciphorPeptide luciphorPeptide) {
        this.luciphorPeptide = luciphorPeptide;
    }
}



