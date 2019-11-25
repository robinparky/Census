package edu.scripps.pms.census.util.dtaselect;

import java.util.*;

/**
 * @author  Robin Park
 * @version $Id: Peptide.java,v 1.1 2014/09/09 19:28:27 rpark Exp $
 */

import java.util.*;

public class Peptide
{

    private boolean unique;
    private String fileName;
    private String xCorr;
    private String deltCN;
    private String mhPlus;
    private String calcMHplus;
    private String totalIntensity;
    private String spRank;
    private String spScore;
    private String ionProportion;
    private String redundancy;
    private String sequence;
    private String[] peptideLine;
    private String tmpStr;
    private String scanNum;
    private String conf;
    private String kd;

    private String ptmIndex="NA";
    private String ptmIndexProtein="NA";

    private static int uniqueIndex = 0;
    private static int scanNumIndex = 1;
    private static int xcorrIndex = 2;
    private static int dcnIndex = 3;
    private static int confIndex = 4;
    private static int mPlusHIndex = 5;
    private static int calcMassIndex = 6;
    private static int totalIntensityIndex = 7;
    private static int spRankIndex = 8;
    private static int spScoreIndex = 9;
    private static int ionProportionIndex = 10;
    private static int redundancyIndex = 11;
    private static int sequenceIndex = 12;
    private static int pIIndex = 9; 
    private static int ppmIndex = -1; 
    private static int kdIndex = -1; 
    private int hashcode = -1;

    private static int ptmIndexIndex = -1;
    private static int ptmIndexProteinIndex = -1;

    private String filePath;
    private Hashtable<String, String> scoreHt = new Hashtable<String, String>();
    
    private String chargeState = null;
    private double[][] bionSample;
    private double[][] bionRef;
    private double[][] yionSample;
    private double[][] yionRef;    

    public Peptide() {
    }

    public int hashCode() {
        if(hashcode == -1) {
            // fileName contains the scan number
            hashcode = (getSequence() + fileName).hashCode();
            //hashcode = (fileName + scanNum).hashCode();
            //System.out.println("file: " + fileName + "\tscannumber: " + scanNum);
        } 
//System.out.println(hashcode);
        return hashcode;
    } 
    public boolean equals(Object o) {
        Peptide p = (Peptide)o;
//System.out.println("in equals");
        return getSequence().equals(p.getSequence()) && fileName.equals(p.fileName);
    }
    public String getInfo() {
        StringBuffer sb = new StringBuffer(1000);
        sb.append(sequence + "\t" + xCorr + "\t" + deltCN);
        sb.append("\t" + mhPlus + "\t" + calcMHplus + "\t" + getDeltaMass());
        sb.append("\t" + spScore + "\t" + conf + "\t" + fileName);
        return sb.toString();
    }
    /* DTASelect 2.0 file */
    private void parseLine2() throws ArrayIndexOutOfBoundsException
    {

        scanNum = peptideLine[scanNumIndex];
        
      //  System.out.println("confIndex====: " + scanNumIndex + " " + confIndex);
        
        scanNum = scanNum.substring(0, scanNum.lastIndexOf("."));
        scanNum = scanNum.substring(scanNum.lastIndexOf(".")+1);
//System.out.println("scanNum: " + scanNum);
        this.setUnique((peptideLine[uniqueIndex]).startsWith("*"));
        this.setFileName(peptideLine[scanNumIndex]);
        this.setXCorr(peptideLine[xcorrIndex]);
        this.setDeltCN(peptideLine[dcnIndex]);
        this.setConf(peptideLine[confIndex]);
        this.setMhPlus(peptideLine[mPlusHIndex]);
        this.setCalcMHplus(peptideLine[calcMassIndex]);
        this.setTotalIntensity(peptideLine[totalIntensityIndex]);
        this.setSpRank(peptideLine[spRankIndex]);
        this.setSpScore(peptideLine[spScoreIndex]);
        this.setIonProportion(peptideLine[ionProportionIndex]);
        this.setRedundancy(peptideLine[redundancyIndex]);
        this.setSequence(peptideLine[sequenceIndex]);

        //robin add here


        if(this.ptmIndexIndex>0 && peptideLine.length>this.ptmIndexIndex) {
            this.setPtmIndex(peptideLine[ptmIndexIndex]);
        }

        if(this.ptmIndexProteinIndex>=0 && peptideLine.length>ptmIndexProteinIndex) {
            this.setPtmIndexProtein(peptideLine[ptmIndexProteinIndex]);
        }



        if(kdIndex>0 && peptideLine.length>kdIndex)
	        this.setKd(peptideLine[kdIndex]);
//System.out.println("\t11: " + peptideLine[11] + "\t12: " + peptideLine[12]+ "\t13: " + peptideLine[13] + "\t14: " + peptideLine[14]);
    }
    
    public static void setFeatureIndices(String features) {

        String [] contents = features.split("\t");
        for(int i = 0; i < contents.length; i++) {
            String s = contents[i].trim();
            uniqueIndex = s.startsWith("Uni")? i : uniqueIndex; 
            scanNumIndex = s.startsWith("File")? i : scanNumIndex;
            xcorrIndex = s.startsWith("XC")? i : xcorrIndex;
            dcnIndex = s.startsWith("DeltCN")? i : dcnIndex;
            confIndex = s.startsWith("Conf%")? i : confIndex;
            mPlusHIndex = s.startsWith("M")? i : mPlusHIndex;
            calcMassIndex = s.startsWith("CalcM")? i : calcMassIndex;
            totalIntensityIndex = s.startsWith("Total")? i : totalIntensityIndex;
            spRankIndex = s.startsWith("SpR")? i : spRankIndex;
            spScoreIndex = s.startsWith("Prob")? i : spScoreIndex;
            ionProportionIndex = s.startsWith("IonP")? i : ionProportionIndex;
            redundancyIndex = s.startsWith("Red")? i : redundancyIndex;
            sequenceIndex = s.startsWith("Seq")? i : sequenceIndex;
            pIIndex = s.startsWith("pI")? i : pIIndex;            
            ppmIndex = s.startsWith("PPM")? i : ppmIndex;            
            kdIndex = s.startsWith("KD")? i : kdIndex;


            ptmIndexIndex = s.equals("PTMIndex")? i :ptmIndexIndex;
            ptmIndexProteinIndex = s.equals("PTMIndex Protein List")? i:ptmIndexProteinIndex;
            
        }

    } 
    //For DTASelect version 2
    public Peptide(String[] peptideLine, boolean isV2)
    {
        this.peptideLine = peptideLine;
        
        //if(isV2)
            parseLine2();
        //else
            //parseLine();
    }

    
    /* DTASelect file */
    private void parseLine() throws ArrayIndexOutOfBoundsException
    {
        scanNum = peptideLine[1];
        scanNum = scanNum.substring(0, scanNum.lastIndexOf("."));
        scanNum = scanNum.substring(scanNum.lastIndexOf(".")+1);

        //this.setUnique(!"".equals(peptideLine[0]));
        this.setUnique((peptideLine[0]).startsWith("*"));
        this.setFileName(peptideLine[1]);
        this.setXCorr(peptideLine[2]);
        this.setDeltCN(peptideLine[3]);
        this.setMhPlus(peptideLine[4]);
        this.setCalcMHplus(peptideLine[5]);
        this.setTotalIntensity(peptideLine[6]);
        this.setSpRank(peptideLine[7]);
        this.setSpScore(peptideLine[8]);
        this.setIonProportion(peptideLine[9]);
        this.setRedundancy(peptideLine[10]);
        this.setSequence(peptideLine[11]);
    }
    public void setDTASelectTxtPeptideLine(String peptideLine) {
	String[] arr = peptideLine.split("\t");

        scanNum = arr[1];
        scanNum = scanNum.substring(0, scanNum.lastIndexOf("."));
        scanNum = scanNum.substring(scanNum.lastIndexOf(".")+1);

        this.setUnique((arr[arr.length-1]).startsWith("U"));
        this.setFileName(arr[1]);
        this.setXCorr(arr[3]);
        this.setDeltCN(arr[4]);
        this.setConf( arr[5] );
        this.setMhPlus(arr[6]);
        this.setCalcMHplus(arr[7]);
        this.setTotalIntensity(arr[8]);
        this.setSpRank(arr[9]);
        this.setSpScore(arr[10]);
        this.setIonProportion( arr[10] );
        this.setRedundancy("-1"); //no value
        this.setSequence( arr[12] );
    }
    public float getPi() {
        return Float.parseFloat(peptideLine[pIIndex]);
    }
    public float getDeltaMass() {
        
        return ppmIndex != -1 ? Float.parseFloat(peptideLine[ppmIndex]) : 1000;
    }

    public String getLoScan()
    {
        tmpStr = this.fileName.substring( this.fileName.indexOf(".") +1 );

        return tmpStr.substring(0, tmpStr.indexOf(".") );
    }

    public String getFileName()
    {
        return fileName.substring(0, fileName.indexOf("."));
    }

    public String getFileNameWithScan() {
        return fileName;
    }

    public String getXCorr()
    {
        return xCorr;
    }

    public double  getXCorrValue()
    {
        return Double.parseDouble(xCorr);
    }
    public String getDeltCN()
    {
        return deltCN;
    }
    public double  getDeltCnValue()
    {
        return Double.parseDouble(deltCN);
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
    public double getSpScoreValue()
    {
        return Double.parseDouble(spScore);
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
    // return the peptide sequence without leading and tailing residues
    public String getMidSeq()
    {
        int lastindex = sequence.length() - 2;
        return sequence.substring(2, lastindex);
    }

    public boolean isUnique()
    {
        return unique;
    }

    public String getCalcMHplus()
    {
        return calcMHplus;
    }

    public String getTotalIntensity()
    {
        return totalIntensity;
    }

    public String getScanNum()
    {
        return scanNum;
    }
    public int getScanNumber()
    {
        return Integer.parseInt(scanNum);
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

    public void setUnique(boolean unique)
    {
        this.unique = unique;
    }

    public void setCalcMHplus(String calcMHplus)
    {
        this.calcMHplus = calcMHplus;
    }

    public void setTotalIntensity(String totalIntensity)
    {
        this.totalIntensity = totalIntensity;
    }

    public void setScanNum(String scanNum)
    {
        this.scanNum = scanNum;
    }

    public String getConf() {
        return (null==conf)?"":conf;
    }
    public double  getConfValue() {
        if(conf == null) {
            return 0;
        }
        return Double.parseDouble(conf);
    }

    public void setConf(String conf) {
        this.conf = conf;
    }

    public String[] getPeptideLine()
    {
        return peptideLine;
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

    public void addScore(String score, String value) {
	this.scoreHt.put(score, value);
    }

    public String getFilePath() {
        return filePath;
    }

    public void setFilePath(String filePath) {
        this.filePath = filePath;
    }

    public void setChargeState(String chargeState) {
        this.chargeState = chargeState;
    }
    public String getChargeState()
    {
        if(null != this.chargeState || "".equals(this.chargeState))
            return chargeState;
        
        return this.fileName.substring( this.fileName.lastIndexOf(".") + 1 );
    }

    public double[][] getBionSample() {
        return bionSample;
    }

    public void setBionSample(double[][] bionSample) {
        this.bionSample = bionSample;
    }

    public double[][] getBionRef() {
        return bionRef;
    }

    public void setBionRef(double[][] bionRef) {
        this.bionRef = bionRef;
    }

    public double[][] getYionSample() {
        return yionSample;
    }

    public void setYionSample(double[][] yionSample) {
        this.yionSample = yionSample;
    }

    public double[][] getYionRef() {
        return yionRef;
    }

    public void setYionRef(double[][] yionRef) {
        this.yionRef = yionRef;
    }

    public void setKd(String kd) {
	this.kd = kd;
    }

    public String getKd() {
	return this.kd;
    }

    public int getSpectralCount()
    {
        return Integer.parseInt(redundancy);
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


    public static int getPtmIndexIndex() {
        return ptmIndexIndex;
    }

    public static void setPtmIndexIndex(int ptmIndexIndex) {
        Peptide.ptmIndexIndex = ptmIndexIndex;
    }

    public static int getPtmIndexProteinIndex() {
        return ptmIndexProteinIndex;
    }

    public static void setPtmIndexProteinIndex(int ptmIndexProteinIndex) {
        Peptide.ptmIndexProteinIndex = ptmIndexProteinIndex;
    }
}
