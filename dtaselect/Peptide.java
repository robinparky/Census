package edu.scripps.pms.util.dtaselect;

/**
 * @author  Robin Park
 * @version $Id: Peptide.java,v 1.2 2006/10/13 05:50:30 rpark Exp $
 */
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
    
    private String filePath;
    

    /*
    public Peptide(String[] peptideLine)
    {
        this.peptideLine = peptideLine;
        parseLine();
    }
    */
    
    //For DTASelect version 2
    public Peptide(String[] peptideLine, boolean isV2)
    {
        this.peptideLine = peptideLine;
        
        if(isV2)
            parseLine2();
        else
            parseLine();
    }

    /* DTASelect 2.0 file */
    private void parseLine2() throws ArrayIndexOutOfBoundsException
    {
        scanNum = peptideLine[1];
        scanNum = scanNum.substring(0, scanNum.lastIndexOf("."));
        scanNum = scanNum.substring(scanNum.lastIndexOf(".")+1);

        this.setUnique((peptideLine[0]).startsWith("*"));
        this.setFileName(peptideLine[1]);
        this.setXCorr(peptideLine[2]);
        this.setDeltCN(peptideLine[3]);
        this.setConf(peptideLine[4]);
        this.setMhPlus(peptideLine[5]);
        this.setCalcMHplus(peptideLine[6]);
        this.setTotalIntensity(peptideLine[7]);
        this.setSpRank(peptideLine[8]);
        this.setSpScore(peptideLine[9]);
        this.setIonProportion(peptideLine[10]);
        this.setRedundancy(peptideLine[11]);
        this.setSequence(peptideLine[12]);
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

    public String getChargeState()
    {
        return this.fileName.substring( this.fileName.lastIndexOf(".") + 1 );
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

    public void setConf(String conf) {
        this.conf = conf;
    }

    public String[] getPeptideLine()
    {
        return peptideLine;
    }

    public String getFilePath() {
        return filePath;
    }

    public void setFilePath(String filePath) {
        this.filePath = filePath;
    }
}
