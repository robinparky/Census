/*
 * ChroData.java
 *
 * Created on May 31, 2005, 2:10 PM
 */

package edu.scripps.pms.census.model;

/**
 *
 * @author rpark
 */
public class ChroData {
    
    /** Creates a new instance of ChroData */
    
    private int scanNum;
    private long sampleIntensity;
    private long refIntensity;
    private long refIntensity2;

    private long[] bsIntensity = null;
    private long[] ysIntensity = null;
    private long[] brIntensity = null;
    private long[] yrIntensity = null;

    private double[] bsStartMass = null;
    private double[] bsEndMass = null;
    private double[] ysStartMass = null;
    private double[] ysEndMass = null;
    
    private long[] intensityArr; //for iTRAQ
    
    private String bsIntensityString = null;
    private String ysIntensityString = null;
    private String brIntensityString = null;
    private String yrIntensityString = null;
    private double massToleranceLight;
    private double massToleranceHeavy;
    private double massToleranceHeavy2;
    
    private int totalIsoPeaksSam;
    private int foundIsoPeaksSam;
    private int totalIsoPeaksRef;
    private int foundIsoPeaksRef;
    private int totalIsoPeaksRef2;
    private int foundIsoPeaksRef2;
    
    private double retentionTime; 
    private double[] dIntensityArr;
 
    
    public ChroData() {        
    }
    
    public ChroData(String[] arr, String quantType) {

        if("3plexMS1Labeling".equals(quantType)) {  //this is outdated.  we started using "MultipleMs1Labeling"

            this.scanNum = Integer.parseInt(arr[0]);
            this.sampleIntensity = (long)Double.parseDouble(arr[1]);
            this.refIntensity = (long)Double.parseDouble(arr[2]);
            this.refIntensity2 = (long)Double.parseDouble(arr[3]);
            this.massToleranceLight = Double.parseDouble(arr[10]);
            this.massToleranceHeavy = Double.parseDouble(arr[11]);
            this.massToleranceHeavy2 = Double.parseDouble(arr[12]);
            this.totalIsoPeaksSam = Integer.parseInt(arr[4]);
            this.foundIsoPeaksSam = Integer.parseInt(arr[7]);
            this.totalIsoPeaksRef = Integer.parseInt(arr[5]);
            this.foundIsoPeaksRef = Integer.parseInt(arr[8]);
            this.totalIsoPeaksRef = Integer.parseInt(arr[6]);
            this.foundIsoPeaksRef = Integer.parseInt(arr[9]);
        } else if ("MultipleMs1Labeling".equals(quantType)) {
            this.scanNum = Integer.parseInt(arr[0]);
            this.retentionTime = Double.parseDouble(arr[1]);

            this.dIntensityArr = new double[arr.length-2];

            for(int i=2;i<arr.length;i++) {
                this.dIntensityArr[i-2] = Double.parseDouble(arr[i]);
            }

           // System.out.println("aaaaaaaaaaaaaa");

        }

    }

    public ChroData(int scanNum, long sampleIntensity, long refIntensity, double mtLight, double mtHeavy, int totalIPSam, int foundIPSam, int totalIPRef, int foundIPRef) {
        this(scanNum, sampleIntensity, refIntensity, mtLight, mtHeavy);
        
        this.totalIsoPeaksSam = totalIPSam;
        this.foundIsoPeaksSam = foundIPSam;
        this.totalIsoPeaksRef = totalIPRef;
        this.foundIsoPeaksRef = foundIPRef;
        
    }


    public ChroData(int scanNum, long sampleIntensity, long refIntensity, double mtLight, double mtHeavy) {
        this(scanNum, sampleIntensity, refIntensity);
        
        this.massToleranceLight = mtLight;
        this.massToleranceHeavy = mtHeavy;
    }


    public ChroData(int scanNum, long sampleIntensity, long refIntensity) {
        this.scanNum = scanNum;
        this.sampleIntensity = sampleIntensity;
        this.refIntensity = refIntensity;
    }

    public int getScanNum() {
        return scanNum;
    }

    public void setFullScanData(String[] arr) {}
    public void setTandemData(String[] arr) {}
    
    public void setScanNum(int scanNum) {
        this.scanNum = scanNum;
    }

    public long getSampleIntensity() {
        return sampleIntensity;
    }

    public void setSampleIntensity(long sampleIntensity) {
        this.sampleIntensity = sampleIntensity;
    }

    public long getRefIntensity() {
        return refIntensity;
    }

    public void setRefIntensity(long refIntensity) {
        this.refIntensity = refIntensity;
    }
    
    public long getRefIntensity2() {
        return refIntensity2;
    }

    public void setRefIntensity2(long refIntensity2) {
        this.refIntensity2 = refIntensity2;
    }

    public void setBsIntensity(String[] arr) {
        bsIntensity = new long[arr.length];

        for(int i=0;i<arr.length;i++)
        {
            bsIntensity[i] = (long)Double.parseDouble(arr[i]);
        }
    }

    public void setYsIntensityReverse(String[] arr) {
        ysIntensity = new long[arr.length];

        for(int i=0;i<arr.length;i++)
        {
            ysIntensity[arr.length-i-1] = (long)Double.parseDouble(arr[i]);
        }
    }
    
    public void setYsIntensity(String[] arr) {
        ysIntensity = new long[arr.length];

        for(int i=0;i<arr.length;i++)
        {
            ysIntensity[i] = (long)Double.parseDouble(arr[i]);
        }
    }

    public void setBrIntensity(String[] arr) {
        brIntensity = new long[arr.length];

        for(int i=0;i<arr.length;i++)
        {
            brIntensity[i] = (long)Double.parseDouble(arr[i]);
        }
    }

    public void setYrIntensity(String[] arr) {
        yrIntensity = new long[arr.length];

        for(int i=0;i<arr.length;i++)
        {
            yrIntensity[i] = (long)Double.parseDouble(arr[i]);
        }
    }

    private long[] parseIntensity(String intStr)
    {
	String[] arr = intStr.split(" ");
	long[] longArr = new long[arr.length];

	for(int i=0;i<arr.length;i++)
	{
	    longArr[i] = (long)Double.parseDouble(arr[i]);
	}

	return longArr;
    }

    public long[] getBsIntensity() {
	if(null != bsIntensity)
	    return bsIntensity;

	return parseIntensity(bsIntensityString);
    }

    public long[] getYsIntensity() {
	if(null != ysIntensity)
	    return ysIntensity;

	return parseIntensity(ysIntensityString);
    }

    public long[] getBrIntensity() {
	if(null != brIntensity)
	    return brIntensity;

	return parseIntensity(brIntensityString);
    }

    public long[] getYrIntensity() {
	if(null != yrIntensity)
	    return yrIntensity;

	return parseIntensity(yrIntensityString);
    }

    public int getResidueLength()
    {
	if(null != bsIntensity)
	    return bsIntensity.length;

	String[] arr = bsIntensityString.split(" ");
	return arr.length;
    }

    public double[] getBsStartMass() {
        return bsStartMass;
    }

    public void setBsStartMass(double[] bsStartMass) {
        this.bsStartMass = bsStartMass;
    }

    public double[] getBsEndMass() {
        return bsEndMass;
    }

    public void setBsEndMass(double[] bsEndMass) {
        this.bsEndMass = bsEndMass;
    }

    public double[] getYsStartMass() {
        return ysStartMass;
    }

    public void setYsStartMass(double[] ysStartMass) {
        this.ysStartMass = ysStartMass;
    }

    public double[] getYsEndMass() {
        return ysEndMass;
    }

    public void setYsEndMass(double[] ysEndMass) {
        this.ysEndMass = ysEndMass;
    }

    public long[] getIntensityArr() {
        return intensityArr;
    }

    public void setIntensityArr(long[] intensityArr) {
        this.intensityArr = intensityArr;
    }

    public long getIntensityArrTotal() {
	long sum =0;
	for(long l:intensityArr)
	    sum+=l;

	return sum;
    }

    public String getBsIntensityString() {
        return bsIntensityString;
    }

    public void setBsIntensity(String bsIntensityString) {
        this.bsIntensityString = bsIntensityString;
    }

    public String getYsIntensityString() {
        return ysIntensityString;
    }

    public void setYsIntensity(String ysIntensityString) {
        this.ysIntensityString = ysIntensityString;
    }

    public String getBrIntensityString() {
        return brIntensityString;
    }

    public void setBrIntensity(String brIntensityString) {
        this.brIntensityString = brIntensityString;
    }

    public String getYrIntensityString() {
        return yrIntensityString;
    }

    public void setYrIntensity(String yrIntensityString) {
        this.yrIntensityString = yrIntensityString;
    }

    public double getMassToleranceLight() {
        return massToleranceLight;
    }

    public void setMassToleranceLight(double massToleranceLight) {
        this.massToleranceLight = massToleranceLight;
    }

    public double getMassToleranceHeavy() {
        return massToleranceHeavy;
    }

    public void setMassToleranceHeavy(double massToleranceHeavy) {
        this.massToleranceHeavy = massToleranceHeavy;
    }

    public double getMassToleranceHeavy2() {
        return massToleranceHeavy2;
    }

    public void setMassToleranceHeavy2(double massToleranceHeavy2) {
        this.massToleranceHeavy2 = massToleranceHeavy2;
    }


    public int getTotalIsoPeaksSam() {
        return totalIsoPeaksSam;
    }

    public void setTotalIsoPeaksSam(int totalIsoPeaksSam) {
        this.totalIsoPeaksSam = totalIsoPeaksSam;
    }

    public int getFoundIsoPeaksSam() {
        return foundIsoPeaksSam;
    }

    public void setFoundIsoPeaksSam(int foundIsoPeaksSam) {
        this.foundIsoPeaksSam = foundIsoPeaksSam;
    }

    public int getTotalIsoPeaksRef() {
        return totalIsoPeaksRef;
    }

    public void setTotalIsoPeaksRef(int totalIsoPeaksRef) {
        this.totalIsoPeaksRef = totalIsoPeaksRef;
    }

    public int getFoundIsoPeaksRef() {
        return foundIsoPeaksRef;
    }

    public void setFoundIsoPeaksRef(int foundIsoPeaksRef) {
        this.foundIsoPeaksRef = foundIsoPeaksRef;
    }
    
    public int getTotalIsoPeaksRef2() {
        return totalIsoPeaksRef2;
    }

    public void setTotalIsoPeaksRef2(int totalIsoPeaksRef2) {
        this.totalIsoPeaksRef2 = totalIsoPeaksRef2;
    }

    public int getFoundIsoPeaksRef2() {
        return foundIsoPeaksRef2;
    }

    public void setFoundIsoPeaksRef2(int foundIsoPeaksRef2) {
        this.foundIsoPeaksRef2 = foundIsoPeaksRef2;
    } 
    
    public String getJSONString()
    {
//        JSONObject obj = new JSONObject();
       StringBuffer sb = new StringBuffer();

       sb.append(scanNum).append(" ");
       sb.append(sampleIntensity).append(" ");
       sb.append(refIntensity).append(" ");
       sb.append(massToleranceLight).append(" ");
       sb.append(massToleranceHeavy).append(" ");
       sb.append(totalIsoPeaksSam).append(" ");
       sb.append(foundIsoPeaksSam).append(" ");
       sb.append(totalIsoPeaksRef).append(" ");
       sb.append(foundIsoPeaksRef).append(";");


        return sb.toString();
    }

    public double getRetentionTime() {
        return retentionTime;
    }

    public void setRetentionTime(double retentionTime) {
        this.retentionTime = retentionTime;
    }

    public double[] getdIntensityArr() {
        return dIntensityArr;
    }

    public void setdIntensityArr(double[] dIntensityArr) {
        this.dIntensityArr = dIntensityArr;
    }

    
    
}
