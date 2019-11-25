/*
 * ChroProtein.java
 *
 * Created on May 23, 2005, 11:38 AM
 */

package edu.scripps.pms.census.model;

import edu.scripps.pms.census.util.dtaselect.Peptide;
import edu.scripps.pms.census.util.dtaselect.Protein;
import java.util.*;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;

/**
 *
 * @author rpark
 * @version $Id: ChroProtein.java,v 1.9 2014/08/27 18:00:35 rpark Exp $
 */
public class ChroProtein {
    private static int locusIndex = 0;
    private static int seqCountIndex = 1;
    private static int spectrumCountIndex = 2;
    private static int seqCoverageIndex = 3;
    private static int lengthIndex = 4;
    private static int molWtIndex = 5;
    private static int pIIndex = 6;
    private static int validationIndex = 7;
    private static int descriptionIndex = 8;

    private String locus;
    private String seqCount;
    private String spectrumCount;
    private String seqCoverage;
    private String length;
    private String molWt;
    private String pI;
    private String validation;
    private String description;
    private List<ChroPeptide> peptideList = new ArrayList<ChroPeptide>();
    private String proteinLine;
    private String[] strArr=null;
    private boolean redundant=false;

    private Double averageRatio;
    private Double intensity;

    private List<ChroProtein> redunList = new ArrayList<ChroProtein>();
    private String lspectrumCount="NA";
    private String hspectrumCount="NA";
    private int pepCount = -1;
    private List<Double> pepRatioList = new ArrayList<>();
    private List<Integer> specCountList = new ArrayList<>();
    private List<Double> avgIntensityList = new ArrayList<>();// used for labelfree final output.. average of peptide for a protein.
    public void addRedunProtein(ChroProtein protein) {
	redunList.add(protein);
    }

    public List<ChroProtein> getRedunList() {
	return redunList;
    }

    public static ChroProtein convertProtein(Protein p)
    {
        ChroProtein cp = new ChroProtein();
        cp.setLocus(p.getLocus());
        cp.setSeqCount(p.getSeqCount());
        cp.setSpectrumCount(p.getSpectrumCount());
        cp.setSeqCoverage(p.getSeqCoverage());
        cp.setLength(p.getLength());
        cp.setMolWt(p.getMolWt());
        cp.setPI(p.getPI());
        cp.setValidation(p.getValidation());
        cp.setDescription(p.getDescription());
        cp.setProteinLine(p.getProteinLine());

        for(Iterator<Peptide> itr=p.getPeptides(); itr.hasNext(); )
        {
            Peptide pep = itr.next();
            ChroPeptide chPep = new ChroPeptide();
            chPep.setUnique( pep.isUnique() );
            chPep.setFileName( pep.getFileName() );
            chPep.setXCorr( pep.getXCorr() );
            chPep.setDeltCN( pep.getDeltCN() );
            chPep.setMhPlus( pep.getMhPlus() );
            chPep.setCalcMHplus( pep.getCalcMHplus() );
            chPep.setTotalIntensity( Double.parseDouble(pep.getTotalIntensity()) );
            chPep.setSpRank( pep.getSpRank() );
            chPep.setSpScore( pep.getSpScore() );
            chPep.setIonProportion( pep.getIonProportion() );
            chPep.setRedundancy( pep.getRedundancy() );
            chPep.setSequence( pep.getSequence() );
            chPep.setScanNum( Integer.parseInt(pep.getScanNum()) );

            cp.addPeptide( chPep );
        }

        return cp;
    }

    public ChroProtein()
    {
    }

    public ChroProtein(String locus, String spectrumCount, String description, String pLine)
    {
        this.locus = locus;
        this.spectrumCount = spectrumCount;
        this.description = description;
        this.proteinLine= pLine;

        //System.out.println("----------" + this.description + " " + this.spectrumCount);
    }

    public ChroProtein(String proteinLine) throws ArrayIndexOutOfBoundsException
    {
	this( proteinLine.split("\t") );
	this.setProteinLine(proteinLine);

    }

    public ChroProtein(String[] strArr) throws ArrayIndexOutOfBoundsException
    {
	try {
        this.strArr = strArr; //new String[] {strArr[0], strArr[1], strArr[2], strArr[3], strArr[4], strArr[5], strArr[6], strArr[8] };

        init(strArr);
	} catch (ArrayIndexOutOfBoundsException ae) {
	    System.out.println("e" + ae);
	} catch (Exception e) {
	    System.out.println("e" + e);
	}
    }

    public boolean addPeptide(ChroPeptide peptide)
    {
        if(null == peptide)
            return false;

        //System.out.println("." + peptide);
        return peptideList.add(peptide);
    }

    public void deletePeptides()
    {
        this.peptideList.clear();
    }

    public List getPeptideList()
    {
        return peptideList;
    }

    public void setPeptideList(List<ChroPeptide> peptideList)
    {
        this.peptideList = peptideList;
    }

    private void init(String[] strArr) throws ArrayIndexOutOfBoundsException
    {
        try {
            this.setLocus(strArr[locusIndex]);
            this.setSeqCount(strArr[seqCountIndex]);
            this.setSpectrumCount(strArr[spectrumCountIndex]);
            this.setSeqCoverage(strArr[seqCoverageIndex]);
            this.setLength(strArr[lengthIndex]);
            this.setMolWt(strArr[molWtIndex]);
            this.setPI(strArr[pIIndex]);
            this.setValidation(strArr[validationIndex]);
            this.setDescription(strArr[descriptionIndex]);
        }
        catch(ArrayIndexOutOfBoundsException ex) {
            throw new ArrayIndexOutOfBoundsException("Error : Mal formed DTASelect-filter.txt file in protein line");
        }
    }

/*
    public static void setIndex(String str) {
	String [] arr = str.split("\t");
	for(int i = 2; i < arr.length; i++) {
	    int index = i-2;

PLINE   LOCUS   AVERAGE_RATIO   STANDARD_DEVIATION      WEIGHTED_AVERAGE        PEPTIDE_NUM     SPEC_COUNT      DESCRIPTION

	    if(s.startsWith("LOCUS")) {
		locusIndex = index;
	    } else if(s.startsWith("")) {

	    } else if(s.startsWith("")) {

	    } else if(s.startsWith("")) {
	    } else if(s.startsWith("")) {
	    } else if(s.startsWith("")) {
    seqCountIndex = 1;
    spectrumCountIndex = 2;
    seqCoverageIndex = 3;
    lengthIndex = 4;
    molWtIndex = 5;
    pIIndex = 6;
    validationIndex = 7;
    descriptionIndex = 8;
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

	}

    }
*/
    public String[] getProteinData()
    {
        if(this.strArr==null)
        {
//            System.out.println(getSeqCount() + "--" +  getSpectrumCount() + "--" +  getSeqCoverage() + "--" +  getLength() + "--" +  getMolWt() + "--" +  getPI() + "--" +  getDescription() );
            String[] temp = {getLocus(), getSeqCount(), getSpectrumCount(), getSeqCoverage(), getLength(), getMolWt(), getPI(), getDescription(), };
            return temp;
        }

        return this.strArr;
    }

    public String getLocus()
    {
        return locus;
    }

    public String getSeqCount()
    {
        return seqCount;
    }

    public String getSpectrumCount()
    {
        return spectrumCount;
    }

    public int getSpectrumCountAsInt()
    {
        return Integer.parseInt(spectrumCount);
    }

    public String getSeqCoverage()
    {
        return seqCoverage;
    }

    public String getLength()
    {
        return length;
    }

    public String getMolWt()
    {
        return molWt;
    }

    public String getPI()
    {
        return pI;
    }

    public String getValidation()
    {
        return validation;
    }

    public String getDescription()
    {
        return description;
    }

    public void setLocus(String locus)
    {
        this.locus = locus;
    }

    public void setSeqCount(String seqCount)
    {
        this.seqCount = seqCount;
    }

    public void setSpectrumCount(String spectrumCount)
    {
        this.spectrumCount = spectrumCount;
    }

    public void setSeqCoverage(String seqCoverage)
    {
        this.seqCoverage = seqCoverage;
    }

    public void setLength(String length)
    {
        this.length = length;
    }

    public void setMolWt(String molWt)
    {
        this.molWt = molWt;
    }

    public void setPI(String pI)
    {
        this.pI = pI;
    }

    public void setValidation(String validation)
    {
        this.validation = validation;
    }

    public void setDescription(String description)
    {
        this.description = description;
    }

    public String getProteinLine()
    {
        return this.proteinLine;
    }

    public boolean isRedundant() {
        return redundant;
    }

    public void setRedundant(boolean redundant) {
        this.redundant = redundant;
    }

    public void setProteinLine(String proteinLine) {
        this.proteinLine = proteinLine;
    }

    public void setIntensity(Double intensity) {
	this.intensity = intensity;
    }

    public Double getIntensity() {
	return intensity;
    }

    public Double getAverageRatio() {
        return averageRatio;
    }

    public void setAverageRatio(Double averageRatio) {
        this.averageRatio = averageRatio;
    }

    public String getHspectrumCount() {
        return hspectrumCount;
    }

    public void setHspectrumCount(String hspectrumCount) {
        this.hspectrumCount = hspectrumCount;
    }

    public String getLspectrumCount() {
        return lspectrumCount;
    }

    public void setLspectrumCount(String lspectrumCount) {
        this.lspectrumCount = lspectrumCount;
    }

    public JSONObject getJSONObj()
    {
        JSONObject jsonObject = new JSONObject();
        jsonObject.put("accession", locus);
        jsonObject.put("seq_ct",seqCount );
        jsonObject.put("seq_cov", seqCoverage);
        jsonObject.put("length",length );
//        obj.put("molwt", molWt);
//        obj.put("pi", pI);
//        obj.put("val",validation );
        jsonObject.put("desc", description);
        JSONArray peptideArr = new JSONArray();
        for(ChroPeptide peptide : peptideList)
        {
            peptideArr.add(peptide.getJSONObj());
        }
        jsonObject.put("peptides",peptideArr );


        return jsonObject;
    }
    public static ChroProtein readJsonObject(JSONObject proteinObj)
    {
        ChroProtein protein  = new ChroProtein();
        if(proteinObj.containsKey("val"))
            protein.setValidation((String) proteinObj.get("val"));
        if(proteinObj.containsKey("molwt"))
            protein.setMolWt((String) proteinObj.get("molwt"));
        if(proteinObj.containsKey("desc"))
            protein.setDescription((String) proteinObj.get("desc"));
        if(proteinObj.containsKey("length"))
            protein.setLength((String) proteinObj.get("length"));
        if(proteinObj.containsKey("seq_ct"))
            protein.setSeqCount((String) proteinObj.get("seq_ct"));
        if(proteinObj.containsKey("accession"))
            protein.setLocus((String) proteinObj.get("accession"));

//        if(proteinObj.containsKey("peptides"))
//        {
//            JSONArray pepList =  (JSONArray) proteinObj.get("peptides");
//
//            for(Iterator pepItr = pepList.iterator();pepItr.hasNext();)
//            {
//                JSONObject peptideObj = (JSONObject) pepItr.next();
//                ChroPeptide peptide = new ChroPeptide();
//                peptide.readJsonObject(peptideObj);
//                protein.peptideList.add(peptide);
//            }
//        }
        return protein;
    }

    /**
     * @return the pepCount
     */
    public int getPepCount() {
        return pepCount;
    }

    /**
     * @param pepCount the pepCount to set
     */
    public void setPepCount(int pepCount) {
        this.pepCount = pepCount;
    }

    public List<Double> getPepRatioList() {
        return pepRatioList;
    }

    public void setPepRatioList(List<Double> pepRatioList) {
        this.pepRatioList = pepRatioList;
    }
    public void addPepRatioList(double pepRatio) {
        this.pepRatioList.add(pepRatio);
    }
    public List<Integer> getSpecCountList() {
        return specCountList;
    }

    public void setSpecCountList(List<Integer> specCountList) {
        this.specCountList = specCountList;
    }

    public void addSpecCountList(int specCount) {
        this.specCountList.add(specCount);
    }

    public List<Double> getAvgIntensityList() {
        return avgIntensityList;
    }

    public void setAvgIntensityList(List<Double> avgIntensityList) {
        this.avgIntensityList = avgIntensityList;
    }

    public void addAvgIntensityList(double value) {
        this.avgIntensityList.add(value);
    }




}
