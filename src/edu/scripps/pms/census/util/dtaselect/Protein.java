
package edu.scripps.pms.census.util.dtaselect;

/**
 * @author  Robin Park
 * @version $Id: Protein.java,v 1.1 2014/09/09 19:28:27 rpark Exp $
 */

import edu.scripps.pms.census.model.*;
import edu.scripps.pms.util.seq.Fasta;
import java.util.*;

public class Protein
{
    private String locus;
    private String seqCount = "N/A";
    private String spectrumCount = "N/A";
    private String seqCoverage = "N/A";
    private String length = "N/A";
    private String molWt = "N/A";
    private String pI = "N/A";
    private String validation = "N/A";
    private String description = "N/A";
    private List<Peptide> peptideList;
    private String proteinLine="N/A";
    private Hashtable<String, Integer> proteinIndexHt;
    private String lSpectrumCount = "N/A";
    private String hSpectrumCount = "N/A";
    private String nsaf=null;
    private String empai=null;
    
    //for non-labeling union
    private Hashtable<String, Peptide> peptideHt;
	private String heavySpec;
	private String lightSpec;
    private boolean isProblematic = false;
    private ProteinGroup proteinGroup;

    private static int locusIndex = -1;
    private static int seqCountIndex = -1;
    private static int spectrumCountIndex = -1;
    private static int seqCoverageIndex = -1;
    private static int lengthIndex = -1;
    private static int molWtIndex = -1;
    private static int pIIndex = -1;
    private static int validationIndex = -1;
    private static int descriptionIndex = -1;
    private static int heavySpecIndex = -1;
    private static int lightSpecIndex = -1;
    private static int listOfExpNamesIndex = -1;
    private static int listOfSearchNamesIndex = -1;
    private static int listOfSeqCountIndex = -1;
    private static int listOfSpecCountIndex = -1;
    private static int listOfSeqCoverageIndex = -1;
    private static int nsafIndex = -1;
    private static int empaiIndex = -1;

    private boolean valid=true;
        
    public Protein(String proteinLine) throws ArrayIndexOutOfBoundsException
    {
        this( proteinLine.split("\t") );
        this.proteinLine = proteinLine;
    }

    public int addPeptideHt(Iterator<Peptide> pepList, String filePath)
    {
        int count=0;
        
        for(Iterator<Peptide> itr=pepList; itr.hasNext(); )
        {
            Peptide pep = itr.next();
            pep.setFilePath(filePath);
            
            Peptide tempPep = peptideHt.get(pep.getSequence() + pep.getChargeState());
         
            if(null == tempPep)
            {
                peptideHt.put(pep.getSequence() + pep.getChargeState(), pep);
                count++;
                
                continue;
            }
            
            if( Double.parseDouble(tempPep.getXCorr()) < Double.parseDouble(pep.getXCorr()) )
            {
                peptideHt.put(pep.getSequence() + pep.getChargeState(), pep);                                
            }
        }
        
        return count;
        
    }

    public void populatePeptideHt(String filePath)
    {
        for(Iterator<Peptide> itr=this.peptideList.iterator(); itr.hasNext(); )
        {
            Peptide pep = itr.next();
            pep.setFilePath(filePath);
            peptideHt.put(pep.getSequence() + pep.getChargeState(), pep);
        }
    }
   
    public Protein(String[] strArr) throws ArrayIndexOutOfBoundsException
    {
        peptideList = new ArrayList<Peptide>();

        init(strArr);
    }

   
    public void setPeptideList(List<Peptide> list) {
        peptideList = list;
    }
    public List<Peptide> getPeptideList() {
        return peptideList;
    }
    public void setProteinGroup(ProteinGroup pg) {
        proteinGroup = pg;
    }

    public ProteinGroup getProteinGroup() {
        return proteinGroup;
    }

    public boolean equals(Object o) {

        Protein p = (Protein) o;
        return p == null? false : getAccession().equals(p.getAccession());
    }
    public int hashCode() {
        return getAccession().hashCode();
    }

    public String getAccession() {
        return getLocus();
    }
    public String getAccessionWithoutVersion() {
        return Fasta.getAccessionWithNoVersion(getLocus());
    }
    public Protein()
    {
        peptideList = new ArrayList<Peptide>();
    }

    public void setElement(String[] strArr)
    {
	init(strArr);
        //System.out.println("in Protein, seqCoverage: " + seqCoverage);
    }

    public boolean addPeptide(Peptide peptide)
    {
        return peptideList.add(peptide);
    }

    public void removePeptide(Peptide peptide)
    {
        peptideList.remove(peptide);
    }

    public Iterator<Peptide> getPeptides()
    {
        return peptideList.iterator();
    }

    
    private void init(String[] strArr) throws ArrayIndexOutOfBoundsException
    {
        try {
            this.setLocus( Fasta.getAccession(strArr[this.getLocusIndex()]) );
            this.setSeqCount(strArr[this.getSeqCountIndex()]);
            this.setSpectrumCount(strArr[this.getSpectrumCountIndex()]);
            this.setSeqCoverage(strArr[this.getSeqCoverageIndex()]);
            this.setLength(strArr[this.getLengthIndex()]);
            this.setMolWt(strArr[this.getMolWtIndex()]);
            this.setPI(strArr[this.getpIIndex()]);
            this.setValidation(strArr[this.getValidationIndex()]);
            this.setDescription(strArr[this.getDescriptionIndex()]);
            if(this.getNsafIndex()>=0) {
                this.setNsaf(strArr[this.getNsafIndex()]);
                this.setEmpai(strArr[this.getEmpaiIndex()]);
            }
            if(strArr[this.getHeavySpecIndex()]!=null && strArr[this.getLightSpecIndex()]!=null){
                    this.setHeavySpec(strArr[this.getHeavySpecIndex()]);
                    this.setLightSpec(strArr[this.getLightSpecIndex()]);
            }


            /*
            this.setLocus( Fasta.getAccession(strArr[0]) );
            this.setSeqCount(strArr[1]);
            this.setSpectrumCount(strArr[2]);
            this.setSeqCoverage(strArr[3]);
            this.setLength(strArr[4]);
            this.setMolWt(strArr[5]);
            this.setPI(strArr[6]);
            this.setValidation(strArr[7]);
            this.setDescription(strArr[8]);
		if(strArr[9]!=null && strArr[10]!=null){
			this.setHeavySpec(strArr[9]);
			this.setLightSpec(strArr[10]);
		}
                */
        }
        catch(ArrayIndexOutOfBoundsException ex) {
            isProblematic = true;
            //throw new ArrayIndexOutOfBoundsException("Error : Mal formed DTASelect-filter.txt file in protein line: " + strArr[0]);
        }
    }
    public boolean isProblematic() {
        return isProblematic;
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

    public String getSeqCoverage()
    {
        return seqCoverage;
    }

    public String getLength()
    {
        return length;
    }

    public int getLengthValue()
    {
        return Integer.parseInt(length);
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
    	public void setHeavySpec(String hs){
		this.heavySpec = hs;
    	}
	public void setLightSpec(String ls){
		this.lightSpec = ls;
	}
	public String getHeavySpec(){
		return heavySpec;
	}
	public String getLightSpec(){
		return lightSpec;
	}
    public int getPeptideSize()
    {
        return peptideList.size();
    }
    public int getNumPeptides()
    {
        return peptideList.size();
    }

    public String getProteinLine()
    {
        return proteinLine;
    }
    public Hashtable<String, Peptide> getPeptideHt() {
        return peptideHt;
    }

    public void setPeptideHt(Hashtable<String, Peptide> peptideHt) {
        this.peptideHt = peptideHt;
    }

    public String gethSpectrumCount() {
        return hSpectrumCount;
    }

    public void sethSpectrumCount(String hSpectrumCount) {
        this.hSpectrumCount = hSpectrumCount;
    }

    public String getlSpectrumCount() {
        return lSpectrumCount;
    }

    public void setlSpectrumCount(String lSpectrumCount) {
        this.lSpectrumCount = lSpectrumCount;
    }

    public String getpI() {
        return pI;
    }

    public void setpI(String pI) {
        this.pI = pI;
    }

    public Hashtable<String, Integer> getProteinIndexHt() {
        return proteinIndexHt;
    }

    public void setProteinIndexHt(Hashtable<String, Integer> proteinIndexHt) {
        this.proteinIndexHt = proteinIndexHt;
    }

    public boolean isIsProblematic() {
        return isProblematic;
    }

    public void setIsProblematic(boolean isProblematic) {
        this.isProblematic = isProblematic;
    }

    public static int getLocusIndex() {
        return locusIndex;
    }

    public static void setLocusIndex(int locusIndex) {
        Protein.locusIndex = locusIndex;
    }

    public static int getSeqCountIndex() {
        return seqCountIndex;
    }

    public static void setSeqCountIndex(int seqCountIndex) {
        Protein.seqCountIndex = seqCountIndex;
    }

    public static int getSpectrumCountIndex() {
        return spectrumCountIndex;
    }

    public static void setSpectrumCountIndex(int spectrumCountIndex) {
        Protein.spectrumCountIndex = spectrumCountIndex;
    }

    public static int getSeqCoverageIndex() {
        return seqCoverageIndex;
    }

    public static void setSeqCoverageIndex(int seqCoverageIndex) {
        Protein.seqCoverageIndex = seqCoverageIndex;
    }

    public static int getLengthIndex() {
        return lengthIndex;
    }

    public static void setLengthIndex(int lengthIndex) {
        Protein.lengthIndex = lengthIndex;
    }

    public static int getMolWtIndex() {
        return molWtIndex;
    }

    public static void setMolWtIndex(int molWtIndex) {
        Protein.molWtIndex = molWtIndex;
    }

    public static int getpIIndex() {
        return pIIndex;
    }

    public static void setpIIndex(int pIIndex) {
        Protein.pIIndex = pIIndex;
    }

    public static int getValidationIndex() {
        return validationIndex;
    }

    public static void setValidationIndex(int validationIndex) {
        Protein.validationIndex = validationIndex;
    }

    public static int getDescriptionIndex() {
        return descriptionIndex;
    }

    public static void setDescriptionIndex(int descriptionIndex) {
        Protein.descriptionIndex = descriptionIndex;
    }

    public static int getHeavySpecIndex() {
        return heavySpecIndex;
    }

    public static void setHeavySpecIndex(int heavySpecIndex) {
        Protein.heavySpecIndex = heavySpecIndex;
    }

    public static int getLightSpecIndex() {
        return lightSpecIndex;
    }

    public static void setLightSpecIndex(int lightSpecIndex) {
        Protein.lightSpecIndex = lightSpecIndex;
    }

    public static int getListOfExpNamesIndex() {
        return listOfExpNamesIndex;
    }

    public static void setListOfExpNamesIndex(int listOfExpNamesIndex) {
        Protein.listOfExpNamesIndex = listOfExpNamesIndex;
    }

    public static int getListOfSearchNamesIndex() {
        return listOfSearchNamesIndex;
    }

    public static void setListOfSearchNamesIndex(int listOfSearchNamesIndex) {
        Protein.listOfSearchNamesIndex = listOfSearchNamesIndex;
    }

    public static int getListOfSeqCountIndex() {
        return listOfSeqCountIndex;
    }

    public static void setListOfSeqCountIndex(int listOfSeqCountIndex) {
        Protein.listOfSeqCountIndex = listOfSeqCountIndex;
    }

    public static int getListOfSpecCountIndex() {
        return listOfSpecCountIndex;
    }

    public static void setListOfSpecCountIndex(int listOfSpecCountIndex) {
        Protein.listOfSpecCountIndex = listOfSpecCountIndex;
    }

    public static int getListOfSeqCoverageIndex() {
        return listOfSeqCoverageIndex;
    }

    public static void setListOfSeqCoverageIndex(int listOfSeqCoverageIndex) {
        Protein.listOfSeqCoverageIndex = listOfSeqCoverageIndex;
    }

    public static int getNsafIndex() {
        return nsafIndex;
    }

    public static void setNsafIndex(int nsafIndex) {
        Protein.nsafIndex = nsafIndex;
    }

    public static int getEmpaiIndex() {
        return empaiIndex;
    }

    public static void setEmpaiIndex(int empaiIndex) {
        Protein.empaiIndex = empaiIndex;
    }

    public static void setFeatureIndices(String features) {
        clearIndex();
        String [] contents = features.split("\t");
//Locus                                                    
        
        for(int i = 0; i < contents.length; i++) {
            String s = contents[i].trim();
            locusIndex = s.startsWith("Locus")? i :locusIndex;
            seqCountIndex = s.startsWith("Sequence Count")? i :seqCountIndex;
            spectrumCountIndex = s.startsWith("Spectrum Count")? i : spectrumCountIndex;
            seqCoverageIndex = s.startsWith("Sequence Coverage")? i :seqCoverageIndex;
            lengthIndex = s.startsWith("Length")? i :lengthIndex;
            molWtIndex = s.startsWith("MolWt")? i :molWtIndex;
            pIIndex = s.startsWith("pI")? i :pIIndex;
            validationIndex = s.startsWith("Validation Status")? i :validationIndex;
            descriptionIndex = s.startsWith("Descriptive Name")? i :descriptionIndex;
            heavySpecIndex = s.startsWith("HRedundancy")? i :heavySpecIndex;
            lightSpecIndex = s.startsWith("LRedundancy")? i :lightSpecIndex;
            nsafIndex = s.startsWith("NSAF")? i :nsafIndex;
            empaiIndex = s.startsWith("EMPAI")? i :empaiIndex;

            //listOfSeqCountIndex = s.startsWith("")? i :listOfSeqCountIndex;
            //listOfSpecCountIndex = s.startsWith("")? i :listOfSpecCountIndex;
            //listOfSeqCoverageIndex = s.startsWith("")? i :listOfSeqCoverageIndex;

        }
	    //System.out.println("conf: " + confIndex);
    }
    
    public static void clearIndex(){
        seqCountIndex = -1;
         spectrumCountIndex = -1;
         seqCoverageIndex = -1;
         lengthIndex = -1; 
         molWtIndex = -1;
         pIIndex = -1;
         validationIndex = -1;
         descriptionIndex = -1;
         heavySpecIndex = -1;
         lightSpecIndex = -1;
         listOfExpNamesIndex = -1;
         listOfSearchNamesIndex = -1;
         listOfSeqCountIndex = -1;
         listOfSpecCountIndex = -1;
         listOfSeqCoverageIndex = -1;
         nsafIndex = -1;
         empaiIndex = -1;
         //ptmIndexIndex = -1;
         //ptmIndexProteinIndex = -1;

    }        

    public String getNsaf() {
        return nsaf;
    }

    public void setNsaf(String nsaf) {
        this.nsaf = nsaf;
    }

    public String getEmpai() {
        return empai;
    }

    public void setEmpai(String empai) {
        this.empai = empai;
    }

    public boolean isValid() {
        return valid;
    }

    public void setValid(boolean valid) {
        this.valid = valid;
    }


    public void setProteinLine(String proteinLine) {
        this.proteinLine = proteinLine;
    }

    public void setProblematic(boolean problematic) {
        isProblematic = problematic;
    }
}
