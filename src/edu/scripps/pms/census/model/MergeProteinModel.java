/*
 * MergeProteinModel.java
 *
 * Created on February 23, 2006, 2:55 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package edu.scripps.pms.census.model;

import java.util.*;
import edu.scripps.pms.census.util.*;

/**
 *
 * @author rpark
 * @version $Id: MergeProteinModel.java,v 1.8 2012/04/17 23:58:33 rpark Exp $
 */
public class MergeProteinModel {
    
    private Set<String> proteinNames = new HashSet<String>();
    private Hashtable<String, String> proteinDescHt = new Hashtable<String, String>();
    private List<ChroProtein> proList = new Vector<ChroProtein>();
    private List<Peptide> pepList = new Vector<Peptide>();
    private Set<String> peptideSeqFileNames = new HashSet<String>(); //check if there are duplicate peptides
    private double totalRatio=0;
    private List descList = new Vector();
    
    /** Creates a new instance of MergeProteinModel */
    public MergeProteinModel() {
    }
   
    private int totalSpecCount=0;
 
    public void removeOutlier(double pValue)
    {
        if(pepList.size()>3)
            edu.scripps.pms.stats.GrubbsTest.filterMerge(pepList, pValue);                               
        
        int index=0;
        
        //peptideSeqFileNames.clear();

        this.totalRatio = 0;
        
        List<Peptide> tmpList = new Vector<Peptide>();
        
        for(int i=0;i<pepList.size();i++)
        //for(Iterator<Peptide> itr=pepList.iterator(); itr.hasNext(); )
        {            
            Peptide pep = pepList.get(i); //itr.next();
            if(!pep.isFilterOut())
            {
                String fName = pep.getFileName();
                //peptideSeqFileNames.add( pep.getSequence() + fName.substring(0, fName.indexOf(".")) );
                peptideSeqFileNames.add( pep.getSequence() + fName );
                totalRatio += pep.getRatio();
                tmpList.add(pep);
                
            }
            
            index++;
        }
        
        pepList.clear();
        pepList.addAll(tmpList);
            
    }    
           
    public List<String> getDescList()
    {
	return this.descList;
    }

    public void addProteinInfo(String proteinName, String specCount, String description, boolean addSpecCount, String pLine)
    {
	    if(addSpecCount && null != specCount && !"".equals(specCount))
	    {
		this.totalSpecCount += Integer.parseInt(specCount); 
//		System.out.println(totalSpecCount  + "===>>" + addSpecCount + " " + specCount + "<==");
	    }

	descList.add(description);

        if(!proteinNames.contains(proteinName))
        {
            proteinNames.add(proteinName);            
	    proteinDescHt.put(proteinName, description);

            proList.add(new ChroProtein(proteinName, specCount, description, pLine));
        }        
    }

    public Hashtable<String, String> getProteinDescHt()
    {
	return proteinDescHt;
    }
    
    public String getProteinDesc(String proteinName)
    {
	return proteinDescHt.get(proteinName);
    }
    
    public boolean isPeptideContain(String uniqueId)
    {
        for(Iterator<Peptide> itr1 = this.getPeptides().iterator(); itr1.hasNext(); )
        {
            Peptide pep = itr1.next();

            if(uniqueId.equals(pep.getUniqueIdentifier()))
                return true;
        }
        
        return false;
    }
    
    public void addPeptideList(List<Peptide> l)
    {
        for(Iterator<Peptide> itr=l.iterator(); itr.hasNext(); )
        {
            Peptide each = itr.next();
            
            
            //if( !peptideSeqFileNames.contains(each.getUniqueIdentifier()) )
            if( !this.isPeptideContain(each.getUniqueIdentifier()) )
            {
                pepList.add(each);                
                totalRatio += each.getRatio();
            }
        }        
    }
    
    public void addProteinList(List<ChroProtein> l)
    {
	int specCount = 0;
        for(Iterator<ChroProtein> itr=l.iterator(); itr.hasNext(); )
        {
            ChroProtein each = itr.next();
            
            if( !proteinNames.contains(each.getLocus()) )
            {
                proteinNames.add(each.getLocus());
		
                proList.add(each);                
            }
		String specC = each.getSpectrumCount();

		if(null != specC && !"".equals(specC))
		    specCount = Integer.parseInt(specC);

        }        

	this.totalSpecCount += specCount; 
    }
  
    public int getTotalSpecCount()
    {
	return this.totalSpecCount;

/*
	int sum = 0;
        for(Iterator<ChroProtein> itr=proList.iterator(); itr.hasNext(); )
        {
            ChroProtein each = itr.next();
	    sum += Integer.parseInt(each.getSpectrumCount());

	    System.out.println(sum);
	}

	    System.out.println("==>>" + sum);
	return sum;
            
*/
    }


    //weighted average
    public double getAverageRatio()
    {
	WeightedProtein.ProteinModel pModel = new WeightedProtein.ProteinModel();
	//for(Iterator<MergeProteinModel.Peptide> pepItr=protein.getPeptides().iterator(); pepItr.hasNext(); )

	if(pepList.size()<0)
	    return 0.0;
	else if(pepList.size()==1)
	{
	    Peptide pep = pepList.get(0);
	    return pep.getRatio();
	}
	    
        for(Iterator<Peptide> itr=pepList.iterator(); itr.hasNext(); )
	{
	    Peptide each = itr.next();
	    double invStdev = CalcUtil.getWeightedStdev(each.getRegFactor());

	    double devSum=0;
	    //grap only valid ratio
	    if(each.getRatio()>0)
		pModel.add(invStdev, each.getRatio());                                            
	}

	return pModel.getStandardWeightedAverage();

       // return totalRatio/pepList.size();       
    }
   
    //weighted average
    public double getStdev()
    {
        double devSum=0;
        double averageRatio = getAverageRatio();
        
        for(Iterator<Peptide> itr=pepList.iterator(); itr.hasNext(); )
        {
            Peptide pep = itr.next();
            double dev = pep.getRatio()-averageRatio;
            devSum += dev*dev;    
            
        }
        
        return Math.sqrt(devSum/(pepList.size()-1));        
    }
    
    /*
    public void addPeptide(boolean unique, String sequence, double ratio, double regFactor, double samIntensity, double refIntensity, String fileName, String peptideLine)
    {
        String tempName = fileName.substring(0, fileName.indexOf("."));
        //peptideSeqFileNames.add(sequence + tempName);        
        peptideSeqFileNames.add(sequence + fileName);
        pepList.add(new Peptide(unique, sequence, ratio, regFactor, samIntensity, refIntensity, fileName, peptideLine));
        totalRatio += ratio;
    }
    */
    public void addPeptide(boolean unique, String sequence, double ratio, double regFactor, double samIntensity, double refIntensity, double areaRatio, double profileScore, String fileName, String peptideLine)
    {
      //  String tempName = fileName.substring(0, fileName.indexOf("."));
        //peptideSeqFileNames.add(sequence + tempName);        
        peptideSeqFileNames.add(sequence + fileName);
        pepList.add(new Peptide(unique, sequence, ratio, regFactor, samIntensity, refIntensity, areaRatio, profileScore, fileName, peptideLine));
        totalRatio += ratio;
    }
    
    public void addPeptide(boolean unique, String sequence, double ratio, double regFactor, double samIntensity, double refIntensity, double areaRatio, double profileScore, String fileName)
    {
      //  String tempName = fileName.substring(0, fileName.indexOf("."));
        //peptideSeqFileNames.add(sequence + tempName);        
        peptideSeqFileNames.add(sequence + fileName);
        pepList.add(new Peptide(unique, sequence, ratio, regFactor, samIntensity, refIntensity, areaRatio, profileScore, fileName));
        totalRatio += ratio;
    }
    
    public boolean isContain(String proteinName)
    {
        return proteinNames.contains(proteinName);
    }
   
    private boolean mergeCheck(List<MergeProteinModel.Peptide> pepCompareList, Set<String> proteinNameSet, Hashtable<String, String> ht)
        throws Exception
    {
        boolean isMergeable=true;
        
        for(Iterator<Peptide> itr=pepCompareList.iterator();itr.hasNext(); )
        {
            Peptide peptide = itr.next();
            String pepSequence = peptide.getNoModSequence();
            
            pepSequence = pepSequence.substring(2, pepSequence.length()-2);

                       
            for(Iterator<String> proItr=proteinNameSet.iterator(); proItr.hasNext(); )
            {
                String locus = proItr.next();
                
                String seq = ht.get(locus);   
                
                if(seq == null)
                {
                    throw new Exception(locus + " was not found in the database file. It seems wrong database file was used");
                }
    
              //  System.out.println( locus + "\t" + pepSequence + "\t" + seq.indexOf(pepSequence) + "\t" + seq.indexOf(pepSequence));
                
                if( seq.indexOf(pepSequence)<0 )
                {
                    isMergeable=false;
                    break;
                }
                
                
//                System.out.println(pepSequence + "\t" + seq + "\t" + seq.indexOf(pepSequence));
            }
            
            if(!isMergeable)
                break;
        }
        
        return isMergeable;        
        
    }
    
    
    public boolean isMergeable(MergeProteinModel proteinToCompare, Hashtable<String, String> ht) throws Exception
    {
        //boolean isMergeable=true;

        /*
        mergeCheck(proteinToCompare.getPeptides(), this.getProteinNames(), ht);
        System.out.println("==================");
             mergeCheck(this.getPeptides(), proteinToCompare.getProteinNames(), ht);
             System.out.println("==================");
        
             return false;
        */
        
        
        return mergeCheck(proteinToCompare.getPeptides(), this.getProteinNames(), ht) 
            && mergeCheck(this.getPeptides(), proteinToCompare.getProteinNames(), ht);
    }
            
            
            /*
    public boolean isMergeable(List<Peptide> pepCompareList, Hashtable<String, String> ht) throws Exception
    {
        boolean isMergeable=true;
                
        for(Iterator<Peptide> itr=pepCompareList.iterator();itr.hasNext(); )
        {
            Peptide peptide = itr.next();
            //String pepSequence = peptide.getSequence();
            String pepSequence = peptide.getNoModSequence();
            
            pepSequence = pepSequence.substring(2, pepSequence.length()-2);
            
            //System.out.println("==>>pepSequence" + pepSequence);
            
            
            for(Iterator<String> proItr=this.proteinNames.iterator(); proItr.hasNext(); )
            {
                String locus = proItr.next();
                
                String seq = ht.get(locus);   
                
                if(seq == null)
                {
                    throw new Exception(locus + " was not found in the database file. It seems wrong database file was used");
                }
    
                
                if( seq.indexOf(pepSequence)<0 )
                {
                    isMergeable=false;
                    break;
                }
                
                
            }
            
            if(!isMergeable)
                break;
        }

  //      if(true)                
  //          System.exit(0);
        
        return isMergeable;
    }
    */
    
    public Set getProteinNames()
    {
	return proteinNames;
    }

    public List<ChroProtein> getProteins()
    {
        return proList;
    }
    
    public List<Peptide> getPeptides()
    {
        return pepList;
    }
    
    public class Peptide
    {
        private boolean unique;
        private String sequence;
        private double ratio;
        private double regFactor;
        private String fileName;
        private String uniqueIdentifier; // sequence + fileName
        private double probability;
        private boolean filterOut=false;
        private double samIntensity=-1;
        private double refIntensity=-1;
	private String peptideLine;

	private double areaRatio;
	private double profileScore;


/*
        public Peptide(boolean unique, String sequence, double ratio, double regFactor, double samIntensity, double refIntensity, String fileName, String peptideLine)
	{
	    this(unique, sequence, ratio, regFactor, samIntensity, refIntensity, fileName);
	    this.peptideLine = peptideLine;
	}
*/	
        public Peptide(boolean unique, String sequence, double ratio, double regFactor, double samIntensity, double refIntensity, double areaRatio, double profileScore, String fileName)
	{
            this.unique = unique;
            this.sequence = sequence;
            this.ratio = ratio;
            this.regFactor = regFactor;
            this.fileName = fileName;
            this.samIntensity = samIntensity;
            this.refIntensity = refIntensity;
            
            this.areaRatio = areaRatio;
            this.profileScore = profileScore;
            
            this.setUniqueIdentifier(sequence + fileName);


	}


        public Peptide(boolean unique, String sequence, double ratio, double regFactor, double samIntensity, double refIntensity, double areaRatio, double profileScore, String fileName, String peptideLine)
        {
	/*
            this.unique = unique;
            this.sequence = sequence;
            this.ratio = ratio;
            this.regFactor = regFactor;
            this.fileName = fileName;
            this.samIntensity = samIntensity;
            this.refIntensity = refIntensity;
            
            this.areaRatio = areaRatio;
            this.profileScore = profileScore;
            this.setUniqueIdentifier(sequence + fileName);
*/
	    this(unique, sequence, ratio, regFactor, samIntensity, refIntensity, areaRatio, profileScore, fileName);
	    this.peptideLine = peptideLine;
            
        }
       
	public String getPeptideLine() {
	    return this.peptideLine;
	}

	public void setPeptideLine(String peptideLine) {
	    this.peptideLine = peptideLine;
	}
	
        public boolean isUnique() {
            return unique;
        }

        public void setUnique(boolean unique) {
            this.unique = unique;
        }

        public String getSequence() {
            return sequence;
        }

        public String getNoModSequence() {
            return sequence.replaceAll("[*#@]", "");
        }
        
        public void setSequence(String sequence) {
            this.sequence = sequence;
        }

        public double getRatio() {
            return ratio;
        }

        public void setRatio(double ratio) {
            this.ratio = ratio;
        }
        
        public double getRegFactor() {
            return regFactor;
        }

        public void setRegFactor(double regFactor) {
            this.regFactor = regFactor;
        }

        public String getFileName() {
            return fileName;
        }

        public void setFileName(String fileName) {
            this.fileName = fileName;
        }

        public String getUniqueIdentifier() {
            return uniqueIdentifier;
        }

        public void setUniqueIdentifier(String uniqueIdentifier) {
            this.uniqueIdentifier = uniqueIdentifier;
        }

        public double getProbability() {
            return probability;
        }

        public void setProbability(double probability) {
            this.probability = probability;
        }

        public boolean isFilterOut() {
            return filterOut;
        }

        public void setFilterOut(boolean filterOut) {
            this.filterOut = filterOut;
        }

        public double getSamIntensity() {
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

        public double getAreaRatio() {
            return areaRatio;
        }

        public void setAreaRatio(double areaRatio) {
            this.areaRatio = areaRatio;
        }

        public double getProfileScore() {
            return profileScore;
        }

        public void setProfileScore(double profileScore) {
            this.profileScore = profileScore;
        }
        
        
    }

}

