package edu.scripps.pms.census.util.dtaselect;

/**
 * @author Tao Xu 
 * @version $Id: ProteinGroup.java
 */
import java.util.*;
import edu.scripps.pms.util.seq.Fasta;

public class ProteinGroup {
    private List<Peptide> peptideList = new ArrayList();
    private List<Protein> proteinList = new ArrayList<Protein>(5);
    private Protein representative;
    public static final String DECOYLABEL = "Reverse";

    public ProteinGroup(Protein p) throws ArrayIndexOutOfBoundsException {
        representative = p;
        addProtein(p); 
    }

    public void addProtein(Protein p) {
        if(p== null) return;

        if(p.getNumPeptides() > getNumPeptides()) {
            peptideList = p.getPeptideList();
            //peptideList = new ArrayList<Peptide>();
            //for(Iterator<Peptide> it = p.getPeptides(); it.hasNext();) {
            //    peptideList.add(it.next());
            //}
        }
        if(p.getLengthValue() > representative.getLengthValue()) {
            representative = p;
        } 
        p.setProteinGroup(this);
        proteinList.add(p);
    }

   public String getDescription() {

        return representative.getDescription();
   }

    public ArrayList<Peptide> getHighestSpectralCountPeptide(int num) {
        ArrayList<Peptide> maxpeptides = new ArrayList(num);
        PeptideComparator pc = new PeptideComparator(1); // 1 means sort by spectral count
        Collections.sort(peptideList, pc);
        for(int i = 0; i < num && i < peptideList.size(); i++) {
            maxpeptides.add(peptideList.get(i));
        }
//System.out.println("num pepitdes: " + peptideList.size() + "\tnumber of peptide added: " + maxpeptides.size());
        return maxpeptides;
    }
    public Peptide getHighestSpectralCountPeptide() {
        int maxspectralcount = -1;
        double maxZScore = -1;
        Peptide maxpeptide = null;
        for(Iterator<Peptide> it = peptideList.iterator(); it.hasNext();) {
            Peptide p = it.next(); 
            if(it != null) {
                int spectralcount = p.getSpectralCount(); 
                if(spectralcount > maxspectralcount) {
                    maxpeptide = p;
                } else if(spectralcount == maxspectralcount) {

                    double zscore = p.getSpScoreValue();
                    if(zscore > maxZScore) {
                        maxZScore = zscore;
                        maxpeptide = p;
                    }
                }
           
            }
        } 
        return maxpeptide;
    }

    public Protein getRepresentativeProtein() {
        return representative;
    }

  
    public Protein getShortestProtein() {
        int minlength = 10000000;
        Protein shortest = null; 
        for(Iterator<Protein> it = proteinList.iterator(); it.hasNext();) {
            Protein p = it.next();
            if(p != null && p.getLengthValue() < minlength) {
                shortest = p;
                minlength = p.getLengthValue();
            }
        }
        return shortest;
    }

    public Protein getLongestProtein() {
        int maxlength = -1;
        Protein longest = null; 
        for(Iterator<Protein> it = proteinList.iterator(); it.hasNext();) {
            Protein p = it.next();
            if(p != null && p.getLengthValue() > maxlength) {
                longest = p;
                maxlength = p.getLengthValue();
            }
        }
        return longest;
    }
    public boolean isDecoyHit() {
   
        for(Iterator<Protein> it = proteinList.iterator(); it.hasNext();) {
            Protein p = it.next();
            if(p.getAccession().startsWith(DECOYLABEL)) return true;
        }        

        return false;         

    }
    public List<Peptide> getPeptideList() {
        return peptideList;
    }
    public void setPeptideList(List<Peptide> list) {
        peptideList = list;
    }
    public int getNumPeptides() {
        return peptideList.size();
    }
    public String getPeptideCount() {
        return representative.getSeqCount();
    }
    public String getSpectrumCount() {
        return representative.getSpectrumCount();
    }
    public String getSeqCoverage(){
        return representative.getSeqCoverage();
    }
    public Iterator<Peptide> getPeptides() {
        return peptideList.iterator();
    }
    public Iterator<Protein> getProteins() {
        return proteinList.iterator();
    }
    public boolean equals(Object o) {
        ProteinGroup p = (ProteinGroup)o;
        return o == null? false : representative.equals(p.getRepresentativeProtein());
    }
    public int hashCode() {
        return representative.getAccession().hashCode(); 
    }
}
