/*
 * PepXmlReader.java
 *
 * Created on July 30, 2007, 3:07 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package edu.scripps.pms.census.util.io;

import java.io.IOException;
import java.util.*;

import edu.scripps.pms.census.util.dtaselect.Protein;
import edu.scripps.pms.census.util.dtaselect.Peptide;

import org.jdom.Document;
import org.jdom.Element;
import org.jdom.JDOMException;
import org.jdom.input.SAXBuilder;
import org.jdom.Namespace;

/**
 *
 * @author rpark
 */
public class PepXmlReader extends BaseIdentificationReader {
    private String fileName;
    private SAXBuilder builder = new SAXBuilder(); // create XML Builder/Parser
    private Document doc = null;
    private Element root = null;
    private Namespace ns = null;
    private int totalPeptideNumber=0;
    
    public static void main(String[] args) throws Exception
    {
	BaseIdentificationReader idReader = new PepXmlReader(args[0]);
            Iterator<Protein> pitr = idReader.getProteins(); //need to run to calculate redundnat peptides
	//new PepXmlReader("/home/rpark/rpark_on_data/project/census/mzxml/Thibault_Data/test/small_pep.xml");
    }

    /** Creates a new instance of PepXmlReader */
    public PepXmlReader(String fileName) throws JDOMException, IOException {
        this.setFileName(fileName);
        
        doc = builder.build(fileName);
        root = doc.getRootElement(); // get XML's root == <ResultFile> tag
	ns = root.getNamespace();

	reOrganizeProteins();

	calcTotalPeptideNum();
    }

    private void calcTotalPeptideNum()
    {
	List msmsRunList = root.getChildren("msms_run_summary", ns);
	for(Iterator<Element> mitr=msmsRunList.iterator(); mitr.hasNext(); )
	{
	    Element msmsEle = mitr.next();
	    List spectrumList = msmsEle.getChildren("spectrum_query", ns);

	    for(Iterator<Element> itr=spectrumList.iterator(); itr.hasNext(); )
	    {
		Element sQuery = itr.next();
		Element searchResultEle = sQuery.getChild("search_result", ns);

		if(null == searchResultEle)
		    continue;

		totalPeptideNumber++;
	    }
	}
    }


    public void reOrganizeProteins()
    {
	List<Protein> proteinList = new ArrayList<Protein>();

	Hashtable<String, Protein> proteinHt = new Hashtable<String, Protein>();
	List msmsRunList = root.getChild("msms_run_summary", ns).getChildren("spectrum_query", ns);
    }
   
    private class ModResidue implements Comparable {
	private int pos;
	private double mass;

	public ModResidue(int pos, double mass)
	{
	    this.pos = pos;
	    this.mass = mass;
	}

	public int compareTo(Object o) throws ClassCastException {

	    ModResidue mr = (ModResidue)o;

	    //return this.getPos() - mr.getPos();
	    return mr.getPos() - this.getPos();
	}

        public int getPos() {
            return pos;
        }

        public void setPos(int pos) {
            this.pos = pos;
        }

        public double getMass() {
            return mass;
        }

        public void setMass(double mass) {
            this.mass = mass;
        }
    }
    
    public Iterator <Protein> getProteins() throws IOException
    {
	//List<Protein> proteinList = new ArrayList<Protein>();
	Hashtable<String, Protein> proteinHt = new Hashtable<String, Protein>();

	List msmsRunList = root.getChildren("msms_run_summary", ns);

//	List msmsRunList = root.getChild("msms_run_summary", ns).getChildren("spectrum_query", ns);

	for(Iterator<Element> mitr=msmsRunList.iterator(); mitr.hasNext(); )
	{
	    Element msmsEle = mitr.next();

	    List spectrumList = msmsEle.getChildren("spectrum_query", ns);
	    String baseName = msmsEle.getAttributeValue("base_name");

	    if(baseName.startsWith("/"))
		baseName = baseName.substring(baseName.lastIndexOf("/") + 1);
	    else
		baseName = baseName.substring(baseName.lastIndexOf("\\") + 1); //windows

	    baseName += ".mzXML";

	    for(Iterator<Element> itr=spectrumList.iterator(); itr.hasNext(); )
	    {
		Element sQuery = itr.next();
		Element searchResultEle = sQuery.getChild("search_result", ns);


		if(null == searchResultEle)
		    continue;

		Element searchHit = searchResultEle.getChild("search_hit", ns);
		if(null == searchHit)
		    continue;

		String locus = searchHit.getAttributeValue("protein");

		Protein protein = proteinHt.get(locus);
		String seq = searchHit.getAttributeValue("peptide");

		//<modification_info modified_peptide="LAEVEAALEK[239]QR">
		//<mod_aminoacid_mass position="10" mass="239.174100"/>
		//</modification_info>
		Element modInfo = searchHit.getChild("modification_info", ns);

		if(null != modInfo)
		{

		    Element modAminoMassEle = modInfo.getChild("mod_aminoacid_mass", ns);

		    Hashtable<String, String> tmpModHt = new Hashtable<String, String>();
		    

		    List modList =  modInfo.getChildren("mod_aminoacid_mass", ns);
		    List<ModResidue> modPosList = new ArrayList<ModResidue>();

		    for(Iterator<Element> modItr=modList.iterator(); modItr.hasNext(); )
		    {
			Element tmpModEle = modItr.next();
			int pos = Integer.parseInt(tmpModEle.getAttributeValue("position"));
			double mass = Double.parseDouble(tmpModEle.getAttributeValue("mass"));

			modPosList.add( new ModResidue( pos, mass) );
		    }

		    Collections.sort(modPosList);

		    String modNtermMass = modInfo.getAttributeValue("mod_nterm_mass");
		    String modCtermMass = modInfo.getAttributeValue("mod_cterm_mass");


		    for(Iterator<ModResidue> modItr=modPosList.iterator(); modItr.hasNext(); )
		    {
			ModResidue modRes = modItr.next();

		//	System.out.println(modRes.getMass());
		//	System.out.println(modRes.getPos());
		//	System.out.println(seq);
			
			seq = seq.substring(0, modRes.getPos()) + "[" + modRes.getMass() + "]" + seq.substring(modRes.getPos());
		    }
				/*
		    String tmpseq = modInfo.getAttributeValue("modified_peptide");

		    if(null != modAminoMassEle)
		    {
			if(null == tmpseq)
			{
			    int aaPos = Integer.parseInt( modAminoMassEle.getAttributeValue("position") );
			    String modMass = modAminoMassEle.getAttributeValue("mass");

			    seq = seq.substring(0, aaPos) + "[" + modMass + "]" + seq.substring(aaPos);
			}
			else 
			{
			    String mstr = modAminoMassEle.getAttributeValue("mass");
			    if(null != mstr)
			    {
				double d = Double.parseDouble(mstr);
				if(tmpseq.startsWith("n"))
				    tmpseq = tmpseq.substring(1);

				char[] chArr = tmpseq.toCharArray();
				StringBuffer newSeq = new StringBuffer();

				boolean isModValue = false;

				int aaIndex=0;
				char prevRes = '0'; //just add zero

				for( char c : chArr)
				{
				    if(Character.isDigit(c))
					continue;

				    if(c == '[')
				    {
					//System.out.println(d + " " + isoReader.getAAMonoMass(prevRes) + " " + prevRes);
if(seq.equals("KITDRMIPK"))
System.out.println("====>>" + d + " "  + tmpseq);
					if(prevRes == '0')
					    newSeq.append(c).append(d);
					else
					    newSeq.append(c).append(d-isoReader.getAAMonoMass(prevRes));

					continue;
				    }

				    if( c == ']')
				    {
					newSeq.append(c);
					continue;
				    }

				    newSeq.append(c);

				    prevRes = c;

				    aaIndex++;

if(seq.equals("KITDRMIPK"))
    System.out.println("------------->>" + aaIndex + "\t" + c);

				}

if(seq.equals("KITDRMIPK"))
				seq = newSeq.toString();
			    }
			}
		    }
				*/

		    if(null != modNtermMass)
			seq = "[" + modNtermMass + "]" + seq;

		    if(null != modCtermMass)
			seq = seq + "[" + modCtermMass + "]";

		}

		seq = searchHit.getAttributeValue("peptide_prev_aa") + 
		    "." + seq + "." + 
		    searchHit.getAttributeValue("peptide_next_aa");


		if(null==protein)
		{
		    Protein p = new Protein();
		    p.setLocus( edu.scripps.pms.util.seq.Fasta.getAccession(locus) );
		    p.setDescription(locus);

		    p.addPeptide( this.buildPeptide(seq, searchHit, sQuery, baseName) );

		    proteinHt.put(locus, p);
		}
		else
		{
		    protein.addPeptide( this.buildPeptide(seq, searchHit, sQuery, baseName) );                    
		}

	    }
	}

	List<Protein> pList = new ArrayList<Protein>();

	for(Enumeration<Protein> e=proteinHt.elements(); e.hasMoreElements(); )
	{
	    pList.add(e.nextElement());
	}

	return pList.iterator();


	//return Collections.EMPTY_LIST.iterator();
    }
    
    public Peptide buildPeptide(String seq, Element searchHit, Element sQuery, String fileName)
    {
        Peptide pep = new Peptide();
        pep.setSequence(seq);
        //pep.setFileName(sQuery.getAttributeValue("spectrum"));
        pep.setFileName(fileName);
        pep.setRedundancy("0");

	pep.setChargeState(sQuery.getAttributeValue("assumed_charge"));
        
        String scanNum = null;
        
        try {
            scanNum = sQuery.getAttributeValue("start_scan");
            
            if(null == scanNum || "0".equals(scanNum))
            {
                String[] sarr = sQuery.getAttributeValue("spectrum").split("\\.");
                scanNum = sarr[sarr.length-2];
            }
            
        } catch (Exception e) {
            System.out.println("Scan Number in pepXML is invalid");
            System.exit(0);
        }
                
        pep.setScanNum(scanNum);

        for(Iterator<Element> scoreItr=searchHit.getChildren("search_score", ns).iterator(); scoreItr.hasNext(); )
        {
            Element scoreEle = scoreItr.next();
	    
	    String name = scoreEle.getAttributeValue("name");
	    String value = scoreEle.getAttributeValue("value");

	    pep.addScore(name, value);

	    if("xcorr".equals(name))
		pep.setXCorr(value);
	    else if("deltacn".equals(name))
		pep.setDeltCN(value);

/*

	    Hashtable<String, Double> scoreHt = new Hashtable<String, Double>();
	    for(Iterator<String> itrScr=scoreHt.keySet().iterator(); itrScr.hasNext(); )
	    {
		String score = itrScr.next();
		double value = scoreHt.get(scoreHt);

		Element scoreEle = new Score
	    }


	    <search_score name="ionscore" value="20.91"/>
		<search_score name="identityscore" value="48.64"/>
		<search_score name="star" value="0"/>
		<search_score name="homologyscore" value="31.40"/>
		<search_score name="expect" value="29.64"/>
*/

/*
	    
	    
	    
            String xCorr = scoreEle.getAttributeValue("xcorr");
            if(null != xCorr)
                pep.setXCorr(xCorr);

            String deltacn = scoreEle.getAttributeValue("deltacn");                        
            if(null != deltacn)
                pep.setDeltCN(deltacn);                        
*/		
        }

        return pep;
    }
    public int getTotalPeptideNumber() throws IOException
    {
        return totalPeptideNumber;
        
    }

    public String getFileName() {
        return fileName;
    }

    public void setFileName(String fileName) {
        this.fileName = fileName;
    }
    
}
