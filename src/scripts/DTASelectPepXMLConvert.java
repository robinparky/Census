package scripts;

import edu.scripps.pms.census.util.io.DTASelectFilterReader;
import java.util.*;
import java.io.*;


import java.io.BufferedReader;
import java.util.*;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.FileReader;
import java.io.RandomAccessFile;

import edu.scripps.pms.census.util.dtaselect.Protein;
import edu.scripps.pms.census.util.dtaselect.Peptide;

import edu.scripps.pms.census.model.*;

import org.jdom.*;
import org.jdom.output.*;
import org.jdom.input.*;

public class DTASelectPepXMLConvert
{
    public static void main(String args[]) throws Exception
    {
	//Convert DTASelect to pepXML
	if(args.length<1)
	{
	    System.out.println("Usage: java DTASelectPepXMLConvert DTASelect-filter.txt");
	    System.exit(0);
	}

	//build jdom
	Element rootEle = new Element("msms_pipeline_analysis");

        DTASelectFilterReader reader = new DTASelectFilterReader( args[0] );

        Iterator<Protein> pitr = reader.getProteins();

        Protein protein;
	Peptide peptide;
        ArrayList<Protein> aList = new ArrayList<Protein>();

	HashSet set = new HashSet();

	int i=0;
        for (Iterator<Protein> itr = pitr; itr.hasNext(); )
        {
            protein = itr.next();

//            if( protein.getLocus().startsWith("Rever") || protein.getLocus().startsWith("contam") || protein.getLocus().startsWith("Contam") )
  //              continue;

            if(protein.getPeptideSize()<=0)
            {
                aList.add(protein);
                set.add(protein.getLocus());
            }
            else
            {
//		Element proteinEle = new Element("protein");

//		proteinEle.setAttribute("protein_name", protein.getLocus());
/*
		for(Iterator<Protein> proItr=aList.iterator(); proItr.hasNext(); )
		{
		    Protein eachPro = proItr.next(); 
		    Element indistProEle = new Element("indistinguishable_protein");
		    indistProEle.setAttribute("protein_name", eachPro.getLocus());
		    proteinEle.addContent(indistProEle);
		}
*/
		for(Iterator<Peptide> pepItr=protein.getPeptides(); pepItr.hasNext(); )
		{
		    peptide = pepItr.next();
		    //		if(protein.getLocus().equals("IPI00430839.1"))
		    Element spectrumQueryEle = new Element("spectrum_query");
		    Element searchResultEle = new Element("search_result");
		    searchResultEle.setAttribute("assumed_charge", peptide.getChargeState());
		    Element searchHitEle = new Element("search_hit");
		    searchHitEle.setAttribute("hit_rank", "1");
		    searchHitEle.setAttribute("protein", protein.getLocus());
		    String seq = peptide.getSequence();
		    searchHitEle.setAttribute("peptide", seq.substring(2, seq.length()-2));
		    searchHitEle.setAttribute("peptide_prev_aa", seq.substring(0,1));
		    searchHitEle.setAttribute("peptide_next_aa", seq.substring(seq.length()-1, seq.length()));

		    Element score = new Element("search_score");
		    score.setAttribute("xcorr", peptide.getXCorr());
		    searchHitEle.addContent(score);

		    score = new Element("search_score");
		    score.setAttribute("deltacn", peptide.getDeltCN());
		    searchHitEle.addContent(score);
		    searchResultEle.addContent(searchHitEle);
		    spectrumQueryEle.setAttribute("spectrum", peptide.getFileName());
		    spectrumQueryEle.setAttribute("start_scan", peptide.getScanNum());
		    spectrumQueryEle.setAttribute("end_scan", peptide.getScanNum());
		    spectrumQueryEle.addContent(searchResultEle);		    
		    rootEle.addContent(spectrumQueryEle);

		    i++;
		}

                aList.clear();
            }
        }
//	System.out.println(i);

	Document doc = new Document(rootEle);
	OutputStream os = new FileOutputStream(new File("converted_pepxml.xml")); //(filePath + "census_chro.xml");
	XMLOutputter outputter = new XMLOutputter();
	outputter.setFormat(Format.getPrettyFormat());
	outputter.output(doc, os);
	os.close();

    }

}

