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

public class DTASelectProtXMLConvert
{
    public static void main(String args[]) throws Exception
    {
	//Convert DTASelect to pepXML
	if(args.length<1)
	{
	    System.out.println("Usage: java DTASelectProtXMLConvert DTASelect-filter.txt");
	    System.exit(0);
	}

	//build jdom
	Element rootEle = new Element("protein_summary");

	Element proteinGroupEle = new Element("protein_group");


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

            System.out.println(protein.getLocus());
            if( protein.getLocus().startsWith("Rever") || protein.getLocus().startsWith("contam") || protein.getLocus().startsWith("Contam") )
                continue;


            if(protein.getPeptideSize()<=0)
            {
                aList.add(protein);
                set.add(protein.getLocus());
            }
            else
            {
		Element proteinEle = new Element("protein");

		proteinEle.setAttribute("protein_name", protein.getLocus());

		for(Iterator<Protein> proItr=aList.iterator(); proItr.hasNext(); )
		{
		    Protein eachPro = proItr.next(); 
		    Element indistProEle = new Element("indistinguishable_protein");
		    indistProEle.setAttribute("protein_name", eachPro.getLocus());
		    proteinEle.addContent(indistProEle);
		}

		for(Iterator<Peptide> pepItr=protein.getPeptides(); pepItr.hasNext(); )
		{
		    peptide = pepItr.next();
		    //		if(protein.getLocus().equals("IPI00430839.1"))
		    Element peptideEle = new Element("peptide");
		    peptideEle.setAttribute("peptide_sequence", peptide.getSequence());
		    peptideEle.setAttribute("charge", peptide.getChargeState());

		    proteinEle.addContent(peptideEle);

		    System.out.println(peptide.getSequence());
		    i++;

		}

		rootEle.addContent(proteinEle);
		
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

