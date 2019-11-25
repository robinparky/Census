/*
 * ChroReader.java
 *
 * Created on May 17, 2005, 12:00 PM
 */

import java.io.*;
import java.util.*;

import org.jdom.*;
import org.jdom.output.*;
import org.jdom.input.*;

import java.nio.*;

/**
 *
 * @author  Robin Park
 * @version $Id: JDomEnrichment.java,v 1.1 2009/03/13 21:55:58 rpark Exp $
 */

public class JDomEnrichment {

    public static void main(String args[]) throws Exception
    {

        //File file = new File("/home/rpark/rpark_on_data/project/census/n15_enrichment_calc/lliao_non_phospho/p1p45BrainNuc1/census_chro.xml");
        File file = new File(args[0]);
	SAXBuilder builder = new SAXBuilder();

        Document doc = builder.build( file );
        Element rootEle = doc.getRootElement();

	for(Iterator<Element> itr=rootEle.getChildren("protein").iterator(); itr.hasNext(); )
	{
	    Element proEle = itr.next();

	    for(Iterator<Element> pepitr=proEle.getChildren("peptide").iterator(); pepitr.hasNext(); )
	    {
		Element pepEle = pepitr.next();
		System.out.print(pepEle.getAttributeValue("seq"));
		System.out.print("\t");
		System.out.print(pepEle.getAttributeValue("enrichment"));
		System.out.print("\t");
		System.out.print(pepEle.getAttributeValue("enrich_corr"));
		System.out.println("");
	    }
	}

    
    }    


}


