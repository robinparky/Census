/*
 * PepSniffer.java
 *
 * Created on Jan 18, 2006, 3:52 PM
 */

package edu.scripps.pms.pepsniffer;

import edu.scripps.pms.census.util.io.SpectrumReader;

import edu.scripps.pms.util.spectrum.*;
import java.util.*;
import java.io.*;

/**
 *
 * @author  Robin Park
 * @version $Id: Census.java,v 1.2 2007/01/03 06:01:44 rpark Exp $
 */

public class PepSniffer {

    /** Creates new form inRelaxMainFrame */
    public static void main(String args[]) throws IOException, Exception
    {
	System.out.println( System.currentTimeMillis() );

	SpectrumReader reader = new SpectrumReader(args[0], "ms1");

	String line = null;
	SpectrumReader sr = new SpectrumReader(args[0], args[1]);
	Hline h = new Hline(sr.getHlines());

	Iterator<PeakList> it = sr.getSpectra();
	int counter = 0;
	int numPeaks = 0;
	//boolean sortByIntensity = true;
	while (it.hasNext()) {
	    PeakList list = it.next();
	    Peak p;
	    StringBuffer sb = new StringBuffer();

	    for(Iterator<Peak> itr=list.getPeaks(); itr.hasNext(); )
	    {
		p = itr.next();
		sb.append(p.getM2z());
System.out.println(p.getM2z() + " " + p.getIntensity());
		sb.append("\t");
		sb.append(p.getIntensity());
		sb.append("\n");
	    }
	    for(Iterator itr = list.getZlines(); itr.hasNext(); )
	    {
		Zline  zline = (Zline)itr.next();

		System.out.println(zline.getM2z() + " " + zline.getChargeState());

	    }

break;
	}
    }

    public PepSniffer() throws IOException, Exception
    {
    }

}
