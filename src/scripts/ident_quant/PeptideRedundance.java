package scripts.ident_quant;

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

public class PeptideRedundance
{
    public static void main(String[] args) throws IOException
    {
	if(args.length<1)
	{
	    System.out.println("Check redundance peptides on different salt steps");
	    System.out.println("Usage: PeptideRedundance DTASelect-filter.txt");
	    return;
	}


        DTASelectFilterReader reader = new DTASelectFilterReader( args[0] );

        Iterator<Protein> pitr = reader.getProteins();

        Protein protein;
	Peptide peptide;
        ArrayList<Protein> aList = new ArrayList<Protein>();

	Hashtable<String, Set> ht = new Hashtable<String, Set>();

        for (Iterator<Protein> itr = pitr; itr.hasNext(); )
        {
            protein = itr.next();

	    for(Iterator<Peptide> pepItr=protein.getPeptides(); pepItr.hasNext(); )
	    {
		peptide = pepItr.next();

		String seq = peptide.getSequence();
		String fName = peptide.getFileName();

//		System.out.println(fName);
//		fName = fName.substring(0, fName.indexOf("."));
		//System.out.println(fName.substring(0,fName.indexOf(".")));

		Set set = ht.get(seq);

		if(null == set)
		{
		    set = new HashSet();
		    set.add(fName);
		    ht.put(seq, set);
		}
		else
		{
		    set.add(fName);
		}

	    }

        }

	Set keySet = ht.keySet();

	//System.out.println(ht);
	for(Iterator<String> itr=keySet.iterator(); itr.hasNext(); )
	{
	    String key = itr.next();

	    Set set = ht.get(key);

	    if(set.size()>1)
		System.out.println(set);

	}

    }
}
