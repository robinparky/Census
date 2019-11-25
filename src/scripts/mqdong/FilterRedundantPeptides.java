package scripts.mqdong;

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

public class FilterRedundantPeptides
{
    public static void main(String[] args) throws IOException
    {
	if(args.length<1)
	{
	    System.out.println("Remove all redundant peptides from DTASelect-filter.txt for census");
	    System.out.println("Usage: FilterRedundantPeptides DTASelect-filter.txt");
	    return;
	}

	BufferedReader br = new BufferedReader(new FileReader(args[0]));
	String eachLine = "";
	StringBuffer summarySb = new StringBuffer();
	boolean readSummary = false;

	while( null != (eachLine = br.readLine()) )
	{
	    System.out.println(eachLine);

	    if(eachLine.startsWith("Unique"))
		break;
	}

	while( null != (eachLine = br.readLine()) )
	{
	    if(eachLine.startsWith("\tProteins\tPeptide") || readSummary)
	    {
		summarySb.append(eachLine).append("\n");
		readSummary=true;
	    }
	}

        DTASelectFilterReader reader = new DTASelectFilterReader( args[0] );

        Iterator<Protein> pitr = reader.getProteins();

        Protein protein;
	Peptide peptide;
        ArrayList<Protein> aList = new ArrayList<Protein>();

	Hashtable<String, String> ht = new Hashtable<String, String>();

        for (Iterator<Protein> itr = pitr; itr.hasNext(); )
        {
            protein = itr.next();

	    HashSet set = new HashSet();
            System.out.println(protein.getProteinLine());

                for(Iterator<Peptide> pepItr=protein.getPeptides(); pepItr.hasNext(); )
                {
                    peptide = pepItr.next();

		    String[] arr = peptide.getPeptideLine();

		    String tmp = arr[1].substring(0, arr[1].indexOf("."));
		    tmp += arr[11];

		    if(set.contains(tmp))
		    {

		    }
		    else
		    {
			for(String each:arr)
			{
			    System.out.print(each);
			    System.out.print("\t");
			}
			System.out.println("");
			set.add(tmp);
		    }


                }

            if(protein.getPeptideSize()<=0)
            {
                aList.add(protein);
                set.add(protein.getLocus());
            }
            else
            {
                aList.clear();


            }
        }

	System.out.println(summarySb.toString());
    }
}
