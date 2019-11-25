package scripts.emily.transmembrane;

import edu.scripps.pms.census.util.io.FastaReader;
import java.util.*;
import java.io.*;

import edu.scripps.pms.util.seq.Fasta;
import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Iterator;
import java.util.List;
import java.util.LinkedList;
import java.io.IOException;
import java.io.FileReader;
import java.io.FileInputStream;
import java.io.File;
import edu.scripps.pms.census.util.dtaselect.Protein;
import java.util.ArrayList;
import edu.scripps.pms.census.util.dtaselect.Peptide;
import edu.scripps.pms.census.util.io.DTASelectFilterReader;
import edu.scripps.pms.census.util.dtaselect.Protein;
import edu.scripps.pms.census.util.dtaselect.Peptide;


public class FindTransMembraneProteins
{
    public static void main(String[] args) throws Exception
    {
	Fasta fasta;
	String defLine;
	Hashtable<String, String> fastaHt = new Hashtable<String, String>();

	if(args.length<3)
	{
	    System.out.println("java FindTransMembraneProteins fasta_db TMHMM_input DTASelect-filter.txt");
	    System.exit(0);
	}

	//args[0] : fasta database file 
	//args[1] : TMHMM input file
	//args[2] : DTASelect-filter.txt

	//*** Build fasta hash table ***/
	//args[0] = /home/rpark/rpark_on_data/project/emily/trans_membrane/Contrast.fasta
	//String fastaFile = "/home/rpark/rpark_on_data/project/emily/trans_membrane/Contrast.fasta";
	String fastaFile = args[0];
	for (Iterator itr = FastaReader.getFastas(new FileInputStream(fastaFile)); itr.hasNext(); ) {
	    fasta = (Fasta) itr.next();
	    defLine = fasta.getDefline();
	
	    fastaHt.put(fasta.getAccession(), fasta.getSequence());
	}

//	System.out.println(fastaHt);

	//if(true)
	//return;

	//*** transmembrane hashtable ***/
	//String TMHMMInput = "/home/rpark/rpark_on_data/project/emily/trans_membrane/All_rat_mem.txt";
	String TMHMMInput = args[1];
	BufferedReader br = new BufferedReader(new FileReader(TMHMMInput));
	String each;
	Hashtable<String, List> tmHt = new Hashtable<String, List>();

	while( null != (each = br.readLine()) )
	{
	    if(!each.startsWith("gi"))
	    //if(!each.startsWith("IPI"))
		continue;

	    //String acc = each.substring( each.indexOf("_")+1 );
	    //acc = acc.substring(0, acc.indexOf("_"));

	    String acc=null;
	    if(each.contains("NP"))
		acc = each.substring( each.indexOf("NP") );
	    else if(each.contains("YP"))
		acc = each.substring( each.indexOf("YP") );
	    else
	    {
		System.out.println("error");
		System.exit(0);
	    }

	    acc = acc.substring(0, acc.lastIndexOf("_"));
	    String seq = fastaHt.get(acc);

	    if(seq == null)
	    {
		System.out.println("Error: no protein found");
		continue;
	    }

	    String[] arr = each.split("\t");
	  
	    if("TMhelix".equals(arr[2]))
	    {
		List l = tmHt.get(acc);

		if(null == l)
		{
		    l = new ArrayList();
		    l.add(arr[3].trim());
		    tmHt.put(acc, l);
		}
		else
		{
		    l.add(arr[3].trim());
		}
		
		//System.out.println(arr[2] + "\t" + arr[3] );  //+ "\t" + arr[4]);
	    }
	}

	//read dtaselect file
	String dtaSelectFile = args[2]; //"/home/rpark/rpark_on_data/project/emily/trans_membrane/Merged_All-All-filter.txt";
	DTASelectFilterReader reader = new DTASelectFilterReader( dtaSelectFile );

	Iterator<Protein> pitr = reader.getProteins();

	Protein protein;
	Peptide peptide;
	ArrayList<Protein> aList = new ArrayList<Protein>();

	HashSet<String> set = new HashSet<String>();

	//print header
	System.out.println("Accession\tPepSequence\tTMPepSeq\tChargeState\tOverlabSeq\tTMStartIndex\tTMEndIndex\tPeptideStartIndex\tPeptideEndIndex\tKD\tScanNum");
	int i=0;
	for (Iterator<Protein> itr = pitr; itr.hasNext(); )
	{
	    protein = itr.next();

	    if( protein.getLocus().startsWith("Rever") || protein.getLocus().startsWith("contam") || protein.getLocus().startsWith("Contam") )
		continue;

	    set.add(protein.getLocus());

	    for(Iterator<Peptide> pepItr=protein.getPeptides(); pepItr.hasNext(); )
	    {
//		System.out.println(set);
		peptide = pepItr.next();
		String pepSeq = peptide.getSequence();

		pepSeq = pepSeq.replaceAll("#","");
		pepSeq = pepSeq.replaceAll("\\*","");
		pepSeq = pepSeq.replaceAll("@","");

		pepSeq = pepSeq.substring(2, pepSeq.length()-2);
		for(Iterator<String> proAccitr=set.iterator(); proAccitr.hasNext(); )
		{
		    String proAcc = proAccitr.next();
		    String proSeq = fastaHt.get(proAcc);

		    List ll = tmHt.get(proAcc);

		    if(null == ll)
		    {
			continue;
		    }

	//	    System.out.println("");
	//	    System.out.println(proAcc + "\t" + peptide.getSequence() + "\t" + peptide.getFileName());
		    int pepStartIndex = proSeq.indexOf(pepSeq) + 1;
		    int pepEndIndex = pepStartIndex + pepSeq.length() - 1;

		    int tmstart = -1;
		    int tmend = -1;
		    int overlap = -1;
		    for(Iterator<String> indexItr=ll.iterator(); indexItr.hasNext(); )
		    {
			String[] indexArr = indexItr.next().split(" ");

			tmstart = Integer.parseInt(indexArr[0]);
			tmend = Integer.parseInt(indexArr[indexArr.length-1]);

			if(pepEndIndex<tmstart || tmend<pepStartIndex) //no overlap
			    continue;

			if(pepStartIndex>=tmstart && pepEndIndex<=tmend)
			{
			    overlap = pepEndIndex - pepStartIndex + 1;
			    break;
			}

			if(pepStartIndex<=tmstart && pepEndIndex>=tmend)
			    overlap = tmend - tmstart + 1;
			    /*
			{
			    overlap = pepEndIndex - pepStartIndex + 1;
			    break;
			}*/

			
			if(pepStartIndex<tmstart)
			{
			    if(pepEndIndex<tmend)
				overlap = pepEndIndex-tmstart+1;
			    else
				overlap = tmend - tmstart + 1;


			    break;
			}
			else
			{
			    if(pepEndIndex>tmend)
				overlap = tmend - pepStartIndex + 1;
			    else
				overlap = pepEndIndex - pepStartIndex + 1;

			    break;
			}

		    }

		    System.out.println(proAcc + "\t" + pepSeq + "\t"  + proSeq.substring(tmstart-1, tmend) + "\t" + peptide.getChargeState() + "\t" + ((overlap<0)?"0":overlap) + "\t" + tmstart  + "\t" + tmend + "\t" + pepStartIndex + "\t" + pepEndIndex + "\t" + peptide.getFileName() + "\t" + peptide.getKd() + "\t" + peptide.getScanNum());

		    
		}
		
		//          if(protein.getLocus().equals("IPI00430839.1"))
		//             System.out.println(peptide.getSequence());

	    }

	    if(protein.getPeptideSize()<=0)
	    {
	    }
	    else
	    {
		set.clear();
	    }
	}
	//              System.out.println(set.size());
	
    }
}
