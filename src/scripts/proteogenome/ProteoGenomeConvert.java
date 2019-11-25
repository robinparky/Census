package scripts.proteogenome;

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
import edu.scripps.pms.census.util.io.FastaReader;
import edu.scripps.pms.util.seq.Fasta;

public class ProteoGenomeConvert
{
    public static void main(String[] args) throws Exception
    {
	if(args.length<1)
	{
	    System.out.println("Convert protein information to genome representation for genome browser");
	    System.out.println("Usage: ProteoGenomeConvert DTASelect-filter.txt");
	    return;
	}

	Hashtable<String, Integer> chromHt = readGenomeDatabaseIndex();
	Hashtable<String, String> proteinDbHt= readProteinDatabase();

	//args[0] = /data/4/cbamberg/nemas/070108_gemone_analysis/0_100/042208detailed/DTASelect-filter.txt 
	args[0] = "/data/4/cbamberg/nemas/070108_gemone_analysis/0_100/042208detailed/DTASelect-filter.txt";
	BufferedReader br = new BufferedReader(new FileReader(args[0]));
	String eachLine = "";
	StringBuffer summarySb = new StringBuffer();
	boolean readSummary = false;
/*
	while( null != (eachLine = br.readLine()) )
	{
	    System.out.println(eachLine);

	    if(eachLine.startsWith("Unique"))
		break;
	}
*/
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
	    String accession = protein.getLocus();
	    if(accession.startsWith("Conta") || accession.startsWith("Rever"))
		continue;	    

	    String proteinSeq = proteinDbHt.get(accession);

	    int accIndex = accession.indexOf("Frame");
	    String scaffoldName = accession.substring(0, accIndex);
	    
	    scaffoldName = scaffoldName.replace("scfd", "scaffold");

	    int scaffoldSize = chromHt.get(scaffoldName);

	    String strand = accession.substring(accession.indexOf("Frame"));
	    String[] tmpArr = accession.split("_");
	    String[] indexArr = tmpArr[3].split("-");

	    int browSide = 30;

	    int strandFrame = (int)strand.charAt(7) - (int)'0';

	    char strandCh = strand.charAt(6);

  
	    String SequenceCount = protein.getSeqCount();
//	    System.out.println("SequenceCount = " + SequenceCount);

	    int blockCount = Integer.parseInt(SequenceCount); 
 
//	    System.out.println("SeqCount = " + SeqCount);
	    
	    int browStart = 0;
	    int browEnd = 0;

//	    if(strandCh == '+') {
//		browStart = Integer.parseInt(indexArr[0])*3 - 2 + strandFrame;
		browStart = Integer.parseInt(indexArr[0])*3;
		browEnd = Integer.parseInt(indexArr[1])*3;

//	    } else
	    if(strandCh == '-') {
		//System.out.println(indexArr[1]);
	
		int tmp_browStart = browStart;
		int tmp_browEnd = browEnd;
	
		browStart = scaffoldSize - tmp_browEnd;
		browEnd = scaffoldSize - tmp_browStart;

	    }

	    System.out.println("browser position " + scaffoldName + ":" + (browStart - browSide) + "-" + (browEnd + browSide));
	    System.out.println("track name=\"" + accession + "\" description=\"peptides in " + accession + "\" visibility=2 color=240,120,0 useScore=1");
	    
	    int count = -1;
	    int[] blockStarts_arr = new int[blockCount];
	    int[] blockSizes_arr = new int[blockCount];
	    
	    for(Iterator<Peptide> pepItr=protein.getPeptides(); pepItr.hasNext(); )
	    	    	    
	    {

		count ++;
		
//		System.out.println("counter = " + count);
		
		peptide = pepItr.next();

		String[] arr = peptide.getPeptideLine();
		String pepSeq = arr[12].substring(2, arr[12].length()-2);
		
		int pepIndex = proteinSeq.indexOf(pepSeq); 
		int pepLength = pepSeq.length();



		    int chromStart = Integer.parseInt(indexArr[0]) *3 + pepIndex *3 + strandFrame;
		    int chromEnd = Integer.parseInt(indexArr[0]) *3 + pepIndex *3 + strandFrame + pepLength *3;

		    blockStarts_arr[count] = pepIndex *3 + strandFrame;
		    blockSizes_arr[count] = pepLength *3;
		    

 	   

		 if(strandCh == '-') {
		    
		    int strandFrame1 = 0;
	    
		    if (strandFrame == 0) {strandFrame1 = 1;};
		    if (strandFrame == 1) {strandFrame1 = 2;};
		    if (strandFrame == 2) {strandFrame1 = 0;};
  
		    chromStart = scaffoldSize + 2 + strandFrame1 - ( Integer.parseInt(indexArr[0]) *3 + pepIndex *3 + pepLength *3);
		    chromEnd = scaffoldSize + 1 + strandFrame1 - ( Integer.parseInt(indexArr[0]) *3 + pepIndex *3 );

    		}


		System.out.println(scaffoldName + "\t" + chromStart + "\t" + chromEnd + "\t" + pepSeq + "\t500\t" + strandCh); 
	

		//String tmp = arr[1].substring(0, arr[1].indexOf("."));
		//tmp += arr[11];
/*
		if(set.contains(tmp))
		{

		}
		else
		{
		    for(String each:arr)
		    {
			// System.out.print(each);
			// System.out.print("\t");
		    }
		    //System.out.println("");
		    set.add(tmp);
		}

*/
	      }

	    System.out.println(" ");

	    System.out.print(scaffoldName + "\t" + browStart+ "\t" + browEnd + "\t" +  accession+ "\t500\t" + strandCh + "\t" + browStart + "\t" + "\t" + browEnd + "\t" + "\t240,240,240\t" + blockCount + "\t");

		for(int each:blockSizes_arr) 
//		for(int i=0;i<blockCount-1;i++) 
	    	   System.out.print(each + ",");
//			System.out.print(blockSizes_arr[blockCount-1]);
	    
	    System.out.print("\t");
	    
                for(int each:blockStarts_arr)
//		for(int i=0;i<blockCount-1;i++)
		   System.out.print(each + ",");
//		        System.out.print(blockStarts_arr[blockCount-1]);
									    
	    System.out.println(" ");
	    System.out.println(" ");

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

//	System.out.println(summarySb.toString());
    }

    public static Hashtable<String, Integer> readGenomeDatabaseIndex() throws Exception {
	String fileName = "/data/2/cbamberg/databases/Nematostella/scaffold_size.txt";
	BufferedReader br = new BufferedReader(new FileReader(fileName));
	String eachLine;

	Hashtable<String, Integer> chromHt = new Hashtable<String, Integer>();

	while( (eachLine=br.readLine()) != null) {
	    String[] arr = eachLine.split("\t");
	    chromHt.put(arr[0], Integer.parseInt(arr[1]));


	}

	return chromHt;

    }

    public static Hashtable<String, String> readProteinDatabase() throws Exception {
	String fileName = "/data/2/cbamberg/databases/Nematostella/Nemve1_6frame_ORFs_contamin_reverse.fasta";

	Fasta fasta;
	String defLine;
	Hashtable<String, String> genomeDbHt = new Hashtable<String, String>();
	System.out.print("reading genome db...");
	for (Iterator itr = FastaReader.getFastas(new FileInputStream(fileName)); itr.hasNext(); ) {
	    fasta = (Fasta) itr.next();
	    defLine = fasta.getDefline();
//	    System.out.println(fasta.getAccession() + "\t" + fasta.getSequence().length());
	    
	    genomeDbHt.put(fasta.getAccession(), fasta.getSequence());
	}
	System.out.println("done");

	return genomeDbHt;
    }
    public static Hashtable<String, String> readGenomeDatabase() throws Exception {
	String fileName = "/data/2/cbamberg/databases/Nematostella/Nemve1.fasta";

	Fasta fasta;
	String defLine;
	Hashtable<String, String> genomeDbHt = new Hashtable<String, String>();
	System.out.print("reading genome db");
	for (Iterator itr = FastaReader.getFastas(new FileInputStream(fileName)); itr.hasNext(); ) {
	    fasta = (Fasta) itr.next();
	    defLine = fasta.getDefline();

//	    System.out.println(fasta.getAccession() + "\t" + fasta.getSequence().length());
	    System.out.print(".");
	    
	    genomeDbHt.put(fasta.getAccession(), fasta.getSequence());
	}
	System.out.println("done");

	return genomeDbHt;
    }
}
