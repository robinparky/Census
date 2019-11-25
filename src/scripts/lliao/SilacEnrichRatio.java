package scripts.lliao;

import edu.scripps.pms.census.io.*;
import java.io.*;
import java.util.*;
import edu.scripps.pms.census.model.*;
import edu.scripps.pms.census.util.io.FastaReader;
import edu.scripps.pms.util.seq.Fasta;

import java.text.DecimalFormat;
import java.text.NumberFormat;

public class SilacEnrichRatio
{
    private static DecimalFormat format = new DecimalFormat("0.000");

    public static void main(String[] args) throws IOException
    {
	if(args.length<2)
	{
	    System.out.println("java /home/rpark/rpark_on_data/project/lliao/silac_enrichment_correct/test/peptide_merged_01.txt /data/6/lliao/neuron/SILAC/060228SILAC/mixed/mix_Ratio1_May07.txt");
	    return;
	}
	
        BufferedReader br = null;
        String eachLine;
        
	//br = new BufferedReader(new FileReader(args[0]));
	//read reference
	//br = new BufferedReader(new FileReader("/data/6/lliao/neuron/SILAC/enrichment060723/enrich060723-new.txt"));
	//br = new BufferedReader(new FileReader("/home/rpark/rpark_on_data/project/lliao/silac_enrichment_correct/test/peptide_merged_01.txt"));
	br = new BufferedReader(new FileReader(args[0])); //"/home/rpark/rpark_on_data/project/lliao/silac_enrichment_correct/test/peptide_merged_01.txt"));

	Hashtable<String,ChroPeptide> ht = new Hashtable<String, ChroPeptide>();

//	Set pepSet = new HashSet();

	while( (eachLine = br.readLine())!=null )
	{
	    //		if( eachLine.startsWith("H\t") )
	    //		    System.out.println(eachLine);

	    if(eachLine.startsWith("S\t"))
	    {
		String[] arr = eachLine.split("\t");
		String seq = arr[2];
		ChroPeptide pep = new ChroPeptide();
		pep.setSequence(seq);
		pep.setRatio( Double.parseDouble(arr[3]) );
		pep.setCorr( Double.parseDouble(arr[4]) );

		ChroPeptide tmpPep = ht.get(seq);
		if(null == tmpPep)
		{
		    ht.put(arr[2], pep);
		}
		else
		{
		    double corr = tmpPep.getCorr();
		    if(corr<Double.parseDouble(arr[4]))
		    {
			ht.put(arr[2], pep);
		    }
		}
		

//		System.out.println(arr[2] + " " + arr[3]);
	    }
	}

	br.close();


/*
System.out.println(ht.size());

*/
int size=ht.size();
double totalRatio=0;
for(Iterator<String> itr=ht.keySet().iterator(); itr.hasNext(); )
{
    String seq =  itr.next();
    ChroPeptide pep = ht.get(seq);

    totalRatio += pep.getRatio();

    
    //System.out.println(pep.getSequence() + " " + pep.getCorr());
}
double ratioEff = (totalRatio/size);
System.out.println("H\tratio efficiency\t" + ratioEff);


	//read mixed experimental file
	//br = new BufferedReader(new FileReader("/data/6/lliao/neuron/SILAC/060723WT-H_WT-L_lowEnrich/WT-H_WT-L_merge2.txt"));
	br = new BufferedReader(new FileReader(args[1])); //"/data/6/lliao/neuron/SILAC/060228SILAC/mixed/mix_Ratio1_May07.txt"));

	int pepCount=0;
	StringBuffer sb = new StringBuffer();
	while( (eachLine = br.readLine())!=null )
	{
	    //		if( eachLine.startsWith("H\t") )
	    //		    System.out.println(eachLine);

	    if(eachLine.startsWith("S\t"))
	    {
		pepCount++;
		String[] arr = eachLine.split("\t");

		ChroPeptide pep = ht.get(arr[2]);


		if(null == pep){
		    sb.append("S\tN/A\t" + arr[2] + "\t" + arr[3] + "\t");
		    

		    double r1 = ratioEff;
		    double r2 = Double.parseDouble(arr[3]);

		    //r = (r2-r1)/(r1+1)
		    sb.append( format.format( (r2-r1)/(r1+1) ) );
		    //sb.append( "N/A" );
		}
		else
		{
		    sb.append("S\t\t" + arr[2] + "\t" + arr[3] + "\t");
		    double r1 = pep.getRatio();
		    double r2 = Double.parseDouble(arr[3]);

		    //r = (r2-r1)/(r1+1)
		    sb.append( format.format( (r2-r1)/(r1+1) ) );

		}

		sb.append("\t" + arr[4] + "\t" + arr[5] + "\t" + arr[6] + "\t" + arr[7] + "\t" + arr[8] );
		sb.append("\n");
		System.out.print(sb.toString());
		sb = new StringBuffer();
	    }
	    /*
	    if(eachLine.startsWith("P\t"))
	    {
		System.out.println(eachLine);

		pepCount=0;
		sb = new StringBuffer();

	    }*/
	    else
	    {
		System.out.println(eachLine);
	    }

	}

	br.close();

    }
}
