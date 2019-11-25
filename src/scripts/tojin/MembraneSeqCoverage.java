package scripts.tojin;

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

public class MembraneSeqCoverage 
{
    public static void main(String args[]) throws Exception
    {
        DTASelectFilterReader reader = new DTASelectFilterReader( args[0] );
	String proteinSeq = "MKHHHHHHHHGGLVPRGSHGMFCTFFEKHHRKWDILLEKSTGVMEAMKVTSEEKEQLSTAIDRMNEGLDAFIQLYNESEIDEPLIQLDDDTAELMKQARDMYGQEKLNEKLNTIIKQILSISVSEEGEKEGSGSGHMHHHHHHSSGGSSTSLYKKAGSLVPRGSGSVMKTLFSRLITVIACFFIFSAAWFCLWSISLHLVERPDMAVLLFPFGLRLGLMLQCPRGYWPVLLGAEWLLIYWLTQAVGLTHFPLLMIGSLLTLLPVALISRYRHQRDWRTLLLQGAALTAAALLQSLPWLWHGKESWNALLLTLTGGLTLAPICLVFWHYLANNTWLPLGPSLVSQPINWRGRHLVWYLLLFVISLWLQLGLPDELSRFTPFCLALPIIALAWHYGWQGALIATLMNAIALIASQTWRDHPVDLLLSLLVQSLTGLLLGAGIQRLRELNQSLQKELARNQHLAERLLETEESVRRDVARELHDDIGQTITAIRTQAGIVQRLAADNASVKQSGQLIEQLSLGVYDAVRRLLGRLRPRQLDDLTLEQAIRSLMREMELEGRGIVSHLEWRIDESALSENQRVTLFRVCQEGLNNIVKHADASAVTLQGWQQDERLMLVIEDDGSGLPPGSGQQGFGLTGMRERVTALGGTLHISCLHGTRVSVSLPQRYV";

	int[] countArr = new int[proteinSeq.length()];


        Iterator<Protein> pitr = reader.getProteins();

        Protein protein;
	Peptide peptide;
        ArrayList<Protein> aList = new ArrayList<Protein>();

	HashSet set = new HashSet();

        for (Iterator<Protein> itr = pitr; itr.hasNext(); )
        {
            protein = itr.next();

//	    if(!protein.getLocus().equals("NP_418124.4"))
//		continue;

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

		    //peptide.getChargeState();
		    //protein.getLocus();


		    String[] arr = peptide.getPeptideLine();
		    int redundance = Integer.parseInt(arr[arr.length-3]);
		    double d = Double.parseDouble(arr[arr.length-5]);

		    //0 for not using KD
		    //1 for using KD
		    int cut = Integer.parseInt(args[1]);
		    if(cut!=0 && d<15)
			continue;

		    String seq = arr[arr.length-2];
		    seq = seq.substring(2, seq.length()-2);
		    seq = seq.replaceAll("\\*", "");
		    seq = seq.replaceAll("\\#", "");
		    int start = proteinSeq.indexOf(seq);
		    int end = start + seq.length()-1;

		    for(int rc=0;rc<redundance;rc++)
			for(int i=start;i<=end;i++)
			{
			    countArr[i]++;
			}
		    
		   // System.out.println(seq);
		   // System.out.println(start + " " + end);

		}


	    }
	}

	char[] charr = proteinSeq.toCharArray();

	for(int i=0;i<charr.length;i++)
	    System.out.println(charr[i] + "\t" + countArr[i]);
	    

    }


/*
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
*/
}

