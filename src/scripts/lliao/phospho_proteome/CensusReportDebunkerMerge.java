
/*
 * CenSusReportReader.java
 *
 * Created on February 23, 2006, 3:03 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package scripts.lliao.phospho_proteome;

import java.io.*;
import java.util.*;
import edu.scripps.pms.census.model.*;
import edu.scripps.pms.census.util.*;
import edu.scripps.pms.census.util.io.FastaReader;
import edu.scripps.pms.util.seq.Fasta;

import gnu.trove.TDoubleArrayList;
import scripts.lliao.*;
/**
 *
 * @author rpark
 * @version $Id: CensusReportDebunkerMerge.java,v 1.4 2014/08/27 18:00:35 rpark Exp $
 */
public class CensusReportDebunkerMerge {
    
    private int totalPeptideCount=0;
    private static double threshold=0;

    public static Hashtable<String, String> getDebunkerOutput(String fileName) throws Exception
    {
	String lastLine;
	BufferedReader br = new BufferedReader(new FileReader(fileName));

	Hashtable<String, String> ht = new Hashtable<String, String>();
	
	while( null != (lastLine=br.readLine()) )
	{
	    String[] arr = lastLine.split("\t");
	    ht.put(arr[0] + arr[1], arr[2] + "\t" + arr[3] + "\t" + arr[4]);
	}

	return ht;

    }
    
    public static void main(String args[]) throws Exception
    {
	if(args.length<3)
	{
	    
	    System.out.println("Usage: ");
	    System.out.println("1st arg is census-merge.txt");
	    System.out.println("2nd arg is debunker-out.txt");
	    System.out.println("3rd arg is debunker threshold. 0.5 is recommended.");
	    
	    System.exit(0);
	}
	
	threshold = Double.parseDouble(args[2]);

	Hashtable<String, String> ht = getDebunkerOutput(args[1]);

        List<MergeProteinModel> list = Collections.EMPTY_LIST;
	CensusReportReaderModified re = new CensusReportReaderModified();
        list = re.read(args[0]);
       
        for(Iterator<MergeProteinModel> itr=list.iterator(); itr.hasNext(); )
        {
            MergeProteinModel protein = itr.next();
            
            List<ChroProtein> pList = protein.getProteins();

/*
	    double ratioSum = 0;
	    int pepSize = 0;
	    
            for(Iterator<MergeProteinModel.Peptide> pepItr=protein.getPeptides().iterator(); pepItr.hasNext(); )
            {
                MergeProteinModel.Peptide pep = pepItr.next();
		ratioSum += pep.getAverageRatio();
		pepSize++;
            }
	    */
            
	    StringBuffer pepSb = new StringBuffer();

	    int pepSize=0;
	    double ratioSum=0;
	    TDoubleArrayList dlist = new TDoubleArrayList();
            for(Iterator<MergeProteinModel.Peptide> pepItr=protein.getPeptides().iterator(); pepItr.hasNext(); )
            {
                MergeProteinModel.Peptide pep = pepItr.next();
		String key = pep.getFileName() + pep.getSequence();
		String value = ht.get(key);

		if(null != value) // && d>0.5)
		{
		    String[] debunkArr = value.split("\t");
		    double d = Double.parseDouble(debunkArr[2]);

		    if(d<=threshold)
			continue;

		    pepSize++;
		    ratioSum += pep.getRatio();
		    dlist.add(pep.getRatio());
		   
		    pepSb.append(pep.getPeptideLine()).append("\t").append(value).append("\n");
		}
		    

            }

	    if(pepSb.length()>0)
	    {
            for(Iterator<String> pItr=protein.getProteinNames().iterator(); pItr.hasNext(); )
            {
                String cPro = pItr.next();
                
                System.out.print("P\t" + cPro + "\t" + CensusHelper.format.format(ratioSum/pepSize) + "\t");
                System.out.println( edu.scripps.pms.util.stats.StatCalc.getStandardDeviation( dlist.toNativeArray() ) + "\t" + "0.0\t" + pepSize + "\t0\t" + protein.getProteinDesc(cPro));
            }
           
	    System.out.print(pepSb.toString());
	    }
        }
        
//	CensusReportReaderModified reader = new CensusReportReaderModified();
	//reader.calculateDeviationError(args[0]);

    }
}
