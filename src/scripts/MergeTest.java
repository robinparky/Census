package scripts;

import org.jdom.input.SAXBuilder;
import org.jdom.*;
import java.io.*;
import java.util.*;
import edu.scripps.pms.census.model.*;

import edu.scripps.pms.census.util.io.FastaReader;
import edu.scripps.pms.census.io.CenSusReportReader;
import edu.scripps.pms.util.seq.Fasta;

public class MergeTest
{
    public static void main(String args[]) throws Exception
    {
	List<Fasta> proteinDbList = FastaReader.getFastaList(new FileInputStream(args[0]));

	//accession and sequence hashtable
	Hashtable<String, String> ht = new Hashtable<String, String>();
	for(Iterator<Fasta> itr=proteinDbList.iterator(); itr.hasNext(); )
	{
	    Fasta fasta = itr.next();
	    ht.put(fasta.getAccession(), fasta.getSequence());
	}

/*
        List<MergeProteinModel> pList = CenSusReportReader.read(args[1], null);
        
        Vector vec = new Vector();
        for(Iterator<MergeProteinModel> itr=pList.iterator(); itr.hasNext(); )
        {
            MergeProteinModel each = itr.next();
            
            List<ChroProtein> chList = each.getProteins();
            
            List<MergeProteinModel.Peptide> pepList = each.getPeptides();

            for(Iterator<ChroProtein> itr1=chList.iterator(); itr1.hasNext(); )
            {
                String proName = itr1.next().getLocus();
                
                String seq = ht.get(proName);
                
                
                for(Iterator<MergeProteinModel.Peptide> itr2=pepList.iterator(); itr2.hasNext(); )
                {
                    MergeProteinModel.Peptide eachPep = itr2.next();
                    String pepSeq = eachPep.getSequence();
                    pepSeq = pepSeq.substring(2, pepSeq.length()-2);
                    
                    int index = seq.indexOf(pepSeq);
                    if(index<0)
		    {
                        System.out.println(proName + "\t" + pepSeq);
		    }
                    else
                        System.out.print("."); //ln("===>>" + proName);
                }
                
            }
            
            
            
            
        }
*/	
	/*
	   List<MergeProteinModel> allList = Collections.EMPTY_LIST; //new Vector<MergeProteinModel>();
	   if(arr.length<=0)
	   return allList;

	   List expList = new Vector();

	   int totalProteinSize = 0;
	   for(int i=0;i<arr.length;i++)
	   {
	   List<MergeProteinModel> tempList = read(arr[i].toString(), fModel);
	   totalProteinSize += tempList.size();
	   expList.add( tempList ); 
	   }

	   allList = (List<MergeProteinModel>)expList.get(0);

	//allList = read(arr[0].toString());
	double eachSeg = (double)100/totalProteinSize;                
	double percent =0;


	for(int i=1;i<arr.length;i++) //each experiment
	{
	//List<MergeProteinModel> list = read(arr[i].toString());
	List<MergeProteinModel> list = (List<MergeProteinModel>)expList.get(i);

	for(Iterator<MergeProteinModel> itr1=list.iterator(); itr1.hasNext(); )  //each protein
	{
	MergeProteinModel protein = itr1.next();
	Set proSet = protein.getProteinNames();
	List<MergeProteinModel.Peptide> pepCompareList = protein.getPeptides();

	//              System.out.println(i + " " + proSet);		
	boolean findSameProtein=false;

	for(Iterator<MergeProteinModel> allItr=allList.iterator(); allItr.hasNext(); )
	{
	MergeProteinModel globalProtein = allItr.next();

	if( globalProtein.isMergeable(pepCompareList, ht) )
	{
	findSameProtein = true;
	globalProtein.addProteinList(protein.getProteins());
	globalProtein.addPeptideList(protein.getPeptides());

	break;
	}
	}                   

	if(!findSameProtein) //add the new protein to the global protein list
	{
	allList.add(protein);
	}

	percent += eachSeg;
	progress.setValue((int)percent);                
	}
	}

	return allList;        
	}    
	 */
    }
}

