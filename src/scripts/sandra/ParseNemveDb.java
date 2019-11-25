package scripts.sandra;

import org.jdom.input.SAXBuilder;
import org.jdom.*;
import java.io.*;
import java.util.*;
import edu.scripps.pms.census.util.io.FastaReader;
import edu.scripps.pms.util.seq.Fasta;

//check census svm file and separate true and false peptide
public class ParseNemveDb
{
    public static String fastaFile = "/home/rpark/rpark_on_data/project/sandra/database/test1.fasta";

    public static void main(String args[]) throws Exception
    {
	//file 1 : /home/rpark/rpark_on_data/deploy/global/census_svm_test.txt
	//file 2 : /home/rpark/rpark_on_data/deploy/global/census_svm_predict.txt
	List<Fasta> proteinDbList = FastaReader.getFastaList(new FileInputStream(fastaFile));

	for(Iterator<Fasta> itr=proteinDbList.iterator(); itr.hasNext(); )
	{
	    Fasta fasta = itr.next();

	    String acc = fasta.getAccession();
	    String proteinSeq = fasta.getSequence();
	    String defline = fasta.getDefline();
	    String newDef = defline.substring(0);

	    System.out.println(acc);
	}
    }
}
