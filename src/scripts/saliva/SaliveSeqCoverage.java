package scripts.saliva;

import org.jdom.input.SAXBuilder;
import org.jdom.*;
import java.io.*;
import java.sql.*;
import java.util.*;

import edu.scripps.pms.census.util.dtaselect.Protein;
import edu.scripps.pms.census.util.dtaselect.Peptide;

import edu.scripps.pms.census.model.*;
import edu.scripps.pms.census.util.io.FastaReader;
import edu.scripps.pms.util.seq.Fasta;


public class SaliveSeqCoverage
{
    public static Hashtable<String, String> proteinDbHt;
    public static void main(String args[]) throws Exception
    {
        proteinDbHt = readProteinDatabase();

	System.out.println("done.");

        SAXBuilder sb = new SAXBuilder();
        Document doc = sb.build(new File("/home/rpark/rpark_on_data/project/saliva/sequence_coverage/9.xml"));

	Element root = doc.getRootElement();
	Element proteinList = root.getChild("protein_list");
	List<Element> proteinEleList = proteinList.getChildren("protein");

	int count=0;
	for(Iterator<Element> itr=proteinEleList.iterator(); itr.hasNext(); )
	{
	    Element protein = itr.next();
	    String expIds = protein.getAttributeValue("exp_ids");

	    String proId = protein.getAttributeValue("id");
	    String accession = protein.getAttributeValue("accession");
//	    String[] expArr = expIds.split(",");
//	    String[] spArr = sp.split(",");
	    //List list = Arrays.asList(expArr);

	    getSeqCov(proId, proteinDbHt.get(accession), accession);
/*
	    for(int i=0;i<expArr.length;i++)
	    {
		System.out.println(expArr[i]);

	    }
	    */
//	    count++;

//	    if(count>500)
//		break;
	    
	}

	Set keySet = seqCovHt.keySet();
	for(Iterator<String> itr=keySet.iterator(); itr.hasNext(); )
	{
	    String key = itr.next();
	    boolean[] arr = seqCovHt.get(key);

	    int covCount=0;
	    for(int i=0;i<arr.length;i++)
	    {
		if(arr[i]==true)
		    covCount++;
	    }

	    System.out.println(key + "\t" + covCount + "\t" + arr.length + "\t" + (double)covCount/arr.length);
	}




    }

    public static Hashtable<String, boolean[]> seqCovHt = new Hashtable<String, boolean[]>();

    public static void getSeqCov(String proId, String proteinSeq, String accession) throws Exception {
	//Class.forName("com.mysql.jdbc.Driver");
	Class.forName("org.gjt.mm.mysql.Driver");
	//innerConnection = DriverManager.getConnection("jdbc:mysql://localhost/pmsdb", "pmsuser", "myd474");
	
	Connection con = DriverManager.getConnection("jdbc:mysql://homer.scripps.edu/newpmsdb", "pmsuser", "myd474");
	//con = DriverManager.getConnection("jdbc:mysql://skinner.scripps.edu/pmsdb", "pmsuser", "myd474");

//	String sql = "select pepH.sequence from protein_hit ph, dtaselect_best_peptide b, peptide_hit pepH, protein_hit_identified_protein pidp, protein p where p.id=" + proId + " and ph.id=b.protein_hit_id and b.peptide_hit_id=pepH.id and pidp.protein_id=p.id and pidp.protein_hit_id=ph.id";
	String sql = "select pidp.protein_hit_id from protein_hit_identified_protein pidp, protein p where p.id=" + proId + " and p.id=pidp.protein_id";
	Statement stmt = con.createStatement();
	ResultSet rs = stmt.executeQuery(sql);
	rs.next();
	String proteinHitId= rs.getString("protein_hit_id");

	//System.out.println(proId + "\t" + accession + "\t" + proteinSeq);
	boolean[] seqCov = seqCovHt.get(accession);
	if(null == seqCov) {
	    seqCov = new boolean[proteinSeq.length()];
	} 

	sql = "select pepH.sequence from protein_hit ph, dtaselect_best_peptide b, peptide_hit pepH where ph.id=" + proteinHitId + " and ph.id=b.protein_hit_id and b.peptide_hit_id=pepH.id";
	rs = stmt.executeQuery(sql);

	while(rs.next()) {
	    String seq = rs.getString("sequence");
	    seq = seq.substring(2, seq.length()-2);

	    int start = proteinSeq.indexOf(seq);
	    int end = start + seq.length();

	    for(int i=start;i<end;i++)
	    {
		if(i<0 || i>=proteinSeq.length())
		{	System.out.println("out of bound"); break; }
		seqCov[i] = true;
	    }
	}

	if(null == seqCovHt.get(accession))
	    seqCovHt.put(accession, seqCov);

    }

    public static Hashtable<String, String> readProteinDatabase() throws Exception {
	String fileName = "/home/rpark/rpark_on_data/project/saliva/database/forward_only.fasta";

	Fasta fasta;
	String defLine;
	Hashtable<String, String> proteinDbHt = new Hashtable<String, String>();
	System.out.print("reading protein db...");
	for (Iterator itr = FastaReader.getFastas(new FileInputStream(fileName)); itr.hasNext(); ) {
	    fasta = (Fasta) itr.next();
	    defLine = fasta.getDefline();
//	    System.out.println(fasta.getAccession() + "\t" + fasta.getSequence().length());
	    
	    proteinDbHt.put(fasta.getAccession(), fasta.getSequence());
	}

	return proteinDbHt;
    }
}

