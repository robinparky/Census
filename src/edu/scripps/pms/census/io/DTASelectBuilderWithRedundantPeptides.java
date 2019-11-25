package edu.scripps.pms.census.io;
        
import java.io.*;
import java.util.*;
import edu.scripps.pms.util.seq.Fasta;

import edu.scripps.pms.census.model.SpecRange;
import edu.scripps.pms.census.util.io.IdentificationReader;
import edu.scripps.pms.census.util.dtaselect.Peptide;


/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: Yates Lab</p>
 *
 * @author Robin Park
 * @version $Id: DTASelectBuilderWithRedundantPeptides.java,v 1.2 2014/08/27 18:00:35 rpark Exp $
 */

public class DTASelectBuilderWithRedundantPeptides 
{
    public static void main(String args[]) throws IOException, Exception
    {
	String path = args[0] + File.separator;
	System.out.println("parsing DTASelect.txt...");
        DTASelectForLabelfreeParser generator = new DTASelectForLabelfreeParser(path + "DTASelect.txt", 0.05);
	Hashtable<String, Hashtable> pepHt = generator.getUnfilteredPeptides();
	System.out.println("done.");
/*	
	System.out.println(pepHt);
*/
	String lastLine =null;
	BufferedReader br = new BufferedReader( new FileReader(path + "DTASelect-filter.txt") );
	while( !(lastLine=br.readLine()).startsWith("Unique") )
	    System.out.println(lastLine);

	System.out.println(lastLine);
	while( null != (lastLine=br.readLine()))
	{
//	    if( !lastLine.startsWith("D\t") )
	    String[] arr = lastLine.split("\t");
	    System.out.println(lastLine);
	    if(arr.length<11)
		continue;
		

//	    System.out.println(arr[1].substring(arr[1].length()-1));
	    String pepSeq = arr[arr.length-1];
	    pepSeq = pepSeq.substring(2, pepSeq.length()-2);
	    String fileName = arr[1].substring(0, arr[1].indexOf("."));

	    String key = pepSeq + arr[1].substring(arr[1].length()-1);
//	    System.out.println(fileName);
//	    System.out.println(pepHt.get(key));
	    Hashtable<String,String> pepFileHt = pepHt.get(key);

	    for(Iterator<String> pepFileItr = pepFileHt.keySet().iterator(); pepFileItr.hasNext(); ) {

		String pepFileKey = pepFileItr.next();

		if(fileName.equals(pepFileKey))
		    continue;

		String tmpPepLine = pepFileHt.get(pepFileKey);
		Peptide tmpPep = new Peptide();
		tmpPep.setDTASelectTxtPeptideLine(tmpPepLine);

		StringBuffer sb = new StringBuffer();
		sb.append("");
		sb.append("\t");
		sb.append(tmpPep.getFileName());
		sb.append("\t");
		sb.append(tmpPep.getXCorr());
		sb.append("\t");
		sb.append(tmpPep.getDeltCN());
		sb.append("\t");
		sb.append(tmpPep.getConf());
		sb.append("\t");
		sb.append(tmpPep.getMhPlus());
		sb.append("\t");
		sb.append(tmpPep.getCalcMHplus());
		sb.append("\t");
		sb.append(tmpPep.getTotalIntensity());
		sb.append("\t");
		sb.append(tmpPep.getSpRank());
		sb.append("\t");
		sb.append("0");  //prob score
		sb.append("\t");
		sb.append(tmpPep.getIonProportion());
		sb.append("\t");
		sb.append(tmpPep.getRedundancy());
		sb.append("\t");
		sb.append(tmpPep.getSequence());

		System.out.println(sb.toString());

/*
        this.setUnique((arr[arr.length-1]).startsWith("U"));
        this.setFileName(arr[1]);
        this.setXCorr(arr[3]);
        this.setDeltCN(arr[4]);
        this.setConf( Double.parseDouble(arr[5])*100 );
        this.setMhPlus(arr[6]);
        this.setCalcMHplus(arr[7]);
        this.setTotalIntensity(arr[8]);
        this.setSpRank(arr[9]);
        this.setSpScore(arr[10]);
        this.setIonProportion( Double.parseDouble(arr[10])*100 );
        this.setRedundancy("-1"); //no value
        this.setSequence( arr[12] );
	*/
	    }

/*
        for (Iterator<String> itr = tmpHt.keySet().iterator(); itr.hasNext(); ) {
	    String key = itr.next();
	    System.out.println(key + "===" + tmpHt.get(key));
        } 
  */      
//        System.out.println(generator.getUnfilteredPeptides());

        } 
    }

}
