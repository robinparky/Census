package edu.scripps.pms.census.io;
        
import java.io.*;
import java.util.*;
import edu.scripps.pms.util.seq.Fasta;

import edu.scripps.pms.census.model.SpecRange;
import edu.scripps.pms.census.util.io.IdentificationReader;


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
 * @version $Id: DTASelectForLabelfreeParser.java,v 1.2 2014/08/27 18:00:35 rpark Exp $
 */

public class DTASelectForLabelfreeParser 
{
    //private RandomAccessFile file=null;
    private String fileName;

    private BufferedReader br;
    private final int CARRIAGE_RETURN = "\n".getBytes().length;

    private double confidence;

    private Hashtable<String, Hashtable> ht = new Hashtable<String, Hashtable>();

    public static void main(String args[]) throws IOException, Exception
    {
        DTASelectForLabelfreeParser generator = new DTASelectForLabelfreeParser(args[0], 0.05);
	Hashtable<String, Hashtable> tmpHt = generator.getUnfilteredPeptides();

/*
        for (Iterator<String> itr = tmpHt.keySet().iterator(); itr.hasNext(); ) {
	    String key = itr.next();
	    System.out.println(key + "===" + tmpHt.get(key));
        } 
  */      
//        System.out.println(generator.getUnfilteredPeptides());

    }

    public Hashtable<String, Hashtable> getUnfilteredPeptides() {
	return ht;
    }

    public DTASelectForLabelfreeParser(String fileName, double confidence) throws IOException, Exception
    {
        this.fileName = fileName;        
        this.confidence = confidence;

        PrintStream p=null;

	br = new BufferedReader( new FileReader(fileName) );
	int pepLocation = 12;

	String lastLine =null;


	while( null != (lastLine=br.readLine()))
	{
	    if( !lastLine.startsWith("D\t") )
		continue;

	    String[] arr = lastLine.split("\t");
	    String seq = arr[pepLocation];
	    seq = seq.substring(2, seq.length()-2);

	    if(confidence>Double.parseDouble(arr[5]))
		continue;

	    //System.out.println(seq + " " + seq.substring(2, seq.length()-2));

	    String pepFileName = arr[1].substring(0, arr[1].indexOf("."));
	    //System.out.println(pepFileName);
	    String cstate = arr[1].substring(arr[1].length()-1);

	    String seqKey = seq+cstate;

	    Hashtable<String, String> fileHt = ht.get(seqKey);	    
	    if(null == fileHt) {
		fileHt = new Hashtable<String, String>();
		fileHt.put(pepFileName, lastLine);
		ht.put(seqKey, fileHt);
	    } else {
		String pepLine = fileHt.get(pepFileName);

		if(pepLine==null) {
		    pepLine = lastLine;
		    fileHt.put(pepFileName, lastLine);

		} else {
		    double currentConf = Double.parseDouble(arr[5]);
		    double bestConf = Double.parseDouble(pepLine.split("\t")[5]);

		    if(bestConf<currentConf)
		    {
			fileHt.put(pepFileName, lastLine);
			bestConf = currentConf;
		    }

//		System.out.println("==" + arr[5]);
//		System.out.println("==" + currentConf);
//		System.out.println("==>" + lastLine);

		}
	    }
	}

	if(null != br)
	    br.close();
    }
}
