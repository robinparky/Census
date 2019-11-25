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
 * @version $Id: SpecRangeGenerator.java,v 1.7 2014/08/27 18:00:35 rpark Exp $
 */

public class SpecRangeGenerator 
{
    //private RandomAccessFile file=null;
    private String fileName;

    private BufferedReader br;
    private final int CARRIAGE_RETURN = "\n".getBytes().length;

    private final double XCORR_THRESHOLD_1 = 1.8;
    private final double XCORR_THRESHOLD_2 = 2.5;
    private final double XCORR_THRESHOLD_3 = 3.5;
    private final double XCORR_THRESHOLD_4 = 4.5;

    private double confidence;
    private boolean isVersion2;

    private Hashtable<String, SpecRange> ht = new Hashtable<String, SpecRange>();

    private int pepLocation;


    public static void main(String args[]) throws IOException, Exception
    {
        //SpecRangeGenerator generator = new SpecRangeGenerator(args[0]);
        //SpecRangeGenerator generator = new SpecRangeGenerator("/data/1/rpark/temp/DTASelect.txt", true, 0.05);
        //SpecRangeGenerator generator = new SpecRangeGenerator(args[0], false, 0.05);
        boolean b =false;
        
        if(args[1].equals("1"))
            b = true;
        
        SpecRangeGenerator generator = new SpecRangeGenerator(args[0], b, 0.05);

        Hashtable<String, SpecRange> ht1 = generator.getTable();
        for (Enumeration<String> e = ht1.keys() ; e.hasMoreElements() ;) {
            String key = e.nextElement();
            SpecRange range = ht1.get(key);

	    //if(key.startsWith("YFL014W011406_50percent_enriched-04.itms"))
	    //pfa3D7predicted|chr13|chr13.phat_578|FulPARC_PFK1_1+2_08STLESFFTEIK

        } 
        
        //SpecRangeGenerator generator = new SpecRangeGenerator("E:\\relex\\newFinal_data\\data-independent\\DTASelect.txt");        
    }

    public Hashtable<String, SpecRange> getTable()
    {
        return ht;
    }
    
    public static SpecRangeGenerator getSpecRangeGenerator(IdentificationReader idReader) throws IOException, Exception
    {
        String fileName = idReader.getFileName();
        
        if(fileName.endsWith("DTASelect-filter.txt"))
        {
            return new SpecRangeGenerator(idReader.getFileName().replace("DTASelect-filter.txt", "DTASelect.txt"), idReader.isVersion2(), idReader.getConfidence());    
        }
        else
            return new SpecRangeGenerator();        
    }
    
    public SpecRangeGenerator()
    {
    }
    
    public SpecRangeGenerator(String fileName, boolean isVersion2, double confidence) throws IOException, Exception
    {
        this.fileName = fileName;        
        this.setIsVersion2(isVersion2);        
        this.confidence = confidence;
        createIndex();
    }

    public SpecRange getSpecRange(String key)
    {
        return ht.get(key);
    }

    public void createIndex() throws IOException, Exception
    {
        PrintStream p=null;

        try {
            br = new BufferedReader( new FileReader(fileName) );

            String lastLine;
            long pos = 0;

            lastLine = br.readLine();

            if(isVersion2)
                pepLocation = 12;
            else //if (lastLine.startsWith("DTASelect v2.0"))
                pepLocation = 11;
            
            while(null != (lastLine = br.readLine()) && !lastLine.startsWith("L\t") )
            {
                pos += lastLine.getBytes().length;
                pos += CARRIAGE_RETURN;
            }

            String accession="";
            String prevPep="";
            String curPep="";

            int minSpecNum=-1;
            int maxSpecNum=-1;

            int curScan;
            int prevScan=0;

            while( null != lastLine && !lastLine.startsWith("C\t"))
            {
		
                lastLine = lastLine.substring(lastLine.indexOf("\t") + 1);
                lastLine = lastLine.substring(0, lastLine.indexOf("\t"));
                accession = Fasta.getAccession(lastLine);

		prevPep = "";

                while( null != (lastLine = br.readLine()) && !lastLine.startsWith("C\t"))
                {
		    if( lastLine.startsWith("L\t") )
		    {       
	    		if(!"".equals(prevPep) && !accession.startsWith("Rev"))
                        {
                            ht.put(accession + prevPep, new SpecRange(minSpecNum, maxSpecNum));
                        }
			break;
		    }
    
                    pos += lastLine.getBytes().length;
                    pos += CARRIAGE_RETURN;

                    String[] arr = lastLine.split("\t");
                    
                    String msFileName=arr[1];
                    for(int i=0;i<3;i++)
			msFileName = msFileName.substring(0, msFileName.lastIndexOf("."));
                    
                    //int tempIndex = arr[1].indexOf(".");
                    //String msFileName = arr[1].substring(0, tempIndex);
                    //String temp = arr[1].substring( tempIndex+1 );
                                        
                    String temp = arr[1].substring(msFileName.length()+1);                        
                    //for(int i=0;i<arr.length;i++)

		    try {
			curPep = msFileName + arr[pepLocation].substring(2,arr[pepLocation].length()-2);
		    } catch(java.lang.StringIndexOutOfBoundsException sbe)
		    {

			//System.out.println("Peptide sequence is invalid in DTASelect.txt " + arr[pepLocation]);

			continue;
		    }

                    if(isVersion2)
                    {
                        double conValue = Double.parseDouble(arr[5]);

                        if(conValue<confidence)
                            continue;
                    }
                    else //version 1.9
                    {
                        double xcorr = Double.parseDouble(arr[3]);

                        if( (temp.endsWith("1") && xcorr< this.XCORR_THRESHOLD_1) ||
                                (temp.endsWith("2") && xcorr< this.XCORR_THRESHOLD_2) ||
                                (temp.endsWith("3") && xcorr< this.XCORR_THRESHOLD_3) )
                            continue;
                    }

		    temp = temp.substring(0, temp.indexOf("."));
                   
                    curScan = Integer.parseInt(temp); //scan number

		    SpecRange tran = (SpecRange)ht.get(accession + prevPep);
		    if(null != tran)
		    {
			if(minSpecNum>tran.getMin())
			    minSpecNum = tran.getMin();
			if(maxSpecNum<tran.getMax())
			    maxSpecNum = tran.getMax();
		    }

                    if( curPep.equals(prevPep) )
                    {
                        if(minSpecNum>curScan)
                            minSpecNum = curScan;

                        if(maxSpecNum<curScan)
                            maxSpecNum = curScan;

                    }
                    else
                    {		
	    		//if(!"".equals(prevPep))
	    		if(!"".equals(prevPep) && !accession.startsWith("Rev"))
                        {
                            ht.put(accession + prevPep, new SpecRange(minSpecNum, maxSpecNum));
                        }

                        minSpecNum = curScan;
                        maxSpecNum = curScan;

                    }

                    prevPep = curPep;
                }

            }
	    if(!"".equals(prevPep))
            {
                ht.put(accession + prevPep, new SpecRange(minSpecNum, maxSpecNum));
            }
            
        /*
		  for (Enumeration<String> e = ht.keys() ; e.hasMoreElements() ;) {
		String key = e.nextElement();
		SpecRange range = ht.get(key);
         System.out.println(key + "\t" + range.getMin() +"\t" + range.getMax());

     	} */   
	//SpecRange range = ht.get("ORFP:YGR161C-C_R.TPTAAAWDSPESHIGVAK.K");
        //    System.out.println(range.getMin());
        //    System.out.println(range.getMax());
	
        } catch (IOException e) {
            System.out.println("IO Error : " + e);
            e.printStackTrace();
            throw new IOException("Failed to read DTASelect.txt file : " + e.toString());
        } catch (Exception e) {
            System.out.println("Error : " + e);
            e.printStackTrace();

            throw new Exception("Failed to read DTASelect.txt file : " + e.toString());
        }
        finally {
            try {

                br.close();
            } catch (IOException e) {
                System.out.println("Error to close br : " + e);
            }
        }

    }

    public void setIsVersion2(boolean isVersion2) {
        this.isVersion2 = isVersion2;
    }
}
