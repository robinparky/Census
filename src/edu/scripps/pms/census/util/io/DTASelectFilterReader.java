package edu.scripps.pms.census.util.io;

import java.io.BufferedReader;
import java.util.*;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.FileReader;
import java.io.RandomAccessFile;

import edu.scripps.pms.census.util.dtaselect.Protein;
import edu.scripps.pms.census.util.dtaselect.Peptide;

import edu.scripps.pms.census.model.*;
import java.io.File;

/**
 * @author  Robin Park
 * @version $Id: DTASelectFilterReader.java,v 1.1 2014/09/09 19:29:52 rpark Exp $
 */
public class DTASelectFilterReader extends BaseIdentificationReader
{
    private final long READ_FROM_THE_END = 500; //position from the end

    private String dbFileName;
    private String dbFilePathAndName;
    //    private InputStreamReader reader;
    private String fileName;
    private BufferedReader br;
    private String lastLine = "";
    private String criteria;
    private int unfilteredProteinNum;
    private int redundantProteinNum;
    private int nonRedundantProteinNum;
    private int unfilteredPeptideNum;
    private int redundantPeptideNum;
    private int nonRedundantPeptideNum;
    private double proteinFP;
    private double peptideFP;
    private double spectrumFP;

    private double confidence;
    private double CONFIDENCE_THRESHOLD=0.05;
    private boolean version1;  //DTASelect version
    private boolean version2;  //DTASelect version
    private Hashtable<String, Integer> proteinIndexHt = new Hashtable<String, Integer>();
    private String label;

    //read total peptide number
    public int getTotalPeptideNumber() throws IOException
    {
	int totalPeptideCount=0;

	DTASelectFilterReader dtaReader = new DTASelectFilterReader(fileName);

	for (Iterator<Protein> itr1 = dtaReader.getProteins(); itr1.hasNext(); ) {
	    Protein protein = itr1.next();

	    totalPeptideCount += protein.getPeptideSize();
	}

	return totalPeptideCount;

    }

    public static void main(String args[]) throws IOException
    {
	DTASelectFilterReader reader = 
                new DTASelectFilterReader("/home/rpark/titus/1311HumanVsCompil/compil/DTASelect-filter.txt0.01");
	//DTASelectFilterReader reader = new DTASelectFilterReader(args[0]);

	Iterator<Protein> pitr = reader.getProteins();

        //ProteinGroup pGroup = new ProteinGroup();
	List l = new ArrayList();
	Set s = new HashSet();
	int decoynum=0;
	int forwardrednum=0;
	int proteinGroup=0;
	boolean isUnitProt=false;
	for (Iterator<Protein> itr = pitr; itr.hasNext(); ) {
	    Protein protein = itr.next();

       //     System.out.println(protein.getLocus() + " " + protein.getPeptideList().size());
	    //if( protein.getLocus().startsWith("Rever") || protein.getLocus().startsWith("contam") || protein.getLocus().startsWith("Contam") )
	    if( protein.getLocus().startsWith("Rever")) {
		decoynum++; continue; }



	forwardrednum++;
	    //ArrayList<Protein> aList = new ArrayList<Protein>();
	    HashSet set = new HashSet();

	    s.add(protein.getLocus());
//	    s.setConfidence(Float.parseFloat(protein.getConfidence()));

	    if(protein.getPeptideSize()<=0)
	    {
	        System.out.println(protein.getDescription());
	        if(protein.getAccession().startsWith("sp"))
	            isUnitProt=true;
	//	pGroup.add(pHit);
	    }
	    else
	    {

            if(protein.getAccession().startsWith("sp")) {
                isUnitProt = true;
                proteinGroup++;
            }
		l.add(s);
		s = new HashSet();
	//	proteinGroup++;

		isUnitProt=false;
//		aList.clear();
	    }
	}

        System.out.println(proteinGroup);
        /*
	System.out.println(reader.getRedundantProteinNum());
	System.out.println(reader.getNonRedundantProteinNum());
	System.out.println(reader.getProteinFP());
	System.out.println(l.size());
	System.out.println(decoynum);

        */

	reader.close();
//	DTASelectFilterReader reader = new DTASelectFilterReader( args[0] );


/*
	Iterator<Protein> pitr = reader.getProteins();

	for (Iterator<Protein> itr = pitr; itr.hasNext(); )
	{
	    Protein protein = itr.next();
	    System.out.println(protein.getLocus());

	}
*/
    }

    public DTASelectFilterReader(String fileName) throws IOException
    {
        this.label = label;
        this.fileName = fileName;
//        this.reader = new FileReader(fileName);
        init();
    }

    /*
    public Iterator <Protein> getProteins(final InputStream is) throws IOException {
        return getProteins(new InputStreamReader(is));
    }
    */
    private void readSummary() throws IOException
    {
        RandomAccessFile file = null;

	try {
		file = new RandomAccessFile(fileName, "r");

		file.seek(file.length()-500);

		String eachLine;
		eachLine=file.readLine();

		while( (eachLine=file.readLine()) != null && !eachLine.startsWith("Unfiltered") );

		String[] arr = eachLine.split("\t");
		this.unfilteredProteinNum = Integer.parseInt(arr[1]);
		this.unfilteredPeptideNum = Integer.parseInt(arr[2]);

		eachLine=file.readLine(); //Redundant line of DTASelect-filter.txt
		arr = eachLine.split("\t");
		this.redundantProteinNum = Integer.parseInt(arr[1]);
		this.redundantPeptideNum = Integer.parseInt(arr[2]);

		eachLine=file.readLine(); //NonRedundant line of DTASelect-filter.txt
		if(null != eachLine) {
			arr = eachLine.split("\t");
			if(arr.length>3) {
				this.nonRedundantProteinNum = Integer.parseInt(arr[1]);
				this.nonRedundantPeptideNum = Integer.parseInt(arr[2]);
			}
		}

		while( (eachLine=file.readLine()) != null && !eachLine.startsWith("Forward FP") && !eachLine.startsWith("Forward FDR"));

		if(null != eachLine) {
			arr = eachLine.split("\t");
			if(arr.length>3) {
				this.proteinFP = Double.parseDouble(arr[1]);

				if(arr.length>2) {
					this.peptideFP = Double.parseDouble(arr[2]);
					this.spectrumFP = Double.parseDouble(arr[3]);
				}
			}
		}
	} catch (Exception e) {
		System.out.println("Error:" + e);	
	} finally {

		if(null != file) file.close(); 
	}
    }


    private void init() throws IOException
    {
        br = new BufferedReader(new FileReader(fileName));

        readHeader();

        while (!(lastLine = br.readLine()).startsWith("Locus"));
        Protein.setFeatureIndices(lastLine);
        
        //Move line to parse protein
        while (!(lastLine = br.readLine()).startsWith("Unique"));
        //add setFeatureIndices here

        Peptide.setFeatureIndices(lastLine);
        lastLine = br.readLine();

    }

    private void readHeader() throws IOException
    {
        lastLine = br.readLine();

        if(lastLine.startsWith("DTASelect v2.0"))
            version2 = true;
        else if(lastLine.startsWith("DTASelect v1")) {
            version2 = false;
            version1 = true;

        }


        //whie (!(lastLine = br.readLine()).endsWith("fasta")); .startsWith("Unique"));
        for(int i=0; i<2; i++)
            lastLine = br.readLine();

        //Remove directory name
        this.dbFileName = lastLine.substring(lastLine.lastIndexOf("/")+1);
        this.dbFilePathAndName = lastLine;

        if(version2)
        {

            //Move line to parameters.  We assume parameters start after carriage return
            //while( !(lastLine = br.readLine()).startsWith("SEQUEST") && !lastLine.startsWith("ProLuCID") && !lastLine.startsWith("TurboS") && !lastLine.startsWith("?") && !lastLine.startsWith("Mascot") && !lastLine.startsWith("Mascot") );


	    //Move line to parameters.  We assume parameters start after carriage return
	    while ((lastLine = br.readLine()) != null) {
		    if(lastLine.startsWith("ProLuCID") || lastLine.startsWith("SEQUEST") || lastLine.startsWith("BlindPTM") || lastLine.startsWith("?") || lastLine.startsWith("MASCOT") || lastLine.startsWith("Blazmass")) {
			    break;
		    }
	    }

            lastLine = br.readLine();

            String tempLine = lastLine;
            int tempIndex = tempLine.indexOf("fp");
            if(tempIndex<0)
            {
                confidence = CONFIDENCE_THRESHOLD;
            }
            else
            {
                tempLine = tempLine.substring(tempIndex+3);

                tempIndex = tempLine.indexOf(" ");
                if(tempIndex>0)
                    tempLine = tempLine.substring(0, tempIndex);
                confidence = Double.parseDouble(tempLine);
            }
            
            StringBuffer sb = new StringBuffer();
            sb.append(lastLine=br.readLine());  //Read this line, which can be either empty or not
            while (!(lastLine = br.readLine()).equals(""))
            {
                sb.append(lastLine);
                sb.append("\n");
            }

            criteria = sb.toString();
        }

        
	String[] arr = lastLine.split("\t");
	for(int i=0;i<arr.length;i++) 
		proteinIndexHt.put(arr[i], i);

    }


    public Iterator <Protein> getProteins() throws IOException {

        return new Iterator<Protein>() {
            private Protein protein;
            private Peptide peptide;

            public boolean hasNext() {
                return lastLine != null && !lastLine.startsWith("\tProteins\t");
            }

            public Protein next() {

                try {
                    protein = getProtein(lastLine.split("\t"));
                }
                catch (IOException e) {
                    e.printStackTrace();
                }

                return protein;
            }

            public void remove() {
                throw new UnsupportedOperationException("Not supported");
            }

            private Protein getProtein(String[] strArr) throws IOException {

//                String[] strArr = lastLine.split("\t");
//		String[] peptideLine;
                //protein = new Protein(strArr);
                protein = new Protein();
                try {
                    protein.setElement(strArr);
                } catch (ArrayIndexOutOfBoundsException e ) {
                    System.out.println(e.getMessage());
                }
                /*
                The format of the DTASelect-filter.txt file peptide line
                 is the following:

                 - a star if the peptide is unique to that protein (optional)
                 - then an integer greater than 1 (optional)
                 - then one of the characters 'M', 'Y', or 'N' (optional)
                 - then tab (mandatory)
                 - then the rest of the fields...

                 Anything that comes until the first tab is optional. You can
                 have a line that starts with "\t", or "*\t", or "2\t", or
                 "*2\t", or "*2M\t", or "Y\t", or "2N\t", or "*Y\t", etc...

		Peptide line starts like (*)(int)(M||Y||N)\t
                 */

                /*
                 *  Some proteins does not have peptide lines, because those proteins
                 *  are assumed to have identical peptides as following protein has.
                 *
                 **/
                
                while ( ((lastLine = br.readLine()) != null && !lastLine.equals(""))
                    && !lastLine.startsWith("\tProteins\t"))
                {
                    strArr = lastLine.split("\t");
                 //   System.out.println("=="  + lastLine);
                 //  System.out.println("2=="  + strArr[2] + "\t" + strArr[2].indexOf(".")); 
                    // If Spectrum Count position does not have a decimal point,
                    // it is not a peptide line
                    if(strArr[2].indexOf(".")<=0 || strArr.length<=5)
                            break;

                    peptide = new Peptide(strArr, version2);
                    protein.addPeptide(peptide);
                }

                return protein;
            }
        };
    }
    public static String getFilterParams(String fname) {

	BufferedReader br = null;
	try {
	    java.io.File f = new java.io.File(fname);
	    if(!f.exists())
		return "NA";

	    br = new BufferedReader(new FileReader(f));

	    String lastLine = null;

	    while ( null != (lastLine = br.readLine())) {
		if(lastLine.contains("SQT format"))
		    break;
	    }

	    lastLine  = br.readLine();

	    br.close();

	    return lastLine;


	} catch (Exception e) {
	    e.printStackTrace();
	    System.out.println("Error: " + e);
	    try {   if(null != br) br.close(); }
	    catch(IOException ie) { }
	    return "NA";

	}
    }


    public void close() throws IOException
    {
        br.close();
    }

    public void setDbFileName(String dbFileName)
    {
        this.dbFileName = dbFileName;
    }

    public void setUnfilteredProteinNum(int unfilteredProteinNum)
    {
        this.unfilteredProteinNum = unfilteredProteinNum;
    }

    public void setRedundantProteinNum(int redundantProteinNum)
    {
        this.redundantProteinNum = redundantProteinNum;
    }

    public void setNonRedundantProteinNum(int nonRedundantProteinNum)
    {
        this.nonRedundantProteinNum = nonRedundantProteinNum;
    }

    public void setUnfilteredPeptideNum(int unfilteredPeptideNum)
    {
        this.unfilteredPeptideNum = unfilteredPeptideNum;
    }

    public void setRedundantPeptideNum(int redundantPeptideNum)
    {
        this.redundantPeptideNum = redundantPeptideNum;
    }

    public void setNonRedundantPeptideNum(int nonRedundantPeptideNum)
    {
        this.nonRedundantPeptideNum = nonRedundantPeptideNum;
    }

    public String getDbFileName()
    {
        return dbFileName;
    }
    public String getDbFilePathAndName()
    {
        return dbFilePathAndName;
    }

    public String getCriteria()
    {
        return criteria;
    }

    public int getUnfilteredProteinNum()
    {
        return unfilteredProteinNum;
    }

    public int getRedundantProteinNum()
    {
        return redundantProteinNum;
    }

    public int getNonRedundantProteinNum()
    {
        return nonRedundantProteinNum;
    }

    public int getUnfilteredPeptideNum()
    {
        return unfilteredPeptideNum;
    }

    public int getRedundantPeptideNum()
    {
        return redundantPeptideNum;
    }

    public int getNonRedundantPeptideNum()
    {
        return nonRedundantPeptideNum;
    }

    public boolean isVersion2() {
        return version2;
    }

    public void setProteinFP(double proteinFP)
    {
	this.proteinFP = proteinFP;
    }

    public double getProteinFP() {
	return proteinFP;
    }

    public void setPeptideFP(double peptideFP)
    {
	this.peptideFP = peptideFP;
    }

    public double getPeptideFP() {
	return peptideFP;
    }
    public void setSpectrumFP(double spectrumFP)
    {
	this.spectrumFP = spectrumFP;
    }

    public double getSpectrumFP() {
	return spectrumFP;
    }
    public String getFileName() {
        return fileName;
    }

    public void setFileName(String fileName) {
        this.fileName = fileName;
    }
    public ArrayList<ChroProtein> getChroProteinList() throws IOException {
        
        ArrayList list = new ArrayList();
        
        for(Iterator<Protein> itr=getProteins(); itr.hasNext(); )
        {
            Protein protein = itr.next();
            
            list.add( ChroProtein.convertProtein(protein) );
            
        }
        
        return list;
    }
    
}
