
/*
 * ChroReader.java
 *
 * Created on May 17, 2005, 12:00 PM
 */

package edu.scripps.pms.census.io;

import java.io.*;
import java.util.*;

import edu.scripps.pms.census.model.ChroProtein;
import edu.scripps.pms.census.model.ChroPeptide;
import edu.scripps.pms.census.model.ChroData;



/**
 *
 * @author  Robin Park
 * @version $Id: ChroReader.java,v 1.2 2007/01/01 23:20:41 rpark Exp $
 */

public class ChroReader {

//    protected String fileName;
    protected File file;
    protected BufferedReader br;
    protected String lastLine = "";
    protected ArrayList list = new ArrayList();
    protected ChroData data=null;    
    protected int dataDependency=0; //0:data dependent, 1:data indenpendent
    
    public ChroReader(String fileName) throws IOException
    {
        this(new File(fileName));
    }

    public ChroReader(File file) throws IOException
    {
        this.file = file;
        init();
    }
    
    protected void init() throws IOException
    {
        br = new BufferedReader(new FileReader(file));
        
        while( (lastLine=br.readLine()).startsWith("H\t") )
        {
            if(lastLine.startsWith("H\tdata_"))
            {
                String[] strArr = lastLine.split("\t");
                
                if( "1".equals(strArr[2]) )
                    this.dataDependency = 1; //data indenpendent
            }
        }
        
        //Move line to parameters.  We assume parameters start after carriage return
        ChroProtein protein=null;
        ChroPeptide peptide=null;         
    }    

    public Iterator<ChroProtein> getProteins() throws IOException {
        return new Iterator<ChroProtein>() {
            private ChroProtein protein;
            private ChroPeptide peptide;

            public boolean hasNext() {
                //System.out.println("last" + lastLine);
                return lastLine != null; // && !lastLine.startsWith("\tProteins\t");
            }

            public ChroProtein next() {

                try {
                    protein = getProtein();                                            
                }
                catch (IOException e) {
                    e.printStackTrace();
                }

                return protein;
            }

            public void remove() {
                throw new UnsupportedOperationException("Not supported");
            }

            private ChroProtein getProtein() throws IOException {

                //String[] strArr = lastLine.split("\t");
                //protein = new Protein(strArr);
                //System.out.println(lastLine);
                lastLine=br.readLine();
                
                protein = new ChroProtein(lastLine);
                ////System.out.println("----" + lastLine);

                
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
                 */
                
                peptide = null;
                
                while( (lastLine=br.readLine())!=null && !lastLine.startsWith("[PRO") )
                {
                    if(lastLine.startsWith("[PEP"))
                    {
                        if(null != peptide)
                            protein.addPeptide(peptide);
                        
                        peptide = new ChroPeptide(lastLine=br.readLine());
                        
                        lastLine=br.readLine();//Read away [CHROMATOGRAMS] line
                        String[] peakRange = (lastLine=br.readLine()).split("\t"); //read away Peak line;
                        if(peakRange.length==3)
                        {
                            peptide.setStartRange(peakRange[1]);
                            peptide.setEndRange(peakRange[2]);
                        }
                        //Read P line here
                        //lastLine=br.readLine();//Read the first data line
                       
                    }
                    else// if(!lastLine.startsWith("[CHR"))
                    {
                        String[] str=lastLine.split("\t");
                        data = new ChroData(Integer.parseInt(str[0]), Long.parseLong(str[1]), Long.parseLong(str[2]) );
                        peptide.addData(data);
                    }
                }
                
                protein.addPeptide(peptide);

                return protein;
            }
        };
    }
        
    public ArrayList getProteinList() throws IOException    
    {
        ArrayList list = new ArrayList();

        ChroProtein protein;
        
        for(Iterator<ChroProtein> itr=getProteins(); itr.hasNext(); )
        {
            protein = itr.next();
            
            if(protein.getPeptideList().size()!=0)
                list.add(protein);
        }
        
        return list;
    }
    
    public int getDataDependency()
    {
        return this.dataDependency;
    }
    
    public static void main(String args[]) throws IOException
    {
        ChroReader cr = new ChroReader( args[0] );
        
        
        ArrayList<ChroProtein> list = cr.getProteinList();
        
        System.out.println(list.get(0).getProteinLine());

        List pepList = list.get(0).getPeptideList();
        
        ChroPeptide peptide;
        System.out.println(pepList.size());
        for(int i=0;i<pepList.size();i++)
        {
            peptide = (ChroPeptide)pepList.get(i);
            System.out.println( peptide.getPeptideLine() );

        }

        //System.out.println(list.size());
        /*
        DTASelectFilterReader reader = new DTASelectFilterReader( args[0] );

        Iterator<Protein> pitr = reader.getProteins();

        Protein p;
	Peptide peptide;

        for(Iterator<Protein> itr = pitr; itr.hasNext(); )
        {
            p = itr.next();


        }

        StringBuffer sb = new StringBuffer();
        sb.append("SCAN\tTIME\tSAMPLE\tREFERENCE");

        System.out.println( sb.toString() );
*/
        //System.out.println(reader.getDbFileName());
        //System.out.println(reader.getCriteria());
        //System.out.println(reader.getNonRedundantProteinNum());
        //System.out.println(reader.getNonRedundantPeptideNum());
    }

}



/*
 *    private final long READ_FROM_THE_END = 500; //position from the end

    private String dbFileName;
//    private InputStreamReader reader;



    public Iterator <Protein> getProteins() throws IOException {
        return new Iterator<Protein>() {
            private Protein protein;
            private Peptide peptide;

            public boolean hasNext() {
                return lastLine != null && !lastLine.startsWith("\tProteins\t");
            }

            public Protein next() {

                try {
                    protein = getProtein();
                }
                catch (IOException e) {
                    e.printStackTrace();
                }

                return protein;
            }

            public void remove() {
                throw new UnsupportedOperationException("Not supported");
            }

        };
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
}
*/    
