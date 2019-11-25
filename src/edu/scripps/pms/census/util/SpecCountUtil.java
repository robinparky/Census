package edu.scripps.pms.census.util;

import java.util.*;
import java.io.*;

import edu.scripps.pms.census.util.io.DTASelectFilterReader;
import edu.scripps.pms.census.util.dtaselect.Protein;
import edu.scripps.pms.census.util.dtaselect.Peptide;

/**
 *
 * @author rpark
 * @version $Id: SpecCountUtil.java,v 1.3 2014/08/27 18:00:35 rpark Exp $
 */
public class SpecCountUtil {

    String outFileName;
    
    private ArrayList<String> fList = new ArrayList<String>();

    public void addFile(String fileName)
    {
	fList.add(fileName);
    }
    
    public SpecCountUtil(String outFileName)
    {
        this.outFileName = outFileName;
    }
    
    private class ProteinSp implements Comparable {
       public String accession;
       public int[] spArr;
       
       public ProteinSp (int fileSize)
       {
           spArr = new int[fileSize];
       }        
       
       public int getTotalSpecC()
       {
           int total=0;
           
           for(int i:spArr)
               total+=i;
           
           return total;
       }
       
       public int compareTo(Object o) throws ClassCastException {

            ProteinSp p = (ProteinSp)o;

            int temp = p.getTotalSpecC() - this.getTotalSpecC();

            if(temp<0)
                return -1;
            else if(temp>0)
                return 1;
            else
                return 0;
        }
       
    }
    
    public void runSpecCount() throws IOException
    {
        Hashtable<String, ProteinSp> pHt = new Hashtable<String, ProteinSp>();
        
        int fileSize = fList.size();
        
        int index=0;
        
        
        FileOutputStream out = new FileOutputStream(this.outFileName);
        PrintStream ps = new PrintStream( out );

        
        for(Iterator<String> itr=fList.iterator(); itr.hasNext(); )
        {
            String fileName = itr.next();            
            
            DTASelectFilterReader reader = new DTASelectFilterReader( fileName );

            Protein protein;
            Peptide peptide;
            ArrayList<Protein> aList = new ArrayList<Protein>();
            
            for (Iterator<Protein> pitr = reader.getProteins(); pitr.hasNext(); )
            {
                protein = pitr.next();

                if( protein.getLocus().startsWith("Rever") || protein.getLocus().startsWith("contam") || protein.getLocus().startsWith("Contam") )
                    continue;
                
                String desc = protein.getDescription();
                
                
                ProteinSp psp = pHt.get(protein.getLocus());
                
                if(null == psp) {
                    psp = new ProteinSp(fileSize);
                    psp.accession=protein.getLocus();
                    psp.spArr[index] = Integer.parseInt(protein.getSpectrumCount());
                    pHt.put(protein.getLocus(), psp);
                } else {
                    psp.spArr[index] = Integer.parseInt(protein.getSpectrumCount());
                }
            }

            index++;            
        }
        
        Object[] objArr = pHt.values().toArray();
        
        Arrays.sort(objArr);

        ps.print("Locus\t");
        for(int i=1;i<=fileSize;i++)
        {
            ps.print("Sample" + i);
            ps.print("\t");
        }
        ps.print("\n");
        
        for(int i=0;i<objArr.length;i++)
        {
            ProteinSp p = (ProteinSp)objArr[i];
            ps.print(p.accession);
            ps.print("\t");
            
            int[] spArr = p.spArr;
            
            for(int j=0;j<spArr.length;j++)
            {
                ps.print(spArr[j]);
                ps.print("\t");
            }
            
            ps.println("");
        }
        
        ps.close();
        out.close();

    }
}
