
/*
 * CenSusReportReader.java
 *
 * Created on February 23, 2006, 3:03 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package scripts;

import java.io.*;
import java.util.*;
import edu.scripps.pms.census.model.*;
import edu.scripps.pms.census.util.*;
import edu.scripps.pms.census.util.io.FastaReader;
import edu.scripps.pms.util.seq.Fasta;

/**
 *
 * @author rpark
 * @version $Id: CensusReportReader.java,v 1.3 2014/08/27 18:00:35 rpark Exp $
 */
public class CensusReportReader {
    
    private int totalPeptideCount=0;
    
    public static void main(String args[]) throws Exception
    {
        List<MergeProteinModel> list = Collections.EMPTY_LIST;
	    CensusReportReader re = new CensusReportReader();
        list = re.read(args[0]);
        
        int truePepCount=0;
        int falsePepCount=0;
       
        for(Iterator<MergeProteinModel> itr=list.iterator(); itr.hasNext(); )
        {
            MergeProteinModel protein = itr.next();
            
            List<ChroProtein> pList = protein.getProteins();

	    double ratioSum = 0;
	    int pepSize = 0;
	    
            for(Iterator<MergeProteinModel.Peptide> pepItr=protein.getPeptides().iterator(); pepItr.hasNext(); )
            {
                MergeProteinModel.Peptide pep = pepItr.next();
		ratioSum += pep.getRatio();
		pepSize++;
            }
                        
            for(Iterator<String> pItr=protein.getProteinNames().iterator(); pItr.hasNext(); )
            {
                String cPro = pItr.next();
                
                System.out.print("P\t" + cPro + "\t" + CensusHelper.format.format(ratioSum/pepSize) + "\t");
                System.out.print( CensusHelper.format.format(protein.getStdev())  + "\t" + pepSize + "\t" + protein.getTotalSpecCount()  + "\t" + protein.getDescList().get(0));
                System.out.println("");
                
            }
            
            for(Iterator<MergeProteinModel.Peptide> pepItr=protein.getPeptides().iterator(); pepItr.hasNext(); )
            {
                MergeProteinModel.Peptide pep = pepItr.next();
                System.out.println(pep.getPeptideLine());
                
            }
            
        }
        
//	CensusReportReaderModified reader = new CensusReportReaderModified();
	//reader.calculateDeviationError(args[0]);
        
        System.out.println(truePepCount + " " + falsePepCount);

    }

    /** Creates a new instance of CenSusReportReader */
    public CensusReportReader() {
      
    }

    public List<MergeProteinModel> read(String fileName) throws IOException, Exception
    {        
    /*    BufferedReader br = new BufferedReader(new FileReader(fileName));

        String eachLine;
        
        while( !(eachLine = br.readLine()).startsWith("H\tSLINE") )
	{
	    System.out.println(eachLine);
	}
        
	System.out.println(eachLine);

        boolean isOldFormat = false;
        if(eachLine.split("\t").length<=8)
            isOldFormat = true;
        
        while( !(eachLine = br.readLine()).startsWith("P\t") );
        
        List<MergeProteinModel> list = new Vector<MergeProteinModel>();
        String[] arr = eachLine.split("\t");
        MergeProteinModel protein = new MergeProteinModel();
                
        if(arr.length>2)
            protein.addProteinInfo(arr[1], arr[arr.length-2], arr[arr.length-1], true);
        else
            protein.addProteinInfo(arr[1], "", "", true);
        
        list.add(protein);
        
        boolean isNewProtein=false;        
                
        if(isOldFormat)
        {
            arr = eachLine.split("\t");                

            if(!isNewProtein)
            {
                if(arr.length>2)
                    protein.addProteinInfo(arr[1], arr[arr.length-2], arr[arr.length-1], false);
                else
                    protein.addProteinInfo(arr[1], "", "", false);
            }
            else
            {
                protein = new MergeProteinModel();

                if(arr.length>2)
                    protein.addProteinInfo(arr[1], arr[arr.length-2], arr[arr.length-1], true);
                else
                    protein.addProteinInfo(arr[1], "", "", true);


                list.add(protein);
                isNewProtein = false;
            }
            
            while( (eachLine = br.readLine()) != null)
            {       

                if(eachLine.startsWith("P\t"))
                {
                    arr = eachLine.split("\t");                

                    if(!isNewProtein)
                    {
                        if(arr.length>2)
                            protein.addProteinInfo(arr[1], arr[arr.length-2], arr[arr.length-1], false);
                        else
                            protein.addProteinInfo(arr[1], "", "", false);
                    }
                    else
                    {
                        protein = new MergeProteinModel();
                         
                        if(arr.length>2)
                            protein.addProteinInfo(arr[1], arr[arr.length-2], arr[arr.length-1], true);
                        else
                            protein.addProteinInfo(arr[1], "", "", true);


                        list.add(protein);
                        isNewProtein = false;
                    }

                 } else {
                    
                    totalPeptideCount++;
                    
                    isNewProtein = true;                
                    arr = eachLine.split("\t");

                    boolean isUnique = arr[1].equals("U")?true:false;
                    double ratio = Double.parseDouble(arr[3]);
                    double regFactor = Double.parseDouble(arr[4]);
                    
                    protein.addPeptide(isUnique, arr[2], ratio, regFactor, -1.0, -1.0, -1.0, -1.0, eachLine);
                    
                    
                    
                    System.out.print(this.totalPeptideCount + " " );
                }            
            }            
        }
        else
        {
            //System.out.println("pro " + eachLine + " " + isNewProtein);
            arr = eachLine.split("\t");                

            if(!isNewProtein)
            {
                if(arr.length>2)
                    protein.addProteinInfo(arr[1], arr[arr.length-2], arr[arr.length-1], false);
                else
                    protein.addProteinInfo(arr[1], "", "", false);
            }                
            else
            {
                protein = new MergeProteinModel();

                if(arr.length>2)
                    protein.addProteinInfo(arr[1], arr[arr.length-2], arr[arr.length-1], true);
                else
                    protein.addProteinInfo(arr[1], "", "", true);
                

                list.add(protein);
                isNewProtein = false;
            }
                    
            while( (eachLine = br.readLine()) != null)
            {       
                if(eachLine.startsWith("P\t"))
                {
                    //System.out.println("pro " + eachLine + " " + isNewProtein);
                    arr = eachLine.split("\t");                

                    if(!isNewProtein)
                    {
                        if(arr.length>2)
                            protein.addProteinInfo(arr[1], arr[arr.length-2], arr[arr.length-1], false);
                        else
                            protein.addProteinInfo(arr[1], "", "", false);
                        
                    }
                    else
                    {
                        protein = new MergeProteinModel();

                        if(arr.length>2)
                            protein.addProteinInfo(arr[1], arr[arr.length-2], arr[arr.length-1], true);
                        else
                            protein.addProteinInfo(arr[1], "", "", true);


                        list.add(protein);
                        isNewProtein = false;
                    }

                } else {
                    
                    totalPeptideCount++;
                    
                    //System.out.println("pep " + eachLine + " " + isNewProtein);
                    isNewProtein = true;                
                    arr = eachLine.split("\t");
                    boolean isUnique = arr[1].equals("U")?true:false;
                    double ratio = Double.parseDouble(arr[3]);
                    double regFactor = Double.parseDouble(arr[4]);
                    
                    try {
                        //protein.addPeptide(isUnique, arr[2], ratio, regFactor, Double.parseDouble(arr[6]), Double.parseDouble(arr[7]), arr[arr.length-1], eachLine);
			protein.addPeptide(isUnique, arr[2], ratio, regFactor, -1.0, -1.0, -1.0, -1.0, eachLine);
                    }
                    catch (Exception e) {
                        throw new Exception(fileName + " contains malformed data format : " + eachLine);
                    }
                    
                }            
            }

        }
       */ 
        /*
        for(Iterator<MergeProteinModel> itr = list.iterator(); itr.hasNext(); )
        {
            MergeProteinModel each = itr.next();
            System.out.println("=============>>" + each.getProteinNames());
            System.out.println("=============>>" + each.getPeptides());
        }   
           
        */
        
                       /*     
        if(list.size()>3 && fModel.isPValueSelect())
            edu.scripps.pms.stats.GrubbsTest.filterMerge(list, fModel.getPValue());
                               
        if( fModel.isPValueSelect() && pep.isFilterOut() )
            continue;
                     */       
            //                System.out.println("after" + fModel.isPValueSelect() + " " + fModel.getPValue() + " " + pepList.size());
                            
                            
        return null;        
    }
    

    public int getTotalPeptideCount() {
        return totalPeptideCount;
    }

    public void setTotalPeptideCount(int totalPeptideCount) {
        this.totalPeptideCount = totalPeptideCount;
    }
}
