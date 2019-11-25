
/*
 * CenSusReportReader.java
 *
 * Created on February 23, 2006, 3:03 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package edu.scripps.pms.census.io;

import java.io.*;
import java.util.*;
import edu.scripps.pms.census.model.*;
import edu.scripps.pms.census.util.io.FastaReader;
import edu.scripps.pms.util.seq.Fasta;

/**
 *
 * @author rpark
 * @version $Id: CenSusReportReader.java,v 1.10 2014/08/27 18:00:35 rpark Exp $
 */
public class CenSusReportReader {
    
    private int totalPeptideCount=0;
    
    public static void main(String args[]) throws Exception
    {
        //List<MergeProteinModel> list = CenSusReportReader.read(args[0]);
        
        
        //List<MergeProteinModel> list = Collections.EMPTY_LIST;
	

	/*
        List<MergeProteinModel> list = CenSusReportReader.merge(args, null);
       
       //System.out.println(list);
       
        for(Iterator<MergeProteinModel> itr=list.iterator(); itr.hasNext(); )
        {
            MergeProteinModel protein = itr.next();
            
            List<ChroProtein> pList = protein.getProteins();
                        
            for(Iterator<String> pItr=protein.getProteinNames().iterator(); pItr.hasNext(); )
            {
                String cPro = pItr.next();
                
                System.out.println(cPro);
            }
            
            for(Iterator<MergeProteinModel.Peptide> pepItr=protein.getPeptides().iterator(); pepItr.hasNext(); )
            {
                MergeProteinModel.Peptide pep = pepItr.next();
                System.out.println("   " + pep.getSequence());
            }
            
        }
*/	
       
///home/rpark/rpark_on_data/project/quant_compare_tools/RoverCompareAnalysis/./census/524/rerun/census-out.txt 
///home/rpark/rpark_on_data/project/quant_compare_tools/RoverCompareAnalysis/./census/520/rerun/census-out.txt 
	CenSusReportReader reader = new CenSusReportReader();
	List<MergeProteinModel> plist = reader.read(args[0], null, 0);

	for(Iterator<MergeProteinModel> itr = plist.iterator(); itr.hasNext(); ) {
		MergeProteinModel pro = itr.next();
		List<MergeProteinModel.Peptide> pepList = pro.getPeptides();

		int upepCount=0;	
		//System.out.println(pro.getProteinNames());

		Set<String> set = new HashSet();
		for(Iterator<MergeProteinModel.Peptide> pepitr = pepList.iterator(); pepitr.hasNext(); ) {
			MergeProteinModel.Peptide pep = pepitr.next();

			if(pep.isUnique()) 
			{
				set.add(pep.getSequence());
				upepCount++;
			}
//			System.out.println(pep.isUnique() + " " + pep.getSequence());
		}
		
		// PRINT minimum two unique peptides for max quant comparison
		//if(set.size()>1) 
		//	System.out.println(upepCount + "\t" + set.size());
		//if(false) {
		if(set.size()>1) {
		//if(true) {
			ChroProtein cp = pro.getProteins().get(0);
			System.out.println(cp.getProteinLine());

			for(Iterator<MergeProteinModel.Peptide> pepitr = pepList.iterator(); pepitr.hasNext(); ) {
				MergeProteinModel.Peptide pep = pepitr.next();

				if(pep.isUnique()) 
				System.out.println(upepCount);
				//			System.out.println(pep.isUnique() + " " + pep.getSequence());
			}

		}

	}
	//reader.calculateDeviationError(args[0]);

    }

    /*
    public void calculateDeviationError(String fileName) throws IOException, Exception
    {
	
        List<MergeProteinModel> list = read(fileName);

        for(Iterator<MergeProteinModel> itr=list.iterator(); itr.hasNext(); )
        {
            MergeProteinModel protein = itr.next();
	    double proRatio = protein.getAverageRatio();

	    if(protein.getPeptides().size()>3)
            for(Iterator<MergeProteinModel.Peptide> pepItr=protein.getPeptides().iterator(); pepItr.hasNext(); )
            {
                MergeProteinModel.Peptide pep = pepItr.next();
//                System.out.println("   " + pep.getSequence() + " " + pep.getAverageRatio());
		double diff = proRatio - pep.getAverageRatio();
		//if(diff<0) diff = -diff;

                //System.out.println(diff + "\t" +  pep.getSequence() + " " + pep.getAverageRatio() + " " + proRatio + " " + diff);

		System.out.println(diff);
            }
            
        }
    }
*/
    /** Creates a new instance of CenSusReportReader */
    public CenSusReportReader() {
      
    }

    public List<MergeProteinModel> read(String fileName, FilterModel fModel, double correctFactorValue) throws IOException, Exception
    {        
        BufferedReader br = new BufferedReader(new FileReader(fileName));

        String eachLine;
        
        while( !(eachLine = br.readLine()).startsWith("H\tSLINE") );
        
//        boolean isOldFormat = false;
//        if(eachLine.split("\t").length<=8)
  //          isOldFormat = true;

        while( !(eachLine = br.readLine()).startsWith("P\t") );
        
        List<MergeProteinModel> list = new Vector<MergeProteinModel>();
        String[] arr = eachLine.split("\t");
        MergeProteinModel protein = new MergeProteinModel();

	//System.out.println("==" + eachLine);
	//System.out.println("==" + arr[1] + " == " + arr[6] );
            
        if(arr.length>2)
            //protein.addProteinInfo(arr[1], arr[arr.length-2], arr[arr.length-1], true);
            protein.addProteinInfo(arr[1], arr[6], arr[arr.length-1], true, eachLine);
        else
            protein.addProteinInfo(arr[1], "", "", true, eachLine);
        
        list.add(protein);
        
        boolean isNewProtein=false;        
 
	//it does not support old format any more
	       /*
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
                    
                    if(fModel != null && fModel.isRemoveNegative() && ratio<0)
                        continue;
                    if(fModel != null && fModel.isDetSelect() && fModel.getDetValue()>(regFactor*regFactor))
                        continue;
                    if(fModel != null && fModel.isUniquePeptide() && !isUnique)
                        continue;
                   
		    double correctedRatio = Math.exp(Math.log(ratio) + correctFactorValue);

                    protein.addPeptide(isUnique, arr[2], correctedRatio, regFactor, -1, -1, arr[arr.length-1]);
                    
//			System.out.println("===============>>" + arr[2] + " " + arr[0]);
                }            
            }            
        }
	
        else
        {
	*/
	arr = eachLine.split("\t");                

	if(!isNewProtein)
	{
	    if(arr.length>2)
		protein.addProteinInfo(arr[1], arr[6], arr[arr.length-1], false, eachLine);
	    else
		protein.addProteinInfo(arr[1], "", "", false, eachLine);
	}                
	else
	{
	    protein = new MergeProteinModel();

	    if(arr.length>2)
		protein.addProteinInfo(arr[1], arr[6], arr[arr.length-1], true, eachLine);
	    else
		protein.addProteinInfo(arr[1], "", "", true, eachLine);


	    list.add(protein);
	    isNewProtein = false;
	}

	while( (eachLine = br.readLine()) != null)
	{       
	    if(eachLine.startsWith("&")) //skip singleton peptide
		continue;
		
	    if(eachLine.startsWith("P\t"))
	    {
		arr = eachLine.split("\t");                

		if(!isNewProtein)
		{
		    if(arr.length>2)
			protein.addProteinInfo(arr[1], arr[6], arr[arr.length-1], false, eachLine);
		    else
			protein.addProteinInfo(arr[1], "", "", false, eachLine);

		}
		else
		{
		    protein = new MergeProteinModel();

if("NA".equals(arr[6])) continue;

		    if(arr.length>2)
			protein.addProteinInfo(arr[1], arr[6], arr[arr.length-1], true, eachLine);
		    else
			protein.addProteinInfo(arr[1], "", "", true, eachLine);


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

		//filter
		if(fModel != null && fModel.isRemoveNegative() && ratio<0)
		    continue;
		if(fModel != null && fModel.isDetSelect() && fModel.getDetValue()>(regFactor*regFactor))
		    continue;
		if(fModel != null && fModel.isUniquePeptide() && !isUnique)
		    continue;
		try {
		    double correctedRatio = Math.exp(Math.log(ratio) + correctFactorValue);

//		    String tmpValue = arr[10];

//		    if("INF".equals(tmpValue))
//			tmpValue = "
		   
		    protein.addPeptide(isUnique, 
			arr[2], 
			correctedRatio, 
			regFactor, 
			arr[8].equals("")?0:Double.parseDouble(arr[8]), 
			Double.parseDouble(arr[9]), 
			Double.parseDouble(arr[10]), 
			Double.parseDouble(arr[11]), 
			arr[arr.length-1]);
		}
		catch (Exception e) {
		    System.out.println("Error: " + e);
		    e.printStackTrace();
		    throw new Exception(fileName + " contains malformed data format : " + eachLine);
		}

	    }            
	}

        /*
        for(Iterator<MergeProteinModel> itr = list.iterator(); itr.hasNext(); )
        {
            MergeProteinModel each = itr.next();
        }   
           
        */
        
                       /*     
        if(list.size()>3 && fModel.isPValueSelect())
            edu.scripps.pms.stats.GrubbsTest.filterMerge(list, fModel.getPValue());
                               
        if( fModel.isPValueSelect() && pep.isFilterOut() )
            continue;
                     */       
            //                System.out.println("after" + fModel.isPValueSelect() + " " + fModel.getPValue() + " " + pepList.size());
                            
                            
        return list;        
    }
    
    public List<MergeProteinModel> merge(Object[] arr, javax.swing.JProgressBar progress, String dbFilePath, FilterModel fModel, double correctFactorValue) throws IOException, Exception
    {
        //File dbFile = new File(dbFilePath);
        List<Fasta> proteinDbList = FastaReader.getFastaList(new FileInputStream(dbFilePath));
        
        //accession and sequence hashtable
        Hashtable<String, String> ht = new Hashtable<String, String>();
        for(Iterator<Fasta> itr=proteinDbList.iterator(); itr.hasNext(); )
        {
            Fasta fasta = itr.next();
            ht.put(fasta.getAccession(), fasta.getSequence());
        }
        
        List<MergeProteinModel> allList = Collections.EMPTY_LIST; //new Vector<MergeProteinModel>();
        if(arr.length<=0)
	    return allList;
        
        List expList = new Vector();
        
        int totalProteinSize = 0;
        for(int i=0;i<arr.length;i++)
        {
            List<MergeProteinModel> tempList = read(arr[i].toString(), fModel, correctFactorValue);            
            
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
                
		boolean findSameProtein=false;

                for(Iterator<MergeProteinModel> allItr=allList.iterator(); allItr.hasNext(); )
                {
                    MergeProteinModel globalProtein = allItr.next();

                    //if( globalProtein.isMergeable(pepCompareList, ht) )
                    if( pepCompareList.size()>0 && globalProtein.isMergeable(protein, ht) )
                    {
                        findSameProtein = true;
                        globalProtein.addProteinList(protein.getProteins());
                        globalProtein.addPeptideList(protein.getPeptides());

                        break;
                    }
                }                   
                
                if(!findSameProtein && pepCompareList.size()>0) //add the new protein to the global protein list
                {
                    allList.add(protein);
                }
                
                percent += eachSeg;
                progress.setValue((int)percent);                
            }
        }

        return allList;        
    }    

    public int getTotalPeptideCount() {
        return totalPeptideCount;
    }

    public void setTotalPeptideCount(int totalPeptideCount) {
        this.totalPeptideCount = totalPeptideCount;
    }


}
