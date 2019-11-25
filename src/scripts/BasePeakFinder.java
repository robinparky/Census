
/*
 * CenSusReportReader.java
 *
 * Created on February 23, 2006, 3:03 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

import java.io.*;
import java.util.*;

/**
 *
 * @author rpark
 * @version $Id: BasePeakFinder.java,v 1.1 2008/10/30 21:15:15 rpark Exp $ 
 */
public class BasePeakFinder {
    
    private int totalPeptideCount=0;
    
    public static void main(String args[]) throws Exception
    {
	if(args.length<=0)
	{
	    System.out.println("Usage: java BasePeakFinder path option");
	    System.out.println("option 1 : display basepeak for each spectrum (no display in default)");
	    System.exit(0);
	}

	File f = new File(args[0]);
	int option=0;

	if(args.length>=2)
	    option = Integer.parseInt(args[1]);

	String[] arr = f.list();

	for(String each : arr)
	{
	    if(!each.endsWith("ms1"))
		continue;
		
	    BufferedReader br = new BufferedReader(new FileReader(args[0] + "/" + each));

	    String eachLine;

	    double basepeak=0;
	    double sumbasepeak=0;

	    System.out.println(each);

	    while( (eachLine = br.readLine()) != null )
	    {
		if(eachLine.startsWith("S")) {
		    sumbasepeak += basepeak;

		    if(option==1) {
			System.out.println(basepeak);
			System.out.print(eachLine + "\t");
		    }
		    basepeak=0;
		}

		String[] tArr = eachLine.split("\t");
		if(tArr.length>2)
		    continue;

		tArr = eachLine.split(" ");

		double tmpD = Double.parseDouble(tArr[1]);
		if(tmpD>basepeak)
		    basepeak = tmpD;
	    }

	    if(option==1)
		System.out.println(basepeak);

	    System.out.println("Total intensity : " + sumbasepeak);
	}

/*
        for(Iterator<MergeProteinModel> itr=list.iterator(); itr.hasNext(); )
        {
            MergeProteinModel protein = itr.next();
            
            List<ChroProtein> pList = protein.getProteins();

	    double ratioSum = 0;
	    int pepSize = 0;
	    
            for(Iterator<MergeProteinModel.Peptide> pepItr=protein.getPeptides().iterator(); pepItr.hasNext(); )
            {
                MergeProteinModel.Peptide pep = pepItr.next();
		ratioSum += pep.getAverageRatio();
		pepSize++;
            }
                        
            for(Iterator<String> pItr=protein.getProteinNames().iterator(); pItr.hasNext(); )
	}
*/	
    }
}
