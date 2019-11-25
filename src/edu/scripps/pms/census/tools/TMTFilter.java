/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.scripps.pms.census.tools;

import edu.scripps.pms.census.io.CenSusReportReader;
import edu.scripps.pms.census.model.ChroProtein;
import edu.scripps.pms.census.model.MergeProteinModel;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

/**
 *
 * @author rpark
 */
public class TMTFilter {
    
    public static void main(String[] args) throws Exception {
        
        String filename = "/home/rpark/pms/Census/build/classes/tmt_dmcclat/census-out.txt";
        if(args.length>0)
            filename = args[0];
        
	CenSusReportReader reader = new CenSusReportReader();
	List<MergeProteinModel> plist = reader.read(filename, null, 0);

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
    
}
