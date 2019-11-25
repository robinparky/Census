package scripts;

import java.util.*;
import java.io.*;
import edu.scripps.pms.util.stats.StatCalc;

class N15EnrichmentCalc {

	public static void main(String[] args) throws Exception {

		BufferedReader br = null;
		//br = new BufferedReader(new FileReader(args[0]));
                br = new BufferedReader(new FileReader("/home/rpark/rpark_on_data/project/N15/K_all_N15/census-out.txt"));
		String eachLine;
		List<Double> list = new ArrayList<Double>();

                
                while( (eachLine = br.readLine()) != null && !eachLine.startsWith("H\tSLINE\tUNIQUE") );
                
                System.out.println(eachLine);
                String[] arr = eachLine.split("\t");
                
                int lightIntIndex=-1;
                int heavyIntIndex=-1;
                int profileIndex=-1;
                int enrichIndex=-1;
                int bestEnrichCorrIndex=-1;
                int corrPlusIndex=-1;
                int corrMinusIndex=-1;                
                int correctIndex=-1;
                
                for(int i=0;i<arr.length;i++) {
                    if(arr[i].equals("SAM_INT"))
                        lightIntIndex = i + correctIndex;
                    else if(arr[i].equals("REF_INT"))
                        heavyIntIndex = i + correctIndex;
                    else if(arr[i].equals("PROFILE_SCORE"))
                        profileIndex =  i + correctIndex;
                    else if(arr[i].equals("BEST_ENRICH_CORR"))
                        bestEnrichCorrIndex =  i + correctIndex;
                    else if(arr[i].equals("ENRICHMENT"))
                        enrichIndex =  i + correctIndex;
                    else if(arr[i].equals("CORR_ONE_PLUS"))
                        corrPlusIndex =  i + correctIndex;
                    else if(arr[i].equals("CORR_ONE_MINUS"))
                        corrMinusIndex =  i + correctIndex;
                    
                }
                
		//String qType = "rOverR";
		while( (eachLine = br.readLine()) != null ){
			if(eachLine.startsWith("H")) continue;
                        if(eachLine.startsWith("P")) continue;
			arr = eachLine.split("\t");

                        /*
			if(eachLine.startsWith("&S")) {
			//	System.out.println("&" +eachLine);
			//	System.out.println(arr[8] + " " + arr[9] + "\t" + arr[10]);
				Double n14 = Double.parseDouble(arr[8]);
				Double n15 = Double.parseDouble(arr[9]);
				double enrich = n15/(n14+n15);
				list.add(enrich);
			//	System.out.println( n15/(n14+n15) );
			} else {*/
			//	System.out.println(eachLine);
				Double n14 = Double.parseDouble(arr[7]);
				Double n15 = Double.parseDouble(arr[8]);
				double enrich = n15/(n14+n15);
				double isoEnrich = Double.parseDouble(arr[arr.length-1]);
				double diff = Math.abs(enrich-isoEnrich);
				if(diff<0.1)
					list.add(enrich);
				
				//System.out.println(arr[7] + " " + arr[8] + "\t" + arr[arr.length-1]);
			//	System.out.print(enrich+ "\t" + isoEnrich + "\t" + diff + "\t");
				//if(diff<0.1) System.out.println(enrich + "\t" + isoEnrich);
			//	else System.out.println("");

			//}
		}
		
		double[] darr = getDoublArr(list);
		int[] harr = StatCalc.getHistogram(darr, 0.05, 1);
		for(int each:harr) 
			System.out.println(each);
	}


	public static double[] getDoublArr (List<Double> list)
	{
		double[] ret = new double[list.size()];
		Iterator<Double> iterator = list.iterator();
		for (int i = 0; i < ret.length; i++)
		{
			ret[i] = iterator.next().doubleValue();
		}
		return ret;
	}
}
