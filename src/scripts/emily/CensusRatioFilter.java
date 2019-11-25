import java.io.*;
import java.util.*;

public class CensusRatioFilter 
{
	public static void main(String args[]) throws Exception
	{

	if(args.length<4)
	{
	    System.out.println("Usage: CensusRatioFilter lower_bound_ratio upper_bound_ratio R_sqrt_threshold census_out_unfiltered.txt");

	    System.exit(0);
	}

        BufferedReader br = null;
        String eachLine;
        
            //br = new BufferedReader(new FileReader("/home/rpark/pms/RelEx/data/emily/census-outNF.txt"));            


            br = new BufferedReader(new FileReader(args[3]));

	    double lowerBound = Double.parseDouble(args[0]);
	    double upperBound = Double.parseDouble(args[1]);
	    double rsqrt = Double.parseDouble(args[2]);

	    Hashtable<String,Set> ht = new Hashtable<String, Set>();
           
	    int totalPepCount=0;
	    String accession =null;

            while( (eachLine = br.readLine())!=null )
	    {
//		if( eachLine.startsWith("H\t") )
//		    System.out.println(eachLine);

		if(eachLine.startsWith("P\t"))
		{
		    String[] arr = eachLine.split("\t");
		    accession = arr[1];
	//	    System.out.println(arr[1]);
		    
		    continue;
		}

		if(eachLine.startsWith("S\t"))
		{
		    String[] arr = eachLine.split("\t");

		    double ratio = Double.parseDouble(arr[5]);
		    double rValue = Double.parseDouble(arr[4]);

		    
		    //if(ratio<0.5)
		    if(rValue<=0 || ratio<rsqrt)
		    //if(ratio<rsqrt)
		    {
			double d1 = Double.parseDouble(arr[6]);
			double d2 = Double.parseDouble(arr[7]);

			double tempRatio = d1/d2;
//System.out.println(d1 + " " + d2 + " " + tempRatio);
			//if(tempRatio>=2.0 || tempRatio<0.5)
			if(tempRatio>=upperBound || tempRatio<lowerBound)
			{
			    Set set = ht.get(accession);

			    if(set==null) {
				set = new HashSet();
				set.add(eachLine + "\t" + tempRatio);
				ht.put(accession, set);
			    }
			    else
			    {
				set.add(eachLine + "\t" + tempRatio);
			    }

	//		    System.out.println(eachLine);
			}

		    }
		    //System.out.println(arr[6]);
		    //System.out.println(arr[7]);

		}
		
		    

	    }
            
            br.close();


	    //read next DF file
	    
            //br = new BufferedReader(new FileReader("/home/rpark/pms/Census/data/emily/census-outDF.txt"));            
            br = new BufferedReader(new FileReader(args[3]));

            while( (eachLine = br.readLine())!=null )
	    {
		if( eachLine.startsWith("H\t") )
		{
		    System.out.println(eachLine);
		    continue;
		}

		if(eachLine.startsWith("P\t"))
		{
//	System.out.println("==accession" + ht.get(accession));	    
		    String[] arr = eachLine.split("\t");
		    accession = arr[1];
		    Set set = ht.get(accession);

		    if(null != set)
		    {
			System.out.println(eachLine);
			for(Iterator itr=set.iterator(); itr.hasNext(); )
			{
			    System.out.println("&" + itr.next()); 
			}

		    }
		    
		    
		    continue;
		}

/*
		if(eachLine.startsWith("S\t"))
		{
		    String[] arr = eachLine.split("\t");

		}
		*/
	    }
            
            br.close();
	}
}
