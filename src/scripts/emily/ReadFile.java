import java.io.*;
import java.util.*;

public class ReadFile
{
	public static void main(String args[]) throws Exception
	{
        BufferedReader br = null;
        String eachLine;
        
            br = new BufferedReader(new FileReader("/home/rpark/pms/RelEx/data/emily/census-outNF.txt"));            

	    Hashtable<String,Set> ht = new Hashtable<String, Set>();
           
	    int totalPepCount=0;
	    String accession =null;

            while( (eachLine = br.readLine())!=null )
	    {
		if( eachLine.startsWith("H\t") )
		    System.out.println(eachLine);

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

		    if(ratio<0.5)
		    {
			double d1 = Double.parseDouble(arr[6]);
			double d2 = Double.parseDouble(arr[7]);

			double tempRatio = d1/d2;
			if(tempRatio>=2.0 || tempRatio<0.5)
			{
			    Set set = ht.get(accession);

			    if(set==null) {
				set = new HashSet();
				set.add(eachLine);
				ht.put(accession, set);
			    }
			    else
			    {
				set.add(eachLine);
			    }

	//		    System.out.println(eachLine);
			}

		    }
		    //System.out.println(arr[6]);
		    //System.out.println(arr[7]);

		}
		
		    

	    }
            
            br.close();

	    //System.out.println(ht);


	    //read next DF file
	    
            br = new BufferedReader(new FileReader("/home/rpark/pms/RelEx/data/emily/census-outDF.txt"));            

            while( (eachLine = br.readLine())!=null )
	    {
		if( eachLine.startsWith("H\t") )
		    System.out.println(eachLine);

		if(eachLine.startsWith("P\t"))
		{
//	System.out.println("==accession" + ht.get(accession));	    

		    Set set = ht.get(accession);
		    if(null != set)
		    {
			for(Iterator itr=set.iterator(); itr.hasNext(); )
			{
			    System.out.println("&" + itr.next()); 
			}

		    }
		    
		    String[] arr = eachLine.split("\t");
		    accession = arr[1];
		    System.out.println(eachLine);
		    
		    continue;
		}

		if(eachLine.startsWith("S\t"))
		{
		    String[] arr = eachLine.split("\t");


/*
		    if(ratio<0.5)
		    {
			double d1 = Double.parseDouble(arr[6]);
			double d2 = Double.parseDouble(arr[7]);

			double tempRatio = d1/d2;
			if(tempRatio>=2.0 || tempRatio<0.5)
			{
			    Set set = ht.get(accession);

			    if(set==null) {
				set = new HashSet();
				set.add(eachLine);
				ht.put(accession, set);
			    }
			    else
			    {
				set.add(eachLine);
			    }

			    System.out.println(eachLine);
			}

		    }
		    //System.out.println(arr[6]);
		    //System.out.println(arr[7]);

		    */
		    System.out.println(eachLine);

		}
	    }
            
            br.close();
	}
}
