
import java.util.*;
import java.io.*;

public class CensusContrast 
{
    public static Hashtable readFile(String fileName, Hashtable<String, String> pepHt) throws Exception
    {
        BufferedReader br = new BufferedReader(new FileReader(fileName));

        String eachLine;
        
        while( !(eachLine = br.readLine()).startsWith("H\tCorrection") );
        
        boolean isNewProtein=false;        

	while( !(eachLine = br.readLine()).startsWith("P\t") );

	String[] arr = eachLine.split("\t");

	Set s = new HashSet();
	s.add(arr[1]);
	Hashtable<String, Set> ht = new Hashtable<String, Set>();

	while( (eachLine = br.readLine()) != null)
	{       
	    if(eachLine.startsWith("P\t"))
	    {
		arr = eachLine.split("\t");                

		if(!isNewProtein)
		{
		    s.add(arr[1]);
		}
		else
		{
		    s = new HashSet();
		    s.add(arr[1]);
		    //if(arr.length>2)
		//	protein.addProteinInfo(arr[1], arr[arr.length-2], arr[arr.length-1], true);
		  //  else
		//	protein.addProteinInfo(arr[1], "", "", true);


		    //list.add(protein);
		    isNewProtein = false;
		}

	    } else {

		isNewProtein = true;                
		arr = eachLine.split("\t");
		boolean isUnique = arr[1].equals("U")?true:false;
		String ratio = arr[3];
		String regFactor = arr[4];
		String charge = arr[arr.length-1];
		String area1 = arr[arr.length-5];
		String area2 = arr[arr.length-4];

		charge = charge.substring(charge.length()-1);

		Set tmpSet = ht.get(arr[2] + "\t" + charge);

		if(null == tmpSet)
		{
		    ht.put(arr[2] + "\t" + charge, s);
		}
		else
		{
		    tmpSet.addAll(s);

		    ht.put(arr[2] + "\t" + charge, tmpSet);
		}
		
		pepHt.put(arr[2] + "\t" + charge, arr[arr.length-1] + "\t" + ratio + "\t" + regFactor + "\t" + area1 + "\t" + area2);
	    }            
	}

	return ht;
    }

    public static void main(String[] args) throws IOException, Exception
    {       
	List<Hashtable<String, String>> pepHtList = new ArrayList<Hashtable<String, String>>();
	List<Hashtable<String, Set>> setList = new ArrayList<Hashtable<String, Set>>();
    

	for(int i=0;i<args.length;i++)
	{
	    Hashtable<String, String> pepHt = new Hashtable<String, String>();

	    Hashtable<String, Set> ht = readFile(args[i], pepHt);
	    setList.add(ht);
	    pepHtList.add(pepHt);

	}

	Set<String> analyzedPepSet = new HashSet<String>();
	
	for(int i=0;i<setList.size();i++)
	{
	    Hashtable<String, Set> ht = setList.get(i);

	    Set s = ht.keySet();
	    Hashtable<String, String> pepHt = pepHtList.get(i);

	    for(Iterator<String> itr=s.iterator(); itr.hasNext(); )
	    {
		String key = itr.next();
		String value1 = pepHt.get(key);
		Set set1 = ht.get(key);

		if(analyzedPepSet.contains(key))
		    continue;

		analyzedPepSet.add(key);

		System.out.print(key + "\t" + set1 + "\t");
		for(int tmpC=0;tmpC<i;tmpC++)
		{
		    System.out.print("X\tX\tX\tX\tX\t");
		}
/*
		for(int j=0;j<setList.size();j++)
		{
		    if(j==i)
			continue;

		    
		}
*/
		System.out.print(value1 + "\t"); // + "\t" + ((set2==null)?"X\tX\tX\t":pepHt2.get(key)) );
		for(int tmpC=i+1;tmpC<setList.size();tmpC++)
		{
		    Hashtable<String, String> tmpPepHt = pepHtList.get(tmpC);
		    String tValue = tmpPepHt.get(key);
		//    tmpPepHt.remove(key);
		    Hashtable<String, Set> tht = setList.get(tmpC);
		    Set tset = tht.get(key);

		    System.out.print(((tset==null)?"X\tX\tX\tX\tX\t":tmpPepHt.get(key)) + "\t");
		}

		System.out.println("");
	    }
	}

/*
	Hashtable<String, Set> ht1 = readFile(args[0], pepHt1);
	Hashtable<String, Set> ht2 = readFile(args[1], pepHt2);

	s = ht2.keySet();
	for(Iterator<String> itr=s.iterator(); itr.hasNext(); )
	{
	    String key = itr.next();
	    String value = pepHt2.get(key);
	    Set set2 = ht2.get(key);
	    //System.out.println(key + "\tX\t" + set2 + "\t" + "X\tX\tX" + "\t" + value);
	    System.out.println(key + "\t" + set2 + "\t" + set2 + "\t" + "X\tX\tX" + "\t" + value);

	}
	*/
    }
}

