import java.io.*;
import java.util.*;

//count sequest searched spectra
public class SearchSpectraCount
{
    public static void main(String args[]) throws Exception
    {
	File file = new File(args[0]);
	String[] arr = file.list();

	System.out.println("Check ms2 files");
	checkMS2(arr);
	System.out.println("");
	System.out.println("Check sqt files");
	checkSQT(arr);

    }

    public static void checkMS2(String[] arr) throws Exception {
	
	List<String> list = new ArrayList<String>();

	for(String each : arr) {
	    if(!each.endsWith("ms2")) 
		continue;

	    list.add(each);
	}

	Object[] larr = list.toArray();
	Arrays.sort(larr);

	//for(Iterator<String> itr=list.iterator(); itr.hasNext(); ) {
	for(Object each : larr) {
	    BufferedReader br = new BufferedReader(new FileReader(each.toString()));

	    String eachLine;
	    int count=0;
	    while( null != (eachLine=br.readLine()) ) {
		if(eachLine.startsWith("Z")) {
		    String[] zarr = eachLine.split("\t");

		    double mass = Double.parseDouble(zarr[2]);

		    if(mass>=500)
			count++;
		}
	    }

	    System.out.println(each + "\t" + count);

	}




    }
     
    public static void checkSQT(String[] arr) throws Exception {
	List<String> list = new ArrayList<String>();

	for(String each : arr) {
	    if(!each.endsWith("sqt")) 
		continue;

	    list.add(each);
	}

	Object[] larr = list.toArray();
	Arrays.sort(larr);

	//for(Iterator<String> itr=list.iterator(); itr.hasNext(); ) {
	for(Object each : larr) {
	    BufferedReader br = new BufferedReader(new FileReader(each.toString()));

	    String eachLine;
	    int count=0;
	    while( null != (eachLine=br.readLine()) ) {
		if(eachLine.startsWith("S"))
		    count++;
	    }

	    System.out.println(each + "\t" + count);
	    

	}

    }

/*
    public static void writeIndex(String fileName) throws Exception
    {
	InputStream fisIn = new FileInputStream(fileName);


	int size = (int)file.length();

	byte[] byteBuffer = new byte[size];

	fisIn.read(byteBuffer, 0, size);
	StringBuffer sb = new StringBuffer();

	long pos;

    }
    */

}

