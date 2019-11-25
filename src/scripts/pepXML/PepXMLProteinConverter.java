package scripts;

import java.util.*;
import java.io.*;

public class PepXMLProteinConverter
{
    public static void main(String[] args) throws IOException
    {
        BufferedReader br = null;
        String eachLine;
        
	br = new BufferedReader(new FileReader("exp_5.txt"));            

	int totalSpecCount=0;

	while( (eachLine = br.readLine())!=null )
	{
	    String[] arr = eachLine.split(";");

	    System.out.println(arr.length-2);
	    for(int i=1;i<arr.length-1;i++)
	    {
		totalSpecCount += arr[i].split(" ").length/2;

		
//		System.out.println(totalSpecCount + "\t" + (arr[i].split(" ").length/2));
	    }
	    
	}

    }
}
