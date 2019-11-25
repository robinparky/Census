
import java.io.*;
import java.util.*;


public class VennDiagram 
{
    public static void main(String[] args) throws Exception
    {
	if(args.length<3)
	{
	    printUsage();
	    return;
	}
	
	if("-i".equals(args[2]))
	{
	    intersection(args[0], args[1]);
	}
	else if("-e".equals(args[2]))
	{
	    exclusive(args[0], args[1]);
	}
	else
	    printUsage();

    }

    public static void printUsage()
    {
	System.out.println("Usage: java VennDiagram file1 file2 option");
	System.out.println("-i : intersection");
	System.out.println("-e : exclusive for file2");
    }

    public static void intersection(String file1, String file2) throws Exception
    {
	Set set = new HashSet();

	BufferedReader br = new BufferedReader(new FileReader(file1));

	String eachLine="";

	while( null != (eachLine=br.readLine()) )
	{
	    set.add(eachLine.trim());
	}

	br = new BufferedReader(new FileReader(file2));

	while( null != (eachLine=br.readLine()) )
	{
	    if( set.contains(eachLine.trim()) )
	    {
		System.out.println(eachLine);
	    }
	}

    }

    public static void exclusive(String file1, String file2) throws Exception
    {
	Set set = new HashSet();

	BufferedReader br = new BufferedReader(new FileReader(file1));

	String eachLine="";

	while( null != (eachLine=br.readLine()) )
	{
	    set.add(eachLine.trim());
	}

	br = new BufferedReader(new FileReader(file2));

	while( null != (eachLine=br.readLine()) )
	{
	    if( !set.contains(eachLine.trim()) )
	    {
		System.out.println(eachLine);
	    }
	}
    }
}


