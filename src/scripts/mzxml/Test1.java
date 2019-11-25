import java.io.*;
import java.util.*;

import org.jdom.*;
import org.jdom.output.*;
import org.jdom.input.*;

//separate FTMS and ITMS spectra in full scans
public class Test1
{
    public static void main(String[] args) throws Exception
    {

	    BufferedReader br = new BufferedReader(new FileReader(args[0]));
	    String fileName = args[0];

	    System.out.println("fileName==>>" + fileName);

	    String eachLine;
	    long pos=0;

	    String temp2Name = fileName + ".new";

	    FileOutputStream out = new FileOutputStream( temp2Name );
	    PrintStream p = new PrintStream( out );

	    while( (eachLine=br.readLine()) != null )
	    {
		p.println(eachLine);
	    }

	    p.close();
	    out.close();


	System.out.println("done");
    }

}
