import org.jdom.input.SAXBuilder;
import org.jdom.*;
import java.io.*;
import java.util.*;


//:!java ChroSplit /home/rpark/rpark_on_data/project/census/new_scoring_test/./data/census_chro.xml

public class ChroSplit
{
    public static void main(String args[]) throws Exception
    {
	String fileName = args[0];
	BufferedReader br = new BufferedReader( new FileReader(new File(fileName)) );

	String eachLine;
	StringBuffer head = new StringBuffer();
	StringBuffer tail = new StringBuffer("</relex_chro>");

	while( null != (eachLine=br.readLine()) )
	{
	    head.append(eachLine).append("\n");
	    if(eachLine.trim().startsWith("<quantLevel"))
		break;
	}

//	System.out.println(head.toString());

	fileName = fileName.replace(".xml", "");

	FileOutputStream out = null;
	PrintStream p = null;

	int count=0;
	int fileNameCount=1;
	int proteinCount=0;

	String splitFile = fileName + "_split_" + fileNameCount++ + ".xml";
	out = new FileOutputStream( splitFile );
	p = new PrintStream( out );
	p.println(head.toString());

	while( null != (eachLine=br.readLine()) )
	{
	    
	    if(eachLine.trim().startsWith("<protein l"))
	    {
		proteinCount++;
	    }

	    if(proteinCount%1000==0)
	    {
		p.println("</relex_chro>");
		proteinCount++;
		if(null != p)
		    p.close();
		if(null != out)
		    out.close();

		splitFile = fileName + "_split_" + fileNameCount++ + ".xml";
		out = new FileOutputStream( splitFile );
		p = new PrintStream( out );
		p.println(head.toString());
//		return finalFileName;
		System.out.println(splitFile);
		System.out.println(count);
		count++;
	    }

	    p.println(eachLine);

	}

	p.println("</relex_chro>");
	p.close();
	out.close();

    }

}

