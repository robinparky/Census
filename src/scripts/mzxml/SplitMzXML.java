import java.io.*;
import java.util.*;

import org.jdom.*;
import org.jdom.output.*;
import org.jdom.input.*;


//separate FTMS and ITMS spectra in full scans
public class SplitMzXML
{
    private static int KEEP_ITMS=1; //first scan
    private static int KEEP_FTMS=2; //second scan
    private static Namespace ns = Namespace.getNamespace("http://sashimi.sourceforge.net/schema_revision/mzXML_2.0");
    private static int lastIndex=0;

    //split high and low scan (FTMS and ITMS)
    //args[0] : mzxml file name
    //args[1] : 1 keep first scan(ITMS) and 2 second scan(FTMS)
    public static void main(String[] args) throws Exception
    {

	if(args.length<2)
	{
	    System.out.println("split high and low scan (FTMS and ITMS)");
	    System.out.println("args[0] : mzxml file name");
	    System.out.println("args[1] : 1 keep first scan(ITMS) and 2 second scan(FTMS)");
	    System.exit(0);
	}

//	if(true)
//	    addIndex(args[0]);

//	    System.exit(0);
	
	String newFileName = splitFile(args);
	//String newFileName = "test123_itms_temp.mzXML";
	
	//addIndex(newFileName);
	///addIndex("test123.mzXML");
	//String indexStr = addIndex("120205_yeast_1to1_LTQ_FTMSComp-03.mzXML");
	String indexStr = addIndex(newFileName);
	//"120205_yeast_1to1_LTQ_FTMSComp-03.mzXML");

//	String fileName = "test123.mzXML";

	String temp2Name = insertIndex(newFileName, indexStr);
	String finalFileName = replaceInsertOffset(temp2Name);

	System.out.println(finalFileName + " was generated.");
	System.out.println("done");
    }


    public static String replaceInsertOffset(String fileName) throws Exception
    {
	BufferedReader br = new BufferedReader( new FileReader(fileName) );
	String eachLine;

	String finalFileName = fileName.replace("_temp2.", ".");
	FileOutputStream out = new FileOutputStream( finalFileName );
	PrintStream p = new PrintStream( out );

	while( (eachLine=br.readLine()) != null )
	{
	    if(eachLine.trim().startsWith("<indexOffset>"))
	    {
		p.print("<indexOffset>");
		p.print(lastIndex);
		p.println("</indexOffset>");
	    }
	    else
		p.println(eachLine);

	}
     
	p.close();
	out.close();

	return finalFileName;



    /*
	    if(
		(i+4)<byteBuffer.length &&
		(char)byteBuffer[i] == '/' && 
		(char)byteBuffer[i+1] == 'm' &&
		(char)byteBuffer[i+2] == 's' &&
		(char)byteBuffer[i+3] == 'R' &&
		(char)byteBuffer[i+4] == 'u' &&
		(char)byteBuffer[i+5] == 'n' && 
		(char)byteBuffer[i+6] == '>' 
	    )
	    {
		System.out.println("---------------");
		lastIndex = i;
	    }

	while( (eachLine=br.readLine()) != null )
	{
	    if(eachLine.trim().startsWith("<indexOffset>"))
	    {
		p.print("<indexOffset>");
		p.print(lastIndex);
		p.println("</indexOffset>");
	    }
	    else
		p.println(eachLine);

	}

	*/

    }
    

    public static String insertIndex(String fileName, String indexStr) throws Exception
    {
	BufferedReader br = new BufferedReader(new FileReader(fileName));

	String eachLine;
	long pos=0;

	String temp2Name = fileName.replace("_temp.", "_temp2.");

	FileOutputStream out = new FileOutputStream( temp2Name );
	PrintStream p = new PrintStream( out );

	while( (eachLine=br.readLine()) != null )
	{
	    p.println(eachLine);
	    if(eachLine.trim().startsWith("</msRun>"))
		break;
	}
     
	p.println(indexStr);
	
	while( (eachLine=br.readLine()) != null )
	   p.println(eachLine);

	p.close();
	out.close();

	//replace indexOffset
	//String temp3Name = fileName.replace("_temp.", "_temp3.");
	//out = new FileOutputStream( temp3Name );
	//p = new PrintStream( out );

	File file = new File(temp2Name);

	InputStream fisIn = new FileInputStream(temp2Name);

	br = new BufferedReader(new FileReader(temp2Name));

	int size = (int)file.length();
	byte[] byteBuffer = new byte[size];

	fisIn.read(byteBuffer, 0, size);
	StringBuffer sb = new StringBuffer();

	for(int i=0;i<byteBuffer.length;i++)
	{
	    if(
		(i+4)<byteBuffer.length &&
		(char)byteBuffer[i] == '/' && 
		(char)byteBuffer[i+1] == 'm' &&
		(char)byteBuffer[i+2] == 's' &&
		(char)byteBuffer[i+3] == 'R' &&
		(char)byteBuffer[i+4] == 'u' &&
		(char)byteBuffer[i+5] == 'n' && 
		(char)byteBuffer[i+6] == '>' 
	    )
	    {
		lastIndex = i+7;
	    }
	}

	fisIn.close();

	//p.close();
	//out.close();
/*	
	String finalFileName = fileName.replace("_temp2.", ".");
	br = new BufferedReader(new FileReader(finalFileName));
	FileOutputStream out = new FileOutputStream( finalFileName );
	p = new PrintStream( out );

	while( (eachLine=br.readLine()) != null )
	{
	    p.println(eachLine);
	    if(eachLine.trim().startsWith("</msRun>"))
		break;
	}
     
	p.println(indexStr);
	
	while( (eachLine=br.readLine()) != null )
	    p.println(eachLine);

	p.close();
	out.close();
*/
	return temp2Name;
    }

    public static String addIndex(String fileName) throws Exception
    {

	StringBuffer indexSb = new StringBuffer("<index name=\"scan\">");
	indexSb.append("\n");

	File file = new File(fileName);
	InputStream fisIn = new FileInputStream(fileName);

	BufferedReader br = new BufferedReader(new FileReader(fileName));

	String temp;

	int size = (int)file.length();

	byte[] byteBuffer = new byte[size];

	fisIn.read(byteBuffer, 0, size);
	StringBuffer sb = new StringBuffer();

	long pos;
	boolean findScan = true;
	for(int i=0;i<byteBuffer.length;i++)
	{
	/*
	    if(
		(char)byteBuffer[i] == '<' &&
		    (char)byteBuffer[i+1] == 's' &&
		    (char)byteBuffer[i+2] == 'c' &&
		    (char)byteBuffer[i+3] == 'a' &&
		    (char)byteBuffer[i+4] == 'n' &&
		    (char)byteBuffer[i+6] == 'n' &&
		    (char)byteBuffer[i+7] == 'u' &&
		    (char)byteBuffer[i+8] == 'm' &&
		    (char)byteBuffer[i+9] == '=' &&
		    (char)byteBuffer[i+10] == '\"' &&
		    (char)byteBuffer[i+11] == '7' &&
		    (char)byteBuffer[i+12] == '9' &&
		    (char)byteBuffer[i+13] == '5' &&
		    (char)byteBuffer[i+14] == '5'


	      )   //17936315
*/


	    
	    if( 
		findScan &&
		(char)byteBuffer[i] == '<' && 
		(char)byteBuffer[i+1] == 's' && 
		(char)byteBuffer[i+2] == 'c' &&
		(char)byteBuffer[i+3] == 'a' &&
		(char)byteBuffer[i+4] == 'n' &&
		(char)byteBuffer[i+6] == 'n' && 
		(char)byteBuffer[i+7] == 'u' &&
		(char)byteBuffer[i+8] == 'm' &&
		(char)byteBuffer[i+9] == '=' &&
		(char)byteBuffer[i+10] == '\"'
		)
	    {
		//17936305 o
		//17964149 x
	//	System.out.println("==-->>" + (i-2));

		indexSb.append("<offset id=\"");

		int j = i+11;
	//	System.out.println("");
		while( (char)byteBuffer[j] != '\"' )
		{
		    indexSb.append((char)byteBuffer[j]);
	//	System.out.print((char)byteBuffer[j]);
		    j++;
		}

		indexSb.append("\">");
		indexSb.append(i-2);

	//	System.out.println("==>>" + (i-2));
		indexSb.append("</offset>");
		indexSb.append("\n");

		while( true )
		{
		    if(
			(char)byteBuffer[j] == 'm' &&
			(char)byteBuffer[j+1] == 's' &&
			(char)byteBuffer[j+2] == 'L' &&
			(char)byteBuffer[j+3] == 'e' &&
			(char)byteBuffer[j+4] == 'v' &&
			(char)byteBuffer[j+5] == 'e' &&
			(char)byteBuffer[j+6] == 'l' && 
			(char)byteBuffer[j+7] == '=' 
			)  break;
		    
		   j++; 
		}
		    

		findScan = false;

	    }

	    if(
		!findScan &&
		(char)byteBuffer[i] == 'i' && 
		(char)byteBuffer[i+1] == 'n' &&
		(char)byteBuffer[i+2] == 't' &&
		(char)byteBuffer[i+3] == '\"' &&
		(char)byteBuffer[i+4] == '>' 
	    )
	    {

		findScan = true;

	    }


	}

	indexSb.append("</index>");
	fisIn.close();

	return indexSb.toString();


    }

    public static String splitFile(String[] args) throws Exception
    {

	//args[0] mzXML file
	//e.g. /home/rpark/rpark_on_data/project/census/mzxml/orbit_yeast/120205_yeast_1to1_LTQ_FTMSComp-06_itms.mzXML

	int type = Integer.parseInt(args[1]);

	Document doc = new SAXBuilder().build(new FileReader(args[0]));
	Element rootEle = doc.getRootElement();
	Element indexEle = rootEle.getChild("index", ns);
	indexEle.detach();

	Element msRunEle = rootEle.getChild("msRun", ns);
	
//	List<Element> chList = msRunEle.getChildren("scan"); //rootEle.getChildren();
	List<Element> chList = msRunEle.getChildren("scan", ns); //rootEle.getChildren();

//	System.out.println(chList);

	//for(Iterator<Element> itr=chList.iterator(); itr.hasNext(); )
	int prevScan=-1;
	int currentScan;
	int size = chList.size();

	Element prevEle=null;

	for(int i=size-1;i>=0;i--) 
	{
//	    Element e = itr.next();
	    Element e = chList.get(i); 

	    int scan = Integer.parseInt( e.getAttributeValue("num") );
	    int msLevel = Integer.parseInt( e.getAttributeValue("msLevel") );
	    currentScan = scan;

	    int diff = prevScan - currentScan;
	    
	    //System.out.println(scan + " " + msLevel + " " + currentScan + " " + prevScan + " " + diff + " " + type);

	    if(diff == 1)
	    {
	//	msRunEle.removeContent(e);

		if(type == KEEP_FTMS)
		    e.detach(); 
		else if(type == KEEP_ITMS)
		    prevEle.detach();
	    }

	    prevScan = currentScan;
	    prevEle = e;
	}

//	Document doc = new Document(rootEle);

	String newFileName = args[0].replace(".mzXML", "_" + ((type==KEEP_FTMS)?"ftms_temp.mzXML":"itms_temp.mzXML"));
	OutputStream os = new FileOutputStream(new File(newFileName)); //(filePath + "census_chro.xml");
	XMLOutputter outputter = new XMLOutputter();
	outputter.setFormat(Format.getPrettyFormat());
	outputter.output(doc, os);
	os.close();

	Process p = Runtime.getRuntime().exec("dos2unix " + newFileName);
	p.waitFor();

	return newFileName;

    }
}
