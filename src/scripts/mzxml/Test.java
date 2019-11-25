import java.io.*;
import java.util.*;

import org.jdom.*;
import org.jdom.output.*;
import org.jdom.input.*;


//separate FTMS and ITMS spectra in full scans
public class Test
{
    public static void main(String[] args) throws Exception
    {
	//addIndex("input.txt");
	//addIndex("/home/rpark/rpark_on_data/project/census/census_paper/data/input/1_2_4/120205_yeast_1to1_ITMS/120205_yeast_1to1_LTQ_FTMSComp-06_itms.mzXML");
	//addIndex("/home/rpark/rpark_on_data/project/census/census_paper/data/input/1_2_4/120205_yeast_1to1_ITMS/120205_yeast_1to1_LTQ_FTMSComp-06_itms.mzXML");
	addIndex(args[0]);

	System.out.println("done");
    }

    public static void addIndex(String fileName) throws Exception
    {
	File file = new File(fileName);
	InputStream fisIn = new FileInputStream(fileName);

	System.out.println("filename===>>" + fileName);

	BufferedReader br = new BufferedReader(new FileReader(fileName));

	String temp;

	int size = (int)file.length();

	byte[] byteBuffer = new byte[size];

	fisIn.read(byteBuffer, 0, size);
	StringBuffer sb = new StringBuffer();

	long pos=0;
	boolean findScan = true;
	for(int i=0;i<byteBuffer.length;i++)
	{
	    char ch = (char)byteBuffer[i];

	    if(
	    /*
		(char)byteBuffer[i] == '\"' &&
		(char)byteBuffer[i+1] == '7' &&
		(char)byteBuffer[i+2] == '9' &&
		(char)byteBuffer[i+3] == '5' &&
		(char)byteBuffer[i+4] == '5' &&
		(char)byteBuffer[i+5] == '\"')
	    */
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


	    )	//17936315
	    {		
	    //if(ch == '1')
		pos = i;
		break;
	    }

//	    System.out.print(ch);
	}
	System.out.println(pos + " " + size);


	RandomAccessFile rfile = new RandomAccessFile( new File(fileName), "r");
	rfile.seek(pos);
	byte[] bytes = new byte[50];
	rfile.readFully(bytes);
	for(int i=0;i<bytes.length;i++)
	    System.out.print( (char)bytes[i] );

	
    }
}
