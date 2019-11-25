import java.io.*;
import java.util.*;

public class IndexTest
{
    public static void main(String args[]) throws Exception
    {

	if(args.length<3)
	{
	    System.out.println("Usage: java IndexTest mzxml_file_name start_buffer size");
	    System.exit(0);
	}

//	writeIndex(args[0]);
	int start = Integer.parseInt(args[1]);
	int size = Integer.parseInt(args[2]);
    readRandom(args[0], start, start+size);
    //readRandom(args[0], 1262, 16404);
    }

    public static void writeIndex(String fileName) throws Exception
    {
	File file = new File(fileName);
	InputStream fisIn = new FileInputStream(fileName);

	BufferedReader br = new BufferedReader(new FileReader(fileName));

	int size = (int)file.length();

	byte[] byteBuffer = new byte[size];

	fisIn.read(byteBuffer, 0, size);
	StringBuffer sb = new StringBuffer();

	long pos;

	/*
	for(int i=0;i<500;i++)
	{
		pos = i;
	    System.out.print((char)byteBuffer[i] );
	    System.out.print("\t");
	    System.out.println(pos);
	}
*/
	for(int i=0;i<byteBuffer.length;i++)
	{
	    if( (char)byteBuffer[i] == '<' && (char)byteBuffer[i+1] == 's' && (char)byteBuffer[i+2] == 'c' && (char)byteBuffer[i+3] == 'a')
	    //if( (char)byteBuffer[i] == '<' ) // && (char)byteBuffer[i+1] == 'm' && (char)byteBuffer[i+2] == 's' && (char)byteBuffer[i+3] == 'M')
	    {
		pos = i;
	//	int j = i; //skip S char and space
	//	char ch = (char)byteBuffer[j];
	//	int tabCount=0;

		//System.out.print(Integer.parseInt(sb.toString().trim()));
		for(int j=i;j<i+20;j++)		
		    System.out.print((char)byteBuffer[j] );
		System.out.print("\t");

		//sb.delete(0, sb.length());
	    }
	}

	fisIn.close();


    }

    public static void readRandom(String fileName, int start, int end) throws Exception
    {

	RandomAccessFile rfile = new RandomAccessFile( new File(fileName), "r");

	System.out.println(rfile.length());

	rfile.seek(start);

	int size = end-start;

	System.out.println("size==>>" + rfile.length() + " " + size + " " + start + " " + end);

	byte[] bytes = new byte[size];

	rfile.readFully(bytes);

	for(int i=0;i<size;i++)
	    System.out.print( (char)bytes[i] );

    }

}

