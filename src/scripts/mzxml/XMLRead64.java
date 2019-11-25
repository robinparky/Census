package scripts.mzxml;

/*
 * ChroReader.java
 *
 * Created on May 17, 2005, 12:00 PM
 */

import java.io.*;
import java.util.*;

//import org.jdom.*;
//import org.jdom.output.*;
//import org.jdom.input.*;

import java.nio.*;

import org.systemsbiology.jrap.Base64;
/**
 *
 * @author  Robin Park
 * @version $Id: XMLRead64.java,v 1.2 2007/08/07 17:59:42 rpark Exp $
 */

public class XMLRead64 {

    
    public static byte [] floatTobyte(float num)
    {
	ByteBuffer buf = ByteBuffer.allocate(4);
	buf.putFloat(num);
	return buf.array();
    }

    public static void main(String args[]) throws Exception
    {
    /*
        File file = new File(args[0]);
	SAXBuilder builder = new SAXBuilder();

        Document doc = builder.build( file );
        Element rootEle = doc.getRootElement();
	String str = rootEle.getChildText("msRun");

	System.out.println( str.trim() );

	
	System.out.println( Base64.decodeToString(str.trim()) );

	String encStr = Base64.encodeString(str.trim());

	System.out.println(encStr);
	System.out.println( Base64.decodeToString(encStr) );

*/
	//byte[] barr = Base64.decodeToString("Q3sMUEApm5ZDfgu6P8q+4EN/DyhAy7eVQ4Skp0AyYfhDhXLYQF7bsEOGJEs/wKKRQ4aSS0Bl6wZDhxizP7jKp0OINH9AZ8YUQ4iXnEDDoP5DiTSDQRIOnkOKE+xBFVw+Q4rzdEB6PWJDi4VnQJI/tUOMFww/14UlQ41Os0CRz+RDjaGTQQ3DRUOOD+pABossQ46PUUDiv1ZDjzs0QAGU/kOPtr1AKfXgQ5AuG0AYvh5DkLJEP7VBKEORjy9AZPvHQ5KXEkCQD29DkwF6QCUuV0OULUJAaBELQ5SrZEAALS9DlRPEQE9etEOVpDJA1OxPQ5ZK5ECGvG9DlqdSQOw4bkOXI0JAEJUaQ5iXpkEo4tRDmS9CQI0B6UOZiChAFnmxQ5uoKkFM53hDnCnQQMFIXUOcoE5BHznuQ50RZUB6IHtDnaOqQOEZWkOeFDxAdMAkQ57U");
	//byte[] barr = Base64.decode("Q3sMUEApm5ZDfgu6P8q+4EN/DyhAy7eVQ4Skp0AyYfhDhXLYQF7bsEOGJEs/wKKRQ4aSS0Bl6wZDhxizP7jKp0OINH9AZ8YUQ4iXnEDDoP5DiTSDQRIOnkOKE+xBFVw+Q4rzdEB6PWJDi4VnQJI/tUOMFww/14UlQ41Os0CRz+RDjaGTQQ3DRUOOD+pABossQ46PUUDiv1ZDjzs0QAGU/kOPtr1AKfXgQ5AuG0AYvh5DkLJEP7VBKEORjy9AZPvHQ5KXEkCQD29DkwF6QCUuV0OULUJAaBELQ5SrZEAALS9DlRPEQE9etEOVpDJA1OxPQ5ZK5ECGvG9DlqdSQOw4bkOXI0JAEJUaQ5iXpkEo4tRDmS9CQI0B6UOZiChAFnmxQ5uoKkFM53hDnCnQQMFIXUOcoE5BHznuQ50RZUB6IHtDnaOqQOEZWkOeFDxAdMAkQ57U");

	/*
	byte[] arr = Base64.decode("eJw03Hdcz9/3AHCzSGloa8kuQoOspNc5r7cGKSpCyQrRtPfeGQ0qJDvbx95775kQmUVoGFnF73yd8/v883y8ed/3a917");


	String input = "251.0 4812.650 12125 4.0458 4.450121 255.05 921.583958";
	//byte[] arr1 = new byte[4];
	//byte[] arr1 = floatTobyte(402.7995675464f); 
	//byte[] arr2 = floatTobyte(118000f); 

	byte[] arr = new byte[8];
	int	i;
	for (i=0;i<arr1.length;i++) {
                System.out.println(i + "\t" + arr1[i]);
        }
	for(i=0;i<4;i++)
	{
	   arr[i] = arr1[i]; 
	}

	for(i=4;i<8;i++)
	{
	   arr[i] = arr2[i-4]; 
	}

	for(i=0;i<8;i++)
	{
	    System.out.print( (char)arr[i] );
	}

	String encoded = Base64.encodeBytes(arr1);

	System.out.println("====>>" + encoded);
/*
	String encoded = Base64.encodeString(input);
*/
	String peakData = "Q3sMUEApm5ZDfgu6P8q+4EN/DyhAy7eVQ4Skp0AyYfhDhXLYQF7bsEOGJEs/wKKRQ4aSS0Bl6wZDhxizP7jKp0OINH9AZ8YUQ4iXnEDDoP5DiTSDQRIOnkOKE+xBFVw+Q4rzdEB6PWJDi4VnQJI/tUOMFww/14UlQ41Os0CRz+RDjaGTQQ3DRUOOD+pABossQ46PUUDiv1ZDjzs0QAGU/kOPtr1AKfXgQ5AuG0AYvh5DkLJEP7VBKEORjy9AZPvHQ5KXEkCQD29DkwF6QCUuV0OULUJAaBELQ5SrZEAALS9DlRPEQE9etEOVpDJA1OxPQ5ZK5ECGvG9DlqdSQOw4bkOXI0JAEJUaQ5iXpkEo4tRDmS9CQI0B6UOZiChAFnmxQ5uoKkFM53hDnCnQQMFIXUOcoE5BHznuQ50RZUB6IHtDnaOqQOEZWkOeFDxAdMAkQ57U";

	byte[] tmpArr = Base64.decode("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");

	//int floatBytes = 32 / 8;
	int floatBytes = 64 / 8;
	double[][] tmpMassIntensityList =
	    new double[2][tmpArr.length / floatBytes / 2];
	int peakIndex = 0;
	int fieldIndex = 0;
	int i;

	if (floatBytes <= 0)
	    System.err.println("FLOATBYTES <= 0!!!");

	    System.out.println( "-------------------" + tmpArr.length + " " + floatBytes);
	for (i = 0; i < tmpArr.length - floatBytes + 8; i += floatBytes)
	{
	    long intBits = 0;
	    intBits |= (((long) tmpArr[i]) & 0xff);
	    intBits <<= 8;
	    intBits |= (((long) tmpArr[i + 1]) & 0xff);
	    intBits <<= 8;
	    intBits |= (((int) tmpArr[i + 2]) & 0xff);
	    intBits <<= 8;
	    intBits |= (((long) tmpArr[i + 3]) & 0xff);
	    // Must be in IEEE 754 encoding!
	    intBits <<= 8;
	    intBits |= (((long) tmpArr[i + 4]) & 0xff);
	    intBits <<= 8;
	    intBits |= (((long) tmpArr[i + 5]) & 0xff);
	    intBits <<= 8;
	    intBits |= (((long) tmpArr[i + 6]) & 0xff);
	    intBits <<= 8;
	    intBits |= (((long) tmpArr[i + 7]) & 0xff);

	    System.out.println( Double.longBitsToDouble(intBits) );
	    //System.out.print( (char)tmpArr[i] );
	    tmpMassIntensityList[fieldIndex++][peakIndex] =
		Double.longBitsToDouble(intBits);
	    if (fieldIndex == 2)
	    {
		fieldIndex = 0;
		peakIndex++;
	    }
	}
	    System.out.println( "-------------------" );


	    for(i=0;i<tmpMassIntensityList.length;i++)
	    {
		for(int j=0;j<tmpMassIntensityList[i].length;j++)
		    System.out.print(tmpMassIntensityList[i][j] + " " );

		System.out.println(" " );
		System.out.println(" " );
	    }

    }    


}

