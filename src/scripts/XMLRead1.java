package scripts;

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
 * @version $Id: XMLRead1.java,v 1.2 2007/10/31 18:52:19 rpark Exp $
 */

public class XMLRead1 {

    
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

	int floatBytes = 32 / 8;
	float[][] tmpMassIntensityList =
	    new float[2][tmpArr.length / floatBytes / 2];
	int peakIndex = 0;
	int fieldIndex = 0;
	int i;

	if (floatBytes <= 0)
	    System.err.println("FLOATBYTES <= 0!!!");

	    System.out.println( "-------------------" + tmpArr.length + " " + floatBytes);
	for (i = 0; i < tmpArr.length - floatBytes + 4; i += floatBytes)
	{
	    int intBits = 0;
	    intBits |= (((int) tmpArr[i]) & 0xff);
	    intBits <<= 8;
	    intBits |= (((int) tmpArr[i + 1]) & 0xff);
	    intBits <<= 8;
	    intBits |= (((int) tmpArr[i + 2]) & 0xff);
	    intBits <<= 8;
	    intBits |= (((int) tmpArr[i + 3]) & 0xff);
	    // Must be in IEEE 754 encoding!

	    System.out.println( Float.intBitsToFloat(intBits) );
	    //System.out.print( (char)tmpArr[i] );
	    tmpMassIntensityList[fieldIndex++][peakIndex] =
		Float.intBitsToFloat(intBits);
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


