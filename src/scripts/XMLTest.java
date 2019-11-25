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
 * @version $Id: XMLTest.java,v 1.2 2007/10/31 18:52:19 rpark Exp $
 */

public class XMLTest {

    
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
	byte[] barr = Base64.decode("eJw03Hdcz9/3AHCzSGloa8kuQoOspNc5r7cGKSpCyQrRtPfeGQ0qJDvbx95775kQmUVoGFnF73yd8/v883y8ed/3a917");


	String input = "251.0 4812.650 12125 4.0458 4.450121 255.05 921.583958";
	//byte[] arr1 = new byte[4];
	byte[] arr1 = floatTobyte(402.7995675464f); 
	byte[] arr2 = floatTobyte(118000f); 

/*	
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
*/

		
	String encoded = Base64.encodeBytes(arr1);

	System.out.println("====>>" + encoded);
/*
	String encoded = Base64.encodeString(input);
*/
/*
	String peakData = "Q3sMUEApm5ZDfgu6P8q+4EN/DyhAy7eVQ4Skp0AyYfhDhXLYQF7bsEOGJEs/wKKRQ4aSS0Bl6wZDhxizP7jKp0OINH9AZ8YUQ4iXnEDDoP5DiTSDQRIOnkOKE+xBFVw+Q4rzdEB6PWJDi4VnQJI/tUOMFww/14UlQ41Os0CRz+RDjaGTQQ3DRUOOD+pABossQ46PUUDiv1ZDjzs0QAGU/kOPtr1AKfXgQ5AuG0AYvh5DkLJEP7VBKEORjy9AZPvHQ5KXEkCQD29DkwF6QCUuV0OULUJAaBELQ5SrZEAALS9DlRPEQE9etEOVpDJA1OxPQ5ZK5ECGvG9DlqdSQOw4bkOXI0JAEJUaQ5iXpkEo4tRDmS9CQI0B6UOZiChAFnmxQ5uoKkFM53hDnCnQQMFIXUOcoE5BHznuQ50RZUB6IHtDnaOqQOEZWkOeFDxAdMAkQ57U";


	byte[] tmpArr = Base64.decode("eJw03Hdcz9/3AHCzSGloa8kuQoOspNc5r7cGKSpCyQrRtPfeGQ0qJDvbx95775kQmUVoGFnF73yd8/v883y8ed/3a917");
	int floatBytes = 32 / 8;
	float[][] tmpMassIntensityList =
	    new float[2][tmpArr.length / floatBytes / 2];
	int peakIndex = 0;
	int fieldIndex = 0;

	if (floatBytes <= 0)
	    System.err.println("FLOATBYTES <= 0!!!");

	    System.out.println( "-------------------" );
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

*/
    
    }    


}


