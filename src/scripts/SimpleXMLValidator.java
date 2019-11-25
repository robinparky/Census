package scripts;

import java.io.*;
import org.jdom.*;
import org.jdom.input.*;


public class SimpleXMLValidator 
{
    public static void main(String args[])
    {

	try {
            SAXBuilder sb = new SAXBuilder();
	    //Document doc = sb.build( new File("031405_RPal_N14N15_1to1_LTQ_04.mzXML") );
	    Document doc = sb.build( new File(args[0]) );
	    Element ele = doc.getRootElement();

	} catch (Exception e) {
	    // instance document is invalid!

	    System.out.println(e.toString());
	}
    }
}
