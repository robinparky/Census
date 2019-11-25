package scripts;

import javax.xml.transform.dom.*;
import org.w3c.dom.*;
import java.io.*;
import java.util.*;
import javax.xml.validation.*;
import javax.xml.parsers.*;
import javax.xml.*;
import org.xml.sax.*;

public class XMLValidator 
{
    public static void main(String args[]) throws Exception
    {
	SchemaFactory factory =
	    SchemaFactory.newInstance(XMLConstants.W3C_XML_SCHEMA_NS_URI);
	//Schema schema = factory.newSchema(new File("http://sashimi.sourceforge.net/schema_revision/mzXML_2.0/mzXML_2.0.xsd"));
	//Schema schema = factory.newSchema(new File("mzXML_idx_2.0.xsd"));
	Schema schema = factory.newSchema(new File("mzXML_idx_2.1.xsd"));
	Validator validator = schema.newValidator();



	    DocumentBuilder parser =
	    DocumentBuilderFactory.newInstance().newDocumentBuilder();
	//Document document = parser.parse(new File("test031405_RPal_N14N15_1to1_LTQ_04.mzXML"));
	Document document = parser.parse(new File("031405_RPal_N14N15_1to1_LTQ_04.mzXML"));
	//Document document = parser.parse(new File("031405_RPal_N14N15_1to1_LTQ_04.mzXML"));


	            

	// validate the DOM tree
	try {
	    validator.validate(new DOMSource(document));
	} catch (SAXException e) {
	    // instance document is invalid!

	    System.out.println(e.toString());
	}
    }
}
