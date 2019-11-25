package edu.scripps.pms.census.labelFree;

import edu.scripps.pms.census.ChroGenerator;
import edu.scripps.pms.census.conf.Configuration;
import org.jdom.Document;
import org.jdom.Element;
import org.jdom.input.SAXBuilder;
//import sun.security.krb5.Config;

public class LabelfreeOnSingleExperiment {
    public static void main(String[] args) throws Exception {
        System.out.println("aaaaaaaa");


        SAXBuilder sb = new SAXBuilder();

        Document doc = sb.build("/data/2/rpark/ip2_data/hyukim/nomogram/labelfree_quant/labelfree_18910/small_test/test_config.xml");
        final Element root = doc.getRootElement();


        Configuration conf = Configuration.getInstance();
        conf.setConfigFilePath("/data/2/rpark/ip2_data/hyukim/nomogram/labelfree_quant/labelfree_18910/small_test/test_config.xml");

        ChroGenerator chro = new ChroGenerator(
                null, null, null, null, root, null
        );

        chro.createLabelFreeBasedOnDirectId(false, false);
    }
}
