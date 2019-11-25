
/*
* Copyright (c) 2008 Integrated Proteomics Applications.  All rights reserved.  
*/

package scripts.mrm;

/**
 *
 * @author Sung Kyu, Robin, Park
 * @email robinparky@yahoo.com
 * Created on Jan 29, 2010 
 * $Revision:$
 * $Date:$
 */

import edu.scripps.pms.census.util.io.SpectrumReader;
import edu.scripps.pms.util.spectrum.*;
import java.util.*;

import org.jdom.*;
import org.jdom.output.*;
import org.jdom.input.*;

import java.io.File;

public class MRMAnalyzer {

    public MRMAnalyzer() throws Exception{

        //this.parseConfig(configFile);
    }

    public List<PeptideModel> parseConfig(String configFile) throws Exception {
        File file = new File(configFile);
        SAXBuilder builder = new SAXBuilder();

        Document doc = builder.build( file );
        Element paramsEle = doc.getRootElement().getChild("params");

        List<PeptideModel> pepList = new ArrayList<PeptideModel>();
        
        for(Iterator<Element> itr=paramsEle.getChildren("peptide").iterator(); itr.hasNext(); )
        {
            Element pepEle = itr.next();

            //System.out.println(pgroupEle);

            PeptideModel peptide = new PeptideModel(pepEle.getAttributeValue("name"));

            
            for(Iterator<Element> pepitr=pepEle.getChildren("precursor").iterator(); pepitr.hasNext(); )
            {
                Element preEle = pepitr.next();

                double pmass = Double.parseDouble(preEle.getAttributeValue("mass"));

                Precursor p = new Precursor(pmass);

                for(Iterator<Element> ditr=preEle.getChildren("daughter").iterator(); ditr.hasNext(); )
                {
                    Element dauEle = ditr.next();
                    double dmass = Double.parseDouble(dauEle.getAttributeValue("mass"));
                    p.addDaughter(dmass);
                }

                peptide.addPrecursor(p);
            }             
            
            pepList.add(peptide);

        }

        return pepList;
    }


    public static void main(String[] args) throws Exception {

        double precursor1 = 622.220;
        double daughter11 = 830.500;
        double daughter12 = 716.400;


        double precursor2 = 629.720;
        double daughter21 = 495.600;
        double daughter22 = 841.900;



        Precursor lp = new Precursor(precursor1);
        lp.addDaughter(daughter11);
        lp.addDaughter(daughter12);

        Precursor hp = new Precursor(precursor2);
        hp.addDaughter(daughter21);
        hp.addDaughter(daughter22);

        List<Precursor> plist = new ArrayList<Precursor>();
        plist.add(lp);
        plist.add(hp);

        MRMAnalyzer a = new MRMAnalyzer();
        List<PeptideModel> pepList = a.parseConfig("/home/rpark/rpark_on_data/project/jcopping/mrm_012910/mrm_config.xml");
        a.run(pepList);
    }

    public void run(List<PeptideModel> pepList) throws Exception {

        //SpectrumReader sr = new SpectrumReader(args[0], args[1]);
        SpectrumReader sr = new SpectrumReader("/home/rpark/rpark_on_data/project/jcopping/mrm_012910/CF301hrrep1-01.ms2", "ms2");
	Hline h = new Hline(sr.getHlines());

        Iterator<PeakList> it = sr.getSpectra();


        Hashtable<Double, Precursor> pht = new Hashtable<Double, Precursor>();

        for(Iterator<PeptideModel> itr=pepList.iterator(); itr.hasNext(); ) {
            PeptideModel pmodel = itr.next();

            for(Iterator<Precursor> pitr=pmodel.getPlist().iterator(); pitr.hasNext(); ) {
                Precursor precursor = pitr.next();

                pht.put(precursor.getMass(), precursor);
            }

        }

//        int counter = 0;
//        int numPeaks = 0;
        //boolean sortByIntensity = true;
        while (it.hasNext()) {
            PeakList list = it.next();
            
	    Peak p;
	    StringBuffer sb = new StringBuffer();

            Precursor precursor = pht.get(list.getPrecursorMass());

            if(null == precursor)
                continue;

//            System.out.println("===>>" + precursor);
            
            //precursor.addIntensity(p.getM2z(), p.getIntensity());

            for(Iterator<Peak> itr=list.getPeaks(); itr.hasNext(); )
            {
                p = itr.next();
                        
                //precursor.addIntensity(p.getM2z(), p.getIntensity());

                        //double mass = p.getM2z();

                        //if(mass == daughter) {
  //                      System.out.print(p.getM2z());
    //                    System.out.print("\t");
      //                  System.out.print(p.getIntensity());
        //                System.out.print("\t");
                        //}

              
//                System.out.print("\n");

            } 

        }


            


/*
            for(Iterator itr = list.getZlines(); itr.hasNext(); )
            {
                Zline  zline = (Zline)itr.next();

                System.out.println(zline.getM2z() + " " + zline.getChargeState());

            }*/

        
System.out.println(pepList);
        PeptideModel peptide = pepList.get(0);
       // peptide.print(0);
        //peptide.calculateRatio();
        peptide.print();

    }
}
