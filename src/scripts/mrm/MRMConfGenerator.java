
/*
* Copyright (c) 2008 Integrated Proteomics Applications.  All rights reserved.  
*/

package scripts.mrm;

/**
 *
 * @author Sung Kyu, Robin, Park
 * @email robinparky@yahoo.com
 * Created on Apr 5, 2010 
 * $Revision:$
 * $Date:$
 */

import java.io.*;
import java.util.*;
import org.jdom.Element;
import org.jdom.Comment;
import org.jdom.Document;
import org.jdom.output.*;

public class MRMConfGenerator {

    public MRMConfGenerator() {

    }

    public void parseMethodFile(String filename, String path, String ms2fileName) throws Exception {
        BufferedReader br = new BufferedReader(new FileReader(filename));

        String eachLine;

        while( (eachLine = br.readLine()) != null )
        {
            eachLine = eachLine.trim();
            if(eachLine.startsWith("Parent"))
                break;
        }

        Element rootEle = new Element("config");
        Element labelType = new Element("label_type");
        labelType.setAttribute("labeling", "true");
        rootEle.addContent(labelType);
        Element paramEle = new Element("params");
        Element mrmparamEle = new Element("mrm_params");
        Element spectraFile = new Element("ms2_file");
        spectraFile.addContent(ms2fileName);
        mrmparamEle.addContent(spectraFile);

        String prevPrecursorMass="";
        Element peptideEle = null;
        Element precursorEle = null;

        while( (eachLine = br.readLine()) != null )
        {
            eachLine = eachLine.trim();

            if("".equals(eachLine) || eachLine.startsWith("Syringe"))
                break;

            //System.out.println(eachLine);
            String[] arr = eachLine.split(" +");

            String curPrecursorMass = arr[0];
            if(!curPrecursorMass.equals(prevPrecursorMass)) {
                peptideEle = new Element("peptide");
                peptideEle.setAttribute("name", arr[0] + "_" + arr[1]);
                precursorEle = new Element("precursor");
                precursorEle.setAttribute("mass", arr[0]);
                precursorEle.setAttribute("name", "none");
                precursorEle.setAttribute("desc", "none");
            }


            Element transitionEle = new Element("transition");
            transitionEle.setAttribute("mass", arr[1]);            
            precursorEle.addContent(transitionEle);

            if(!curPrecursorMass.equals(prevPrecursorMass)) {
                peptideEle.addContent(precursorEle);
                mrmparamEle.addContent(peptideEle);
            }

            prevPrecursorMass = curPrecursorMass;
        }



        /*
        //rootEle.addContent(new Comment(expTypeSb.toString()));
        //rootEle.addContent(expType);
        
        //rootEle.addContent(eleCompEle);
*/
        paramEle.addContent(mrmparamEle);
        rootEle.addContent(paramEle);
        
        try
        {
            File outputFile = new File(path + File.separator + "mrm_config.xml");
            Document doc = new Document(rootEle);
            OutputStream os = new FileOutputStream(outputFile); //(filePath + "census_chro.xml");
            XMLOutputter outputter = new XMLOutputter();
            outputter.setFormat(Format.getPrettyFormat());
            outputter.output(doc, os);
            os.close();
        } catch (IOException e)
        {
            System.out.println(e.toString());
        }
    }

    public static void main(String[] args) throws Exception {
        //String filename = args[0];
        MRMConfGenerator gen = new MRMConfGenerator();

//        gen.parseMethodFile(filename);
        //String methodFile = "/home/rpark/rpark_on_data/project/jcopping/mrm_012910/tsq_method.txt";
        //String outputPath = "/home/rpark";
        //String ms2fileName = "/home/rpark/rpark_on_data/project/jcopping/mrm_012910/CF301hrrep1-01.ms2";

        if(args.length<3) {
            System.out.println("Usave java MRMConfGenerator methodFile, outputPath, ms2fileName");
            System.exit(0);
        }

        String methodFile = args[0];
        String outputPath = args[1];
        String ms2fileName = args[2];

        gen.parseMethodFile(methodFile, outputPath, ms2fileName);

        System.out.println(outputPath + File.separator + "mrm_config.xml was successfully generated");


        
    }
}
