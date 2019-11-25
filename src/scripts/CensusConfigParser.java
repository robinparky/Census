/*
 * Configuration.java
 *
 * Created on March 18, 2005, 3:51 PM
 */

import java.util.*;

import java.io.*;

import java.text.*;

import org.jdom.*;
import org.jdom.output.*;
import org.jdom.input.*;

import javax.swing.*;

/**
 *
 * @author  Robin Park
 *
 */

//Singleton object
public class CensusConfigParser {
    
   
	public static void main(String[] args)  throws Exception {  

	readXMLParam(args[0]);
}
 
    //read config xml file
    public static void readXMLParam(String fileName) throws Exception { 

        SAXBuilder builder = new SAXBuilder();
       
        Document doc = builder.build(new File(fileName));
        Element rootEle = doc.getRootElement();

	List<Element> sList = rootEle.getChild("element_comp").getChildren("each_sample");

	for(Iterator<Element> itr=sList.iterator(); itr.hasNext(); )
	{
		Element each = itr.next();
		List<Element> rList = each.getChildren("residue");

		for(Iterator<Element> ritr=rList.iterator(); ritr.hasNext(); )
		{
			Element residueEle = ritr.next();
			String rname = residueEle.getAttributeValue("name");
			System.out.print(rname);
			System.out.print("\t");
			System.out.print( residueEle.getChildText("ele_C") );
			System.out.print("\t");
			System.out.print( residueEle.getChildText("ele_H") );
			System.out.print("\t");
			System.out.print( residueEle.getChildText("ele_O") );
			System.out.print("\t");
			System.out.print( residueEle.getChildText("ele_N") );
			System.out.print("\t");
			System.out.print( residueEle.getChildText("ele_S") );
			System.out.print("\t");
			System.out.print( residueEle.getChildText("ele_P") );
			System.out.print("\t");
			System.out.print( residueEle.getChildText("ele_15N") );
			System.out.print("\t");
			System.out.print( residueEle.getChildText("ele_2H") );
			System.out.print("\t");
			System.out.print( residueEle.getChildText("ele_13C") );
			System.out.print("\n");

		}
	}


/*
        Element paramEle = rootEle.getChild("params");        
        Element qtypeEle = rootEle.getChild("quant_type");        
	if(null != qtypeEle)
	else { //for back compatibility we need to check if it is 15N or not



		if(sList.size()>1) {
			Element hEle = sList.get(1);
			
			List<Element> rList = hEle.getChildren("residue");
			boolean n15allvalue=true;
			if(n15allvalue)
				conf.setQuantType("15N");
		}

	}
		

	Element mrmParamEle = paramEle.getChild("mrm_params");
	if(null != mrmParamEle) {
	    this.expType = Configuration.MRM_WITHOUT_ID;
	    this.setSpectrumFormat( Configuration.MS2_FILE_FORMAT );
	}
        
        readMRMParams(mrmParamEle);
//        BufferedReader br = new BufferedReader(new FileReader(filePath + File.separator + fileName));
        Hashtable<String, String> paramTable = new Hashtable<String, String>();
        
        String scanType = paramEle.getChildText("scan_type");
        
        String labeling = rootEle.getChild("label_type").getAttributeValue("labeling");
        
        if("false".equals(labeling))
            this.isLabeling = false;
        else
            this.isLabeling = true;

        String useMassDiff = rootEle.getChild("label_type").getAttributeValue("useMassDiff");
        if("true".equals(useMassDiff))
            this.useMassDiff = true;
        else
            this.useMassDiff = false;

        String expTypeStr = rootEle.getChildText("experiment_type");
        
        if(null != expTypeStr && !"".equals(expTypeStr))
            this.expType = Integer.parseInt(expTypeStr);
       
        if( "MS".equals(scanType) )
        {
            this.quantLevel=1;        
	//    conf.setSpectrumFormat( Configuration.MS_FILE_FORMAT );
        }
        else if( "MS/MS".equals(scanType) )
        {
            if(conf.getSpectrumFormat() != Configuration.MZXML_FILE_FORMAT)
                conf.setSpectrumFormat( Configuration.MS2_FILE_FORMAT );
            
            this.quantLevel=2;
            
            Element msmsParams = paramEle.getChild("msms_params");

	    //Element 
         
            if(null != msmsParams)
            {

                String scanShiftStr = msmsParams.getAttributeValue("scan_shift");

		if(null != scanShiftStr)
			scanShift = Integer.parseInt(scanShiftStr);
		//System.out.println("mass shift===========" + scanShift);
		
                Element msmsFragParams = paramEle.getChild("msms_params").getChild("frag_ext_type");

                if("single".equals(msmsFragParams.getChildText("spectrum")))
                {
                    this.msmsSpectrumNum = this.MSMS_SINGLE_SPECTRUM;
                }
                else if("multiple".equals(msmsFragParams.getChildText("spectrum")) )
                {
                    this.msmsSpectrumNum = this.MSMS_MULTIPLE_SPECTRA;
                }

                if(null != msmsFragParams)
                {
                    String tmpMsmsType = msmsFragParams.getAttributeValue("type");

                    if("auto".equals( tmpMsmsType ))
                        this.setMsmsFragType(this.AUTOMATIC_FRAGMENT_ION);
                    else if("specific".equals(tmpMsmsType))
                    {
                        this.setMsmsFragType(this.SPECIFIC_FRAGMENT_ION);

                        List l = msmsFragParams.getChildren("specific_mass");

                        for(Iterator<Element> massItr=l.iterator(); massItr.hasNext(); )
                        {
                            Element massEach = massItr.next();
                            getMsmsMassArr().add( massEach.getText() );
                        }

                        setMsmsSpecificTolerance(Double.parseDouble( msmsFragParams.getChildText("msms_tolerance")  ));                                                
                    }
                }                
            }
            
        }
        
        Element alignElement = paramEle.getChild("alignment");
        
        if(null != alignElement)
        {
            if( "false".equals(alignElement.getAttributeValue("align")) )
                this.align=false;
            else
                this.align=true;

	    Element peakCountEle = alignElement.getChild("peak_count");

	    if(null != peakCountEle) {
		this.peakCount = peakCountEle.getText();
	    }

	    Element basedOnIdEle = alignElement.getChild("based_on_id");
	    if(null != basedOnIdEle)
	    {
		String basedOnIdStr = basedOnIdEle.getText();
		if("true".equals(basedOnIdStr))
		    conf.setBasedOnId(true);
                
                this.databaseFile = basedOnIdEle.getAttributeValue("database_file");

	    }
        }
        
        //if(null == ht.get(CensusConstants.DATA_TYPE))
        //    throw new Exception("Data type must be defined.");
        
        //int tempInt = Integer.parseInt( (String)ht.get(CensusConstants.DATA_TYPE) );        
        //if(tempInt==0) this.isDataIndependent=false;
        
        String tempIsoWin = paramEle.getChildText("iso_window");
        
        if(null != tempIsoWin)
	{
            this.isolationWindow = Double.parseDouble( tempIsoWin );
	    
	    if(this.expType == CensusConstants.MSMS_DATA_INDEPENDENT)
		this.isolationWindow -= 0.03; //correction factor for data independent
	}
       
	String tmpExtMethod = paramEle.getChildText("extract_method");

	if(null != tmpExtMethod)
	    this.extMethod = Integer.parseInt( tmpExtMethod );
        
        //this.startMassRange = Double.parseDouble( paramEle.getChildText("start_mass") );
        //this.endMassRange = Double.parseDouble( paramEle.getChildText("end_mass") );

	String tmpMassTol = paramEle.getChildText("mass_accuracy");

	if(null != tmpMassTol)
	    this.massTolerance = Double.parseDouble( tmpMassTol );

	String tmpIsodistThreshold = paramEle.getChildText("isodist_threshold");

	if(null != tmpIsodistThreshold)
	    this.isodistThreshold = Double.parseDouble( tmpIsodistThreshold );


        //if(extMethod==1) //m/z            


	String unitValue = null;
	if(null == paramEle.getChild("mass_accuracy"))
		unitValue = "mz";
	else 
		unitValue = paramEle.getChild("mass_accuracy").getAttributeValue("unit");

        if("ppm".equals(unitValue)) //ppm
	{
            this.massTolerance *= 0.001;         
	    this.isHighRes=true;        
	}
        
	Element enrichEle = paramEle.getChild("enrich");
        
        if(null!=enrichEle) 
	{
	    String tempEnrich = paramEle.getChildText("enrich");
	    if(null != tempEnrich && !"".equals(tempEnrich)) {
		this.enrichment = Double.parseDouble( tempEnrich );
	    }

	    String sEnrich = enrichEle.getAttributeValue("start");

	    if(null != sEnrich) {
		this.calculateEnrich = true;
		this.startEnrich = Double.parseDouble(enrichEle.getAttributeValue("start"));
		this.endEnrich = Double.parseDouble(enrichEle.getAttributeValue("end"));
		this.enrichmentMaxDeviation = Double.parseDouble(enrichEle.getAttributeValue("max_deviation"));
	    }
	}
        
        String maxWin = paramEle.getChildText("max_win");
	if(null != maxWin)
	    this.maxWindow = Integer.parseInt( maxWin );

        String isotopePeakStr = paramEle.getChildText("isotope_peak");
	if(null != isotopePeakStr)
	    this.isotopePeak = Integer.parseInt( isotopePeakStr );

	Element prolineEle = paramEle.getChild("proline_use");
	if(null != prolineEle)
	{
	    String proValue = prolineEle.getText();

	    if("true".equals(proValue))
		useProline = true;

	    String prolineC = prolineEle.getAttributeValue("num");
	    if(null != prolineC)
		prolineCount = Integer.parseInt(prolineC);

	    String prolineS = prolineEle.getAttributeValue("shift");
	    if(null != prolineS)
		PROLINE_SHIFT = Double.parseDouble(prolineS);
	}


	String specShift = paramEle.getChildText("max_spectrum_shift");
	if(null != specShift)
	    this.maxSpectrumShift = Integer.parseInt(specShift);

	String marginValue = paramEle.getChildText("margin");
	if(marginValue != null)
	    this.margin = Integer.parseInt(marginValue); //  This is hard coded for now

*/
    }
   
    
}
