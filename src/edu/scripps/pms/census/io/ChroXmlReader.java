
/*
 * ChroReader.java
 *
 * Created on May 17, 2005, 12:00 PM
 */

package edu.scripps.pms.census.io;

import java.io.*;
import java.util.*;

import edu.scripps.pms.census.model.*;
import edu.scripps.pms.census.*;
import edu.scripps.pms.census.conf.Configuration;
import edu.scripps.pms.census.util.IsoData;

import org.jdom.*;
import org.jdom.output.*;
import org.jdom.input.*;

import java.net.URI;
import java.util.logging.Level;
import java.util.logging.Logger;
/**
 *
 * @author  Robin Park
 * @version $Id: ChroXmlReader.java,v 1.43 2014/08/06 05:29:53 rpark Exp $
 */

public class ChroXmlReader {

    private BufferedReader br;
    private String lastLine = "";
    private ArrayList list = new ArrayList();
    private ChroData data=null;
    private URI uri=null;
    private boolean isDataDependent;

    private SAXBuilder builder = new SAXBuilder();
    private String fileName;
    private File file=null;
    private Element rootEle;
    private Configuration conf;
    private int quantLevel=1; //default quant is from MS
    private boolean labeled=true;
    private String quantType=null;
    private int expType = -1;

    private ArrayList<String> fileList = new ArrayList<String>();
    private ArrayList<String> sampleList = new ArrayList<String>(); //sample list
    private Hashtable<String, String> fileSampleHt = new Hashtable<String, String>();  //filename(key) sample(value)
    private Hashtable<String, Sample> sampleObjList = new Hashtable<String, Sample>(); //sample array
    private ArrayList<Sample> sampleExpList = new ArrayList<Sample>();

    public Sample getSampleObj(String sampleName)
    {
        return this.sampleObjList.get(sampleName);
    }

    public void addSampleObj(Sample sample)
    {
        this.sampleObjList.put(sample.getSampleName(), sample);
    }

    public static class Sample
    {
        private ArrayList intensityArr = new ArrayList();
        private ArrayList expList = new ArrayList();

        private String sampleName;

        private ArrayList<Double> normalizedSpecC = new ArrayList<Double>();

        public void addNormalizedSpecC(double d)
        {
            this.normalizedSpecC.add(d);
        }

        //average of normalized values
        public double getNormalizedSpecC()
        {
            double sum=0;

            for(Iterator<Double> itr=normalizedSpecC.iterator(); itr.hasNext(); )
            {
                double each = itr.next().doubleValue();
                sum += each;
            }

            return sum/normalizedSpecC.size();
        }

        public void addExperiment(String expName)
        {
            this.expList.add(expName);
        }

        public ArrayList getExpList()
        {
            return this.expList;
        }

        public Sample(String sampleName)
        {
            this.sampleName = sampleName;
        }


        public void addIntensity(long intensity)
        {
            this.intensityArr.add(intensity);
        }

        public long getSumIntensity()
        {
            long total=0;

            for(Iterator<Long> itr=intensityArr.iterator(); itr.hasNext(); )
            {
                total += itr.next();
            }

            return total;
        }

        public long getAverage()
        {
            if(this.intensityArr.size() <=0 )
                return 0;

            return this.getSumIntensity() / this.intensityArr.size();
        }

        public String getSampleName() {
            return sampleName;
        }

        public void setSampleName(String sampleName) {
            this.sampleName = sampleName;
        }


        public ArrayList getIntensityArr() {
            return intensityArr;
        }

        public void setIntensityArr(ArrayList intensityArr) {
            this.intensityArr = intensityArr;
        }
    }

    public ChroXmlReader(String fileName) throws IOException, JDOMException, Exception
    {
        this(new File(fileName), true);
    }

    public ChroXmlReader(String fileName, boolean doInit) throws IOException, JDOMException, Exception
    {
        this(new File(fileName), doInit);
    }

    public ChroXmlReader(File file) throws IOException, JDOMException, Exception
    {
	this(file, true);
    }

    public ChroXmlReader(File file, boolean doInit) throws IOException, JDOMException, Exception
    {
        this.file = file;

	if(doInit)
	    init();
    }

    public ChroXmlReader(java.net.URI uri) throws IOException, JDOMException, Exception
    {
	this.uri = uri;
	this.file = null;

	init();
    }

    private void init() throws IOException, JDOMException, Exception
    {
        Document doc = null;

	if( null != this.uri) {
	    InputStream in = uri.toURL().openStream();
	    doc = builder.build( in );
	} else {
	    doc = builder.build( file );
	}

        rootEle = doc.getRootElement();
        String expTypeStr = rootEle.getChildText("experiment_type");

        if( !"relex_chro".equals(rootEle.getName()) )
            throw new JDOMException();

        if(null == conf)
            conf = Configuration.getInstance();

        if( "1".equals(rootEle.getChildText("data_dependency")) )
        {
            this.isDataDependent= true; //data denpendent
            conf.setIsDataIndependent(false);
        }
        else if( "0".equals(rootEle.getChildText("data_dependency")) )
        {
            this.isDataDependent= false; //data indenpendent
            conf.setIsDataIndependent(true);
        }


        if(null != expTypeStr)
            setExpType(Integer.parseInt(expTypeStr));

        String labelType = rootEle.getChildText("label_type");

        if(null != labelType)
            if("false".endsWith(labelType))
                this.labeled = false;
            else
                this.labeled = true;
        else if("census_chro_temp.xml".equals(file.getName()))
            this.labeled = false;

	conf.setLabeling(this.labeled);
        String quantLevelString = rootEle.getChildText("quantLevel");
        quantType = rootEle.getChildText("quantType");

        if(null != quantLevelString )
            this.quantLevel = Integer.parseInt( quantLevelString );

	conf.setQuantLevel(this.quantLevel);
	Element each = null;
	Element eachFile = null;
	Sample s = null;

        for(Iterator<Element> itr=rootEle.getChildren("sample").iterator(); itr.hasNext(); )
        {
            each = itr.next();

            this.sampleList.add( each.getChildText("name") );
            s = new Sample(each.getChildText("name"));

            for(Iterator<Element> fitr=each.getChild("ms_files").getChildren("file").iterator(); fitr.hasNext(); )
            {
                eachFile = fitr.next();
                fileList.add(eachFile.getText());
                s.addExperiment(eachFile.getText());

                fileSampleHt.put(eachFile.getText(), each.getChildText("name"));

            }

            this.sampleExpList.add(s);
        }

        this.getProteinList(); //initialize all parameters, especially for stand-alone report generation
        //Move line to parameters.  We assume parameters start after carriage return
    }

    public ArrayList<ChroProtein> getMrmCrvProteinList() throws IOException, JDOMException
    {
        Document doc = builder.build( file );
        rootEle = doc.getRootElement();

        ArrayList<ChroProtein> list = new ArrayList<ChroProtein>();

        ChroProtein protein;
        ChroPeptide peptide;
        ChroData data;

        List proList = rootEle.getChildren("protein");
        for(Iterator<Element> itrPro=proList.iterator(); itrPro.hasNext(); )
        {
            Element eachProtein = itrPro.next();

            protein = new ChroProtein();
            protein.setLocus( eachProtein.getAttributeValue("locus") );
            protein.setSeqCount( eachProtein.getAttributeValue("seq_ct") );
            protein.setSpectrumCount( eachProtein.getAttributeValue("spec_ct") );
            protein.setSeqCoverage( eachProtein.getAttributeValue("seq_cov") );
            protein.setLength( eachProtein.getAttributeValue("length") );
            protein.setMolWt( eachProtein.getAttributeValue("molwt") );
            protein.setPI( eachProtein.getAttributeValue("pi") );
            protein.setValidation( eachProtein.getAttributeValue("val") );
            protein.setDescription( eachProtein.getAttributeValue("desc") );

            List pepList = eachProtein.getChildren("peptide");

            for(Iterator<Element> itrPep=pepList.iterator(); itrPep.hasNext(); )
            {
                Element eachPeptide = itrPep.next();

                peptide = new ChroPeptide();
                peptide.setUnique( !"".equals(eachPeptide.getAttributeValue("unique")) );
                peptide.setFileName( eachPeptide.getAttributeValue("file") );
                peptide.setScanNum( Integer.parseInt(eachPeptide.getAttributeValue("scan")) );
                peptide.setSequence( eachPeptide.getAttributeValue("seq") );
                peptide.setXCorr( eachPeptide.getAttributeValue("xcorr") );
                peptide.setChargeState( eachPeptide.getAttributeValue("charge") );
                peptide.setDeltCN( eachPeptide.getAttributeValue("deltaCN") );
                peptide.setMhPlus( eachPeptide.getAttributeValue("MHPlus") );
                peptide.setDeltMass( eachPeptide.getAttributeValue("deltaMass") );
                peptide.setSpRank(eachPeptide.getAttributeValue("spRank") );
                peptide.setSpScore(eachPeptide.getAttributeValue("spScore") );

                peptide.setDtaStartRange( Integer.parseInt(eachPeptide.getAttributeValue("start_scan")) ); //For dtaselect range
                peptide.setDtaEndRange( Integer.parseInt(eachPeptide.getAttributeValue("end_scan")) ); //for dtaselect range

		addScores(peptide, eachPeptide);

                String[] dataArr = eachPeptide.getChildText("chro").split(";");
                String[] peakArr = dataArr[0].split(" ");

                data = new ChroData(); //i, Long.parseLong(intArr[1]), Long.parseLong(intArr[2]));

                data.setBsIntensity( dataArr[1].split(" ") );
                data.setYsIntensity( dataArr[2].split(" ") );
                //data.setYsIntensityReverse( dataArr[2].split(" ") );

                Element fragEle = eachPeptide.getChild("frag");
                String bsText = fragEle.getChildText("bs");
                String ysText = fragEle.getChildText("ys");
                String brText = fragEle.getChildText("br");
                String yrText = fragEle.getChildText("yr");

                peptide.setBsText(bsText);
                peptide.setYsText(ysText);
                peptide.setBrText(brText);
                peptide.setYrText(yrText);

                String[] bsArr = bsText.split(",");
                String[] ysArr = ysText.split(",");

                double[] bsStartMassArr = new double[bsArr.length];
                double[] bsEndMassArr = new double[bsArr.length];
                double[] ysStartMassArr = new double[ysArr.length];
                double[] ysEndMassArr = new double[ysArr.length];

                for(int i=0;i<bsArr.length;i++)
                {
                    String[] strArr = bsArr[i].split(" ");
                    bsStartMassArr[i] = Double.parseDouble(strArr[1]);
                    bsEndMassArr[i] = Double.parseDouble(strArr[2]);

                    strArr = ysArr[i].split(" ");
                    ysStartMassArr[i] = Double.parseDouble(strArr[1]);
                    ysEndMassArr[i] = Double.parseDouble(strArr[2]);
                }

                data.setBsStartMass(bsStartMassArr);
                data.setBsEndMass(bsEndMassArr);
                data.setYsStartMass(ysStartMassArr);
                data.setYsEndMass(ysEndMassArr);

                peptide.addData(data);


                protein.addPeptide(peptide);
            }

	    if(protein.getPeptideList().size() == 0)
		protein.setRedundant(true);

            list.add(protein);
        }

	//populate all redundant proteins, which has no peptides
	int listSize = list.size();
	for(int i=0;i<listSize;i++)
	{
	    if(i>=listSize-1)
		break;

	    ChroProtein p = list.get(i);

	    if(p.getPeptideList().size()<=0)
		p.setPeptideList( list.get(i+1).getPeptideList() );
	}

        return list;

    }

    public ArrayList<ChroProtein> getProteinList() throws IOException, Exception
    {
        return getProteinList(null, 0);
    }

    public ArrayList<ChroProtein> getProteinList(ProgressTask task, int initSize) throws IOException, Exception
    {

        ArrayList<ChroProtein> list = new ArrayList<ChroProtein>();

        ChroProtein protein;
        ChroPeptide peptide;
        ChroData data;

        List proList = rootEle.getChildren("protein");

        int proSize = proList.size();
        int leftSize = 100 - initSize;
        float seg = leftSize/(float)proSize;
        float progress = (float)initSize;

        for(Iterator<Element> itrPro=proList.iterator(); itrPro.hasNext(); )
        {
            Element eachProtein = itrPro.next();

            if(null != task)
            {
                progress += seg;
                task.updateProgress((int)progress);
            }

            protein = new ChroProtein();
            protein.setLocus( eachProtein.getAttributeValue("locus") );
            protein.setSeqCount( eachProtein.getAttributeValue("seq_ct") );
            protein.setSpectrumCount( eachProtein.getAttributeValue("spec_ct") );


            String ltmpSpec = eachProtein.getAttributeValue("lspec_ct");
            //System.out.println("======" + ltmpSpec);

            if(null != ltmpSpec)
                protein.setLspectrumCount(ltmpSpec);

            String htmpSpec = eachProtein.getAttributeValue("hspec_ct");
            if(null != htmpSpec)
                protein.setHspectrumCount(htmpSpec);

            protein.setSeqCoverage( eachProtein.getAttributeValue("seq_cov") );
            protein.setLength( eachProtein.getAttributeValue("length") );
            protein.setMolWt( eachProtein.getAttributeValue("molwt") );
            protein.setPI( eachProtein.getAttributeValue("pi") );
            protein.setValidation( eachProtein.getAttributeValue("val") );
            protein.setDescription( eachProtein.getAttributeValue("desc") );


	    Element redEle = eachProtein.getChild("redundant");
	    if(null != redEle) {
		List redunList = redEle.getChildren("protein");

		for(Iterator<Element> itrRPro=redunList.iterator(); itrRPro.hasNext(); ) {

		    Element eachRPro = itrRPro.next();

		    ChroProtein rprotein = new ChroProtein();
		    rprotein.setLocus( eachRPro.getAttributeValue("locus") );
		    rprotein.setSeqCount( eachRPro.getAttributeValue("seq_ct") );
		    rprotein.setSpectrumCount( eachRPro.getAttributeValue("spec_ct") );
		    rprotein.setSeqCoverage( eachRPro.getAttributeValue("seq_cov") );
		    rprotein.setLength( eachRPro.getAttributeValue("length") );
		    rprotein.setMolWt( eachRPro.getAttributeValue("molwt") );
		    rprotein.setPI( eachRPro.getAttributeValue("pi") );
		    rprotein.setValidation( eachRPro.getAttributeValue("val") );
		    rprotein.setDescription( eachRPro.getAttributeValue("desc") );

		    protein.addRedunProtein(rprotein);
		}
	    }

            List pepList = eachProtein.getChildren("peptide");

	    long freeMemory = Runtime.getRuntime().freeMemory();

            //System.out.println(freeMemory);
	    if(freeMemory<10000) //virtual memory is available less than 1M
		throw new Exception("Memory is not enough. Increase virtual memory in the bat file in Windows or command in Linux. " + freeMemory + " memory left");

            for(Iterator<Element> itrPep=pepList.iterator(); itrPep.hasNext(); )
            {
                Element eachPeptide = itrPep.next();

                peptide = new ChroPeptide();
                peptide.setUnique( !"".equals(eachPeptide.getAttributeValue("unique")) );
                peptide.setFileName( eachPeptide.getAttributeValue("file") );
                peptide.setScanNum( Integer.parseInt( (null!=eachPeptide.getAttributeValue("scan"))?eachPeptide.getAttributeValue("scan"):eachPeptide.getAttributeValue("start_scan") ) );
                peptide.setSequence( eachPeptide.getAttributeValue("seq") );
                peptide.setXCorr( eachPeptide.getAttributeValue("xcorr") );
                peptide.setChargeState( eachPeptide.getAttributeValue("charge") );
                peptide.setDeltCN( eachPeptide.getAttributeValue("deltaCN") );
                peptide.setSpecCount( eachPeptide.getAttributeValue("spC") );
                peptide.setDeltMass( eachPeptide.getAttributeValue("deltaMass") );
                peptide.setSpRank(eachPeptide.getAttributeValue("spRank") );
                peptide.setSpScore(eachPeptide.getAttributeValue("spScore") );
                String purity = eachPeptide.getAttributeValue("tmt_purity");
                String signalNoise = eachPeptide.getAttributeValue("signal-noise");
                String reportIonSignalNoise = eachPeptide.getAttributeValue("report-ion-signal-noise");
                if(purity!=null)peptide.setTmtPurity(Double.parseDouble(purity));
             //   if(signalNoise!=null)peptide.setSignalNoise(Double.parseDouble(signalNoise));
                if(reportIonSignalNoise!=null)peptide.setSignalNoise(Double.parseDouble(reportIonSignalNoise));
                String enrichment = eachPeptide.getAttributeValue("enrichment");
                if(null != enrichment)
                    peptide.setEnrichment(Double.parseDouble(enrichment));
                String calcEnrichment = eachPeptide.getAttributeValue("calc_enrich");
                if(null != calcEnrichment)
                    peptide.setEnrichment(Double.parseDouble(calcEnrichment));

                String enrichTmp = eachPeptide.getAttributeValue("enrich_corr");
                if(null != enrichTmp)
                    peptide.setBestEnrichCorr(Double.parseDouble(enrichTmp));

                        enrichTmp = eachPeptide.getAttributeValue("bestEnrichDelCN");
                if(null != enrichTmp)
                    peptide.setBestEnrichDelCN(Double.parseDouble(enrichTmp));
                        enrichTmp = eachPeptide.getAttributeValue("corrOnePlus");

                if(null != enrichTmp)
                    peptide.setCorrOnePlus(Double.parseDouble(enrichTmp));
                        enrichTmp = eachPeptide.getAttributeValue("corrOneMinus");
                if(null != enrichTmp)
                    peptide.setCorrOneMinus(Double.parseDouble(enrichTmp));
                        //System.out.println(eachPeptide.getAttributeValue("spC") );


                String ptmIndexValue = eachPeptide.getAttributeValue("ptmIndex");
                if(null != ptmIndexValue)
                    peptide.setPtmIndex(ptmIndexValue);

                ptmIndexValue = eachPeptide.getAttributeValue("ptmIndexProtein");
                if(null != ptmIndexValue)
                    peptide.setPtmIndexProtein(ptmIndexValue);

                enrichTmp = eachPeptide.getAttributeValue("corrOneMinus");
                if(null != enrichTmp)
                    peptide.setCorrOneMinus(Double.parseDouble(enrichTmp));



                addScores(peptide, eachPeptide);

                String lMass = eachPeptide.getAttributeValue("lightAvgMass");
                String hMass = eachPeptide.getAttributeValue("heavyAvgMass");
                if(null != lMass)
                    peptide.setLightMass(Double.parseDouble(lMass));
                if(null != hMass)
                    peptide.setHeavyMass(Double.parseDouble(hMass));


                //we will use exp type only in the future.  No more many if else.. checking quantlevel or labeled check.

                if(this.expType>0)
                {
                    String[] dataArr;
                    String[] massMonitorArr;
                    String[] tempArr;
                    String[] peakArr;

//		if(true)
//		    continue;
                    switch(this.expType)
                    {

                        case CensusConstants.MSMS_DATA_INDEPENDENT : //data independent
                            peptide.setDtaStartRange( Integer.parseInt(eachPeptide.getAttributeValue("start_scan")) ); //For dtaselect range
                            peptide.setDtaEndRange( Integer.parseInt(eachPeptide.getAttributeValue("end_scan")) ); //for dtaselect range

                            dataArr = eachPeptide.getChildText("chro").split(";");
                            peakArr = dataArr[0].split(" ");
                            peptide.setStartRange( peakArr[1] );
                            peptide.setEndRange( peakArr[2] );


                            for(int i=1;i<dataArr.length;i++)
                            {
                                tempArr = dataArr[i].split(",");
                                String[] intArr = tempArr[0].split(" ");

                                data = new ChroData(Integer.parseInt(intArr[0]), Long.parseLong(intArr[1]), Long.parseLong(intArr[2]));

                                //data.setBsIntensity( tempArr[1].split(" ") );
                                //data.setYsIntensity( tempArr[2].split(" ") );
                                //data.setBrIntensity( tempArr[3].split(" ") );
                                //data.setYrIntensity( tempArr[4].split(" ") );

                                data.setBsIntensity( tempArr[1] );
                                data.setYsIntensity( tempArr[2] );
                                data.setBrIntensity( tempArr[3] );
                                data.setYrIntensity( tempArr[4] );

                                peptide.addData(data);
                            }


                            break;

                        case CensusConstants.MSMS_QUANT: //ms2 quant cbamberg

                            Element fragEle = eachPeptide.getChild("frag");

                            String bsText = fragEle.getChildText("bs");
                            String ysText = fragEle.getChildText("ys");
                            String brText = fragEle.getChildText("br");
                            String yrText = fragEle.getChildText("yr");

                            peptide.setBsText(bsText);
                            peptide.setYsText(ysText);
                            peptide.setBrText(brText);
                            peptide.setYrText(yrText);

                            break;

                        case CensusConstants.MSMS_SPECIFIC_SINGLE_MASS : //iTRAQ
                            peptide.setStartRange(eachPeptide.getAttributeValue("start_scan") );
                            peptide.setEndRange(eachPeptide.getAttributeValue("end_scan")) ;
                            /*  example
                             <protein locus="Reverse_P53365|ARFP2_HUMAN" seq_ct="1" spec_ct="1" seq_cov="4.4%" length="341" molwt="37856" pi="6.0"
                                val="U" desc="Arfaptin-2 (ADP-ribosylation factor-interacting protein 2) (Partner of RAC1) (Protein POR1) - Homo sapiens (Human) ">
                                <peptide unique="" file="040607_angio_glufib_1to1-03_itms.01926.01926.2" scan="01926"
                                <chro>114 117,1926.0 0.0 0.0</chro>
                                </peptide>
                            </protein> */

                            dataArr = eachPeptide.getChildText("chro").split(",");
                            massMonitorArr = dataArr[0].split(" "); //iTRAQ user defined masses to monitor

                            peptide.setMassMonitorArr(massMonitorArr);
                            List<ReportIon> ionList = new ArrayList<ReportIon>();

                            for(String each:massMonitorArr) {
                                ReportIon ri = new ReportIon();
                                ri.setMass(Double.parseDouble(each));
                                ionList.add(ri);

                            }
                            conf.setReportIonList(ionList);
//			    System.out.println("==" + Arrays.asList(massMonitorArr));

                            tempArr = dataArr[1].split(" ");
                            data = new ChroiTRAQLabelData();
                            data.setFullScanData(tempArr);
                            peptide.addData(data);

                            break;

                        case CensusConstants.MSMS_SPECIFIC_MULTIPLE_MASS : //iTRAQ
                            peptide.setDtaStartRange( Integer.parseInt(eachPeptide.getAttributeValue("start_scan")) ); //For dtaselect range
                            peptide.setDtaEndRange( Integer.parseInt(eachPeptide.getAttributeValue("end_scan")) ); //for dtaselect range

                            dataArr = eachPeptide.getChildText("chro").split(";");
                            peakArr = dataArr[0].split(" ");
                            peptide.setStartRange( peakArr[1] );
                            peptide.setEndRange( peakArr[2] );

                            massMonitorArr = dataArr[1].split(" "); //iTRAQ user defined masses to monitor
                            peptide.setMassMonitorArr(massMonitorArr);
                            //conf.setMsmsMassArr(Arrays.asList(massMonitorArr));

                            List<ReportIon> ionList2 = new ArrayList<ReportIon>();

                            for(String each:massMonitorArr) {
                                ReportIon ri = new ReportIon();
                                ri.setMass(Double.parseDouble(each));
                                ionList2.add(ri);
                            }

                            conf.setReportIonList(ionList2);


                            for(int i=2;i<dataArr.length;i++)
                            {
                                tempArr = dataArr[i].split(" ");

                                data = new ChroiTRAQLabelData();
                                data.setFullScanData(tempArr);
                                peptide.addData(data);
                            }

                            break;

                        default :
                            break;
                    }

                }
                else if(quantLevel==1 && !labeled)
                {
                    peptide.setDtaStartRange( (int)Double.parseDouble(eachPeptide.getAttributeValue("start_scan")) ); //For dtaselect range
                    peptide.setDtaEndRange( (int)Double.parseDouble(eachPeptide.getAttributeValue("end_scan")) ); //for dtaselect range
                    if(null != eachPeptide.getAttributeValue("start_rt"))
                        peptide.setStartRt( Double.parseDouble(eachPeptide.getAttributeValue("start_rt")) ); //for dtaselect range
                    else
                        peptide.setStartRt(0.0); //for dtaselect range

                    if(null != eachPeptide.getAttributeValue("end_rt"))
                        peptide.setEndRt( Double.parseDouble(eachPeptide.getAttributeValue("end_rt")) ); //for dtaselect range
                    else
                        peptide.setEndRt(0.0); //for dtaselect range

                    if(null != eachPeptide.getAttributeValue("peak_area"))
                        peptide.setPeakArea(Double.parseDouble(eachPeptide.getAttributeValue("peak_area")));

                  //  peptide.setPeakSigma( Double.parseDouble(eachPeptide.getAttributeValue("peak_sigma")) );
                  //  peptide.setPeakx( Double.parseDouble(eachPeptide.getAttributeValue("peak_x")) );
                  //  peptide.setPeaky( Double.parseDouble(eachPeptide.getAttributeValue("peak_y")) );

                   peptide.setChroData(eachPeptide.getChildText("chro"));
                   peptide.setGaussianPeakString(eachPeptide.getChildText("peaks"));

                    String[] dataArr = eachPeptide.getChildText("chro").split(";");
                    String[] peakArr = dataArr[0].split(" ");

                    peptide.setStartRange( peakArr[1] );
                    peptide.setEndRange( peakArr[2] );


                    /*
                    for(int i=1;i<dataArr.length;i++)
                    {

                        String[] tempArr = dataArr[i].split(" ");
                        //String[] intArr = tempArr[0].split(" ");

                        //data = new ChroNonLabelData();
                       // data = new ChroData();
			//data.setFullScanData(tempArr);
                       // peptide.addData(data);


                    } */


                }
                else if(quantLevel==2 && !labeled) //Typical MRM data supposed to be, but in this case it is for msms quantification of multiple samples with alignment
                {
                    peptide.setDtaStartRange( Integer.parseInt(eachPeptide.getAttributeValue("start_scan")) ); //For dtaselect range
                    peptide.setDtaEndRange( Integer.parseInt(eachPeptide.getAttributeValue("end_scan")) ); //for dtaselect range

                    String[] dataArr = eachPeptide.getChildText("chro").split(";");
                    String[] peakArr = dataArr[0].split(" "); //peak range

                    peptide.setStartRange( peakArr[1] );
                    peptide.setEndRange( peakArr[2] );

                    for(int i=1;i<dataArr.length;i++)
                    {
                        String[] eachScanArr = dataArr[i].split(",");

//                        String[] intArr = tempArr[0].split(" ");

                        //data = new ChroNonLabelData(Integer.parseInt(intArr[0]), Long.parseLong(intArr[1]), Long.parseLong(intArr[2]));
                        data = new ChroNonLabelMSMSData(); //Integer.parseInt(intArr[0]), Long.parseLong(intArr[1]), Long.parseLong(intArr[2]));

			data.setTandemData(eachScanArr);

                        peptide.addData(data);
                    }
                }
                else if(quantLevel==2 && labeled) //MRM data with labeling
                {

                    peptide.setDtaStartRange( Integer.parseInt(eachPeptide.getAttributeValue("start_scan")) ); //For dtaselect range
                    peptide.setDtaEndRange( Integer.parseInt(eachPeptide.getAttributeValue("end_scan")) ); //for dtaselect range
                    String[] dataArr = eachPeptide.getChildText("chro").split(";");
                    String[] peakArr = dataArr[0].split(" ");
                    peptide.setStartRange( peakArr[1] );
                    peptide.setEndRange( peakArr[2] );
                        for(int i=1;i<dataArr.length;i++)
                        {
                            String[] tempArr = dataArr[i].split(",");
                            String[] intArr = tempArr[0].split(" ");

                            data = new ChroData(Integer.parseInt(intArr[0]), Long.parseLong(intArr[1]), Long.parseLong(intArr[2]));

                            if(tempArr.length<4)  //new census_chro.xml format
                            {

                            }
                            else //old census_chro.xml format
                            {
                                data.setBsIntensity( tempArr[1].split(" ") );
                                data.setYsIntensity( tempArr[2].split(" ") );
                                data.setBrIntensity( tempArr[3].split(" ") );
                                data.setYrIntensity( tempArr[4].split(" ") );
                            }

                            peptide.addData(data);
                        }

                }
                else
                {
                    peptide.setDtaStartRange( (int)Double.parseDouble(eachPeptide.getAttributeValue("start_scan")) ); //For dtaselect range
                    peptide.setDtaEndRange( (int)Double.parseDouble(eachPeptide.getAttributeValue("end_scan")) ); //for dtaselect range
                    String[] dataArr = eachPeptide.getChildText("chro").split(";");
                    String[] peakArr = dataArr[0].split(" ");


//--added Harshil Shah
                    String[] listOfEntries = null;

                    String theoMassStr = eachPeptide.getChildText("theo-mass");

// ISO tag generation...
                    if(null != theoMassStr) {
                        listOfEntries = theoMassStr.split(":");// lowMass:highmass
                        double[] tempData;
                        TheoryData theData=new TheoryData();
                        tempData = new double[listOfEntries[0].split(",").length];
                        for(int i=0;i<listOfEntries[0].split(",").length;i++) // lightmass
                            tempData[i]=Double.parseDouble(listOfEntries[0].split(",")[i]);
                        theData.setLightMass(tempData);
                        tempData = new double[listOfEntries[1].split(",").length];
                        for(int i=0;i<listOfEntries[1].split(",").length;i++) //heavy mass.........
                            tempData[i]=Double.parseDouble(listOfEntries[1].split(",")[i]);
                        theData.setHeavyMass(tempData);

                        listOfEntries = eachPeptide.getChildText("theo-int").split(":");// llightInternsity : hightintebsity...
                        tempData = new double[listOfEntries[0].split(",").length];
                        for(int i=0;i<listOfEntries[0].split(",").length;i++) // lightmass
                            tempData[i]=Double.parseDouble(listOfEntries[0].split(",")[i]);
                        theData.setLightIntensity(tempData);
                        tempData = new double[listOfEntries[1].split(",").length];
                        for(int i=0;i<listOfEntries[1].split(",").length;i++) //heavy mass.........
                            tempData[i]=Double.parseDouble(listOfEntries[1].split(",")[i]);
                        theData.setHeavyIntensity(tempData);

                        peptide.setTheoryData(theData);
                        double[] normLightData = peptide.getTheoryData().getNormLightIntensity();
                        double[] normHeavyData = peptide.getTheoryData().getNormHeavyIntensity();

                        listOfEntries = eachPeptide.getChildText("iso").split(";");
                        List<IsoData> listIsoData= new ArrayList();
                        List<IsoData> listOrigianlIsoData= new ArrayList();
                        boolean isFirstLine = true;
                        for(String currrentEntry : listOfEntries)
                        {
                            if(isFirstLine)
                            {
                                isFirstLine = false;
                                continue;
                            }
                            String[] currentIsoData = currrentEntry.split(":");//scanNum:data1:data2.........
                            int scanNumber = Integer.parseInt(currentIsoData[0]);
                            String[] data1 = currentIsoData[1].split(" ");
                            String[] data2 = currentIsoData[2].split(" ");
                            double[] lightData = new double[data1.length];
                            double[] heavyData = new double[data2.length];
                            double[] lightOrigData = new double[data1.length];
                            double[] heavyOrigData = new double[data2.length];

                            for(int i =0 ;i<data1.length;i++)
                            {
                                lightData[i] = Double.parseDouble(data1[i]) * normLightData[i];
                                lightOrigData[i] = Double.parseDouble(data1[i]);

                            }
                            for(int i =0 ;i<data2.length;i++)
                            {
                                heavyData[i] = Double.parseDouble(data2[i]) * normHeavyData[i];
                                heavyOrigData[i] = Double.parseDouble(data2[i]);
                            }
//>>>>>>> 1.41



                            IsoData currentData = new IsoData(scanNumber,lightData,heavyData);
                            IsoData currentOriginalData = new IsoData(scanNumber,lightOrigData, heavyOrigData);
                            listIsoData.add(currentData);
                            listOrigianlIsoData.add(currentOriginalData);
                        }
                        peptide.setIsoDataList(listIsoData);
                        peptide.setIsoOrigDataList(listOrigianlIsoData);

                    }
//                        }
                    /*
// this is only for the array which will have a list of array each ith value of each scan NUmber..


                     listOfEntries = eachPeptide.getChildText("iso").split(";");
                      //scan:lightInt : heavyInt ; scan:lightInt : heavyInt
                    double[][] isoData = new double[listOfEntries[0].split(":")[1].split(" ").length*2][listOfEntries.length];

                    double[][] isoLightData = new double[listOfEntries[0].split(":")[1].split(" ").length][listOfEntries.length];
                    double[][] isoHeavyData = new double[listOfEntries[0].split(":")[2].split(" ").length][listOfEntries.length];
                    int counter=0;

                    for(String currrentEntry : listOfEntries)
                    {
                        String[] currentIsoData = currrentEntry.split(":");//scanNum:data1:data2.........
                        int scanNumber = Integer.parseInt(currentIsoData[0]);
                        int length = currentIsoData[1].split(" ").length;
                        String[] data1 = currentIsoData[1].split(" ");
                        String[] data2 = currentIsoData[2].split(" ");

                        for(int i=0;i<length-1;i++)
                        {
                            isoLightData[i][counter] = Double.parseDouble(data1[i]);
                            isoHeavyData[i][counter] = Double.parseDouble(data2[i]);
//                            isoData[i][counter] = Double.parseDouble(data1[i]);
//                            isoData[length+i][counter] = Double.parseDouble(data2[i]);
                        }
                        counter++;
                    }

          */


//-- till here.....
                    peptide.setStartRange( peakArr[1] );
                    peptide.setEndRange( peakArr[2] );

                    //if(isDataDependent)
                    {
                        for(int i=1;i<dataArr.length;i++)
                       {
                            String[] eachArr = dataArr[i].split(" ");

                             if("3plexMS1Labeling".equals(this.quantType) || "MultipleMs1Labeling".equals(this.quantType)) {
				data = new ChroData(eachArr, this.quantType);

			    } else if(eachArr.length>5)
			    {
				data = new ChroData(
					Integer.parseInt(eachArr[0]),
					//Long.parseLong(eachArr[1]),
					//Long.parseLong(eachArr[2]),
					(long)Double.parseDouble(eachArr[1]),
					(long)Double.parseDouble(eachArr[2]),
					Double.parseDouble(eachArr[7]),
					Double.parseDouble(eachArr[8]),
					Integer.parseInt(eachArr[3]), //total iso sam
					Integer.parseInt(eachArr[5]), //found iso sam
					Integer.parseInt(eachArr[4]), //total iso ref
					Integer.parseInt(eachArr[6]) //found iso ref
					);
			    }
			    else
			    {
				data = new ChroData(
					Integer.parseInt(eachArr[0]),
					//Long.parseLong(eachArr[1]),
					//Long.parseLong(eachArr[2])
					(long)Double.parseDouble(eachArr[1]),
					(long)Double.parseDouble(eachArr[2])
					);


			    }




                            peptide.addData(data);
                        }

                    }
                    /*
                    else
                    {
                        for(int i=1;i<dataArr.length;i++)
                        {
                            String[] tempArr = dataArr[i].split(",");
                            String[] intArr = tempArr[0].split(" ");

                            data = new ChroData(Integer.parseInt(intArr[0]), Long.parseLong(intArr[1]), Long.parseLong(intArr[2]));

                            if(tempArr.length<4)  //new census_chro.xml format
                            {

                            }
                            else //old census_chro.xml format
                            {
                                data.setBsIntensity( tempArr[1].split(" ") );
                                data.setYsIntensity( tempArr[2].split(" ") );
                                data.setBrIntensity( tempArr[3].split(" ") );
                                data.setYrIntensity( tempArr[4].split(" ") );
                            }

                            peptide.addData(data);
                        }
                    }*/
                }

                protein.addPeptide(peptide);

            }

	    if(protein.getPeptideList().size() == 0)
		protein.setRedundant(true);

            list.add(protein);
        }

	//populate all redundant proteins, which has no peptides
	int listSize = list.size();
	for(int i=0;i<listSize;i++)
	{
	    if(i>=listSize-1)
		break;

	    ChroProtein p = list.get(i);

	    if(p.getPeptideList().size()<=0)
		p.setPeptideList( list.get(i+1).getPeptideList() );
	}

        return list;
    }

    public boolean isDataDependent()
    {
        return this.isDataDependent;
    }
    /*
    public Iterator<ChroProtein> getProteins() throws IOException {
        return new Iterator<ChroProtein>() {
            private ChroProtein protein;
            private ChroPeptide peptide;

            public boolean hasNext() {
                return lastLine != null; // && !lastLine.startsWith("\tProteins\t");
            }

            public ChroProtein next() {

                try {
                    protein = getProtein();
                }
                catch (IOException e) {
                    e.printStackTrace();
                }

                return protein;
            }

            public void remove() {
                throw new UnsupportedOperationException("Not supported");
            }

            private ChroProtein getProtein() throws IOException {

                //String[] strArr = lastLine.split("\t");
                //protein = new Protein(strArr);
                //System.out.println(lastLine);
                lastLine=br.readLine();

                protein = new ChroProtein(lastLine);
                peptide = null;

                while( (lastLine=br.readLine())!=null && !lastLine.startsWith("[PRO") )
                {
                    if(lastLine.startsWith("[PEP"))
                    {
                        if(null != peptide)
                            protein.addPeptide(peptide);

                        peptide = new ChroPeptide(lastLine=br.readLine());

                        lastLine=br.readLine();//Read away [CHROMATOGRAMS] line
                        String[] peakRange = (lastLine=br.readLine()).split("\t"); //read away Peak line;
                        if(peakRange.length==3)
                        {
                            peptide.setStartRange(peakRange[1]);
                            peptide.setEndRange(peakRange[2]);
                        }
                        //Read P line here
                        //lastLine=br.readLine();//Read the first data line

                    }
                    else// if(!lastLine.startsWith("[CHR"))
                    {
                        String[] str=lastLine.split("\t");
                        data = new ChroData(Integer.parseInt(str[0]), Integer.parseInt(str[1]), Integer.parseInt(str[2]) );
                        peptide.addData(data);
                    }
                }

                protein.addPeptide(peptide);

                return protein;
            }
        };
    }
       */

    public static void main(String args[]) throws Exception
    {
        ChroXmlReader cr = new ChroXmlReader( args[0] );

        ArrayList<ChroProtein> list = cr.getProteinList();

        System.out.println(list.size());
        System.out.println(list.get(0).getProteinLine());

//        List pepList = list.get(0).getPeptideList();

        ChroPeptide peptide;
        for(int i=0;i<list.size();i++)
        {
		ChroProtein protein = (ChroProtein)list.get(i);
		List pepList = protein.getPeptideList();
		for(int j=0;j<pepList.size();j++) {
			peptide = (ChroPeptide)pepList.get(j);
                        //peptide.getF
			//System.out.println( peptide.getDataList().get(0) );
			System.out.println( peptide );


		}
        }

	ArrayList sampleList = cr.getSampleList();
	for(Iterator<String> itr=sampleList.iterator(); itr.hasNext(); )
	{
	    String each = itr.next();
	}
    }

    public int getQuantLevel() {
        return quantLevel;
    }

    public void setQuantLevel(int quantLevel) {
        this.quantLevel = quantLevel;
    }

    public boolean isLabeled() {
        return labeled;
    }

    public void setLabeled(boolean labeled) {
        this.labeled = labeled;
    }

    public ArrayList<String> getFileList() {
        return fileList;
    }

    public void setFileList(ArrayList<String> fileList) {
        this.fileList = fileList;
    }

    public String getSampleName(String fileName)
    {
        return this.fileSampleHt.get(fileName);
    }

    public String getSampleName(int fileIndex)
    {
        String curDir = getFileList().get(fileIndex);
        return getSampleName(curDir);
    }


    public String getFileName(int index)
    {
        String tmpFileName = this.fileList.get(index);

        return tmpFileName.substring( tmpFileName.lastIndexOf(File.separator)+1 );
    }

    public ArrayList getSampleList() {
        return sampleList;
    }

    public void setSampleList(ArrayList sampleList) {
        this.sampleList = sampleList;
    }

    public Hashtable<String, Sample> getSampleObjList() {
        return sampleObjList;
    }

    public void setSampleObjList(Hashtable<String, Sample> sampleObjList) {
        this.sampleObjList = sampleObjList;
    }

    public ArrayList<Sample> getSampleExpList() {
        return sampleExpList;
    }

    public void setSampleExpList(ArrayList<Sample> sampleExpList) {
        this.sampleExpList = sampleExpList;
    }

    public int getExpType() {
        return expType;
    }

    public void setExpType(int expType) {
        this.expType = expType;
    }

    public void addScores(ChroPeptide chroPeptide, Element peptide) {

	List scrList = peptide.getChildren("search_score");
	for(Iterator<Element> scrItr=scrList.iterator(); scrItr.hasNext(); )
	{
	    Element eachScr = scrItr.next();
	    chroPeptide.addScore( eachScr.getAttributeValue("name"), eachScr.getAttributeValue("value") );
	}
    }

    public String getQuantType() {
        return quantType;
    }

    public void setQuantType(String quantType) {
        this.quantType = quantType;
    }
    public static HashMap<String, ChroPeptide> getPeptideMap(String fileName)
    {//key is the Sequence_cs
        ChroXmlReader cr = null;
        HashMap<String, ChroPeptide> peptideMap = new HashMap<>();
        try {
             cr = new ChroXmlReader( fileName );
             for(ChroProtein currentPRotein : cr.getProteinList())
             {
                 for(ChroPeptide currentPeptide :(List<ChroPeptide>) currentPRotein.getPeptideList())
                 {
                     peptideMap.put(currentPeptide.getSequence()+"_" + currentPeptide.getChargeState(), currentPeptide);

                 }
             }

        } catch (JDOMException ex) {
            Logger.getLogger(ChroXmlReader.class.getName()).log(Level.SEVERE, null, ex);
        } catch (Exception ex) {
            Logger.getLogger(ChroXmlReader.class.getName()).log(Level.SEVERE, null, ex);
        }
        return peptideMap;
    }


}


