/*
 * Configuration.java
 *
 * Created on March 18, 2005, 3:51 PM
 */

package edu.scripps.pms.census.conf;

import java.util.*;
import edu.scripps.pms.census.CensusConstants;

import java.io.*;

import gnu.trove.TDoubleArrayList;
import edu.scripps.pms.census.util.StringUtil;
import edu.scripps.pms.census.model.*;
import edu.scripps.pms.census.model.mrm.*;

import java.text.*;

import org.jdom.*;
import org.jdom.input.*;

import edu.scripps.pms.census.ChroProgressDialog;
import edu.scripps.pms.census.hash.IndexedFile;
import edu.scripps.pms.census.io.IsotopeReader;

/**
 *
 * @author  Robin Park
 * @version $Id: Configuration.java,v 1.57 2014/06/06 03:41:41 rpark Exp $
 *
 */


/*

version 2.51 For generating TMT census-out.txt file, we made memory efficient
 *
 *
 */
//Singleton object
public class Configuration {
    private static Configuration conf;
    private double isolationWindow;
    private double startMassRange;
    private double endMassRange;
    private int numOfIsolationWindow=1;
    private double[] precursorArr;
    private double[] windowArr=null;;
    private double massRange;
    private double enrichment;
    private double enrichmentMaxDeviation;
    private String ms2ScanTypeValue ="";
    private int margin=20;
    private int steepArea=3;
    private float steepRatioThreshold=0.95f;
    private double isobaricIsolationWindow = 0.8;
//    private float thresholdIntensity=0.5f;
    //private int maxWindow=2000;
    private int maxWindow=50;
    private boolean isHighRes=false; //high resolution
    private boolean printLog=false;
    public final int HIGH_RES_THRESHOLD=5000;


    private boolean isSimpleIndexGenerator=false;
    private boolean chargeColumn=false;
    private double massTolerance; //=0.03; //for high resolution
    private String filePath;
    private int resolution=0;
    private long startTime;
    private boolean luciphorMode = false;

    private boolean isDataIndependent=false;
    private boolean isLabeling=true;
    private boolean readConfigFile = false;
    private String errorMessage = null;
    private String outputFilename = null;
    private String quantType= null;
    private double retentionTimeWindow = -1;
    //private JLabel progressLabel;

    private TDoubleArrayList precursorList = new TDoubleArrayList();

    private String[] fileList;
    private List<SampleGroup> sampleGroupList = new ArrayList<>();

    private int numFileSize;

    private static SimpleDateFormat timerFormat = new SimpleDateFormat("HH:mm:ss");

    private ArrayList<Sample> expList = new ArrayList<Sample>();  //list of sample classes


    private String refFileName;  //reference file name for non-labeling quantification
    private int quantLevel=1;

    private Set nonlabelFilePaths;
    private List nonlabelFilenameList;
    private Map<String,List<String>> nonlabelFilenameGroupMap = new HashMap<String,List<String>>();
    private String configFilePath;

    //private Vector<String> pathFileNameList;
    private NonLabelMappingModel mapModel;

    //private String dtaSelectFile;
    private String idFileName = "DTASelect-filter.txt"; //default file name is DTASelect-filter.txt

    //extraction method
    private int extMethod; //1 for whole extraction, and 2 for individual peaks
    private Element rootConfEle; //root config element
    private double calcSamAvgMass; //average mass
    private double calcRefAvgMass; //average mass
    private double calcRef2AvgMass; //average mass

    private ChroProgressDialog progressDialog;

    //file number, data size
    private double retArr[][] = null;
    private boolean align=true;
    private boolean xmlConf=false;

    private String version = "2.54.2";
    private String tmtOutlierLevel = "NO_OUTLIER";


    /*******************  EXPORT PARAMETERS *********************/
    private double outlierPValue = 0.1; //outlier p value. default is 0.1
    private int specCountNormal = 0;

    private Hashtable<String, int[]> spHt = new Hashtable<String, int[]>(); //spec count hashtable

    public static final int SPEC_COUNT_NO_NORMALIZATION = 0;
    public static final int SPEC_COUNT_NORMALIZATION1 = 1;

    private int spectrumFormat = 0;
    private String mzXMLFilePath;
    public static final int MS_FILE_FORMAT = 0; //ms1 file format
    public static final int MS2_FILE_FORMAT = 2; //ms2 file format
    public static final int MZXML_FILE_FORMAT = 1;

//    public static int SPEC_COUNT_NO_NORMALIZATION = 0;
    /*******************  END OF EXPORT PARAMETERS **************/


    /*******************  MSMS SPECIFIC PARAMETERS (e.g. iTRAQ) *********************/
    public static final int AUTOMATIC_FRAGMENT_ION=1;
    public static final int SPECIFIC_FRAGMENT_ION=2;

    private int msmsFragType=-1;
    //private List msmsMassArr = new ArrayList();
    private double msmsSpecificTolerance;
    private boolean msmsTolerancePpm=false;

    public static final int MSMS_SINGLE_SPECTRUM=1;
    public static final int MSMS_MULTIPLE_SPECTRA=2;

    private int msmsSpectrumNum;
    public static final double MRM_PRECURSOR_TOLERANCE=0.05; //This works for both data dependent and MRM experiments
//Added by harshil
    private List experimentGroup = new ArrayList();
    private List experimentGroupAttribute = new ArrayList();
    private IsotopeReader isotopeReader = null;
    //Experiment Types
    public static final int MRM_WITHOUT_ID=20;

    public int expType=-1;

    private List<MRMPeptideGroup> mrmPeptideGroupList = null;
    private int mrmPeptideSize=0;

    /*******************  END OF MSMS SPECIFIC PARAMETERS (e.g. iTRAQ) **************/


    private boolean useProline = false; // for proline from SILAC
    private int prolineCount= 1; // # of prolines to consider
    public static double PROLINE_SHIFT = 6.013804;

    private int scoreType=DTASELECT;
    public static final int DTASELECT=1;
    public static final int PEPXML=2;

    private ArrayList<String> scoreList = new ArrayList<String>();

    private boolean basedOnId = false; //label free analysis

    private List<SampleModel> sampleList = new ArrayList<SampleModel>();
    private String databaseFile;

    //N15 enrichment calculation
    private boolean calculateEnrich=false;

    private double startEnrich; //15N enrichment range
    private double endEnrich; //15N enrichment range

    //private boolean searchAdditionalRefMass = false; //see isuseProline
    private double[] addtionalRefMassArr = null;

    private boolean useMassDiff=false;

    //multiple peaks or single peak for label free based on id
    private String peakCount;
    private double isodistThreshold=0.1; // remove masses below 10%
    private int isotopePeak=0;
    private int scanShift=0;
    private boolean ms3ScanRandom=false;
    private String scanShiftType="HCD";
    private int sampleNum=0;

    private String fileShift;
    private List<ReportIon> reportIonList = new ArrayList<ReportIon>();  //list of sample classes

    //cbamberg ms2quant
    private double ms2WideTolerance;
    private double ms2NarrowTolerance;
    private double[] ms2Standards = null;


    private double intensityThreshold;
    private int maxSpectrumShift=0;
    private boolean pulseLabeling = false;
    private String pulseLightMass;
    private String pulseHeavyMass;
    private String pulseResidue;

    private boolean labelfree=false;
    private String ms2Label;
    private String allIDSline ="true";
    //for labelfree
    Hashtable<String, IndexedFile> indexHt = null;
    private boolean smooth=false;
    private String labelfreeFillPeptide = "true";
    private boolean labelfreeCheckChargeState=false;
    private boolean targeted=false;

    private int targetedStartCharge=2;
    private int targetedEndCharge=5;

    private boolean isTimstofXicMode = false;


    public String getAllIDSline() {
        return allIDSline;
    }
    public void setAllIDSline(String allIDSline) {
        this.allIDSline = allIDSline;
    }
    public String getLabelfreeFillPeptide() {
        return labelfreeFillPeptide;
    }

    public void setLabelfreeFillPeptide(String labelfreeFillPeptide) {
        this.labelfreeFillPeptide = labelfreeFillPeptide;
    }
    public void setOutlierPValue(double outlierPValue)
    {
	this.outlierPValue = outlierPValue;
    }

    public double getOutlierPValue()
    {
	return outlierPValue;
    }

    public void addExp(Sample sam)
    {
        this.expList.add(sam);
    }

    public boolean isBasedOnId() {
        return basedOnId;
    }

    public void setBasedOnId(boolean basedOnId) {
        this.basedOnId = basedOnId;
    }

    public List<SampleModel> getSampleList() {
        return sampleList;
    }

    public void setSampleList(List<SampleModel> sampleList) {
        this.sampleList = sampleList;
    }

    public void addSample(SampleModel sModel) {
        this.sampleList.add(sModel);
    }

    public String getDatabaseFile() {
        return databaseFile;
    }

    public void setDatabaseFile(String databaseFile) {
        this.databaseFile = databaseFile;
    }

    public boolean isCalculateEnrich() {
        return calculateEnrich;
    }

    public void setCalculateEnrich(boolean calculateEnrich) {
        this.calculateEnrich = calculateEnrich;
    }

    public double getStartEnrich() {
        return startEnrich;
    }

    public void setStartEnrich(double startEnrich) {
        this.startEnrich = startEnrich;
    }

    public double getEndEnrich() {
        return endEnrich;
    }

    public void setEndEnrich(double endEnrich) {
        this.endEnrich = endEnrich;
    }

    public double getEnrichmentMaxDeviation() {
        return enrichmentMaxDeviation;
    }

    public void setEnrichmentMaxDeviation(double enrichmentMaxDeviation) {
        this.enrichmentMaxDeviation = enrichmentMaxDeviation;
    }

    //public //15N enrichment range
    /*
    boolean isSearchAdditionalRefMass() {
        return searchAdditionalRefMass;
    }

    public void setSearchAdditionalRefMass(boolean searchAdditionalRefMass) {
        this.searchAdditionalRefMass = searchAdditionalRefMass;
    }
*/
    public double[] getAddtionalRefMassArr() {
        return addtionalRefMassArr;
    }

    public void setAddtionalRefMassArr(double[] addtionalRefMassArr) {
        this.addtionalRefMassArr = addtionalRefMassArr;
    }

    /**
     * @return the experimentGroup
     */
    public List getExperimentGroup() {
        return experimentGroup;
    }

    /**
     * @param experimentGroup the experimentGroup to set
     */
    public void setExperimentGroup(String experimentGroup) {
        String[] s= experimentGroup.split(",");
        double[] d = new double[s.length];
        for(int i=0;i<s.length ; i++)
        {
            d[i]=(double) Double.parseDouble(s[i]);
        }

        this.experimentGroup.add(d);
    }

    /**
     * @return the experimentGroupAttribute
     */
    public List getExperimentGroupAttribute() {
        return experimentGroupAttribute;
    }

    /**
     */
    public void setExperimentGroupAttribute(String name) {
        experimentGroupAttribute.add(name);
    }

    /**
     * checks whether we have read teh cofig file or not
     * @return true: if read the cofig file
     *         false: if we haven't read the  config file
     */
    public boolean isReadConfigFile() {
        return readConfigFile;
    }

    public void setReadConfigFile(boolean readConfigFile) {
        this.readConfigFile = readConfigFile;
    }

    /**
     * @return the isotopeReader
     */
    public IsotopeReader getIsotopeReader() {
        return isotopeReader;
    }

    /**
     * @param isotopeReader the isotopeReader to set
     */
    public void setIsotopeReader(IsotopeReader isotopeReader) {
        this.isotopeReader = isotopeReader;
    }

    public double getRetentionTimeWindow() {
        return retentionTimeWindow;
    }

    public void setRetentionTimeWindow(double retentionTimeWindow) {
        this.retentionTimeWindow = retentionTimeWindow;
    }

    public static class Sample {
        private String name;
        private ArrayList msFileList = new ArrayList();

        public String getName() {
            return name;
        }

        public void setName(String name) {
            this.name = name;
        }

        public ArrayList getMsFileList() {
            return msFileList;
        }

        public void setMsFileList(ArrayList msFileList) {
            this.msFileList = msFileList;
        }

        public void addFile(String fileName)
        {
            this.msFileList.add(fileName);
        }

    }


    public boolean isSimpleIndexGenerator() {
        return this.isSimpleIndexGenerator;
    }

    public void setSimpleIndexGenerator(boolean isSimpleIndexGenerator) {
        this.isSimpleIndexGenerator = isSimpleIndexGenerator;
    }

    private int colNum; //column number of mass intensity line

    public void setColNum(int colNum) {
        this.colNum = colNum;
    }

    public int getColNum() {
        return this.colNum;
    }

    public void setStartTime()
    {
	this.startTime = System.currentTimeMillis();
    }

    public long getStartTime()
    {
	return this.startTime;
    }

    public String getRunningTime()
    {
	return timerFormat.format(new java.util.Date(System.currentTimeMillis()-startTime));
    }

    public void setIsHighRes(boolean isHighRes) {
        this.isHighRes = isHighRes;
    }

    public boolean isHighRes() {
        return this.isHighRes;
    }

    private Configuration() {
    }

    public static Configuration getInstance() {
        if (conf == null)
            conf = new Configuration();

        return conf;
    }

    public double getMassRange() {
        return endMassRange - startMassRange;
    }

    public void setIsolationWindow(double isolationWindow) {
        this.isolationWindow = isolationWindow;
    }

    public double getIsolationWindow() {
        return this.isolationWindow;
    }

    public boolean isDataIndependent()
    {
        return this.isDataIndependent;
    }

    public void setIsDataIndependent(boolean isDataIndependent)
    {
        this.isDataIndependent = isDataIndependent;
    }


    public void setStartMassRange(double startMassRange) {
        this.startMassRange = startMassRange;
    }

    public void setEndMassRange(double endMassRange) {
        this.endMassRange = endMassRange;
    }

    public double getStartMassRange() {
        return this.startMassRange;
    }

    public double getEndMassRange() {
        return this.endMassRange;
    }

    public int getNumOfIsolationWindow() {
        return this.numOfIsolationWindow;
    }

    //never called from outside
    private void setNumOfIsolationWindow(int numOfIsolationWindow) {
        this.numOfIsolationWindow = numOfIsolationWindow;
    }

    public double[] getPrecursorArr() {
        return this.precursorArr;
    }

    public void setPrecursorArr(double[] precursorArr) {
        this.precursorArr = precursorArr;

        setNumOfIsolationWindow(precursorArr.length);
    }

    public double getPrecursor(double avgMass) {
        if(null == windowArr)
            calculateWindowArr();

        int i;

        for(i=0;i<this.windowArr.length;i++) {
            if(this.windowArr[i]>avgMass)
                break;
        }



        return this.precursorArr[(i-1)>0?i-1:0];// + this.isolationWindow;
        //return this.precursorArr[i]; //(i-1)>0?i-1:0] + this.isolationWindow;
    }

    public double[] getWindowArr() {
        if(null == windowArr)
            calculateWindowArr();

        return this.windowArr;
    }

    public void setWindowArr(double[] windowArr) {
        this.windowArr = windowArr;
    }

    public void calculateWindowArr() {
        this.windowArr = new double[this.precursorArr.length];

        for(int i=0;i<this.precursorArr.length;i++)
            this.windowArr[i] = this.precursorArr[i] - this.isolationWindow/2;
    }

    public double getEnrichment() {
        /*
        if(enrichment>0)
            return enrichment;

        this.enrichment = Double.parseDouble( ht.get( CensusConstants.ENRICHMENT ).toString() );
        */
        return this.enrichment;
    }

    public void setEnrichment(double enrichment) {
        this.enrichment = enrichment;
    }

    public int getMargin() {
        return margin;
    }

    public void setMargin(int margin) {
        this.margin = margin;
    }

    public int getSteepArea() {
        return steepArea;
    }

    public void setSteepArea(int steepArea) {
        this.steepArea = steepArea;
    }

    public float getSteepRatioThreshold() {
        return steepRatioThreshold;
    }

    public void setSteepRatioThreshold(float steepRatioThreshold) {
        this.steepRatioThreshold = steepRatioThreshold;
    }

   /*
    public float getThresholdIntensity() {
        return thresholdIntensity;
    }

    public void setThresholdIntensity(float thresholdIntensity) {
        this.thresholdIntensity = thresholdIntensity;
    }
    */

    public int getMaxWindow() {
        return maxWindow;
    }

    public void setMaxWindow(int maxWindow) {
        this.maxWindow = maxWindow;
    }

    public boolean isChargeColumn() {
        return chargeColumn;
    }

    public void setChargeColumn(boolean chargeColumn) {
        this.chargeColumn = chargeColumn;
    }

    public double getMassTolerance() {
        return massTolerance;
    }

    public double getMassToleranceInPPM() {
      return massTolerance*1000;
    }
    public void setMassTolerance(double massTolerance) {
        this.massTolerance = massTolerance;
    }

    public void readSimpleXml(String fileName) throws Exception {

        SAXBuilder builder = new SAXBuilder();

        Document doc = builder.build(new File(fileName));
        Element rootEle = doc.getRootElement();

        this.rootConfEle = rootEle;

    }

    public void readMRMParams(Element mrmParamEle)
    {
        if(null == mrmParamEle)
            return;

        List mrmPeptideList = mrmParamEle.getChildren("peptide"); //null;

	/*
        if(null != mrmParamEle && "false".equals(mrmParamEle.getAttributeValue("id")) )
	{
	    this.expType = MRM_WITHOUT_ID;
	    this.spectrumFormat = MS2_FILE_FORMAT;
            mrmPeptideList = mrmParamEle.getChildren("peptide_group");
	}*/

        mrmPeptideGroupList = new ArrayList<MRMPeptideGroup>();

	for(Iterator<Element> gitr=mrmPeptideList.iterator(); gitr.hasNext(); )
	{
	    Element pepGroupEle = gitr.next();
	    MRMPeptideGroup mrmGroup = new MRMPeptideGroup();

	    for(Iterator<Element> pitr=pepGroupEle.getChildren("precursor").iterator(); pitr.hasNext(); )
	    {
		Element pepEach = pitr.next();

		MRMPeptideModel peptide = new MRMPeptideModel();
                mrmPeptideSize++;

		//String value = pepEach.getAttributeValue("sequence");
		//if(null != value && !"".equals(value))
		//    peptide.setSequence(value);

		String value = pepEach.getAttributeValue("mass");
		if(null != value && !"".equals(value))
		    peptide.setParentMass( Double.parseDouble(value) );

		value = pepEach.getAttributeValue("name");
		if(null != value && !"".equals(value))
		    peptide.setName(value);

		value = pepEach.getAttributeValue("desc");
		if(null != value && !"".equals(value))
		    peptide.setDesc(value);

		/*
                value = pepEach.getAttributeValue("label");
                if(null != value && "true".equals(value))
                    peptide.setLabeled(true);
                else if(null != value && "false".equals(value))
                    peptide.setLabeled(false);
		*/
                /*
                peptide.setScanNum( pepEach.getAttributeValue("scan_num") );
                peptide.setRt( pepEach.getAttributeValue("rt") );
                peptide.setFileName( pepEach.getAttributeValue("file_name") );

                value = mrmParamEle.getAttributeValue("sn_tolerance");

                if(null != value)
                    peptide.setSnTolerance( Double.parseDouble(value) );

                value = mrmParamEle.getAttributeValue("rt_tolerance");

                if(null != value)
                peptide.setRtTolerance( Double.parseDouble(value) );
                */
                Element daughterEle = pepEach.getChild("daughter");
                List<Element> daughterList = Collections.EMPTY_LIST;

                if(null != daughterEle)
                   daughterList = daughterEle.getChildren("mass");

		for(Iterator<Element> ditr=daughterList.iterator(); ditr.hasNext(); )
		{
		    Element massEle = ditr.next();
		    peptide.addDaughter( Double.parseDouble(massEle.getText()), Double.parseDouble(massEle.getAttributeValue("rt")) );
		}

                mrmGroup.addPeptide(peptide);
	    }

	    mrmPeptideGroupList.add(mrmGroup);


	}

    }

    /**
     * Reads config xml file
     * @param fileName
     * @throws Exception
     */
    public void readXMLParam(String fileName) throws Exception {

        xmlConf = true;
        SAXBuilder builder = new SAXBuilder();

        Document doc = builder.build(new File(fileName));
        Element rootEle = doc.getRootElement();
        this.rootConfEle = rootEle;

        Element paramEle = rootEle.getChild("params");
        Element qtypeEle = rootEle.getChild("quant_type");
        Element labelfreeEle = paramEle.getChild("labelfree");

        Element luciphorMode = rootEle.getChild("luciphor_mode");
        if(luciphorMode !=null)
        {
            this.luciphorMode = "true".equalsIgnoreCase(luciphorMode.getAttributeValue("readMode"));
        }

	this.sampleNum = rootEle.getChild("element_comp").getChildren("each_sample").size();
	if(null != qtypeEle)
		conf.setQuantType(qtypeEle.getText() );
	else { //for back compatibility we need to check if it is 15N or not

		List<Element> sList = rootEle.getChild("element_comp").getChildren("each_sample");

		if(sList.size()>1) {
			Element hEle = sList.get(1);

			List<Element> rList = hEle.getChildren("residue");
			boolean n15allvalue=true;

      for (Iterator<Element> rItr = rList.iterator(); rItr.hasNext(); ) {
        Element eachR = rItr.next();

        String tstr = eachR.getAttributeValue("name");

        if ("NTERM".equals(tstr) ||
          "CTERM".equals(tstr) ||
          "*".equals(tstr) ||
          "#".equals(tstr) ||
          "@".equals(tstr))
          continue;

        if ("0".equals(eachR.getChildText("ele_15N")))
          n15allvalue = false;
      }
      if (n15allvalue)
        conf.setQuantType("15N");
    }

	}

        //labelfree parser
         List<Element> sampleEleList = rootEle.getChildren("sample");
         if(labelfree)
            readLabelfreeSampleModel(sampleEleList);
        String timstofXicMode = paramEle.getChildText("timstof_xic");
        if(timstofXicMode !=null)
        {
            this.isTimstofXicMode = timstofXicMode.equalsIgnoreCase("true");
        }

        String ms2LabelTxt = paramEle.getChildText("ms2_label");
        if(null != ms2LabelTxt) {
          Element msLabelEle = paramEle.getChild("ms2_label");
          String ms2ScanTypeValue = msLabelEle.getAttributeValue("scan_type");
          if(null != ms2ScanTypeValue)
          {
            //assign this.scanTypeValue =
              this.ms2ScanTypeValue = ms2ScanTypeValue;
              this.scanShiftType = ms2ScanTypeValue;
          }
            this.ms2Label = ms2LabelTxt;
        }
        String smoothtxt = paramEle.getChildText("smooth");
        if(null != smoothtxt) {
            this.smooth = Boolean.valueOf(smoothtxt);
        }

        if (labelfreeEle != null) {
            String labelfree_missing_peptide = labelfreeEle.getChildText("find_missing_peptide");
            if (null != labelfree_missing_peptide) {
                this.labelfreeFillPeptide = labelfree_missing_peptide;
            }

            String allSlineOnly = labelfreeEle.getChildText("all_id_sline_only");
            if (null != labelfree_missing_peptide) {
                this.allIDSline = allSlineOnly;
            }

            String cs_check = labelfreeEle.getChildText("charge_state_check");
            if (null != cs_check && !"false".equals(cs_check) && "true".equals(cs_check)) {
                this.labelfreeCheckChargeState = true;
            }

        }

        Element targetEle = paramEle.getChild("targeted");
        if(null != targetEle) {
            setTargeted(true);
            String chargeStart = targetEle.getAttributeValue("charge_start");
            String chargeEnd = targetEle.getAttributeValue("charge_end");
            if(chargeStart != null)
                targetedStartCharge = Integer.parseInt(chargeStart);
            if(chargeEnd != null)
                targetedEndCharge = Integer.parseInt(chargeEnd);

        }
        String massTol = paramEle.getChildText("isobaric_isolation_window");
        if(massTol!=null) isobaricIsolationWindow = Double.parseDouble(massTol);

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

        String labelingType = rootEle.getChild("label_type").getAttributeValue("type");
        if(null != labelingType && "pulse".equals(labelingType)) {
            //add pulse reader
            Element pulseLabelingEle = rootEle.getChild("label_type").getChild("pulse_label");
            if(null != pulseLabelingEle) {
		pulseLabeling = true;
                pulseLightMass = pulseLabelingEle.getAttributeValue("light");
                pulseHeavyMass = pulseLabelingEle.getAttributeValue("heavy");
                pulseResidue = pulseLabelingEle.getAttributeValue("residue");
            }
        }


        String useMassDiff = rootEle.getChild("label_type").getAttributeValue("useMassDiff");
        if("true".equals(useMassDiff))
            this.useMassDiff = true;
        else
            this.useMassDiff = false;

        String expTypeStr = rootEle.getChildText("experiment_type");

        if(null != expTypeStr  && !"".equals(expTypeStr))
            this.expType = Integer.parseInt(expTypeStr);


      if(null != expTypeStr && expTypeStr.equals("17")) // DIA
        quantType = "DIA_LF";



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
		String scanShiftTypeStr = msmsParams.getAttributeValue("scan_shift_type");
               // fileShift = msmsParams.getAttributeValue("ms3");  //file_shift tag was removed.  this ms3 is redundnat, but keep them for back compatibility

                String ms3ScanDefined = msmsParams.getAttributeValue("ms3");
                if(null != ms3ScanDefined && "random".equals(ms3ScanDefined)) {
                        this.ms3ScanRandom = true;
                        this.fileShift ="ms3";
                }

		if(null != scanShiftStr)
			scanShift = Integer.parseInt(scanShiftStr);
		if(null != scanShiftTypeStr)
			scanShiftType = scanShiftTypeStr;

                Element msmsFragParams = paramEle.getChild("msms_params").getChild("frag_ext_type");
                Element tmtOutlierLevel = paramEle.getChild("msms_params").getChild("outlier_level");

                if(null != tmtOutlierLevel)
                    this.tmtOutlierLevel = tmtOutlierLevel.getValue();

		if(null != msmsFragParams) {
			if("single".equals(msmsFragParams.getChildText("spectrum")))
			{
				this.msmsSpectrumNum = this.MSMS_SINGLE_SPECTRUM;
			}
			else if("multiple".equals(msmsFragParams.getChildText("spectrum")) )
			{
				this.msmsSpectrumNum = this.MSMS_MULTIPLE_SPECTRA;
			}
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

                            ReportIon ri = new ReportIon();
                            ri.setMass( Double.parseDouble(massEach.getText()) );
                            String cpList = massEach.getAttributeValue("cp_list");
                            String cmList = massEach.getAttributeValue("cm_list");

                            if(null != cpList) ri.setcPlusList(cpList);
                            if(null != cmList) ri.setcMinusList(cmList);
			    reportIonList.add(ri);

                            //getMsmsMassArr().add( massEach.getText() );
                        }
//Added by Harshil
                        if(msmsFragParams.getChildren("group")!= null)
                        {
                            List l1= msmsFragParams.getChildren("group");

                            int counter =0;
                            for(Iterator<Element> itr=l1.iterator(); itr.hasNext(); )
                            {
                                Element currentGroup = itr.next();
                                if((String)currentGroup.getAttributeValue("name") == null)
                                {
                                    setExperimentGroupAttribute(Integer.toString(counter));
                                    counter++;
                                }
                                else
                                {
                                    setExperimentGroupAttribute(currentGroup.getAttributeValue("name"));
                                }

                                setExperimentGroup(currentGroup.getText());

                            }
                        }
 //Till here- Harshil
                        Element msmsToleranceEle = msmsFragParams.getChild("msms_tolerance");
                        String unit = msmsToleranceEle.getAttributeValue("unit");

                        if("ppm".equals(unit))
                            this.msmsTolerancePpm=true;
                        else
                            this.msmsTolerancePpm=false;

                        setMsmsSpecificTolerance(Double.parseDouble(msmsToleranceEle.getValue()));
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
        if(paramEle.getChild("retentionTime_window") != null)
            this.retentionTimeWindow = Double.parseDouble(paramEle.getChild("retentionTime_window").getValue());


	tmpMassTol = paramEle.getChildText("ms2_wide_masstolerance");

	if(null != tmpMassTol)
	    this.ms2WideTolerance = Double.parseDouble( tmpMassTol );

        tmpMassTol = paramEle.getChildText("ms2_narrow_masstolerance");
        if(null != tmpMassTol)
	    this.ms2NarrowTolerance = Double.parseDouble( tmpMassTol );

      String tmpMs2Standard = paramEle.getChildText("ms2_standards");
      if(null != tmpMs2Standard) {
        String[] arr = tmpMs2Standard.split(",");
        this.ms2Standards = new double[arr.length];
        for(int i=0;i<arr.length;i++) {
          this.ms2Standards[i] = Double.parseDouble(arr[i]);
        }
      }


        String tmpInt = paramEle.getChildText("intensity_threshold");
        if(null != tmpInt)
	    this.intensityThreshold= Double.parseDouble( tmpInt );

	String tmpIsodistThreshold = paramEle.getChildText("isodist_threshold");

	if(null != tmpIsodistThreshold)
	    this.isodistThreshold = Double.parseDouble( tmpIsodistThreshold );


/*

            if("ppm".equals(pReader.getExtMethod())) {
                this.highRes = true;
                this.massAccuracy = pReader.getMassTolerance()*1000;
            } else {
                this.highRes = false;
                this.envelopMargin = pReader.getMassTolerance();
            }

*/
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


      //read header
      this.setIsotopeReader(new IsotopeReader(rootEle));

        this.readConfigFile = true;
    }

    /**
     * //this method is deprecated.  read param as non-xml format.
     * @param filePath
     * @param fileName
     * @throws Exception
     */

    public void readParam(String filePath, String fileName) throws Exception {
        //ht.addParamLine(paramLine);
        xmlConf = false;

        this.filePath = filePath;

        File f = new File(filePath);
	fileList = f.list();

        numFileSize = this.fileList.length;

        BufferedReader br = new BufferedReader(new FileReader(filePath + File.separator + fileName));

        Hashtable<String, String> paramTable = new Hashtable<String, String>();

        //read header
        String lastLine = br.readLine();

        while( (lastLine=br.readLine())!=null )
        {
            ht.addParamLine(lastLine);
        }

        if(null == ht.get(CensusConstants.DATA_TYPE))
            throw new Exception("Data type must be defined.");

        int tempInt = Integer.parseInt( (String)ht.get(CensusConstants.DATA_TYPE) );
        if(tempInt==0) this.isDataIndependent=false;
        else this.isDataIndependent=true;
        if(null != ht.get(CensusConstants.ISOLATION_WINDOW))
            this.isolationWindow = Double.parseDouble( (String)ht.get(CensusConstants.ISOLATION_WINDOW) );

        if(null != ht.get(CensusConstants.START_RANGE))
            this.startMassRange = Double.parseDouble( (String)ht.get(CensusConstants.START_RANGE) );

        if(null != ht.get(CensusConstants.END_RANGE))
            this.endMassRange = Double.parseDouble( (String)ht.get(CensusConstants.END_RANGE) );

        if(null != ht.get(CensusConstants.MASS_ACCURACY))
            this.massTolerance = Double.parseDouble((String)ht.get(CensusConstants.MASS_ACCURACY))*0.001;

        if(null != ht.get(CensusConstants.ENRICHMENT))
            this.enrichment = Double.parseDouble((String)ht.get(CensusConstants.ENRICHMENT));

        if(null != ht.get(CensusConstants.MAX_WINDOW))
            this.maxWindow = Integer.parseInt((String)ht.get(CensusConstants.MAX_WINDOW));

        if(null != ht.get(CensusConstants.WINDOW_MARGIN))
            this.margin = Integer.parseInt((String)ht.get(CensusConstants.WINDOW_MARGIN));

        this.resolution = Integer.parseInt( (String)ht.get(CensusConstants.RESOLUTION ) );
	if( conf.HIGH_RES_THRESHOLD <= resolution )
	    this.isHighRes=true;

        if(this.isDataIndependent)
        {

            ArrayList<Double> aList = readHlines();
	    double[] winArr = new double[aList.size()];
	    double eachWidth = this.isolationWindow/2;
	    double[] precurArr = new double[winArr.length];

            int index=0;

	    for(Iterator<Double> itr=aList.iterator(); itr.hasNext(); )
	    {
		precurArr[index] = itr.next().doubleValue();
		winArr[index] = precurArr[index] - eachWidth;

		index++;
	    }

	    if(precurArr.length<=0)
	    {
                double startPrecursor = this.startMassRange + this.isolationWindow/2;

		for(double i=startPrecursor; i<=this.endMassRange; i+=this.isolationWindow)
		    precursorList.add(i);

		double[] arr = precursorList.toNativeArray();
		conf.setPrecursorArr(arr);
	    } else {

		conf.setPrecursorArr(precurArr);
		conf.setWindowArr(winArr);
	    }

	    haveChargeStateColumn("ms2");
	}
	else
        {
            haveChargeStateColumn("ms1");
        }

        br.close();
    }

    private ParamTable ht = new ParamTable();

    //Read isolation window from H lines
    //return true, if it contains isolation window ranges
    private ArrayList<Double> readHlines() throws IOException
    {
	String eachLine;

	ArrayList<Double> aList = new ArrayList<Double>();

	int size=0;

	BufferedReader br=null;

	for(int i=0;i<fileList.length;i++)
	{
	    if(fileList[i].endsWith("ms2"))
	    {
		br = new BufferedReader(new FileReader(filePath + File.separator + fileList[i]));

		while ( (eachLine = br.readLine()) != null && eachLine.startsWith("H\t"))// some H lines are followed by space, not tab. So, don't use H\t. :(
		{
		    if(!eachLine.startsWith("H\t"))
			break;

		    if(eachLine.startsWith("H\tScanFilter"))
		    {
			int tempIndex = eachLine.indexOf('@');
			if(tempIndex<0)
			    continue;
			eachLine = eachLine.substring(0, eachLine.indexOf('@') );
			eachLine = eachLine.substring( eachLine.lastIndexOf(' ')+1 );

			aList.add(new Double(eachLine));
			size++;
		    }
		}

		while ( (eachLine = br.readLine()) != null && eachLine.startsWith("I\t"));
		while ( (eachLine = br.readLine()) != null && eachLine.startsWith("Z\t"));

		conf.setColNum( eachLine.split(" ").length );

		break;
	    }
	}

	return aList;
    }


    private void haveChargeStateColumn(String extension) throws IOException
    {
	String eachLine=null;

	    for(int i=0;i<fileList.length;i++)
	    {
		if(fileList[i].endsWith(extension))
		{
		    BufferedReader br = new BufferedReader(new FileReader(filePath + File.separator + fileList[i]));
		    while ( (eachLine = br.readLine()) != null && !StringUtil.startsWithDigit(eachLine) );

		    if(eachLine.split(" ").length==3)
			conf.setChargeColumn(true);
		    else
			conf.setChargeColumn(false);

		    break;
		}

	}
    }


    class ParamTable extends Hashtable {
        public void addParamLine(String paramLine) {
            if(null == paramLine || paramLine.startsWith("H\t") || "".equals(paramLine.trim()))
                return;

            String[] arr = paramLine.split("\t");
            int index = arr[1].indexOf("#");

            if(index<=0)
                ht.put(arr[0], arr[1].trim());
            else
                ht.put(arr[0], arr[1].substring(0, index).trim());

        }
    }

    public String getFilePath() {
	if(null == filePath)
	    return null;

        //file could be generated from either linux or window
        if( !this.filePath.endsWith("/") && !this.filePath.endsWith("\\") ) {
            //for linux
            if(this.filePath.startsWith("/") || this.filePath.equals("."))
                this.filePath += "/";
            else //for window
                this.filePath += "\\";
        }

        return filePath;
    }

    public void setFilePath(String filePath) {
        this.filePath = filePath;
    }

    public String getElementCompFile() {
        return (String)ht.get( CensusConstants.ELEMENT_COMPOSITION_FILE );
    }

    public int getResolution() {
        return resolution;
    }

    public void setResolution(int resolution) {
        this.resolution = resolution;
    }

    public String getRefFileName() {
        return refFileName;
    }

    public void setRefFileName(String refFileName) {
        this.refFileName = refFileName;
    }

    public int getQuantLevel() {
        return quantLevel;
    }

    public void setQuantLevel(int quantLevel) {
        this.quantLevel = quantLevel;
    }

    public Set getNonlabelFilePaths() {
        return nonlabelFilePaths;
    }

    public void setNonlabelFilePaths(Set nonlabelFilePaths) {
        this.nonlabelFilePaths = nonlabelFilePaths;
    }

    public boolean isLabeling() {
        return isLabeling;
    }

    public void setLabeling(boolean isLabeling) {
        this.isLabeling = isLabeling;
    }

    public NonLabelMappingModel getMapModel() {
        return mapModel;
    }

    public void setMapModel(NonLabelMappingModel mapModel) {
        this.mapModel = mapModel;
    }

    public int getNumFileSize() {
        return numFileSize;
    }

    public void setNumFileSize(int numFileSize) {
        this.numFileSize = numFileSize;
    }

    public int getExtMethod() {
        return extMethod;
    }

    public void setExtMethod(int extMethod) {
        this.extMethod = extMethod;
    }

        /*
    public String getDtaSelectFile() {
        return dtaSelectFile;
    }


    public void setDtaSelectFile(String dtaSelectFile) {
        this.dtaSelectFile = dtaSelectFile;
    }
*/

    public Element getRootConfEle() {
        return rootConfEle;
    }

    public void setRootConfEle(Element rootConfEle) {
        this.rootConfEle = rootConfEle;
    }

    public double getCalcSamAvgMass() {
        return calcSamAvgMass;
    }

    public void setCalcSamAvgMass(double calcSamAvgMass) {
        this.calcSamAvgMass = calcSamAvgMass;
    }

    public double getCalcRefAvgMass() {
        return calcRefAvgMass;
    }

    public void setCalcRefAvgMass(double calcRefAvgMass) {
        this.calcRefAvgMass = calcRefAvgMass;
    }

    public double getCalcRef2AvgMass() {
        return calcRef2AvgMass;
    }

    public void setCalcRef2AvgMass(double calcRef2AvgMass) {
        this.calcRef2AvgMass = calcRef2AvgMass;
    }

    /*
    public JLabel getProgressLabel() {
        return progressLabel;
    }

    public void setProgressLabel(JLabel progressLabel) {
        this.progressLabel = progressLabel;
    }*/

    public ChroProgressDialog getProgressDialog() {
        return progressDialog;
    }

    public void setProgressDialog(ChroProgressDialog progressDialog) {
        this.progressDialog = progressDialog;
    }

    public double[][] getRetArr() {
        return retArr;
    }

    public void setRetArr(double[][] retArr) {
        this.retArr = retArr;
    }

    public boolean isAlign() {
        return align;
    }

    public void setAlign(boolean align) {
        this.align = align;
    }

    public String getErrorMessage() {
        return errorMessage;
    }

    public void setErrorMessage(String errorMessage) {
        this.errorMessage = errorMessage;
    }

    public boolean isXmlConf() {
        return xmlConf;
    }

    public void setXmlConf(boolean xmlConf) {
        this.xmlConf = xmlConf;
    }

    public String getVersion() {
        return version;
    }

    public void setVersion(String version) {
        this.version = version;
    }

    public int getSpecCountNormal() {
        return specCountNormal;
    }

    public void setSpecCountNormal(int specCountNormal) {
        this.specCountNormal = specCountNormal;
    }

    public List getNonlabelFilenameList() {
        return nonlabelFilenameList;
    }

    public void setNonlabelFilenameList(List nonlabelFilenameList) {
        this.nonlabelFilenameList = nonlabelFilenameList;
    }

    public Hashtable<String, int[]> getSpHt() {
        return spHt;
    }

    public void setSpHt(Hashtable<String, int[]> spHt) {
        this.spHt = spHt;
    }

    public int getMsmsFragType() {
        return msmsFragType;
    }

    public void setMsmsFragType(int msmsFragType) {
        this.msmsFragType = msmsFragType;
    }

    public double getMsmsSpecificTolerance() {
        return msmsSpecificTolerance;
    }

    public void setMsmsSpecificTolerance(double msmsSpecificTolerance) {
        this.msmsSpecificTolerance = msmsSpecificTolerance;
    }

    public int getMsmsSpectrumNum() {
        return msmsSpectrumNum;
    }

    public void setMsmsSpectrumNum(int msmsSpectrumNum) {
        this.msmsSpectrumNum = msmsSpectrumNum;
    }

    public int getExpType() {
        return expType;
    }

    public void setExpType(int expType) {
        this.expType = expType;
    }

    public int getSpectrumFormat() {
        return spectrumFormat;
    }

/*
    public String getSpectrumFormatAsString() {

	String format=null;

	switch(spectrumFormat)
	{
	    case Configuration.MS_FILE_FORMAT:
		format = CensusConstants.MS1_FILE;
		break;

	    case Configuration.MS2_FILE_FORMAT:
		format = CensusConstants.MS2_FILE;
		break;

	    case Configuration.MZXML_FILE_FORMAT:
		format = CensusConstants.MZXML;
		break;

	    default:
		break;
	}

        return format;
    }
*/

    public void setSpectrumFormat(int spectrumFormat) {
        this.spectrumFormat = spectrumFormat;
    }

    public String getMzXMLFilePath() {
        return mzXMLFilePath;
    }

    public void setMzXMLFilePath(String mzXMLFilePath) {
        this.mzXMLFilePath = mzXMLFilePath;
    }

    public String getIdFileName() {
        return idFileName;
    }

    public void setIdFileName(String idFileName) {
        this.idFileName = idFileName;
    }

    public List<MRMPeptideGroup> getMrmPeptideGroupList() {
        return mrmPeptideGroupList;
    }

    public void setMrmPeptideGroupList(List<MRMPeptideGroup> mrmPeptideGroupList) {
        this.mrmPeptideGroupList = mrmPeptideGroupList;
    }

    public int getMrmPeptideSize() {
        return mrmPeptideSize;
    }

    public void setMrmPeptideSize(int mrmPeptideSize) {
        this.mrmPeptideSize = mrmPeptideSize;
    }

    public boolean isUseProline() {
        return useProline;
    }

    public void setUseProline(boolean useProline) {
        this.useProline = useProline;
    }

    public boolean isPrintLog() {
        return printLog;
    }

    public void setPrintLog(boolean printLog) {
        this.printLog = printLog;
    }

    public void setScoreType(int scoreType) {
	this.scoreType = scoreType;
    }

    public int getScoreType() {
	return scoreType;
    }

    public void setScoreList(ArrayList scoreList) {
	this.scoreList = scoreList;
    }

    public ArrayList<String> getScoreList() {
	return scoreList;
    }

    public void addScoreName(String score) {
	this.scoreList.add(score);
    }

    public int getProlineCount() {
	return this.prolineCount;
    }


    public String getPeakCount() {
        return peakCount;
    }

    public void setPeakCount(String peakCount) {
        this.peakCount = peakCount;
    }

    public String getOutputFilename() {
        return outputFilename;
    }

    public void setOutputFilename(String outputFilename) {
        this.outputFilename = outputFilename;
    }

    public boolean isUseMassDiff() {
        return useMassDiff;
    }

    public void setUseMassDiff(boolean useMassDiff) {
        this.useMassDiff = useMassDiff;
    }

    public void setIsodistThreshold(double isodistThreshold) {
        this.isodistThreshold = isodistThreshold;
    }

    public double getIsodistThreshold() {
        return isodistThreshold;
    }

    public int getIsotopePeak() {
        return isotopePeak;
    }

    public String getQuantType() {
        return quantType;
    }

    public void setQuantType(String quantType) {
        this.quantType = quantType;
    }


    /**
     * @return the scanShift
     */
    public int getScanShift() {
        return scanShift;
    }

    /**
     * @param scanShift the scanShift to set
     */
    public void setScanShift(int scanShift) {
        this.scanShift = scanShift;
    }

    /**
     * @return the sampleNum
     */
    public int getSampleNum() {
        return sampleNum;
    }

    /**
     * @param sampleNum the sampleNum to set
     */
    public void setSampleNum(int sampleNum) {
        this.sampleNum = sampleNum;
    }

    public String getFileShift() {
        return fileShift;
    }

    public void setFileShift(String fileShift) {
        this.fileShift = fileShift;
    }

    public String getScanShiftType() {
        return scanShiftType;
    }

    public void setScanShiftType(String scanShiftType) {
        this.scanShiftType = scanShiftType;
    }

    public List<ReportIon> getReportIonList() {
        return reportIonList;
    }

    public void setReportIonList(List<ReportIon> reportIonList) {
        this.reportIonList = reportIonList;
    }

    public String getTmtOutlierLevel() {
        return tmtOutlierLevel;
    }

    public void setTmtOutlierLevel(String tmtOutlierLevel) {
        this.tmtOutlierLevel = tmtOutlierLevel;
    }

    public double getMs2WideTolerance() {
        return ms2WideTolerance;
    }

    public void setMs2WideTolerance(double ms2WideTolerance) {
        this.ms2WideTolerance = ms2WideTolerance;
    }

    public double getMs2NarrowTolerance() {
        return ms2NarrowTolerance;
    }

    public void setMs2NarrowTolerance(double ms2NarrowTolerance) {
        this.ms2NarrowTolerance = ms2NarrowTolerance;
    }

    public double getIntensityThreshold() {
        return intensityThreshold;
    }

    public void setIntensityThreshold(double intensityThreshold) {
        this.intensityThreshold = intensityThreshold;
    }

    public boolean isMs3ScanRandom() {
        return ms3ScanRandom;
    }

    public void setMs3ScanRandom(boolean ms3ScanRandom) {
        this.ms3ScanRandom = ms3ScanRandom;
    }


    public int getMaxSpectrumShift() {
        return maxSpectrumShift;
    }

    public void setMaxSpectrumShift(int maxSpectrumShift) {
        this.maxSpectrumShift = maxSpectrumShift;
    }

    public String getPulseLightMass() {
        return pulseLightMass;
    }

    public void setPulseLightMass(String pulseLightMass) {
        this.pulseLightMass = pulseLightMass;
    }

    public String getPulseHeavyMass() {
        return pulseHeavyMass;
    }

    public void setPulseHeavyMass(String pulseHeavyMass) {
        this.pulseHeavyMass = pulseHeavyMass;
    }

    public String getPulseResidue() {
        return pulseResidue;
    }

    public void setPulseResidue(String pulseResidue) {
        this.pulseResidue = pulseResidue;
    }

    public boolean isPulseLabeling() {
        return pulseLabeling;
    }

    public void setPulseLabeling(boolean pulseLabeling) {
        this.pulseLabeling = pulseLabeling;
    }

        public void readLabelfreeSampleModel(List<Element> sampleEleList) {
        for(Iterator<Element> itr = sampleEleList.iterator(); itr.hasNext(); )
            {
                Element sam = itr.next();
                //String samName = sam.getChildText("name");
                String groupName = sam.getAttributeValue("group");
                //populate configuration class
                //Configuration.Sample confSam = new Configuration.Sample();
                //confSam.setName(samName);
                SampleGroup sampleGroup = new SampleGroup();
                sampleGroup.setName(groupName);

                List<Element> tsampleList = sam.getChildren("each_sample");
                for(Iterator<Element> eachSampleItr=tsampleList.iterator(); eachSampleItr.hasNext(); ) {
                    Element eachSample = eachSampleItr.next();

                    List<Element> fileList = eachSample.getChild("ms_files").getChildren("file");
                    String samName = eachSample.getAttributeValue("name");

                    SampleModel sampleModel = new SampleModel(samName);

                    if(fileList.size()>0)
                    {
                        String firstFileName = fileList.get(0).getText();
                        String filePath = firstFileName.substring(0, firstFileName.lastIndexOf(File.separator));
                        firstFileName = firstFileName.substring(firstFileName.lastIndexOf(File.separator)+1);

                        if(firstFileName.startsWith("*"))
                        {
                            String extension = firstFileName.substring(firstFileName.lastIndexOf(".")+1);

                            edu.scripps.pms.census.util.RelExFileFilter fFilter = new edu.scripps.pms.census.util.RelExFileFilter(extension);


                            File specFile = new File(filePath);
                            String[] splist = specFile.list(new edu.scripps.pms.census.util.RelExFileFilter(extension));


                            for(String eachFileName : splist)
                            {

                                //eachFileName = filePath + File.separator + eachFileName;
                                //set.add(filePath);
                                sampleModel.addPath(filePath);
                                sampleModel.addLabelfreeFilename(eachFileName);


                                //confSam.addFile(eachFileName);
                               // pathFileNameList.add(eachFileName);
                                //fileNameList.add(eachFileName);
                                //sampleNameList.add(samName);
                            }
                        }
                        else {
                            for(Iterator<Element> itr1 = fileList.iterator(); itr1.hasNext(); )
                            {
                                Element eachFile = itr1.next();
                                String fileName = eachFile.getText();

                                //confSam.addFile(fileName);
                               // pathFileNameList.add(fileName);

                               // set.add(fileName.substring(0, fileName.lastIndexOf(File.separator)));
                                sampleModel.addPath(fileName.substring(0, fileName.lastIndexOf(File.separator)));

                                if(fileName.endsWith("ms2"))
                                {
                                    fileName = fileName.substring(0, fileName.length()-3);
                                    fileName += "ms1";
                                }

                               // fileNameList.add(fileName);
                               // sampleNameList.add(samName);
//
                            }
                        }

                        sampleGroup.addSample(sampleModel);
                     //   conf.addExp(confSam);
                        addSample(sampleModel);
                    //List<Element> fileList = sam.getChild("ms_files").getChildren("file");
                    }



                }

         //     System.out.println("=====>>>==================\t" + sampleGroup.getName());
              addSampleGroup(sampleGroup);
            }

        System.out.println("conf done");
    }


    public Hashtable<String, IndexedFile> getIndexHt() {
        return indexHt;
    }

    public void setIndexHt(Hashtable<String, IndexedFile> indexHt) {
        this.indexHt = indexHt;
    }

    public boolean isLabelfree() {
        return labelfree;
    }

    public void setLabelfree(boolean labelfree) {
        this.labelfree = labelfree;
    }


    public Map<String, List<String>> getNonlabelFilenameGroupMap() {
        return nonlabelFilenameGroupMap;
    }

    public void setNonlabelFilenameGroupMap(Map<String, List<String>> nonlabelFilenameGroupMap) {
        this.nonlabelFilenameGroupMap = nonlabelFilenameGroupMap;
    }

    public String getConfigFilePath() {
        return configFilePath;
    }

    public void setConfigFilePath(String configFilePath) {
        this.configFilePath = configFilePath;
    }

    public boolean isMsmsTolerancePpm() {
        return msmsTolerancePpm;
    }

    public void setMsmsTolerancePpm(boolean msmsTolerancePpm) {
        this.msmsTolerancePpm = msmsTolerancePpm;
    }

    public String getMs2Label() {
        return ms2Label;
    }

    public void setMs2Label(String ms2Label) {
        this.ms2Label = ms2Label;
    }

    public boolean isSmooth() {
        return smooth;
    }

    public void setSmooth(boolean smooth) {
        this.smooth = smooth;
    }

    public boolean isLabelfreeCheckChargeState() {
        return labelfreeCheckChargeState;
    }

    public void setLabelfreeCheckChargeState(boolean labelfreeCheckChargeState) {
        this.labelfreeCheckChargeState = labelfreeCheckChargeState;
    }

    public static void main(String[] args) throws Exception {
	        Configuration conf = Configuration.getInstance();

     // conf.readXMLParam("/data/2/rpark/ip2_data/negraosf/Label_free_experiments/labelfree_quant/labelfree_16441/census_config_labelfree_16441.xml");
        conf.readXMLParam("/home/rpark/test_data/pfizer/lfree_john/rerun/config.xml");

      Element rootEle = conf.getRootConfEle();
      List<Element> sampleList = rootEle.getChildren("sample");
      for(Element sampleEle:sampleList) {
        List<Element> eachSamList = sampleEle.getChildren("each_sample");
        for(Element eachSample:eachSamList) {
          List<Element> msfileList = eachSample.getChildren("ms_files");
          for(Element msfile:msfileList) {
            List<Element> fileList = msfile.getChildren("file");
            for(Element file:fileList) {
              String path = file.getText();
              String searchId = path.substring(path.lastIndexOf("_")+1);
              searchId = searchId.substring(0, searchId.indexOf("/"));
              System.out.println(searchId);
            }
          }

        }
      }

      if(true) return;


        if (!conf.isReadConfigFile()) {
            conf.setLabelfree(true);
            //conf.readXMLParam("/data/2/rpark/ip2_data/rpark/Megan/default_params/census_config_targeted.xml");
          conf.readXMLParam("/home/rpark/test_data/targeted_lfree/census_config_targeted_21.xml");
            //conf.readXMLParam("/home/rpark/test_data/targeted_lfree/census_config_labelfree_target.xml");
        }

        System.out.println(conf.getSampleGroupList());
      System.out.println(conf.getSampleGroupList().size());


        //List<SampleModel> sampleList = conf.getSampleList();

        System.out.println("charge state check:\t" + conf.isLabelfreeCheckChargeState());


    }

    public boolean isTargeted() {
        return targeted;
    }

    public void setTargeted(boolean targeted) {
        this.targeted = targeted;
    }

    public void addSampleGroup(SampleGroup sampleGroup) {
        sampleGroupList.add(sampleGroup);
    }
    public List<SampleGroup> getSampleGroupList() {
        return sampleGroupList;
    }

    public void setSampleGroupList(List<SampleGroup> sampleGroupList) {
        this.sampleGroupList = sampleGroupList;
    }

    public int getTargetedStartCharge() {
        return targetedStartCharge;
    }

    public void setTargetedStartCharge(int targetedStartCharge) {
        this.targetedStartCharge = targetedStartCharge;
    }

    public int getTargetedEndCharge() {
        return targetedEndCharge;
    }

    public void setTargetedEndCharge(int targetedEndCharge) {
        this.targetedEndCharge = targetedEndCharge;
    }
    /**
     * Reads config xml file but does not open directories
     * @param fileName
     * @throws Exception
     */
    public void softReadXMLParam(String fileName) throws Exception {

        xmlConf = true;
        SAXBuilder builder = new SAXBuilder();

        Document doc = builder.build(new File(fileName));
        Element rootEle = doc.getRootElement();
        this.rootConfEle = rootEle;

        //read header
        this.setIsotopeReader(new IsotopeReader(rootEle));
        Element paramEle = rootEle.getChild("params");
        Element qtypeEle = rootEle.getChild("quant_type");
        Element labelfreeEle = paramEle.getChild("labelfree");

        this.sampleNum = rootEle.getChild("element_comp").getChildren("each_sample").size();
        if(null != qtypeEle)
            conf.setQuantType(qtypeEle.getText() );
        else { //for back compatibility we need to check if it is 15N or not

            List<Element> sList = rootEle.getChild("element_comp").getChildren("each_sample");

            if(sList.size()>1) {
                Element hEle = sList.get(1);

                List<Element> rList = hEle.getChildren("residue");
                boolean n15allvalue=true;
                for(Iterator<Element> rItr=rList.iterator(); rItr.hasNext(); )
                {
                    Element eachR = rItr.next();

                    String tstr = eachR.getAttributeValue("name");

                    if("NTERM".equals(tstr) ||
                            "CTERM".equals(tstr) ||
                            "*".equals(tstr) ||
                            "#".equals(tstr) ||
                            "@".equals(tstr) )
                        continue;

                    if("0".equals(eachR.getChildText("ele_15N")))
                        n15allvalue=false;
                }
                if(n15allvalue)
                    conf.setQuantType("15N");
            }

        }

        //labelfree parser
        List<Element> sampleEleList = rootEle.getChildren("sample");
      /*  if(labelfree)
            readLabelfreeSampleModel(sampleEleList);
*/

        String ms2LabelTxt = paramEle.getChildText("ms2_label");
        if(null != ms2LabelTxt) {
            this.ms2Label = ms2LabelTxt;
        }
        String smoothtxt = paramEle.getChildText("smooth");
        if(null != smoothtxt) {
            this.smooth = Boolean.valueOf(smoothtxt);
        }

        if (labelfreeEle != null) {
            String labelfree_missing_peptide = labelfreeEle.getChildText("find_missing_peptide");
            if (null != labelfree_missing_peptide) {
                this.labelfreeFillPeptide = labelfree_missing_peptide;
            }

            String allSlineOnly = labelfreeEle.getChildText("all_id_sline_only");
            if (null != labelfree_missing_peptide) {
                this.allIDSline = allSlineOnly;
            }

            String cs_check = labelfreeEle.getChildText("charge_state_check");
            if (null != cs_check && !"false".equals(cs_check) &&  "true".equals(cs_check)) {
                this.labelfreeCheckChargeState = true;
            }

        }

        Element targetEle = paramEle.getChild("targeted");
        if(null != targetEle) {
            setTargeted(true);
            String chargeStart = targetEle.getAttributeValue("charge_start");
            String chargeEnd = targetEle.getAttributeValue("charge_end");
            if(chargeStart != null)
                targetedStartCharge = Integer.parseInt(chargeStart);
            if(chargeEnd != null)
                targetedEndCharge = Integer.parseInt(chargeEnd);

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

        String labelingType = rootEle.getChild("label_type").getAttributeValue("type");
        if(null != labelingType && "pulse".equals(labelingType)) {
            //add pulse reader
            Element pulseLabelingEle = rootEle.getChild("label_type").getChild("pulse_label");
            if(null != pulseLabelingEle) {
                pulseLabeling = true;
                pulseLightMass = pulseLabelingEle.getAttributeValue("light");
                pulseHeavyMass = pulseLabelingEle.getAttributeValue("heavy");
                pulseResidue = pulseLabelingEle.getAttributeValue("residue");
            }
        }


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
                String scanShiftTypeStr = msmsParams.getAttributeValue("scan_shift_type");
                // fileShift = msmsParams.getAttributeValue("ms3");  //file_shift tag was removed.  this ms3 is redundnat, but keep them for back compatibility

                String ms3ScanDefined = msmsParams.getAttributeValue("ms3");
                if(null != ms3ScanDefined && "random".equals(ms3ScanDefined)) {
                    this.ms3ScanRandom = true;
                    this.fileShift ="ms3";
                }

                if(null != scanShiftStr)
                    scanShift = Integer.parseInt(scanShiftStr);
                if(null != scanShiftTypeStr)
                    scanShiftType = scanShiftTypeStr;

                Element msmsFragParams = paramEle.getChild("msms_params").getChild("frag_ext_type");
                Element tmtOutlierLevel = paramEle.getChild("msms_params").getChild("outlier_level");

                if(null != tmtOutlierLevel)
                    this.tmtOutlierLevel = tmtOutlierLevel.getValue();

                if(null != msmsFragParams) {
                    if("single".equals(msmsFragParams.getChildText("spectrum")))
                    {
                        this.msmsSpectrumNum = this.MSMS_SINGLE_SPECTRUM;
                    }
                    else if("multiple".equals(msmsFragParams.getChildText("spectrum")) )
                    {
                        this.msmsSpectrumNum = this.MSMS_MULTIPLE_SPECTRA;
                    }
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

                            ReportIon ri = new ReportIon();
                            ri.setMass( Double.parseDouble(massEach.getText()) );
                            String cpList = massEach.getAttributeValue("cp_list");
                            String cmList = massEach.getAttributeValue("cm_list");

                            if(null != cpList) ri.setcPlusList(cpList);
                            if(null != cmList) ri.setcMinusList(cmList);
                            reportIonList.add(ri);

                            //getMsmsMassArr().add( massEach.getText() );
                        }
//Added by Harshil
                        if(msmsFragParams.getChildren("group")!= null)
                        {
                            List l1= msmsFragParams.getChildren("group");

                            int counter =0;
                            for(Iterator<Element> itr=l1.iterator(); itr.hasNext(); )
                            {
                                Element currentGroup = itr.next();
                                if((String)currentGroup.getAttributeValue("name") == null)
                                {
                                    setExperimentGroupAttribute(Integer.toString(counter));
                                    counter++;
                                }
                                else
                                {
                                    setExperimentGroupAttribute(currentGroup.getAttributeValue("name"));
                                }

                                setExperimentGroup(currentGroup.getText());

                            }
                        }
                        //Till here- Harshil
                        Element msmsToleranceEle = msmsFragParams.getChild("msms_tolerance");
                        String unit = msmsToleranceEle.getAttributeValue("unit");

                        if("ppm".equals(unit))
                            this.msmsTolerancePpm=true;
                        else
                            this.msmsTolerancePpm=false;

                        setMsmsSpecificTolerance(Double.parseDouble(msmsToleranceEle.getValue()));
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
        if(paramEle.getChild("retentionTime_window") != null)
            this.retentionTimeWindow = Double.parseDouble(paramEle.getChild("retentionTime_window").getValue());


        tmpMassTol = paramEle.getChildText("ms2_wide_masstolerance");

        if(null != tmpMassTol)
            this.ms2WideTolerance = Double.parseDouble( tmpMassTol );

        tmpMassTol = paramEle.getChildText("ms2_narrow_masstolerance");
        if(null != tmpMassTol)
            this.ms2NarrowTolerance = Double.parseDouble( tmpMassTol );

        String tmpInt = paramEle.getChildText("intensity_threshold");
        if(null != tmpInt)
            this.intensityThreshold= Double.parseDouble( tmpInt );

        String tmpIsodistThreshold = paramEle.getChildText("isodist_threshold");

        if(null != tmpIsodistThreshold)
            this.isodistThreshold = Double.parseDouble( tmpIsodistThreshold );


/*

            if("ppm".equals(pReader.getExtMethod())) {
                this.highRes = true;
                this.massAccuracy = pReader.getMassTolerance()*1000;
            } else {
                this.highRes = false;
                this.envelopMargin = pReader.getMassTolerance();
            }

*/
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


        this.readConfigFile = true;
    }

    public String getMs2ScanTypeValue() {
        return ms2ScanTypeValue;
    }
    public boolean isMs2ScanTypeValueUnassigned() {return conf.getMs2ScanTypeValue().equals("");}

  public double[] getMs2Standards() {
    return ms2Standards;
  }

  public void setMs2Standards(double[] ms2Standards) {
    this.ms2Standards = ms2Standards;
  }

    public double getIsobaricIsolationWindow() {
        return isobaricIsolationWindow;
    }

    public void setIsobaricIsolationWindow(double isobaricIsolationWindow) {
        this.isobaricIsolationWindow = isobaricIsolationWindow;
    }

    public boolean isTimstofXicMode() {
        return isTimstofXicMode;
    }


    public boolean isLuciphorMode() {
        return luciphorMode;
    }
}
