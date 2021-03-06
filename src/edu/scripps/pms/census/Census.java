/*
 * RelaxMainFrame.java
 *
 * Created on March 18, 2005, 3:52 PM
 */

package edu.scripps.pms.census;

import com.TmtFilter.DAN.TmtfilterDan;
import edu.scripps.pms.census.labelFree.LabelfreeTargeted;
import edu.scripps.pms.census.tmtFilter.TmtfilterPeptide;
import javax.swing.*;
import java.io.*;
import java.awt.Color;
import java.util.*;
import java.text.DecimalFormat;
import java.awt.event.*;
import java.awt.Toolkit;

import edu.scripps.pms.census.util.SimpleFileNameFilter;

import edu.scripps.pms.census.io.*;
import edu.scripps.pms.census.io.parse.ChroXMLParser;

import edu.scripps.pms.census.util.RelExFileFilter;
import edu.scripps.pms.census.util.PostCalculation;
import edu.scripps.pms.census.util.CalcUtil;

import edu.scripps.pms.census.model.ChroProtein;
import edu.scripps.pms.census.model.ChroPeptide;
import edu.scripps.pms.census.model.ChroData;
import edu.scripps.pms.census.model.SampleModel;
import edu.scripps.pms.census.dialog.*;

import javax.swing.table.DefaultTableModel;

import ptolemy.plot.*;
import ptolemy.plot.plotml.PlotBoxMLParser;
import ptolemy.plot.plotml.PlotMLParser;
import edu.scripps.pms.census.plot.*;
import edu.scripps.pms.census.conf.*;

import edu.scripps.pms.census.util.LinearRegression;
import edu.scripps.pms.census.util.*;
import edu.scripps.pms.census.exception.*;
import org.jdom.*;
import org.jdom.input.*;

import edu.scripps.pms.census.chroalign.*;
import edu.scripps.pms.census.labelFree.ChroJSONReader;
import edu.scripps.pms.census.model.NonLabelMappingModel;

import org.apache.commons.cli.*;

import edu.scripps.pms.util.Mzxml2Ms;
/**
 *
 * @author  Robin Park
 * @version $Id: Census.java,v 1.39 2014/08/08 17:07:37 rpark Exp $
 */

public class Census {

    /** Creates new form inRelaxMainFrame */

    private float massTolerance;
    private final String MS1_FILE = "ms1";
    private final String MS2_FILE = "ms2";
    private Configuration conf;
    private String filePath;
    private String elementFile;
    private DecimalFormat twoDigitFormat = new DecimalFormat("0.00");

    //private final String PARAM_FILE="relex.param";
    private final String PARAM_FILE="census.param";

    public static void printHelp()
    {
    /*
	System.out.println("Welcome to Census");
        System.out.println("identification file such as DTASelect-filter.txt, pepXML or census_id.txt");
	System.out.println("-c\tconfig file name");
	System.out.println("-x\tuse mzXML files instead of MS files.");
	System.out.println("-g\tGUI mode run");
	System.out.println("-f\tspectra folder for labeling analysis");
	System.out.println("-h\tprint help");
	System.out.println("-e\texport report file");
	System.out.println("-chro\tchro xml file");
        System.out.println("-p\tAdd proline shift for SILAC");
        System.out.println("-etan\texport tandem tag analysis (TMT, iTRAQ, etc)");
        System.out.println("-int\tintensity threshold for tandem tag analysis (TMT, iTRAQ, etc)");



        //System.out.println("-e\texport file.  Specify location of chro file.  census-out.txt will be generated");

        //System.out.println("example java -jar census.jar -e -ch /data/1/census_chro.xml");
	*/
	System.out.println("example1. java -jar census.jar -c /my_folder/census_config.xml -f folder_name -i DTASelect-filter.txt");
	System.out.println("example2. java -jar census.jar -c /my_folder/census_config.xml -f folder_name");
	System.out.println("example3. java -jar census.jar -e census_report.params -f folder_name -chro census_chro.xml");
	System.out.println("example4. java -jar census.jar  -int 50 -etan census_chro.xml -f /temp");
	System.out.println("example4. java -jar census.jar  -int 50 -etan census_chro.xml -f /temp -t mean");
	System.out.println("example5. java -jar census.jar  noise -etan census_chro.xml -f /temp");
        System.out.println("-int value : This option is for isobaric labeling.  value can be hard threshold or z-score e.g. 'z2' will discard low intensity values lower than z-score 2.");
    }

    public static void main(String args[]) throws IOException, CensusGeneralException, Exception
    {
//	this.filePath = args[0];
//        System.out.println("dfsf");
        long start = System.currentTimeMillis();
        //System.out.println( start );
        Census census;

		Options opt = new Options();
		opt.addOption("c", true, "config file name with path");
		//opt.addOption("x", true, "The data source to use");
		opt.addOption("g", false, "run in GUI mode");
		opt.addOption("noise", false, "export iTRAQ/TMT based on 90% baseline");
		opt.addOption("ref", false, "normalize using 1st report ion as reference");
		opt.addOption("f", true, "path to working directory. Default is current folder.");
		opt.addOption("h", false, "print help");
        opt.addOption("e", true, "export report file");
        opt.addOption("chro", true, "chro xml file");
        opt.addOption("x", false, "mzXML files instead of MS file format.");
        opt.addOption("i", true, "identification file such as DTASelect-filter.txt, pepXML or census_id.txt");
        opt.addOption("m", false, "MRM experiment without identification");
        opt.addOption("v", false, "Census version");
        opt.addOption("p", false, "Add proline shift for SILAC");
        opt.addOption("log", false, "print log file");
        opt.addOption("etan", true, "export tandem tag analysis (TMT, iTRAQ, etc)");
        opt.addOption("int", true, "intensity threshold for tandem tag analysis (TMT, iTRAQ, etc)");
        opt.addOption("t", true, "type of normalization (TMT, iTRAQ, etc).  mean or mode");
        opt.addOption("spc", true, "Compare spec count from multiple DTASelect results.  Usage: java census.jar -spc DTASelect-filter.txt(1st) -spc DTASelect-filter.txt(2nd) ...");
	opt.addOption("aa", false, "select default values automatically for commandline interactive questions (support currently only for label free analysis)");
	opt.addOption("d", false, "read existing census_chro_temp.xml for label-free analysis");
	opt.addOption("of", true, "output file name (support currently only for label free analysis)");
        opt.addOption("lf", true, "labelfree Analysis run");

        //export file example java -jar census.jar -e -f folder_name -ch census_chro.xml
        //This will generate census-out.txt file in the same folder

	//opt.addOption("example", false, "java -jar census.jar -c /my_folder/census_config.xml -f folder_path -i DTASelect-filter.txt");
	//opt.addOption("example_to_generate_report", false, "java -jar census.jar -e census.param -f folder_path -chro census_chro.xml");

        BasicParser parser = null;
        CommandLine cl = null;

        try {
            parser = new BasicParser();
            cl = parser.parse(opt, args);
            if ( cl.hasOption('h') || !cl.iterator().hasNext() ) {
                HelpFormatter f = new HelpFormatter();
                f.printHelp("OptionsTip", opt);
		printHelp();

                return;
            }
        } catch (org.apache.commons.cli.MissingArgumentException me)
        {
            System.out.println(me.getMessage());
            System.out.println();
            System.out.println("See Census help");
            printHelp();

            System.exit(0);
        }
	/*
	else {


	}
*/
    try {
       census = new Census(cl);

       }
      catch (IOException ioe)
	{
            ioe.printStackTrace();
	    System.out.println(ioe.getMessage());
	}
      catch (CensusGeneralException ce)
	{
	    System.out.println(ce.getMessage());
	}
	catch (Exception e)
	{
	    System.out.println(e.getMessage());
	    e.printStackTrace();

	}

       // if(args.length>0)
	 //   census = new Census(args);
        //else

//        System.out.println( System.currentTimeMillis() -start );
    }

    public Census(CommandLine cl) throws IOException, CensusGeneralException, Exception
    {
        init(cl);
    //    run();
    }

    public void init(CommandLine cl) throws IOException, CensusGeneralException, Exception
    {

	conf = Configuration.getInstance();
	conf.setFilePath(filePath);

	conf.setStartTime();
        System.out.println("Census new version " + conf.getVersion());

	//System.out.println(args.length + " " + args[1]);

	String configFile = null;
	boolean isXmlConfig = true; //alway true from this version
//	isXmlConfig = true;

        if( cl.hasOption('v') )
            return;

        if( cl.hasOption("spc") )
	{
	    String[] arr = cl.getOptionValues("spc");
	    SpecCountUtil util = new SpecCountUtil("spec_count_out.txt");

	    for(int i=0;i<arr.length;i++)
	    {
		util.addFile(arr[i]);
	    }

	    util.runSpecCount();

	    System.out.println("spec_count_out.txt was successfully generated");
            return;
	}

        if( cl.hasOption("etan") ) {

	    String chroFile = null;
	    if( cl.hasOption("etan") ) {
		chroFile = cl.getOptionValue("etan");
	    } else {
		System.out.println("etan option needs parameter");
		return;
	    }


	    String outfile = "census-out.txt";
	    if( cl.hasOption("of") )
		    conf.setOutputFilename( cl.getOptionValue("of") );


		String path = "";
	    if( cl.hasOption("f") ) {
		outfile = cl.getOptionValue("f") + File.separator + outfile;
		path = cl.getOptionValue("f");

		conf = Configuration.getInstance();
		conf.readXMLParam(path + File.separator + "census_config.xml");

	    } else {

		outfile = "." + File.separator + outfile;
                path = ".";
                                conf = Configuration.getInstance();
		conf.readXMLParam(path + File.separator + "census_config.xml");

	    }

	    PrintStream p = new PrintStream( new BufferedOutputStream(new FileOutputStream(outfile)));

	    double intensityThreshold = 0;
            int intThresholdType = 0; //0 don't use, 1 statistical approach (z score), 2 hard threshold
	    //String type = "mode";
            String type = "mean";
	    if( cl.hasOption("int") ) {

                String tmpValue = cl.getOptionValue("int");

                if(tmpValue.startsWith("z")) {
                    intThresholdType = 1;
                    intensityThreshold = Double.parseDouble(tmpValue.substring(1));
                } else {
                    intThresholdType = 2;
                    intensityThreshold = Double.parseDouble(tmpValue);
                }
	    }

	    if( cl.hasOption("noise") ) {
		intensityThreshold = -1;
	    }

	    if( cl.hasOption("t") ) {
		type = cl.getOptionValue("t");
	    }

	    ChroXmlReader cr = new ChroXmlReader(chroFile);
//	    if( cl.hasOption("ref") ) {
//		    result = RelExMainFrame.exportITRAQSingleReportCoreRef(p, intensityThreshold, null, "", true, cr.getProteinList(), path);
//	    } else {
//            String outlier="PROTEIN_LEVEL";




			ReportParam rParam = ReportParamReader.readParam(path + File.separator+"census.param");
			rParam.setConf( Configuration.getInstance() );
			rParam.setIsGui(false);

			RelExMainFrame.exportITRAQSingleReportCore(p, intensityThreshold, intThresholdType, null, "", true, cr.getProteinList(), path, type,rParam);

    //p.print(result);

	    if(null != p)
		p.close();

            System.out.println("census-out.txt was successfully generated");

	    String outlierLevel = Configuration.getInstance().getTmtOutlierLevel();
	    if("DMCCLAT".equals(outlierLevel)) {
//		    TmtfilterPeptide tmtfilter= new TmtfilterPeptide();
                TmtfilterDan tmtfilter = new TmtfilterDan();
		    tmtfilter.parseCensusFile(path);
		    System.out.println("census-out-Tmtfilter.txt was successfully generated");
	    }

	    return;
	}

        if( cl.hasOption('e') )
        {

            conf = Configuration.getInstance();

            String path = "";
	    if( cl.hasOption("f") ) path = cl.getOptionValue("f");
            else path = ".";

            conf.readXMLParam(path + File.separator + "census_config.xml");


            String chroFileName = cl.getOptionValue("chro");

            String param = cl.getOptionValue("e");


            if(null == chroFileName || null == param)
            {
                System.out.println("Error.  Usage example export file example: java -jar census.jar -e census_report.params -f folder_name -ch census_chro.xml");
                return;
            }


	    if(null == path)
	    {
		path = ".";
		System.out.println("Working folder is set to the current folder");
	    }

	    //            ReportGenerator.exportITRAQMultipleReport(chroFileName, new File(folder + File.separator + "census-out.txt"));
	    if(!path.endsWith(File.separator))
		path += File.separator;

	    BufferedReader br = new BufferedReader(new FileReader(param));
	    String eachLine=null;

	    String reportFile = "census-out.txt";

	    if( cl.hasOption("of") )
		    reportFile = cl.getOptionValue("of");

	    //File file = new File(path + File.separator + "census-out.txt");
	    File file = new File(path + File.separator + reportFile);

	    File singleFile = new File(path + File.separator + "census-out_singleton.txt");
	    File chroFile = new File(path + File.separator + chroFileName);

	    boolean isSuccessful=true;
            String errorMessage="";
            PrintStream p = new PrintStream( new BufferedOutputStream(new FileOutputStream(file)));
            PrintStream singleP = new PrintStream( new BufferedOutputStream(new FileOutputStream(singleFile)));
	    ReportResult rResult = null;

	    try {

		ChroXmlReader cr = new ChroXmlReader(chroFileName);
	        ArrayList<ChroProtein> proteinList = cr.getProteinList();

		ReportParam rParam = ReportParamReader.readParam(param);
		rParam.setProteinList(proteinList);
		rParam.setConf( Configuration.getInstance() );
		rParam.setIsGui(false);

               // if(((ChroPeptide)rParam.getProteinList().get(0).getPeptideList().get(0)).getIsoDataList() == null)
               // {

                if(null != conf.getQuantType() && conf.getQuantType().equals("15N enrichment")) {
                    rParam.addN15Params(param);
                    rResult = RelExMainFrame.runReportN15(rParam, p, singleP);
                }
                else {
                    //RelExMainFrame.runReport(rParam, p, singleP)
                    if("MultipleMs1Labeling".equals(conf.getQuantType()))
                        rResult = RelExMainFrame.runReportMultipleMs1labeling(rParam, p);
                    else if("DIA_LF".equals(conf.getQuantType())) {
                      rParam.setFilterFragmentIons(true);

                      ChroXmlReaderDIALF crdia = new ChroXmlReaderDIALF(chroFileName);
                      proteinList = crdia.getProteinList();
                      rParam.setProteinList(proteinList);

                      rResult = RelExMainFrame.runReportDIALF(rParam, p, singleP);
                    }
                    else
                        rResult = RelExMainFrame.runReport(rParam, p, singleP);
                }
                //else
//                    rResult = RelExMainFrame.runReportIsotops(rParam, p, singleP);
// added By Harshil Shah
//
//-- till here...
	    }
	    catch(Exception e)
	    {
		isSuccessful=false;
		e.printStackTrace();
	    }

	    int totalCount=rResult.getTotalCount();
	    int quantifiedCount=rResult.getQuantifiedCount();
	    int redunProteinCount=rResult.getRedunProteinCount();
	    int uniqueProteinCount=rResult.getUniqueProteinCount();
	    int proteinGroupCount=rResult.getProteinGroupCount();

	    if(isSuccessful)
	    {
		System.out.println("Report file was successfully created");
		System.out.println("Total peptides : " + totalCount);
		System.out.println("Quantified peptides : " + quantifiedCount);
		System.out.println("Quantified peptides including singletons : " + rResult.getQuantifiedCountWithSingleton());

		System.out.println("Quantification efficiency : " + CensusHelper.format.format( (double)quantifiedCount/totalCount*100 ) + " %");

      //  System.out.println("==" + rResult.getQuantifiedCount() + " "  + rResult.getQuantifiedCountWithSingleton()  + " " + totalCount);

        System.out.println("Quantification efficiency including singletons : " + CensusHelper.format.format( (double)(rResult.getQuantifiedCountWithSingleton())/totalCount*100 ) + " %");
	    }
	    else
	    {
		System.out.println("Failed to generate the report file");
	    }

	    if(null != p)
		    p.close();

	    if(null != singleP)
		singleP.close();

            return;
        }

	if ( cl.hasOption('g') ) {

	    UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
	    RelExMainFrame census = new RelExMainFrame();

	    census.addWindowListener(new WindowAdapter() {
		    public void windowClosing(WindowEvent e) {System.exit(0);}    });

	    census.setVisible(true);
	    census.validate();
	    census.pack();
	    census.setSize(Toolkit.getDefaultToolkit().getScreenSize());

	    //also for the web start
            String chroFileName = cl.getOptionValue("chro");

	    System.out.println("chro name  ==" + chroFileName);
	    if(null != chroFileName && !"".equals(chroFileName)) {

		if(chroFileName.startsWith("http"))
		    census.openChroFile(new java.net.URI(chroFileName));
		else
		    census.openChroFile(chroFileName);
	    }

	    return;
	}

	if ( cl.hasOption('c') ) {
	    configFile = cl.getOptionValue("c");
	}

	if(null == configFile)
	{
	    System.out.println("Error : config file required. Use -c option to specify location of the config file.");
	    System.exit(0);
	}

	int fileIndex = -1;
	if(File.separator.equals("/"))  //unix
	{
		fileIndex = configFile.lastIndexOf("/");
	}else{ //windows
		fileIndex = configFile.lastIndexOf("\\");
	}
	if( fileIndex <0 ) //current folder
			conf.setConfigFilePath(".");
	else
			conf.setConfigFilePath(configFile.substring(0, fileIndex));

	String path = cl.getOptionValue("f");

        String idFileName = cl.getOptionValue("i");

	if(null == path)
	{
	    path = ".";
	    System.out.println("Working folder is set to the current folder");
	}

	if(!path.endsWith(File.separator))
		path += File.separator;

	if( cl.hasOption('x') ) {
	    String mzXMLFile = cl.getOptionValue("x");
	    conf.setSpectrumFormat( Configuration.MZXML_FILE_FORMAT );
//	    conf.setMzXMLFilePath( mzXMLFile );
	}

        if( cl.hasOption("log") ) {
	    conf.setPrintLog(true);
	}

        if( cl.hasOption('p') ) {
	    conf.setUseProline(true);
	}

/*
	//file could be generated from either linux or window
	if( !path.endsWith("/") && !path.endsWith("\\") ) {
	    //for linux
	    if(path.startsWith("/"))
		path += "/";
	    else //for window
		path += "\\";
	}
*/
	conf.setFilePath(path);

	SAXBuilder sb = new SAXBuilder();

	Document doc = sb.build(configFile);
	final Element root = doc.getRootElement();
	boolean isLabeled = "true".equals( root.getChild("label_type").getAttributeValue("labeling"))?true:false;

	Element expTypeEle = root.getChild("experiment_type");
	if(null != expTypeEle)
	{
	    String expType = expTypeEle.getText();
	    if(null != expType && !"".equals(expType))
	    {
		conf.setExpType( Integer.parseInt(expType) );

	    }

	}

	if(isXmlConfig)
	    conf.readXMLParam(configFile);
	else
	    conf.readParam(path, "census.param");
	//conf.readParam(path, "census_config.txt");
         String label = conf.getMs2Label();

	//Generate MS1 files
	File ff = new File(path);
        String[] flist = ff.list(new RelExFileFilter(CensusConstants.MS1_FILE));
        String[] ms2flist = ff.list(new RelExFileFilter(CensusConstants.MS2_FILE));
        String[] sqlite_list = ff.list(new RelExFileFilter(CensusConstants.SQLITE_FILE));
	String[] mzflist = ff.list(new RelExFileFilter(CensusConstants.MZXML));

	/***************  MRM EXPERIMENT WITHOUT IDENTIFICATION ***********************/
	if( conf.getExpType() == Configuration.MRM_WITHOUT_ID )
	{
	    ChroGenerator chro = new ChroGenerator();
            OpenConfigDialog.runMRMWithoutId(chro);

	    return;
	}

	if(null != flist && isLabeled && flist.length<=0 && conf.getQuantLevel()==1) // && mzflist.length<=0)
	{
//	    System.out.println("Error: Spectral files are not found.");
	    //throw new CensusGeneralException("Error: Spectral files are ot found. If you use mzXML, please use option '-x'");
	    System.out.println("Start converting mzXML to MS1 files...  Census does this only once.");
	    Mzxml2Ms.converMzXML2MS(path);
	    System.out.println("Converting mzXML to MS1 files is completed.");
	}

	if( ms2flist.length<=0 && sqlite_list.length<=0 && conf.getQuantLevel() == 2)
	{
	    System.out.println("Start converting mzXML to MS2 files... Census does this only once.");
	    Mzxml2Ms.converMzXML2MS2(path);
	    System.out.println("Converting mzXML to MS2 files is completed.");
	}


	/***************  END OF MRM EXPERIMENT WITHOUT IDENTIFICATION ****************/

	//conf.setDtaSelectFile(path + "DTASelect-filter.txt");
        conf.setIdFileName(idFileName);


	//conf.setLabeling(true);


	//final File file = choose.getSelectedFile();

	//if(isXmlConfig)
	//{

	//if( "true".equals( root.getChild("label_type").getAttributeValue("labeling") ) )
	if( isLabeled || (conf.getQuantType() != null && conf.getQuantType().equals("DIA_LF")) )
	{
	    ChroGenerator chro = new ChroGenerator(
		    //null, null //
		    null//.getProgressBar(),
		    //chroProgress.getProgressText()
		    //configFile,
		    //dtaFile
		    //massTolerance
		    );

            OpenConfigDialog.runLabeledAnalysis(chro, conf);

	}
	else
	{
	    if(null == conf.getFilePath() || "".equals(conf.getFilePath()))
	    {
		int fIndex = -1;

		if(File.separator.equals("/"))  //unix
		{
		    fIndex = configFile.lastIndexOf("/");
		}
		else //windows
		{
		    fIndex = configFile.lastIndexOf("\\");
		}

		if( fIndex <0 ) //current folder
		    conf.setFilePath(".");
		else
		    conf.setFilePath(configFile.substring(0, fIndex));
	    }

	    List sampleEleList = root.getChildren("sample");

	    //Configuration.
	    String refFileName = "";
	    Element refEle = root.getChild("ref");
	    if(null != refEle)
	    {
		//String refSamName = refEle.getChildText("sample_name");
		refFileName = refEle.getChildText("file_name");

		conf.setRefFileName(refFileName);
	    }


	    int index=0;
	    int refIndex=0;

	    Vector<String> fileNameList = new Vector<String>();
	    Vector<String> sampleNameList = new Vector<String>();
	    Vector<String> pathFileNameList = new Vector<String>();

	    //confSam.setRefFileName(refFileName);
	    //Hashtable<String, Hashtable> masterHt = new Hashtable<String, Hashtable>();
	    HashSet set = new HashSet();

	    for(Iterator<Element> itr = sampleEleList.iterator(); itr.hasNext(); )
	    {
		Element sam = itr.next();
		//String samName = sam.getChildText("name");
                String samName = sam.getAttributeValue("group");
		//populate configuration class
		Configuration.Sample confSam = new Configuration.Sample();
		confSam.setName(samName);

			List<String> filesList  = new ArrayList<String>();
			List<String> pathList  = new ArrayList<String>();

                List<Element> sampleList = sam.getChildren("each_sample");
                for(Iterator<Element> eachSampleItr=sampleList.iterator(); eachSampleItr.hasNext(); ) {
                    Element eachSample = eachSampleItr.next();
                    List<Element> fileList = eachSample.getChild("ms_files").getChildren("file");



					String eachPath = fileList.get(0).getText();
					eachPath = eachPath.substring(0, eachPath.lastIndexOf(File.separator));
					pathList.add(eachPath);
                    SampleModel sampleModel = new SampleModel(samName);

                    if(fileList.size()>0)
                    {
        String firstFileName = fileList.get(0).getText();
                        String filePath = firstFileName.substring(0, firstFileName.lastIndexOf(File.separator));

                        firstFileName = firstFileName.substring(firstFileName.lastIndexOf(File.separator)+1);

                        if(firstFileName.startsWith("*"))
                        {
                            String extension = firstFileName.substring(firstFileName.lastIndexOf(".")+1);

                    //        edu.scripps.pms.census.util.RelExFileFilter fFilter = new edu.scripps.pms.census.util.RelExFileFilter(extension);

                            File specFile = new File(filePath);
                            String[] splist = specFile.list(new edu.scripps.pms.census.util.RelExFileFilter(extension));

                            for(String eachFileName : splist)
                            {
                                eachFileName = filePath + File.separator + eachFileName;
                                set.add(filePath);
								filesList.add(filePath);
                                sampleModel.addPath(filePath);

                                confSam.addFile(eachFileName);
                                pathFileNameList.add(eachFileName);
                                fileNameList.add(eachFileName);
                                sampleNameList.add(samName);
                            }
                        }
                        else {
                            for(Iterator<Element> itr1 = fileList.iterator(); itr1.hasNext(); )
                            {
                                Element eachFile = itr1.next();
                                String fileName = eachFile.getText();

                                confSam.addFile(fileName);
                                pathFileNameList.add(fileName);

                                set.add(fileName.substring(0, fileName.lastIndexOf(File.separator)));
								filesList.add(fileName.substring(0, fileName.lastIndexOf(File.separator)));
                                sampleModel.addPath(fileName.substring(0, fileName.lastIndexOf(File.separator)));

                                if(fileName.endsWith("ms2"))
                                {
                                    fileName = fileName.substring(0, fileName.length()-3);
                                    fileName += "ms1";
                                }

                                fileNameList.add(fileName);
                                sampleNameList.add(samName);

                                if( null != refEle && refFileName.equals(fileName) )
                                    refIndex = index;

                                index++;
                            }
                        }

                        conf.addExp(confSam);
                        conf.addSample(sampleModel);
                    //List<Element> fileList = sam.getChild("ms_files").getChildren("file");
                    }
		}
			conf.getNonlabelFilenameGroupMap().put(samName, pathList);
	    }

	    conf.setNonlabelFilenameList(pathFileNameList);
	    conf.setNonlabelFilePaths(set);

	    long start = System.currentTimeMillis();

	    if(cl.hasOption("lf"))//labelfree Run
        {
            String jsonFile = cl.getOptionValue("lf");

            ChroJSONReader cr = new ChroJSONReader(configFile);
           // cr.parse(configFile, jsonFile);
             cr.parse(jsonFile);
        } else
            {
	    if( conf.isBasedOnId() )
	    {

			if( cl.hasOption("of") )
				conf.setOutputFilename( cl.getOptionValue("of") );

			runLabelFreeBasedOnId(root, cl.hasOption("aa"), cl.hasOption("d"));
			//	runLabelFreeBasedOnId(root, cl.hasOption("aa"), false);

			return;
	    } else if(conf.isTargeted()) {

			System.out.println("targeted...");
			String peptideFile = configFile.substring(0, configFile.lastIndexOf(File.separator)) + File.separator + "peptides.txt";
			if(!new File(peptideFile).exists()) {
				System.out.println("peptides.txt input file is required in same config file folder");
				return;
			}

			//LabelfreeTargeted.xxx
			return;
		}


	    chroalign align = new chroalign();

	    int[][][] pathArray = null;

	    boolean isAlign = conf.isAlign();

	    //********  check if alignment was done already ***//
	    File alignFile = new File(conf.getFilePath() + "aligned_out.xml");
	    File pathoutFile = new File(conf.getFilePath() + "path_out.xml");
	    File chrooutFile = new File(conf.getFilePath() + "chro_out.xml");

	    int selection = -1;


	    if(alignFile.exists() && pathoutFile.exists() && chrooutFile.exists() )
	    {
		System.out.print("Alignment out files are found.  Do you want to read them? (y|n) : ");

		while(true)
		{

		    try {
			BufferedReader ibr = new BufferedReader(new InputStreamReader(System.in));
			String input = ibr.readLine();

			if("y".equals(input))
			{
			    selection = 0;
			    break;
			}
			else if("n".equals(input))
			{
			    selection = -1;
			    break;
			}

			System.out.print("y or n?");
			//    userName = br.readLine();
		    } catch (IOException ioe) {
			System.out.println("IO error trying to read your input!");
			System.exit(1);
		    }

		}
	    }

	    if(selection != 0) //new align
	    {
		if(isAlign) //align based on chromatogram profile
		{
		    pathArray = align.alignChro(null, fileNameList.toArray(), sampleNameList.toArray(), refIndex, false, 500, conf.getFilePath());

		    NonLabelMappingModel mapModel = new NonLabelMappingModel(pathArray, pathFileNameList, refIndex);
		    conf.setMapModel(mapModel);
		}
		else //align based on retention time
		{
		    pathArray = align.noAlignChro(null, fileNameList.toArray(), sampleNameList.toArray(), conf.getFilePath());

		    NonLabelMappingModel mapModel = new NonLabelMappingModel(pathArray, pathFileNameList, refIndex);
		    mapModel.reinitializeMaxIndexByRet(); //this method is called only for the aligning by retention time.  Not nice idea.
		    conf.setMapModel(mapModel);
		}

	    }
	    else //just read aligned data
	    {
		System.out.println("Reading pre-aligned xml file...");

		SAXBuilder builder = new SAXBuilder();
		Document tdoc = builder.build( pathoutFile );
		Element rootEle = tdoc.getRootElement();

		int numberOfAlignment = 1;  //number of alignment is one bigger than actual alignment number :-(
		int biggestNum=0;

		for(Iterator<Element> itr=rootEle.getChildren("dataset").iterator(); itr.hasNext(); )
		{
		    Element ele = itr.next();

		    if(biggestNum< ele.getChildren().size())
			biggestNum = ele.getChildren().size();
		    numberOfAlignment++;
		}

		pathArray = new int[numberOfAlignment][2][biggestNum+1];

		int firstIndex=0;

		for(Iterator<Element> itr=rootEle.getChildren("dataset").iterator(); itr.hasNext(); )
		{
		    Element ele = itr.next();

		    int thirdIndex=0;
		    for(Iterator<Element> pitr = ele.getChildren("p").iterator(); pitr.hasNext(); )
		    {
			Element pele = pitr.next();

			pathArray[firstIndex][0][thirdIndex] = Integer.parseInt(pele.getAttributeValue("x"));
			pathArray[firstIndex][1][thirdIndex] = Integer.parseInt(pele.getAttributeValue("y"));
			//System.out.println( pathArray[firstIndex][0][thirdIndex]  + "\t" + pathArray[firstIndex][1][thirdIndex] );
			thirdIndex++;
		    }

		    firstIndex++;
		}

		System.out.println("Populating mapping model...");

		NonLabelMappingModel mapModel = new NonLabelMappingModel(pathArray, pathFileNameList, refIndex);
		conf.setMapModel(mapModel);

		System.out.println("Pre-alignment file loading done.");


	    }

	    String[] targetMS1Files = align.getTargetMS1Files();
	    String referenceMS1File = align.getReferenceMS1File();

	    conf.setNonlabelFilePaths(set);

	    //chroProgress.addMessage("Reading Spectra and Calculating Quantification...");

	    ChroGenerator chro = new ChroGenerator(
		    //aJProgressBar, //commented out for testing
		    null, //.getProgressBar(),
		    //null,
		    null,
		    referenceMS1File,
		    targetMS1Files,
		    root,
		    pathArray
		    //isotopeFileField.getText().trim(),
		    );


	    //if( 1 == conf.getQuantLevel() )
	    chro.createNonlabelXmlChro();

	}
        }

        //elementFile = filePath + File.separator + conf.getElementCompFile();

    }

    private void runLabelFreeBasedOnId(Element root, boolean autoAnswer, boolean regenerate) throws IOException, Exception
    {
	ChroGenerator chro = new ChroGenerator(
		//aJProgressBar, //commented out for testing
		null, //.getProgressBar(),
		//null,
		null,
//		referenceMS1File,
		null,
//		targetMS1Files,
		null,
		root,
//		pathArray
		null
		//isotopeFileField.getText().trim(),
		);

	chro.createLabelFreeBasedOnDirectId(autoAnswer, regenerate);
    }

    private void printMode(boolean isDataIndependent, boolean isHighRes)
    {
	System.out.println("\n****************************************************************");
	System.out.println("*                                                              *");
	    System.out.println("*  Running Census...                                           *");

	if(isDataIndependent)
	    System.out.println("*  in Data independent mode                                    *");
	else
	    System.out.println("*  in Data dependent mode                                      *");

	if( conf.isHighRes() )
	    System.out.println("*  in High Resolution mode                                     *");
	else
	    System.out.println("*  in Low Resolution mode                                      *");

	System.out.println("*                                                              *");
	System.out.println("****************************************************************");

    }
}
