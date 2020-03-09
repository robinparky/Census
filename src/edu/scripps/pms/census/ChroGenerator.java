/*
 * ChroGenerator.java
 *
 * Created on March 21, 2005, 10:22 AM
 */
package edu.scripps.pms.census;

import com.google.gson.Gson;
import edu.scripps.pms.census.tandem.MZValues;
import edu.scripps.pms.census.util.SpectrumUtil;
import edu.scripps.pms.census.util.XYPoint;
import edu.scripps.pms.census.util.io.*;

import java.sql.SQLException;
import java.util.*;
import java.text.DecimalFormat;
import java.io.*;
import javax.swing.JProgressBar;
import javax.swing.JTextArea;

import edu.scripps.pms.mspid.Modification;
import edu.scripps.pms.util.FileFilterUtil;
import edu.scripps.pms.util.spectrum.*;

import edu.scripps.pms.census.util.dtaselect.Protein;
import edu.scripps.pms.census.util.dtaselect.Peptide;
import edu.scripps.pms.census.model.IsotopeTable;

import edu.scripps.pms.census.io.*;
import edu.scripps.pms.census.io.SpecRangeGenerator;
import edu.scripps.pms.census.model.SpecRange;
import edu.scripps.pms.census.tools.Formatter;

import edu.scripps.pms.census.util.*;

import edu.scripps.pms.census.util.RelExFileFilter;
import edu.scripps.pms.census.hash.*;
import edu.scripps.pms.census.conf.Configuration;
import edu.scripps.pms.census.exception.*;
import edu.scripps.pms.census.labelFree.util.LabelfreeChroUtil;

import edu.scripps.pms.util.sqlite.spectra.CreateDb;
import edu.scripps.pms.util.sqlite.spectra.SpectraDB;
import gnu.trove.*;

import org.apache.commons.lang3.StringUtils;
import org.jdom.*;
import org.jdom.output.*;

import edu.scripps.pms.census.model.*;
import edu.scripps.pms.census.tmtFilter.TMTUtil;
import org.apache.commons.lang3.ArrayUtils;
import scripts.MSSplitFolderCreation;
import scripts.histogram.Histogram;
import scripts.mrm.*;

/**
 *
 * @author rpark
 * @author rohan
 * @version $Id: ChroGenerator.java,v 1.73 2014/08/27 18:00:34 rpark Exp $
 */
public class ChroGenerator {

    public final static String MASTER_JSON = "JSON_OBJ";
    public final static String XY_VALUES = "xyvalues";
    private int scanBefore;
    private int scanAfter;
    private String isotopeFile;
    public final static String PURITY_GRAPH_PATH = "purity.json";
    public final static String SIGNAL_GRAPH_PATH="signal_to_noise.json";
    //private String dtaSelectFilterFile = "DTASelect-filter.txt";
    //private String dtaSelectFile = "DTASelect.txt";

    private String filePath = null;
    //private JProgressBar bar;
    //private JTextArea text;
    //private String newline = "\n";
    public static final String CENSUS_CONFIG = "census_config.xml";
    //private final String INDEX_FILE = "index";
    public static final double PROTON_MASS = 1.00728;

    private double massTolerance;
    private static Configuration conf;
    private DecimalFormat formatter = new DecimalFormat("0.000");

    private Element confRootEle = null;

    private String refMS1File;
    private String[] targetMS1Files;
    private int[][][] pathArray;

    private static ChroProgressDialog progress;

    private long startTime;

    public ChroGenerator() {
        startTime = System.currentTimeMillis();
    }

    /**
     * Creates a new instance of ChroGenerator. To be deprecated.
     */
    public ChroGenerator(
            JProgressBar bar,
            JTextArea text,
            int scanBefore,
            int scanAfter,
            String isotopeFile,
            double massTolerance
    ) throws IOException {
        startTime = System.currentTimeMillis();
        //this.text = text;
        this.scanBefore = scanBefore;
        this.scanAfter = scanAfter;
        this.isotopeFile = isotopeFile;
        this.filePath = isotopeFile.substring(0, isotopeFile.lastIndexOf(File.separator) + 1);
        //this.dtaSelectFilterFile = filePath + dtaSelectFilterFile;
        this.massTolerance = massTolerance;

        conf = Configuration.getInstance();

    }

    /**
     * new ChroGenerator constractor
     */
    public ChroGenerator(
            ChroProgressDialog progress
    //JTextArea text
    ) throws IOException {
        //this.bar = bar;
        this.progress = progress;
        startTime = System.currentTimeMillis();

        // this.text = text;
//        this.isotopeFile = isotopeFile;
//        this.filePath = isotopeFile.substring(0, isotopeFile.lastIndexOf(File.separator) + 1);
//        this.dtaSelectFilterFile = filePath + dtaSelectFilterFile;
//        this.massTolerance = massTolerance;
        conf = Configuration.getInstance();
        this.filePath = conf.getFilePath();

        if (conf.isXmlConf()) {
            this.confRootEle = conf.getRootConfEle();
        }

        this.massTolerance = conf.getMassTolerance();

    }

    //* non-labeling quantification */
    public ChroGenerator(
            ChroProgressDialog progress,
            JTextArea text,
            String refMS1File,
            String[] targetMS1Files,
            Element confRootEle,
            int[][][] pathArray
    ) {
        this.progress = progress;
        //this.text = text;

        this.refMS1File = refMS1File;
        this.targetMS1Files = targetMS1Files;
        this.pathArray = pathArray;

        //this.isotopeFile = isotopeFile;
        //if(!this.filePath.endsWith("/"))
        //    this.filePath += File.separator;
        //this.dtaSelectFilterFile = dtaFile;
        //SAXBuilder builder = new SAXBuilder();
        //Document doc = builder.build(new File(configFile));
        this.confRootEle = confRootEle; //doc.getRootElement();

        //this.dtaSelectFilterFile = filePath + dtaSelectFilterFile;
        //this.massTolerance = massTolerance;
        conf = Configuration.getInstance();
        this.filePath = conf.getFilePath();

        startTime = System.currentTimeMillis();

    }

    //non labeling //label free based on direct id comparison new approach
    public void createLabelFreeBasedOnDirectId(boolean autoAnswer, boolean readExisting) throws IOException, Exception {
        conf = Configuration.getInstance();
        //List<String> pathFileNameList = conf.getNonlabelFilenameList();
        Map<String, List<String>> nonlabelFilenameGroupMap = conf.getNonlabelFilenameGroupMap();
        Map<String, Integer> groupCompletionProcess = new HashMap<String, Integer>();
        Set<String> groupNameSet = nonlabelFilenameGroupMap.keySet();
        //List tmpFileNameList = conf.getNonlabelFilenameList();

        String chroTmpFileName = "census_chro_temp.xml";

        int fileCount = 0;
        boolean isReadAll = false;
        if(readExisting) isReadAll=true;

        if (autoAnswer) {
            isReadAll = true;
        }

        for (String groupName : groupNameSet) {
            groupCompletionProcess.put(groupName, 0);
        }
        String processFile = conf.getConfigFilePath() + File.separator + "progress.properties";
        CensusChroProgress censusChroProgress = new CensusChroProgress(processFile, groupCompletionProcess);
        for (String groupName : groupNameSet) {
            List<String> pathList = nonlabelFilenameGroupMap.get(groupName);
            if (pathList.size() > 0) {
                int eachFileProcess = Math.abs(100 / pathList.size());
                for (int pathCounter = 0; pathCounter < pathList.size(); pathCounter++) {
                    String eachPath = pathList.get(pathCounter);
                    if (!eachPath.endsWith(File.separator)) {
                        eachPath += File.separator;
                    }

                    File tmpF = new File(eachPath + chroTmpFileName);

                    boolean isRead = false;

                    if (readExisting) {
                        if (tmpF.exists() && !isReadAll) {

                            if (autoAnswer) {
                                isRead = true;
                            } else {
                                System.out.print(eachPath + chroTmpFileName + " are found.  Do you want to read them? (y|n|ya(yes all)) : ");
                                while (true) {
                                    try {

                                        BufferedReader ibr = new BufferedReader(new InputStreamReader(System.in));
                                        String input = ibr.readLine();


                                        if ("y".equals(input)) {
                                            isRead = true;
                                            break;
                                        } else if ("n".equals(input)) {
                                            isRead = false;
                                            break;
                                        } else if ("ya".equals(input)) {
                                            isRead = true;
                                            isReadAll = true;
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
                        }

                        if (isRead || (isReadAll && tmpF.exists())) {
                            continue;
                        }
                    }

                    String eachIdPath = eachPath + CensusConstants.SEARCH_OUTPUT;
                    String spectraPath = eachPath + "../../spectra/";
                    String splitSpectraPath = eachPath + "../../spectra/split/";
                    //System.out.println("====================\t" + spectraPath);
                    MSSplitFolderCreation msp = new MSSplitFolderCreation();

                    // if version2 file exists in split folder, don't index again.
                  //  HashMap<String, String> splitSpectraMap = msp.splitMS1Files(spectraPath, CensusConstants.LABELFREE_MS1_SPLIT_SCAN_NUM, true);  //if exist, don't index again
//
                    //String splitSpectraPath = spectraPath+ File.separator + "split";
                 //   Hashtable<String, IndexedFile> splitMs1FileHt = null;
                    Hashtable<String, IndexedFile> origMs1FileHt = null;
                    HashMap<String, HashMap<Integer, Integer>> ms2ToMs1Map = null;

                    switch (conf.getSpectrumFormat()) {

                        case Configuration.MS_FILE_FORMAT:

                            origMs1FileHt = createIndexedFiles(spectraPath, CensusConstants.MS1_FILE);
                          //  splitMs1FileHt = createIndexedFiles(splitSpectraPath, CensusConstants.MS1_FILE);
                            ms2ToMs1Map = IndexUtil.buildMS2toMS1ScanMapFiles(spectraPath);
                            break;

                        case Configuration.MZXML_FILE_FORMAT:
                        //    splitMs1FileHt = createIndexedFiles(splitSpectraPath, CensusConstants.MZXML);
                            break;

                        default:
                            break;
                    }

           //         int[] keys;

                    IndexedFile iFile = null;
                    BufferedOutputStream out = null;
                    PrintStream p = null;

                    try {
                        /**
                         * *****************************************************************
                         * Read DTASelect.txt file to find spectrum range for
                         * each peptide
						 *****************************************************************
                         */
                        //TIntLongHashMap index;
                        System.out.print("Parsing " + eachIdPath + "...");
                        conf.setIdFileName(eachIdPath);

                        if (null != progress) {
                            progress.addMessage("\nParsing Identified Peptides...");
                        }

                        if (null != progress) {
                            progress.addMessage("\ndone.");
                        }

                        IsotopeReader isoReader = null;

                        if (null != isotopeFile) {
                            isoReader = new IsotopeReader(isotopeFile);
                        } else {
                            isoReader = new IsotopeReader(this.confRootEle);
                        }

                        IdentificationReader idReader = BaseIdentificationReader.getIdentificationInst(isoReader);
                        System.out.println("done.");

                        SpecRangeGenerator rangeGen = null;

//		File dtaFile = new File(filePath + "DTASelect.txt");

                        /*
		System.out.print("Parsing DTASelect.txt...");
		if(null != progress)
		    progress.addMessage("\nParsing DTASelect.txt...");

		DTASelectForLabelfreeParser generator = new DTASelectForLabelfreeParser(eachPath + "DTASelect.txt", 0.05);
		Hashtable<String, Hashtable> dtaHt = generator.getUnfilteredPeptides();

		System.out.println("\nParsing DTASelect.txt complete.");
                         */
                        rangeGen = new SpecRangeGenerator();

                        Element rootEle = this.createXmlChroHeader(1);
                        ElementComposition totalElement;

                        Element proteinEle = null;

                        //	Element filteredEle = new Element("filtered_peptides");
                        Iterator<Protein> pitr = idReader.getProteins(); //need to run to calculate redundnat peptides
                        int redundantPeptideNum = idReader.getTotalPeptideNumber();
                        double percent = 0.0;
                        double eachSeg = 100.0 / redundantPeptideNum;
                        int pepCount = 0;

                        List<Element> proteinEleList = new ArrayList<Element>();

                        for (Iterator<Protein> itr = pitr; itr.hasNext();) {
                            Protein protein = itr.next();

                            proteinEle = new Element("protein");
                            proteinEle.setAttribute("locus", protein.getLocus());
                            proteinEle.setAttribute("seq_ct", protein.getSeqCount());
                            proteinEle.setAttribute("spec_ct", protein.getSpectrumCount());
                            proteinEle.setAttribute("seq_cov", protein.getSeqCoverage());
                            proteinEle.setAttribute("length", protein.getLength());
                            proteinEle.setAttribute("molwt", protein.getMolWt());
                            proteinEle.setAttribute("pi", protein.getPI());
                            proteinEle.setAttribute("val", protein.getValidation());

                            try {
                                proteinEle.setAttribute("desc", protein.getDescription());
                            } catch (org.jdom.IllegalDataException ide) {
                                proteinEle.setAttribute("desc", StringUtil.removeIsoControlChar(protein.getDescription()));
                            }

//		    List<Peptide> additionPepList = new ArrayList<Peptide>();
                            for (Iterator<Peptide> pepItr = protein.getPeptides(); pepItr.hasNext();) {
                                Peptide peptide = pepItr.next();
                                pepCount++;
                                String pepSequence = peptide.getSequence();
                                String pepKey = pepSequence.substring(2, pepSequence.length() - 2) + peptide.getChargeState();
//			System.out.println(pepSequence.substring(2, pepSequence.length()-2) + peptide.getChargeState() + " " + peptide.getFileName());
                                //Hashtable<String,String> pepFileHt = dtaHt.get(pepKey);

                                /*
			for(Iterator<String> pepFileItr = pepFileHt.keySet().iterator(); pepFileItr.hasNext(); ) {

			    String pepFileKey = pepFileItr.next();

			    if(peptide.getFileName().equals(pepFileKey))
				continue;

			    String tmpPepLine = pepFileHt.get(pepFileKey);
			    Peptide tmpPep = new Peptide();
			    tmpPep.setDTASelectTxtPeptideLine(tmpPepLine);
			    additionPepList.add(tmpPep);

			}
                                 */
                                percent += eachSeg;

                                System.out.print(pepCount);
                                System.out.print("/");
                                System.out.print(redundantPeptideNum);
                                System.out.print(" peptides, ");
                                System.out.print((int) percent);
                                System.out.print(" % is complete\r");

                                if (null != progress) {
                                    progress.setProgress((int) percent);
                                }

                                try {
                                   // Element peptideEle = LabelfreeChroUtil.getPeptideDomElement(peptide, isoReader, spectraPath, origMs1FileHt, splitMs1FileHt, splitSpectraMap, ms2ToMs1Map);
                                      Element peptideEle = LabelfreeChroUtil.getPeptideDomElement(peptide, isoReader, spectraPath, origMs1FileHt,ms2ToMs1Map);
                                    if (null != peptideEle) {
                                        proteinEle.addContent(peptideEle);
                                    }
                                } catch (IOException e) {

                                    e.printStackTrace();

                                    System.out.println("not quantifiable\t" + peptide.getSequence());
                                    continue;
                                } catch (Exception e) {
                                    System.out.println("not quantifiable\t" + peptide.getSequence());

                                    e.printStackTrace();

//                                    System.out.println("uncccccccccccccccomment below");
  //                                  System.exit(0);
                                    continue;
                                }

                            }

//		    System.out.println("------" + protein.getLocus() + " " + proteinEle.getChildren().size());
                            if (proteinEle.getChildren().size() > 0) {

                                for (Iterator<Element> proEleItr = proteinEleList.iterator(); proEleItr.hasNext();) {
                                    Element protempEle = proEleItr.next();

                                    List<Element> pepchildList = proteinEle.getChildren();

                                    for (Iterator<Element> pepItr = pepchildList.iterator(); pepItr.hasNext();) {
                                        Element pepcopy = pepItr.next();
                                        protempEle.addContent((Element) pepcopy.clone());
                                    }

                                    //protempEle.addContent(proteinEle.getChildren().);
                                    rootEle.addContent(protempEle);

                                }

                                rootEle.addContent(proteinEle);
                                proteinEleList.clear();

                            } else {
                                proteinEleList.add(proteinEle);
                            }

                        }

                        Document doc = new Document(rootEle);
                        OutputStream os = new FileOutputStream(eachPath + chroTmpFileName);
                        XMLOutputter outputter = new XMLOutputter();
                        outputter.setFormat(Format.getPrettyFormat());
                        outputter.output(doc, os);
                        os.close();

                        if (pathCounter == (pathList.size() - 1)) {
                            groupCompletionProcess.put(groupName, 100);
                            censusChroProgress.updateProcess(groupName, 100);
                        } else {
                            Integer progress = groupCompletionProcess.get(groupName) + eachFileProcess;
                            groupCompletionProcess.put(groupName, progress);
                            censusChroProgress.updateProcess(groupName, progress);
                        }
                        System.out.println(groupCompletionProcess.get(groupName) + " % completed");
                    } catch (IOException e) {
                        System.out.println("IO Error while generating msms chro file : " + e);
                        e.printStackTrace();
                        throw new IOException(e.toString());
                    } catch (Exception e) {
                        System.out.println("Error while generating msms chro file : " + e);
                        e.printStackTrace();
                        throw new Exception(e.toString());
                    } finally {
                        if (null != p) {
                            p.close();
                        }

                        if (null != out) {
                            out.close();
                        }

                        //Close all random files
                      /*  for (Enumeration e = splitMs1FileHt.keys(); e.hasMoreElements();) {
                            iFile = splitMs1FileHt.get(e.nextElement());

                            if (null != iFile) {
                                iFile.close();
                            }
                        }*/
                        for(Map.Entry<String,IndexedFile> entry : origMs1FileHt.entrySet())
                        {
                            IndexedFile indexedFile = entry.getValue();
                            if(indexedFile!=null)
                            {
                                SpectraDB db = indexedFile.getSpectraDB(false);
                                if(db!=null)
                                {
                                    db.close();
                                }
                            }
                        }

                    }
                    System.out.println((System.currentTimeMillis() - startTime) * 0.001 + " seconds taken");
                }
            }
        }

        MergeLabelFreeChro.mergeLabelFreeChro(conf.getSampleList(), chroTmpFileName);

    }

    private Element getPeptideDomElement(Peptide peptide, IsotopeReader isoReader, String eachPath, Hashtable<String, IndexedFile> ht) throws IOException, Exception {
        String pepSequence = peptide.getSequence();
        String fileName = peptide.getFileName();
        fileName = cleanFileName(fileName);

        //trim additional characters from peptide sequence at both ends
        char[] ch = pepSequence.substring(2, peptide.getSequence().length() - 2).toCharArray();

        ElementComposition element = new ElementComposition(ch, 0, ch.length, isoReader.getIsotope());

        try {
            element.calculate();
        } catch (InvalidAAException invE) {
            System.out.println("Not Quantifiable peptide : " + pepSequence);
            return null;
        }

        if (!element.isQuantifiable()) {
            System.out.print("\nError : ");
            System.out.println(pepSequence + " is not quantifiable.");
            return null;
        }

        Configuration conf = Configuration.getInstance();

        IsotopeDist sampleDist = null;
        IsotopeDist refDist = null;
        if ("15N".equals(conf.getQuantType())) {
            sampleDist = new IsotopeDist15N(element.getElementSampleArr(), element.getModShift(), true);
            refDist = new IsotopeDist15N(element.getElementRefArr(), element.getModShift(), false);
        } else {
            sampleDist = new IsotopeDist(element.getElementSampleArr(), element.getModShift(), true);
            refDist = new IsotopeDist(element.getElementRefArr(), element.getModShift(), false);
        }

        IndexedFile iFile = ht.get(eachPath + fileName + ".ms1");

        TIntDoubleHashMap retentionTimeMap = iFile.getRetentionTimeMap();
        TIntDoubleHashMap ionInjectionMap = iFile.getIonInjectionMap();

        int scanNumber = peptide.getScanNumber();
        double retTime = retentionTimeMap.get(scanNumber);
        double ionInjectionTime = ionInjectionMap.get(scanNumber);

        int tmpScanNumber = scanNumber;
        int tmpCount = 0;
        while (retTime <= 0) {
            retTime = retentionTimeMap.get(--tmpScanNumber);
            tmpCount++;
            if (tmpCount > 200) { //check if ms1.index file contains retention time
                System.out.println("==retention time is required in ms1.index file.");
                System.exit(0);
            }
        }

        ionInjectionTime = ionInjectionMap.get(tmpScanNumber);

        // System.out.println( retentionTimeMap.get( peptide.getScanNumber()) );
        //  System.out.println( retTime + " == " + tmpScanNumber );
        Element peptideEle = this.createXmlChroPeptideTitle(true, peptide); //true is for full scan
        peptideEle.setAttribute("rt", String.valueOf(retTime));
        peptideEle.setAttribute("iit", String.valueOf(ionInjectionTime));
        conf.setCalcSamAvgMass(sampleDist.getAvgMass());
        conf.setCalcRefAvgMass(refDist.getAvgMass());
        peptideEle.setAttribute("lightStartMass", String.valueOf(sampleDist.getStartMass()));
        peptideEle.setAttribute("heavyStartMass", String.valueOf(refDist.getStartMass()));
        peptideEle.setAttribute("lightAvgMass", String.valueOf(conf.getCalcSamAvgMass()));
        peptideEle.setAttribute("heavyAvgMass", String.valueOf(conf.getCalcRefAvgMass()));

        //String fileName = peptide.getFileName().substring(0, peptide.getFileName().indexOf(".")+1) + "ms1";
        SpecRange range = null;

        int tmpScanNum = Integer.parseInt(peptide.getScanNum());
        range = new SpecRange(tmpScanNum, tmpScanNum);

        //peptideEle.setAttribute("start_scan", peptide.getScanNum());
        //peptideEle.setAttribute("end_scan", peptide.getScanNum());
        int scanNum = Integer.parseInt(peptide.getScanNum());

        switch (conf.getSpectrumFormat()) {
            case Configuration.MS_FILE_FORMAT:
                fileName += ".ms1";
                break;

            case Configuration.MZXML_FILE_FORMAT:
                fileName += ".mzXML";
                break;

            default:
                break;
        }

        //iFile = ht.get(fileName);
        //
        if (null == iFile) {
            iFile = ht.get(eachPath + fileName.substring(1, fileName.length()));
            if (null == iFile) {
                System.out.println("Error : cannot find the file " + eachPath + fileName.substring(1, fileName.length()));
                throw new IOException("Error : cannot find the file " + eachPath + fileName);
            }
            //System.exit(0);
        }

        int[] keys = iFile.getKeys();
        int keyIndex = -1;
        keyIndex = Arrays.binarySearch(keys, scanNum);

        if (keyIndex < 0) //Cannot find index
        {
            keyIndex = -(++keyIndex); //Math.abs(++keyIndex);
        }
        if (keyIndex >= keys.length) {
            keyIndex--;
        }

        Element chro = new Element("chro");

        int chargeState = Integer.parseInt(peptide.getChargeState());

        String chroText = "";

        if (conf.isHighRes()) {
            double[] samIsoArr = sampleDist.getHighMassList();
            double[] refIsoArr = refDist.getHighMassList();

            for (int i = 0; i < samIsoArr.length; i++) {
                samIsoArr[i] = (samIsoArr[i] + chargeState * PROTON_MASS) / chargeState;
            }

            for (int i = 0; i < refIsoArr.length; i++) {
                refIsoArr[i] = (refIsoArr[i] + chargeState * PROTON_MASS) / chargeState;
            }

            if (conf.isUseProline() && pepSequence.substring(2, pepSequence.length() - 2).contains("P")) {
                double[] refProlineIsoArr = new double[refIsoArr.length * 2];

                for (int i = 0; i < refIsoArr.length; i++) {
                    refProlineIsoArr[i] = refIsoArr[i];
                }

                for (int i = refIsoArr.length; i < refProlineIsoArr.length; i++) {
                    refProlineIsoArr[i] = refIsoArr[i - refIsoArr.length] + Configuration.PROLINE_SHIFT / chargeState;
                }

                chroText = CalcUtil.calculateFullMS(keyIndex, iFile, samIsoArr, refProlineIsoArr, range, null);
                chro.setText(chroText);
            } else {
                //chro.setText( CalcUtil.calculateFullMS( keyIndex, iFile, samIsoArr, refIsoArr, range) );
                //chro.setText( CalcUtil.calculateFullMS( keyIndex, iFile, samIsoArr, null, range) );
                chroText = CalcUtil.calculateFullMS(keyIndex, iFile, samIsoArr, null, range, null);
                chro.setText(chroText);
            }

        } else {

            massTolerance = conf.getMassTolerance();

            double sampleStartMass = (sampleDist.getStartMass() + chargeState * PROTON_MASS) / chargeState - massTolerance;
            double sampleEndMass = (sampleDist.getEndMass() + chargeState * PROTON_MASS) / chargeState + massTolerance;

            double refStartMass = (refDist.getStartMass() + chargeState * PROTON_MASS) / chargeState - massTolerance;
            double refEndMass = (refDist.getEndMass() + chargeState * PROTON_MASS) / chargeState + massTolerance;

            //chro.setText( CalcUtil.calculateFullMS( keyIndex, iFile, sampleStartMass, sampleEndMass, refStartMass, refEndMass, range) );
            chroText = CalcUtil.calculateFullMS(keyIndex, iFile, sampleStartMass, sampleEndMass, refStartMass, refEndMass, range);
            chro.setText(chroText);
        }

        String[] tmpStrArr = chroText.substring(0, chroText.indexOf(";")).split(" ");
        peptideEle.setAttribute("start_scan", String.valueOf(Integer.parseInt(tmpStrArr[1])));
        peptideEle.setAttribute("end_scan", String.valueOf(Integer.parseInt(tmpStrArr[2])));

        peptideEle.addContent(chro);

        return peptideEle;
    }

    public void getUnfilteredPeptideChro() {
//	chro.setText( CalcUtil.calculateFullMS( keyIndex, iFile, samIsoArr, null, range) );
    }

    //non labeling //label free
    public void createNonlabelXmlChro() throws IOException, Exception {

        //Element outRootEle = createXmlChroHeader(conf.getQuantLevel());
        ElementComposition element;
        ElementComposition totalElement;

        IndexedFile iFile;

        /**
         * ****************************************************************
         * Read DTASelect.txt file to find spectrum range for each peptide
         *****************************************************************
         */
        Hashtable<String, IndexedFile> ht = new Hashtable<String, IndexedFile>();
        Hashtable<String, IndexedFile> ht2 = new Hashtable<String, IndexedFile>(); //for tandem mass

        //DTASelectFilterReader dtaReader;
        IdentificationReader idReader;
        Set pathSet = conf.getNonlabelFilePaths();
        List fileNameList = conf.getNonlabelFilenameList();

        IsotopeReader isoReader = new IsotopeReader(this.confRootEle);
        //for tandem scans
        IsotopeTable<String, int[]> isoTable = isoReader.getIsotope();
        int[] sampleNterm = isoTable.get("sampleNTERM");
        int[] sampleCterm = isoTable.get("sampleCTERM");
        int[] refNterm = isoTable.get("refNTERM");
        int[] refCterm = isoTable.get("refCTERM");

        //build sp hashtable
        List tmpFileNameList = conf.getNonlabelFilenameList();

        Hashtable<String, int[]> spHt = new Hashtable<String, int[]>();

        int fileCount = 0;
        for (Iterator<String> itr = tmpFileNameList.iterator(); itr.hasNext();) {
            String eachPath = itr.next();
            eachPath = eachPath.substring(0, eachPath.lastIndexOf(File.separator) + 1);
            //eachPath += "DTASelect-filter.txt";
            eachPath += CensusConstants.SEARCH_OUTPUT;

            DTASelectFilterReader idReadert = new DTASelectFilterReader(eachPath);

            for (Iterator<Protein> itr1 = idReadert.getProteins(); itr1.hasNext(); ) {
                Protein protein = itr1.next();
                String accession = protein.getLocus();
                int[] tmpArr = spHt.get(accession);

                if (null == tmpArr) {
                    tmpArr = new int[tmpFileNameList.size()];
                    tmpArr[fileCount] = Integer.parseInt(protein.getSpectrumCount());

                    spHt.put(accession, tmpArr);
                } else {
                    tmpArr[fileCount] = Integer.parseInt(protein.getSpectrumCount());
                }

            }

            fileCount++;
        }

        conf.setSpHt(spHt);

        Hashtable<String, SpecRangeGenerator> specRangeHt = new Hashtable<String, SpecRangeGenerator>();
        Hashtable<String, Protein> proteinHt = new Hashtable<String, Protein>();

        System.out.println("Checking duplicate proteins and peptides...");
        //populate index hashtable
        int totalPeptideCount = 0;
        for (Iterator<String> itr = pathSet.iterator(); itr.hasNext();) {
            String path = itr.next();

            //file could be generated from either linux or window
            if (!path.endsWith("/") && !path.endsWith("\\")) {
                //for linux
                if (path.startsWith("/")) {
                    path += "/";
                } else //for window
                {
                    path += "\\";
                }
            }

            //if(1==conf.getQuantLevel())
            ht.putAll(ChroGenerator.createIndexedFiles(path, CensusConstants.MS1_FILE));

            if (2 == conf.getQuantLevel()) {
                ht2.putAll(ChroGenerator.createIndexedFiles(path, CensusConstants.MS2_FILE));
            }

            idReader = new DTASelectFilterReader(path + CensusConstants.SEARCH_OUTPUT);

            //SpecRangeGenerator rangeGen = new SpecRangeGenerator(path + "DTASelect.txt", dtaReader.isVersion2(), dtaReader.getConfidence());
            SpecRangeGenerator rangeGen = SpecRangeGenerator.getSpecRangeGenerator(idReader);

            specRangeHt.put(path, rangeGen);

            for (Iterator<Protein> itr1 = idReader.getProteins(); itr1.hasNext();) {
                Protein protein = itr1.next();

                Protein tempPro = proteinHt.get(protein.getLocus());

                if (null == tempPro) {
                    protein.populatePeptideHt(path);
                    proteinHt.put(protein.getLocus(), protein);

                    totalPeptideCount += protein.getPeptideSize();

                    continue;
                }

                totalPeptideCount += tempPro.addPeptideHt(protein.getPeptides(), path);
            }
        }

        NonLabelMappingModel mapModel = conf.getMapModel();
        if (conf.getQuantLevel() == 2) {
            mapModel.setMsmsMap(ht, ht2);
        }

        //increase status bar
        //int redundantPeptideNum = dtaReader.getRedundantPeptideNum();
        double percent = 30.0;
        double eachSeg = 70.0 / totalPeptideCount;
        int pepCount = 0;

        Element proteinEle = null;
        Element peptideEle = null;
        TIntLongHashMap index;

        Element outRootEle = createXmlChroHeader(conf.getQuantLevel());

        Element labelTypeEle = new Element("label_type");
        labelTypeEle.setText("false");

        outRootEle.addContent(labelTypeEle);

        Element alignEle = new Element("align");

        if (conf.isAlign()) {
            alignEle.setText("true");
        } else {
            alignEle.setText("false");
        }

        outRootEle.addContent(alignEle);

        for (Iterator<Element> itr2 = this.confRootEle.getChildren("sample").iterator(); itr2.hasNext();) {
            Element each = itr2.next();

            Element sample = new Element("sample");
            sample.addContent(new Element("name").addContent(each.getChildText("name")));
            Element msFile = new Element("ms_files");

            for (Iterator<Element> fileItr = each.getChild("ms_files").getChildren("file").iterator(); fileItr.hasNext();) {
                Element fileEach = fileItr.next();
                //Element fileEle = new Element("file");

                msFile.addContent(new Element("file").addContent(fileEach.getText()));
            }

            sample.addContent(msFile);

            outRootEle.addContent(sample);
        }

        if (conf.isAlign()) {
            Element ref = new Element("ref");
            Element confRefEle = this.confRootEle.getChild("ref");

            ref.addContent(new Element("sample_name").addContent(confRefEle.getChildText("sample_name")));
            ref.addContent(new Element("file_name").addContent(confRefEle.getChildText("file_name")));
            outRootEle.addContent(ref);
        }

        System.exit(0);

        /**
         * ** read spec count info ***
         */
        for (Iterator<Protein> itr = proteinHt.values().iterator(); itr.hasNext();) {
            Protein protein = itr.next();

            int[] spcArr = spHt.get(protein.getLocus());

            StringBuffer spcSb = new StringBuffer();

            for (int i = 0; i < spcArr.length; i++) {
                spcSb.append(spcArr[i]).append(",");
            }

            proteinEle = new Element("protein");
            proteinEle.setAttribute("locus", protein.getLocus());
            proteinEle.setAttribute("seq_ct", protein.getSeqCount());
            //proteinEle.setAttribute("spec_ct", protein.getSpectrumCount());
            proteinEle.setAttribute("spec_ct", spcSb.toString());
            proteinEle.setAttribute("seq_cov", protein.getSeqCoverage());
            proteinEle.setAttribute("length", protein.getLength());
            proteinEle.setAttribute("molwt", protein.getMolWt());
            proteinEle.setAttribute("pi", protein.getPI());
            proteinEle.setAttribute("val", protein.getValidation());

            try {
                proteinEle.setAttribute("desc", protein.getDescription());
            } catch (org.jdom.IllegalDataException ide) {
                proteinEle.setAttribute("desc", StringUtil.removeIsoControlChar(protein.getDescription()));
            }

            //for (Iterator<Peptide> pepItr = protein.getPeptides(); pepItr.hasNext(); )
            for (Iterator<Peptide> pepItr = protein.getPeptideHt().values().iterator(); pepItr.hasNext();) {
                Peptide peptide = pepItr.next();
                pepCount++;

                String pepSequence = peptide.getSequence();

                char[] ch = pepSequence.substring(2, peptide.getSequence().length() - 2).toCharArray();
                //element = new ElementComposition(peptide.getSequence().substring(2, peptide.getSequence().length()-2), isoReader.getIsotope());

                try {
                    //element = new ElementComposition(ch, 0, ch.length, isoReader.getIsotope());
                    totalElement = new ElementComposition(ch, 0, ch.length, isoTable);
                    totalElement.calculate();

                } catch (InvalidAAException aaEx) {

                    System.out.println("Not Quantifiable peptide : " + pepSequence);
                    percent += eachSeg;

                    if (null != progress) {
                        ChroProgressDialog.addMessageWithLine(progress, "Not Quantifiable peptide : " + pepSequence);
                        progress.setProgress((int) percent);
                    }
                    continue;

                }

                if (!totalElement.isQuantifiable()) {
                    System.out.print("\nError : ");
                    System.out.println(pepSequence + " is not quantifiable.");

                    percent += eachSeg;
                    progress.setProgress((int) percent);
                    System.out.print(pepCount);
                    System.out.print("/");
                    System.out.print(totalPeptideCount);
                    System.out.print(" peptides, ");
                    System.out.print((int) percent);
                    System.out.print(" % is complete\r");

                    if (null != this.progress) {
                        progress.addMessage("\nError : ");
                        progress.addMessage(pepSequence);
                        progress.addMessage(" is not quantifiable.\n");

                        System.out.print(String.valueOf(pepCount));
                        System.out.print("/");
                        System.out.print(String.valueOf(totalPeptideCount));
                        System.out.print(" peptides, ");
                        System.out.print(String.valueOf((int) percent));
                        System.out.print(" % is complete\r");
                    }

                    continue;
                }

                int chargeState = Integer.parseInt(peptide.getChargeState());

                IsotopeDist sampleDist = null;
                //IsotopeDist refDist = null;

                String fileName = peptide.getFileName();

                for (int i = 0; i < 3; i++) {
                    fileName = fileName.substring(0, fileName.lastIndexOf("."));
                }

                String fpath = peptide.getFilePath();
                SpecRangeGenerator rangeGen = specRangeHt.get(fpath);

                StringBuffer rangeKey = new StringBuffer();
                rangeKey.append(protein.getLocus());
                rangeKey.append(fileName);
                rangeKey.append(peptide.getSequence().substring(2, peptide.getSequence().length() - 2));

                SpecRange range = rangeGen.getSpecRange(rangeKey.toString());
                //sb.append("\tStartScan\tEndScan\tDTAPeakStart\tDTAPeakEnd\n");

                peptideEle = this.createXmlChroPeptideTitle(false, peptide);

                if (null == range) {
                    int tmpScanNum = Integer.parseInt(peptide.getScanNum());
                    range = new SpecRange(tmpScanNum, tmpScanNum);
                    peptideEle.setAttribute("start_scan", peptide.getScanNum());
                    peptideEle.setAttribute("end_scan", peptide.getScanNum());
                } else {
                    peptideEle.setAttribute("start_scan", String.valueOf(range.getMin()));
                    peptideEle.setAttribute("end_scan", String.valueOf(range.getMax()));
                }

                Element chro = new Element("chro");

                if (1 == conf.getQuantLevel()) {
                    peptideEle = this.createXmlChroPeptideTitle(true, peptide); //true for full scan
                    fileName += ".";
                    fileName += CensusConstants.MS1_FILE;
                    iFile = ht.get(fpath + fileName);
                } else //level is 2
                {
                    peptideEle = this.createXmlChroPeptideTitle(false, peptide);
                    fileName += ".";
                    fileName += CensusConstants.MS2_FILE;
                    iFile = ht2.get(fpath + fileName);
                }

                peptideEle.setAttribute("start_scan", String.valueOf(range.getMin()));
                peptideEle.setAttribute("end_scan", String.valueOf(range.getMax()));

                if (null == iFile) {
                    System.out.println("Error : cannot find the file " + fpath + fileName);
                    throw new IOException("Error : cannot find the file " + fpath + fileName);
                    //System.exit(0);
                }

                int[] keys = iFile.getKeys();
                int keyIndex = Arrays.binarySearch(keys, Integer.parseInt(peptide.getScanNum()));

                if (keyIndex < 0) //Cannot find index
                {
                    keyIndex = -(++keyIndex); //Math.abs(++keyIndex);
                }
                if (keyIndex >= keys.length) {
                    keyIndex--;
                }

                if (conf.getQuantLevel() == 1) {

                    sampleDist = new IsotopeDist(totalElement.getElementSampleArr(), totalElement.getModShift(), true);

                    if (conf.isHighRes()) {
                        double[] samIsoArr = sampleDist.getHighMassList();
                        //double[] refIsoArr = refDist.getHighMassList();

                        for (int i = 0; i < samIsoArr.length; i++) {
                            samIsoArr[i] = (samIsoArr[i] + chargeState * PROTON_MASS) / chargeState;
                        }

                        try {

                            chro.setText(CalcUtil.calculateNonlabelMS(keyIndex, fileName, ht, pathArray, samIsoArr, range, fpath));
                            //System.out.println( CalcUtil.calculateNonlabelMS(keyIndex, fileName, ht, pathArray, samIsoArr, range, path) );
                        } catch (CensusIndexOutOfBoundException ciob) {
                            percent += eachSeg;
                            System.out.print("Warning: the peptide is not quantifiable - the peptide seems to be outside of chromatogram alignment range.\r");
                            continue;
                        }

                    } else {
                        double sampleStartMass = sampleStartMass = (sampleDist.getStartMass() + chargeState * PROTON_MASS) / chargeState - massTolerance;
                        double sampleEndMass = sampleEndMass = (sampleDist.getEndMass() + chargeState * PROTON_MASS) / chargeState + massTolerance;

                        chro.setText(CalcUtil.calculateNonlabelMS(keyIndex, fileName, ht, pathArray, sampleStartMass, sampleEndMass, range, fpath));

                        //      System.exit(0);
//                        throw new Exception("not supporting low resolution yet. Robin.");
                        //refStartMass = (refDist.getStartMass() + chargeState*PROTON_MASS) / chargeState - massTolerance;
                        //refEndMass = (refDist.getEndMass() + chargeState*PROTON_MASS) / chargeState + massTolerance;
                        //chro.setText( CalcUtil.calculateFullMS( keyIndex, iFile, sampleStartMass, sampleEndMass, refStartMass, refEndMass, conf, range) );
                    }

                    peptideEle.addContent(chro);
                    proteinEle.addContent(peptideEle);

                } else if (conf.getQuantLevel() == 2) {

                    double[][] bionSample;
                    double[][] bionRef;
                    double[][] yionSample;
                    double[][] yionRef;

                    int pepLength = 0;

                    for (int i = 0; i < ch.length; i++) {
                        if (ch[i] == '*' || ch[i] == '@' || ch[i] == '#') {
                            continue;
                        }

                        pepLength++;
                    }

                    bionSample = new double[pepLength][chargeState * 3];
                    bionRef = new double[pepLength][chargeState * 3];
                    //Yions
                    yionSample = new double[pepLength][chargeState * 3];
                    yionRef = new double[pepLength][chargeState * 3];

                    int pepIndex = 0;

                    for (int i = 0; i < ch.length; i++) {
                        if (ch[i] == '*' || ch[i] == '@' || ch[i] == '#') {
                            continue;
                        }

                        try {
                            element = new ElementComposition(ch, 0, i + 1, isoTable);
                            element.calculate();
                        } catch (InvalidAAException ive) {
                            System.out.println("Not Quantifiable peptide : " + pepSequence);

                            percent += eachSeg;
                            if (null != progress) {
                                ChroProgressDialog.addMessageWithLine(progress, "Not Quantifiable peptide : " + pepSequence);
                                progress.setProgress((int) percent);

                            }
                            continue;
                        }

                        //Y ions
                        sampleDist = new IsotopeDist(
                                getComplementaryComposition(totalElement.getElementSampleArr(), element.getElementSampleArr(), sampleNterm, sampleCterm), element.getModShift(), true);

                        //refDist = new IsotopeDist(
                        //        getComplementaryComposition(totalElement.getElementRefArr(), element.getElementRefArr(), refNterm, refCterm), false);   //fix this
                        switch (chargeState) {
                            case 3:
                                yionSample[(pepIndex + 1) % pepLength][8] = (sampleDist.getEndMass() + 3 * PROTON_MASS) / 3 + massTolerance;
                                yionSample[(pepIndex + 1) % pepLength][7] = (sampleDist.getStartMass() + 3 * PROTON_MASS) / 3 - massTolerance;
                                yionSample[(pepIndex + 1) % pepLength][6] = (sampleDist.getAvgMass() + 3 * PROTON_MASS) / 3;
                            //      yionRef[(pepIndex+1)%pepLength][8] = (refDist.getEndMass()+3*PROTON_MASS)/3+massTolerance;
                            //      yionRef[(pepIndex+1)%pepLength][7] = (refDist.getStartMass()+3*PROTON_MASS)/3-massTolerance;
                            //      yionRef[(pepIndex+1)%pepLength][6] = (refDist.getAvgMass()+3*PROTON_MASS)/3;

                            case 2:
                                yionSample[(pepIndex + 1) % pepLength][5] = (sampleDist.getEndMass() + 2 * PROTON_MASS) / 2 + massTolerance;
                                yionSample[(pepIndex + 1) % pepLength][4] = (sampleDist.getStartMass() + 2 * PROTON_MASS) / 2 - massTolerance;
                                yionSample[(pepIndex + 1) % pepLength][3] = (sampleDist.getAvgMass() + 2 * PROTON_MASS) / 2;
                            //      yionRef[(pepIndex+1)%pepLength][5] = (refDist.getEndMass()+2*PROTON_MASS)/2+massTolerance;
                            //      yionRef[(pepIndex+1)%pepLength][4] = (refDist.getStartMass()+2*PROTON_MASS)/2-massTolerance;
                            //      yionRef[(pepIndex+1)%pepLength][3] = (refDist.getAvgMass()+2*PROTON_MASS)/2;

                            case 1:
                                yionSample[(pepIndex + 1) % pepLength][2] = sampleDist.getEndMass() + 1 * PROTON_MASS + massTolerance; //add proton to give b fragment ion
                                yionSample[(pepIndex + 1) % pepLength][1] = sampleDist.getStartMass() + 1 * PROTON_MASS - massTolerance;
                                yionSample[(pepIndex + 1) % pepLength][0] = sampleDist.getAvgMass() + 1 * PROTON_MASS;
                            //     yionRef[(pepIndex+1)%pepLength][2] = refDist.getEndMass()+1*PROTON_MASS+massTolerance;
                            //     yionRef[(pepIndex+1)%pepLength][1] = refDist.getStartMass()+1*PROTON_MASS-massTolerance;
                            //     yionRef[(pepIndex+1)%pepLength][0] = refDist.getAvgMass()+1*PROTON_MASS;

                            default:
                                break;
                        }

                        element.calculateBion();
                        //element.printComposition();

                        sampleDist = new IsotopeDist(element.getElementSampleArr(), element.getModShift(), true);

                        switch (chargeState) {
                            case 3:
                                bionSample[pepIndex][8] = (sampleDist.getEndMass() + 2 * PROTON_MASS) / 3 + massTolerance;
                                bionSample[pepIndex][7] = (sampleDist.getStartMass() + 2 * PROTON_MASS) / 3 - massTolerance;
                                bionSample[pepIndex][6] = (sampleDist.getAvgMass() + 2 * PROTON_MASS) / 3;
                            //      bionRef[pepIndex][8] = (refDist.getEndMass()+2*PROTON_MASS)/3+massTolerance;
                            //      bionRef[pepIndex][7] = (refDist.getStartMass()+2*PROTON_MASS)/3-massTolerance;
                            //      bionRef[pepIndex][6] = (refDist.getAvgMass()+2*PROTON_MASS)/3;

                            case 2:
                                bionSample[pepIndex][5] = (sampleDist.getEndMass() + 1 * PROTON_MASS) / 2 + massTolerance;
                                bionSample[pepIndex][4] = (sampleDist.getStartMass() + 1 * PROTON_MASS) / 2 - massTolerance;
                                bionSample[pepIndex][3] = (sampleDist.getAvgMass() + 1 * PROTON_MASS) / 2;
                            //      bionRef[pepIndex][5] = (refDist.getEndMass()+1*PROTON_MASS)/2+massTolerance;
                            //      bionRef[pepIndex][4] = (refDist.getStartMass()+1*PROTON_MASS)/2-massTolerance;
                            //      bionRef[pepIndex][3] = (refDist.getAvgMass()+1*PROTON_MASS)/2;

                            case 1:
                                bionSample[pepIndex][2] = sampleDist.getEndMass() + massTolerance;
                                bionSample[pepIndex][1] = sampleDist.getStartMass() - massTolerance;
                                bionSample[pepIndex][0] = sampleDist.getAvgMass();
                            //     bionRef[pepIndex][2] = refDist.getEndMass()+massTolerance;
                            //     bionRef[pepIndex][1] = refDist.getStartMass()-massTolerance;
                            //     bionRef[pepIndex][0] = refDist.getAvgMass();

                            default:
                                break;

                        }

                        pepIndex++;
                    }

                    try {
                        element = new ElementComposition(ch, 0, ch.length, isoTable);
                        element.calculate();
                    } catch (InvalidAAException ive) {
                        System.out.println("Not Quantifiable peptide : " + pepSequence);

                        percent += eachSeg;
                        if (null != progress) {
                            ChroProgressDialog.addMessageWithLine(progress, "Not Quantifiable peptide : " + pepSequence);
                            progress.setProgress((int) percent);

                        }
                        continue;
                    }
                    sampleDist = new IsotopeDist(element.getElementSampleArr(), element.getModShift(), true);

                    int[] tempA = element.getElementSampleArr();

//                    for(int i=0;i<tempA.length;i++)
                    switch (chargeState) {
                        case 3:
                            yionSample[0][8] = (sampleDist.getEndMass() + 3 * PROTON_MASS) / 3 + massTolerance;
                            yionSample[0][7] = (sampleDist.getStartMass() + 3 * PROTON_MASS) / 3 - massTolerance;
                            yionSample[0][6] = (sampleDist.getAvgMass() + 3 * PROTON_MASS) / 3;
                        //      yionRef[0][8] = (refDist.getEndMass()+3*PROTON_MASS)/3+massTolerance;
                        //      yionRef[0][7] = (refDist.getStartMass()+3*PROTON_MASS)/3-massTolerance;
                        //      yionRef[0][6] = (refDist.getAvgMass()+3*PROTON_MASS)/3;

                        case 2:
                            yionSample[0][5] = (sampleDist.getEndMass() + 2 * PROTON_MASS) / 2 + massTolerance;
                            yionSample[0][4] = (sampleDist.getStartMass() + 2 * PROTON_MASS) / 2 - massTolerance;
                            yionSample[0][3] = (sampleDist.getAvgMass() + 2 * PROTON_MASS) / 2;
                        //      yionRef[0][5] = (refDist.getEndMass()+2*PROTON_MASS)/2+massTolerance;
                        //      yionRef[0][4] = (refDist.getStartMass()+2*PROTON_MASS)/2-massTolerance;
                        //      yionRef[0][3] = (refDist.getAvgMass()+2*PROTON_MASS)/2;

                        case 1:
                            yionSample[0][2] = sampleDist.getEndMass() + 1 * PROTON_MASS + massTolerance; //add proton to give b fragment ion
                            yionSample[0][1] = sampleDist.getStartMass() + 1 * PROTON_MASS - massTolerance;
                            yionSample[0][0] = sampleDist.getAvgMass() + 1 * PROTON_MASS;
                        //      yionRef[0][2] = refDist.getEndMass()+1*PROTON_MASS+massTolerance;
                        //      yionRef[0][1] = refDist.getStartMass()+1*PROTON_MASS-massTolerance;
                        //      yionRef[0][0] = refDist.getAvgMass()+1*PROTON_MASS;

                        default:
                            break;

                    }

//                    peptideEle = this.createXmlChroPeptideTitle(false, peptide);
                    //String ms2FileName = peptide.getFileName().substring(0, peptide.getFileName().indexOf("."));
                    // iFile = ht.get( this.filePath + fileName + "." + "ms2");
                    //keys = iFile.getKeys();
                    // keyIndex = Arrays.binarySearch(keys, Integer.parseInt(peptide.getScanNum()));
                    // if(keyIndex<0) //Cannot find index
                    //     keyIndex=-(++keyIndex); //Math.abs(++keyIndex);
                    //System.out.println("key" + keyIndex);
                    /**
                     * Find start scan number same as a precursor of
                     * sampleAvgMass Then, later the program will find following
                     * scan # quicker.
                     *
                     */
                    /* End of finding start spectrum number */
 /*
                    samplePrecursor=1;
                    refPrecursor=1;

                    TIntDoubleHashMap precursorMap = iFile.getPrecursorMap();
                     */
                    double calcSamMass = (sampleDist.getAvgMass() + chargeState * PROTON_MASS) / chargeState;
                    // double calcRefMass = (refDist.getAvgMass()+chargeState*PROTON_MASS)/chargeState;

                    conf.setCalcSamAvgMass(calcSamMass);
                    peptideEle.setAttribute("lightStartMass", String.valueOf(sampleDist.getStartMass()));
                    peptideEle.setAttribute("lightAvgMass", String.valueOf(conf.getCalcSamAvgMass()));
                    //conf.setCalcRefAvgMass(calcRefMass);

                    if (conf.isHighRes()) {
                        System.out.println("high resolution for tandem spectra is not supported yet.  Please change it to low resolution or contact Robin by rpark@scripps.edu");
                        System.exit(0);
                    } else { //low resolution
                        double[] samIsoArr = sampleDist.getHighMassList();
                        //double[] refIsoArr = refDist.getHighMassList();

                        for (int i = 0; i < samIsoArr.length; i++) {
                            samIsoArr[i] = (samIsoArr[i] + chargeState * PROTON_MASS) / chargeState;
                            //System.out.println("\t\t" + samIsoArr[i]);
                        }

                        try {
                            String outStr = CalcUtil.calculateNonlabelMS(keyIndex, fileName, ht2, pathArray, samIsoArr, range, fpath, bionSample, yionSample, iFile);

                            if (outStr == null) {
                                percent += eachSeg;
                                continue;
                            }

                            chro.setText(outStr);
                        } catch (CensusIndexOutOfBoundException ciob) {
                            percent += eachSeg;
                            System.out.print("Warning: the peptide is not quantifiable - the peptide seems to be outside of chromatogram alignment range.\r");

                            continue;
                        } catch (Exception e) {
                            percent += eachSeg;
                            System.out.print("Warning: the peptide is not quantifiable.\r");

                            continue;
                        }

                    }

                    peptideEle.addContent(chro);

                    Element fragEle = new Element("frag");

                    Element bSample = new Element("bs");
                    StringBuffer tempSb = new StringBuffer();
                    for (int i = 0; i < bionSample.length; i++) {
                        int j;
                        for (j = 0; j < bionSample[i].length - 1; j++) {
                            tempSb.append(formatter.format(bionSample[i][j])).append(" ");
                        }

                        tempSb.append(formatter.format(bionSample[i][j])).append(",");
                    }

                    bSample.setText(tempSb.toString());
                    fragEle.addContent(bSample);

                    Element bRef = new Element("br");
                    tempSb.delete(0, tempSb.length());
                    for (int i = 0; i < bionRef.length; i++) {
                        int j;
                        for (j = 0; j < bionRef[i].length - 1; j++) {
                            tempSb.append(formatter.format(bionRef[i][j])).append(" ");
                        }

                        tempSb.append(formatter.format(bionRef[i][j])).append(",");

                    }

                    bRef.setText(tempSb.toString());
                    fragEle.addContent(bRef);

                    Element ySample = new Element("ys");
                    tempSb.delete(0, tempSb.length());
                    for (int i = 0; i < yionSample.length; i++) {
                        int j;
                        for (j = 0; j < yionSample[i].length - 1; j++) {
                            tempSb.append(formatter.format(yionSample[i][j])).append(" ");
                        }

                        tempSb.append(formatter.format(yionSample[i][j])).append(",");
                    }

                    ySample.setText(tempSb.toString());
                    fragEle.addContent(ySample);

                    Element yRef = new Element("yr");
                    tempSb.delete(0, tempSb.length());
                    for (int i = 0; i < yionRef.length; i++) {
                        int j;
                        for (j = 0; j < yionRef[i].length - 1; j++) {
                            tempSb.append(formatter.format(yionRef[i][j])).append(" ");
                        }

                        tempSb.append(formatter.format(yionRef[i][j])).append(",");
                    }

                    yRef.setText(tempSb.toString());
                    fragEle.addContent(yRef);

                    peptideEle.addContent(fragEle);

                    //peptideEle.addContent(chro);
                    proteinEle.addContent(peptideEle);

                }

                percent += eachSeg;

                System.out.print(pepCount);

                System.out.print("/");
                System.out.print(totalPeptideCount);
                System.out.print(" peptides, ");
                System.out.print((int) percent);
                System.out.print(" % is complete\r");

                if (null != progress) {
                    progress.setProgress((int) percent);
                    //progress.addMessageOffset(pepCount + "/" + totalPeptideCount + " Peptide, " + (int)percent +  "% is complete");
                }
            }

            if (proteinEle.getChildren().size() > 0) {
                outRootEle.addContent(proteinEle);
            }

        }

        Document doc = new Document(outRootEle);
        OutputStream os = null;

        if (conf.getFilePath().endsWith(File.separator)) {
            os = new FileOutputStream(conf.getFilePath() + "census_chro.xml");
        } else {
            os = new FileOutputStream(conf.getFilePath() + File.separator + "census_chro.xml");
        }

        XMLOutputter outputter = new XMLOutputter();
        outputter.setFormat(Format.getPrettyFormat());
        outputter.output(doc, os);
        os.close();

        System.out.println("\n100% complete");
        ChroProgressDialog.addMessageWithLine(progress, "\n100% complete");
        System.out.println("done.");
        ChroProgressDialog.addMessageWithLine(progress, "done.");
    }


    public static SpectraDB connectCreateSpectraDB(String filePath, File spectraDir, String path, String name) throws SQLException, IOException {
        String sqliteDBPath = path+".sqlite";
        File sqliteDB = new File(sqliteDBPath);
        if(sqliteDB.exists()  && sqliteDB.length()>0)
        {
            SpectraDB db = SpectraDB.connectToDBReadOnly(sqliteDBPath);
            return db;
        }
        else if(spectraDir.exists()){
            String spectraSqliteDBPath = spectraDir.getAbsolutePath() +File.separatorChar + name + ".sqlite";
            File spectraSqliteDB = new File(spectraSqliteDBPath);
            if(spectraSqliteDB.exists()  && spectraSqliteDB.length()>0)
            {
                SpectraDB spectraDB = SpectraDB.connectToDBReadOnly(spectraSqliteDB.getAbsolutePath());
                return spectraDB;
            }
            else
            {
                CreateDb.createNewDatabase(spectraDir.getAbsolutePath(),name+".sqlite",name);
                SpectraDB spectraDB = SpectraDB.connectToDBReadOnly(spectraSqliteDB.getAbsolutePath());
                return spectraDB;
            }
        }
        else if(!spectraDir.exists())
        {
            CreateDb.createNewDatabase(filePath,sqliteDB.getName(),name);
            SpectraDB spectraDB = SpectraDB.connectToDBReadOnly(sqliteDBPath);
            return spectraDB;
        }
        return null;
    }

    public static SpectraDB connectCreateSpectraDB(String filePath, File spectraDir, File ms2File) throws SQLException, IOException {
        String parentPath = ms2File.getParentFile().getAbsolutePath();
        if(isHeavyFile(ms2File.getName(),parentPath) || isMediumFile(ms2File.getName(), parentPath)
                || isLightFile(ms2File.getName(), parentPath))
        {
            ms2File = new File(parentPath + File.separator + ms2File.getName().substring(1));
        }
        String sqliteDBPath = ms2File.getAbsolutePath()+".sqlite";
        File sqliteDB = new File(sqliteDBPath);
        if(sqliteDB.exists()  && sqliteDB.length()>0)
        {
            SpectraDB db = SpectraDB.connectToDBReadOnly(sqliteDBPath);
            return db;
        }
        else if(spectraDir.exists()){
            String spectraSqliteDBPath = spectraDir.getAbsolutePath() +File.separatorChar + ms2File.getName() + ".sqlite";
            File spectraSqliteDB = new File(spectraSqliteDBPath);
            if(spectraSqliteDB.exists()  && spectraSqliteDB.length()>0)
            {
                SpectraDB spectraDB = SpectraDB.connectToDBReadOnly(spectraSqliteDB.getAbsolutePath());
                return spectraDB;
            }
            else
            {
                CreateDb.createNewDatabase(spectraDir.getAbsolutePath(),ms2File.getName()+".sqlite",ms2File.getName());
                SpectraDB spectraDB = SpectraDB.connectToDBReadOnly(spectraSqliteDB.getAbsolutePath());
                return spectraDB;
            }
        }
        else if(!spectraDir.exists())
        {
            CreateDb.createNewDatabase(filePath,sqliteDB.getName(),ms2File.getName());
            SpectraDB spectraDB = SpectraDB.connectToDBReadOnly(sqliteDBPath);
            return spectraDB;
        }
        return null;
    }

    public static Map<String, SpectraDB> connectCreateSpectraDBIndex(String filePath, String extension) throws Exception {
        Map<String,SpectraDB> result = new HashMap<>();
        File spectraDir = new File(filePath+"/../../spectra/");
        File currentDir = new File(filePath);
        File [] arr = currentDir.listFiles(new RelExFileFilter(extension+".index"));
        for(File ms2File: arr)
        {
            String path  = ms2File.getAbsolutePath().replaceAll("\\.index", "");
            String name = ms2File.getName().replaceAll("\\.index","");
            SpectraDB db = connectCreateSpectraDB(filePath,spectraDir,path,name);
            result.put( name, db);
        }
        return result;
    }



    public static Map<String, SpectraDB> connectCreateSpectraDB(String filePath, String extension) throws Exception {
        Map<String,SpectraDB> result = new HashMap<>();
        File spectraDir = new File(filePath+"/../../spectra/");
        File currentDir = new File(filePath);
        File [] arr = currentDir.listFiles(new RelExFileFilter(extension));
        if(arr.length==0 && spectraDir !=null)
            arr = spectraDir.listFiles(new RelExFileFilter(extension));

        if(arr !=null)
        {
            for(File ms2File: arr)
            {

              //  System.out.println("<<>> " + ms2File.getAbsolutePath());
                SpectraDB db = connectCreateSpectraDB(filePath,spectraDir,ms2File);
                result.put( ms2File.getName(), db);
            }
        }
       return result;
    }



    public static Map<String, SpectraDB> connectSpectraDB(String filePath, String extension) throws SQLException {
        return connectSpectraDB(filePath,extension, new HashMap<>());
    }


    public static Map<String, SpectraDB> connectSpectraDB(String filePath, String extension, Map<String,SpectraDB> dbMap)
                throws SQLException {
        File f = new File(filePath);
        if(f.exists() && f.isDirectory())
        {
            File[] list = f.listFiles(new RelExFileFilter(extension+".sqlite"));
            String indexFileName;
            File indexFile;

            if (null == list || list.length <= 0) {

                return dbMap;
                //throw new CensusGeneralException("Error: Spectral files are ot found at " + filePath + " If you use mzXML, please use option '-x'");
            }

            if (null == conf) {
                conf = Configuration.getInstance();
            }
            for(File msFile : list)
            {
                String path =  msFile.getAbsolutePath();
                String name = msFile.getName().replace(".sqlite","");
                if(!dbMap.containsKey(name))
                {
                    SpectraDB db = SpectraDB.connectToDBReadOnly(path);
                    dbMap.put(name,db);
                }
            }
        }
        return dbMap;
    }

    public static void createSpectraDB(String filePath, String extension) throws Exception {
        File f = new File(filePath);


        File[] list = f.listFiles(new RelExFileFilter(extension));
        String indexFileName;
        File indexFile;

        Map<String, SpectraDB> dbMap = new HashMap<>();
        if (null == list || list.length <= 0) {

            return ;
            //throw new CensusGeneralException("Error: Spectral files are ot found at " + filePath + " If you use mzXML, please use option '-x'");
        }

        if (null == conf) {
            conf = Configuration.getInstance();
        }
        for(File msFile : list)
        {
            String name  = msFile.getAbsolutePath()+".sqlite";
            File sqliteDB = new File(name);
            if(!sqliteDB.exists())
            {
                CreateDb.createNewDatabase("",name,msFile.getName());
            }

        }
        return;
    }




    public static Hashtable spectraSqliteDB(String filePath, String extension) throws IOException, CensusGeneralException {
        return createIndexedFilesHelper(filePath, extension, true, false,false);
    }
    public static Hashtable createIndexedFilesNoMs(String filePath, String extension) throws IOException, CensusGeneralException {
        return createIndexedFilesHelper(filePath, extension, true, false,true);
    }
    public static Hashtable createIndexedFiles(String filePath, String extension) throws IOException, CensusGeneralException {
        return createIndexedFilesHelper(filePath, extension, true, false,false);
    }

    public static Hashtable createIndexedFiles(String filePath, String extension, boolean includePath, boolean indexCheck) throws IOException, CensusGeneralException {
            return createIndexedFilesHelper(filePath, extension, true, false,false);
    }
        //robin
    public static Hashtable createIndexedFilesHelper(String filePath, String extension, boolean includePath, boolean indexCheck, boolean indexOnly) throws IOException, CensusGeneralException {

        //for non-labeling quantification, same file names may appear in different folders.  That means
        // it is dangerous to use file names as keys.  Fix me later.

        if (!filePath.endsWith(File.separator)) {
            filePath += File.separator;
        }

        File f = new File(filePath);


        String[] list = f.list(new RelExFileFilter(indexOnly ? extension+".index" : extension));

        Hashtable<String, IndexedFile> ht = new Hashtable<String, IndexedFile>();
        IndexedFile iFile = null;

        String indexFileName;
        File indexFile;

        if (null == list || list.length <= 0) {

            return new Hashtable<String, IndexedFile>();
            //throw new CensusGeneralException("Error: Spectral files are ot found at " + filePath + " If you use mzXML, please use option '-x'");
        }

        if (null == conf) {
            conf = Configuration.getInstance();
        }

        switch (conf.getSpectrumFormat()) {
            case Configuration.MS_FILE_FORMAT:
                if(indexCheck)
                  System.out.print("creating index file");
                for (int i = 0; i < list.length; i++) {
                    indexFileName = filePath + list[i] + (indexOnly ? "" : ".index");
                    indexFile = new File(indexFileName);

                    //System.out.println(indexFileName + " " + indexFile.exists());
                    //Create index file
                    //if (indexCheck && (!indexFile.exists() || indexFile.length() <= 0)) {
                  if(indexCheck) {
                    System.out.print(".");
                    ChroProgressDialog.addMessageWithLine(progress, "creating index file " + indexFileName);
                  }

                    //System.out.println("uncomment below... robin");
                     //   MSIndexFileCreator.createIndexFile(filePath + list[i]); //text.append("Index files are required");


                    //System.out.println("Reading index file for " +  indexFile + " & " + list[i] + " & " + filePath);
                    //if index file is corrupted, delete the file and re-try to read it.
                    try {
                        iFile = new IndexedFile(indexFile, filePath + list[i]);
                    } catch (FileNotFoundException fnfe) {
                        System.out.println("Error: Spectral files are not found.  If you use mzXML, please use option '-x'");
                        fnfe.printStackTrace();
                        System.exit(0);
                    } catch (Exception ioe) {
                        System.out.println("re-creating index file " + indexFileName);

                        ChroProgressDialog.addMessageWithLine(progress, "re-creating index file " + indexFileName);
                        MSIndexFileCreator.createIndexFile(filePath + list[i]); //text.append("Index files are required");
                        iFile = new IndexedFile(indexFile, filePath + list[i]);
                    }

                    if (includePath) {
                        ht.put(filePath + list[i], iFile);
                    } else {
                        ht.put(list[i], iFile);
                    }

                }

                System.out.println("");
                break;

            case Configuration.MZXML_FILE_FORMAT:
                for (int i = 0; i < list.length; i++) {
                    try {
                        System.out.print("Indexing on " + list[i] + "...");
                        ChroProgressDialog.addMessageWithLine(progress, "Indexing on " + list[i] + "...");
                        iFile = new IndexedFile(filePath + list[i]);
                        System.out.println("done");
                        ChroProgressDialog.addMessageWithLine(progress, "done");

                    } catch (FileNotFoundException fnfe) {
                        System.out.println("Error: Spectral file is not found");

                        System.exit(0);
                    } catch (Exception ioe) {

                        System.out.println("re-creating index file "); // + indexFileName);

                        //ChroProgressDialog.addMessageWithLine(progress, "re-creating index file " + indexFileName);
                        //MSIndexFileCreator.createIndexFile(filePath + list[i]); //text.append("Index files are required");
                    }

                    //System.out.println(ht + " " + iFile + " " + list[i]);
//		    ht.put(filePath + list[i], iFile);
                    if (includePath) {
                        ht.put(filePath + list[i], iFile);
                    } else {
                        ht.put(list[i], iFile);
                    }
                }

                break;

            case Configuration.MS2_FILE_FORMAT:
                for (int i = 0; i < list.length; i++) {
                    indexFileName = filePath + list[i] + (indexOnly ? "" : ".index");
                    indexFile = new File(indexFileName);
                    String path = filePath + (indexOnly? list[i].replace(".index","") : list[i]);

                    if (!indexFile.exists() || indexFile.length() <= 0) {
                        System.out.println("creating index file " + indexFileName);

                        ChroProgressDialog.addMessageWithLine(progress, "creating index file " + indexFileName);
                        MSIndexFileCreator.createIndexFile(filePath + list[i]); //text.append("Index files are required");
                    }

                    try {
                        iFile = new IndexedFile(indexFile,
                                path,
                                false);
                    } catch (FileNotFoundException fnfe) {
                        System.out.println("Error: Spectral file is not found");

                        System.exit(0);
                    } catch (Exception ioe) {
                        ioe.printStackTrace();
                        System.out.println(ioe);
                        System.out.println("re-creating index file " + indexFileName);

                        ChroProgressDialog.addMessageWithLine(progress, "re-creating index file " + indexFileName);
                        MSIndexFileCreator.createIndexFile(filePath + list[i]); //text.append("Index files are required");
                    }

                    if (includePath) {
                        ht.put(path, iFile);
                    } else {
                        ht.put(list[i], iFile);
                    }
                    //ht.put(filePath + list[i], iFile);
                }

                break;
            default:
                break;
        }

        return ht;
    }

    private int[] getComplementaryComposition(int[] totalEle, int[] element, int[] nterm, int[] cterm) {
        int[] temp = new int[element.length];

        for (int i = 0; i < element.length; i++) {
            temp[i] = totalEle[i] - element[i] + nterm[i] + cterm[i];
        }

        return temp;
    }

    /*
    private static Element createXmlChroHeader(boolean isDataDependent)
    {
        return createXmlChroHeader( isDataDependent?0:1 );
    }
     */
    private static Element createXmlChroHeader(int quantLevel) {
        return createXmlChroHeader(quantLevel, -1, null);
    }

    private static Element createXmlChroHeader(int quantLevel, String quantType) {
        return createXmlChroHeader(quantLevel, -1, quantType);
    }

    private static Element createXmlChroHeader(int quantLevel, int expType) {
        return createXmlChroHeader(quantLevel, expType, null);
    }

    //quantLevel 1 for MS1, 2 for MS2, and so on
    private static Element createXmlChroHeader(int quantLevel, int expType, String quantType) {
        //expType 13: iTRAQ single spectrum
        //expType 14: iTRAQ multiple spectra
        Element root = new Element("relex_chro");

        Element version = new Element("version");
        version.setText("Census v. " + conf.getVersion() + " Chro file");
        root.addContent(version);

        Element author = new Element("author");
        author.setAttribute("name", "Robin, Sung Kyu Park");
        author.setAttribute("email", "rpark@scripps.edu");
        root.addContent(author);

        author = new Element("author");
        author.setAttribute("name", "John Venable");
        author.setAttribute("email", "jvenable@gnf.org");
        root.addContent(author);

        Element createdDate = new Element("created_date");
        createdDate.setText(new Date().toString());
        root.addContent(createdDate);

        Element dataDependency = new Element("data_dependency");
        dataDependency.setText(String.valueOf(quantLevel));
        /*
        if(isDataDependent)
            dataDependency.setText("0");
        else
            dataDependency.setText("1");
         */
        Element quantLevelEle = new Element("quantLevel");
        quantLevelEle.setText(String.valueOf(quantLevel));

        if (null != quantType) {
            Element quantTypeEle = new Element("quantType");
            quantTypeEle.setText(quantType);
            root.addContent(quantTypeEle);

        }


        //System.out.println("exp type==============" + expType);
        if (expType > 0) {
            Element expTypeEle = new Element("experiment_type");
            expTypeEle.setText(String.valueOf(expType));
            root.addContent(expTypeEle);
        }

        root.addContent(dataDependency);
        root.addContent(quantLevelEle);

        Element smoothEle = new Element("smooth");
        smoothEle.setText("" + conf.isSmooth());
        root.addContent(smoothEle);


        //0:data dependent, 1:data independent\n");
        return root;
    }

    private String createChroHeader(boolean isDataDependent) {
        StringBuffer sb = new StringBuffer();
        sb.append("H\t").append("Census v. 2.0 Chro file\n");
        sb.append("H\t").append("Create by John Venable jvenable@scripps.edu\n");
        sb.append("H\t").append("Robin, Sung Kyu Park rpark@scripps.edu\n");
//        sb.append("H\t").append("Michael J. MacCoss \n");
        sb.append("H\t").append("The Scripps Research Institute, La Jolla, CA\n");
        sb.append("H\t").append("created date\t").append(new Date()).append("\n");
        sb.append("H\t").append("data_dependency\t");

        if (isDataDependent) {
            sb.append("0");
        } else {
            sb.append("1");
        }

        sb.append("\t; 0:data dependent, 1:data independent\n");
        sb.append("H\t\n");

        return sb.toString();

    }

    private Element createXmlChroPeptideTitle(boolean isDataDependent, Peptide peptide) {
        Element peptideEle = new Element("peptide");
        peptideEle.setAttribute("unique", peptide.isUnique() ? "*" : "");
        peptideEle.setAttribute("file", peptide.getFileName());
        peptideEle.setAttribute("scan", peptide.getScanNum());
        peptideEle.setAttribute("seq", peptide.getSequence());
        peptideEle.setAttribute("xcorr", peptide.getXCorr());
        peptideEle.setAttribute("calcMHplus", peptide.getCalcMHplus());
        peptideEle.setAttribute("MHplus", peptide.getMhPlus());
        peptideEle.setAttribute("totalIntensity", peptide.getTotalIntensity());
        peptideEle.setAttribute("spRank", peptide.getSpRank());
        peptideEle.setAttribute("spScore", peptide.getSpScore());
        peptideEle.setAttribute("redundancy", peptide.getRedundancy());
        peptideEle.setAttribute("deltaCN", peptide.getDeltCN());
        peptideEle.setAttribute("deltaMass", String.valueOf(peptide.getDeltaMass()));
        if (null == peptide.getDeltCN()) {
            peptideEle.setAttribute("deltaCN", "");
        }

        peptideEle.setAttribute("charge", peptide.getChargeState());
        peptideEle.setAttribute("spC", peptide.getRedundancy());

        Hashtable<String, String> scoreHt = peptide.getScoreHt();
        for (Iterator<String> itrScr = scoreHt.keySet().iterator(); itrScr.hasNext();) {
            String score = itrScr.next();
            String value = scoreHt.get(score);

            Element scoreEle = new Element("search_score");
            scoreEle.setAttribute("name", score);
            scoreEle.setAttribute("value", value);
            peptideEle.addContent(scoreEle);
        }

        return peptideEle;
    }

    private Element createXmlMrmChroPeptideTitle(boolean isDataDependent, Peptide peptide) {
        Element peptideEle = new Element("peptide");
        peptideEle.setAttribute("unique", peptide.isUnique() ? "*" : "");
        peptideEle.setAttribute("file", peptide.getFileName());
//        peptideEle.setAttribute("scan", peptide.getScanNum());
        peptideEle.setAttribute("seq", peptide.getSequence());
        peptideEle.setAttribute("xcorr", peptide.getXCorr());
        peptideEle.setAttribute("deltaCN", peptide.getDeltCN());
        peptideEle.setAttribute("charge", peptide.getChargeState());

        return peptideEle;
    }

    /*
    private String createChroPeptideTitle(boolean isDataDependent, Peptide peptide)
    {
        StringBuffer sb = new StringBuffer();

        sb.append("[PEPTIDE]Unique\tFile Name\tScan Num\tSequence\tXCorr\tDentaCN\tCharge");

        if(isDataDependent)
            sb.append("\n");
        else
            sb.append("\tStartScan\tEndScan\n");
            //sb.append("\tStartScan\tEndScan\tDTAPeakStart\tDTAPeakEnd\n");

        sb.append(peptide.isUnique()?"*":"").append("\t");
        sb.append(peptide.getFileName()).append("\t");
        sb.append(peptide.getScanNum()).append("\t");
        sb.append(peptide.getSequence()).append("\t");
        sb.append(peptide.getXCorr()).append("\t");
        sb.append(peptide.getDeltCN()).append("\t");
        sb.append(peptide.getChargeState());

        if(isDataDependent)
            sb.append("\n");
        else
            sb.append("\t");

        return sb.toString();
    }
     */
    public void createMRMFragmentIons(ChroProgressDialog progress) throws IOException, Exception {
        this.progress = progress;
        int[] keys;

        if (null == conf) {
            conf = Configuration.getInstance();
        }
        this.filePath = conf.getFilePath();

        conf.readSimpleXml(filePath + "census_config.xml");
        conf.setIdFileName(filePath + "DTASelect-filter.txt");

        Hashtable<String, IndexedFile> ht = createIndexedFiles(filePath, CensusConstants.MS2_FILE);

        IndexedFile iFile;
        BufferedOutputStream out = null;
        PrintStream p = null;

        long startTime = System.currentTimeMillis();
        String pepSequence = null;

        try {

            /**
             * ****************************************************************
             * Read DTASelect.txt file to find spectrum range for each peptide
             *****************************************************************
             */
            ChroProgressDialog.addMessageWithLine(progress, "");

            IsotopeReader isoReader = new IsotopeReader(conf.getRootConfEle());

            SpecRangeGenerator rangeGen = null;
            File dtaFile = new File(filePath + "DTASelect.txt");
            IdentificationReader idReader = BaseIdentificationReader.getIdentificationInst(isoReader);

            if (dtaFile.exists()) {
                System.out.print("Parsing DTASelect.txt...");
                ChroProgressDialog.addMessageWithLine(progress, "Parsing DTASelect.txt...");
                rangeGen = SpecRangeGenerator.getSpecRangeGenerator(idReader);

                System.out.println("done.");

            } else {
                rangeGen = new SpecRangeGenerator();
            }

//            ChroProgressDialog.addMessageWithLine(progress, "done");
            TIntLongHashMap index;

            IsotopeTable<String, int[]> isoTable = isoReader.getIsotope();
            //IsotopeTable<String, int[]> isoTable = IsotopeReader.getStandardIsotopeTable();

            int[] sampleNterm = isoTable.get("sampleNTERM");
            int[] sampleCterm = isoTable.get("sampleCTERM");

            //int redundantPeptideNum = idReader.getRedundantPeptideNum();
            int redundantPeptideNum = idReader.getTotalPeptideNumber();
            //increase status bar
            double percent = 0.0;
            double eachSeg = 100.0 / redundantPeptideNum;
            int pepCount = 0;

            IsotopeDist sampleDist;
            Protein protein;
            Peptide peptide;

            Element rootEle = this.createXmlChroHeader(2);

            ElementComposition element;
            ElementComposition totalElement;

            //int[] elementSampleArr;
            int keyIndex;
            int start;
            int last;

            double samplePrecursor;
            double refPrecursor;

            double[][] bionSample;
            double[][] bionRef;
            double[][] yionSample;
            double[][] yionRef;

            Element proteinEle = null;
            Element peptideEle = null;

            for (Iterator<Protein> itr = idReader.getProteins(); itr.hasNext();) {
                protein = itr.next();

                proteinEle = new Element("protein");
                proteinEle.setAttribute("locus", protein.getLocus());
                proteinEle.setAttribute("seq_ct", protein.getSeqCount());
                proteinEle.setAttribute("spec_ct", protein.getSpectrumCount());
                proteinEle.setAttribute("seq_cov", protein.getSeqCoverage());
                proteinEle.setAttribute("length", protein.getLength());
                proteinEle.setAttribute("molwt", protein.getMolWt());
                proteinEle.setAttribute("pi", protein.getPI());
                proteinEle.setAttribute("val", protein.getValidation());

                try {
                    proteinEle.setAttribute("desc", protein.getDescription());
                } catch (org.jdom.IllegalDataException ide) {
                    proteinEle.setAttribute("desc", StringUtil.removeIsoControlChar(protein.getDescription()));
                }

                for (Iterator<Peptide> pepItr = protein.getPeptides(); pepItr.hasNext();) {
                    peptide = pepItr.next();
                    pepCount++;

                    pepSequence = peptide.getSequence();

                    //System.out.println(pepSequence);
                    //System.out.println("");
                    //What is the purpose of hs???
                    char[] ch = pepSequence.substring(2, peptide.getSequence().length() - 2).toCharArray();

                    try {
                        totalElement = new ElementComposition(ch, 0, ch.length, isoTable);
                        totalElement.lightCalculate();

                    } catch (InvalidAAException ive) {
                        System.out.println("Not Quantifiable peptide : " + pepSequence);

                        percent += eachSeg;
                        if (null != progress) {
                            ChroProgressDialog.addMessageWithLine(progress, "Not Quantifiable peptide : " + pepSequence);
                            progress.setProgress((int) percent);
                        }
                        continue;
                    }

                    int chargeState = Integer.parseInt(peptide.getChargeState());

                    IsotopeDist totalDist = new IsotopeDist(totalElement.getElementSampleArr(), totalElement.getModShift(), true);

                    int pepLength = 0;

                    for (int i = 0; i < ch.length; i++) {
                        if (ch[i] == '*' || ch[i] == '@' || ch[i] == '#') {
                            continue;
                        }

                        pepLength++;
                    }

                    bionSample = new double[pepLength][chargeState * 3];
                    //  bionRef = new double[pepLength][chargeState*3];
                    //Yions
                    yionSample = new double[pepLength][chargeState * 3];
                    //yionRef = new double[pepLength][chargeState*3];

                    int pepIndex = 0;

                    //System.out.println("aamass" + massTolerance + " " + conf.getMassTolerance());
                    for (int i = 0; i < ch.length; i++) {
                        if (ch[i] == '*' || ch[i] == '@' || ch[i] == '#') {
                            continue;
                        }

                        element = new ElementComposition(ch, 0, i + 1, isoTable);
                        element.lightCalculate();

                        //Y ions
                        sampleDist = new IsotopeDist(
                                getComplementaryComposition(totalElement.getElementSampleArr(), element.getElementSampleArr(), sampleNterm, sampleCterm), element.getModShift(), true);

                        switch (chargeState) {
                            case 3:
                                yionSample[(pepIndex + 1) % pepLength][8] = (sampleDist.getEndMass() + 3 * PROTON_MASS) / 3 + massTolerance;
                                yionSample[(pepIndex + 1) % pepLength][7] = (sampleDist.getStartMass() + 3 * PROTON_MASS) / 3 - massTolerance;
                                yionSample[(pepIndex + 1) % pepLength][6] = (sampleDist.getAvgMass() + 3 * PROTON_MASS) / 3;

                            case 2:
                                yionSample[(pepIndex + 1) % pepLength][5] = (sampleDist.getEndMass() + 2 * PROTON_MASS) / 2 + massTolerance;
                                yionSample[(pepIndex + 1) % pepLength][4] = (sampleDist.getStartMass() + 2 * PROTON_MASS) / 2 - massTolerance;
                                yionSample[(pepIndex + 1) % pepLength][3] = (sampleDist.getAvgMass() + 2 * PROTON_MASS) / 2;

                            case 1:
                                yionSample[(pepIndex + 1) % pepLength][2] = sampleDist.getEndMass() + 1 * PROTON_MASS + massTolerance; //add proton to give b fragment ion
                                yionSample[(pepIndex + 1) % pepLength][1] = sampleDist.getStartMass() + 1 * PROTON_MASS - massTolerance;
                                yionSample[(pepIndex + 1) % pepLength][0] = sampleDist.getAvgMass() + 1 * PROTON_MASS;

                            default:
                                break;
                        }

                        element.calculateBion();
                        sampleDist = new IsotopeDist(element.getElementSampleArr(), element.getModShift(), true);
                        //refDist = new IsotopeDist(element.getElementRefArr(), false);

                        switch (chargeState) {
                            case 3:
                                bionSample[pepIndex][8] = (sampleDist.getEndMass() + 2 * PROTON_MASS) / 3 + massTolerance;
                                bionSample[pepIndex][7] = (sampleDist.getStartMass() + 2 * PROTON_MASS) / 3 - massTolerance;
                                bionSample[pepIndex][6] = (sampleDist.getAvgMass() + 2 * PROTON_MASS) / 3;

                            case 2:
                                bionSample[pepIndex][5] = (sampleDist.getEndMass() + 1 * PROTON_MASS) / 2 + massTolerance;
                                bionSample[pepIndex][4] = (sampleDist.getStartMass() + 1 * PROTON_MASS) / 2 - massTolerance;
                                bionSample[pepIndex][3] = (sampleDist.getAvgMass() + 1 * PROTON_MASS) / 2;

                            case 1:
                                bionSample[pepIndex][2] = sampleDist.getEndMass() + massTolerance;
                                bionSample[pepIndex][1] = sampleDist.getStartMass() - massTolerance;
                                bionSample[pepIndex][0] = sampleDist.getAvgMass();

                            default:
                                break;

                        }

                        pepIndex++;

                    }

                    element = new ElementComposition(ch, 0, ch.length, isoTable);
                    element.lightCalculate();

                    sampleDist = new IsotopeDist(element.getElementSampleArr(), element.getModShift(), true);

                    switch (chargeState) {
                        case 3:
                            yionSample[0][8] = (sampleDist.getEndMass() + 3 * PROTON_MASS) / 3 + massTolerance;
                            yionSample[0][7] = (sampleDist.getStartMass() + 3 * PROTON_MASS) / 3 - massTolerance;
                            yionSample[0][6] = (sampleDist.getAvgMass() + 3 * PROTON_MASS) / 3;

                        case 2:
                            yionSample[0][5] = (sampleDist.getEndMass() + 2 * PROTON_MASS) / 2 + massTolerance;
                            yionSample[0][4] = (sampleDist.getStartMass() + 2 * PROTON_MASS) / 2 - massTolerance;
                            yionSample[0][3] = (sampleDist.getAvgMass() + 2 * PROTON_MASS) / 2;
                        case 1:
                            yionSample[0][2] = sampleDist.getEndMass() + 1 * PROTON_MASS + massTolerance; //add proton to give b fragment ion
                            yionSample[0][1] = sampleDist.getStartMass() + 1 * PROTON_MASS - massTolerance;
                            yionSample[0][0] = sampleDist.getAvgMass() + 1 * PROTON_MASS;

                        //System.out.println(sampleDist.getAvgMass() + " " + yionSample[0][0] + " " + yionSample[0][1] + " " + yionSample[0][2]);
                        default:
                            break;

                    }

                    peptideEle = this.createXmlChroPeptideTitle(false, peptide);

                    String ms2FileName = peptide.getFileName().substring(0, peptide.getFileName().indexOf("."));

                    StringBuffer rangeKey = new StringBuffer();
                    rangeKey.append(protein.getLocus());
                    rangeKey.append(ms2FileName);
                    rangeKey.append(peptide.getSequence().substring(2, peptide.getSequence().length() - 2));
                    SpecRange range = rangeGen.getSpecRange(rangeKey.toString());

                    if (null == range) {
                        int tmpScanNum = Integer.parseInt(peptide.getScanNum());
                        range = new SpecRange(tmpScanNum, tmpScanNum);
                        peptideEle.setAttribute("start_scan", peptide.getScanNum());
                        peptideEle.setAttribute("end_scan", peptide.getScanNum());
                    } else {
                        peptideEle.setAttribute("start_scan", String.valueOf(range.getMin()));
                        peptideEle.setAttribute("end_scan", String.valueOf(range.getMax()));
                    }

                    peptideEle.setAttribute("MHPlus", String.valueOf(Formatter.formatDecimal(totalDist.getAvgMass() + PROTON_MASS)));
                    peptideEle.setAttribute("CalcMHplus", peptide.getCalcMHplus());

                    Element chro = new Element("chro"); //scan # and intensity
                    //output.append("[CHROMATOGRAMS]\tSCAN\tSAMPLE\tREFERENCE\n");

                    iFile = ht.get(this.filePath + ms2FileName + "." + "ms2");
                    keys = iFile.getKeys();
                    keyIndex = Arrays.binarySearch(keys, Integer.parseInt(peptide.getScanNum()));

                    if (keyIndex < 0) //Cannot find index
                    {
                        keyIndex = -(++keyIndex); //Math.abs(++keyIndex);
                    }
                    //System.out.println("key" + keyIndex);

                    samplePrecursor = 1;
                    refPrecursor = 1;

                    TIntDoubleHashMap precursorMap = iFile.getPrecursorMap();

                    double calcSamMass = (sampleDist.getAvgMass() + chargeState * PROTON_MASS) / chargeState;

                    conf.setCalcSamAvgMass(calcSamMass);
                    peptideEle.setAttribute("lightStartMass", String.valueOf(sampleDist.getStartMass()));
                    peptideEle.setAttribute("lightAvgMass", String.valueOf(conf.getCalcSamAvgMass()));

                    try {

                        double[] d = CalcUtil.readSpectrum(keys, keyIndex, iFile, bionSample, yionSample);
                        StringBuffer sb = new StringBuffer();

                        for (int i = 0; i < 2; i++) {
                            sb.append(d[i]).append(" ");
                        }

                        sb.append(";");
                        for (int i = 2; i < d.length; i++) {
                            if (i == (d.length / 2 + 1)) {
                                sb.append(";");
                            }

                            sb.append(d[i]).append(" ");
                        }
//	System.out.println(d[i]);
                        chro.setText(sb.toString()); //CalcUtil.readSpectrum(keys, keyIndex, iFile, bionSample, yionSample ));

                        //} catch (PrecursorNotFoundException ive)
                    } catch (Exception ive) {
                        System.out.println("Precursor not found for " + pepSequence);
                        ChroProgressDialog.addMessageWithLine(progress, "Error : Precursor not found for " + pepSequence);
                        percent += eachSeg;

                        if (null != progress) {
                            progress.setProgress((int) percent);
                        }

                        continue;
                    }

                    peptideEle.addContent(chro);

                    Element fragEle = new Element("frag");

                    Element bSample = new Element("bs");
                    StringBuffer tempSb = new StringBuffer();
                    for (int i = 0; i < bionSample.length; i++) {
                        int j;
                        for (j = 0; j < bionSample[i].length - 1; j++) {
                            tempSb.append(formatter.format(bionSample[i][j])).append(" ");
                        }

                        tempSb.append(formatter.format(bionSample[i][j])).append(",");
                    }

                    bSample.setText(tempSb.toString());
                    fragEle.addContent(bSample);

//                    Element bRef = new Element("br");
                    tempSb.delete(0, tempSb.length());

                    Element ySample = new Element("ys");
                    tempSb.delete(0, tempSb.length());
                    for (int i = 0; i < yionSample.length; i++) {
                        int j;
                        for (j = 0; j < yionSample[i].length - 1; j++) {
                            tempSb.append(formatter.format(yionSample[i][j])).append(" ");

                        }

                        tempSb.append(formatter.format(yionSample[i][j])).append(",");
                    }

                    ySample.setText(tempSb.toString());
                    fragEle.addContent(ySample);
                    tempSb.delete(0, tempSb.length());

                    peptideEle.addContent(fragEle);

                    percent += eachSeg;

                    if (null != progress) {
                        progress.setProgress((int) percent);
                    }

                    System.out.print(pepCount);
                    System.out.print("/");
                    System.out.print(redundantPeptideNum);
                    System.out.print(" peptides, ");
                    System.out.print((int) percent);
                    //System.out.print(" % is complete\n");
                    System.out.print(" % is complete\r");

                    proteinEle.addContent(peptideEle);

                }

                if (proteinEle.getChildren().size() > 0) {
                    rootEle.addContent(proteinEle);
                }

            }

            Document doc = new Document(rootEle);
            OutputStream os = new FileOutputStream(filePath + "mrm_frags.xml");
            XMLOutputter outputter = new XMLOutputter();
            outputter.setFormat(Format.getPrettyFormat());
            outputter.output(doc, os);
            os.close();
            System.out.println("\n100% complete");
            //System.out.println( System.out.println"\n100% complete");

        } catch (IOException e) {
            System.out.println("IO Error while generating chro file : " + pepSequence + e);
            e.printStackTrace();
            throw new IOException(e.toString());
        } catch (Exception e) {
            System.out.println("Error while generating chro file : " + pepSequence + e);
            e.printStackTrace();
            throw new Exception(e.toString());
        } finally {
            if (null != p) {
                p.close();
            }

            if (null != out) {
                out.close();
            }

            //Close all random files
            for (Enumeration e = ht.keys(); e.hasMoreElements();) {
                iFile = ht.get(e.nextElement());

                if (null != iFile) {
                    iFile.close();
                }
            }

        }
    }

    //iTRAQ data tmt
    public void createMsmsSpecificChro() throws IOException, Exception {
        createMsmsSpecificChroQuick(null);
    }

    public void createMsmsSpecificChro(ChroProgressDialog progress) throws IOException, Exception {

        this.progress = progress;
        int[] keys;

        this.filePath = conf.getFilePath();


        //check if the metabolic heavy labeling or not.  if yes, create soft link with H
        List<String> sqtFilelist = FileFilterUtil.getFilesBySuffix(this.filePath, "sqt");
        for(String each:sqtFilelist) {

            String ms2file = each.substring(0, each.lastIndexOf(".")) + ".ms2";


            File ms2f = new File(this.filePath + File.separator + ms2file);

            boolean isHeavy = FileFilterUtil.isHeavyFile(ms2file, this.filePath);

            //System.out.println(each + " " + ms2file + " " + ms2f.exists() + " " + isHeavy);

            if(!ms2f.exists() && isHeavy) {
                FileFilterUtil.makeSymbolicLink(this.filePath, this.filePath, ms2file.substring(1), ms2file);
                FileFilterUtil.makeSymbolicLink(this.filePath, this.filePath, ms2file.substring(1) + ".index", ms2file + ".index");
            }
        }

        Hashtable<String, IndexedFile> htMS3 = createIndexedFiles(filePath, CensusConstants.MS3_FILE);
        Hashtable<String, IndexedFile> ht = createIndexedFiles(filePath, CensusConstants.MS2_FILE);
        Hashtable<String, IndexedFile> htMs1 = createIndexedFiles(filePath, CensusConstants.MS1_FILE);
        String ms1path = null;

        boolean usingOldPath = false;
        if(htMs1.isEmpty())
        {
            ms1path = filePath.concat("/../../spectra");
            htMs1 = createIndexedFiles(ms1path,CensusConstants.MS1_FILE);
            usingOldPath = !htMs1.isEmpty();
        }

        if (conf.isMs3ScanRandom() || "ms3".equals(conf.getFileShift()) || conf.getScanShift() > 0) {
            ht.putAll(createIndexedFiles(filePath, CensusConstants.MS3_FILE));
        }

        //////////////////// remove this later robin
//	Hashtable<String, IndexedFile> hcdht = createIndexedFiles("/home/rpark/rpark_on_data/project/xmhan/pqd_hcd_comparison/hcd_cid_comparison/CIDHCD-091808OT3/hcd", CensusConstants.MS2_FILE);
        IndexedFile iFile;
        IndexedFile iFileMs1 = null;
        BufferedOutputStream out = null;
        PrintStream p = null;

        long startTime = System.currentTimeMillis();
        String pepSequence = null;
        String tempPath = filePath;
        double spTolerance = 10;

        try {

            /**
             * ****************************************************************
             * Red DTASelect.txt file to find spectrum range for each peptide
             *****************************************************************
             *
             */


            ChroProgressDialog.addMessageWithLine(progress, "");

            IsotopeReader isoReader = null;

            double lkMass =0;
            double lnMass =0;

            double kMass =0;
            double nMass =0;

            double hKMass =0;
            double hNMass =0;

            double mKMass =0;
            double mNMass =0;

            int searchCode =0;

            String searchXMlFile = this.filePath+File.separatorChar+"search.xml";
            String hSearchXMlFile = this.filePath+File.separatorChar+"Hsearch.xml";
            String mSearchXMlFile = this.filePath+File.separatorChar+"Msearch.xml";

            File searchFile = new File(searchXMlFile);
            edu.scripps.pms.mspid.SearchParams sp =null;
            edu.scripps.pms.mspid.SearchParams lSp =null;
            edu.scripps.pms.mspid.SearchParams hSp =null;
            edu.scripps.pms.mspid.SearchParams mSp =null;
            if(searchFile.exists())
                    lSp= new edu.scripps.pms.mspid.SearchParams(searchXMlFile);



            if(lSp!=null)
            {
                for(Iterator<Modification> mItr = lSp.getStaticMods(); mItr.hasNext();)
                {
                    Modification mod = mItr.next();
                    if(mod.getResidue()=='K')
                    {
                        lkMass = mod.getMassShift();
                    }
                }
                 lnMass = lSp.getStaticNTermMod();
                spTolerance = lSp.getPrecursorTolerance();
            }


            if (conf.isXmlConf()) {
                isoReader = new IsotopeReader(conf.getRootConfEle());
            } else {
                isoReader = new IsotopeReader(isotopeFile);
            }

            SpecRangeGenerator rangeGen = null;
            //File dtaFile = new File(filePath + "DTASelect.txt");
            File JSONfile =new File(filePath+MASTER_JSON);
            String xyvaluePath = filePath+MASTER_JSON+File.separatorChar+ XY_VALUES+File.separatorChar;
            File xyValuesDir = new File(filePath+MASTER_JSON+File.separatorChar+ XY_VALUES);
            if(!JSONfile.exists())
            {
                JSONfile.mkdir();
            }
            if(!xyValuesDir.exists())
            {
                xyValuesDir.mkdir();
            }

            IdentificationReader idReader = BaseIdentificationReader.getIdentificationInst(isoReader);
            int redundantPeptideNum = idReader.getTotalPeptideNumber();


            //if (dtaFile.exists()) {
            if (false) {  //we don't use DTASelect.txt parsing any more
                System.out.print("Parsing DTASelect.txt...");
                ChroProgressDialog.addMessageWithLine(progress, "Parsing DTASelect.txt...");
                rangeGen = SpecRangeGenerator.getSpecRangeGenerator(idReader);
                System.out.println("done.");
            } else {
                rangeGen = new SpecRangeGenerator();
            }





            if (conf.isPrintLog()) {
                p = new PrintStream(new FileOutputStream(filePath + "progress.log"));
            }




            //DTASel
            //
            // ectFilterReader dtaReader = new DTASelectFilterReader(conf.getIdFileName());
            //SpecRangeGenerator rangeGen = new SpecRangeGenerator(this.filePath + this.dtaSelectFile, idReader.isVersion2(), idReader.getConfidence());
            Hashtable tempht = rangeGen.getTable();

            ChroProgressDialog.addMessageWithLine(progress, "done");
            TIntLongHashMap index;

            IsotopeTable<String, int[]> isoTable = isoReader.getIsotope();

            //increase status bar
            double percent = 0.0;
            double eachSeg = 100.0 / redundantPeptideNum;
            int pepCount = 0;

            Protein protein;
            Peptide peptide;

            Element rootEle = null;

            if (conf.getMsmsSpectrumNum() == conf.MSMS_SINGLE_SPECTRUM) {
                rootEle = this.createXmlChroHeader(2, CensusConstants.MSMS_SPECIFIC_SINGLE_MASS);
            } else if (conf.getMsmsSpectrumNum() == conf.MSMS_MULTIPLE_SPECTRA) {
                rootEle = this.createXmlChroHeader(2, CensusConstants.MSMS_SPECIFIC_MULTIPLE_MASS);
            } else {
                System.out.println("Cannot find spectral extraction type");
                throw new Exception("Cannot find spectral extraction type");
            }

            ElementComposition element=null;
            ElementComposition totalElement;

            //int[] elementSampleArr;
            int keyIndex;
            int start;
            int last;

            double samplePrecursor;
            double refPrecursor;

//            double[][] bionSample;
//            double[][] bionRef;
//            double[][] yionSample;
//            double[][] yionRef;
            Element proteinEle = null;
            Element peptideEle = null;
            Gson gson = new Gson();
            List<Double> signalList = new ArrayList<>();
            List<Double> purityList = new ArrayList<>();
            double maxSignalToNoise = Double.MIN_VALUE;
            double minSignalTONoise = Double.MAX_VALUE;
            List<Protein> aList = new ArrayList<>();
            for (Iterator<Protein> itr = idReader.getProteins(); itr.hasNext();) {
                protein = itr.next();

                proteinEle = new Element("protein");
                proteinEle.setAttribute("locus", protein.getLocus());
                proteinEle.setAttribute("seq_ct", protein.getSeqCount());
                proteinEle.setAttribute("spec_ct", protein.getSpectrumCount());
                proteinEle.setAttribute("seq_cov", protein.getSeqCoverage());
                proteinEle.setAttribute("length", protein.getLength());
                proteinEle.setAttribute("molwt", protein.getMolWt());
                proteinEle.setAttribute("pi", protein.getPI());
                proteinEle.setAttribute("val", protein.getValidation());
                List<MZValues> mzValuesList =null;
                try {
                    proteinEle.setAttribute("desc", protein.getDescription());
                } catch (org.jdom.IllegalDataException ide) {
                    proteinEle.setAttribute("desc", StringUtil.removeIsoControlChar(protein.getDescription()));
                }
              //  List<TandemTagPeptide> tandemTagPeptideList = new ArrayList<>();
                for (Iterator<Peptide> pepItr = protein.getPeptides(); pepItr.hasNext(); ) {
                    peptide = pepItr.next();
                   // double theorMass = Double.parseDouble(peptide.getCalcMHplus());
                    pepCount++;
               //     TandemTagPeptide tandemPeptide = new TandemTagPeptide();
                    pepSequence = peptide.getSequence();
             /*       tandemPeptide.setSequence(pepSequence);
                    tandemPeptide.setFilename(peptide.getFileName());
                    tandemPeptide.setCstate(peptide.getChargeState());
                    tandemPeptide.setScanNum(peptide.getScanNum());
                    tandemPeptide.setSpc(peptide.getRedundancy());*/
                    //System.out.println(pepSequence);
                    //System.out.println("");
                    //System.out.println("");
                    char[] ch = pepSequence.substring(2, peptide.getSequence().length() - 2).toCharArray();

                    try {
                        totalElement = new ElementComposition(ch, 0, ch.length, isoTable);
                        totalElement.lightCalculate();
                    } catch (InvalidAAException ive) {
                        System.out.println("Not Quantifiable peptide - invalid amino acid : " + pepSequence);

                        percent += eachSeg;
                        if (null != progress) {
                            ChroProgressDialog.addMessageWithLine(progress, "Not Quantifiable peptide - invalid amino acid : " + pepSequence);
                            progress.setProgress((int) percent);
                        }

                        continue;
                    }

                    int chargeState = Integer.parseInt(peptide.getChargeState());

                    int pepLength = 0;

                    for (int i = 0; i < ch.length; i++) {
                        if (ch[i] == '*' || ch[i] == '@' || ch[i] == '#') {
                            continue;
                        }

                        pepLength++;
                    }

                    //bionSample = new double[pepLength][chargeState*3];
                    //bionRef = new double[pepLength][chargeState*3];
                    //Yions
                    //yionSample = new double[pepLength][chargeState*3];
                    //yionRef = new double[pepLength][chargeState*3];
                    int pepIndex = 0;

                    //System.out.println("aamass" + massTolerance + " " + conf.getMassTolerance());
                    for (int i = 0; i < ch.length; i++) {
                        if (ch[i] == '*' || ch[i] == '@' || ch[i] == '#') {
                            continue;
                        }

                        try {
                            element = new ElementComposition(ch, 0, i + 1, isoTable);
                            element.lightCalculate();

                        } catch (InvalidAAException ive) {
                            System.out.println("Not Quantifiable peptide : " + pepSequence);

                            percent += eachSeg;
                            if (null != progress) {
                                ChroProgressDialog.addMessageWithLine(progress, "Not Quantifiable peptide : " + pepSequence);
                                progress.setProgress((int) percent);

                            }
                            continue;
                        }

                        pepIndex++;
                    }
                    DisplayData.DisplayPeptide displayPeptide = new DisplayData.DisplayPeptide(peptide);

                    peptideEle = this.createXmlChroPeptideTitle(false, peptide);

                    //String ms2FileName = peptide.getFileName().substring(0, peptide.getFileName().indexOf("."));
                    String ms2FileName = peptide.getFileName();

                    StringBuffer rangeKey = new StringBuffer();
                    rangeKey.append(protein.getLocus());
                    rangeKey.append(ms2FileName);
                    rangeKey.append(peptide.getSequence().substring(2, peptide.getSequence().length() - 2));
                    SpecRange range = rangeGen.getSpecRange(rangeKey.toString());
                    DisplayData.DisplayChroData chroData = new DisplayData.DisplayChroData();

                    if (null == range) {
                        int tmpScanNum = Integer.parseInt(peptide.getScanNum());
                        range = new SpecRange(tmpScanNum, tmpScanNum);
                        peptideEle.setAttribute("start_scan", peptide.getScanNum());
                        peptideEle.setAttribute("end_scan", peptide.getScanNum());
                    } else {
                        peptideEle.setAttribute("start_scan", String.valueOf(range.getMin()));
                        peptideEle.setAttribute("end_scan", String.valueOf(range.getMax()));
                    }

                    peptideEle.setAttribute("ptmIndex", peptide.getPtmIndex());
                    peptideEle.setAttribute("ptmIndexProtein", peptide.getPtmIndexProtein());


                    //sb.append("\tStartScan\tEndScan\tDTAPeakStart\tDTAPeakEnd\n");
                    Element chro = new Element("chro");
                    //output.append("[CHROMATOGRAMS]\tSCAN\tSAMPLE\tREFERENCE\n");
                    iFile = ht.get(this.filePath + ms2FileName + "." + "ms2");

                    if (null == iFile) {
                        iFile = ht.get(filePath + ms2FileName.substring(1) + "." + "ms2");

                        if (null == iFile) {

                            System.out.println("Error : cannot find the file " + ms2FileName + "." + "ms2");
                            System.exit(0);
                        }
                    }

                    keys = iFile.getKeys();
                    String msFile = "";

                    if (filePath.endsWith("/")) {
                        tempPath = filePath.substring(0, filePath.length() - 1);
                    } else {
                        tempPath = filePath;
                    }
                    if (null != htMs1) {
                        String testName = ms2FileName+".ms2";

                        if(isHeavyFile(testName,this.filePath))
                        {
                            searchCode = 2;
                            if(hSp==null)
                            {
                                if( new File(hSearchXMlFile).exists())
                                {
                                    hSp= new edu.scripps.pms.mspid.SearchParams(hSearchXMlFile);
                                    for(Iterator<Modification> mItr = hSp.getStaticMods(); mItr.hasNext();)
                                    {
                                        Modification mod = mItr.next();
                                        if(mod.getResidue()=='K')
                                        {
                                            hKMass = mod.getMassShift();
                                        }
                                    }
                                    hNMass = hSp.getStaticNTermMod();
                                }
                            }
                            sp = hSp;
                            kMass = hKMass;
                            nMass = hNMass;

                        }
                        else if(isMediumFile(testName,this.filePath))
                        {
                            searchCode = 1;
                            if(mSp==null)
                            {
                                if( new File(mSearchXMlFile).exists())
                                {
                                    mSp= new edu.scripps.pms.mspid.SearchParams(mSearchXMlFile);
                                    for(Iterator<Modification> mItr = mSp.getStaticMods(); mItr.hasNext();)
                                    {
                                        Modification mod = mItr.next();
                                        if(mod.getResidue()=='K')
                                        {
                                            mKMass = mod.getMassShift();
                                        }
                                    }
                                    mNMass = mSp.getStaticNTermMod();
                                }
                            }
                            sp = hSp;
                            kMass = mKMass;
                            nMass = mNMass;
                        }
                        else
                        {
                            searchCode = 0;
                            kMass = lkMass;
                            nMass = lnMass;
                            sp = lSp;
                        }


                     //   System.out.println(">> "+testName);
                        String fname = searchCode>0?
                                ms2FileName.substring(1,ms2FileName.length()) : ms2FileName;
                        msFile = fname + ".ms1";

                        iFileMs1 = htMs1.get(this.filePath + fname + "." + "ms1");
                        //System.out.println(">>looking for "+iFileMs1.getFileName());
                        if (iFileMs1 == null) {
                       //     System.out.println(">> "+this.filePath + "/../../spectra/" + fname + "." + "ms1");
                            iFileMs1 = htMs1.get(this.filePath + "/../../spectra/" + fname + "." + "ms1");
                    //        System.out.println(">>looking for "+iFileMs1.getFileName());
                      //     // msFile = this.filePath + "/../../spectra/" + ms2FileName + "." + "ms1";
                        }
                        else
                        {
                     //       System.out.println(">>looking for "+iFileMs1.getFileName());
                        }
                    }



                    double isolationWindow = conf.getIsobaricIsolationWindow();
                    boolean contaminanted = false;
                    int purityCount = 0;

                    IsotopeDist sampleDist = new IsotopeDist(element.getElementSampleArr(), element.getModShift(), true);
                    double[] samIsoArr = sampleDist.getHighMassList();
                    double[] samIsoIntArr = sampleDist.getRelabun();

                    if(sp!=null)
                    {
                        int kCount = StringUtils.countMatches(pepSequence.substring(1,pepSequence.length()-1), "K");


                        double totalKMass = kCount *kMass;

                        for (int i = 0; i < samIsoArr.length; i++) {
                            samIsoArr[i] = (samIsoArr[i]+totalKMass+nMass + chargeState * PROTON_MASS) / chargeState;
                        }
                        double massTol = spTolerance ;
//                    spReader.init();

                        if (null != iFileMs1) {
                            String path = usingOldPath ? ms1path : tempPath;
                            double slineMass = iFile.getPrecursorMap().get(peptide.getScanNumber());
//MS_READ
                         // int[] arr = iFile.getPrecursorMap().keys();
                         // for(int iii:arr)
                         //   System.out.println("====\t" + iii + " " + peptide.getScanNumber());




                            double totalIntensity = 0;
                            double precursorSum =0;
                            int prevScanFromMs2 = iFile.getIdScanToMs1ScanMap().get(peptide.getScanNumber());
                            TDoubleArrayList massArrList = new TDoubleArrayList();
                            TDoubleArrayList intArrList = new TDoubleArrayList();
                            SpectrumUtil.getMs1ScanPeaksArr(path, msFile, peptide.getScanNumber(),
                                    slineMass - isolationWindow, slineMass + isolationWindow, prevScanFromMs2,massArrList,intArrList);
                            //MS_READ

                            if(intArrList == null)
                            {
                                purityList.add(0.0);
                                peptideEle.setAttribute("tmt_purity", Double.toString(0.0));
                             //   System.out.println("Cannot find in isolation window ");
                            }
                            else
                            {
                                double [] intArr = intArrList.toNativeArray();
                                double [] massArr = massArrList.toNativeArray();
                                for(double d: intArr)
                                {
                                    totalIntensity+=d;
                                }

                                for(int i=0; i<samIsoArr.length; i++)
                                {
                                    double mass = samIsoArr[i];
                                    double tempDiff = mass-slineMass;
                                    //   if(tempDiff>3)System.out.println("..... "+tempDiff + " " + chargeState);
                                    precursorSum+=  CalcUtil.intensitySumForSinglePeakNewMassTolerance(massArr, intArr, mass, massTol);
                                }
                                double purity = precursorSum/totalIntensity;
                       /* if(purity>1 || purity< 0 || Double.isNaN(purity)|| Double.isInfinite(purity))
                        {
                            System.out.println("''''''''Bad result");
                        }*/
                                purityList.add(purity);
                                peptideEle.setAttribute("tmt_purity", Double.toString(purity));
                            }


                        }
                    }
                    else
                    {
                        peptideEle.setAttribute("tmt_purity", Double.toString(-1.0));

                    }







                    //System.out.println("====" + peptide.getScanNum() + " " + peptide.getScanNumber());
                   /*
                    if (null != iFileMs1) {
                        int cs = Integer.parseInt(peptide.getChargeState());
                        double d = (Double.parseDouble(peptide.getMhPlus()) - cs - 1) / cs;
                        //System.out.println(d + "  ------------ " + peptide.getMhPlus());
                        //System.out.println(d + "  ------------ " + peptide.getMhPlus());
//                        System.out.println("We are in the  Chrogenerator...............");
                        //contaminanted = CalcUtilGeneric.isIsolateContaminanted(iFileMs1, peptide.getScanNumber(), 1.0, d, cs, peptide.getSequence());
                        // String ms2File = ms2FileName+".ms2";
                        double mass = iFile.getPrecursorMap().get(peptide.getScanNumber());
                        String path = usingOldPath ? ms1path : tempPath;
                        int prevScanFromMs2 = iFile.getIdScanToMs1ScanMap().get(peptide.getScanNumber());
                        List<XYPoint> points = SpectrumUtil.getMs1ScanPeaksArr(path, msFile, peptide.getScanNumber(), mass - mass_tol, mass + mass_tol, prevScanFromMs2);
                        purityCount = points.size();
                        peptideEle.setAttribute("tmt_purity", Integer.toString(purityCount));
                        if (!purityFrequencyMap.adjustValue(purityCount, 1))
                            purityFrequencyMap.put(purityCount, 1);

                    }*/

                    keyIndex = Arrays.binarySearch(keys, Integer.parseInt(peptide.getScanNum()));
                    if (keyIndex < 0) //Cannot find index
                    {
                        keyIndex = -(++keyIndex); //Math.abs(++keyIndex);
                    }

                    if (conf.isMs3ScanRandom()) {  //ms3 based search
                        String ms3FileName = this.filePath + ms2FileName + ".ms3";
                        IndexedFile ms3Ifile = ht.get(ms3FileName);

                        if (null == ms3Ifile) {  //for heavy file
                            ms3Ifile = ht.get(filePath + ms2FileName.substring(1) + "." + "ms3");

                            if (null == ms3Ifile) {  //for replaced ext

                                ms3Ifile = ht.get(ms3FileName.replace(".ms3", "_ms3.ms2"));
                                if (null == ms3Ifile) {
                                        ms3Ifile = ht.get(ms3FileName.replace(".ms3", ".ms2"));
                                }

                                if (null == ms3Ifile) {

                                    //check if ms3 file can be heavy one
                                    if (ms2FileName.startsWith("H") || ms2FileName.startsWith("M")) {
                                        String heavyMs3File = this.filePath + ms2FileName.substring(1) + "_ms3.ms2";
                                        ms3Ifile = ht.get(heavyMs3File);
                                    }
                                    if (null == ms3Ifile) {

/*
                                        //check one more time, if it is heavy search or not.
                                        boolean isHeavy = FileFilterUtil.isHeavyFile(ms2FileName + ".ms2", this.filePath);
                                        if(isHeavy) {
                                            ms3Ifile = ht.get(filePath + ms2FileName.substring(1) + "." + "ms2");
                                        }

                                        if(null == iFile) {
                                        } */

                                        System.out.println("Error : cannot find the file " + ms2FileName + "." + "ms3");
                                        System.exit(0);
                                    }
                                }
                            }
                        }

                        int ms3Scan = TMTUtil.getCorrespondingMs3Scan(ms3Ifile, keys[keyIndex]);

                        if (ms3Scan < 0) {
                            continue;
                        }

                        iFile = ms3Ifile;
                        //keyIndex

                        keyIndex = Arrays.binarySearch(iFile.getKeys(), ms3Scan);
                        if (keyIndex < 0) //Cannot find index
                        {
                            keyIndex = -(++keyIndex); //Math.abs(++keyIndex);
                        }
                    } else if (conf.getScanShift() >= 0) {  //ms2 scan shift
                        keyIndex += conf.getScanShift();
                    }// else {  //regular search


                    try {
//MS_READ
                        String outStr = CalcUtil.calculateMS2Mass(iFile, range, keyIndex, null, null, null, null, conf, chargeState);
                        peptideEle.setAttribute("signal-noise", CensusHelper.format.format(CalcUtil.SignalToNoise));
                        peptideEle.setAttribute("report-ion-signal-noise", CensusHelper.format.format(CalcUtil.ReportIonSignalToNoise));

                        peptideEle.setAttribute("noise_level", CensusHelper.format.format(CalcUtil.NoiseLevel));

                        signalList.add(CalcUtil.SignalToNoise);
                        maxSignalToNoise = CalcUtil.SignalToNoise> maxSignalToNoise ? CalcUtil.SignalToNoise : maxSignalToNoise;
                        minSignalTONoise = CalcUtil.SignalToNoise< minSignalTONoise ? CalcUtil.SignalToNoise : minSignalTONoise;
                  //      mzValuesList = CalcUtil.MZValuesList;
                    //    tandemPeptide.setMzValues(mzValuesList);
                        if (null == outStr) {
                            chro = null;
                        } else {
                            chro.setText(outStr);
                        }
                    } catch (PrecursorNotFoundException ive) {
                        System.out.println("Precursor not found for " + pepSequence);
                        ChroProgressDialog.addMessageWithLine(progress, "Error : Precursor not found for " + pepSequence);
                        percent += eachSeg;
                        progress.setProgress((int) percent);

                        continue;
                    }

                    if (null != chro) {
                        peptideEle.addContent(chro);
                    }

                    percent += eachSeg;

                    if (null != progress) {
                        progress.setProgress((int) percent);
                    }
                    int ms3Scan;
                    // if(!ms3ScanMap.containsKey(fileName)) {


                    List<ReportIon> reportIonList = conf.getReportIonList();
                    double startmass = reportIonList.get(0).getMass();
                    double endMass = reportIonList.get(reportIonList.size()-1).getMass();



                    IndexedFile iFile2 = ht.get(filePath+ peptide.getFileName() + "." + "ms2");
                    if(null == iFile2) {

                        if(FileFilterUtil.isHeavyFile(peptide.getFileName() + ".ms2", filePath)) {
                            iFile2 = ht.get(filePath + peptide.getFileName().substring(1) + "." + "ms2");
                            if(null != iFile2)
                                peptide.setFileName(peptide.getFileNameWithScan().substring(1));
                        }
                        else if(null == iFile2) {
                            iFile2 = ht.get(filePath+ peptide.getFileName() + "." + "ms3");

                            if(null == iFile2 && FileFilterUtil.isHeavyFile(peptide.getFileName() + ".ms3", filePath)) {
                                iFile2 = ht.get(filePath + peptide.getFileName().substring(1) + "." + "ms3");
                                if(null != iFile2)
                                    peptide.setFileName(peptide.getFileNameWithScan().substring(1));

                            }
                        }



                    }
                    //check if it is heavy search or not.
                    /*if(null == iFile2) {
                        boolean isHeavy = FileFilterUtil.isHeavyFile(ms2FileName + ".ms2", this.filePath);
                        if (isHeavy) {
                            iFile2 = ht.get(filePath + ms2FileName.substring(1) + "." + "ms2");
                         //   peptide.setFileName(ms2FileName.substring(1));
                        }
                    } */


                    int [] keys2 = iFile2.getKeys();

                    //  long start = System.currentTimeMillis();

                    String ms3FilePattern="orig";
                    List<edu.scripps.pms.census.model.XYPoint> xyPoints;
                    if (conf.isMs3ScanRandom()) {

                        String ms3FileName = filePath + peptide.getFileName() + ".ms3";

                        IndexedFile ms3Ifile = ht.get(ms3FileName);




                        if (null == ms3Ifile) { // for heavy file
                            ms3Ifile = ht.get(filePath+ peptide.getFileName().substring(1) + "." + "ms3");

                            if (null == ms3Ifile) { // for replaced ext

                                ms3Ifile = ht.get(ms3FileName.replace(".ms3", "_ms3.ms2"));
                                ms3FilePattern="_ms3.ms2";

                                if (null == ms3Ifile) {

                                    if (peptide.getFileName().startsWith("H")
                                            || peptide.getFileName().startsWith("M")) {
                                        String heavyMs3File =filePath + peptide.getFileName().substring(1)
                                                + "_ms3.ms2";
                                        ms3Ifile = ht.get(heavyMs3File);
                                    }
                                    if (null == ms3Ifile) {
                                        System.out.println("Error : cannot find the file " + peptide.getFileName()
                                                + "." + "ms3");
                                    }
                                }
                            }
                        }


                        int keyIndex2 = Arrays.binarySearch(keys2, Integer.parseInt(peptide.getScanNum()));

                        ms3Scan = TMTUtil.getCorrespondingMs3Scan(ms3Ifile, keys2[keyIndex2]);

                        if(ms3FilePattern.equals("orig")) {
                            //MS_READ

                            File f = new File(filePath + peptide.getFileName() + ".ms3");
			                if(f.exists())
				                xyPoints = SpectrumUtil.getSpectrum(filePath, peptide.getFileName() + ".ms3", ms3Scan,
						    startmass-1, endMass+1 );
			                else
				                xyPoints = SpectrumUtil.getSpectrum(filePath, peptide.getFileName().substring(1) + ".ms3", ms3Scan,
						    startmass-1, endMass+1 );
                        }
                        else {
                            //MS_READ

                            File f = new File(filePath + peptide.getFileName() + "_ms3.ms2");
                            if(f.exists())
                                xyPoints = SpectrumUtil.getSpectrum(filePath, peptide.getFileName() + "_ms3.ms2", ms3Scan,
                                        startmass-1, endMass+1 );
                            else
                                xyPoints = SpectrumUtil.getSpectrum(filePath, peptide.getFileName().substring(1) + "_ms3.ms2", ms3Scan,
                                        startmass-1, endMass+1 );
                        }
                    }
                    else{
                        //MS_READ

                        xyPoints = SpectrumUtil.getSpectrum(filePath, peptide.getFileName() + ".ms2",
                                Integer.parseInt(peptide.getScanNum()),startmass -1 , endMass + 1);
                    }
                    String pepKey = peptide.getSequence() + peptide.getChargeState() + peptide.getScanNum();
                    FileUtil.writeJSON(gson.toJson(xyPoints), xyvaluePath +  pepKey + "_xy.JSON");


                    System.out.print(pepCount);
                    System.out.print("/");
                    System.out.print(redundantPeptideNum);
                    System.out.print(" peptides, ");
                    System.out.print((int) percent);
                    //System.out.print(" % is complete\n");
                    System.out.print(" % is complete\r");



                    if (conf.isPrintLog()) {
                        p.println((int) percent + "\t" + pepCount + "\t" + redundantPeptideNum);
                    }



                    //   tandemTagPeptideList.add(tandemPeptide);
                    if (null != chro) {
                        proteinEle.addContent(peptideEle);
                    }
                }

                if (proteinEle.getChildren().size() > 0) {

                    Element redunEle = new Element("redundant");

                    int rPepCount = 0;

                    for (Iterator<Protein> rItr = aList.iterator(); rItr.hasNext();) {
                        rPepCount++;
                        Protein rp = rItr.next();

                        Element rpEle = new Element("protein");
                        rpEle.setAttribute("locus", rp.getLocus());
                        rpEle.setAttribute("seq_ct", rp.getSeqCount());
                        rpEle.setAttribute("spec_ct", rp.getSpectrumCount());
                        rpEle.setAttribute("seq_cov", rp.getSeqCoverage());
                        rpEle.setAttribute("length", rp.getLength());
                        rpEle.setAttribute("molwt", rp.getMolWt());
                        rpEle.setAttribute("pi", rp.getPI());
                        rpEle.setAttribute("val", rp.getValidation());

                        try {
                            rpEle.setAttribute("desc", rp.getDescription());
                        } catch (org.jdom.IllegalDataException ide) {
                            rpEle.setAttribute("desc", StringUtil.removeIsoControlChar(rp.getDescription()));
                        }

                        redunEle.addContent(rpEle);
                    }

                    if (rPepCount > 0) {
                        proteinEle.addContent(redunEle);
                    }

                    rootEle.addContent(proteinEle);
                    aList.clear();
                } else {
                    aList.add(protein);
                    //	    set.add(protein.getLocus());
                }

                //FileUtil.writeJSON(gson.toJson(tandemTagPeptideList), filePath +  File.separator + MASTER_JSON + File.separator +protein.getLocus() + ".JSON");
            }
            List<XYPoint> purityFrequencyPointList = new ArrayList<>();

            Histogram purityHistogram = new Histogram(11,0,1.1);
            purityHistogram.loadData(purityList);
            double [] purityBins = purityHistogram.getBins();
            int[] purityFrequency = purityHistogram.getFreqArr();
            for(int i=0; i<purityBins.length; i++)
            {
                purityFrequencyPointList.add(new XYPoint(purityBins[i],purityFrequency[i],Double.toString(purityBins[i])));
            }



            List<XYPoint> signalHistPointList = new ArrayList<>();
            //Histogram signalHist= new Histogram(20,minSignalTONoise,maxSignalToNoise);
            Histogram signalHist= new Histogram(20,0,200);  //max value is 200
            signalHist.loadData(signalList);
            double [] signalBins = signalHist.getBins();
            int[] signalFreq = signalHist.getFreqArr();
            for(int i=0; i<signalBins.length; i++)
            {
                signalHistPointList.add(new XYPoint(signalBins[i],signalFreq[i],Double.toString(signalBins[i])));
            }
            FileUtil.writeJSON(gson.toJson(purityFrequencyPointList), filePath+PURITY_GRAPH_PATH);
            FileUtil.writeJSON(gson.toJson(signalHistPointList), filePath+SIGNAL_GRAPH_PATH);




            Document doc = new Document(rootEle);
            OutputStream os = new FileOutputStream(filePath + "census_chro.xml");
            XMLOutputter outputter = new XMLOutputter();
            outputter.setFormat(Format.getPrettyFormat());
            outputter.output(doc, os);
            os.close();

            System.out.println("\n100% complete");
            //System.out.println( System.out.println"\n100% complete");

        } catch (IOException e) {
            System.out.println("IO Error while generating msms chro file : " + pepSequence + e);
            e.printStackTrace();
            throw new IOException(e.toString());
        } catch (java.lang.IndexOutOfBoundsException e) {
            System.out.println("Error while generating msms chro file : " + pepSequence + e);
            e.printStackTrace();
            throw new Exception(e.toString());
        } catch (Exception e) {
            System.out.println("Error while generating msms chro file : " + pepSequence + e);
            e.printStackTrace();
            throw new Exception(e.toString());
        } finally {
            if (null != p) {
                p.close();
            }

            if (null != out) {
                out.close();
            }

            //Close all random files
            for (Enumeration e = ht.keys(); e.hasMoreElements();) {
                iFile = ht.get(e.nextElement());

                if (null != iFile) {
                    iFile.close();
                }
            }

        }

    }




    public void createMsmsSpecificChroQuick(ChroProgressDialog progress) throws IOException, Exception {

        this.progress = progress;
        int[] keys;

        this.filePath = conf.getFilePath();


        //check if the metabolic heavy labeling or not.  if yes, create soft link with H
        List<String> sqtFilelist = FileFilterUtil.getFilesBySuffix(this.filePath, "sqt");
        for(String each:sqtFilelist) {

            String ms2file = each.substring(0, each.lastIndexOf(".")) + ".ms2";


            File ms2f = new File(this.filePath + File.separator + ms2file);

            boolean isHeavy = FileFilterUtil.isHeavyFile(ms2file, this.filePath);

            if(!ms2f.exists() && isHeavy) {
                FileFilterUtil.makeSymbolicLink(this.filePath, this.filePath, ms2file.substring(1), ms2file);
                FileFilterUtil.makeSymbolicLink(this.filePath, this.filePath, ms2file.substring(1)
                        + ".index", ms2file + ".index");
            }
        }

       // createSpectraDB(filePath,CensusConstants.MS3_FILE);
        Map<String, SpectraDB> ms3Map = connectCreateSpectraDB(filePath, CensusConstants.MS3_FILE);

    //    createSpectraDB(filePath,CensusConstants.MS2_FILE);
        Map<String, IndexedFile> ms2Ht = createIndexedFiles(filePath, CensusConstants.MS2_FILE);
        Map<String, SpectraDB> ms2Map = connectCreateSpectraDB(filePath, CensusConstants.MS2_FILE);

       // createSpectraDB(filePath,CensusConstants.MS1_FILE);
        Map<String, SpectraDB> ms1Map = connectCreateSpectraDB(filePath, CensusConstants.MS1_FILE);

        String ms1path = null;

     //   boolean usingOldPath = false;
        if(ms1Map.isEmpty())
        {
         //   ms1path = filePath.concat("/../../spectra");
       //     ms1Map = connectCreateSpectraDB(ms1path,CensusConstants.MS1_FILE);
   //         usingOldPath = !ms1Map.isEmpty();
        }

        if (conf.isMs3ScanRandom() || "ms3".equals(conf.getFileShift()) || conf.getScanShift() > 0) {
            ms2Ht.putAll(createIndexedFilesNoMs(filePath, CensusConstants.MS3_FILE));
            ms2Map.putAll(connectCreateSpectraDB(filePath, CensusConstants.MS3_FILE));
        }

        //////////////////// remove this later robin
//	Hashtable<String, IndexedFile> hcdht = createIndexedFiles("/home/rpark/rpark_on_data/project/xmhan/pqd_hcd_comparison/hcd_cid_comparison/CIDHCD-091808OT3/hcd", CensusConstants.MS2_FILE);
        IndexedFile iFile;
        SpectraDB iFileMs1 = null;
        BufferedOutputStream out = null;
        PrintStream p = null;

        long startTime = System.currentTimeMillis();
        String pepSequence = null;
        String tempPath = filePath;
        double spTolerance = 10;

        try {

            /**
             * ****************************************************************
             * Red DTASelect.txt file to find spectrum range for each peptide
             *****************************************************************
             *
             */


            ChroProgressDialog.addMessageWithLine(progress, "");

            IsotopeReader isoReader = null;

            double lkMass =0;
            double lnMass =0;

            double kMass =0;
            double nMass =0;

            double hKMass =0;
            double hNMass =0;

            double mKMass =0;
            double mNMass =0;

            int searchCode =0;

            String searchXMlFile = this.filePath+File.separatorChar+"search.xml";
            String hSearchXMlFile = this.filePath+File.separatorChar+"Hsearch.xml";
            String mSearchXMlFile = this.filePath+File.separatorChar+"Msearch.xml";

            File searchFile = new File(searchXMlFile);
            edu.scripps.pms.mspid.SearchParams sp =null;
            edu.scripps.pms.mspid.SearchParams lSp =null;
            edu.scripps.pms.mspid.SearchParams hSp =null;
            edu.scripps.pms.mspid.SearchParams mSp =null;
            if(searchFile.exists())
                lSp= new edu.scripps.pms.mspid.SearchParams(searchXMlFile);



            if(lSp!=null)
            {
                for(Iterator<Modification> mItr = lSp.getStaticMods(); mItr.hasNext();)
                {
                    Modification mod = mItr.next();
                    if(mod.getResidue()=='K')
                    {
                        lkMass = mod.getMassShift();
                    }
                }
                lnMass = lSp.getStaticNTermMod();
                spTolerance = lSp.getPrecursorTolerance();
            }


            if (conf.isXmlConf()) {
                isoReader = new IsotopeReader(conf.getRootConfEle());
            } else {
                isoReader = new IsotopeReader(isotopeFile);
            }

            SpecRangeGenerator rangeGen = null;
            //File dtaFile = new File(filePath + "DTASelect.txt");
            File JSONfile =new File(filePath+MASTER_JSON);
            String xyvaluePath = filePath+MASTER_JSON+File.separatorChar+ XY_VALUES+File.separatorChar;
            File xyValuesDir = new File(filePath+MASTER_JSON+File.separatorChar+ XY_VALUES);
            if(!JSONfile.exists())
            {
                JSONfile.mkdir();
            }
            if(!xyValuesDir.exists())
            {
                xyValuesDir.mkdir();
            }

            IdentificationReader idReader = BaseIdentificationReader.getIdentificationInst(isoReader);
            int redundantPeptideNum = idReader.getTotalPeptideNumber();
            System.out.println(">><<<>> "+redundantPeptideNum);

            //if (dtaFile.exists()) {
            if (false) {  //we don't use DTASelect.txt parsing any more
                System.out.print("Parsing DTASelect.txt...");
                ChroProgressDialog.addMessageWithLine(progress, "Parsing DTASelect.txt...");
                rangeGen = SpecRangeGenerator.getSpecRangeGenerator(idReader);
                System.out.println("done.");
            } else {
                rangeGen = new SpecRangeGenerator();
            }





            if (conf.isPrintLog()) {
                p = new PrintStream(new FileOutputStream(filePath + "progress.log"));
            }




            //DTASel
            //
            // ectFilterReader dtaReader = new DTASelectFilterReader(conf.getIdFileName());
            //SpecRangeGenerator rangeGen = new SpecRangeGenerator(this.filePath + this.dtaSelectFile, idReader.isVersion2(), idReader.getConfidence());
            Hashtable tempht = rangeGen.getTable();

            ChroProgressDialog.addMessageWithLine(progress, "done");
            TIntLongHashMap index;

            IsotopeTable<String, int[]> isoTable = isoReader.getIsotope();

            //increase status bar
            double percent = 0.0;
            double eachSeg = 100.0 / redundantPeptideNum;
            int pepCount = 0;

            Protein protein;
            Peptide peptide;

            Element rootEle = null;

            if (conf.getMsmsSpectrumNum() == conf.MSMS_SINGLE_SPECTRUM) {
                rootEle = this.createXmlChroHeader(2, CensusConstants.MSMS_SPECIFIC_SINGLE_MASS);
            } else if (conf.getMsmsSpectrumNum() == conf.MSMS_MULTIPLE_SPECTRA) {
                rootEle = this.createXmlChroHeader(2, CensusConstants.MSMS_SPECIFIC_MULTIPLE_MASS);
            } else {
                System.out.println("Cannot find spectral extraction type");
                throw new Exception("Cannot find spectral extraction type");
            }

            ElementComposition element=null;
            ElementComposition totalElement;

            //int[] elementSampleArr;
            int keyIndex;
            int start;
            int last;

            double samplePrecursor;
            double refPrecursor;

//            double[][] bionSample;
//            double[][] bionRef;
//            double[][] yionSample;
//            double[][] yionRef;
            Element proteinEle = null;
            Element peptideEle = null;
            Gson gson = new Gson();
            List<Double> signalList = new ArrayList<>();
            List<Double> purityList = new ArrayList<>();
            double maxSignalToNoise = Double.MIN_VALUE;
            double minSignalTONoise = Double.MAX_VALUE;
            List<Protein> aList = new ArrayList<>();
            SpectraDB ispectraDB = null;

            for (Iterator<Protein> itr = idReader.getProteins(); itr.hasNext();) {
                protein = itr.next();

                proteinEle = new Element("protein");
                proteinEle.setAttribute("locus", protein.getLocus());
                proteinEle.setAttribute("seq_ct", protein.getSeqCount());
                proteinEle.setAttribute("spec_ct", protein.getSpectrumCount());
                proteinEle.setAttribute("seq_cov", protein.getSeqCoverage());
                proteinEle.setAttribute("length", protein.getLength());
                proteinEle.setAttribute("molwt", protein.getMolWt());
                proteinEle.setAttribute("pi", protein.getPI());
                proteinEle.setAttribute("val", protein.getValidation());
                List<MZValues> mzValuesList =null;
                try {
                    proteinEle.setAttribute("desc", protein.getDescription());
                } catch (org.jdom.IllegalDataException ide) {
                    proteinEle.setAttribute("desc", StringUtil.removeIsoControlChar(protein.getDescription()));
                }
                //  List<TandemTagPeptide> tandemTagPeptideList = new ArrayList<>();
                for (Iterator<Peptide> pepItr = protein.getPeptides(); pepItr.hasNext(); ) {
                    peptide = pepItr.next();
                    // double theorMass = Double.parseDouble(peptide.getCalcMHplus());
                    pepCount++;
                    //     TandemTagPeptide tandemPeptide = new TandemTagPeptide();
                    pepSequence = peptide.getSequence();
             /*       tandemPeptide.setSequence(pepSequence);
                    tandemPeptide.setFilename(peptide.getFileName());
                    tandemPeptide.setCstate(peptide.getChargeState());
                    tandemPeptide.setScanNum(peptide.getScanNum());
                    tandemPeptide.setSpc(peptide.getRedundancy());*/
                    //System.out.println(pepSequence);
                    //System.out.println("");
                    //System.out.println("");
                    char[] ch = pepSequence.substring(2, peptide.getSequence().length() - 2).toCharArray();

                    try {
                        totalElement = new ElementComposition(ch, 0, ch.length, isoTable);
                        totalElement.lightCalculate();
                    } catch (InvalidAAException ive) {
                        System.out.println("Not Quantifiable peptide - invalid amino acid : " + pepSequence);

                        percent += eachSeg;
                        if (null != progress) {
                            ChroProgressDialog.addMessageWithLine(progress, "Not Quantifiable peptide - invalid amino acid : " + pepSequence);
                            progress.setProgress((int) percent);
                        }

                        continue;
                    }

                    int chargeState = Integer.parseInt(peptide.getChargeState());

                    int pepLength = 0;

                    for (int i = 0; i < ch.length; i++) {
                        if (ch[i] == '*' || ch[i] == '@' || ch[i] == '#') {
                            continue;
                        }

                        pepLength++;
                    }

                    //bionSample = new double[pepLength][chargeState*3];
                    //bionRef = new double[pepLength][chargeState*3];
                    //Yions
                    //yionSample = new double[pepLength][chargeState*3];
                    //yionRef = new double[pepLength][chargeState*3];
                    int pepIndex = 0;

                    //System.out.println("aamass" + massTolerance + " " + conf.getMassTolerance());
                    for (int i = 0; i < ch.length; i++) {
                        if (ch[i] == '*' || ch[i] == '@' || ch[i] == '#') {
                            continue;
                        }

                        try {
                            element = new ElementComposition(ch, 0, i + 1, isoTable);
                            element.lightCalculate();

                        } catch (InvalidAAException ive) {
                            System.out.println("Not Quantifiable peptide : " + pepSequence);

                            percent += eachSeg;
                            if (null != progress) {
                                ChroProgressDialog.addMessageWithLine(progress, "Not Quantifiable peptide : " + pepSequence);
                                progress.setProgress((int) percent);

                            }
                            continue;
                        }

                        pepIndex++;
                    }
                    DisplayData.DisplayPeptide displayPeptide = new DisplayData.DisplayPeptide(peptide);

                    peptideEle = this.createXmlChroPeptideTitle(false, peptide);

                    //String ms2FileName = peptide.getFileName().substring(0, peptide.getFileName().indexOf("."));
                    String ms2FileName = peptide.getFileName();

                    StringBuffer rangeKey = new StringBuffer();
                    rangeKey.append(protein.getLocus());
                    rangeKey.append(ms2FileName);
                    rangeKey.append(peptide.getSequence().substring(2, peptide.getSequence().length() - 2));
                    SpecRange range = rangeGen.getSpecRange(rangeKey.toString());
                    DisplayData.DisplayChroData chroData = new DisplayData.DisplayChroData();

                    if (null == range) {
                        int tmpScanNum = Integer.parseInt(peptide.getScanNum());
                        range = new SpecRange(tmpScanNum, tmpScanNum);
                        peptideEle.setAttribute("start_scan", peptide.getScanNum());
                        peptideEle.setAttribute("end_scan", peptide.getScanNum());
                    } else {
                        peptideEle.setAttribute("start_scan", String.valueOf(range.getMin()));
                        peptideEle.setAttribute("end_scan", String.valueOf(range.getMax()));
                    }

                    peptideEle.setAttribute("ptmIndex", peptide.getPtmIndex());
                    peptideEle.setAttribute("ptmIndexProtein", peptide.getPtmIndexProtein());


                    //sb.append("\tStartScan\tEndScan\tDTAPeakStart\tDTAPeakEnd\n");
                    Element chro = new Element("chro");
                    //output.append("[CHROMATOGRAMS]\tSCAN\tSAMPLE\tREFERENCE\n");
                    iFile = ms2Ht.get(this.filePath + ms2FileName + "." + "ms2");
                    ispectraDB = iFile.getSpectraDB();
                    //spectraDBMs3 = ms2Map.get(ms2FileName + "." + "ms2");
                    if (null == iFile) {
                        iFile = ms2Ht.get(filePath + ms2FileName.substring(1) + "." + "ms2");
                        ispectraDB =  iFile.getSpectraDB();

                        if (null == iFile) {

                            System.out.println("Error : cannot find the file " + ms2FileName + "." + "ms2");
                            System.exit(0);
                        }
                    }

                    keys = iFile.getKeys();
                    String msFile = "";

                    if (filePath.endsWith("/")) {
                        tempPath = filePath.substring(0, filePath.length() - 1);
                    } else {
                        tempPath = filePath;
                    }
                    if (null != ms1Map) {
                        String testName = ms2FileName+".ms2";

                        if(isHeavyFile(testName,this.filePath))
                        {
                            searchCode = 2;
                            if(hSp==null)
                            {
                                if( new File(hSearchXMlFile).exists())
                                {
                                    hSp= new edu.scripps.pms.mspid.SearchParams(hSearchXMlFile);
                                    for(Iterator<Modification> mItr = hSp.getStaticMods(); mItr.hasNext();)
                                    {
                                        Modification mod = mItr.next();
                                        if(mod.getResidue()=='K')
                                        {
                                            hKMass = mod.getMassShift();
                                        }
                                    }
                                    hNMass = hSp.getStaticNTermMod();
                                }
                            }
                            sp = hSp;
                            kMass = hKMass;
                            nMass = hNMass;

                        }
                        else if(isMediumFile(testName,this.filePath))
                        {
                            searchCode = 1;
                            if(mSp==null)
                            {
                                if( new File(mSearchXMlFile).exists())
                                {
                                    mSp= new edu.scripps.pms.mspid.SearchParams(mSearchXMlFile);
                                    for(Iterator<Modification> mItr = mSp.getStaticMods(); mItr.hasNext();)
                                    {
                                        Modification mod = mItr.next();
                                        if(mod.getResidue()=='K')
                                        {
                                            mKMass = mod.getMassShift();
                                        }
                                    }
                                    mNMass = mSp.getStaticNTermMod();
                                }
                            }
                            sp = hSp;
                            kMass = mKMass;
                            nMass = mNMass;
                        }
                        else
                        {
                            searchCode = 0;
                            kMass = lkMass;
                            nMass = lnMass;
                            sp = lSp;
                        }


                        //   System.out.println(">> "+testName);
                        String fname = searchCode>0?
                                ms2FileName.substring(1,ms2FileName.length()) : ms2FileName;
                        msFile = fname + ".ms1";

                        iFileMs1 = ms1Map.get( fname + "." + "ms1");
                        //System.out.println(">>looking for "+iFileMs1.getFileName());
                        if (iFileMs1 == null) {
                            //     System.out.println(">> "+this.filePath + "/../../spectra/" + fname + "." + "ms1");
                            iFileMs1 = ms1Map.get(fname + "." + "ms1");
                            //        System.out.println(">>looking for "+iFileMs1.getFileName());
                            //     // msFile = this.filePath + "/../../spectra/" + ms2FileName + "." + "ms1";
                        }
                        else
                        {
                            //       System.out.println(">>looking for "+iFileMs1.getFileName());
                        }
                    }



                    double isolationWindow = conf.getIsobaricIsolationWindow();
                    boolean contaminanted = false;
                    int purityCount = 0;

                    IsotopeDist sampleDist = new IsotopeDist(element.getElementSampleArr(), element.getModShift(), true);
                    double[] samIsoArr = sampleDist.getHighMassList();
                    double[] samIsoIntArr = sampleDist.getRelabun();

                    if(sp!=null)
                    {
                        int kCount = StringUtils.countMatches(pepSequence.substring(1,pepSequence.length()-1), "K");


                        double totalKMass = kCount *kMass;

                        for (int i = 0; i < samIsoArr.length; i++) {
                            samIsoArr[i] = (samIsoArr[i]+totalKMass+nMass + chargeState * PROTON_MASS) / chargeState;
                        }
                        double massTol = spTolerance ;
//                    spReader.init();

                        if (null != iFileMs1) {
                          //  String path = usingOldPath ? ms1path : tempPath;
                            double slineMass = iFile.getPrecursorMap().get(peptide.getScanNumber());
//MS_READ
                            // int[] arr = iFile.getPrecursorMap().keys();
                            // for(int iii:arr)
                            //   System.out.println("====\t" + iii + " " + peptide.getScanNumber());




                            double totalIntensity = 0;
                            double precursorSum =0;
                            int prevScanFromMs2 = iFile.getIdScanToMs1ScanMap().get(peptide.getScanNumber());
                            TDoubleArrayList massArrList = new TDoubleArrayList();
                            TDoubleArrayList intArrList = new TDoubleArrayList();
                            SpectraDBCalcUtil.getSpectrumMS1Arr(iFileMs1, peptide.getScanNumber(),
                                    slineMass - isolationWindow, slineMass + isolationWindow, prevScanFromMs2,massArrList,intArrList);
                            //MS_READ

                            if(intArrList == null)
                            {
                                purityList.add(0.0);
                                peptideEle.setAttribute("tmt_purity", Double.toString(0.0));
                                //   System.out.println("Cannot find in isolation window ");
                            }
                            else
                            {
                                double [] intArr = intArrList.toNativeArray();
                                double [] massArr = massArrList.toNativeArray();
                                for(double d: intArr)
                                {
                                    totalIntensity+=d;
                                }

                                for(int i=0; i<samIsoArr.length; i++)
                                {
                                    double mass = samIsoArr[i];
                                    double tempDiff = mass-slineMass;
                                    //   if(tempDiff>3)System.out.println("..... "+tempDiff + " " + chargeState);
                                    precursorSum+=  CalcUtil.intensitySumForSinglePeakNewMassTolerance(massArr, intArr, mass, massTol);
                                }
                                double purity = precursorSum/totalIntensity;
                       /* if(purity>1 || purity< 0 || Double.isNaN(purity)|| Double.isInfinite(purity))
                        {
                            System.out.println("''''''''Bad result");
                        }*/
                                purityList.add(purity);
                                peptideEle.setAttribute("tmt_purity", Double.toString(purity));
                            }


                        }
                    }
                    else
                    {
                        peptideEle.setAttribute("tmt_purity", Double.toString(-1.0));

                    }




                    keyIndex = Arrays.binarySearch(keys, Integer.parseInt(peptide.getScanNum()));
                    if (keyIndex < 0) //Cannot find index
                    {
                        keyIndex = -(++keyIndex); //Math.abs(++keyIndex);
                    }

                    if (conf.isMs3ScanRandom()) {  //ms3 based search
                        String ms3FileName;
                        String ms3FileName2;
                        if(ms2FileName.endsWith("_ms3"))
                        {
                            ms3FileName = this.filePath + ms2FileName + ".ms2";
                            ms3FileName2 =  ms2FileName + ".ms2";
                        }
                        else
                        {
                             ms3FileName = this.filePath + ms2FileName + ".ms3";
                             ms3FileName2 =  ms2FileName + ".ms3";
                        }


                        IndexedFile ms3Ifile = ms2Ht.get(ms3FileName);
                        SpectraDB ispectraDB3 = ms2Map.get(ms3FileName2);
                        if (null == ms3Ifile) {  //for heavy file
                            ms3Ifile = ms2Ht.get(filePath + ms2FileName.substring(1) + "." + "ms3");
                            ispectraDB3 = ms2Map.get( ms2FileName.substring(1) + "." + "ms3");
                            if (null == ms3Ifile) {  //for replaced ext

                                ms3Ifile = ms2Ht.get(ms3FileName.replace(".ms3", "_ms3.ms2"));

                                ispectraDB3 = ms2Map.get(ms3FileName2.replace(".ms3", "_ms3.ms2"));
                                if (null == ms3Ifile) {
                                    ispectraDB3 = ms2Map.get(ms3FileName2.replace(".ms3", ".ms2"));
                                    ms3Ifile = ms2Ht.get(ms3FileName.replace(".ms3", ".ms2"));
                                }

                                if (null == ms3Ifile) {

                                    //check if ms3 file can be heavy one
                                    if (ms2FileName.startsWith("H") || ms2FileName.startsWith("M")) {
                                        String heavyMs3File = this.filePath + ms2FileName.substring(1) + "_ms3.ms2";
                                        ispectraDB3 = ms2Map.get(ms3FileName2.substring(1) + "_ms3.ms2");

                                        ms3Ifile = ms2Ht.get(heavyMs3File);
                                    }
                                    if (null == ms3Ifile) {

/*
                                        //check one more time, if it is heavy search or not.
                                        boolean isHeavy = FileFilterUtil.isHeavyFile(ms2FileName + ".ms2", this.filePath);
                                        if(isHeavy) {
                                            ms3Ifile = ms2Map.get(filePath + ms2FileName.substring(1) + "." + "ms2");
                                        }

                                        if(null == iFile) {
                                        } */

                                        System.out.println("Error : cannot find the file " + ms2FileName + "." + "ms3");
                                        System.exit(0);
                                    }
                                }
                            }
                        }

                        int ms3Scan = TMTUtil.getCorrespondingMs3Scan(ms3Ifile, keys[keyIndex]);
                        //MS_READ

                        if (ms3Scan < 0) {
                            System.out.println("MISSING MS3SCAN ");
                            continue;
                        }

                        iFile = ms3Ifile;
                        ispectraDB = ispectraDB3;
                        //keyIndex

                        keyIndex = Arrays.binarySearch(iFile.getKeys(), ms3Scan);
                        if (keyIndex < 0) //Cannot find index
                        {
                            keyIndex = -(++keyIndex); //Math.abs(++keyIndex);
                        }
                    } else if (conf.getScanShift() >= 0) {  //ms2 scan shift
                        keyIndex += conf.getScanShift();
                    }// else {  //regular search


                    try {
            //            System.out.println(">><<>> ");
//
             //           System.out.println(">><<>> "+iFile.getFileName());
            //            System.out.println(">><<>> "+ispectraDB.path);
                        //MS_READ
                        String outStr = SpectraDBCalcUtil.calculateMS2Mass(iFile, range, keyIndex, null,
                                null, null, null, conf, chargeState, ispectraDB);
                        peptideEle.setAttribute("signal-noise", CensusHelper.format.format(CalcUtil.SignalToNoise));
                        peptideEle.setAttribute("report-ion-signal-noise", CensusHelper.format.format(CalcUtil.ReportIonSignalToNoise));

                        peptideEle.setAttribute("noise_level", CensusHelper.format.format(CalcUtil.NoiseLevel));

                        signalList.add(CalcUtil.SignalToNoise);
                        maxSignalToNoise = CalcUtil.SignalToNoise> maxSignalToNoise ? CalcUtil.SignalToNoise : maxSignalToNoise;
                        minSignalTONoise = CalcUtil.SignalToNoise< minSignalTONoise ? CalcUtil.SignalToNoise : minSignalTONoise;
                        //      mzValuesList = CalcUtil.MZValuesList;
                        //    tandemPeptide.setMzValues(mzValuesList);
                        if (null == outStr) {
                            chro = null;
                        } else {
                            chro.setText(outStr);
                        }
                    } catch (PrecursorNotFoundException ive) {
                        System.out.println("Precursor not found for " + pepSequence);
                        ChroProgressDialog.addMessageWithLine(progress, "Error : Precursor not found for " + pepSequence);
                        percent += eachSeg;
                        progress.setProgress((int) percent);

                        continue;
                    }

                    if (null != chro) {
                        peptideEle.addContent(chro);
                    }

                    percent += eachSeg;

                    if (null != progress) {
                        progress.setProgress((int) percent);
                    }
                    int ms3Scan;
                    // if(!ms3ScanMap.containsKey(fileName)) {


                    List<ReportIon> reportIonList = conf.getReportIonList();
                    double startmass = reportIonList.get(0).getMass();
                    double endMass = reportIonList.get(reportIonList.size()-1).getMass();



                    IndexedFile iFile2 = ms2Ht.get(filePath+ peptide.getFileName() + "." + "ms2");
                    if(null == iFile2) {

                        if(FileFilterUtil.isHeavyFile(peptide.getFileName() + ".ms2", filePath)) {
                            iFile2 = ms2Ht.get(filePath + peptide.getFileName().substring(1) + "." + "ms2");
                            if(null != iFile2)
                                peptide.setFileName(peptide.getFileNameWithScan().substring(1));
                        }
                        else if(null == iFile2) {
                            iFile2 = ms2Ht.get(filePath+ peptide.getFileName() + "." + "ms3");

                            if(null == iFile2 && FileFilterUtil.isHeavyFile(peptide.getFileName() + ".ms3", filePath)) {
                                iFile2 = ms2Ht.get(filePath + peptide.getFileName().substring(1) + "." + "ms3");
                                if(null != iFile2)
                                    peptide.setFileName(peptide.getFileNameWithScan().substring(1));

                            }
                        }



                    }
                    //check if it is heavy search or not.
                    /*if(null == iFile2) {
                        boolean isHeavy = FileFilterUtil.isHeavyFile(ms2FileName + ".ms2", this.filePath);
                        if (isHeavy) {
                            iFile2 = ms2Map.get(filePath + ms2FileName.substring(1) + "." + "ms2");
                         //   peptide.setFileName(ms2FileName.substring(1));
                        }
                    } */


                    int [] keys2 = iFile2.getKeys();

                    //  long start = System.currentTimeMillis();

                    String ms3FilePattern="orig";
                    List<edu.scripps.pms.census.model.XYPoint> xyPoints;
                    if (conf.isMs3ScanRandom()) {

                        String ms3FileName = filePath + peptide.getFileName() + ".ms3";

                        IndexedFile ms3Ifile = ms2Ht.get(ms3FileName);




                        if (null == ms3Ifile) { // for heavy file
                            ms3Ifile = ms2Ht.get(filePath+ peptide.getFileName().substring(1) + "." + "ms3");

                            if (null == ms3Ifile) { // for replaced ext

                                ms3Ifile = ms2Ht.get(ms3FileName.replace(".ms3", "_ms3.ms2"));
                                ms3FilePattern="_ms3.ms2";

                                if (null == ms3Ifile) {

                                    if (peptide.getFileName().startsWith("H")
                                            || peptide.getFileName().startsWith("M")) {
                                        String heavyMs3File =filePath + peptide.getFileName().substring(1)
                                                + "_ms3.ms2";
                                        ms3Ifile = ms2Ht.get(heavyMs3File);
                                    }
                                    if (null == ms3Ifile) {
                                        System.out.println("Error : cannot find the file " + peptide.getFileName()
                                                + "." + "ms3");
                                    }
                                }
                            }
                        }


                        int keyIndex2 = Arrays.binarySearch(keys2, Integer.parseInt(peptide.getScanNum()));

                        ms3Scan = TMTUtil.getCorrespondingMs3Scan(ms3Ifile, keys2[keyIndex2]);

                        if(ms3FilePattern.equals("orig")) {
                            //MS_READ

                            File f = new File(filePath + peptide.getFileName() + ".ms3");
                            SpectraDB ms3SpectraDB = ms3Map.get(peptide.getFileName() + ".ms3");
                            if(f.exists())
                                xyPoints = SpectrumUtil.getSpectrum(filePath, ms3SpectraDB, ms3Scan,
                                        startmass-1, endMass+1 );
                            else
                            {
                                String lms3Key = peptide.getFileName().substring(1)+".ms3";
                                SpectraDB lms3SpectraDB = ms3Map.get(lms3Key);
                                xyPoints = SpectrumUtil.getSpectrum(filePath, lms3SpectraDB, ms3Scan,
                                        startmass-1, endMass+1 );
                            }

                        }
                        else {
                            //MS_READ

                            File f = new File(filePath + peptide.getFileName() + "_ms3.ms2");
                            String ms3Key = f.getName();

                            if(f.exists())
                            {
                                SpectraDB ms3SpectraDB = ms2Map.get(ms3Key);
                                xyPoints = SpectrumUtil.getSpectrum(filePath, ms3SpectraDB, ms3Scan,
                                        startmass-1, endMass+1 );
                            }
                            else
                            {
                                SpectraDB ms3SpectraDB = ms2Map.get(ms3Key.substring(1));
                                xyPoints = SpectrumUtil.getSpectrum(filePath, ms3SpectraDB, ms3Scan,
                                        startmass-1, endMass+1 );
                            }
                        }
                    }
                    else{
                        //MS_READ
                        String ms2Key =peptide.getFileName()+".ms2";
                        SpectraDB ms2SpectraDB = ms2Map.get(ms2Key);

                        xyPoints = SpectrumUtil.getSpectrum(filePath, ms2SpectraDB,
                                Integer.parseInt(peptide.getScanNum()),startmass -1 , endMass + 1);
                    }
                    String pepKey = peptide.getSequence() + peptide.getChargeState() + peptide.getScanNum();
                    FileUtil.writeJSON(gson.toJson(xyPoints), xyvaluePath +  pepKey + "_xy.JSON");


                    System.out.print(pepCount);
                    System.out.print("/");
                    System.out.print(redundantPeptideNum);
                    System.out.print(" peptides, ");
                    System.out.print((int) percent);
                    //System.out.print(" % is complete\n");
                    System.out.print(" % is complete\r");



                    if (conf.isPrintLog()) {
                        p.println((int) percent + "\t" + pepCount + "\t" + redundantPeptideNum);
                    }



                    //   tandemTagPeptideList.add(tandemPeptide);
                    if (null != chro) {
                        proteinEle.addContent(peptideEle);
                    }
                }

                if (proteinEle.getChildren().size() > 0) {

                    Element redunEle = new Element("redundant");

                    int rPepCount = 0;

                    for (Iterator<Protein> rItr = aList.iterator(); rItr.hasNext();) {
                        rPepCount++;
                        Protein rp = rItr.next();

                        Element rpEle = new Element("protein");
                        rpEle.setAttribute("locus", rp.getLocus());
                        rpEle.setAttribute("seq_ct", rp.getSeqCount());
                        rpEle.setAttribute("spec_ct", rp.getSpectrumCount());
                        rpEle.setAttribute("seq_cov", rp.getSeqCoverage());
                        rpEle.setAttribute("length", rp.getLength());
                        rpEle.setAttribute("molwt", rp.getMolWt());
                        rpEle.setAttribute("pi", rp.getPI());
                        rpEle.setAttribute("val", rp.getValidation());

                        try {
                            rpEle.setAttribute("desc", rp.getDescription());
                        } catch (org.jdom.IllegalDataException ide) {
                            rpEle.setAttribute("desc", StringUtil.removeIsoControlChar(rp.getDescription()));
                        }

                        redunEle.addContent(rpEle);
                    }

                    if (rPepCount > 0) {
                        proteinEle.addContent(redunEle);
                    }

                    rootEle.addContent(proteinEle);
                    aList.clear();
                } else {
                    aList.add(protein);
                    //	    set.add(protein.getLocus());
                }

                //FileUtil.writeJSON(gson.toJson(tandemTagPeptideList), filePath +  File.separator + MASTER_JSON + File.separator +protein.getLocus() + ".JSON");
            }
            List<XYPoint> purityFrequencyPointList = new ArrayList<>();

            Histogram purityHistogram = new Histogram(11,0,1.1);
            purityHistogram.loadData(purityList);
            double [] purityBins = purityHistogram.getBins();
            int[] purityFrequency = purityHistogram.getFreqArr();
            for(int i=0; i<purityBins.length; i++)
            {
                purityFrequencyPointList.add(new XYPoint(purityBins[i],purityFrequency[i],Double.toString(purityBins[i])));
            }



            List<XYPoint> signalHistPointList = new ArrayList<>();
            //Histogram signalHist= new Histogram(20,minSignalTONoise,maxSignalToNoise);
            Histogram signalHist= new Histogram(20,0,200);  //max value is 200
            signalHist.loadData(signalList);
            double [] signalBins = signalHist.getBins();
            int[] signalFreq = signalHist.getFreqArr();
            for(int i=0; i<signalBins.length; i++)
            {
                signalHistPointList.add(new XYPoint(signalBins[i],signalFreq[i],Double.toString(signalBins[i])));
            }
            FileUtil.writeJSON(gson.toJson(purityFrequencyPointList), filePath+PURITY_GRAPH_PATH);
            FileUtil.writeJSON(gson.toJson(signalHistPointList), filePath+SIGNAL_GRAPH_PATH);




            Document doc = new Document(rootEle);
            OutputStream os = new FileOutputStream(filePath + "census_chro.xml");
            XMLOutputter outputter = new XMLOutputter();
            outputter.setFormat(Format.getPrettyFormat());
            outputter.output(doc, os);
            os.close();

            System.out.println("\n100% complete");
            //System.out.println( System.out.println"\n100% complete");

        } catch (IOException e) {
            System.out.println("IO Error while generating msms chro file : " + pepSequence + e);
            e.printStackTrace();
            throw new IOException(e.toString());
        } catch (java.lang.IndexOutOfBoundsException e) {
            System.out.println("Error while generating msms chro file : " + pepSequence + e);
            e.printStackTrace();
            throw new Exception(e.toString());
        } catch (Exception e) {
            System.out.println("Error while generating msms chro file : " + pepSequence + e);
            e.printStackTrace();
            throw new Exception(e.toString());
        } finally {
            if (null != p) {
                p.close();
            }

            if (null != out) {
                out.close();
            }

            //Close all random files
            for (String  e : ms2Ht.keySet()) {
                iFile = ms2Ht.get(e);

                if (null != iFile) {
                    iFile.close();
                }
            }

            for(Map.Entry<String,SpectraDB> entry: ms2Map.entrySet())
            {
                if(entry.getValue() !=null)
                {
                    SpectraDB db = entry.getValue();
                    db.close();
                }
            }
            for(Map.Entry<String,SpectraDB> entry: ms1Map.entrySet())
            {
                if(entry.getValue() !=null)
                {
                    SpectraDB db = entry.getValue();
                    db.close();
                }
            }
            for(Map.Entry<String,SpectraDB> entry: ms3Map.entrySet())
            {
                if(entry.getValue() !=null)
                {
                    SpectraDB db = entry.getValue();
                    db.close();
                }
            }
        }

    }




    public void createDIALabelfree() throws IOException, Exception {

      int[] keys;

      this.filePath = conf.getFilePath();

      //if( !filePath.endsWith(File.separator) )
      //    filePath += File.separator;
        /*
        //file could be generated from either linux or window
        if( !this.filePath.endsWith("/") && !this.filePath.endsWith("\\") ) {
            //for linux
            if(this.filePath.startsWith("/"))
                this.filePath += "/";
            else //for window
                this.filePath += "\\";
        }
         */
      Hashtable<String, IndexedFile> ht = createIndexedFiles(filePath, CensusConstants.MS2_FILE);

      IndexedFile iFile;
      BufferedOutputStream out = null;
      PrintStream p = null;

      long startTime = System.currentTimeMillis();
      String pepSequence = null;

      try {

        /**
         * ****************************************************************
         * Read DTASelect.txt file to find spectrum range for each peptide
         *****************************************************************
         */
        ChroProgressDialog.addMessageWithLine(progress, "");

        IsotopeReader isoReader = null;

        if (conf.isXmlConf()) {
          isoReader = new IsotopeReader(conf.getRootConfEle());
        } else {
          isoReader = new IsotopeReader(isotopeFile);
        }

        IdentificationReader idReader = BaseIdentificationReader.getIdentificationInst(isoReader);
        SpecRangeGenerator rangeGen = new SpecRangeGenerator();
        int redundantPeptideNum = idReader.getTotalPeptideNumber();

        TIntLongHashMap index;

        IsotopeTable<String, int[]> isoTable = isoReader.getIsotope();
        int[] sampleNterm = isoTable.get("sampleNTERM");
        int[] sampleCterm = isoTable.get("sampleCTERM");

        //int redundantPeptideNum = -1; //dtaReader.getRedundantPeptideNum();
        //int redundantPeptideNum = DTASelectFilterReader.getTotalPeptideNumber(conf.getIdFileName());
        //increase status bar
        double percent = 0.0;
        double eachSeg = 100.0 / redundantPeptideNum;
        int pepCount = 0;

        //IsotopeDist sampleDist;
        //IsotopeDist refDist;
        IsotopeDist sampleDist;
        Protein protein;
        Peptide peptide;

        Element rootEle = this.createXmlChroHeader(2, conf.getExpType());

        ElementComposition element;
        ElementComposition totalElement;

        //int[] elementSampleArr;
        int keyIndex;
        int start;
        int last;

        double samplePrecursor;

        double[][] bionSample;
        double[][] yionSample;

        Element proteinEle = null;
        Element peptideEle = null;

        for (Iterator<Protein> itr = idReader.getProteins(); itr.hasNext();) {
          protein = itr.next();

          proteinEle = new Element("protein");
          proteinEle.setAttribute("locus", protein.getLocus());
          proteinEle.setAttribute("seq_ct", protein.getSeqCount());
          proteinEle.setAttribute("spec_ct", protein.getSpectrumCount());
          proteinEle.setAttribute("seq_cov", protein.getSeqCoverage());
          proteinEle.setAttribute("length", protein.getLength());
          proteinEle.setAttribute("molwt", protein.getMolWt());
          proteinEle.setAttribute("pi", protein.getPI());
          proteinEle.setAttribute("val", protein.getValidation());

          try {
            proteinEle.setAttribute("desc", protein.getDescription());
          } catch (org.jdom.IllegalDataException ide) {
            proteinEle.setAttribute("desc", StringUtil.removeIsoControlChar(protein.getDescription()));
          }

          for (Iterator<Peptide> pepItr = protein.getPeptides(); pepItr.hasNext();) {
            peptide = pepItr.next();
            pepCount++;

            pepSequence = peptide.getSequence();

            //What is the purpose of hs???
            char[] ch = pepSequence.substring(2, peptide.getSequence().length() - 2).toCharArray();

            try {
              totalElement = new ElementComposition(ch, 0, ch.length, isoTable);
              totalElement.calculate();
            } catch (InvalidAAException ive) {
              System.out.println("Not Quantifiable peptide : " + pepSequence);

              percent += eachSeg;
              if (null != progress) {
                ChroProgressDialog.addMessageWithLine(progress, "Not Quantifiable peptide : " + pepSequence);
                progress.setProgress((int) percent);

              }
              continue;
            }

            int chargeState = Integer.parseInt(peptide.getChargeState());

            int pepLength = 0;

            for (int i = 0; i < ch.length; i++) {
              if (ch[i] == '*' || ch[i] == '@' || ch[i] == '#') {
                continue;
              }

              pepLength++;
            }

            bionSample = new double[pepLength-1][chargeState * 3];
            //Yions
            yionSample = new double[pepLength-1][chargeState * 3];

            int pepIndex = 0;

            //System.out.println("aamass" + massTolerance + " " + conf.getMassTolerance());
            for (int i = 0; i < ch.length-1; i++) {
              if (ch[i] == '*' || ch[i] == '@' || ch[i] == '#') {
                continue;
              }

              try {
                element = new ElementComposition(ch, 0, i + 1, isoTable);
                element.calculate();

              } catch (InvalidAAException ive) {
                System.out.println("Not Quantifiable peptide : " + pepSequence);

                percent += eachSeg;
                if (null != progress) {
                  ChroProgressDialog.addMessageWithLine(progress, "Not Quantifiable peptide : " + pepSequence);
                  progress.setProgress((int) percent);

                }
                continue;
              }

              //Y ions
              sampleDist = new IsotopeDist(
                getComplementaryComposition(totalElement.getElementSampleArr(), element.getElementSampleArr(), sampleNterm, sampleCterm), element.getModShift(), true);

              //System.out.println("s>||>==" + yionSample.length + " " + pepIndex + "\t" + chargeState + "\t" + sampleDist.getStartMass());
/*
              for(int csIndex=0;csIndex<chargeState;csIndex++) {

                int tempCs = csIndex + 1;
              //  System.out.println("s==" + csIndex + "\t" + (csIndex*3+2) + "\t" + sampleDist.getStartMass() + " " +  ((pepIndex + 1) % pepLength));
                yionSample[(pepIndex) % pepLength][csIndex*3+2] = (sampleDist.getEndMass() + 3 * PROTON_MASS) / tempCs + massTolerance;
                yionSample[(pepIndex) % pepLength][csIndex*3+1] = (sampleDist.getStartMass() + 3 * PROTON_MASS) / tempCs - massTolerance;
                yionSample[(pepIndex) % pepLength][csIndex*3+0] = (sampleDist.getAvgMass() + 3 * PROTON_MASS) / tempCs;
              }
*/

//                yionSample[(pepIndex) % pepLength][(chargeState-1)*3+2] = (sampleDist.getEndMass() + 3 * PROTON_MASS) / chargeState + massTolerance;
                //yionSample[(pepIndex) % pepLength][(chargeState-1)*3+1] = (sampleDist.getStartMass() + 3 * PROTON_MASS) / chargeState - massTolerance;
//                yionSample[(pepIndex) % pepLength][(chargeState-1)*3+0] = (sampleDist.getAvgMass() + 3 * PROTON_MASS) / chargeState;

                yionSample[(pepIndex) % pepLength][(chargeState-1)*3+1] = (sampleDist.getStartMass() + 1 * PROTON_MASS);
//System.out.print(yionSample[(pepIndex) % pepLength][(chargeState-1)*3+2] + " ");
//System.out.print(yionSample[(pepIndex) % pepLength][(chargeState-1)*3+1] + " ");
//System.out.println(yionSample[(pepIndex) % pepLength][(chargeState-1)*3+0]);
//System.out.println(pepIndex + " " + sampleDist.getStartMass());

              /*
              switch (chargeState) {


                case 3:
                  yionSample[(pepIndex + 1) % pepLength][8] = (sampleDist.getEndMass() + 3 * PROTON_MASS) / 3 + massTolerance;
                  yionSample[(pepIndex + 1) % pepLength][7] = (sampleDist.getStartMass() + 3 * PROTON_MASS) / 3 - massTolerance;
                  yionSample[(pepIndex + 1) % pepLength][6] = (sampleDist.getAvgMass() + 3 * PROTON_MASS) / 3;

                case 2:
                  yionSample[(pepIndex + 1) % pepLength][5] = (sampleDist.getEndMass() + 2 * PROTON_MASS) / 2 + massTolerance;
                  yionSample[(pepIndex + 1) % pepLength][4] = (sampleDist.getStartMass() + 2 * PROTON_MASS) / 2 - massTolerance;
                  yionSample[(pepIndex + 1) % pepLength][3] = (sampleDist.getAvgMass() + 2 * PROTON_MASS) / 2;

                case 1:
                  yionSample[(pepIndex + 1) % pepLength][2] = sampleDist.getEndMass() + 1 * PROTON_MASS + massTolerance; //add proton to give b fragment ion
                  yionSample[(pepIndex + 1) % pepLength][1] = sampleDist.getStartMass() + 1 * PROTON_MASS - massTolerance;
                  yionSample[(pepIndex + 1) % pepLength][0] = sampleDist.getAvgMass() + 1 * PROTON_MASS;

                default:
                  break;
              }  */

              element.calculateBion();
              //element.printComposition();

              sampleDist = new IsotopeDist(element.getElementSampleArr(), element.getModShift(), true);


           //   System.out.println("s>>==" + pepIndex + "\t" + chargeState + "\t" + sampleDist.getStartMass());
/*
              for(int csIndex=0;csIndex<chargeState;csIndex++) {

                int tempCs = csIndex + 1;
                bionSample[pepIndex][csIndex*3+2] = (sampleDist.getEndMass() + 2 * PROTON_MASS) / tempCs + massTolerance;
                bionSample[pepIndex][csIndex*3+1] = (sampleDist.getStartMass() + 2 * PROTON_MASS) / tempCs - massTolerance;
                bionSample[pepIndex][csIndex*3] = (sampleDist.getAvgMass() + 2 * PROTON_MASS) / tempCs;
              }
*/
                bionSample[pepIndex][(chargeState-1)*3+1] = (sampleDist.getStartMass()-0.0005);
              //  bionSample[pepIndex][(chargeState-1)*3+2] = (sampleDist.getEndMass() + 2 * PROTON_MASS) / chargeState + massTolerance;
               // bionSample[pepIndex][(chargeState-1)*3+1] = (sampleDist.getStartMass() + 2 * PROTON_MASS) / chargeState - massTolerance;
             //   bionSample[pepIndex][(chargeState-1)*3] = (sampleDist.getAvgMass() + 2 * PROTON_MASS) / chargeState;


              pepIndex++;
            }



            try {
              element = new ElementComposition(ch, 0, ch.length, isoTable);
              element.calculate();
            } catch (InvalidAAException ive) {
              System.out.println("Not Quantifiable peptide : " + pepSequence);

              percent += eachSeg;
              if (null != progress) {
                ChroProgressDialog.addMessageWithLine(progress, "Not Quantifiable peptide : " + pepSequence);
                progress.setProgress((int) percent);

              }
              continue;
            }

            sampleDist = new IsotopeDist(element.getElementSampleArr(), element.getModShift(), true);


/*
            for(int i=0;i<bionSample.length;i++)
            {
              for(int j=0;j<bionSample[i].length;j++)
                System.out.print("b + " +i + " " + j + " " + bionSample[i][j] + " ");

              System.out.println(" ");
            }
            System.out.println(" ");

            for(int i=0;i<yionSample.length;i++)
            {
              for(int j=0;j<yionSample[i].length;j++)
                System.out.print("y + " +i + " " + j + " " + yionSample[i][j] + " ");

              System.out.println(" ");
            }
            System.out.println(" ");
*/


            double calcSamMass = (sampleDist.getAvgMass() + chargeState * PROTON_MASS) / chargeState;

            conf.setCalcSamAvgMass(calcSamMass);

            //peptideEle.setAttribute("lightStartMass", String.valueOf(sampleDist.getStartMass()));
            //peptideEle.setAttribute("lightAvgMass", String.valueOf(conf.getCalcSamAvgMass()));

            peptideEle = new Element("peptide");

            peptideEle.setAttribute("lightStartMass", String.valueOf(0.0));
            peptideEle.setAttribute("lightAvgMass", String.valueOf(0.0));


            peptideEle = this.createXmlChroPeptideTitle(false, peptide);

            //String ms2FileName = peptide.getFileName().substring(0, peptide.getFileName().indexOf("."));
            String ms2FileName = peptide.getFileName();

            StringBuffer rangeKey = new StringBuffer();
            rangeKey.append(protein.getLocus());
            rangeKey.append(ms2FileName);
            rangeKey.append(peptide.getSequence().substring(2, peptide.getSequence().length() - 2));
            SpecRange range = rangeGen.getSpecRange(rangeKey.toString());

            if (null == range) {
              int tmpScanNum = Integer.parseInt(peptide.getScanNum());
              range = new SpecRange(tmpScanNum, tmpScanNum);
              peptideEle.setAttribute("start_scan", peptide.getScanNum());
              peptideEle.setAttribute("end_scan", peptide.getScanNum());
            } else {
              peptideEle.setAttribute("start_scan", String.valueOf(range.getMin()));
              peptideEle.setAttribute("end_scan", String.valueOf(range.getMax()));
            }

            //sb.append("\tStartScan\tEndScan\tDTAPeakStart\tDTAPeakEnd\n");
            Element chro = new Element("chro");
            //output.append("[CHROMATOGRAMS]\tSCAN\tSAMPLE\tREFERENCE\n");

            iFile = ht.get(this.filePath + ms2FileName + "." + "ms2");
            keys = iFile.getKeys();

            //samplePrecursor = conf.getPrecursor( (sampleDist.getAvgMass()+chargeState)/chargeState );
            //refPrecursor = conf.getPrecursor( (refDist.getAvgMass()+chargeState)/chargeState );
            keyIndex = Arrays.binarySearch(keys, Integer.parseInt(peptide.getScanNum()));

            if (keyIndex < 0) //Cannot find index
            {
              keyIndex = -(++keyIndex); //Math.abs(++keyIndex);
            }
            /**
             * Find start scan number same as a precursor of
             * sampleAvgMass Then, later the program will find following
             * scan # quicker.
             *
             */
                    /* End of finding start spectrum number */

            samplePrecursor = 1;

            TIntDoubleHashMap precursorMap = iFile.getPrecursorMap();

            // {
            //   System.out.println("Outside isowindow for sample precursor " + keys[keyIndex] + " " + (sampleDist.getAvgMass()+chargeState)/chargeState + " " + precursorMap.get(keys[keyIndex]) + " " + samplePrecursor);
            //}
            //if(precursorMap.get(keys[keyIndex+diff]) != refPrecursor)
            //   System.out.println("Outside isowindow for ref precursor " + keys[keyIndex+diff] + " " + (refDist.getAvgMass()+chargeState)/chargeState  + " " + precursorMap.get(keys[keyIndex+diff]) + " " + refPrecursor);
            //System.out.println("set text" + iFile + " " + range  + " " + keyIndex + " " + diff + " " + bionSample + " " + bionRef + " " + yionSample + " " + yionRef + " " + samplePrecursor + " " + refPrecursor + " " + conf + " " + chargeState);
            //chro.setText( CalcUtil.calculateMS2Mass( iFile, range, keyIndex, diff, bionSample, bionRef, yionSample, yionRef, samplePrecursor, refPrecursor, conf, chargeState) );
            try {
              //String outStr = CalcUtil.calculateMS2Mass(iFile, range, keyIndex, bionSample, bionRef, yionSample, yionRef, conf, chargeState);
              //String outStr = CalcUtil.calculateMS2Mass(iFile, range, keyIndex, bionSample, null, yionSample, null, conf, chargeState);



/*
for(int i=0;i<yionSample.length;i++) {
for(int j=0;j<yionSample[i].length;j++) {

System.out.print(yionSample[i][j] + " ");
}
System.out.println("==>>" + i);
}
              System.out.println("--------------------");
              for(int i=0;i<bionSample.length;i++) {
                for(int j=0;j<bionSample[i].length;j++) {

                  System.out.print(bionSample[i][j] + " ");
                }
                System.out.println("==>>" + i);
              }
System.out.println(chargeState);
System.exit(0);
              */
              String outStr = CalcUtil.calculateMS2Labelfree(iFile, range, keyIndex, bionSample, yionSample, conf, chargeState);

              if (outStr == null) {
                continue;
              }

              chro.setText(outStr);
            } catch (PrecursorNotFoundException ive) {
              System.out.println("Precursor not found for " + pepSequence);

              if (null != progress) {
                ChroProgressDialog.addMessageWithLine(progress, "Error : Precursor not found for " + pepSequence);
                percent += eachSeg;
                progress.setProgress((int) percent);
              }

              ive.printStackTrace();

              continue;
            } catch (Exception e) {
              System.out.println("Precursor not found for " + pepSequence);

              if (null != progress) {
                ChroProgressDialog.addMessageWithLine(progress, "Error : Precursor not found for " + pepSequence);
                percent += eachSeg;
                progress.setProgress((int) percent);
              }
              e.printStackTrace();

              continue;
            }

            peptideEle.addContent(chro);

            Element fragEle = new Element("frag");

            Element bSample = new Element("bs");
            StringBuffer tempSb = new StringBuffer();
            for (int i = 0; i < bionSample.length; i++) {
              int j;
              for (j = 0; j < bionSample[i].length - 1; j++) {
                tempSb.append(formatter.format(bionSample[i][j])).append(" ");
              }

              tempSb.append(formatter.format(bionSample[i][j])).append(",");
            }

            bSample.setText(tempSb.toString());
            fragEle.addContent(bSample);

            tempSb.delete(0, tempSb.length());


            Element ySample = new Element("ys");
            tempSb.delete(0, tempSb.length());
            for (int i = 0; i < yionSample.length; i++) {
              int j;
              for (j = 0; j < yionSample[i].length - 1; j++) {
                tempSb.append(formatter.format(yionSample[i][j])).append(" ");
              }

              tempSb.append(formatter.format(yionSample[i][j])).append(",");
            }

            ySample.setText(tempSb.toString());
            fragEle.addContent(ySample);

            tempSb.delete(0, tempSb.length());

            peptideEle.addContent(fragEle);

            percent += eachSeg;

            if (null != progress) {
              progress.setProgress((int) percent);
            }

            System.out.print(pepCount);
            System.out.print("/");
            System.out.print(redundantPeptideNum);
            System.out.print(" peptides, ");
            System.out.print((int) percent);
            //System.out.print(" % is complete\n");
            System.out.print(" % is complete\r");

            proteinEle.addContent(peptideEle);
          }

          if (proteinEle.getChildren().size() > 0) {
            rootEle.addContent(proteinEle);
          }
        }

        Document doc = new Document(rootEle);
        OutputStream os = new FileOutputStream(filePath + "census_chro.xml");
        XMLOutputter outputter = new XMLOutputter();
        outputter.setFormat(Format.getPrettyFormat());
        outputter.output(doc, os);
        os.close();

        System.out.println("\n100% complete");
        //System.out.println( System.out.println"\n100% complete");

        //Create relex.chro file
            /*
            out = new BufferedOutputStream(new FileOutputStream(filePath + "relex.chro"));
            p = new PrintStream(out);
            p.print(output.toString());
             */
      } catch (IOException e) {
        System.out.println("IO Error while generating msms chro file : " + pepSequence + e);
        e.printStackTrace();
        throw new IOException(e.toString());
      } catch (Exception e) {
        System.out.println("Error while generating msms chro file : " + pepSequence + e);
        e.printStackTrace();
        throw new Exception(e.toString());
      } finally {
        if (null != p) {
          p.close();
        }

        if (null != out) {
          out.close();
        }

        //Close all random files
        for (Enumeration e = ht.keys(); e.hasMoreElements();) {
          iFile = ht.get(e.nextElement());

          if (null != iFile) {
            iFile.close();
          }
        }

      }

    }


    //include data independent
    public void createMsmsXmlChro() throws IOException, Exception {
        createMsmsXmlChro(null);
    }

    public void createMsmsXmlChro(ChroProgressDialog progress) throws IOException, Exception {

        this.progress = progress;
        int[] keys;

        this.filePath = conf.getFilePath();

        //if( !filePath.endsWith(File.separator) )
        //    filePath += File.separator;
        /*
        //file could be generated from either linux or window
        if( !this.filePath.endsWith("/") && !this.filePath.endsWith("\\") ) {
            //for linux
            if(this.filePath.startsWith("/"))
                this.filePath += "/";
            else //for window
                this.filePath += "\\";
        }
         */
        Hashtable<String, IndexedFile> ht = createIndexedFiles(filePath, CensusConstants.MS2_FILE);

        IndexedFile iFile;
        BufferedOutputStream out = null;
        PrintStream p = null;

        long startTime = System.currentTimeMillis();
        String pepSequence = null;

        try {

            /**
             * ****************************************************************
             * Read DTASelect.txt file to find spectrum range for each peptide
             *****************************************************************
             */
            ChroProgressDialog.addMessageWithLine(progress, "");

            IsotopeReader isoReader = null;

            if (conf.isXmlConf()) {
                isoReader = new IsotopeReader(conf.getRootConfEle());
            } else {
                isoReader = new IsotopeReader(isotopeFile);
            }

            IdentificationReader idReader = BaseIdentificationReader.getIdentificationInst(isoReader);
            SpecRangeGenerator rangeGen = new SpecRangeGenerator();
            int redundantPeptideNum = idReader.getTotalPeptideNumber();

            /*
            System.out.print("Parsing DTASelect.txt...");
            ChroProgressDialog.addMessageWithLine(progress, "Parsing DTASelect.txt...");
            IdentificationReader idReader = BaseIdentificationReader.getIdentificationInst(filePath);
            int redundantPeptideNum = idReader.getTotalPeptideNumber();
            SpecRangeGenerator rangeGen = SpecRangeGenerator.getSpecRangeGenerator(idReader);
            Hashtable tempht = rangeGen.getTable();
            System.out.println("done.");
            ChroProgressDialog.addMessageWithLine(progress, "done");
             */
            TIntLongHashMap index;

            //IsotopeReader isoReader = new IsotopeReader(isotopeFile);
            IsotopeTable<String, int[]> isoTable = isoReader.getIsotope();
            int[] sampleNterm = isoTable.get("sampleNTERM");
            int[] sampleCterm = isoTable.get("sampleCTERM");
            int[] refNterm = isoTable.get("refNTERM");
            int[] refCterm = isoTable.get("refCTERM");

            //int redundantPeptideNum = -1; //dtaReader.getRedundantPeptideNum();
            //int redundantPeptideNum = DTASelectFilterReader.getTotalPeptideNumber(conf.getIdFileName());
            //increase status bar
            double percent = 0.0;
            double eachSeg = 100.0 / redundantPeptideNum;
            int pepCount = 0;

            //IsotopeDist sampleDist;
            //IsotopeDist refDist;
            IsotopeDist sampleDist;
            IsotopeDist refDist;
            Protein protein;
            Peptide peptide;

            Element rootEle = this.createXmlChroHeader(2, conf.getExpType());

            ElementComposition element;
            ElementComposition totalElement;

            //int[] elementSampleArr;
            int keyIndex;
            int start;
            int last;

            double samplePrecursor;
            double refPrecursor;

            double[][] bionSample;
            double[][] bionRef;
            double[][] yionSample;
            double[][] yionRef;

            Element proteinEle = null;
            Element peptideEle = null;

            for (Iterator<Protein> itr = idReader.getProteins(); itr.hasNext();) {
                protein = itr.next();

                proteinEle = new Element("protein");
                proteinEle.setAttribute("locus", protein.getLocus());
                proteinEle.setAttribute("seq_ct", protein.getSeqCount());
                proteinEle.setAttribute("spec_ct", protein.getSpectrumCount());
                proteinEle.setAttribute("seq_cov", protein.getSeqCoverage());
                proteinEle.setAttribute("length", protein.getLength());
                proteinEle.setAttribute("molwt", protein.getMolWt());
                proteinEle.setAttribute("pi", protein.getPI());
                proteinEle.setAttribute("val", protein.getValidation());

                try {
                    proteinEle.setAttribute("desc", protein.getDescription());
                } catch (org.jdom.IllegalDataException ide) {
                    proteinEle.setAttribute("desc", StringUtil.removeIsoControlChar(protein.getDescription()));
                }

                for (Iterator<Peptide> pepItr = protein.getPeptides(); pepItr.hasNext();) {
                    peptide = pepItr.next();
                    pepCount++;

                    pepSequence = peptide.getSequence();

                    //What is the purpose of hs???
                    char[] ch = pepSequence.substring(2, peptide.getSequence().length() - 2).toCharArray();

                    try {
                        totalElement = new ElementComposition(ch, 0, ch.length, isoTable);
                        totalElement.calculate();
                    } catch (InvalidAAException ive) {
                        System.out.println("Not Quantifiable peptide : " + pepSequence);

                        percent += eachSeg;
                        if (null != progress) {
                            ChroProgressDialog.addMessageWithLine(progress, "Not Quantifiable peptide : " + pepSequence);
                            progress.setProgress((int) percent);

                        }
                        continue;
                    }

                    int chargeState = Integer.parseInt(peptide.getChargeState());

                    int pepLength = 0;

                    for (int i = 0; i < ch.length; i++) {
                        if (ch[i] == '*' || ch[i] == '@' || ch[i] == '#') {
                            continue;
                        }

                        pepLength++;
                    }

                    bionSample = new double[pepLength][chargeState * 3];
                    bionRef = new double[pepLength][chargeState * 3];
                    //Yions
                    yionSample = new double[pepLength][chargeState * 3];
                    yionRef = new double[pepLength][chargeState * 3];

                    int pepIndex = 0;

                    //System.out.println("aamass" + massTolerance + " " + conf.getMassTolerance());
                    for (int i = 0; i < ch.length; i++) {
                        if (ch[i] == '*' || ch[i] == '@' || ch[i] == '#') {
                            continue;
                        }

                        try {
                            element = new ElementComposition(ch, 0, i + 1, isoTable);
                            element.calculate();

                        } catch (InvalidAAException ive) {
                            System.out.println("Not Quantifiable peptide : " + pepSequence);

                            percent += eachSeg;
                            if (null != progress) {
                                ChroProgressDialog.addMessageWithLine(progress, "Not Quantifiable peptide : " + pepSequence);
                                progress.setProgress((int) percent);

                            }
                            continue;
                        }

                        //Y ions
                        sampleDist = new IsotopeDist(
                                getComplementaryComposition(totalElement.getElementSampleArr(), element.getElementSampleArr(), sampleNterm, sampleCterm), element.getModShift(), true);

                        refDist = new IsotopeDist(
                                getComplementaryComposition(totalElement.getElementRefArr(), element.getElementRefArr(), refNterm, refCterm), element.getModShift(), false);   //fix this

                        switch (chargeState) {
                            case 3:
                                yionSample[(pepIndex + 1) % pepLength][8] = (sampleDist.getEndMass() + 3 * PROTON_MASS) / 3 + massTolerance;
                                yionSample[(pepIndex + 1) % pepLength][7] = (sampleDist.getStartMass() + 3 * PROTON_MASS) / 3 - massTolerance;
                                yionSample[(pepIndex + 1) % pepLength][6] = (sampleDist.getAvgMass() + 3 * PROTON_MASS) / 3;
                                yionRef[(pepIndex + 1) % pepLength][8] = (refDist.getEndMass() + 3 * PROTON_MASS) / 3 + massTolerance;
                                yionRef[(pepIndex + 1) % pepLength][7] = (refDist.getStartMass() + 3 * PROTON_MASS) / 3 - massTolerance;
                                yionRef[(pepIndex + 1) % pepLength][6] = (refDist.getAvgMass() + 3 * PROTON_MASS) / 3;

                            case 2:
                                yionSample[(pepIndex + 1) % pepLength][5] = (sampleDist.getEndMass() + 2 * PROTON_MASS) / 2 + massTolerance;
                                yionSample[(pepIndex + 1) % pepLength][4] = (sampleDist.getStartMass() + 2 * PROTON_MASS) / 2 - massTolerance;
                                yionSample[(pepIndex + 1) % pepLength][3] = (sampleDist.getAvgMass() + 2 * PROTON_MASS) / 2;
                                yionRef[(pepIndex + 1) % pepLength][5] = (refDist.getEndMass() + 2 * PROTON_MASS) / 2 + massTolerance;
                                yionRef[(pepIndex + 1) % pepLength][4] = (refDist.getStartMass() + 2 * PROTON_MASS) / 2 - massTolerance;
                                yionRef[(pepIndex + 1) % pepLength][3] = (refDist.getAvgMass() + 2 * PROTON_MASS) / 2;

                            case 1:
                                yionSample[(pepIndex + 1) % pepLength][2] = sampleDist.getEndMass() + 1 * PROTON_MASS + massTolerance; //add proton to give b fragment ion
                                yionSample[(pepIndex + 1) % pepLength][1] = sampleDist.getStartMass() + 1 * PROTON_MASS - massTolerance;
                                yionSample[(pepIndex + 1) % pepLength][0] = sampleDist.getAvgMass() + 1 * PROTON_MASS;
                                yionRef[(pepIndex + 1) % pepLength][2] = refDist.getEndMass() + 1 * PROTON_MASS + massTolerance;
                                yionRef[(pepIndex + 1) % pepLength][1] = refDist.getStartMass() + 1 * PROTON_MASS - massTolerance;
                                yionRef[(pepIndex + 1) % pepLength][0] = refDist.getAvgMass() + 1 * PROTON_MASS;

                            default:
                                break;
                        }

                        element.calculateBion();
                        //element.printComposition();

                        sampleDist = new IsotopeDist(element.getElementSampleArr(), element.getModShift(), true);
                        refDist = new IsotopeDist(element.getElementRefArr(), element.getModShift(), false);

                        switch (chargeState) {
                            case 3:
                                bionSample[pepIndex][8] = (sampleDist.getEndMass() + 2 * PROTON_MASS) / 3 + massTolerance;
                                bionSample[pepIndex][7] = (sampleDist.getStartMass() + 2 * PROTON_MASS) / 3 - massTolerance;
                                bionSample[pepIndex][6] = (sampleDist.getAvgMass() + 2 * PROTON_MASS) / 3;
                                bionRef[pepIndex][8] = (refDist.getEndMass() + 2 * PROTON_MASS) / 3 + massTolerance;
                                bionRef[pepIndex][7] = (refDist.getStartMass() + 2 * PROTON_MASS) / 3 - massTolerance;
                                bionRef[pepIndex][6] = (refDist.getAvgMass() + 2 * PROTON_MASS) / 3;

                            case 2:
                                bionSample[pepIndex][5] = (sampleDist.getEndMass() + 1 * PROTON_MASS) / 2 + massTolerance;
                                bionSample[pepIndex][4] = (sampleDist.getStartMass() + 1 * PROTON_MASS) / 2 - massTolerance;
                                bionSample[pepIndex][3] = (sampleDist.getAvgMass() + 1 * PROTON_MASS) / 2;
                                bionRef[pepIndex][5] = (refDist.getEndMass() + 1 * PROTON_MASS) / 2 + massTolerance;
                                bionRef[pepIndex][4] = (refDist.getStartMass() + 1 * PROTON_MASS) / 2 - massTolerance;
                                bionRef[pepIndex][3] = (refDist.getAvgMass() + 1 * PROTON_MASS) / 2;

                            case 1:
                                bionSample[pepIndex][2] = sampleDist.getEndMass() + massTolerance;
                                bionSample[pepIndex][1] = sampleDist.getStartMass() - massTolerance;
                                bionSample[pepIndex][0] = sampleDist.getAvgMass();
                                bionRef[pepIndex][2] = refDist.getEndMass() + massTolerance;
                                bionRef[pepIndex][1] = refDist.getStartMass() - massTolerance;
                                bionRef[pepIndex][0] = refDist.getAvgMass();

                            default:
                                break;

                        }

                        pepIndex++;
                    }

                    try {
                        element = new ElementComposition(ch, 0, ch.length, isoTable);
                        element.calculate();
                    } catch (InvalidAAException ive) {
                        System.out.println("Not Quantifiable peptide : " + pepSequence);

                        percent += eachSeg;
                        if (null != progress) {
                            ChroProgressDialog.addMessageWithLine(progress, "Not Quantifiable peptide : " + pepSequence);
                            progress.setProgress((int) percent);

                        }
                        continue;
                    }

                    sampleDist = new IsotopeDist(element.getElementSampleArr(), element.getModShift(), true);

                    int[] tempA = element.getElementSampleArr();

                    refDist = new IsotopeDist(element.getElementRefArr(), element.getModShift(), false);

                    switch (chargeState) {
                        case 3:
                            yionSample[0][8] = (sampleDist.getEndMass() + 3 * PROTON_MASS) / 3 + massTolerance;
                            yionSample[0][7] = (sampleDist.getStartMass() + 3 * PROTON_MASS) / 3 - massTolerance;
                            yionSample[0][6] = (sampleDist.getAvgMass() + 3 * PROTON_MASS) / 3;
                            yionRef[0][8] = (refDist.getEndMass() + 3 * PROTON_MASS) / 3 + massTolerance;
                            yionRef[0][7] = (refDist.getStartMass() + 3 * PROTON_MASS) / 3 - massTolerance;
                            yionRef[0][6] = (refDist.getAvgMass() + 3 * PROTON_MASS) / 3;

                        case 2:
                            yionSample[0][5] = (sampleDist.getEndMass() + 2 * PROTON_MASS) / 2 + massTolerance;
                            yionSample[0][4] = (sampleDist.getStartMass() + 2 * PROTON_MASS) / 2 - massTolerance;
                            yionSample[0][3] = (sampleDist.getAvgMass() + 2 * PROTON_MASS) / 2;
                            yionRef[0][5] = (refDist.getEndMass() + 2 * PROTON_MASS) / 2 + massTolerance;
                            yionRef[0][4] = (refDist.getStartMass() + 2 * PROTON_MASS) / 2 - massTolerance;
                            yionRef[0][3] = (refDist.getAvgMass() + 2 * PROTON_MASS) / 2;

                        case 1:
                            yionSample[0][2] = sampleDist.getEndMass() + 1 * PROTON_MASS + massTolerance; //add proton to give b fragment ion
                            yionSample[0][1] = sampleDist.getStartMass() + 1 * PROTON_MASS - massTolerance;
                            yionSample[0][0] = sampleDist.getAvgMass() + 1 * PROTON_MASS;
                            yionRef[0][2] = refDist.getEndMass() + 1 * PROTON_MASS + massTolerance;
                            yionRef[0][1] = refDist.getStartMass() + 1 * PROTON_MASS - massTolerance;
                            yionRef[0][0] = refDist.getAvgMass() + 1 * PROTON_MASS;

                        default:
                            break;

                    }

                    /*
                    for(int i=0;i<bionSample.length;i++)
                    {
                        for(int j=0;j<bionSample[i].length;j++)
                            System.out.print(bionSample[i][j] + " ");

                        System.out.println(" ");
                    }
                    System.out.println(" ");

                    for(int i=0;i<yionSample.length;i++)
                    {
                        for(int j=0;j<yionSample[i].length;j++)
                            System.out.print(yionSample[i][j] + " ");

                        System.out.println(" ");
                    }
                    System.out.println(" ");

                    for(int i=0;i<bionRef.length;i++)
                    {
                        for(int j=0;j<bionRef[i].length;j++)
                            System.out.print(bionRef[i][j] + " ");

                        System.out.println(" ");
                    }
                    System.out.println(" ");

                    for(int i=0;i<yionRef.length;i++)
                    {
                        for(int j=0;j<yionRef[i].length;j++)
                            System.out.print(yionRef[i][j] + " ");

                        System.out.println(" ");
                    }
                    System.out.println(" ");


                    System.exit(1);
                     */
                    double calcSamMass = (sampleDist.getAvgMass() + chargeState * PROTON_MASS) / chargeState;
                    double calcRefMass = (refDist.getAvgMass() + chargeState * PROTON_MASS) / chargeState;

                    conf.setCalcSamAvgMass(calcSamMass);
                    conf.setCalcRefAvgMass(calcRefMass);

                    peptideEle.setAttribute("lightStartMass", String.valueOf(sampleDist.getStartMass()));
                    peptideEle.setAttribute("heavyStartMass", String.valueOf(refDist.getStartMass()));
                    peptideEle.setAttribute("lightAvgMass", String.valueOf(conf.getCalcSamAvgMass()));
                    peptideEle.setAttribute("heavyAvgMass", String.valueOf(conf.getCalcRefAvgMass()));

                    peptideEle = this.createXmlChroPeptideTitle(false, peptide);

                    String ms2FileName = peptide.getFileName().substring(0, peptide.getFileName().indexOf("."));

                    StringBuffer rangeKey = new StringBuffer();
                    rangeKey.append(protein.getLocus());
                    rangeKey.append(ms2FileName);
                    rangeKey.append(peptide.getSequence().substring(2, peptide.getSequence().length() - 2));
                    SpecRange range = rangeGen.getSpecRange(rangeKey.toString());

                    if (null == range) {
                        int tmpScanNum = Integer.parseInt(peptide.getScanNum());
                        range = new SpecRange(tmpScanNum, tmpScanNum);
                        peptideEle.setAttribute("start_scan", peptide.getScanNum());
                        peptideEle.setAttribute("end_scan", peptide.getScanNum());
                    } else {
                        peptideEle.setAttribute("start_scan", String.valueOf(range.getMin()));
                        peptideEle.setAttribute("end_scan", String.valueOf(range.getMax()));
                    }

                    //sb.append("\tStartScan\tEndScan\tDTAPeakStart\tDTAPeakEnd\n");
                    Element chro = new Element("chro");
                    //output.append("[CHROMATOGRAMS]\tSCAN\tSAMPLE\tREFERENCE\n");

                    iFile = ht.get(this.filePath + ms2FileName + "." + "ms2");
                    keys = iFile.getKeys();

                    //samplePrecursor = conf.getPrecursor( (sampleDist.getAvgMass()+chargeState)/chargeState );
                    //refPrecursor = conf.getPrecursor( (refDist.getAvgMass()+chargeState)/chargeState );
                    keyIndex = Arrays.binarySearch(keys, Integer.parseInt(peptide.getScanNum()));

                    if (keyIndex < 0) //Cannot find index
                    {
                        keyIndex = -(++keyIndex); //Math.abs(++keyIndex);
                    }
                    /**
                     * Find start scan number same as a precursor of
                     * sampleAvgMass Then, later the program will find following
                     * scan # quicker.
                     *
                     */
                    /* End of finding start spectrum number */

                    samplePrecursor = 1;
                    refPrecursor = 1;

                    TIntDoubleHashMap precursorMap = iFile.getPrecursorMap();

                    // {
                    //   System.out.println("Outside isowindow for sample precursor " + keys[keyIndex] + " " + (sampleDist.getAvgMass()+chargeState)/chargeState + " " + precursorMap.get(keys[keyIndex]) + " " + samplePrecursor);
                    //}
                    //if(precursorMap.get(keys[keyIndex+diff]) != refPrecursor)
                    //   System.out.println("Outside isowindow for ref precursor " + keys[keyIndex+diff] + " " + (refDist.getAvgMass()+chargeState)/chargeState  + " " + precursorMap.get(keys[keyIndex+diff]) + " " + refPrecursor);
                    //System.out.println("set text" + iFile + " " + range  + " " + keyIndex + " " + diff + " " + bionSample + " " + bionRef + " " + yionSample + " " + yionRef + " " + samplePrecursor + " " + refPrecursor + " " + conf + " " + chargeState);
                    //chro.setText( CalcUtil.calculateMS2Mass( iFile, range, keyIndex, diff, bionSample, bionRef, yionSample, yionRef, samplePrecursor, refPrecursor, conf, chargeState) );
                    try {
                        String outStr = CalcUtil.calculateMS2Mass(iFile, range, keyIndex, bionSample, bionRef, yionSample, yionRef, conf, chargeState);

                        if (outStr == null) {
                            continue;
                        }

                        chro.setText(outStr);
                    } catch (PrecursorNotFoundException ive) {
                        System.out.println("Precursor not found for " + pepSequence);

                        if (null != progress) {
                            ChroProgressDialog.addMessageWithLine(progress, "Error : Precursor not found for " + pepSequence);
                            percent += eachSeg;
                            progress.setProgress((int) percent);
                        }

                        ive.printStackTrace();

                        continue;
                    } catch (Exception e) {
                        System.out.println("Precursor not found for " + pepSequence);

                        if (null != progress) {
                            ChroProgressDialog.addMessageWithLine(progress, "Error : Precursor not found for " + pepSequence);
                            percent += eachSeg;
                            progress.setProgress((int) percent);
                        }
                        e.printStackTrace();

                        continue;
                    }

                    peptideEle.addContent(chro);

                    Element fragEle = new Element("frag");

                    Element bSample = new Element("bs");
                    StringBuffer tempSb = new StringBuffer();
                    for (int i = 0; i < bionSample.length; i++) {
                        int j;
                        for (j = 0; j < bionSample[i].length - 1; j++) {
                            tempSb.append(formatter.format(bionSample[i][j])).append(" ");
                        }

                        tempSb.append(formatter.format(bionSample[i][j])).append(",");
                    }

                    bSample.setText(tempSb.toString());
                    fragEle.addContent(bSample);

                    Element bRef = new Element("br");
                    tempSb.delete(0, tempSb.length());
                    for (int i = 0; i < bionRef.length; i++) {
                        int j;
                        for (j = 0; j < bionRef[i].length - 1; j++) {
                            tempSb.append(formatter.format(bionRef[i][j])).append(" ");
                        }

                        tempSb.append(formatter.format(bionRef[i][j])).append(",");

                    }

                    bRef.setText(tempSb.toString());
                    fragEle.addContent(bRef);

                    Element ySample = new Element("ys");
                    tempSb.delete(0, tempSb.length());
                    for (int i = 0; i < yionSample.length; i++) {
                        int j;
                        for (j = 0; j < yionSample[i].length - 1; j++) {
                            tempSb.append(formatter.format(yionSample[i][j])).append(" ");
                        }

                        tempSb.append(formatter.format(yionSample[i][j])).append(",");
                    }

                    ySample.setText(tempSb.toString());
                    fragEle.addContent(ySample);

                    Element yRef = new Element("yr");
                    tempSb.delete(0, tempSb.length());
                    for (int i = 0; i < yionRef.length; i++) {
                        int j;
                        for (j = 0; j < yionRef[i].length - 1; j++) {
                            tempSb.append(formatter.format(yionRef[i][j])).append(" ");
                        }

                        tempSb.append(formatter.format(yionRef[i][j])).append(",");
                    }

                    yRef.setText(tempSb.toString());
                    fragEle.addContent(yRef);

                    peptideEle.addContent(fragEle);

                    percent += eachSeg;

                    if (null != progress) {
                        progress.setProgress((int) percent);
                    }

                    System.out.print(pepCount);
                    System.out.print("/");
                    System.out.print(redundantPeptideNum);
                    System.out.print(" peptides, ");
                    System.out.print((int) percent);
                    //System.out.print(" % is complete\n");
                    System.out.print(" % is complete\r");

                    proteinEle.addContent(peptideEle);
                }

                if (proteinEle.getChildren().size() > 0) {
                    rootEle.addContent(proteinEle);
                }
            }

            Document doc = new Document(rootEle);
            OutputStream os = new FileOutputStream(filePath + "census_chro.xml");
            XMLOutputter outputter = new XMLOutputter();
            outputter.setFormat(Format.getPrettyFormat());
            outputter.output(doc, os);
            os.close();

            System.out.println("\n100% complete");
            //System.out.println( System.out.println"\n100% complete");

            //Create relex.chro file
            /*
            out = new BufferedOutputStream(new FileOutputStream(filePath + "relex.chro"));
            p = new PrintStream(out);
            p.print(output.toString());
             */
        } catch (IOException e) {
            System.out.println("IO Error while generating msms chro file : " + pepSequence + e);
            e.printStackTrace();
            throw new IOException(e.toString());
        } catch (Exception e) {
            System.out.println("Error while generating msms chro file : " + pepSequence + e);
            e.printStackTrace();
            throw new Exception(e.toString());
        } finally {
            if (null != p) {
                p.close();
            }

            if (null != out) {
                out.close();
            }

            //Close all random files
            for (Enumeration e = ht.keys(); e.hasMoreElements();) {
                iFile = ht.get(e.nextElement());

                if (null != iFile) {
                    iFile.close();
                }
            }

        }

    }

    //cbamberg approach
    public void createMsmsQautnXmlChro(ChroProgressDialog progress) throws IOException, Exception {

        //  System.out.println("addddddddd"); index
        this.progress = progress;
        int[] keys;

        this.filePath = conf.getFilePath();

        Hashtable<String, IndexedFile> ht = createIndexedFiles(filePath, CensusConstants.MS2_FILE);

        IndexedFile iFile;
        BufferedOutputStream out = null;
        PrintStream p = null;

        long startTime = System.currentTimeMillis();
        String pepSequence = null;

        try {

            /**
             * ****************************************************************
             * Read DTASelect.txt file to find spectrum range for each peptide
             /*****************************************************************
             */
            ChroProgressDialog.addMessageWithLine(progress, "");

            IsotopeReader isoReader = null;

            if (conf.isXmlConf()) {
                isoReader = new IsotopeReader(conf.getRootConfEle());
            } else {
                isoReader = new IsotopeReader(isotopeFile);
            }

            IdentificationReader idReader = BaseIdentificationReader.getIdentificationInst(isoReader);
            SpecRangeGenerator rangeGen = new SpecRangeGenerator();
            int redundantPeptideNum = idReader.getTotalPeptideNumber();

            IsotopeTable<String, int[]> isoTable = isoReader.getIsotope();
            int[] sampleNterm = isoTable.get("sampleNTERM");
            int[] sampleCterm = isoTable.get("sampleCTERM");
            int[] refNterm = isoTable.get("refNTERM");
            int[] refCterm = isoTable.get("refCTERM");

            //increase status bar
            double percent = 0.0;
            double eachSeg = 100.0 / redundantPeptideNum;
            int pepCount = 0;

            //IsotopeDist sampleDist;
            //IsotopeDist refDist;
            IsotopeDist sampleDist;
            IsotopeDist refDist;
            Protein protein;
            Peptide peptide;

            Element rootEle = this.createXmlChroHeader(2, conf.getExpType());

            ElementComposition belement;
            ElementComposition yelement;
            ElementComposition totalElement;

            //int[] elementSampleArr;
            int keyIndex;
            int start;
            int last;

            double samplePrecursor;
            double refPrecursor;

            double[][] bionSample;
            double[][] bionRef;
            double[][] yionSample;
            double[][] yionRef;

            Element proteinEle = null;

            HashSet<String> hcdScanPool = new HashSet<String>();

            for (Iterator<Protein> itr = idReader.getProteins(); itr.hasNext();) {
                protein = itr.next();

                for (Iterator<Peptide> pepItr = protein.getPeptides(); pepItr.hasNext();) {
                    peptide = pepItr.next();
                    iFile = ht.get(this.filePath + peptide.getFileName() + "." + "ms2");

                    if (null != iFile) {
                        String scanType = iFile.getScanType(peptide.getScanNumber());
                        if ("HCD".equals(scanType)) {
                            hcdScanPool.add(peptide.getFileName() + peptide.getScanNumber());
                        }

                    }
                }
            }

            idReader = BaseIdentificationReader.getIdentificationInst(isoReader);
            ArrayList<Protein> aList = new ArrayList<Protein>();

            for (Iterator<Protein> itr = idReader.getProteins(); itr.hasNext();) {
                protein = itr.next();

                proteinEle = new Element("protein");
                proteinEle.setAttribute("locus", protein.getLocus());
                proteinEle.setAttribute("seq_ct", protein.getSeqCount());
                proteinEle.setAttribute("spec_ct", protein.getSpectrumCount());
                proteinEle.setAttribute("seq_cov", protein.getSeqCoverage());
                proteinEle.setAttribute("length", protein.getLength());
                proteinEle.setAttribute("molwt", protein.getMolWt());
                proteinEle.setAttribute("pi", protein.getPI());
                proteinEle.setAttribute("val", protein.getValidation());

                try {
                    proteinEle.setAttribute("desc", protein.getDescription());
                } catch (org.jdom.IllegalDataException ide) {
                    proteinEle.setAttribute("desc", StringUtil.removeIsoControlChar(protein.getDescription()));
                }
                String label = conf.getMs2Label();
                for (Iterator<Peptide> pepItr = protein.getPeptides(); pepItr.hasNext();) {
                    peptide = pepItr.next();
                    pepCount++;

                    //if (!"".equals(label) && !(peptide.getFileName().startsWith(label))) {
                    if (null != label && !"".equals(label) && !(peptide.getFileName().startsWith(label))) {
                        continue;
                    } else {
                        String ms2File = conf.getFilePath()+ peptide.getFileName().split("\\.")[0] + ".ms2";
                        //String lightfile = path + "/" + file.substring(1, file.length());
                        File lf = new File(ms2File);

                        if (!(lf.exists())) {
                            System.out.println(ms2File+" not present in the folder");
                            continue;

                        }
                        pepSequence = peptide.getSequence();

                        Element peptideEle = new Element("peptide");
                        peptideEle.setAttribute("file", peptide.getFileName());
                        peptideEle.setAttribute("scan", peptide.getScanNum());
                        peptideEle.setAttribute("seq", peptide.getSequence());

                        //What is the purpose of hs???
                        char[] ch = pepSequence.substring(2, peptide.getSequence().length() - 2).toCharArray();
                        char[] revCh = pepSequence.substring(2, peptide.getSequence().length() - 2).toCharArray();
                        ArrayUtils.reverse(revCh);

                        try {
                            totalElement = new ElementComposition(ch, 0, ch.length, isoTable);
                            totalElement.calculate();
                        } catch (InvalidAAException ive) {
                            System.out.println("Not Quantifiable peptide : " + pepSequence);

                            percent += eachSeg;
                            if (null != progress) {
                                ChroProgressDialog.addMessageWithLine(progress, "Not Quantifiable peptide : " + pepSequence);
                                progress.setProgress((int) percent);

                            }
                            continue;
                        }

                        int chargeState = Integer.parseInt(peptide.getChargeState());

                        int pepLength = 0;

                        for (int i = 0; i < ch.length; i++) {
                            if (ch[i] == '*' || ch[i] == '@' || ch[i] == '#') {
                                continue;
                            }

                            pepLength++;
                        }

                        bionSample = new double[pepLength][1];
                        bionRef = new double[pepLength][1];
                        //Yions
                        yionSample = new double[pepLength][1];
                        yionRef = new double[pepLength][1];

                        int pepIndex = 0;

                        //System.out.println("aamass" + massTolerance + " " + conf.getMassTolerance());
                        for (int i = 0; i < ch.length; i++) {
                            if (ch[i] == '*' || ch[i] == '@' || ch[i] == '#') {
                                continue;
                            }

                            try {
                                belement = new ElementComposition(ch, 0, i + 1, isoTable);
                                belement.calculate();
                                belement.calculateBion();

                                //int[] aa = belement.getElementSampleArr();
                                //for(int ii:aa) System.out.print(ii + " ");
                                //System.out.println("");
                                //aa = belement.getElementRefArr();
                                // for(int ii:aa) System.out.print(ii + " ");
                                //System.out.println("");
                                yelement = new ElementComposition(revCh, 0, i + 1, isoTable);
                                yelement.calculateYion();

                            } catch (InvalidAAException ive) {
                                System.out.println("Not Quantifiable peptide : " + pepSequence);

                                percent += eachSeg;
                                if (null != progress) {
                                    ChroProgressDialog.addMessageWithLine(progress, "Not Quantifiable peptide : " + pepSequence);
                                    progress.setProgress((int) percent);

                                }
                                continue;
                            }

                            //System.exit(0);
                            //  System.out.println("1==============");
                            //sampleDist = new IsotopeDist(yelement.getElementSampleArr(), yelement.getModShift(), true);
                            sampleDist = new IsotopeDist(); //yelement.getElementSampleArr(), yelement.getModShift(), true);
                            sampleDist.setElement(yelement.getElementSampleArr());
                            sampleDist.setModShift(yelement.getModShift());
                            sampleDist.calculateLightSimple();
                            double[] marr = sampleDist.getHighMassList();

                            //refDist = new IsotopeDist(yelement.getElementRefArr(), yelement.getModShift(), false);
                            refDist = new IsotopeDist(); //yelement.getElementRefArr(), yelement.getModShift(), false);
                            refDist.setElement(yelement.getElementRefArr());
                            refDist.setModShift(yelement.getModShift());
                            refDist.calculateHeavySimple();

                            double[] rmarr = refDist.getHighMassList();

                            for (int j = 0; j < 1; j++) {
                                if (j >= marr.length) {
                                    break;
                                }
                                //yionSample[(pepIndex+1)%pepLength][j] = marr[j];
                                yionSample[pepIndex][j] = marr[j];
                            }

                            for (int j = 0; j < 1; j++) {
                                if (j >= rmarr.length) {
                                    break;
                                }
                                yionRef[pepIndex][j] = rmarr[j];
                                //System.out.println("== " + rmarr[j]);
                            }

                            //sampleDist = new IsotopeDist(belement.getElementSampleArr(), belement.getModShift(), true);
                            //refDist = new IsotopeDist(belement.getElementRefArr(), belement.getModShift(), false);
                            sampleDist = new IsotopeDist(); //yelement.getElementSampleArr(), yelement.getModShift(), true);
                            sampleDist.setElement(belement.getElementSampleArr());
                            sampleDist.setModShift(belement.getModShift());
                            sampleDist.calculateLightSimple();
                            marr = sampleDist.getHighMassList();

//                        System.out.println("bbbbbbbbb");
                            //refDist = new IsotopeDist(yelement.getElementRefArr(), yelement.getModShift(), false);
                            refDist = new IsotopeDist(); //yelement.getElementRefArr(), yelement.getModShift(), false);
                            refDist.setElement(belement.getElementRefArr());
                            refDist.setModShift(belement.getModShift());
                            refDist.calculateHeavySimple();

                            rmarr = refDist.getHighMassList();

                            /*
                        for(double a:marr)
                            System.out.print(a + " ");
                        System.out.println(" ");
                        System.out.println("1==============");


                        for(double a:rmarr)
                            System.out.print(a + " ");
                        System.out.println(" ");
                        System.out.println("1==============");
/*
                         //System.exit(0);

                       // System.out.println(" ");
                       // System.out.println("1==============");
                        //arr1 = yelement.getElementRefArr();
                        //for(int a:arr1)
                        //for(double a:rmarr)
                        //    System.out.print(a + " ");

                        //System.out.println(" ");
                        //System.out.println("1==============");

                             */
                            //for(int j=0;j<3;j++) {
                            for (int j = 0; j < 1; j++) {
                                if (j >= marr.length) {
                                    break;
                                }
//                            System.out.println(marr.length + " " + j + " " + marr[j]);
                                bionSample[pepIndex][j] = marr[j]; // - 0.00054857;
                                //System.out.println("mm" + marr[j]);
                            }

                            for (int j = 0; j < 1; j++) {
                                if (j >= rmarr.length) {
                                    break;
                                }
                                bionRef[pepIndex][j] = rmarr[j]; // - 0.00054857;
                                //System.out.println("rm" + rmarr[j]);
                            }
                            pepIndex++;
                        }

                        peptideEle = this.createXmlChroPeptideTitle(false, peptide);

                        String ms2FileName = peptide.getFileName();

                        StringBuffer rangeKey = new StringBuffer();
                        rangeKey.append(protein.getLocus());
                        rangeKey.append(ms2FileName);
                        rangeKey.append(peptide.getSequence().substring(2, peptide.getSequence().length() - 2));

                        Element chro = new Element("chro");
                        //output.append("[CHROMATOGRAMS]\tSCAN\tSAMPLE\tREFERENCE\n");

                        iFile = ht.get(this.filePath + ms2FileName + "." + "ms2");

                        if (iFile == null) {
                            System.out.println("Error: cannot find spectral file " + this.filePath + ms2FileName + "." + "ms2");
                        }

                        String sType = iFile.getScanType(peptide.getScanNumber());
                        keys = iFile.getKeys();

                        if ("CID".equals(sType)) {
                            boolean samePrec = TMTUtil.isSamePrecursorWithFilename(iFile, conf.getScanShift(), peptide.getScanNumber(), hcdScanPool, peptide.getFileName());
                            if (samePrec) {
                                continue;
                            }
                        }

                        keyIndex = Arrays.binarySearch(keys, Integer.parseInt(peptide.getScanNum()));

                        if (keyIndex < 0) //Cannot find index
                        {
                            keyIndex = -(++keyIndex); //Math.abs(++keyIndex);
                        }
                        samplePrecursor = 1;
                        refPrecursor = 1;

                        //System.exit(0);
                        try {
                            chro.setText("NA");
                        } catch (Exception e) {
                            System.out.println("Precursor not found for " + pepSequence);

                            if (null != progress) {
                                ChroProgressDialog.addMessageWithLine(progress, "Error : Precursor not found for " + pepSequence);
                                percent += eachSeg;
                                progress.setProgress((int) percent);
                            }
                            e.printStackTrace();

                            continue;
                        }

                        peptideEle.addContent(chro);

                        Element fragEle = new Element("frag");

                        Element bSample = new Element("bs");
                        StringBuffer tempSb = new StringBuffer();
                        //keyIndex += conf.getScanShift();

                        String bStr = CalcUtil.calculateSingleMS2Mass(iFile, keyIndex, bionSample, conf, chargeState, bionRef, peptide.getScanNumber());
                        if (null != bStr) {
                            tempSb.append(bStr);
                        }

                        bSample.setText(tempSb.toString());
                        fragEle.addContent(bSample);

                        Element bRef = new Element("br");
                        tempSb.delete(0, tempSb.length());
                        String brStr = CalcUtil.calculateSingleMS2Mass(iFile, keyIndex, bionRef, conf, chargeState, bionSample, peptide.getScanNumber());
                        if (null != brStr) {
                            tempSb.append(brStr);
                        }

                        bRef.setText(tempSb.toString());
                        fragEle.addContent(bRef);

                        Element ySample = new Element("ys");
                        tempSb.delete(0, tempSb.length());

                        String yStr = CalcUtil.calculateSingleMS2Mass(iFile, keyIndex, yionSample, conf, chargeState, yionRef, peptide.getScanNumber());
                        if (null != yStr) {
                            tempSb.append(yStr);
                        }

                        // System.out.println("==>>" + tempSb.toString());
                        ySample.setText(tempSb.toString());
                        fragEle.addContent(ySample);

                        Element yRef = new Element("yr");
                        tempSb.delete(0, tempSb.length());

                        String yrStr = CalcUtil.calculateSingleMS2Mass(iFile, keyIndex, yionRef, conf, chargeState, yionSample, peptide.getScanNumber());
                        if (null != yrStr) {
                            tempSb.append(yrStr);
                        }

                        //System.out.println(tempSb.toString());
                        yRef.setText(tempSb.toString());
                        fragEle.addContent(yRef);

                        peptideEle.addContent(fragEle);

                        percent += eachSeg;

                        if (null != progress) {
                            progress.setProgress((int) percent);
                        }

                        System.out.print(pepCount);
                        System.out.print("/");
                        System.out.print(redundantPeptideNum);
                        System.out.print(" peptides, ");
                        System.out.print((int) percent);
                        //System.out.print(" % is complete\n");
                        System.out.print(" % is complete\r");

                        proteinEle.addContent(peptideEle);
                    }
                }

                //if(proteinEle.getChildren().size()>0)
                //    rootEle.addContent(proteinEle);
                if (proteinEle.getChildren().size() > 0) {

                    Element redunEle = new Element("redundant");

                    int rPepCount = 0;

                    for (Iterator<Protein> rItr = aList.iterator(); rItr.hasNext();) {
                        rPepCount++;
                        Protein rp = rItr.next();

                        Element rpEle = new Element("protein");
                        rpEle.setAttribute("locus", rp.getLocus());
                        rpEle.setAttribute("seq_ct", rp.getSeqCount());
                        rpEle.setAttribute("spec_ct", rp.getSpectrumCount());
                        rpEle.setAttribute("seq_cov", rp.getSeqCoverage());
                        rpEle.setAttribute("length", rp.getLength());
                        rpEle.setAttribute("molwt", rp.getMolWt());
                        rpEle.setAttribute("pi", rp.getPI());
                        rpEle.setAttribute("val", rp.getValidation());

                        try {
                            rpEle.setAttribute("desc", rp.getDescription());
                        } catch (org.jdom.IllegalDataException ide) {
                            rpEle.setAttribute("desc", StringUtil.removeIsoControlChar(rp.getDescription()));
                        }

                        redunEle.addContent(rpEle);
                    }

                    if (rPepCount > 0) {
                        proteinEle.addContent(redunEle);
                    }

                    rootEle.addContent(proteinEle);
                    aList.clear();
                } else {
                    aList.add(protein);
                    //	    set.add(protein.getLocus());
                }

            }

            Document doc = new Document(rootEle);
            OutputStream os = new FileOutputStream(filePath + "census_chro.xml");
            XMLOutputter outputter = new XMLOutputter();
            outputter.setFormat(Format.getPrettyFormat());
            outputter.output(doc, os);
            os.close();

            System.out.println("\n100% complete");
            //System.out.println( System.out.println"\n100% complete");

            //Create relex.chro file
            /*
            out = new BufferedOutputStream(new FileOutputStream(filePath + "relex.chro"));
            p = new PrintStream(out);
            p.print(output.toString());
             */
        } catch (IOException e) {
            System.out.println("IO Error while generating msms chro file : " + pepSequence + e);
            e.printStackTrace();
            throw new IOException(e.toString());
        } catch (Exception e) {
            System.out.println("Error while generating msms chro file : " + pepSequence + e);
            e.printStackTrace();
            throw new Exception(e.toString());
        } finally {
            if (null != p) {
                p.close();
            }

            if (null != out) {
                out.close();
            }

            //Close all random files
            for (Enumeration e = ht.keys(); e.hasMoreElements();) {
                iFile = ht.get(e.nextElement());

                if (null != iFile) {
                    iFile.close();
                }
            }

        }

    }

    /*
     * Find first scan number matching avgMass
     */
    private int findStartScanNum(int scanNum, RandomAccessFile file, int[] keys, double precursor) throws IOException {
        System.out.println("scan" + scanNum);
        System.out.println("scan" + keys[scanNum]);
        file.seek(keys[scanNum]);
        String str = file.readLine();
        System.out.println("sss" + str);

        double d = Double.parseDouble(str.substring(str.lastIndexOf("\t")));

        System.out.println(d + "\t" + precursor);

        while (d != precursor) {
            System.out.println(d + "\t" + precursor);

            file.seek(keys[++scanNum]);
            str = file.readLine();
            d = Double.parseDouble(str.substring(str.lastIndexOf("\t")));
        }

        return 0;
    }

    private double findScanNum(RandomAccessFile file, int scanNum, long pos, double avgMass) throws IOException {
        file.seek(pos);
        String str = file.readLine();

        double d = Double.parseDouble(str.substring(str.lastIndexOf("\t")));

        return 0;
    }

    public static String cleanFileName(String fileName) {
        // this is a problem of mzxml data format inconsistency
        if (null != fileName && fileName.contains("rawFile")) {
            fileName = fileName.substring(fileName.indexOf("rawFile")).trim();
            fileName = fileName.substring(fileName.indexOf(":") + 1).trim();
            fileName = fileName.substring(0, fileName.lastIndexOf("."));
        } else {
            int index = fileName.indexOf(".");

            if (index > 0) {
                fileName = fileName.substring(0, index);
            }

            //for(int i=0;i<3;i++)
            //	fileName = fileName.substring(0, fileName.lastIndexOf("."));
//	    fileName = fileName.substring(0, fileName.lastIndexOf("."));
        }

        return fileName;
    }

    public void runMRMWithoutId() throws IOException, CensusGeneralException, Exception {

        this.conf = Configuration.getInstance();
        this.filePath = conf.getFilePath();
        Element mrmParamsEle = conf.getRootConfEle().getChild("params").getChild("mrm_params");

        List<scripts.mrm.PeptideModel> pepList = new ArrayList<scripts.mrm.PeptideModel>();

        for (Iterator<Element> itr = mrmParamsEle.getChildren("peptide").iterator(); itr.hasNext();) {
            Element pepEle = itr.next();

            //System.out.println(pgroupEle);
            scripts.mrm.PeptideModel peptide = new scripts.mrm.PeptideModel(pepEle.getAttributeValue("name"));

            for (Iterator<Element> pepitr = pepEle.getChildren("precursor").iterator(); pepitr.hasNext();) {
                Element preEle = pepitr.next();

                double pmass = Double.parseDouble(preEle.getAttributeValue("mass"));

                scripts.mrm.Precursor p = new scripts.mrm.Precursor(pmass);

                for (Iterator<Element> ditr = preEle.getChildren("transition").iterator(); ditr.hasNext();) {
                    Element dauEle = ditr.next();
                    double dmass = Double.parseDouble(dauEle.getAttributeValue("mass"));

                    p.addDaughter(dmass);
                }

                peptide.addPrecursor(p);
            }

            pepList.add(peptide);

        }

        String ms2File = mrmParamsEle.getChildText("ms2_file");
        //SpectrumReader sr = new SpectrumReader("/home/rpark/rpark_on_data/project/jcopping/mrm_012910/CF301hrrep1-01.ms2", "ms2");
        SpectrumReader sr = new SpectrumReader(ms2File, "ms2");
        Hline h = new Hline(sr.getHlines());

        Iterator<PeakList> it = sr.getSpectra();
        Hashtable<Double, scripts.mrm.Precursor> pht = new Hashtable<Double, scripts.mrm.Precursor>();

        for (Iterator<scripts.mrm.PeptideModel> itr = pepList.iterator(); itr.hasNext();) {
            scripts.mrm.PeptideModel pmodel = itr.next();

            for (Iterator<scripts.mrm.Precursor> pitr = pmodel.getPlist().iterator(); pitr.hasNext();) {
                scripts.mrm.Precursor precursor = pitr.next();

                pht.put(precursor.getMass(), precursor);
            }

        }

//        int counter = 0;
//        int numPeaks = 0;
        //boolean sortByIntensity = true;
        while (it.hasNext()) {
            PeakList list = it.next();

            StringBuffer sb = new StringBuffer();

            scripts.mrm.Precursor precursor = pht.get(list.getPrecursorMass());

            if (null == precursor) {
                continue;
            }

            //precursor.addIntensity(p.getM2z(), p.getIntensity());
            precursor.addIntensity(list.getPeakList(), list.getLoscan(), list.getRetentionTime());
            /*
            for(Iterator<Peak> itr=list.getPeaks(); itr.hasNext(); )
            {
                Peak p = itr.next();

                precursor.addIntensity(p.getM2z(), p.getIntensity(), list.getLoscan());

		//System.out.println("====" + list.getPrecursorMass() + " " + list.getLoscan() + " " + p.getM2z() + " " +  p.getIntensity());
            } */

        }

        // peptide.print(0);
        //peptide.calculateRatio();
        //peptide.print();
        Element rootEle = this.createXmlChroHeader(1);

        for (Iterator<PeptideModel> itr = pepList.iterator(); itr.hasNext();) {
            PeptideModel peptide = itr.next(); //pepList.get(0);

            Element proteinEle = new Element("protein");
            proteinEle.setAttribute("locus", peptide.getName());
            proteinEle.setAttribute("seq_ct", "");
            proteinEle.setAttribute("spec_ct", "");
            proteinEle.setAttribute("seq_cov", "");
            proteinEle.setAttribute("length", "");
            proteinEle.setAttribute("molwt", "");
            proteinEle.setAttribute("pi", "");
            proteinEle.setAttribute("val", "");

            List<Precursor> preList = peptide.getPlist();
            for (Iterator<Precursor> pitr = preList.iterator(); pitr.hasNext();) {
                Precursor precursor = pitr.next();

                int peakIndex = precursor.findpeak();
                //daughterSize += pre.getDmassList().size();

                for (int i = 0; i < precursor.getDmassList().size(); i++) {
                    Daughter d = precursor.getDaughter(i);
                    d.findPeakArea(peakIndex);

                    Element peptideEle = new Element("peptide");
                    peptideEle.setAttribute("unique", "*");
                    peptideEle.setAttribute("file", String.valueOf(d.getRt(peakIndex)));
                    peptideEle.setAttribute("scan", String.valueOf(d.getScan(peakIndex)));
                    peptideEle.setAttribute("seq", String.valueOf(d.getMass()));
                    peptideEle.setAttribute("xcorr", "0");
                    peptideEle.setAttribute("deltaCN", "0");
                    peptideEle.setAttribute("charge", "0");
                    peptideEle.setAttribute("spC", "0");
                    peptideEle.setAttribute("start_scan", String.valueOf(d.getStartScan()));
                    peptideEle.setAttribute("end_scan", String.valueOf(d.getEndScan()));
                    proteinEle.addContent(peptideEle);

//		    System.out.println("--------" + d.getScan(peakIndex) + " " + d.getStartScan() + " " + d.getEndScan());
                    Element chro = new Element("chro");

///		    System.out.println("== " + peakIndex  + " " + d.getMass() + " " + d.getScan(peakIndex) + " " + d.getIntensity(peakIndex));
                    StringBuffer sb = new StringBuffer();
                    int range = 150;
                    sb.append("P ").append(d.getStartScan()).append(" ").append(d.getEndScan()).append(";");
                    for (int ii = peakIndex - range; ii < peakIndex + range; ii++) {
                        sb.append(d.getScan(ii)).append(" ").append((int) d.getIntensity(ii));
                        for (int j = 0; j < 7; j++) {
                            sb.append(" 0");
                        }

                        sb.append(";");
                    }

                    chro.setText(sb.toString());

                    //chro.setText( CalcUtil.calculateFullMS( keyIndex, iFile, sampleStartMass, sampleEndMass, refStartMass, refEndMass, range) );
                    //<chro>P 7651 7779;7532 301798 0 6 7 2 0 0.02081426666666175 -1.0;7538 405049 0 6 7 1 0 0.01874126666666598 -1.0;7544 0 0 6 7 0 0 -1.0 -1.0;7550 143416 0 6 7 1 0 0.019541266666692536 -1.0;7556 101896 0 6 7 1 0 0.01864126666669108 -1.0;7561 0 0 6 7 0 0 -1.0 -1.0;7567 0 0 6 7 0 0 -1.0 -1.0;7573 0 0 6 7 0 0 -1.0 -1.0;7579 0 0 6 7 0 0 -1.0 -1.0;7585 0 0 6 7 0 0 -1.0 -1.0;7591 88080 0 6 7 1 0 0.023107399999958034 -1.0;7597 0 0 6 7 0 0 -1.0 -1.0;7603 0 0 6 7 0 0 -1.0 -1.0;7609 0 0 6 7 0 0 -1.0 -1.0;7615 0 0 6 7 0 0 -1.0 -1.0;7621 0 0 6 7 0 0 -1.0 -1.0;7627 0 0 6 7 0 0 -1.0 -1.0;7633 0 0 6 7 0 0 -1.0 -1.0;7639 0 0 6 7 0 0 -1.0 -1.0;7645 71517 138077 6 7 1 1 0.024858733333303462 4.707333332589769E-4;7651 0 0 6 7 0 0 -1.0 -1.0;7657 0 0 6 7 0 0 -1.0 -1.0;7663 0 0 6 7 0 0 -1.0 -1.0;7669 0 0 6 7 0 0 -1.0 -1.0;7675 0 0 6 7 0 0 -1.0 -1.0;7681 149706 0 6 7 2 0 0.004214266666679123 -1.0;7687 85787 99826 6 7 1 1 0.0016439333332982642 9.707333332471535E-4;7692 0 567381 6 7 0 3 -1.0 8.415777777296777E-4;7698 821825 744934 6 7 3 4 0.0023592666667203352 0.0033786000000759486;7704 1235302 1168124 6 7 6 5 0.00200722222219459 0.004089266666687763;7710 1664041 977082 6 7 4 4 0.0018719666666697776 0.0021536000000423883;7716 1199666 962760 6 7 5 3 0.0018212666667295706 0.0019292666667449037;7722 1139038 243499 6 7 4 2
                    peptideEle.addContent(chro);

                }

            }

            rootEle.addContent(proteinEle);

            //Element peptideEle=null;
            //peptide.findpeak();
        }

        Document doc = new Document(rootEle);
        OutputStream os = new FileOutputStream(filePath + File.separator + "census_chro.xml");
        XMLOutputter outputter = new XMLOutputter();
        outputter.setFormat(Format.getPrettyFormat());
        outputter.output(doc, os);
        os.close();
        System.out.println("\n100% complete");

    }

    public double getMedianN15Ratio(IdentificationReader idReader, IsotopeReader isoReader, Hashtable<String, IndexedFile> ht)
            throws IOException, CensusGeneralException, Exception {

        System.out.println("Calculating 15N enrichment...");
        int totalPepCount = idReader.getTotalPeptideNumber();

        Iterator<Protein> pitr = idReader.getProteins(); //need to run to calculate redundnat peptides
        int[] beanArr = new int[101];

        int count = 0;

        for (Iterator<Protein> itr = pitr; itr.hasNext();) {
            Protein protein = itr.next();

            for (Iterator<Peptide> pepItr = protein.getPeptides(); pepItr.hasNext();) {
                Peptide peptide = pepItr.next();
                count++;
                double progress = (double) count / totalPepCount * 100;
                System.out.print(Formatter.formatDecimal(progress));
                System.out.print(" % is complete\r");

                String pepSequence = peptide.getSequence();
                char[] ch = pepSequence.substring(2, peptide.getSequence().length() - 2).toCharArray();
                ElementComposition element = new ElementComposition(ch, 0, ch.length, isoReader.getIsotope());

                try {
                    element.calculate();
                } catch (InvalidAAException invE) {
                    continue;
                }

                String fileName = peptide.getFileName();
                fileName = cleanFileName(fileName);

                switch (conf.getSpectrumFormat()) {
                    case Configuration.MS_FILE_FORMAT:
                        fileName += ".ms1";
                        break;

                    case Configuration.MZXML_FILE_FORMAT:
                        fileName += ".mzXML";
                        break;

                    default:
                        break;
                }

                int scanNum = Integer.parseInt(peptide.getScanNum());
                IndexedFile iFile = ht.get(filePath + fileName);

                if (null == iFile) {
                    continue;
                }

                int chargeState = Integer.parseInt(peptide.getChargeState());

                IsotopeDist refDist = null;

                double startEnrich = conf.getStartEnrich();
                double endEnrich = conf.getEndEnrich();
                double enrichmentMaxDeviation = conf.getEnrichmentMaxDeviation();

                refDist = new edu.scripps.pms.census.util.N15EnrichmentCalc(element.getElementRefArr(), element.getModShift(), -1, startEnrich, endEnrich, enrichmentMaxDeviation, iFile, scanNum, chargeState);
                //System.out.println("=====" + refDist.getEnrichCorr());
                //System.out.println("=====" + refDist.getFixedEnrichRatio());

                int beanIndex = (int) (refDist.getBestEnrichRatio() * 100);
                if (beanIndex < 0) {
                    beanIndex = 0;
                }

                beanArr[beanIndex] += 1;

            }
        }

        int medianIndex = -1;
        int medianValue = -1;
        for (int i = 0; i < beanArr.length; i++) {
            int temp = beanArr[i];
            if (temp > medianValue) {
                medianIndex = i;
                medianValue = temp;
            }
        }

        System.out.println("");
        System.out.println("done.");

        return (double) medianIndex * 0.01;
    }

    //pulse labeling labeled label
    //aha tev tag
    public void createPulseChro() throws IOException, CensusGeneralException, Exception {

        System.out.println("running pulse labeling...");
       // int pulseCount = conf.getPulseLabelingCount();
        Gson gson = new Gson();

        final int NO_LABEL = 0;
        final int LIGHT_LABEL=1;
        final int HEAVY_LABEL=2;
        int labelingType = NO_LABEL;

        double pulseLightMassDouble = Double.parseDouble(conf.getPulseLightMass());
        double pulseHeavyMassDouble = Double.parseDouble(conf.getPulseHeavyMass());
        double massDiff = pulseHeavyMassDouble - pulseLightMassDouble;

        String pulseLightResidueMass = conf.getPulseResidue();
        String pulseLightResidueMassThreeDigit = conf.getPulseResidue();

        if(!"0".equals(conf.getPulseLightMass())) {
            pulseLightResidueMass += "(" + conf.getPulseLightMass() + ")";
            pulseLightResidueMassThreeDigit += "(" + Formatter.halfRoundUpThreeDigitFormat(conf.getPulseLightMass())  + ")";
        }

        String pulseHeavyResidueMass = conf.getPulseResidue() + "(" + conf.getPulseHeavyMass() + ")";
        String pulseHeavyResidueMassThreeDigit = conf.getPulseResidue() + "(" + Formatter.halfRoundUpThreeDigitFormat(conf.getPulseHeavyMass()) + ")";


        int[] keys;
        this.filePath = conf.getFilePath();

        Hashtable<String, IndexedFile> ht = null;

        switch (conf.getSpectrumFormat()) {
            case Configuration.MS_FILE_FORMAT:
                ht = createIndexedFiles(filePath, CensusConstants.MS1_FILE);
                break;

            case Configuration.MZXML_FILE_FORMAT:
                ht = createIndexedFiles(filePath, CensusConstants.MZXML);
                break;

            default:
                break;
        }

        if (!filePath.endsWith(File.separator)) {
            filePath += File.separator;
        }

        IndexedFile iFile = null;
        BufferedOutputStream out = null;
        PrintStream p = null;

        long startTime = System.currentTimeMillis();

        try {
            /**
             * ****************************************************************
             * Read DTASelect.txt file to find spectrum range for each peptide
		     *****************************************************************
             */
            //TIntLongHashMap index;
            System.out.print("Parsing " + conf.getIdFileName() + "...");

            if (null != progress) {
                progress.addMessage("\nParsing Identified Peptides...");
            }

            if (null != progress) {
                progress.addMessage("\ndone.");
            }

            IsotopeReader isoReader = null;

            if (null != isotopeFile) {
                isoReader = new IsotopeReader(isotopeFile);
            } else {
                isoReader = new IsotopeReader(this.confRootEle);
            }

            File JSONfile =new File(filePath+MASTER_JSON);
            String xyvaluePath = filePath+MASTER_JSON+File.separatorChar+ XY_VALUES+File.separatorChar;
            File xyValuesDir = new File(filePath+MASTER_JSON+File.separatorChar+ XY_VALUES);
            if(!JSONfile.exists())
            {
                JSONfile.mkdir();
            }
            if(!xyValuesDir.exists())
            {
                xyValuesDir.mkdir();
            }

            IdentificationReader idReader = BaseIdentificationReader.getIdentificationInst(isoReader);
            System.out.println("done.");

            SpecRangeGenerator rangeGen = null;

            File dtaFile = new File(filePath + "DTASelect.txt");

            if (conf.isPrintLog()) {
                p = new PrintStream(new FileOutputStream(filePath + "progress.log"));
            }

            //System.out.println(filePath + " " + dtaFile + " " + dtaFile.exists());
            if (dtaFile.exists()) {
                if (conf.getIdFileName().endsWith("txt")) {
                    System.out.print("Parsing DTASelect.txt...");
                    if (null != progress) {
                        progress.addMessage("\nParsing DTASelect.txt...");
                    }
                }

                rangeGen = SpecRangeGenerator.getSpecRangeGenerator(idReader);
                //            SpecRangeGenerator rangeGen = new SpecRangeGenerator(this.filePath + this.dtaSelectFile, idReader.isVersion2(), idReader.getConfidence());
                //SpecRangeGenerator rangeGen = new SpecRangeGenerator(this.filePath + this.dtaSelectFile, idReader);
                //SpecRangeGenerator rangeGen = null;

                if (conf.getIdFileName().endsWith("txt")) {
                    System.out.println("done.");
                    if (null != progress) {
                        progress.addMessage("\ndone.");
                    }
                }
            } else {
                rangeGen = new SpecRangeGenerator();

            }
            Element rootEle = this.createXmlChroHeader(1, conf.getQuantType());
            ElementComposition element;
            ElementComposition totalElement;

            int keyIndex = -1;
            //int start;
            //int last;

            double sampleStartMass;
            double sampleEndMass;
            double refStartMass;
            double refEndMass;

            Element proteinEle = null;
            Element peptideEle = null;

            Iterator<Protein> pitr = idReader.getProteins(); //need to run to calculate redundnat peptides
            int redundantPeptideNum = idReader.getTotalPeptideNumber();
            double percent = 0.0;
            double eachSeg = 100.0 / redundantPeptideNum;
            int pepCount = 0;

            double medianRatio = -1;
            if (conf.isCalculateEnrich()) {
                IdentificationReader tempReader = BaseIdentificationReader.getIdentificationInst(isoReader);
                medianRatio = getMedianN15Ratio(tempReader, isoReader, ht);

                Element enrichmentEle = new Element("n15_enrich_median");
                enrichmentEle.setText(String.valueOf(medianRatio));
                rootEle.addContent(enrichmentEle);
            }

            ArrayList<Protein> aList = new ArrayList<Protein>();
            List<edu.scripps.pms.census.model.XYPoint> xyPoints;

            for (Iterator<Protein> itr = pitr; itr.hasNext();) {
                Protein protein = itr.next();


               // System.out.println("===\t" + protein.getLocus());

                proteinEle = new Element("protein");
                proteinEle.setAttribute("locus", protein.getLocus());
                proteinEle.setAttribute("seq_ct", protein.getSeqCount());
                proteinEle.setAttribute("spec_ct", protein.getSpectrumCount());
                proteinEle.setAttribute("seq_cov", protein.getSeqCoverage());
                proteinEle.setAttribute("length", protein.getLength());
                proteinEle.setAttribute("molwt", protein.getMolWt());
                proteinEle.setAttribute("pi", protein.getPI());
                proteinEle.setAttribute("val", protein.getValidation());

                String tmpStr = protein.getlSpectrumCount();
                if (null != tmpStr && !"N/A".equals(tmpStr)) {
                    proteinEle.setAttribute("lspec_ct", tmpStr);
                }

                tmpStr = protein.gethSpectrumCount();
                if (null != tmpStr && !"N/A".equals(tmpStr)) {
                    proteinEle.setAttribute("hspec_ct", tmpStr);
                }

                try {
                    proteinEle.setAttribute("desc", protein.getDescription());
                } catch (org.jdom.IllegalDataException ide) {
                    proteinEle.setAttribute("desc", StringUtil.removeIsoControlChar(protein.getDescription()));
                }

                int peptideSize = protein.getPeptideSize();
                int tmpPepCount=0;
                for (Iterator<Peptide> pepItr = protein.getPeptides(); pepItr.hasNext();) {
                    Peptide peptide = pepItr.next();
                    pepCount++;
                    tmpPepCount++;

                    //trim additional characters from peptide sequence at both ends
                    String pepSequence = peptide.getSequence();
                    //  System.out.println("---------------------------");

                    //element = new ElementComposition(peptide.getSequence().substring(2, peptide.getSequence().length()-2), isoReader.getIsotope());


                    int lightCount = org.apache.commons.lang3.StringUtils.countMatches(pepSequence, pulseLightResidueMass);
                    if(lightCount<=0)
                        lightCount = org.apache.commons.lang3.StringUtils.countMatches(pepSequence, pulseLightResidueMassThreeDigit);

                  if(!pulseLightResidueMass.contains("(")) lightCount=0;

                    // if(conf.getPulseLabelingCount()==0)
                    //    pulseCount = org.apache.commons.lang3.StringUtils.countMatches(pepSequence, conf.getPulseResidue()); //.getPulseLightMass());
                    int heavyCount = org.apache.commons.lang3.StringUtils.countMatches(pepSequence, pulseHeavyResidueMass);
                    if(heavyCount<=0)
                        heavyCount = org.apache.commons.lang3.StringUtils.countMatches(pepSequence, pulseHeavyResidueMassThreeDigit);//int totalCount = lightCount + heavyCount;

                    //if(pulseCount>1 && conf.) pulseCount=1;

                    String cleanedSeq = pepSequence.substring(2, peptide.getSequence().length() - 2);
                 //   cleanedSeq = cleanedSeq.replaceAll("\\(" + conf.getPulseLightMass() + "\\)", "");
                 //   cleanedSeq = cleanedSeq.replaceAll("\\(" + conf.getPulseHeavyMass() + "\\)", "");


                    if(cleanedSeq.contains(pulseHeavyResidueMass) || cleanedSeq.contains(pulseHeavyResidueMassThreeDigit))
                        labelingType = HEAVY_LABEL;
                    else if(cleanedSeq.contains(pulseLightResidueMass) || cleanedSeq.contains(pulseLightResidueMassThreeDigit))
                        labelingType = LIGHT_LABEL;
                    else
                        labelingType = NO_LABEL;

                    char[] ch = cleanedSeq.toCharArray();

                    boolean isLight = true;

			    element = new ElementComposition(ch, 0, ch.length, isoReader.getIsotope());

                try {
                    element.calculate();
                } catch (InvalidAAException invE) {
                    System.out.println("Not Quantifiable peptide : " + pepSequence);
                    percent += eachSeg;

                    if (null != progress) {
                        ChroProgressDialog.addMessageWithLine(progress, "Not Quantifiable peptide : " + pepSequence);
                        progress.setProgress((int) percent);
                    }


                    if( proteinEle.getChildren().size() <= 0 && (pepCount == peptideSize)) aList.clear();

				    continue;
			    }

			    if (!element.isQuantifiable()) {
                    System.out.print("\nError : ");
                    System.out.println(pepSequence + " is not quantifiable.");

                    percent += eachSeg;
                    System.out.print(pepCount);
                    System.out.print("/");
                    System.out.print(redundantPeptideNum);
                    System.out.print(" peptides, ");
                    System.out.print((int) percent);
                    System.out.print(" % is complete\r");

                    if (null != progress) {
                        progress.addMessage("\nError : ");
                        progress.addMessage(pepSequence);
                        progress.addMessage(" is not quantifiable.\n");
                        progress.setProgress((int) percent);
                    }


                    if( proteinEle.getChildren().size() <= 0 && (pepCount == peptideSize)) aList.clear();

                    continue;
			    }

			    //System.out.println( keyIndex, iFile, samIsoArr, refIsoArr, range) );
			    String fileName = peptide.getFileName();
			    fileName = cleanFileName(fileName);

			    switch (conf.getSpectrumFormat()) {
				case Configuration.MS_FILE_FORMAT:
				    fileName += ".ms1";
				    break;

				case Configuration.MZXML_FILE_FORMAT:
				    fileName += ".mzXML";
				    break;

				default:
				    break;
			    }

			    int scanNum = Integer.parseInt(peptide.getScanNum());
			    iFile = ht.get(filePath + fileName);

			    if (null == iFile) {
				iFile = ht.get(filePath + fileName.substring(1));

				if (null == iFile) {

				    int myIndex = fileName.indexOf("_my_");
				    if (myIndex <= 0) {
					System.out.println("Error : cannot find the file " + fileName);
					System.exit(0);
				    }

				    String unSplitMs1File = fileName.substring(0, myIndex);
				    iFile = ht.get(filePath + unSplitMs1File + ".ms1");
				    if (null == iFile) {
					System.out.println("Error : cannot find the file " + fileName);
					System.exit(0);
				    }
				}
			    }

			    int chargeState = Integer.parseInt(peptide.getChargeState());

			    IsotopeDist sampleDist = null;
			    //IsotopeDist refDist = null;
			   // IsotopeDist ref2Dist = null;

			    sampleDist = new IsotopeDist(element.getElementSampleArr(), element.getModShift(), true);
			    /*
	int[] tarrr = element.getElementSampleArr();
	for(int ii=0;ii<tarrr.length;ii++)
	System.out.println("=========" + tarrr[ii]);
	System.out.println("=========---------");
	tarrr = element.getElementRefArr();
	for(int ii=0;ii<tarrr.length;ii++)
	System.out.println("=========" + tarrr[ii]);
	System.out.println("=========---------");
	tarrr = element.getHeavyArr();
	for(int ii=0;ii<tarrr.length;ii++)
	System.out.println("=========" + tarrr[ii]);
	element.printComposition();
			     */

			    peptideEle = this.createXmlChroPeptideTitle(true, peptide); //true is for full scan

			    //  refDist = new IsotopeDist(element.getElementRefArr(), element.getModShift(), false);
			    peptideEle.setAttribute("enrichment", String.valueOf(conf.getEnrichment()));


			    //String fileName = peptide.getFileName().substring(0, peptide.getFileName().indexOf(".")+1) + "ms1";
			    StringBuffer rangeKey = new StringBuffer();
			    rangeKey.append(protein.getLocus());
			    rangeKey.append(fileName);
			    rangeKey.append(peptide.getSequence().substring(2, peptide.getSequence().length() - 2));

			    SpecRange range = rangeGen.getSpecRange(rangeKey.toString());
			    //sb.append("\tStartScan\tEndScan\tDTAPeakStart\tDTAPeakEnd\n");

			    if (null == range) {
				int tmpScanNum = Integer.parseInt(peptide.getScanNum());
				range = new SpecRange(tmpScanNum, tmpScanNum);
				peptideEle.setAttribute("start_scan", peptide.getScanNum());
				peptideEle.setAttribute("end_scan", peptide.getScanNum());
			    } else {
				peptideEle.setAttribute("start_scan", String.valueOf(range.getMin()));
				peptideEle.setAttribute("end_scan", String.valueOf(range.getMax()));
			    }

			    peptideEle.setAttribute("start_scan", String.valueOf(range.getMin()));
			    peptideEle.setAttribute("end_scan", String.valueOf(range.getMax()));
	//		    int scanNum = Integer.parseInt(peptide.getScanNum());

			    keys = iFile.getKeys();
			    keyIndex = Arrays.binarySearch(keys, scanNum);

			    if (keyIndex < 0) //Cannot find index
			    {
				keyIndex = -(++keyIndex); //Math.abs(++keyIndex);
			    }
			    if (keyIndex >= keys.length) {
				keyIndex--;
			    }

			    Element chro = new Element("chro");
	//harshil shah
			    Element iso = null;
			    Element theoMass = null;
			    Element theoIntensity = null;
//till here....

                    iso = new Element("iso");
                    theoMass = new Element("theo-mass");
                    theoIntensity = new Element("theo-int");

                    double[] samIsoArr = sampleDist.getHighMassList();
                    double[] refIsoArr = sampleDist.getHighMassList();
                    boolean noLabel=false;
                    for (int i = 0; i < samIsoArr.length; i++) {
                       // samIsoArr[i] += pulseLightMassDouble * pulseCount; //pulse labeling
                        if(labelingType == LIGHT_LABEL) {
                            samIsoArr[i] = (samIsoArr[i] + chargeState * PROTON_MASS) / chargeState;
                            //  refIsoArr[i] += pulseHeavyMassDouble * pulseCount;
                            refIsoArr[i] = (refIsoArr[i] + massDiff*lightCount + chargeState * PROTON_MASS) / chargeState;
                        } else if(labelingType == HEAVY_LABEL) {
                            samIsoArr[i] = (samIsoArr[i] - massDiff*heavyCount + chargeState * PROTON_MASS) / chargeState;
                            refIsoArr[i] = (refIsoArr[i] + chargeState * PROTON_MASS) / chargeState;
                        } else {
                            noLabel = true;
                            break;
                        }
                    }

                    if(noLabel) {

                        if( proteinEle.getChildren().size() <= 0 && (tmpPepCount == peptideSize)) {
                            aList.clear();
                            protein.setValid(false);
                            continue;
                        }
                        continue;
                    }

                    double[] samIsoIntArr = sampleDist.getRelabun();    //light intensities
                    double[] refIsoIntArr = sampleDist.getRelabun();       //heavy intensities

                    //double[] refIsoArr = null; //refDist.getHighMassList();
                    //double[] ref2IsoArr = null; //refDist.getHighMassList();
                    chro.setText(CalcUtil.calculateFullMS(keyIndex, iFile, samIsoArr, refIsoArr, range, null));

                    StringBuffer isoText = new StringBuffer();
                    for (double value : samIsoArr) {
                        isoText.append(value + ",");
                    }
                    isoText.append(":");
                    for (double value : refIsoArr) {
                        isoText.append(value + ",");
                    }
                    iso.setText(isoText.toString());

                    iso.setText(IsoCalcUtil.calculateFullMS(keyIndex, iFile, samIsoArr, refIsoArr, range));
                    theoMass.setText(IsoCalcUtil.massWrite(samIsoArr, refIsoArr));
                    theoIntensity.setText(IsoCalcUtil.intensityWrite(samIsoIntArr, refIsoIntArr, samIsoArr.length, refIsoArr.length));

                    /*
System.out.println( " " + peptide.getSequence());
System.out.println( " " + peptide.getChargeState());
System.out.println( " " + peptide.getScanNum());
System.exit(0);
                     */
                    peptideEle.addContent(chro);

//added Harshil Shah--------------------
                    if (iso != null) {
                        peptideEle.addContent(iso);
                    }
                    if (theoMass != null) {
                        peptideEle.addContent(theoMass);
                    }
                    if (theoIntensity != null) {
                        peptideEle.addContent(theoIntensity);
                    }
//till here-- harshil SHah-----------------
                    proteinEle.addContent(peptideEle);

                    percent += eachSeg;

                    System.out.print(pepCount);
                    System.out.print("/");
                    System.out.print(redundantPeptideNum);
                    System.out.print(" peptides, ");
                    System.out.print((int) percent);
                    System.out.print(" % is complete\r");

                    if (conf.isPrintLog()) {
                        p.println((int) percent + "\t" + pepCount + "\t" + redundantPeptideNum);
                    }

                    if (null != progress) {
                        progress.setProgress((int) percent);
                    }


                    /*
                    IndexedFile iFile2 = ht.get(filePath+ peptide.getFileName() + "." + "ms2");
                    List<ReportIon> reportIonList = conf.getReportIonList();
                    double startmass = reportIonList.get(0).getMass();
                    double endMass = reportIonList.get(reportIonList.size()-1).getMass();
                    int [] keys2 = iFile2.getKeys();

                    //  long start = System.currentTimeMillis();
                    int ms3Scan;

                    String ms3FilePattern="orig";
                    if (conf.isMs3ScanRandom()) {

                        String ms3FileName = filePath + peptide.getFileName() + ".ms3";

                        IndexedFile ms3Ifile = ht.get(ms3FileName);




                        if (null == ms3Ifile) { // for heavy file
                            ms3Ifile = ht.get(filePath+ peptide.getFileName().substring(1) + "." + "ms3");

                            if (null == ms3Ifile) { // for replaced ext

                                ms3Ifile = ht.get(ms3FileName.replace(".ms3", "_ms3.ms2"));
                                ms3FilePattern="_ms3.ms2";

                                if (null == ms3Ifile) {

                                    if (peptide.getFileName().startsWith("H")
                                            || peptide.getFileName().startsWith("M")) {
                                        String heavyMs3File =filePath + peptide.getFileName().substring(1)
                                                + "_ms3.ms2";
                                        ms3Ifile = ht.get(heavyMs3File);
                                    }
                                    if (null == ms3Ifile) {
                                        System.out.println("Error : cannot find the file " + peptide.getFileName()
                                                + "." + "ms3");
                                    }
                                }
                            }
                        }


                        int keyIndex2 = Arrays.binarySearch(keys2, Integer.parseInt(peptide.getScanNum()));

                        ms3Scan = TMTUtil.getCorrespondingMs3Scan(ms3Ifile, keys2[keyIndex2]);

                        if(ms3FilePattern.equals("orig")) {
                            xyPoints = SpectrumUtil.getSpectrum(filePath, peptide.getFileName() + ".ms3", ms3Scan,
                                    startmass-1, endMass+1 );
                        }
                        else {
                            xyPoints = SpectrumUtil.getSpectrum(filePath, peptide.getFileName() + "_ms3.ms2", ms3Scan,
                                    startmass-1, endMass+1 );
                        }
                    }
                    else{
                        xyPoints = SpectrumUtil.getSpectrum(filePath, peptide.getFileName() + ".ms2",
                                Integer.parseInt(peptide.getScanNum()),startmass -1 , endMass + 1);
                    }
                    String pepKey = peptide.getSequence() + peptide.getChargeState() + peptide.getScanNum();
                    FileUtil.writeJSON(gson.toJson(xyPoints), xyvaluePath +  pepKey + "_xy.JSON");

                    */



                }


                if (proteinEle.getChildren().size() > 0) {

                    Element redunEle = new Element("redundant");

                    int rPepCount = 0;

                    for (Iterator<Protein> rItr = aList.iterator(); rItr.hasNext();) {
                        rPepCount++;
                        Protein rp = rItr.next();

                        Element rpEle = new Element("protein");
                        rpEle.setAttribute("locus", rp.getLocus());
                        rpEle.setAttribute("seq_ct", rp.getSeqCount());
                        rpEle.setAttribute("spec_ct", rp.getSpectrumCount());
                        rpEle.setAttribute("seq_cov", rp.getSeqCoverage());
                        rpEle.setAttribute("length", rp.getLength());
                        rpEle.setAttribute("molwt", rp.getMolWt());
                        rpEle.setAttribute("pi", rp.getPI());
                        rpEle.setAttribute("val", rp.getValidation());

                        try {
                            rpEle.setAttribute("desc", rp.getDescription());
                        } catch (org.jdom.IllegalDataException ide) {
                            rpEle.setAttribute("desc", StringUtil.removeIsoControlChar(rp.getDescription()));
                        }

                        redunEle.addContent(rpEle);
                    }

                    if (rPepCount > 0) {
                        proteinEle.addContent(redunEle);
                    }

                    rootEle.addContent(proteinEle);
                    aList.clear();
                } else {
                    if(protein.isValid())
                        aList.add(protein);
                    //	    set.add(protein.getLocus());
                }
            }

            Document doc = new Document(rootEle);
            OutputStream os = new FileOutputStream(filePath + "census_chro.xml");
            XMLOutputter outputter = new XMLOutputter();
            outputter.setFormat(Format.getPrettyFormat());
            outputter.output(doc, os);
            os.close();

            System.out.println("\n100% complete");
        } catch (IOException e) {
            System.out.println("IO Error while generating msms chro file : " + e);
            e.printStackTrace();
            throw new IOException(e.toString());
        } catch (Exception e) {
            System.out.println("Error while generating msms chro file : " + e);
            e.printStackTrace();
            throw new Exception(e.toString());
        } finally {
            if (null != p) {
                p.close();
            }

            if (null != out) {
                out.close();
            }

            //Close all random files
            for (Enumeration e = ht.keys(); e.hasMoreElements();) {
                iFile = ht.get(e.nextElement());

                if (null != iFile) {
                    iFile.close();
                }
            }

        }

        System.out.println((System.currentTimeMillis() - startTime) * 0.001 + " seconds taken");
    }

    //labeling labeled label
    // MultipleMs1Labeling  dimethyl 3plex
    public void createFullscanXmlChroTriple() throws IOException, CensusGeneralException, Exception {

        int[] keys;
        this.filePath = conf.getFilePath();

        Hashtable<String, IndexedFile> ht = null;

        switch (conf.getSpectrumFormat()) {
            case Configuration.MS_FILE_FORMAT:
                ht = createIndexedFiles(filePath, CensusConstants.MS1_FILE);
                break;

            case Configuration.MZXML_FILE_FORMAT:
                ht = createIndexedFiles(filePath, CensusConstants.MZXML);
                break;

            default:
                break;
        }

        if (!filePath.endsWith(File.separator)) {
            filePath += File.separator;
        }

        IndexedFile iFile = null;
        BufferedOutputStream out = null;
        PrintStream p = null;

        long startTime = System.currentTimeMillis();

        try {
            /**
             * ****************************************************************
             * Read DTASelect.txt file to find spectrum range for each peptide
		     *****************************************************************
             */
            //TIntLongHashMap index;
            System.out.print("Parsing " + conf.getIdFileName() + "...");

            if (null != progress) {
                progress.addMessage("\nParsing Identified Peptides...");
            }

            if (null != progress) {
                progress.addMessage("\ndone.");
            }

            IsotopeReader isoReader = null;
            isoReader = new IsotopeReader(this.confRootEle);

            IdentificationReader idReader = BaseIdentificationReader.getIdentificationInst(isoReader);
            System.out.println("done.");

            // SpecRangeGenerator rangeGen = null;
            //File dtaFile = new File(filePath + "DTASelect.txt");
            if (conf.isPrintLog()) {
                p = new PrintStream(new FileOutputStream(filePath + "progress.log"));
            }

            //System.out.println(filePath + " " + dtaFile + " " + dtaFile.exists());
            Element rootEle = this.createXmlChroHeader(1, conf.getQuantType());
            ElementComposition element;
            // ElementComposition totalElement;

            int keyIndex = -1;
            //int start;
            //int last;

            Element proteinEle = null;
            Element peptideEle = null;

            Iterator<Protein> pitr = idReader.getProteins(); //need to run to calculate redundnat peptides
            int redundantPeptideNum = idReader.getTotalPeptideNumber();
            double percent = 0.0;
            double eachSeg = 100.0 / redundantPeptideNum;
            int pepCount = 0;

            double medianRatio = -1;

            ArrayList<Protein> aList = new ArrayList<Protein>();

            for (Iterator<Protein> itr = pitr; itr.hasNext();) {
                Protein protein = itr.next();

                proteinEle = new Element("protein");
                proteinEle.setAttribute("locus", protein.getLocus());
                proteinEle.setAttribute("seq_ct", protein.getSeqCount());
                proteinEle.setAttribute("spec_ct", protein.getSpectrumCount());
                proteinEle.setAttribute("seq_cov", protein.getSeqCoverage());
                proteinEle.setAttribute("length", protein.getLength());
                proteinEle.setAttribute("molwt", protein.getMolWt());
                proteinEle.setAttribute("pi", protein.getPI());
                proteinEle.setAttribute("val", protein.getValidation());

                String tmpStr = protein.getlSpectrumCount();
                if (null != tmpStr && !"N/A".equals(tmpStr)) {
                    proteinEle.setAttribute("lspec_ct", tmpStr);
                }

                tmpStr = protein.gethSpectrumCount();
                if (null != tmpStr && !"N/A".equals(tmpStr)) {
                    proteinEle.setAttribute("hspec_ct", tmpStr);
                }

                try {
                    proteinEle.setAttribute("desc", protein.getDescription());
                } catch (org.jdom.IllegalDataException ide) {
                    proteinEle.setAttribute("desc", StringUtil.removeIsoControlChar(protein.getDescription()));
                }

                for (Iterator<Peptide> pepItr = protein.getPeptides(); pepItr.hasNext();) {
                    Peptide peptide = pepItr.next();
                    pepCount++;

                    //trim additional characters from peptide sequence at both ends
                    String pepSequence = peptide.getSequence();
                    /*System.out.println("===" + protein.getLocus());
				    System.out.println("\n\n");
				    System.out.println(pepSequence);
				    System.out.println("\n\n");
				    System.out.println(peptide.getFileName());
				    System.out.println("\n\n");*/

                    char[] ch = pepSequence.substring(2, peptide.getSequence().length() - 2).toCharArray();
                    //element = new ElementComposition(peptide.getSequence().substring(2, peptide.getSequence().length()-2), isoReader.getIsotope());

                    element = new ElementComposition(ch, 0, ch.length, isoReader.getIsotope());

                    try {
                        element.calculate();
                    } catch (InvalidAAException invE) {
                        System.out.println("Not Quantifiable peptide : " + pepSequence);
                        percent += eachSeg;

                        if (null != progress) {
                            ChroProgressDialog.addMessageWithLine(progress, "Not Quantifiable peptide : " + pepSequence);
                            progress.setProgress((int) percent);
                        }
                        continue;
                    }

                    /*
				       IsotopeTable<String, int[]> iit2 = isoReader.getIsotope();
				       int[] iiii= iit2.get("sampleK");
				       for(int ii:iiii)
				       System.out.println("before================" + ii);
				       System.out.println("================");
				       iiii= iit2.get("refK");
				       for(int ii:iiii)
				       System.out.println("before================" + ii);


				       for(int ii:element.getElementSampleArr())
				       System.out.println("ss after------------------================" + ii);
				       for(int ii:element.getElementRefArr())
				       System.out.println("rr after------------------================" + ii);
                     */
                    if (!element.isQuantifiable()) {
                        System.out.print("\nError : ");
                        System.out.println(pepSequence + " is not quantifiable.");

                        percent += eachSeg;
                        System.out.print(pepCount);
                        System.out.print("/");
                        System.out.print(redundantPeptideNum);
                        System.out.print(" peptides, ");
                        System.out.print((int) percent);
                        System.out.print(" % is complete\r");

                        if (null != progress) {
                            progress.addMessage("\nError : ");
                            progress.addMessage(pepSequence);
                            progress.addMessage(" is not quantifiable.\n");
                            progress.setProgress((int) percent);
                        }

                        continue;
                    }

                    //System.out.println( keyIndex, iFile, samIsoArr, refIsoArr, range) );
                    String fileName = peptide.getFileName();
                    fileName = cleanFileName(fileName);

                    switch (conf.getSpectrumFormat()) {
                        case Configuration.MS_FILE_FORMAT:
                            fileName += ".ms1";
                            break;

                        case Configuration.MZXML_FILE_FORMAT:
                            fileName += ".mzXML";
                            break;

                        default:
                            break;
                    }

                    int scanNum = Integer.parseInt(peptide.getScanNum());
                    iFile = ht.get(filePath + fileName);

                    if (null == iFile) {
                        iFile = ht.get(filePath + fileName.substring(1));

                        if (null == iFile) {

                            int myIndex = fileName.indexOf("_my_");
                            if (myIndex <= 0) {
                                System.out.println("Error11 : cannot find the file " + fileName);
                                System.exit(0);
                            }

                            String unSplitMs1File = fileName.substring(0, myIndex);
                            iFile = ht.get(filePath + unSplitMs1File + ".ms1");
                            if (null == iFile) {
                                System.out.println("Error22 : cannot find the file " + fileName);
                                System.exit(0);
                            }
                        }
                    }

                    int chargeState = Integer.parseInt(peptide.getChargeState());

                    IsotopeDist sampleDist = null;
                    IsotopeDist refDist = null;
                    IsotopeDist ref2Dist = null;

                    sampleDist = new IsotopeDist(element.getElementSampleArr(), element.getModShift(), true);
                    /*
int[] tarrr = element.getElementSampleArr();
for(int ii=0;ii<tarrr.length;ii++)
System.out.println("=========" + tarrr[ii]);
System.out.println("=========---------");
tarrr = element.getElementRefArr();
for(int ii=0;ii<tarrr.length;ii++)
System.out.println("=========" + tarrr[ii]);
System.out.println("=========---------");
tarrr = element.getHeavyArr();
for(int ii=0;ii<tarrr.length;ii++)
System.out.println("=========" + tarrr[ii]);
element.printComposition();
                     */

//		    System.out.println("===========" + conf.getQuantType());
/*
		    System.out.println("===========" + element.getElementSampleArr());

		    //print element composition
		    System.out.println("sample===========");
		    int[] tempArr = element.getElementSampleArr();
		    for(int t : tempArr)
			System.out.println(t);
		    System.out.println("ref===========");
		    tempArr = element.getElementRefArr();
		    for(int t : tempArr)
                     */
                    peptideEle = this.createXmlChroPeptideTitle(true, peptide); //true is for full scan

                    /*
			    System.out.println("rarr===========");
			    System.out.println("rarr===========");
			    System.out.println("rarr===========");
			    System.out.println("rarr===========");
			int[] eleArr = element.getElementSampleArr();
			for(int ii:eleArr)
			    System.out.println("rarr===========" + ii);

			int[] eleArr = element.getElementSampleArr();
			for(int ii:eleArr)
			    System.out.println("sarr===========" + ii);
			eleArr = element.getElementRefArr();
			    System.out.println("===========");
                     */
                    try {

                        refDist = new IsotopeDist(element.getElementRefArr(), element.getModShift(), false);
                        ref2Dist = new IsotopeDist(element.getHeavyArr(), element.getModShift(), false);

                    } catch (Exception invE) {
                        System.out.println("Isotope Distribution is out of boundary: " + pepSequence);
                        invE.printStackTrace();
                        percent += eachSeg;

                        if (null != progress) {
                            ChroProgressDialog.addMessageWithLine(progress, "Not Quantifiable peptide : " + pepSequence);
                            progress.setProgress((int) percent);
                        }
                        continue;
                    }

                    /*
			double[] dmlist = sampleDist.getHighMassList();
			for(double dd:dmlist) {
			    if(dd<=0)
				continue;
			    System.out.println("sam mass :\t" + dd);
			}

			    System.out.println("===========");
                            System.out.println("==========="+ refDist);
			dmlist = refDist.getHighMassList();
			for(double dd:dmlist) {
			    if(dd<=0)
				continue;
			    System.out.println("ref mass :\t" + dd);
			}*/
                    conf.setCalcSamAvgMass(sampleDist.getAvgMass());
                    conf.setCalcRefAvgMass(refDist.getAvgMass());
                    conf.setCalcRef2AvgMass(ref2Dist.getAvgMass());
                    peptideEle.setAttribute("heavy2StartMass", String.valueOf(ref2Dist.getStartMass()));
                    peptideEle.setAttribute("heavy2AvgMass", String.valueOf(conf.getCalcRef2AvgMass()));

                    //String fileName = peptide.getFileName().substring(0, peptide.getFileName().indexOf(".")+1) + "ms1";
                    StringBuffer rangeKey = new StringBuffer();
                    rangeKey.append(protein.getLocus());
                    rangeKey.append(fileName);
                    rangeKey.append(peptide.getSequence().substring(2, peptide.getSequence().length() - 2));

                    //  SpecRange range = rangeGen.getSpecRange( rangeKey.toString() );
                    //sb.append("\tStartScan\tEndScan\tDTAPeakStart\tDTAPeakEnd\n");

                    peptideEle.setAttribute("start_scan", peptide.getScanNum());
                    peptideEle.setAttribute("end_scan", peptide.getScanNum());

//		    int scanNum = Integer.parseInt(peptide.getScanNum());
                    keys = iFile.getKeys();
                    keyIndex = Arrays.binarySearch(keys, scanNum);

                    if (keyIndex < 0) //Cannot find index
                    {
                        keyIndex = -(++keyIndex); //Math.abs(++keyIndex);
                    }
                    if (keyIndex >= keys.length) {
                        keyIndex--;
                    }

                    Element chro = new Element("chro");
//harshil shah
                    Element iso = null;
                    Element theoMass = null;
                    Element theoIntensity = null;

                    iso = new Element("iso");
                    theoMass = new Element("theo-mass");
                    theoIntensity = new Element("theo-int");

                    double[] samIsoArr = sampleDist.getHighMassList();
                    double[] samIsoIntArr = sampleDist.getRelabun();    //light intensities
                    double[] refIsoIntArr = refDist.getRelabun();       //heavy intensities
                    double[] ref2IsoIntArr = ref2Dist.getRelabun();


                    double[] refIsoArr = null; //refDist.getHighMassList();
                    double[] ref2IsoArr = null; //refDist.getHighMassList();

                    for (int i = 0; i < samIsoArr.length; i++) {
                        samIsoArr[i] = (samIsoArr[i] + chargeState * PROTON_MASS) / chargeState;
                        //   System.out.println("sam m/z\t =======================" + samIsoArr[i]);

                    }

                    if (conf.isUseMassDiff()) {
                        refIsoArr = IsotopeDist.getHeavySilacDist(sampleDist.getHighMassList(), element.getElementRefArr());
                        for (int i = 0; i < refIsoArr.length; i++) {
                            //					    					    System.out.println("r22 before=======================" + refIsoArr[i]);
                            refIsoArr[i] = (refIsoArr[i] + chargeState * PROTON_MASS) / chargeState;
                            //					    System.out.println("r33 before=======================" + refIsoArr[i]);
                        }
                    } else {

                        refIsoArr = refDist.getHighMassList();
                        for (int i = 0; i < refIsoArr.length; i++) {
                            //					    					    System.out.println("r22 before=======================" + refIsoArr[i]);
                            refIsoArr[i] = (refIsoArr[i] + chargeState * PROTON_MASS) / chargeState;
                            //System.out.println("ref1 m/z\t =======================" + refIsoArr[i]);
                            //					    				System.out.println("ref :\t " + refIsoArr[i] + "\t" + refDist.getRelabun()[i]);
                        }

                        ref2IsoArr = ref2Dist.getHighMassList();
                        for (int i = 0; i < ref2IsoArr.length; i++) {

                            ref2IsoArr[i] = (ref2IsoArr[i] + chargeState * PROTON_MASS) / chargeState;
                            //System.out.println("ref2 m/z\t=======================" + ref2IsoArr[i]);
                        }
                    }

                 //   double tolerance = 0.2;
                    int stopIndex = 0;

                    for (int i = 0; i < samIsoArr.length; i++) {
                        //if (samIsoArr[i] > (refIsoArr[0] - 4.2)) {
                        if (samIsoArr[i] > (refIsoArr[0] - conf.getMassTolerance())) {
                            stopIndex = i - 1;
                            break;
                        }
                        stopIndex++;
                    }

                    if (stopIndex < 0) {
                        stopIndex = 0;
                    }
                    double[] tmpArr1 = new double[stopIndex + 1];
                    double[] tmpArr2 = new double[stopIndex + 1];
                    double[] tmpArr3 = new double[stopIndex + 1];
                    if(tmpArr1.length<samIsoArr.length) {
                      for (int i = 0; i < tmpArr1.length; i++) {
                        tmpArr1[i] = samIsoArr[i];
                        tmpArr2[i] = refIsoArr[i];
                        tmpArr3[i] = ref2IsoArr[i];
                      }
                      samIsoArr = tmpArr1;
                      refIsoArr = tmpArr2;
                      ref2IsoArr = tmpArr3;
                    }


                    System.out.println("<<>>>> "+ peptide.getScanNum());
			       for(int i=0;i<samIsoArr.length;i++)
			       {
			       System.out.println("s before=======================" + samIsoArr[i]);
			       }
			       for(int i=0;i<refIsoArr.length;i++)
			       {
			       System.out.println("rr before=======================" + refIsoArr[i]);
			       }

                    //boolean isProline = true;
                    //double prolineMassDiff = 6.0138;
                    //                        System.out.println(pepSequence);
                    //       System.out.println("=======================================");
                    String chroTxt = CalcUtil.calculateFullMS3Plex(keyIndex, iFile, samIsoArr, refIsoArr, ref2IsoArr);
                    if (null == chroTxt) {
                        continue;
                    }

                    chro.setText(chroTxt);
                    //chro.setText( CalcUtil.calculateFullMS3Plex( keyIndex, iFile, samIsoArr, refIsoArr, ref2IsoArr, range) );

                    StringBuffer isoText = new StringBuffer();
                    for (double value : samIsoArr) {
                        isoText.append(value + ",");
                    }
                    isoText.append(":");
                    for (double value : refIsoArr) {
                        isoText.append(value + ",");
                    }
                    iso.setText(isoText.toString());

                    iso.setText(IsoCalcUtil.calculateFullMS(keyIndex, iFile, samIsoArr, refIsoArr, null));
                    //iso.setText(IsoCalcUtil.calculateFullMS( keyIndex, iFile, samIsoArr, refIsoArr, range));
                    theoMass.setText(IsoCalcUtil.massWrite(samIsoArr, refIsoArr, ref2IsoArr));
                    theoIntensity.setText(IsoCalcUtil.intensityWrite(samIsoIntArr, refIsoIntArr, ref2IsoIntArr));
                    //theoIntensity.setText(IsoCalcUtil.intensityWrite(samIsoIntArr, refIsoIntArr, refsamIsoArr.length, refIsoArr.length));



// Till here-- harshil Shah

                    /*
System.out.println( " " + peptide.getSequence());
System.out.println( " " + peptide.getChargeState());
System.out.println( " " + peptide.getScanNum());
System.exit(0);
                     */
                    peptideEle.addContent(chro);

//added Harshil Shah--------------------
                    if (iso != null) {
                        peptideEle.addContent(iso);
                    }
                    if (theoMass != null) {
                        peptideEle.addContent(theoMass);
                    }
                    if (theoIntensity != null) {
                        peptideEle.addContent(theoIntensity);
                    }
//till here-- harshil SHah-----------------
                    proteinEle.addContent(peptideEle);

                    percent += eachSeg;

                    System.out.print(pepCount);
                    System.out.print("/");
                    System.out.print(redundantPeptideNum);
                    System.out.print(" peptides, ");
                    System.out.print((int) percent);
                    System.out.print(" % is complete\r");

                    if (conf.isPrintLog()) {
                        p.println((int) percent + "\t" + pepCount + "\t" + redundantPeptideNum);
                    }

                    if (null != progress) {
                        progress.setProgress((int) percent);
                    }
                }

                if (proteinEle.getChildren().size() > 0) {

                    Element redunEle = new Element("redundant");

                    int rPepCount = 0;

                    for (Iterator<Protein> rItr = aList.iterator(); rItr.hasNext();) {
                        rPepCount++;
                        Protein rp = rItr.next();

                        Element rpEle = new Element("protein");
                        rpEle.setAttribute("locus", rp.getLocus());
                        rpEle.setAttribute("seq_ct", rp.getSeqCount());
                        rpEle.setAttribute("spec_ct", rp.getSpectrumCount());
                        rpEle.setAttribute("seq_cov", rp.getSeqCoverage());
                        rpEle.setAttribute("length", rp.getLength());
                        rpEle.setAttribute("molwt", rp.getMolWt());
                        rpEle.setAttribute("pi", rp.getPI());
                        rpEle.setAttribute("val", rp.getValidation());

                        try {
                            rpEle.setAttribute("desc", rp.getDescription());
                        } catch (org.jdom.IllegalDataException ide) {
                            rpEle.setAttribute("desc", StringUtil.removeIsoControlChar(rp.getDescription()));
                        }

                        redunEle.addContent(rpEle);
                    }

                    if (rPepCount > 0) {
                        proteinEle.addContent(redunEle);
                    }

                    rootEle.addContent(proteinEle);
                    aList.clear();
                } else {
                    aList.add(protein);
                    //	    set.add(protein.getLocus());
                }
            }

            Document doc = new Document(rootEle);
            OutputStream os = new FileOutputStream(filePath + "census_chro.xml");
            XMLOutputter outputter = new XMLOutputter();
            outputter.setFormat(Format.getPrettyFormat());
            outputter.output(doc, os);
            os.close();

            System.out.println("\n100% complete");
        } catch (IOException e) {
            System.out.println("IO Error while generating msms chro file : " + e);
            e.printStackTrace();
            throw new IOException(e.toString());
        } catch (Exception e) {
            System.out.println("Error while generating msms chro file : " + e);
            e.printStackTrace();
            throw new Exception(e.toString());
        } finally {
            if (null != p) {
                p.close();
            }

            if (null != out) {
                out.close();
            }

            //Close all random files
            for (Enumeration e = ht.keys(); e.hasMoreElements();) {
                iFile = ht.get(e.nextElement());

                if (null != iFile) {
                    iFile.close();
                }
            }

        }

        System.out.println((System.currentTimeMillis() - startTime) * 0.001 + " seconds taken");
    }

    //labeling labeled label
    public void createFullscanXmlChro() throws IOException, CensusGeneralException, Exception {

        int[] keys;
        this.filePath = conf.getFilePath();

        Hashtable<String, IndexedFile> ht = null;

        switch (conf.getSpectrumFormat()) {
            case Configuration.MS_FILE_FORMAT:
                ht = createIndexedFiles(filePath, CensusConstants.MS1_FILE);
                break;

            case Configuration.MZXML_FILE_FORMAT:
                ht = createIndexedFiles(filePath, CensusConstants.MZXML);
                break;

            default:
                break;
        }

        if (!filePath.endsWith(File.separator)) {
            filePath += File.separator;
        }

        IndexedFile iFile = null;
        BufferedOutputStream out = null;
        PrintStream p = null;

        long startTime = System.currentTimeMillis();


        try {
            /**
             * ****************************************************************
             * Read DTASelect.txt file to find spectrum range for each peptide
		     *****************************************************************
             */
            //TIntLongHashMap index;
            Gson gson = new Gson();
            System.out.print("Parsing " + conf.getIdFileName() + "...");

            if (null != progress) {
                progress.addMessage("\nParsing Identified Peptides...");
            }

            if (null != progress) {
                progress.addMessage("\ndone.");
            }

            IsotopeReader isoReader = null;

            if (null != isotopeFile) {
                isoReader = new IsotopeReader(isotopeFile);
            } else {
                isoReader = new IsotopeReader(this.confRootEle);
            }

            File jsonDir = new File(filePath+MASTER_JSON);
            if(!jsonDir.exists())
            {
                jsonDir.mkdir();
            }


            IdentificationReader idReader = BaseIdentificationReader.getIdentificationInst(isoReader);
            System.out.println("done.");

            SpecRangeGenerator rangeGen = null;

            File dtaFile = new File(filePath + "DTASelect.txt");

            if (conf.isPrintLog()) {
                p = new PrintStream(new FileOutputStream(filePath + "progress.log"));
            }

            //System.out.println(filePath + " " + dtaFile + " " + dtaFile.exists());
            //if (dtaFile.exists()) {
             if (false) {  //skip DTASelect.txt parsing
                if (conf.getIdFileName().endsWith("txt")) {
                    System.out.print("Parsing DTASelect.txt...");
                    if (null != progress) {
                        progress.addMessage("\nParsing DTASelect.txt...");
                    }
                }

                rangeGen = SpecRangeGenerator.getSpecRangeGenerator(idReader);
                //            SpecRangeGenerator rangeGen = new SpecRangeGenerator(this.filePath + this.dtaSelectFile, idReader.isVersion2(), idReader.getConfidence());
                //SpecRangeGenerator rangeGen = new SpecRangeGenerator(this.filePath + this.dtaSelectFile, idReader);
                //SpecRangeGenerator rangeGen = null;

                if (conf.getIdFileName().endsWith("txt")) {
                    System.out.println("done.");
                    if (null != progress) {
                        progress.addMessage("\ndone.");
                    }
                }
            } else {
                rangeGen = new SpecRangeGenerator();

            }
            Element rootEle = this.createXmlChroHeader(1, conf.getQuantType());
            ElementComposition element;
            ElementComposition totalElement;

            int keyIndex = -1;
            //int start;
            //int last;

            double sampleStartMass;
            double sampleEndMass;
            double refStartMass;
            double refEndMass;

            Element proteinEle = null;
            Element peptideEle = null;

            Iterator<Protein> pitr = idReader.getProteins(); //need to run to calculate redundnat peptides
            int redundantPeptideNum = idReader.getTotalPeptideNumber();
            double percent = 0.0;
            double eachSeg = 100.0 / redundantPeptideNum;
            int pepCount = 0;

            double medianRatio = -1;
            if (conf.isCalculateEnrich()) {
                IdentificationReader tempReader = BaseIdentificationReader.getIdentificationInst(isoReader);
                medianRatio = getMedianN15Ratio(tempReader, isoReader, ht);

                Element enrichmentEle = new Element("n15_enrich_median");
                enrichmentEle.setText(String.valueOf(medianRatio));
                rootEle.addContent(enrichmentEle);
            }

            ArrayList<Protein> aList = new ArrayList<Protein>();

            for (Iterator<Protein> itr = pitr; itr.hasNext();) {
                Protein protein = itr.next();

                proteinEle = new Element("protein");
                proteinEle.setAttribute("locus", protein.getLocus());
                proteinEle.setAttribute("seq_ct", protein.getSeqCount());
                proteinEle.setAttribute("spec_ct", protein.getSpectrumCount());
                proteinEle.setAttribute("seq_cov", protein.getSeqCoverage());
                proteinEle.setAttribute("length", protein.getLength());
                proteinEle.setAttribute("molwt", protein.getMolWt());
                proteinEle.setAttribute("pi", protein.getPI());
                proteinEle.setAttribute("val", protein.getValidation());

                String tmpStr = protein.getlSpectrumCount();
                if (null != tmpStr && !"N/A".equals(tmpStr)) {
                    proteinEle.setAttribute("lspec_ct", tmpStr);
                }

                tmpStr = protein.gethSpectrumCount();
                if (null != tmpStr && !"N/A".equals(tmpStr)) {
                    proteinEle.setAttribute("hspec_ct", tmpStr);
                }

                try {
                    proteinEle.setAttribute("desc", protein.getDescription());
                } catch (org.jdom.IllegalDataException ide) {
                    proteinEle.setAttribute("desc", StringUtil.removeIsoControlChar(protein.getDescription()));
                }

                for (Iterator<Peptide> pepItr = protein.getPeptides(); pepItr.hasNext();) {
                    Peptide peptide = pepItr.next();
                    pepCount++;

                    //trim additional characters from peptide sequence at both ends
                    String pepSequence = peptide.getSequence();
                    /*
				    System.out.println("===" + protein.getLocus());
				    System.out.println("\n\n");
				    System.out.println(pepSequence);
				    System.out.println("\n\n");
				    System.out.println(peptide.getFileName());
				    System.out.println("\n\n");
                    */

                    char[] ch = pepSequence.substring(2, peptide.getSequence().length() - 2).toCharArray();
                    //element = new ElementComposition(peptide.getSequence().substring(2, peptide.getSequence().length()-2), isoReader.getIsotope());

                    element = new ElementComposition(ch, 0, ch.length, isoReader.getIsotope());

                    try {
                        element.calculate();
                    } catch (InvalidAAException invE) {
                        System.out.println("Not Quantifiable peptide : " + pepSequence);
                        percent += eachSeg;

                        if (null != progress) {
                            ChroProgressDialog.addMessageWithLine(progress, "Not Quantifiable peptide : " + pepSequence);
                            progress.setProgress((int) percent);
                        }
                        continue;
                    }

                    /*
				       IsotopeTable<String, int[]> iit2 = isoReader.getIsotope();
				       int[] iiii= iit2.get("sampleK");
				       for(int ii:iiii)
				       System.out.println("before================" + ii);
				       System.out.println("================");
				       iiii= iit2.get("refK");
				       for(int ii:iiii)
				       System.out.println("before================" + ii);


				       for(int ii:element.getElementSampleArr())
				       System.out.println("ss after------------------================" + ii);
				       for(int ii:element.getElementRefArr())
				       System.out.println("rr after------------------================" + ii);
                     */
                    if (!element.isQuantifiable()) {
                        System.out.print("\nError : ");
                        System.out.println(pepSequence + " is not quantifiable.");

                        percent += eachSeg;
                        System.out.print(pepCount);
                        System.out.print("/");
                        System.out.print(redundantPeptideNum);
                        System.out.print(" peptides, ");
                        System.out.print((int) percent);
                        System.out.print(" % is complete\r");

                        if (null != progress) {
                            progress.addMessage("\nError : ");
                            progress.addMessage(pepSequence);
                            progress.addMessage(" is not quantifiable.\n");
                            progress.setProgress((int) percent);
                        }

                        continue;
                    }

                    //System.out.println( keyIndex, iFile, samIsoArr, refIsoArr, range) );
                    String fileName = peptide.getFileName();
                    fileName = cleanFileName(fileName);

                    switch (conf.getSpectrumFormat()) {
                        case Configuration.MS_FILE_FORMAT:
                            fileName += ".ms1";
                            break;

                        case Configuration.MZXML_FILE_FORMAT:
                            fileName += ".mzXML";
                            break;

                        default:
                            break;
                    }

                    int scanNum = Integer.parseInt(peptide.getScanNum());
                    iFile = ht.get(filePath + fileName);
                    String msFile = fileName;
                    if (null == iFile) {
                        iFile = ht.get(filePath + fileName.substring(1));
                        msFile = fileName.substring(1);
                        if (null == iFile) {

                            int myIndex = fileName.indexOf("_my_");
                            if (myIndex <= 0) {
                                System.out.println("Error1 : cannot find the file " + fileName);

                                System.exit(0);
                            }

                            String unSplitMs1File = fileName.substring(0, myIndex);

                            msFile =unSplitMs1File + ".ms1";
                            iFile = ht.get(filePath + unSplitMs1File + ".ms1");
                            if (null == iFile) {
                                System.out.println("Error2 : cannot find the file " + fileName);
                                System.exit(0);
                            }
                        }
                    }
/*
                    TODELETE rereads ms1 files so is inefficient
                    List<XYPoint> points = SpectrumUtil.getSpectrumMS1Arr(filePath,msFile,peptide.getScanNumber(),ht);
                    DisplayData.DisplayChroData displayChroData = new DisplayData.DisplayChroData();
                    displayChroData.setXyPoints(points);
                  //  displayChroData.setMassAccuracy(mass_tol);

                    FileUtil.writeJSON(gson.toJson(displayChroData), filePath + "JSON_OBJ/" + peptide.getSequence() + peptide.getFileName() + peptide.getScanNum() + peptide.getChargeState() + ".JSON");

*/
                    int chargeState = Integer.parseInt(peptide.getChargeState());

                    IsotopeDist sampleDist = null;
                    IsotopeDist refDist = null;
                    IsotopeDist ref2Dist = null;

                    if ("15N".equals(conf.getQuantType())) {
                        sampleDist = new IsotopeDist15N(element.getElementSampleArr(), element.getModShift(), true);
                    } else {
                        sampleDist = new IsotopeDist(element.getElementSampleArr(), element.getModShift(), true);
                    }

                    /*
int[] tarrr = element.getElementSampleArr();
for(int ii=0;ii<tarrr.length;ii++)
System.out.println("=========" + tarrr[ii]);
System.out.println("=========---------");
tarrr = element.getElementRefArr();
for(int ii=0;ii<tarrr.length;ii++)
System.out.println("=========" + tarrr[ii]);
System.out.println("=========---------");
tarrr = element.getHeavyArr();
for(int ii=0;ii<tarrr.length;ii++)
System.out.println("=========" + tarrr[ii]);
element.printComposition();


//		    System.out.println("===========" + conf.getQuantType());
/*
		    System.out.println("===========" + element.getElementSampleArr());

		    //print element composition
		    System.out.println("sample===========");
		    int[] tempArr = element.getElementSampleArr();
		    for(int t : tempArr)
			System.out.println(t);
		    System.out.println("ref===========");
		    tempArr = element.getElementRefArr();
		    for(int t : tempArr)
                     */
                    peptideEle = this.createXmlChroPeptideTitle(true, peptide); //true is for full scan

                    /*
			    System.out.println("rarr===========");
			    System.out.println("rarr===========");
			    System.out.println("rarr===========");
			    System.out.println("rarr===========");
			int[] eleArr = element.getElementSampleArr();
			for(int ii:eleArr)
			    System.out.println("rarr===========" + ii);

			int[] eleArr = element.getElementSampleArr();
			for(int ii:eleArr)
			    System.out.println("sarr===========" + ii);
			eleArr = element.getElementRefArr();
			    System.out.println("===========");
                     */
                    if (conf.isCalculateEnrich()) {
                        double startEnrich = conf.getStartEnrich();
                        double endEnrich = conf.getEndEnrich();
                        double enrichmentMaxDeviation = conf.getEnrichmentMaxDeviation();

                        refDist = new edu.scripps.pms.census.util.N15EnrichmentCalc(element.getElementRefArr(), element.getModShift(), medianRatio, startEnrich, endEnrich, enrichmentMaxDeviation, iFile, scanNum, chargeState);
                        peptideEle.setAttribute("enrichment", String.valueOf(refDist.getFixedEnrichRatio()));
                        peptideEle.setAttribute("calc_enrich", String.valueOf(refDist.getBestEnrichRatio()));
                        peptideEle.setAttribute("enrich_corr", String.valueOf(refDist.getEnrichCorr()));
                        peptideEle.setAttribute("bestEnrichDelCN", String.valueOf(refDist.getBestEnrichDelCN()));
                        peptideEle.setAttribute("corrOnePlus", String.valueOf(refDist.getCorrOnePlus()));
                        peptideEle.setAttribute("corrOneMinus", String.valueOf(refDist.getCorrOneMinus()));

                    } else {
                        try {
                            if ("15N".equals(conf.getQuantType())) {
                                refDist = new IsotopeDist15N(element.getElementRefArr(), element.getModShift(), false);
                            } else {
                                refDist = new IsotopeDist(element.getElementRefArr(), element.getModShift(), false);
                            }

                            if ("3plexMS1Labeling".equals(conf.getQuantType())) {
                                ref2Dist = new IsotopeDist(element.getHeavyArr(), element.getModShift(), false);
                            }

                        } catch (Exception invE) {
                            System.out.println("Isotope Distribution is out of boundary: " + pepSequence);
                            invE.printStackTrace();
                            percent += eachSeg;

                            if (null != progress) {
                                ChroProgressDialog.addMessageWithLine(progress, "Not Quantifiable peptide : " + pepSequence);
                                progress.setProgress((int) percent);
                            }
                            continue;
                        }

                        peptideEle.setAttribute("enrichment", String.valueOf(conf.getEnrichment()));
                    }

                    /*
			double[] dmlist = sampleDist.getHighMassList();
			for(double dd:dmlist) {
			    if(dd<=0)
				continue;
			    System.out.println("sam mass :\t" + dd);
			}

			    System.out.println("===========");
                            System.out.println("==========="+ refDist);
			dmlist = refDist.getHighMassList();
			for(double dd:dmlist) {
			    if(dd<=0)
				continue;
			    System.out.println("ref mass :\t" + dd);
			}*/
                    conf.setCalcSamAvgMass(sampleDist.getAvgMass());
                    conf.setCalcRefAvgMass(refDist.getAvgMass());
                    peptideEle.setAttribute("lightStartMass", String.valueOf(sampleDist.getStartMass()));
                    peptideEle.setAttribute("heavyStartMass", String.valueOf(refDist.getStartMass()));
                    peptideEle.setAttribute("lightAvgMass", String.valueOf(conf.getCalcSamAvgMass()));
                    peptideEle.setAttribute("heavyAvgMass", String.valueOf(conf.getCalcRefAvgMass()));
                    if ("3plexMS1Labeling".equals(conf.getQuantType())) {
                        conf.setCalcRef2AvgMass(ref2Dist.getAvgMass());
                        peptideEle.setAttribute("heavy2StartMass", String.valueOf(ref2Dist.getStartMass()));
                        peptideEle.setAttribute("heavy2AvgMass", String.valueOf(conf.getCalcRef2AvgMass()));
                    }

                    //String fileName = peptide.getFileName().substring(0, peptide.getFileName().indexOf(".")+1) + "ms1";
                    StringBuffer rangeKey = new StringBuffer();
                    rangeKey.append(protein.getLocus());
                    rangeKey.append(fileName);
                    rangeKey.append(peptide.getSequence().substring(2, peptide.getSequence().length() - 2));

                    SpecRange range = rangeGen.getSpecRange(rangeKey.toString());


                    if (null == range) {
                        int tmpScanNum = Integer.parseInt(peptide.getScanNum());
                        range = new SpecRange(tmpScanNum, tmpScanNum);
                        peptideEle.setAttribute("start_scan", peptide.getScanNum());
                        peptideEle.setAttribute("end_scan", peptide.getScanNum());
                    } else {
                        peptideEle.setAttribute("start_scan", String.valueOf(range.getMin()));
                        peptideEle.setAttribute("end_scan", String.valueOf(range.getMax()));
                    }

                    peptideEle.setAttribute("start_scan", String.valueOf(range.getMin()));
                    peptideEle.setAttribute("end_scan", String.valueOf(range.getMax()));
//		    int scanNum = Integer.parseInt(peptide.getScanNum());

                    keys = iFile.getKeys();
                    keyIndex = Arrays.binarySearch(keys, scanNum);

                    if (keyIndex < 0) //Cannot find index
                    {
                        keyIndex = -(++keyIndex); //Math.abs(++keyIndex);
                    }
                    if (keyIndex >= keys.length) {
                        keyIndex--;
                    }

                    Element chro = new Element("chro");
//harshil shah
                    Element iso = null;
                    Element theoMass = null;
                    Element theoIntensity = null;
//till here....

                    if (conf.isHighRes()) {

                        iso = new Element("iso");
                        theoMass = new Element("theo-mass");
                        theoIntensity = new Element("theo-int");

                        double[] samIsoArr = sampleDist.getHighMassList();
                        double[] samIsoIntArr = sampleDist.getRelabun();    //light intensities
                        double[] refIsoIntArr = refDist.getRelabun();       //heavy intensities

                        double[] refIsoArr = null; //refDist.getHighMassList();
                        double[] ref2IsoArr = null; //refDist.getHighMassList();

                        for (int i = 0; i < samIsoArr.length; i++) {
                            samIsoArr[i] = (samIsoArr[i] + chargeState * PROTON_MASS) / chargeState;
                            //System.out.println("sample :\t" + samIsoArr[i] + "\t" + sampleDist.getRelabun()[i]);
                            //			    System.out.println("s222 before=======================" + samIsoArr[i]);
                            //   System.out.println("sample :\t" + samIsoArr[i] + "\t" + sampleDist.getRelabun()[i]);
                        }

                        if (refDist instanceof N15EnrichmentCalc) {
                            refIsoArr = refDist.getBestEnrichMassArr();

                            if (null == refIsoArr) {
                                continue; //cannot find isotope peaks
                            }
                        } else if (refDist instanceof IsotopeDist) {

                            if (conf.isUseMassDiff()) {
                                refIsoArr = IsotopeDist.getHeavySilacDist(sampleDist.getHighMassList(), element.getElementRefArr());
                                for (int i = 0; i < refIsoArr.length; i++) {
                                    //					    					    System.out.println("r22 before=======================" + refIsoArr[i]);
                                    refIsoArr[i] = (refIsoArr[i] + chargeState * PROTON_MASS) / chargeState;
                                    // System.out.println("r33 before=======================" + refIsoArr[i]);
                                    //					    					    System.out.println("r33 before=======================" + refIsoArr[i]);
                                }
                            } else {

                                refIsoArr = refDist.getHighMassList();
                                for (int i = 0; i < refIsoArr.length; i++) {
                                    //					    					    System.out.println("r22 before=======================" + refIsoArr[i]);
                                    refIsoArr[i] = (refIsoArr[i] + chargeState * PROTON_MASS) / chargeState;
                                    //  System.out.println("ref m/z\t =======================" + refIsoArr[i]);
                                    //					    				System.out.println("ref :\t " + refIsoArr[i] + "\t" + refDist.getRelabun()[i]);
                                }

                                if ("3plexMS1Labeling".equals(conf.getQuantType())) {
                                    ref2IsoArr = ref2Dist.getHighMassList();
                                    for (int i = 0; i < ref2IsoArr.length; i++) {
                                        //				    					    System.out.println("r22 before=======================" + ref2IsoArr[i]);
                                        ref2IsoArr[i] = (ref2IsoArr[i] + chargeState * PROTON_MASS) / chargeState;
                                        //				    					    System.out.println("r22 before=======================" + ref2IsoArr[i]);
                                        //				    				System.out.println("ref :\t " + ref2IsoArr[i] + "\t" + ref2Dist.getRelabun()[i]);
                                    }

                                }
                            }

                          //make ref and sameple arr same length, if enrichment is now low.
                          if(conf.getEnrichment()==0 || conf.getEnrichment()>0.9) {
                            if(samIsoArr.length<refIsoArr.length) {
                              double[] tempArr = new double[samIsoArr.length];
                              for(int i=0;i<tempArr.length;i++)
                                 tempArr[i] = refIsoArr[i];

                              refIsoArr = tempArr;
                            }

                            if ("3plexMS1Labeling".equals(conf.getQuantType())) {
                              if(samIsoArr.length<ref2IsoArr.length) {
                                double[] tempArr = new double[samIsoArr.length];
                                for(int i=0;i<tempArr.length;i++)
                                  tempArr[i] = ref2IsoArr[i];

                                ref2IsoArr = tempArr;
                             }

                            }

                          }



                        }

                        /*
			       for(int i=0;i<samIsoArr.length;i++)
			       {
			       System.out.println("s before=======================" + samIsoArr[i]);
			       }
			       for(int i=0;i<refIsoArr.length;i++)
			       {
			       System.out.println("rr before=======================" + refIsoArr[i]);
			       }
                         */
                        //boolean isProline = true;
                        //double prolineMassDiff = 6.0138;
                        //                        System.out.println(pepSequence);
                        //       System.out.println("=======================================");
                        String cleanName = fileName;
                        if(fileName.startsWith("H") && FileFilterUtil.isHeavyFile(fileName, filePath))
                        {
                            cleanName = fileName.substring(1);
                        } else if(fileName.startsWith("M") && FileFilterUtil.isMediumFile(fileName, filePath)) {
                            cleanName = fileName.substring(1);
                        }



                        /*if(peptide.getSequence().equals("R.AVFVDLEPTVIDEVR.T") && peptide.getScanNumber() ==73284 )
                        {
                        }*/
                        String chroLine = "";
                        int ms1Scan= SpectrumUtil.getMS1Scan(ht,filePath+cleanName,peptide.getScanNumber());
                        CalcUtil.setScanToSearch(ms1Scan);
                        DisplayData.DisplayChroData chroData = new DisplayData.DisplayChroData();
                        if (conf.isUseProline() && pepSequence.substring(2, pepSequence.length() - 2).contains("P")) {
                            //if(conf.isUseProline() && pepSequence.substring(2, pepSequence.length()-2).contains("E") && !pepSequence.substring(2, pepSequence.length()-2).contains("P"){
                            int prolineCount = conf.getProlineCount();

                            if (prolineCount > 1) {
                                int tmpCount = 0;
                                char[] chArr = pepSequence.substring(2, pepSequence.length() - 2).toCharArray();

                                for (char c : chArr) {
                                    if (c == 'P') {
                                        tmpCount++;
                                    }
                                }

                                if (tmpCount > 1) {
                                    prolineCount = 2;
                                } else {
                                    prolineCount = 1;
                                }
                            }
                            double[] refProlineIsoArr = new double[refIsoArr.length * (prolineCount + 1)];

                            for (int i = 0; i < refIsoArr.length; i++) {
                                refProlineIsoArr[i] = refIsoArr[i];
                            }

                            switch (prolineCount) {
                                case 1:
                                    for (int i = refIsoArr.length; i < refProlineIsoArr.length; i++) {
                                        refProlineIsoArr[i] = refIsoArr[i - refIsoArr.length] + Configuration.PROLINE_SHIFT / chargeState;
                                    }
                                    break;
                                case 2:
                                    for (int i = refIsoArr.length; i < refIsoArr.length * 2; i++) {
                                        refProlineIsoArr[i] = refIsoArr[i - refIsoArr.length] + Configuration.PROLINE_SHIFT / chargeState;
                                        refProlineIsoArr[i + refIsoArr.length - 1] = refIsoArr[i - refIsoArr.length] + 2 * Configuration.PROLINE_SHIFT / chargeState;
                                    }

                                    break;
                            }
                            chroLine = CalcUtil.calculateFullMS(keyIndex, iFile, samIsoArr, refProlineIsoArr, range,chroData);
                            chro.setText(chroLine);
                        } else if ("3plexMS1Labeling".equals(conf.getQuantType())) {
                            chroLine =CalcUtil.calculateFullMS3Plex(keyIndex, iFile, samIsoArr, refIsoArr, ref2IsoArr);
                            chro.setText(chroLine);
                        } else {
                            chroLine =CalcUtil.calculateFullMS(keyIndex, iFile, samIsoArr, refIsoArr, range, chroData);
                            chro.setText(chroLine);
                        }
                       /* double low = samIsoArr[0];
                        double high = refIsoArr[refIsoArr.length-1];*/

                        double [] massArr = CalcUtil.getXyMassArr();
                        double [] intensityArr = CalcUtil.getXyIntArr();
                       // TDoubleArrayList massList = new TDoubleArrayList();
                       // TDoubleArrayList intList = new TDoubleArrayList();

                        double maxIntensity = Double.MIN_VALUE;


                        List<edu.scripps.pms.census.model.XYPoint> massInt=new ArrayList<>();
                        for(int i=0; i<massArr.length; i++)
                        {
                            double d = massArr[i];
                            double di = intensityArr[i];
                            edu.scripps.pms.census.model.XYPoint xypoint = new edu.scripps.pms.census.model.XYPoint();
                           // System.out.print(d+" "+di+" ");
                            xypoint.setX(d);
                            xypoint.setY((int)di);
                            massInt.add(xypoint);
                        }
                       // System.out.println("");

                        //long threshold = (long)(maxIntensity*MZTHESHOLD);
                        for(DisplayData.DisplayChroDataXY data: chroData.getData1()){
                            if(data.getY()>maxIntensity) maxIntensity = data.getY();
                        }
                        for(DisplayData.DisplayChroDataXY data: chroData.getData2()){
                            if(data.getY()>maxIntensity) maxIntensity = data.getY();
                        }

                        chroData.setStartRange(CalcUtil.getLowbound());
                        chroData.setEndRange(CalcUtil.getHighbound());
                        chroData.setThoMass(IsoCalcUtil.massWrite(samIsoArr, refIsoArr));
                        chroData.setXyPoints(massInt);
                        //chroData.setMaxIntensity((long)maxIntensity);
                        chroData.setMassAccuracy(conf.getMassTolerance());
                        //  displayChroData.setMassAccuracy(mass_tol);

                        FileUtil.writeJSON(gson.toJson(chroData), filePath + MASTER_JSON +File.separatorChar+ peptide.getSequence() + peptide.getFileName() + peptide.getScanNum() + peptide.getChargeState() + ".JSON");



                        StringBuffer isoText = new StringBuffer();
                        for (double value : samIsoArr) {
                            isoText.append(value + ",");
                        }
                        isoText.append(":");
                        for (double value : refIsoArr) {
                            isoText.append(value + ",");
                        }
                        iso.setText(isoText.toString());

                        iso.setText(IsoCalcUtil.calculateFullMS(keyIndex, iFile, samIsoArr, refIsoArr, range));
                        theoMass.setText(IsoCalcUtil.massWrite(samIsoArr, refIsoArr));
                        theoIntensity.setText(IsoCalcUtil.intensityWrite(samIsoIntArr, refIsoIntArr, samIsoArr.length, refIsoArr.length));
// Till here-- harshil Shah
                    } else {
                        //System.out.println("==" + pepSequence);
                        sampleStartMass = (sampleDist.getStartMass() + chargeState * PROTON_MASS) / chargeState - massTolerance;
                        sampleEndMass = (sampleDist.getEndMass() + chargeState * PROTON_MASS) / chargeState + massTolerance;
                        refStartMass = (refDist.getStartMass() + chargeState * PROTON_MASS) / chargeState - massTolerance;
                        refEndMass = (refDist.getEndMass() + chargeState * PROTON_MASS) / chargeState + massTolerance;

                        TDoubleArrayList refAdditionalMassArr = new TDoubleArrayList();

                        if (conf.isUseProline() && pepSequence.substring(2, pepSequence.length() - 2).contains("P")) {
                            int prolineCount = conf.getProlineCount();

                            if (prolineCount > 1) {
                                int tmpCount = 0;
                                char[] chArr = pepSequence.substring(2, pepSequence.length() - 2).toCharArray();

                                for (char c : chArr) {
                                    if (c == 'P') {
                                        tmpCount++;
                                    }
                                }

                                if (tmpCount > 1) {
                                    prolineCount = 2;
                                } else {
                                    prolineCount = 1;
                                }
                            }

                            switch (prolineCount) {
                                case 1:
                                    refAdditionalMassArr.add(refStartMass + Configuration.PROLINE_SHIFT / chargeState);
                                    refAdditionalMassArr.add(refEndMass + Configuration.PROLINE_SHIFT / chargeState);
                                    break;

                                case 2: //2 or more
                                    refAdditionalMassArr.add(refStartMass + Configuration.PROLINE_SHIFT / chargeState);
                                    refAdditionalMassArr.add(refEndMass + Configuration.PROLINE_SHIFT / chargeState);
                                    refAdditionalMassArr.add(refStartMass + 2 * Configuration.PROLINE_SHIFT / chargeState);
                                    refAdditionalMassArr.add(refEndMass + 2 * Configuration.PROLINE_SHIFT / chargeState);

                                    break;
                            }

                            conf.setAddtionalRefMassArr(refAdditionalMassArr.toNativeArray());
                            // System.out.println(prolineCount + " " + Configuration.PROLINE_SHIFT + " " +chargeState + " " + (Configuration.PROLINE_SHIFT/chargeState));
                            // System.out.println("=" + refStartMass + " " + refEndMass + " " + massTolerance);

                            // for(double d : refAdditionalMassArr.toNativeArray())
                            //     System.out.println(d);
                            chro.setText(CalcUtil.calculateFullMS(keyIndex, iFile, sampleStartMass, sampleEndMass, refStartMass, refEndMass, range));

                        } else {
                            chro.setText(CalcUtil.calculateFullMS(keyIndex, iFile, sampleStartMass, sampleEndMass, refStartMass, refEndMass, range));
                        }
                    }
                    /*
System.out.println( " " + peptide.getSequence());
System.out.println( " " + peptide.getChargeState());
System.out.println( " " + peptide.getScanNum());
System.exit(0);
                     */
                    peptideEle.addContent(chro);

//added Harshil Shah--------------------
                    if (iso != null) {
                        peptideEle.addContent(iso);
                    }
                    if (theoMass != null) {
                        peptideEle.addContent(theoMass);
                    }
                    if (theoIntensity != null) {
                        peptideEle.addContent(theoIntensity);
                    }
//till here-- harshil SHah-----------------
                    proteinEle.addContent(peptideEle);

                    percent += eachSeg;

                    System.out.print(pepCount);
                    System.out.print("/");
                    System.out.print(redundantPeptideNum);
                    System.out.print(" peptides, ");
                    System.out.print((int) percent);
                    System.out.print(" % is complete\r");

                    if (conf.isPrintLog()) {
                        p.println((int) percent + "\t" + pepCount + "\t" + redundantPeptideNum);
                    }

                    if (null != progress) {
                        progress.setProgress((int) percent);
                    }
                }

                if (proteinEle.getChildren().size() > 0) {

                    Element redunEle = new Element("redundant");

                    int rPepCount = 0;

                    for (Iterator<Protein> rItr = aList.iterator(); rItr.hasNext();) {
                        rPepCount++;
                        Protein rp = rItr.next();

                        Element rpEle = new Element("protein");
                        rpEle.setAttribute("locus", rp.getLocus());
                        rpEle.setAttribute("seq_ct", rp.getSeqCount());
                        rpEle.setAttribute("spec_ct", rp.getSpectrumCount());
                        rpEle.setAttribute("seq_cov", rp.getSeqCoverage());
                        rpEle.setAttribute("length", rp.getLength());
                        rpEle.setAttribute("molwt", rp.getMolWt());
                        rpEle.setAttribute("pi", rp.getPI());
                        rpEle.setAttribute("val", rp.getValidation());

                        try {
                            rpEle.setAttribute("desc", rp.getDescription());
                        } catch (org.jdom.IllegalDataException ide) {
                            rpEle.setAttribute("desc", StringUtil.removeIsoControlChar(rp.getDescription()));
                        }

                        redunEle.addContent(rpEle);
                    }

                    if (rPepCount > 0) {
                        proteinEle.addContent(redunEle);
                    }

                    rootEle.addContent(proteinEle);
                    aList.clear();
                } else {
                    aList.add(protein);
                    //	    set.add(protein.getLocus());
                }
            }

            Document doc = new Document(rootEle);
            OutputStream os = new FileOutputStream(filePath + "census_chro.xml");
            XMLOutputter outputter = new XMLOutputter();
            outputter.setFormat(Format.getPrettyFormat());
            outputter.output(doc, os);
            os.close();

            System.out.println("\n100% complete");
        } catch (IOException e) {
            System.out.println("IO Error while generating msms chro file : " + e);
            e.printStackTrace();
            throw new IOException(e.toString());
        } catch (Exception e) {
            System.out.println("Error while generating msms chro file : " + e);
            e.printStackTrace();
            throw new Exception(e.toString());
        } finally {
            if (null != p) {
                p.close();
            }

            if (null != out) {
                out.close();
            }

            //Close all random files
            for (Enumeration e = ht.keys(); e.hasMoreElements();) {
                iFile = ht.get(e.nextElement());

                if (null != iFile) {
                    iFile.close();
                }
            }

        }

        System.out.println((System.currentTimeMillis() - startTime) * 0.001 + " seconds taken");
    }

    public static void setConfiguration(Configuration con) {
        conf = con;
    }

    public static boolean isHeavyFile(String file, String path) {

        if (!file.startsWith("H")) {
            return false;
        }

        String lightfile = path + "/" + file.split("\\.")[0] + ".ms2";
        //String lightfile = path + "/" + file.substring(1, file.length());
        File lf = new File(lightfile);

        if (lf.exists()) {
            return true;
        } else {
            return false;
        }

    }

    public static boolean isMediumFile(String file, String path) {

        if (!file.startsWith("M")) {
            return false;
        }

        // String lightfile = path + "/" + file.substring(1, file.length());
        String lightfile = path + "/" + file.split("\\.")[0] + ".ms2";

        File lf = new File(lightfile);

        if (lf.exists()) {
            return true;
        } else {
            return false;
        }

    }

    public static boolean isLightFile(String file, String path) {

        if (!file.startsWith("L")) {
            return false;
        }

        //String lightfile = path + "/" + file.substring(1, file.length());
        String lightfile = path + "/" + file.split("\\.")[0] + ".ms2";

        File lf = new File(lightfile);

        if (lf.exists()) {
            return true;
        } else {
            return false;
        }

    }

}
