/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.scripps.pms.census.labelFree;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Hashtable;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.jdom.Element;
import org.jdom.JDOMException;

import edu.scripps.pms.census.conf.Configuration;
import edu.scripps.pms.census.hash.IndexedFile;
import edu.scripps.pms.census.model.ChroPeptide;
import edu.scripps.pms.census.model.ChroProtein;
import edu.scripps.pms.census.tools.IsotopeModel;
import edu.scripps.pms.census.tools.IsotopeTool;
import edu.scripps.pms.util.spectrum.Range;
import edu.scripps.pms.util.stats.GaussianFitting;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Set;

/**
 * Step 1:- Run the Census.java to genearte the labelfreeout_id.txttmp file
 * using below command -c
 * /home/harshil/data/jolene/labelfree_quant/census_config_labelfree_7841.xml
 * -aa -of
 * /home/harshil/data/jolene/labelfree_quant//census_labelfree_out_7841.txt
 * ignore the -aa option....
 * 
 * OR
 * 
 * java -Xmx16G -cp "Census.jar:lib/*" edu.scripps.pms.census.Census -c
 * /home/harshil/on_data/jolene/labelfree_quant/census_config_labelfree_7841.xml
 * -of
 * /home/harshil/on_data/jolene/labelfree_quant//census_labelfree_out_7841.txt
 * 
 * Step 2:- Run LabelFreeParser.java with below configuration. Where the 1st
 * argument should be the xml file and 2nd should be the txttmp file. RUN
 * Arguments:
 * /home/harshil/data/jolene/labelfree_quant/census_config_labelfree_7841.xml
 * /home/harshil/data/jolene/labelfree_quant/census_labelfree_out_7841.txttmp
 * 
 * @author Harshil
 */
public class LabelFreeParser {

//	BufferedReader br = null;

	String txtTmpFile = null;
	List<List<ChroPeptide>> experiemntPeptideList = new ArrayList<>();
	Hashtable<String, List<Integer>> seqToIndex = new Hashtable<>();
	Hashtable<String, Double> seqToStartTime = new Hashtable<>();
	Hashtable<String, Double> seqToEndTime = new Hashtable<>();
	Hashtable<String, MissedPeptide> seqToMissedPeptide = new Hashtable<>();// key->
																			// seq_cs_experimentIndex
																			// where
																			// experiemtnIndex
																			// start
																			// at
																			// 0
	// List<String> msFileList = new ArrayList<>();

	public LabelFreeParser(String txtTmpFile) {
		try {
			this.txtTmpFile = txtTmpFile;
//			br = new BufferedReader(new FileReader(txtTmpFile));

		} catch (Exception e) {
			Logger.getLogger(LabelFreeParser.class.getName()).log(Level.SEVERE,
					null, e);
		}
		
		/*catch (FileNotFoundException ex) {
			Logger.getLogger(LabelFreeParser.class.getName()).log(Level.SEVERE,
					null, ex);
		}*/
	}

	public LabelFreeParser() {

	}

	public static void printUsage() {
		System.out
				.println("Usage: LabelFreeParser census_config_labelfree_7841.xml census_labelfree_out_7841.txttmp jsonFilePath");
		System.out.println("Note: all the files are with their full path.");

	}

	public static void main(String args[]) {
		
	

		if (args.length < 3) {
			printUsage();
			return;
		}

		String configFile = args[0];
		String tmpFile = args[1];
		String jsonFile = args[2];
		
		
//		String configFile = "/Users/OfirSoft/Ofir/documents/10-IP2/labelfree_data/labelfree_589/census_config_labelfree_589.xml";
//		String tmpFile = "/Users/OfirSoft/Ofir/documents/10-IP2/labelfree_data/labelfree_589/census_labelfree_out_589.txttmp";
////		String tmpFile = "/Users/OfirSoft/Ofir/documents/10-IP2/labelfree_data/labelfree_589/test_589.txttmp";
//		String jsonFile = "/Users/OfirSoft/Ofir/documents/10-IP2/labelfree_data/labelfree_589/census_config_labelfree_589.json";
//		
		LabelFreeParser labelFreeParser = new LabelFreeParser();
		
		if(!labelFreeParser.checkInputFileExist(tmpFile,configFile)){
			System.out.println("Input File Not Found ");
			return;
		}
		
		labelFreeParser.setTxtTmpFile(tmpFile);
                
                labelFreeParser.getsetRange(configFile);
                
		labelFreeParser.regenerateTxtTmpFile();
		labelFreeParser.preProcess();

//		List<MissedPeptide> missedPeptideList = reader.getMissingEntries();

		System.out.println("");

		List<String> msFileList = labelFreeParser.getmsFiles(configFile);
		// MISSED PEPTIDE IS CONSISTANT WITH THE MAP TO THE LIST.....

		int total_peptideCount = labelFreeParser.getExperiemntPeptideList().get(0).size();

		ChroJSONGenerator jsonGenerator = new ChroJSONGenerator(configFile,
				tmpFile, msFileList, labelFreeParser.getSeqToMissedPeptide(), jsonFile);
                
		jsonGenerator.parse(total_peptideCount);
		
		System.out.println("");

	}
	
	private boolean checkInputFileExist(String tmpFile, String configFile){
		
		boolean isResult = false;
		
		try {
			
			File sTmpFile = new File(tmpFile); 
			File sConfigFile = new File(configFile);
			
			isResult = sTmpFile.exists() && sConfigFile.exists();
			
		} catch (Exception e) {
			
			System.out.println("checkInputFileExist Exception->"+e.getMessage());
			
		}
		
		return isResult;
	}
	
	private boolean regenerateTxtTmpFile(){
		
		boolean isResult = true;
		
		try {
			TxttmpReader txtReader = new TxttmpReader(getTxtTmpFile());
			txtReader.regenerateFile();
		} catch (Exception e) {
			isResult = false;
			System.out.println("regenerateTxtTmpFile->"+e.getMessage());
		}
		
		return isResult;
		
	}
	
	private void preProcess(){
		
		try {
			
			init();
			getIndexMap();
			getRangeMap();
			getMissingEntries();
			
		} catch (Exception e) {
			System.out.println("preProcess->"+e.getMessage());
		}
		
	}

	private void showMemoryDetails() {
		int mb = 1024 * 1024;
		Runtime runtime = Runtime.getRuntime();
		System.out.println("Used Memory: "
				+ (runtime.totalMemory() - runtime.freeMemory()) / mb
				+ " ,Free Memory:" + runtime.freeMemory() / mb
				+ " ,Total Memory: " + runtime.totalMemory() / mb
				+ " ,Max Memory: " + runtime.maxMemory() / mb);
	}

	private List<String> getmsFiles(String configFile) {
		// this will get the MsFilePath and name by combining the config file
		// path and global msfilesName...
		List<String> msFiles = new ArrayList<>();
		List<String> msFilePath = null;
		try {
			msFilePath = new ArrayList<>();
			// Document doc = new SAXBuilder().build(new File(configFile));
			// Element root = doc.getRootElement();
			
			Configuration conf = Configuration.getInstance();
			if (!conf.isReadConfigFile())
				conf.readXMLParam(configFile);
			Element root = conf.getRootConfEle();
			List<Element> sampleList = root.getChildren("sample");
			for (Element sample : sampleList) {
				List<Element> eachSampleList = sample
						.getChildren("each_sample");
				for (Element eachSample : eachSampleList) {
					List<Element> eachmsFilesList = eachSample
							.getChildren("ms_files");
					for (Element eachmsFile : eachmsFilesList) {
						msFilePath.add(eachmsFile.getChildText("file"));
					}
				}
			}

			for (int i = 0; i < msFilePath.size(); i++) {
				StringBuffer sb = new StringBuffer();
				String directory = msFilePath.get(i).split("[*]")[0];
				sb.append(directory);
				File folder = new File(directory);
				String fileName = null;

				for (File eachFile : folder.listFiles()) {
					if (eachFile.getName().endsWith("ms1")
							|| eachFile.getName().endsWith("MS1")) {
						fileName = eachFile.getName().split(".ms1")[0];
						break;
					}
				}

				sb.append(fileName);
				// sb.append(msFileList.get(i));
				msFiles.add(sb.toString());

			}
			System.out.println("msFile List Generated.....");
		} catch (JDOMException ex) {
			Logger.getLogger(LabelFreeParser.class.getName()).log(Level.SEVERE,
					null, ex);
		} catch (IOException ex) {
			Logger.getLogger(LabelFreeParser.class.getName()).log(Level.SEVERE,
					null, ex);
		} catch (Exception ex) {
			Logger.getLogger(LabelFreeParser.class.getName()).log(Level.SEVERE,
					null, ex);
		}
		return msFiles;

	}


	/**
	 * 
	 * @param keys
	 *            : keys has to sorted....
	 * @param searchValue
	 * @return
	 */
	private double getActualValue(double[] keys, double searchValue) {
		// Arrays.sort(keys);

		int keyIndex = Arrays.binarySearch(keys, searchValue);

		if (keyIndex < 0) // Cannot find index
		{
			keyIndex = -(++keyIndex); // Math.abs(++keyIndex);
		}
		if (keyIndex >= keys.length) {
			keyIndex--;
		}
		return keys[keyIndex];
	}

	private List<MissedPeptide> getMissingEntries() {
		System.out.println("Generating missing peptideList");
		for (int i = 0; i < experiemntPeptideList.size(); i++) {

			List<ChroPeptide> currentExperiment = experiemntPeptideList.get(i);
			for (int j = i + 1; j < experiemntPeptideList.size(); j++) {
				List<ChroPeptide> nextExperiemnt = experiemntPeptideList.get(j);
				String index = i + "_" + j;
				compare(currentExperiment, nextExperiemnt, index);

			}

		}
		System.out.println("MissedPeptide list generated...."
				+ seqToMissedPeptide.size());
		return new ArrayList<>(seqToMissedPeptide.values());
	}

	private void compare(List<ChroPeptide> currentExperiment,
			List<ChroPeptide> nextExperiemnt, String index) {
		List<String> currentSequenceList = new ArrayList<>();
		List<String> nextSequenceList = new ArrayList<>();
		// genearte the List of peptideSequence+ChargerState for each peptide
		// Group....
		for (ChroPeptide peptide : currentExperiment) {
			currentSequenceList.add(peptide.getSequence() + "_"
					+ peptide.getChargeState());
		}
		for (ChroPeptide peptide : nextExperiemnt) {
			nextSequenceList.add(peptide.getSequence() + "_"
					+ peptide.getChargeState());
		}

		generateRange(currentExperiment, nextExperiemnt, currentSequenceList,
				nextSequenceList, index.split("_")[1]);
		generateRange(nextExperiemnt, currentExperiment, nextSequenceList,
				currentSequenceList, index.split("_")[0]);

	}

	/**
	 * Checks each element of the currentExperiemnt is present in the
	 * nextExperiment or not... always it will find the missing one in the
	 * nextExperiment,.......
	 * 
	 * @param currentExperiment
	 * @param nextExperiment
	 * @param currentSequenceList
	 * @param nextSequenceList
	 * @param experimentIndex
	 */
	private void generateRange(List<ChroPeptide> currentExperiment,
			List<ChroPeptide> nextExperiment, List<String> currentSequenceList,
			List<String> nextSequenceList, String experimentIndex) {
		System.out.println("generateRange start Loop");
		for (int i = 0; i < currentSequenceList.size(); i++) {

			// System.out.print(i + "\t00\t");

			/*
			 * DEBUG purpose if
			 * (currentSequenceList.get(i).contains("K.NFGGGNTAWEEKTLSKYESSEIR.L"
			 * )) { System.err.println(""); }
			 * 
			 * System.out.println("Debug - 1");
			 */
			if (currentSequenceList.get(i).startsWith("null")) {
				continue;
			}
                        
			// System.out.print("11\t");
                        

			if (!nextSequenceList.contains(currentSequenceList.get(i))) {
				int previousIndex = -1;
				int nextIndex = -1;
				for (int j = 1;; j++) {
					if (previousIndex == -1	&& i - j >= 0 && nextSequenceList.contains(currentSequenceList.get(i - j))
							&& !currentSequenceList.get(i - j).startsWith("null")) {
						previousIndex = nextSequenceList
								.indexOf(currentSequenceList.get(i - j));
                                                
					}
					if (nextIndex == -1
							&& i + j < currentSequenceList.size()
							&& nextSequenceList.contains(currentSequenceList
									.get(i + j))) {
                                                
						nextIndex = nextSequenceList
								.indexOf(currentSequenceList.get(i + j));
					}
					// if(i-j<0 && previousIndex == -1)
					// {//we cant find the previous entry in the search.......
					// previousIndex =0;
					// }
					// if(i+j> currentSequenceList.size() && nextIndex == -1)
					// {//we cant find the next entry in the sequence search
					// nextIndex=currentSequenceList.size()-1;
					// }
					if ((previousIndex != -1 && nextIndex != -1)
							|| (i - j < 0 && i + j > currentSequenceList.size())) {
						break;
					}
				}

				// System.out.print("22\t");
				if (previousIndex == -1 && nextIndex != -1) {
					setRange(currentSequenceList.get(i), -1.0, nextExperiment
							.get(nextIndex).getRetentionTime(),
							experimentIndex, currentExperiment.get(i)
									.getRetentionTime());
				} else if (nextIndex == -1 && previousIndex != -1) {
					setRange(currentSequenceList.get(i),
							nextExperiment.get(previousIndex)
									.getRetentionTime(), -1.0, experimentIndex,
							currentExperiment.get(i).getRetentionTime());
				} else if (nextIndex == -1 && previousIndex == -1)
					setRange(currentSequenceList.get(i), -1.0, -1.0,
							experimentIndex, currentExperiment.get(i)
									.getRetentionTime());
				else {
					setRange(currentSequenceList.get(i),
							nextExperiment.get(previousIndex)
									.getRetentionTime(),
							nextExperiment.get(nextIndex).getRetentionTime(),
							experimentIndex, currentExperiment.get(i)
									.getRetentionTime());
				}

				// System.out.println("33");
				// DEBUG purpose
			/*	System.out.print("=====\t"
						+ (i * 100 / (double) currentSequenceList.size())
						+ "%\t\t" + currentSequenceList.get(i) + "\t"
						+ previousIndex + "\t" + nextIndex + "\r");
				// System.out.print(".");
                                */
			}
		}

		System.out.println("\ncomplete. End Loop");
	}

	private void setRange(String sequenceKey, double startTime, double endTime,
			String experimentIndex, double reteinTime) {

		if (seqToStartTime.containsKey(sequenceKey)) {
			seqToStartTime.put(sequenceKey,
					Math.max(seqToStartTime.get(sequenceKey), startTime));

			seqToEndTime.put(sequenceKey,
					Math.min(seqToEndTime.get(sequenceKey), endTime));
		} else {
			seqToStartTime.put(sequenceKey, startTime);
			seqToEndTime.put(sequenceKey, endTime);
		}
		// if(sequenceKey.contains("K.DLYANTVLSGGTTMYPGIADR.M"))
		// System.err.println("asdsda------K.DLYANTVLSGGTTMYPGIADR.M");
		MissedPeptide missedPeptide = new MissedPeptide(sequenceKey, startTime,
				endTime, Integer.parseInt(experimentIndex), reteinTime);
		// DEBUG purpose labelfree
		// System.out.print(endTime + " : " + startTime +"->"+ (endTime -
		// startTime));

		if (seqToMissedPeptide.containsKey(sequenceKey + "_" + experimentIndex)) {
			seqToMissedPeptide.get(sequenceKey + "_" + experimentIndex).merge(
					missedPeptide);
		} else {
			seqToMissedPeptide.put(sequenceKey + "_" + experimentIndex,
					missedPeptide);
		}
	}

	private Hashtable<String, String> getRangeMap() {
		Hashtable<String, String> seqToRange = new Hashtable<>();
		for (String key : seqToIndex.keySet()) {
			List<Integer> indexList = seqToIndex.get(key);
//			List<Double> startIndex = new ArrayList<>();
//			List<Double> endIndex = new ArrayList<>();
                        double max=-Integer.MIN_VALUE;
                        double min= Integer.MAX_VALUE;
			for (int i = 0; i < indexList.size(); i++) {
				int index = indexList.get(i);
				/*if (index != -1) {
					startIndex.add(Double.parseDouble(experiemntPeptideList
							.get(i).get(index).getStartRange()));
					endIndex.add(Double.parseDouble(experiemntPeptideList
							.get(i).get(index).getEndRange()));
				}*/
                                //System.out.println(""+Double.parseDouble(experiemntPeptideList.get(i).get(index).getStartRange()));
                                if (index != -1){
                                    if (Double.parseDouble(experiemntPeptideList.get(i).get(index).getStartRange())>max){
                                        max = Double.parseDouble(experiemntPeptideList.get(i).get(index).getStartRange());
                                    }
                                    if (Double.parseDouble(experiemntPeptideList.get(i).get(index).getEndRange())<min){
                                        min = Double.parseDouble(experiemntPeptideList.get(i).get(index).getEndRange());
                                    }
                                }
                                
                                
                                
			}
			seqToRange.put(key,max+"_"+min);
		}
		return seqToRange;
	}

	private Hashtable<String, List<Integer>> getIndexMap() {
		for (int i = 0; i < experiemntPeptideList.size(); i++) {
			List<ChroPeptide> currentExperiment = experiemntPeptideList.get(i);

			// List<Integer> indexList = new ArrayList<>();
			for (int j = 0; j < currentExperiment.size(); j++) {
				ChroPeptide peptide = currentExperiment.get(j);
				List<Integer> tempList = new ArrayList<>();
				if (peptide.getSequence() == null) {
					continue;
				}
				String key = peptide.getSequence() + "_"
						+ peptide.getChargeState();
				// if(seqToIndex.contains(key))
				// {
				tempList = seqToIndex.get(key);
				tempList.set(i, j);
				// }
			}
		}
		return seqToIndex;
	}

	public List<List<ChroPeptide>> init() {

		TxttmpReader reader = new TxttmpReader();
		
		List<ChroProtein> proteinList = new ArrayList<ChroProtein>();
		
		BufferedReader br = null;
		
		try {
			
			br = new BufferedReader(new FileReader(txtTmpFile));
			
			String currentLine = br.readLine();
			// readHeader();

			while (currentLine != null) {
				if (currentLine.startsWith("H")) {
					// parseHeader();
				} // else if (currentLine.startsWith("PLINE")) {
					// parseProteinHeader();
					// }
				else if (currentLine.startsWith("SLINE")) {
					reader.parsePeptideHeader(currentLine);
				} else if (currentLine.startsWith("PLINE")) {
					reader.parseProteinHeader(currentLine);
				}

				if (currentLine.startsWith("S\t")) {
					List<ChroPeptide> peptideExpetient = reader
							.parsePeptideLine(currentLine);
					List<Integer> tempIndex = new ArrayList<>();
					for (int i = 0; i < reader.totalExperiments; i++) {
						tempIndex.add(-1);
					}
					for (int i = 0; i < reader.totalExperiments; i++) {
						for (int j = 0; j < peptideExpetient.size(); j++) {
							ChroPeptide peptide = peptideExpetient.get(j);
							// if(peptide.getSequence()== null)
							// continue;
							seqToIndex.put(peptide.getSequence() + "_"
									+ peptide.getChargeState(), tempIndex);
						}
						if (experiemntPeptideList.size() == 0) {
							// List<Integer> index = new ArrayList<>();
							for (int j = 0; j < peptideExpetient.size(); j++) {
								List<ChroPeptide> temp = new ArrayList<>();
								ChroPeptide peptide = peptideExpetient.get(j);
								// if(peptide.getSequence()== null)
								// continue;
								temp.add(peptide);
								experiemntPeptideList.add(temp);
							}
							break;
						} else {
							experiemntPeptideList.get(i).add(
									peptideExpetient.get(i));
						}

					}
				} else if (currentLine.startsWith("P\t")) {
					reader.parseProteinLine(currentLine);
				}
				currentLine = br.readLine();
			}
			// this.msFileList = reader.getMsFileList();

		}catch (FileNotFoundException ex) {
			Logger.getLogger(LabelFreeParser.class.getName()).log(Level.SEVERE,
					null, ex);
		} catch (IOException ex) {
			Logger.getLogger(LabelFreeParser.class.getName()).log(Level.SEVERE,
					null, ex);
		}
		
		try {
			br.close();
		} catch (IOException ex) {
			Logger.getLogger(LabelFreeParser.class.getName()).log(Level.SEVERE,
					null, ex);
		}
		experiemntPeptideList = sortByRetentionTime(experiemntPeptideList);
		
		return experiemntPeptideList;
	}

	// public ChroProtein parseProteinLine(String currentLine)
	// {
	// ChroProtein protein = new ChroProtein();
	// String words[] = currentLine.split("\t");
	// if(accessionIndex != -1)
	// protein.setLocus(words[accessionIndex]);
	// if(descriptionIndex != -1)
	// protein.setDescription(words[descriptionIndex]);
	// if(pepCountIndex != -1)
	// protein.setPepCount(Integer.parseInt(words[pepCountIndex]));
	// return protein;
	// }

	// public void parseProteinHeader(String currentLine)
	// {
	// // ACCESSION DESCRIPTION SCOUNT_1 SCOUNT_2 SCOUNT_3 SCOUNT_4 PEP_COUNT
	// NORM_INTENSITY_1 NORM_INTENSITY_2 NORM_INTENSITY_3 NORM_INTENSITY_4
	// String words[] = currentLine.split("\t");
	//
	// for (int i = 0; i < words.length; i++)
	// {
	// if (words[i].equalsIgnoreCase("ACCESSION"))
	// accessionIndex= i;
	// else if (words[i].equalsIgnoreCase("DESCRIPTION"))
	// descriptionIndex= i;
	// else if (words[i].contains("SCOUNT_"))
	// scountIndexList.add(i);
	// else if (words[i].equalsIgnoreCase("PEP_COUNT"))
	// pepCountIndex = i;
	//
	// }
	// }

	private List<List<ChroPeptide>> sortByRetentionTime(
			List<List<ChroPeptide>> peptideList) {
		List<List<ChroPeptide>> newList = new ArrayList<>();
		for (List<ChroPeptide> currentExperiment : peptideList) {
			Collections.sort(currentExperiment, new Comparator<ChroPeptide>() {

				@Override
				public int compare(ChroPeptide o1, ChroPeptide o2) {
					if (o1.getRetentionTime() < o2.getRetentionTime())
						return -1;
					else if (o1.getRetentionTime() > o2.getRetentionTime())
						return 1;
					else
						return 0;
					// int value = (int) (o1.getRetentionTime() -
					// o2.getRetentionTime());
					// return value;
				}
			});
			newList.add(currentExperiment);
		}
		return newList;
	}

	// public List<ChroPeptide> parsePeptideLine(String currentLine) {
	// String words[] = currentLine.split("\t");
	// int wordCounter = 0;
	//
	// List<ChroPeptide> cPeptideList = new ArrayList<>();
	// for (int i = 0; i < totalExperiments; i++) {
	//
	// ChroPeptide currentPeptide = new ChroPeptide();
	// if (sequenceIndexList.size() > 0) {
	// if (words[sequenceIndexList.get(i)].equalsIgnoreCase("NA")) {
	// cPeptideList.add(currentPeptide);
	// continue;
	// } else {
	// currentPeptide.setSequence(words[sequenceIndexList.get(i)]);
	// }
	// }
	// if (fileNameIndexList.size() > 0) {
	// currentPeptide.setFileName(words[fileNameIndexList.get(i)]);
	// if ( !words[fileNameIndexList.get(i)].equalsIgnoreCase("na")) {
	// msFileList.set(i,words[fileNameIndexList.get(i)]);
	// }
	// }
	// if (scanIndexList.size() > 0) {
	// currentPeptide.setScanNum(Integer.parseInt(words[scanIndexList.get(i)]));
	// }
	// if (csIndexList.size() > 0) {
	// currentPeptide.setChargeState(words[csIndexList.get(i)]);
	// }
	// if (intensityIndexList.size() > 0) {
	// currentPeptide.setAverageIntensity(Double.parseDouble(words[intensityIndexList.get(i)]));
	// }
	// if (profileScoreIndexList.size() > 0) {
	// currentPeptide.setAnCompositeScore(Double.parseDouble(words[profileScoreIndexList.get(i)]));
	// }
	// if (mhPlusIndexList.size() > 0) {
	// currentPeptide.setMhPlus(words[mhPlusIndexList.get(i)]);
	// }
	// if (calcMHPlusIndexList.size() > 0) {
	// currentPeptide.setCalcMHplus(words[calcMHPlusIndexList.get(i)]);
	// }
	// if (totalIntensityIndexList.size() > 0) {
	// currentPeptide.setTotalIntensity(Double.parseDouble(words[totalIntensityIndexList.get(i)]));
	// }
	// if (xCorIndexList.size() > 0) {
	// currentPeptide.setXCorr(words[xCorIndexList.get(i)]);
	// }
	// if (dcnIndexList.size() > 0) {
	// currentPeptide.setDeltCN(words[dcnIndexList.get(i)]);
	// }
	// if (dMassIndexList.size() > 0) {
	// currentPeptide.setDeltMass(words[dMassIndexList.get(i)]);
	// }
	// if (sprankIndexList.size() > 0) {
	// currentPeptide.setSpRank(words[sprankIndexList.get(i)]);
	// }
	// if (spScoreIndexList.size() > 0) {
	// currentPeptide.setSpScore(words[spScoreIndexList.get(i)]);
	// }
	// if (redundancyIndexList.size() > 0) {
	// currentPeptide.setRedundancy(words[redundancyIndexList.get(i)]);
	// }
	// if (startIndexList.size() > 0) {
	// currentPeptide.setStartRange(words[startIndexList.get(i)]);
	// }
	// if (endIndexList.size() > 0) {
	// currentPeptide.setEndRange(words[endIndexList.get(i)]);
	// }
	// if (retentionIndexList.size() > 0) {
	// currentPeptide.setRetentionTime(Double.parseDouble(words[retentionIndexList.get(i)]));
	// }
	// if (ionInjectionIndexList.size() > 0) {
	// currentPeptide.setIonInjectionTime(Double.parseDouble(words[ionInjectionIndexList.get(i)]));
	// }
	//
	// cPeptideList.add(currentPeptide);
	// }
	//
	// return cPeptideList;
	// }

	public void getMissingIntensity(List<String> msFileList,
			String configFileName, MissedPeptide missedPeptide,
			Hashtable<String, IndexedFile> ms1ToIndexFile) {
            
                

		double retentionTimeWindow = -1;
		// gets the chro tag details from the ms1 file range using the
		// missedPeptide detiails.....
		try {
			Configuration conf = Configuration.getInstance();
                        
                       
			if (!conf.isReadConfigFile())
				conf.readXMLParam(configFileName);
                     
			retentionTimeWindow = conf.getRetentionTimeWindow();
		} catch (Exception e) {
			Logger.getLogger(LabelFreeParser.class.getName()).log(Level.SEVERE,
					null, e);
		}
                

                
               

		double searchDifference = missedPeptide.getEndRange()
				- missedPeptide.getStartRange();
		double startTime = missedPeptide.getStartRange();
		double endTime = missedPeptide.getEndRange();
		double searchWindow = missedPeptide
				.getRetentionEndTime(retentionTimeWindow)
				- missedPeptide.getRetentionStartTime(retentionTimeWindow);
		if (searchDifference > searchWindow || searchDifference < 0) {
			startTime = missedPeptide
					.getRetentionStartTime(retentionTimeWindow);
			endTime = missedPeptide.getRetentionEndTime(retentionTimeWindow);

		}
		// if(currentPeptide.getSequence().equals("K.VVLAYEPVWAIGTGK.T"))
		// System.err.println("");
		String msFile = msFileList.get(missedPeptide.getExperimentIndex());
		missedPeptide.setFileName(new File(msFile).getName());
		IndexedFile iFile = ms1ToIndexFile.get(msFile);
		if (startTime != -1 || endTime != -1) {
			// startTime =
			// getActualValue(iFile.getRetentonToScanMap().keys(),startTime);
			// endTime =
			// getActualValue(iFile.getRetentonToScanMap().keys(),endTime);
			startTime = getActualValue(iFile.getRetentionTimeSortedArr(),
					startTime);
			endTime = getActualValue(iFile.getRetentionTimeSortedArr(), endTime);
		}

                
		try {
			double[] distArr = IsotopeTool.getIsotopePeaks(configFileName,
					missedPeptide.getSequenceOnly(), missedPeptide.getCs());
                        System.out.println(System.currentTimeMillis());
			 
			/*
			 * if(currentPeptide.getStartRange() == -1) isotopeModelList =
			 * IsotopeTool
			 * .getIntensityChro(currentPeptide.getRetentionStartTime(
			 * RETTOLERANCE), endTime, distArr, 0,massTolerance,iFile); else if
			 * ( currentPeptide.getEndRange() == -1) isotopeModelList =
			 * IsotopeTool.getIntensityChro(startTime,
			 * currentPeptide.getRetentionEndTime(RETTOLERANCE), distArr, 0,
			 * massTolerance, iFile); else if (currentPeptide.getStartRange() >
			 * currentPeptide.getEndRange()) isotopeModelList =
			 * IsotopeTool.getIntensityChro
			 * (currentPeptide.getRetentionStartTime(RETTOLERANCE),
			 * currentPeptide.getRetentionEndTime(RETTOLERANCE), distArr, 0,
			 * massTolerance, iFile); else isotopeModelList =
			 * IsotopeTool.getIntensityChro(startTime,endTime , distArr, 0,
			 * massTolerance, iFile);
			 */

			List<IsotopeModel> isotopeModelList = IsotopeTool.getIntensityChro(startTime, endTime,
					distArr, 0, iFile);
                       System.out.println("End getIntensityChro :"+System.currentTimeMillis());
			// DEBUG purpose labelfree
			// String outputPath = "E:\\jolene\\output\\" +
			// missedPeptide.getSequence()+"_"+missedPeptide.getCs()+"_"+missedPeptide.getExperimentIndex()+".txt";
			// if(isotopeModelList.size() == 0)
			// System.out.println("ZERO entry found in spectra ms1 file");
			// System.out.println(outputPath+"\t"+isotopeModelList.size());

			// DEBUG purpose labelfree
			// System.out.println( " --" + counter++ + "/" + totalMissedPep +
			// "--" +missedPeptide.getFileName());

			// debug purpose
                       
                       
			missedPeptide.setIsoTopeList(isotopeModelList);
			missedPeptide.generateBasePeakScan();
                        
                        System.out.println("Start other part :"+System.currentTimeMillis());
                        
			// Gaussian Peak Finding approach
			// PeakRange range = getPeakRange(isotopeModelList,iFile);
			// Below is old approach for the PeakFinding...
			PeakRange range = LabelFreeCalcUtil.generateRange(
					ms1ToIndexFile.get(msFile), distArr,
					missedPeptide.getBasePeakScanNumber(), isotopeModelList);
                        

			missedPeptide.setStartRange(range.getStart());
			missedPeptide.setEndRange(range.getEnd());
			addNewIsotopeIntensity(range, isotopeModelList,
					ms1ToIndexFile.get(msFile), distArr);

			generateTotalIntensity(missedPeptide);
                        
                        System.out.println("End other part:"+System.currentTimeMillis());

		} catch (Exception ex) {
			Logger.getLogger(LabelFreeParser.class.getName()).log(Level.SEVERE,
					null, ex);
		}
                
                
               // long taskTimeMs  = System.currentTimeMillis( );
               // System.out.println("\n"+taskTimeMs);

	}

	private PeakRange getPeakRange(List<IsotopeModel> isoModelList,
			IndexedFile iFile) {
		double[] rtTimeArr = new double[isoModelList.size()];
		double[] intensityArr = new double[isoModelList.size()];
		for (int i = 0; i < isoModelList.size(); i++) {
			rtTimeArr[i] = isoModelList.get(i).getRetentionTime();
			intensityArr[i] = isoModelList.get(i).getIntensitySum();
		}

		Range range = GaussianFitting.getGaussianPeakRange(rtTimeArr,
				intensityArr);
		int startScan = iFile.getRetentonToScanMap().get(range.getLowBound());
		int endScan = iFile.getRetentonToScanMap().get(range.getHighBound());

		return new PeakRange(startScan, endScan);
	}

	private void generateTotalIntensity(MissedPeptide missedPeptide) {
		double total = 0;
		for (IsotopeModel isotopeModel : missedPeptide.getIsoTopeList()) {
			total += isotopeModel.getIntensitySum();
		}
		missedPeptide.setTotalIntensitySum(total);
	}

	/**
	 * This will add the intensity to the passed isotopeModelList, depending on
	 * the range we pass in the parameter.
	 * 
	 * @param range
	 * @param isotopeModelList
	 * @param iFile
	 * @param distArr
	 * @return the passed isotopeModelList is updated with new start and end
	 *         range.
	 */
	private void addNewIsotopeIntensity(PeakRange range,
			List<IsotopeModel> isotopeModelList, IndexedFile iFile,
			double[] distArr) {
		int actualStartRange = isotopeModelList.get(0).getScanNumber();
		int actualEndRange = isotopeModelList.get(isotopeModelList.size() - 1)
				.getScanNumber();

		if (actualStartRange > range.getStart())// means we dont have info about
												// start cut.....
		{
			List<IsotopeModel> startIsotopeList = IsotopeTool
					.getIntensityChroByScanNumber(range.getStart(),
							actualStartRange, distArr, 0, iFile);
			isotopeModelList.addAll(0, startIsotopeList);
		}
		if (actualEndRange < range.getEnd())// means we dont have info about end
											// cut.....
		{
			List<IsotopeModel> endIsotopeList = IsotopeTool
					.getIntensityChroByScanNumber(actualEndRange,
							range.getEnd(), distArr, 0, iFile);
			isotopeModelList
					.addAll(isotopeModelList.size() - 1, endIsotopeList);
		}

	}
     
        
        private void getsetRange(String configFileName){
            TxttmpReader reader = new TxttmpReader(txtTmpFile);
            
            
            BufferedReader br = null;
            

            try{
                Configuration conf = Configuration.getInstance();
            if (!conf.isReadConfigFile())
				conf.readXMLParam(configFileName);
                double retentionTimeWindow = conf.getRetentionTimeWindow();
                System.out.println(""+retentionTimeWindow);
                ProteinModel cProteinList = new ProteinModel();
                ChroPeptide cPeptide = new ChroPeptide();
                List<ChroPeptide>  cPeptideList = new ArrayList<>();
                ArrayList<Double> startRange = new ArrayList<>();
                ArrayList<Double> endRange = new ArrayList<>();
                br = new BufferedReader(new FileReader(txtTmpFile));
                String currentLine = null;
                reader.readWholeFile();
                while((currentLine = br.readLine()) != null){
                //    System.out.println(""+currentLine);
                    if(currentLine.startsWith("S\t")){
                        cPeptideList = reader.parsePeptideLine(currentLine);
                    }
                    cProteinList.addPeptideList(cPeptideList);
                }
                cPeptideList.clear();
                HashMap<String, List<ChroPeptide>> peptideMap = new LinkedHashMap<>();// Key is accession+chargeState
                peptideMap = cProteinList.getPeptideMap();
                Set<String> keys = peptideMap.keySet();
                for(String key: keys){
                    cPeptideList = (List<ChroPeptide>) peptideMap.get(key);
                    for (int i = 0; i < cPeptideList.size();i++)
                    {
                        cPeptide = cPeptideList.get(i);
                        startRange.add(cPeptide.getRetentionTime());
                        endRange.add(cPeptide.getRetentionTime());
                    }
                    
                        double max=-Integer.MIN_VALUE;
                        double min= Integer.MAX_VALUE;
			for (int i = 0; i < startRange.size(); i++) {
                                    if(startRange.get(i) == 0.0){
                                     continue;
                                    }
                               
                                    else { 
                                        if (startRange.get(i) > max){
                                        max = startRange.get(i);
                                        max+=retentionTimeWindow;
                                    }
                                        if (endRange.get(i) < min){
                                        min = endRange.get(i);
                                        min-= retentionTimeWindow;
                                    }
                                    }
                        }
                        System.out.println(":\t"+min+"\t"+max);
                        startRange.clear();
                        endRange.clear();
                }

                    
                
            }
        
            catch (Exception ex){
                  System.out.println(""+ex);
		}
            
        }
        
        
        
        
        
        

	public String getTxtTmpFile() {
		return txtTmpFile;
	}

	public void setTxtTmpFile(String txtTmpFile) {
		this.txtTmpFile = txtTmpFile;
	}

	public List<List<ChroPeptide>> getExperiemntPeptideList() {
		return experiemntPeptideList;
	}

	public void setExperiemntPeptideList(
			List<List<ChroPeptide>> experiemntPeptideList) {
		this.experiemntPeptideList = experiemntPeptideList;
	}

	public Hashtable<String, List<Integer>> getSeqToIndex() {
		return seqToIndex;
	}

	public void setSeqToIndex(Hashtable<String, List<Integer>> seqToIndex) {
		this.seqToIndex = seqToIndex;
	}

	public Hashtable<String, Double> getSeqToStartTime() {
		return seqToStartTime;
	}

	public void setSeqToStartTime(Hashtable<String, Double> seqToStartTime) {
		this.seqToStartTime = seqToStartTime;
	}

	public Hashtable<String, Double> getSeqToEndTime() {
		return seqToEndTime;
	}

	public void setSeqToEndTime(Hashtable<String, Double> seqToEndTime) {
		this.seqToEndTime = seqToEndTime;
	}

	public Hashtable<String, MissedPeptide> getSeqToMissedPeptide() {
		return seqToMissedPeptide;
	}

	public void setSeqToMissedPeptide(
			Hashtable<String, MissedPeptide> seqToMissedPeptide) {
		this.seqToMissedPeptide = seqToMissedPeptide;
	}
}
