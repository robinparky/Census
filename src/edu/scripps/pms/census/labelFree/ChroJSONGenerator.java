/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package edu.scripps.pms.census.labelFree;

import JSONWriter.JsonWriter;
import edu.scripps.pms.census.conf.Configuration;
import edu.scripps.pms.census.hash.IndexedFile;
import edu.scripps.pms.census.hash.MSIndexFileCreator;
import edu.scripps.pms.census.io.ChroXmlReader;
import edu.scripps.pms.census.model.ChroData;
import edu.scripps.pms.census.model.ChroPeptide;
import edu.scripps.pms.census.model.ChroProtein;
import edu.scripps.pms.census.tools.IsotopeModel;
import edu.scripps.pms.census.tools.IsotopeTool;

import java.awt.BufferCapabilities.FlipContents;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Date;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;

/**
 * 
 * @author Harshil
 * @version $Id:
 */
public class ChroJSONGenerator {
	private String txttmpFile = null;
	private String configFile = "/cheezer_share/Harshil/data/jolene/labelfree_quant/census_config_labelfree_7841.xml";
	private String jsonFile = "";
	Hashtable<String, edu.scripps.pms.census.labelFree.MissedPeptide> seqToMissedPeptide = new Hashtable<>();
	private List<HashMap<String, ChroPeptide>> peptideMap = new ArrayList<>();// key
																				// is
																				// seqquence_chargeState
																				// --
																				// Its
	private List<String> msFileList = null;
	JSONArray proteinList = new JSONArray();
	private BufferedWriter bw = null;
	private Hashtable<String, IndexedFile> ms1ToIndexFile = new Hashtable<>();

	/**
	 * 
	 * @param configFile
	 * @param txttmpFile
	 * @param msFileList
	 * @param seqToMissedPeptide
	 * @param jsonFile
	 */
	public ChroJSONGenerator(String configFile, String txttmpFile,
			List<String> msFileList,
			Hashtable<String, MissedPeptide> seqToMissedPeptide, String jsonFile) {

		this.txttmpFile = txttmpFile;
		this.msFileList = msFileList;
		for (String msFile : msFileList) {
			File f = new File(msFile);
			peptideMap.add(ChroXmlReader.getPeptideMap(f.getParent()
					+ File.separator + "census_chro_temp.xml"));
		}

		for (String file : msFileList) {
			try {
				System.out.println("Reading index file..." + file);
				generateIndexFile(file + ".ms1");
				IndexedFile iFile = new IndexedFile(new File(file
						+ ".ms1.index"), file + ".ms1");
				ms1ToIndexFile.put(file, iFile);
			} catch (IOException ex) {
				Logger.getLogger(LabelFreeParser.class.getName()).log(
						Level.SEVERE, null, ex);
			}
		}
		this.seqToMissedPeptide = seqToMissedPeptide;
		this.configFile = configFile;
		this.jsonFile = jsonFile;

	}

	private void generateIndexFile(String ms1FIle) {
		File f = new File(ms1FIle + ".index");
		if (f.exists())
			return;
		else {

			try {
				MSIndexFileCreator.createIndexFile(ms1FIle);
			} catch (IOException ex) {
				Logger.getLogger(ChroJSONGenerator.class.getName()).log(
						Level.SEVERE, null, ex);
			}
		}
	}

	public static void main(String args[]) {
		List<String> msFileList = new ArrayList<>();

		msFileList
				.add("/cheezer_share/Harshil/data/jolene/projects2014_01_26_22_56781/20140122_MP23_HEK_PBS_1ug_BEH_35cm_140min_35ms_CID_1.ms1");
		// msFileList.add("/cheezer_share/Harshil/data/jolene/projects2014_01_26_22_56782/20140122_MP23_HEK_PBS_1ug_BEH_35cm_140min_35ms_CID_2.ms1");
		msFileList
				.add("/cheezer_share/Harshil/data/jolene/projects2014_03_03_23_58406/20140227_MP23_HEK_PBS_1ug_BEH_50cm_140min_35ms_CID_1.ms1");
		// msFileList.add("/cheezer_share/Harshil/data/jolene/projects2014_03_03_23_58407/20140227_MP23_HEK_PBS_1ug_BEH_50cm_140min_35ms_CID_2.ms1");
               
               
	}

	public void parse(int total_peptideCount) {
		BufferedReader br = null;

		FileWriter bwr = null;
		// boolean isFound = false;
		int counter = 0;
		System.out.println("Adding missed Peptide Informaiton..This might take some time.");
		TxttmpReader reader = new TxttmpReader();
		try {
			br = new BufferedReader(new FileReader(txttmpFile));
			bwr = new FileWriter(new File(txttmpFile).getParent()
					+ File.separator + "process.log");
			String currentLine = null;
			ChroProtein chroProtein = null;
			// JSONObject proteinJsonObj = null;
			JSONArray proteinArrayObj = new JSONArray();
			ChroJSONProteinModel proteinModel = new ChroJSONProteinModel();
			boolean isPreviousLineProtein = false;// True: if previous line is
													// proteinLine
			while ((currentLine = br.readLine()) != null) {
				List<ChroPeptide> chroPepList = new ArrayList<>();

				if (currentLine.startsWith("PLINE"))
					reader.parseProteinHeader(currentLine);
				else if (currentLine.startsWith("SLINE"))
					reader.parsePeptideHeader(currentLine);
				else if (currentLine.startsWith("P\t")) {

					if (chroProtein != null && !isPreviousLineProtein) {
						// proteinList.add(proteinJsonObj);
						// writeJSON(proteinModel);
//						proteinArrayObj.add(proteinModel.getJSONObject());
						proteinArrayObj.add(writeJSON(proteinModel,getNewFileName()));
					}

					chroProtein = reader.parseProteinLine(currentLine);
					if (isPreviousLineProtein == true) {
						proteinModel.addRedundantProtein(chroProtein);
					} else {
						proteinModel = new ChroJSONProteinModel();
						proteinModel.setProtein(chroProtein);
					}

					isPreviousLineProtein = true;
					// if(chroProtein.getLocus().equals("P04350"))
					// System.out.println("");
					// proteinJsonObj = chroProtein.getJSONObj();
				}

				else if (currentLine.startsWith("S\t")) {
					isPreviousLineProtein = false;
					chroPepList = reader.parsePeptideLine(currentLine);

					// if(chroPepList.get(0).getSequence().equals("R.IMNTFSVVPSPK.V"))
					// System.out.println("");
					fillMissedPeptides(chroPepList);

					// generatePeptidesJSONTag(chroPepList,proteinJsonObj);
					generatePeptidesJSONTag(chroPepList);
					proteinModel.addPeptideGroup(chroPepList);
					System.out.print("------" + counter++ + "/"
							+ total_peptideCount + "----"
							+ chroPepList.get(0).getSequence() + "\r");
					bwr.write(counter + "\t" + total_peptideCount + "\n");
				}
			}

//			proteinArrayObj.add(proteinModel.getJSONObject());
			proteinArrayObj.add(writeJSON(proteinModel,getNewFileName()));
			
			writeJSON(proteinArrayObj);
			
			bwr.write("DONE");

			System.out.println(proteinList.toJSONString());

		} catch (Exception ex) {
			Logger.getLogger(ChroJSONGenerator.class.getName()).log(
					Level.SEVERE, null, ex);
			ex.printStackTrace();
		} finally {
			try {
				if (null != bw)
					bw.close();
				if (null != bwr)
					bwr.close();
				if (null != br)
					br.close();
			} catch (IOException ex) {
				ex.printStackTrace();
				Logger.getLogger(ChroJSONGenerator.class.getName()).log(
						Level.SEVERE, null, ex);
			}
		}
		// String fname =new File(this.txttmpFile).getParentFile()+
		// File.separator + "census_labelfree_out_7841.json";

		System.out.println("check file at " + this.jsonFile);
		System.out.print("Generating the IndexFile for json.\r");

//		JsonIndexer.createJsonIndex(this.jsonFile);
		// System.out.println("Check the index file at "+
		// this.jsonFile+".index");
	}

	/**
	 * Write JSON Array Object
	 * 
	 * @param proteinJsonObj
	 */
	public void writeJSON(JSONArray proteinJsonObj) {
		String fname = null;
		System.out.println("Write Protein Info JSON -  Start");
		try {
			BufferedWriter writer = null;
			fname = this.jsonFile;
			File f = new File(fname);
			if (f.exists()){
				f.delete();
			}
			writer = new BufferedWriter(new FileWriter(fname, true));
			writer.write(JsonWriter.formatJson(proteinJsonObj.toJSONString()));
			writer.write("\n");
			writer.close();
		} catch (IOException ex) {
			Logger.getLogger(ChroJSONGenerator.class.getName()).log(
					Level.SEVERE, null, ex);
		}
		System.out.println("Write Protein End");
	}
	
	/*

	public void writeJSON(ChroJSONProteinModel proteinModel) {
		JSONObject proteinJsonObj = proteinModel.getJSONObject();
		String fname = null;

		try {
			if (bw == null) {
				// fname =new File(this.txttmpFile).getParentFile()+
				// File.separator + "census_labelfree_out_7841.json";
				fname = this.jsonFile;
				File f = new File(fname);
				if (f.exists())
					f.delete();
				bw = new BufferedWriter(new FileWriter(fname, true));
			}

			bw.write(JsonWriter.formatJson(proteinJsonObj.toJSONString()));
			bw.write("\n");
		} catch (IOException ex) {
			Logger.getLogger(ChroJSONGenerator.class.getName()).log(
					Level.SEVERE, null, ex);
		}

	}
	
	*/
	
	@SuppressWarnings("unchecked")
	public JSONObject writeJSON(ChroJSONProteinModel proteinModel,String fileName) {
		
		System.out.println("writeJSON : ->"+ fileName + " ");
		
		this.showMemoryDetails();
		JSONObject accessionInfo = new JSONObject();
		JSONObject proteinJsonObj = proteinModel.getJSONObject();
		@SuppressWarnings("static-access")
		String accession = (String) ((JSONObject)proteinJsonObj.get(proteinModel.PROTEIN)).get("accession");
		String desc = (String) ((JSONObject)proteinJsonObj.get(proteinModel.PROTEIN)).get("desc");
		accessionInfo.put("accession", accession);
		accessionInfo.put("desc", desc);
		String fname = null;

		try {
			
			BufferedWriter writer = null;
				
			fname = new File(this.txttmpFile).getParent() +  File.separator + "JSON_OBJ";
			
			makeDir(fname);
			
			fname = fname+  File.separator + getNewFileName(fname);
			
			System.out.println(accession + " : FileName: " + fname);
			
			accessionInfo.put("jsonPath", fname);
			
			File f = new File(fname);
			if (f.exists())
				f.delete();
			writer = new BufferedWriter(new FileWriter(fname, true));

			writer.write(JsonWriter.formatJson(proteinJsonObj.toJSONString()));
			
			writer.close();
			
		} catch (IOException ex) {
			Logger.getLogger(ChroJSONGenerator.class.getName()).log(
					Level.SEVERE, null, ex);
		}
		
		this.showMemoryDetails();
		
		System.out.println("Write Peptide End");
		
		return accessionInfo;

	}
	
	public static void makeDir(String targetDir) {
		File dataFile = new File(targetDir);
		if (!dataFile.exists()) {
			dataFile.mkdirs();
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

	
	Date date = new Date();
	Calendar calendar = Calendar.getInstance();
	
	public String getNewFileName(){
		
		calendar.setTime(date);
		
		String newName = String.valueOf(calendar.getTimeInMillis());
		return newName;
	}
	
	public static String getNewFileName(String filePath) {
		String dir = filePath;
		int i = 0;
		String leadingZeros = "00000000";
		String fileName = "";
		File file = null;
		do {
			fileName = "" + i++;
			fileName = leadingZeros.substring(Math.min(fileName.length(),
					leadingZeros.length() - 1))
					+ fileName
					+ "."
					+ "JSON";
			file = new File(dir + File.separator + fileName);
		} while (file.exists());
		return fileName;
	}

	private void fillMissedPeptides(List<ChroPeptide> chroPepList) {
		// Find the NA entries in the chroPepList and fill it with the details
		// that we get from the previous task...
		String key = null;
                
		for (ChroPeptide chroPeptide : chroPepList) {
			//if (key == null && chroPeptide.getSequence() != null) {
                        if (chroPeptide.getSequence() != null) {
				key = chroPeptide.getSequence() + "_"
						+ chroPeptide.getChargeState();
				break;
			}
		}

		// DEGUB purpose
		// if(key.contains("K.TVLLLADQMISR.I"))
		// System.out.println("");

		for (int i = 0; i < chroPepList.size(); i++) {
			ChroPeptide currentPeptide = chroPepList.get(i);
			if (currentPeptide.getSequence() != null)// if not missed skip..
				continue;
			String tempkey = key + "_" + i;
			if (seqToMissedPeptide.containsKey(tempkey)) {
				generateMissedPeptideInfo(seqToMissedPeptide.get(tempkey));
				mergePeptide(seqToMissedPeptide.get(tempkey), currentPeptide);
				// seqToMissedPeptide.remove(tempkey);
			}
			// System.out.print("--Missed: " +
			// currentPeptide.getIsoTopeModelList().size() + "...");
		}
		// System.out.println("");

	}

	private void generateMissedPeptideInfo(MissedPeptide missedPeptide) {
		LabelFreeParser labelfreeParser = new LabelFreeParser();
		// DEBUG purpose
		// System.out.print("\nMissed-> ");
		labelfreeParser.getMissingIntensity(msFileList, configFile,
				missedPeptide, ms1ToIndexFile);
	}

	private void mergePeptide(MissedPeptide oldPeptide, ChroPeptide newPeptide) {
		newPeptide.setSequence(oldPeptide.getSequence());
		newPeptide.setChargeState(oldPeptide.getCs());
		newPeptide.setStartRange(String.valueOf(oldPeptide.getStartRange()));
		newPeptide.setEndRange(String.valueOf(oldPeptide.getEndRange()));
		newPeptide.setIsoTopeModelList(oldPeptide.getIsoTopeList());
		newPeptide.setMissedPeptide(true);
		newPeptide.setFileName(oldPeptide.getFileName());
		// newPeptide.setTotalIntensity(oldPeptide.getTotalIntensitySum());
		newPeptide.setAverageIntensity(oldPeptide.getTotalIntensitySum());
		oldPeptide.setIsoTopeList(new ArrayList<IsotopeModel>());
	}

	/**
	 * fill chro tag details for the Json Object...
	 * 
	 * @param list
	 *            of ChroPeptide
	 * @param proteinJsonObj
	 */
	private void generatePeptidesJSONTag(List<ChroPeptide> chroPepList) {
		// JSONArray peptideArray = (JSONArray) proteinJsonObj.get("peptides");

		// JSONArray experimentArray = new JSONArray();
		// List<ChroPeptide> experimentArray = new ArrayList<>();
		for (int expetimentIndex = 0; expetimentIndex < chroPepList.size(); expetimentIndex++) {
			ChroPeptide currentPeptide = chroPepList.get(expetimentIndex);

			if (currentPeptide.getSequence() != null)
				fillInChroData(currentPeptide, expetimentIndex);

			// experimentArray.add(currentPeptide);
		}
		// proteinJsonObj.addPeptideGroup(experimentArray);
		// peptideArray.add(experimentArray);

	}

	/**
	 * Reads the peptide information from the census_chro_temp.xml (for each
	 * experiment) different file... in a specific folder Generate the add all
	 * peaks ("isotopmeodelList" ) to that nonMissed Peptide.
	 * 
	 * @param chroPeptide
	 * @param experimentIndex
	 */
	private void fillInChroData(ChroPeptide currChroPeptide, int experimentIndex) {

		String key = currChroPeptide.getSequence() + "_"
				+ currChroPeptide.getChargeState();
		if (peptideMap.get(experimentIndex).get(key) == null)
			return;
		ChroPeptide newPeptide = peptideMap.get(experimentIndex).get(key);
		currChroPeptide.setDataList(newPeptide.getDataList());
		// DEBUG purpose
		// System.out.print("\nNot missed-> " );
		generatePeaks(currChroPeptide);

	}

	private void generatePeaks(ChroPeptide currentPeptide) {
		// /*Debug
		// if(currentPeptide.getSequence().equals("R.VYNWDVK.R"))
		// System.err.println("");
		// */

		try {
			double[] distArr = IsotopeTool.getIsotopePeaks(this.configFile,
					currentPeptide.getSequenceOnly(),
					currentPeptide.getChargeState());
			String msFile = currentPeptide.getFileName();
			for (String filePath : this.msFileList) {
				if (filePath.contains(msFile)) {
					msFile = filePath;
					break;
				}
			}

			if (!Configuration.getInstance().isReadConfigFile()) {
				Configuration.getInstance().readXMLParam(this.configFile);
			}

			// IndexedFile iFile = new IndexedFile(new
			// File(msFile+".ms1.index"), msFile+".ms1");
			if (!ms1ToIndexFile.containsKey(msFile))
				System.err.println("Error............");
			IndexedFile iFile = ms1ToIndexFile.get(msFile);
			// IndexedFile iFile = new IndexedFile(new File(msFile+".index"),
			// msFile);
			int startScan = ((ChroData) currentPeptide.getDataList().get(0))
					.getScanNum();
			int endScan = ((ChroData) currentPeptide.getDataList().get(
					currentPeptide.getDataList().size() - 1)).getScanNum();

			// List<IsotopeModel> isotopeModelList =
			// IsotopeTool.getIntensityChroByScanNumber(Integer.parseInt(currentPeptide.getStartRange()),
			// Integer.parseInt(currentPeptide.getEndRange()), distArr,
			// 0,massTolerance,iFile);
			List<IsotopeModel> isotopeModelList = IsotopeTool
					.getIntensityChroByScanNumber(startScan, endScan, distArr,
							0, iFile);
			// DEGUB puropose
			// if(isotopeModelList.size() == 0)
			// System.out.println("errrrrr");
			currentPeptide.setIsoTopeModelList(isotopeModelList);

		} catch (Exception ex) {
			Logger.getLogger(ChroJSONGenerator.class.getName()).log(
					Level.SEVERE, null, ex);
		}
	}

	/**
	 * DEBUG purpose
	 */
	public void writeJSON() {
		BufferedWriter bw = null;
		String fname = null;
		try {
			// fname= new File(this.txttmpFile).getParentFile()+ File.separator
			// + "census_labelfree_out_7841.json";
			fname = this.jsonFile;
			bw = new BufferedWriter(new FileWriter(fname));
			bw.write(JsonWriter.formatJson(proteinList.toJSONString()));
		} catch (IOException ex) {
			Logger.getLogger(ChroJSONGenerator.class.getName()).log(
					Level.SEVERE, null, ex);
		} finally {
			try {
				bw.close();
			} catch (IOException ex) {
				Logger.getLogger(ChroJSONGenerator.class.getName()).log(
						Level.SEVERE, null, ex);
			}
		}

	}
}
