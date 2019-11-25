package edu.scripps.pms.census.labelFree.json;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

import JSONWriter.JsonWriter;
import edu.scripps.pms.census.labelFree.ChroJSONGenerator;

public class ChroJSONCreator {

	private GsonBuilder builder = null;

	private Gson gson = null;

	public ChroJSONCreator() {
		builder = new GsonBuilder();
		setGson(builder.create());
	}

	public void createJsonForProteinIndex(List<LabelFreeJSONProtein> proteinIndex, String jsonFileName) {

		String jsonData = "";

		jsonData = gson.toJson(proteinIndex);

		writeJSON(jsonData, jsonFileName);

	}

	public void createJsonForPeptideIndex(LabelFreeJSONPeptide peptideIndex, BufferedWriter bw) {

		String jsonData = "";

		jsonData = gson.toJson(peptideIndex);

		writeJSON(jsonData, bw);

	}

	public void writeJSON(String jsonData, String jsonFileName) {

		try {
			BufferedWriter writer = null;

			File f = new File(jsonFileName);
			if (f.exists()) {
				f.delete();
			}
			writer = new BufferedWriter(new FileWriter(jsonFileName, true));
			writer.write(JsonWriter.formatJson(jsonData));
			writer.write("\n");
			writer.close();
		} catch (IOException ex) {
			Logger.getLogger(ChroJSONGenerator.class.getName()).log(Level.SEVERE, null, ex);
		}
	}

	public void writeJSON(String jsonData, BufferedWriter writer) {

		try {


			writer.write(JsonWriter.formatJson(jsonData));
			writer.write("\n");
			//writer.close();
		} catch (IOException ex) {
			Logger.getLogger(ChroJSONGenerator.class.getName()).log(Level.SEVERE, null, ex);
		}
	}

	
	public void createJsonForPeptideList(List<List<LabelFreeJSONPeptide>> peptideList, String jsonFileName) {
		AccessionJSON accessionJSON = new AccessionJSON();
		accessionJSON.setPeptideList(peptideList);
		String jsonData = gson.toJson(accessionJSON);

		writeJSON(jsonData, jsonFileName);
	}
	
	public static void main(String[] args) {

		String jsonData = "";

		GsonBuilder builder = new GsonBuilder();
		Gson gson = builder.create();

		AccessionJSON accessionJSON = new AccessionJSON();

		LabelFreeJSONProtein protein = new LabelFreeJSONProtein();

		protein.setAccession("P60174");
		protein.setDesc("Triosephosphate isomerase OS=Homo sapiens GN=TPI1 PE=1 SV=2 ");

		accessionJSON.setProtein(protein);

		List<List<LabelFreeJSONPeptide>> peptideList = new ArrayList<>();

		List<LabelFreeJSONPeptide> sampleGroup = new ArrayList<>();

		LabelFreeJSONPeptide peptide = new LabelFreeJSONPeptide();

		peptide.setPROFILE_SCORE("PROFILE_SCORE");
		peptide.setChro_iso("Chro-data");

		sampleGroup.add(peptide);

		peptide = new LabelFreeJSONPeptide();

		peptide.setPROFILE_SCORE("PROFILE_SCORE1");
		peptide.setChro_iso("Chro-data1");

		sampleGroup.add(peptide);

		peptideList.add(sampleGroup);

		sampleGroup = new ArrayList<>();

		peptide = new LabelFreeJSONPeptide();

		peptide.setPROFILE_SCORE("PROFILE_SCORE");
		peptide.setChro_iso("Chro-data");

		sampleGroup.add(peptide);

		peptide = new LabelFreeJSONPeptide();

		peptide.setPROFILE_SCORE("PROFILE_SCORE1");
		peptide.setChro_iso("Chro-data1");

		sampleGroup.add(peptide);

		peptideList.add(sampleGroup);

		accessionJSON.setPeptideList(peptideList);

		jsonData = gson.toJson(accessionJSON);

		try {
			JsonWriter.formatJson(jsonData);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public Gson getGson() {
		return gson;
	}

	public void setGson(Gson gson) {
		this.gson = gson;
	}

}
