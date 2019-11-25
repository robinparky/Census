package edu.scripps.pms.census.labelFree;

import edu.scripps.pms.census.conf.Configuration;
import edu.scripps.pms.census.model.SampleGroup;
import edu.scripps.pms.census.model.SampleModel;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;

import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static edu.scripps.pms.census.labelFree.LabelfreeTargetedUtil.createAUCTable;

/**
 * Created by Titus Jung titusj@scripps.edu on 11/8/18.
 */
public class LabelfreeTargetedTableBuilder {


    public static void main(String [] args) throws Exception {
        String configFile = args[0];
        String jsonPeptideListFileName = args[1];
        String outJSONFile = args[2];

        Configuration conf = Configuration.getInstance();
        // conf.setLabelfreeCheckChargeState(true);

        if (!conf.isReadConfigFile()) {
            conf.setLabelfree(true);
            conf.readXMLParam(configFile);
        }


        Map<String, Integer> sampleFileNameMap = new HashMap<>();
        List<SampleGroup> sampleGroupList = conf.getSampleGroupList();

        for (SampleGroup sampleGroup : sampleGroupList) {
            for (SampleModel sampleModel : sampleGroup.getSampleModelList()) {
                List<String> fnameList = sampleModel.getLabelfreeFilenameList();
                List<String> pathList = sampleModel.getPathList();


                for (int i = 0; i < fnameList.size(); i++) {
                    String eachFile = fnameList.get(i);
                    String sampleFileNameKey = sampleGroup.ID + eachFile;
                    sampleFileNameMap.put(sampleFileNameKey, sampleGroup.ID);
                }
            }
        }


        createAUCTable(jsonPeptideListFileName,outJSONFile,sampleFileNameMap);
    }

    public static void buildSampleFileMap(String path) throws IOException, ParseException {
        JSONParser parser=new JSONParser();
        Object obj = parser.parse(new FileReader(path));
        JSONObject jsonObject =  (JSONObject) obj;
        JSONArray peptideLst = (JSONArray) jsonObject.get("peptideList");

        for(int j=0; j<peptideLst.size(); j++) {

            JSONArray lst = (JSONArray) (peptideLst.get(j));
            String sequence = "";
            Map<String, List<String>> aucMap = new HashMap<>();
            Map<String, List<String>> fileNameMap = new HashMap<>();
            Map<String, List<String>> rtMap = new HashMap<>();
            Map<String, List<String>> peakMap = new HashMap<>();
            Map<String, List<String>> isoMap = new HashMap<>();


            String cs = "";
            for (int k = 0; k < lst.size(); k++) {
                JSONObject lstObj = (JSONObject) lst.get(k);
                JSONArray isotopeArray = (JSONArray) lstObj.get("isotopeArr");
                String auc = (String) lstObj.get("auc");

                sequence = (String) lstObj.get("seq");

                String fileName = (String) lstObj.get("file");
                String sampleName = (String) lstObj.get("sampleName");
                String groupName = (String) lstObj.get("rowId");
                String key = groupName + fileName;
            }
        }
    }



}
