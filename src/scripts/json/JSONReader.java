package scripts.json;

/**
 * Created by rpark on 1/5/17.
 */
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;

import java.io.FileReader;
import java.util.Iterator;
import java.util.List;

public class JSONReader {

    public static void main(String[] args) {

        JSONParser parser = new JSONParser();
        String file = "/home/rpark/test_data/targeted_lfree/JSON_OBJ/TARGETED_PEPTIDE.JSON";
        try {

            JSONObject jsonObject =  (JSONObject)parser.parse(new FileReader(file));

            JSONArray jsonArray = (JSONArray)jsonObject.get("peptideList");
            Iterator<JSONArray> eachItr = jsonArray.iterator();
            while (eachItr.hasNext()) {
                JSONArray pepList = (JSONArray)eachItr.next();
                for(Iterator<JSONObject> pepItr=pepList.iterator(); pepItr.hasNext(); ) {
                    JSONObject pep = pepItr.next();
                    System.out.println(pep.get("charge"));
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

}
