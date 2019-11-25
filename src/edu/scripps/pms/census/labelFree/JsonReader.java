/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.scripps.pms.census.labelFree;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.json.simple.JSONObject;
import org.json.simple.JSONArray;
import org.json.simple.parser.ParseException;
import org.json.simple.parser.JSONParser;

public class JsonReader {
        public static void main(String[] args) {
}

    public JSONObject parser(String dir,Object p) {
            try {
                String filePath = dir+ File.separator +p;
                FileReader reader;
                reader = new FileReader(filePath);
                JSONParser parser = new JSONParser();
                JSONObject array = (JSONObject) parser.parse(reader);
                return array;
                
                }   
             
            catch (IOException | ParseException ex) {
                Logger.getLogger(JsonReader.class.getName()).log(Level.SEVERE, null, ex);
            }
            return null;           
}


public static ArrayList<String> getFilesBySuffix(String dir, String suffix) {
         ArrayList<String> suffixFiles = new ArrayList<String>();
         File currentDir = new File(dir);
         String [] files = currentDir.list();
     if(files==null) return suffixFiles;
       for(int i = 0; i < files.length; i++) {
           String s = files[i];
             if(s.endsWith(suffix)) {
                 suffixFiles.add(s);
             }
         }
            return suffixFiles;
    }
}
