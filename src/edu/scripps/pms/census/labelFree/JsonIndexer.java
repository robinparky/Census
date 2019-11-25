/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.scripps.pms.census.labelFree;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;

/**
 *
 * @author Harshil
 */
public class JsonIndexer {
    
    public static final char TAB = '\t';
    private static final String START_TOKEN = "{";
    private static final String END_TOKEN = "}";
    private static final int
            CR  = 13,
            LF  = 10;
    
    /**
     * 
     * @param jsonIndexFile
     * @return Hashmap with key as protein Accession and value as IndexModel(Laocation)
     */
    public static HashMap<String,IndexModel> readJsonIndexFile(String jsonIndexFile)
    {
        HashMap<String,IndexModel> proteinToLocation = new LinkedHashMap<>();
        BufferedReader br = null;
        try {
             br = new BufferedReader(new FileReader(jsonIndexFile));
             String line =null;
             while((line = br.readLine()) != null)
             {
                 String words[] = line.split("\t");
                 if(words.length==0)
                     continue;
                 IndexModel indexModel = new IndexModel(Long.parseLong(words[1]), Long.parseLong(words[2]));
                 proteinToLocation.put(words[0], indexModel);
             }
            
            
        } catch (FileNotFoundException ex) {
            Logger.getLogger(JsonIndexer.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(JsonIndexer.class.getName()).log(Level.SEVERE, null, ex);
        }
        finally{
            try {
                br.close();
            } catch (IOException ex) {
                Logger.getLogger(JsonIndexer.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        
        return proteinToLocation;
    }
    
    /**
     * Creates the index file for the json File you pass, in the same directory.
     * @param jsonFile 
     */
    public static void createJsonIndex (String jsonFile)
    {
        
        BufferedReader br = null;
        BufferedWriter bw = null;
        try {
            String indexFile = jsonFile+".index";
            br = new BufferedReader(new FileReader(jsonFile));
            bw = new BufferedWriter(new FileWriter(indexFile));
            String line =null;
            long endIndex =0;
            int multiplier = 1;// check this one...
            long startindex =0;
            long newlineCount =0;
            boolean loopControl = false;
            boolean peptideListCheck = false;
            boolean singleCntrl = false;
//            List<Long> newLineIndex = new ArrayList<Object>
            while((line = br.readLine() )!= null)
            {
                newlineCount =0;
//                int tempIdnex =0;
//                tempIdnex = (line).getBytes().length;
                newlineCount++;
                StringBuffer proteinString = new StringBuffer();
                boolean newProtein = false;
                String protein = null;
//                if(line.startsWith("  {"))
                if(line.trim().startsWith(START_TOKEN))
                {

                    String newLine = br.readLine();
                    String prevLine = "";
                    proteinString.append(line);
//                    tempIdnex =0;
                    
                    while((newLine!= null)  && !singleCntrl)
                    {
                        if(newLine.contains("accession"))
                        {
                            protein = newLine.split(":")[1].substring(1,newLine.split(":")[1].length()-2);
                        }
                        proteinString.append(newLine);
                        newlineCount++;
                        prevLine = newLine;
                        
                         if(!peptideListCheck && newLine.trim().startsWith("\"peptideList\"")){
                        	 peptideListCheck =true;
                         }
                        
                        if(peptideListCheck && newLine.trim().equals("]")){
                        	loopControl = true;
                        }
                        
                        newLine=br.readLine();
                        
                        if(loopControl && newLine.trim().startsWith(END_TOKEN)){
                        	singleCntrl = true;
                        }
                        
                        
                    }
                    proteinString.append(newLine);
                    newlineCount++;
                    loopControl = false;
                    peptideListCheck = false;
                    singleCntrl = false;
                }
                
//                endIndex += proteinString.toString().getBytes().length + newlineCount*multiplier + tempIdnex;
                endIndex += proteinString.toString().getBytes().length + newlineCount*multiplier;
                
                if(protein != null)
                    bw.write(protein + "\t" + startindex + "\t" + endIndex+"\n");
                    //System.out.println(protein + "----"+ endIndex  + "----" + startindex + " == " + (endIndex-startindex));
                startindex =endIndex;
                
                
            }
            System.out.println("See the index file at " + indexFile);
        } catch (FileNotFoundException ex) {
            Logger.getLogger(JsonIndexer.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(JsonIndexer.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            try {
                br.close();
                bw.close();
                
            } catch (IOException ex) {
                Logger.getLogger(JsonIndexer.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }
    
    
    public static void main(String args[]) throws Exception
    {
    	String fileName = "/Users/OfirSoft/Ofir/sandbox/ip2/app/data/02/census_labelfree_out_10986.JSON";

//    	String fileName = "/Users/OfirSoft/Ofir/sandbox/ip2/app/data/01/org.JSON";
//    	createJsonIndex(fileName);
    	
    	String proteinString = getProteinString("P05387", fileName);
    	proteinString = proteinString.charAt(proteinString.length()-1) == ',' ? proteinString.substring(0, proteinString.length()-1) : proteinString;
    	witerFile(proteinString);
    	System.out.println("--------------------");
        
    }
    
	public static void witerFile(String content) {

		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new FileWriter(
					"/Users/OfirSoft/Ofir/sandbox/ip2/app/data/output.txt"));
			writer.write(content);

		} catch (IOException e) {
		} finally {
			try {
				if (writer != null)
					writer.close();
			} catch (IOException e) {
			}
		}

	}
    
    
    /**
     * DEBUG purpose
     * @param proteinName
     * @param fileName
     * @return 
     */
    public static String getProteinString(String proteinName,String fileName) 
    {
        HashMap<String ,IndexModel> proteinTable = readJsonIndexFile(fileName+".index");
        String proteinString =null;
        if(proteinTable.containsKey(proteinName))
        {
           proteinString= getProteinData(fileName,proteinTable.get(proteinName));
        }
        return proteinString;
    }
           
    private static String getProteinData(String inputFile,IndexModel location) 
    {
        StringBuffer sb = new StringBuffer();
        RandomAccessFile random = null;
        try {
            
            String indexFiel = inputFile+".index";
            random = new RandomAccessFile(inputFile, "r");
            random.seek(location.getStartIndex());
            int diff = (int) (location.getEndIndex()-location.getStartIndex());
            byte[] b = new byte[diff];
            random.readFully(b);
            for(byte byt : b)
            {
                sb.append((char)byt);
            }
            if(sb.charAt(sb.length()-2) == ',')
                sb.deleteCharAt(sb.length()-2);
//            System.out.println(sb.toString());
        } catch (IOException ex) {
            Logger.getLogger(JsonIndexer.class.getName()).log(Level.SEVERE, null, ex);
        }
        finally{
            try {
                random.close();
            } catch (IOException ex) {
                Logger.getLogger(JsonIndexer.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        return sb.toString();
    }
    
    
    /**
     * Debug purpose
     * 
     */
    private static void jsonIndexReader(String fileName)
    {
        HashMap<String ,IndexModel> proteinTable = readJsonIndexFile(fileName+".index");
        String proteinString =null;
        StringBuffer sb = new StringBuffer();
        try {
            IndexModel location = proteinTable.get("P98196-6");
            String indexFiel = fileName+".index";
            RandomAccessFile random = new RandomAccessFile(fileName, "r");
            random.seek(location.getStartIndex());
            int diff = (int) (location.getEndIndex()-location.getStartIndex());
            byte[] b = new byte[diff];
            random.readFully(b);
            for(byte byt : b)
            {
                sb.append((char)byt);
            }
            if(sb.charAt(sb.length()-2) == ',')
                sb.deleteCharAt(sb.length()-2);
//            System.out.println(sb.toString());
        } catch (IOException ex) {
            Logger.getLogger(JsonIndexer.class.getName()).log(Level.SEVERE, null, ex);
        }
        System.out.println(sb.toString());
              
    }
    
    /**
     * If you dont know the location of the protein in the JSON file ... 
     * simply pass the proteinName and the jsonFile it will give you the protein,
     * @param proteinName
     * @param jsonFile
     * @return 
     */
    public static ChroJSONProteinModel getProteinObject (String proteinName,String jsonFile)
    {
        HashMap<String ,IndexModel> proteinTable = readJsonIndexFile(jsonFile+".index");
        String proteinString =null;
        if(proteinTable.containsKey(proteinName))
        {
           proteinString= getProteinData(jsonFile,proteinTable.get(proteinName));
        }
        
        JSONParser parser  = new JSONParser();
        ChroJSONProteinModel jsonModel = null;
        try {
            JSONObject proteinJsonObj = (JSONObject) parser.parse(proteinString);
//            System.out.println(proteinString);
             jsonModel = new ChroJSONProteinModel();
            jsonModel.readJsonObject(proteinJsonObj);
        } catch (ParseException ex) {
            Logger.getLogger(JsonIndexer.class.getName()).log(Level.SEVERE, null, ex);
        }
        return jsonModel;
    }
    
    
    /**
     * This will give you the protein from the json file you passed.
     * @param jsonFile
     * @param indexModel (location of the protein in the jsonFile
     * @return ChroJSONProteinModel
     */
    public static ChroJSONProteinModel getProteinObject (String jsonFile,IndexModel indexModel)
    {
        String proteinString =getProteinData(jsonFile,indexModel);
        JSONParser parser  = new JSONParser();
        ChroJSONProteinModel jsonModel = null;
        try {
        	proteinString = proteinString.charAt(proteinString.length()-1) == ',' ? proteinString.substring(0, proteinString.length()-1) : proteinString;
            JSONObject proteinJsonObj = (JSONObject) parser.parse(proteinString);
//            System.out.println(proteinString);
             jsonModel = new ChroJSONProteinModel();
            jsonModel.readJsonObject(proteinJsonObj);
        } catch (ParseException ex) {
            Logger.getLogger(JsonIndexer.class.getName()).log(Level.SEVERE, null, ex);
        }
        return jsonModel;
    }
    
}
