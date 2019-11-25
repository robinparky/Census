/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.scripps.pms.census;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.HashMap;
import org.apache.commons.io.FilenameUtils;

/**
 *
 * @author rampuria
 */
public class SQTRegeneration {
    public static void main(String [] args){
        
        if(args.length < 1){
            System.out.println("USAGE: input the sqt file path");
            return;
        }
        
        
        String filename = args[0];
        SQTRegeneration sqt = new SQTRegeneration();
        sqt.regenerate(filename);
        
    }

    public void regenerate(String filename) {
        File file = new File(filename);
        String fileName = file.getParent() + File.separator + file.getName();
        String extension = FilenameUtils.getExtension(filename);

        File outfile = new File(file.getParent() + File.separator + "temp.txt");
        HashMap<String,Integer> scanMap = new HashMap<>();
        try{
            
            BufferedReader br2 = new BufferedReader(new FileReader(file));
            BufferedWriter bw  = new BufferedWriter(new FileWriter(outfile));
            String eachLine = br2.readLine();
            while(eachLine != null){
                
                if(eachLine.startsWith("S\t")){
                    String [] words = eachLine.split("\t");
                    scanMap.put(words[0]+"\t"+words[1]+"\t"+words[2],0);
                    eachLine=br2.readLine();
                }
                else{
                    eachLine=br2.readLine();
                }
            }
            br2.close();
            
            BufferedReader br = new BufferedReader(new FileReader(file));
            while((eachLine=br.readLine()) != null){
                
                if(eachLine.startsWith("S\t")){
                    String [] words = eachLine.split("\t");
                    int count = scanMap.get(words[0]+"\t"+words[1]+"\t"+words[2]);
                    //((Integer.parseInt(words[1])*10000)+counter)
                    bw.write("S\t"+((Integer.parseInt(words[1])*100000)+count)+"\t"
                            +((Integer.parseInt(words[2])*100000)+count)+"\t"+
                            words[3]+"\t"+words[4]+"\t"+words[5]+"\t"+words[6]+"\t"+words[7]+"\t"+words[8]+"\t"+words[9]+"\n");
                            count++;
                            scanMap.put(words[0]+"\t"+words[1]+"\t"+words[2], count);
                }
                
                else{
                   // if(eachLine.equals("L\tsp|Q06512|NOC4_YEAST\t278\tSIL.LILHKR.IIP")){
                  //      System.out.println("");
                  //  }
                    
                    bw.write(eachLine+"\n");
                }
            }
            bw.close();
            file.renameTo(new File(file.getParent()+File.separator+file.getName()+"_orig"));
            outfile.renameTo(new File(fileName));
            
            
            
            
        }catch(Exception e){
            e.printStackTrace();
        }

    }

    
    
}
