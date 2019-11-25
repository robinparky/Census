/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.scripps.dia;

import edu.scripps.pms.census.ChroGenerator;
import edu.scripps.pms.census.ChroProgressDialog;
import edu.scripps.pms.census.conf.Configuration;
import edu.scripps.pms.census.exception.CensusGeneralException;
import edu.scripps.pms.census.exception.CensusIndexOutOfBoundException;
import edu.scripps.pms.census.hash.IndexedFile;
import edu.scripps.pms.census.hash.MSIndexFileCreator;
import edu.scripps.pms.census.util.CalcUtilGeneric;
import edu.scripps.pms.census.util.RelExFileFilter;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.List;

/**
 *
 * @author deep
 */
public class DIAAnalysisGetIndexFromMS1 {
    
     static BufferedReader brMS2=null;
    static BufferedWriter bw=null;
    
    public static void main(String[] args) throws CensusGeneralException, IOException, CensusIndexOutOfBoundException{
        
        if(args.length<1) {
            //usage error
            System.err.append("Not enough Argument");
            return;
        }
        
      /*  String ms2="/data/2/rpark/ip2_data/carolfc/DIA_Optimization2_July2016/DIA_400_1200_4mz_1e5_120ms_Rep1_2016_08_08_16_83980/spectra/data.ms2";
        String path=ms2.substring(0, ms2.lastIndexOf(File.separator));
        String output=path.substring(0, path.lastIndexOf(File.separator));*/
             
        String inputMS2=args[0];
        
         String path = inputMS2.substring(0, inputMS2.lastIndexOf(File.separator));
        
        String filename = inputMS2.substring(inputMS2.lastIndexOf(File.separator)+1);
        
         String filenameWOExtension = filename.substring(filename.lastIndexOf("."));
         String filenameFromMS = filename.substring(0,filename.indexOf(filenameWOExtension));
        String newMS2=path+File.separator+filenameFromMS+"_mp.ms2";
        
        String inputMS1=path+File.separator+filenameFromMS+".ms1";
        
         
         String pathMS1 = inputMS1.substring(0, inputMS1.lastIndexOf(File.separator))+File.separator;
        
         String filenameMS1 = inputMS1.substring(inputMS1.lastIndexOf(File.separator)+1);
       
     //   String inputMS2="/home/rampuria/data/DIA/DIAwithoutPrecur/xCT_400_1200_4mz_1e5_120ms_rep4.ms2";
      //  String newMS2="/home/rampuria/data/DIA/DIAwithoutPrecur/xCT_400_1200_4mz_1e5_120ms_rep4_new.ms2";
        
        brMS2=new BufferedReader(new FileReader(inputMS2));
        Hashtable<String, IndexedFile> ht=ChroGenerator.createIndexedFiles(pathMS1, filenameMS1);
      
         bw  = new BufferedWriter(new FileWriter(newMS2,false));
        
        DIAAnalysisGetIndexFromMS1.getMS1Scan(ht,pathMS1+filenameMS1);
         // System.out.println(ht);
       
          
    }
    
    public static void getMS1Scan(Hashtable<String, IndexedFile> ht, String filename) throws IOException, CensusIndexOutOfBoundException {
        
              
              IndexedFile file=ht.get(filename);
              
            
              String eachLine;
              int scan = 0,keyIndex = 0,prevscan=0;
              int[] keys = null;
              int count;
              double mass=0;
              List<String> zvalues=null;
              List<String> massList=null;
             
              while((eachLine=brMS2.readLine()) != null){
                  
                  
                if(eachLine.startsWith("S\t")){
                    
                  /*  if(zvalues!=null && massList!=null && scan!=0){
                        writeIntoFile(massList, zvalues, scan);
                    }*/
                    String[] words=eachLine.split("\t");
                    scan=Integer.parseInt(words[1]);
                    keys = file.getKeys();
                //    massList=new ArrayList<String>();
                    mass=getMS2ByMass(Double.parseDouble(words[3]));
              
                keyIndex = Arrays.binarySearch(keys, scan);

                if (keyIndex < 0) //Cannot find index
                {
                    keyIndex = -(++keyIndex); //Math.abs(++keyIndex);
                }
                if (keyIndex >= keys.length) {
                    keyIndex--;
                }
                
                  if(keyIndex>0) keyIndex--;
                  
                  prevscan=keys[keyIndex];
                  
                  bw.write(eachLine+"\n");

                }
              
               /*else if(eachLine.startsWith("I\t") || eachLine.startsWith("H\t")){
                    bw.write(eachLine+"\n");
                } */
               // System.out.println(keyIndex + "\t" + keys[keyIndex]);
               
                else if(eachLine.startsWith("Z\t")){
                
                //  zvalues=new ArrayList<String>();  
                    
                 String str=CalcUtilGeneric.getSpectrumString(file, keys[keyIndex]);
               //  System.out.println(str);
                 String[] data=str.split("\n");
              //   count=0;
                 for(String s:data){
                    // System.out.println(s);
                    if(!s.startsWith("S") && !s.startsWith("I")){
                     String[] massData=s.split("\t");
                  //   System.out.println(massData[0].split(" ")[0]);
                     
                     double massVal=Double.parseDouble(massData[0].split(" ")[0]);
                     if(Math.abs(mass-massVal)<=2){
                         int precursor=Integer.parseInt(massData[0].split(" ")[2]);
                      //   count++;
                         if(precursor==0){
                          //   zvalues.add("Z\t"+massVal+"\t"+2+"\n"+"Z\t"+massVal+"\t"+3+"\n");
                             bw.write("Z\t"+2+"\t"+massVal+"\n");
                             bw.write("Z\t"+3+"\t"+massVal+"\n");
                         }else{
                            // zvalues.add("Z\t"+massVal+"\t"+precursor+"\n");
                             bw.write("Z\t"+precursor+"\t"+massVal+"\n");
                         }
                     }
                 }
                  }
                 
               }
                
                else{
                   // System.out.println(eachLine+"\n");
                  //  massList.add(eachLine+"\n");
                    /*for(int i=0;i<zvalues.size();i++){
                        bw.write("S\t"+scan+"00001");
                        bw.write(zvalues.get(i));
                        for(int j=0;j<massList.size();j++){
                            bw.write(massList.get(j));
                        }
                    }*/
                    bw.write(eachLine+"\n");
                }
                 
                 
                 
              }
             
          
    }
    
  /*  public static void writeIntoFile(List<String> massList,List<String> zvalues,int scan) throws IOException{
        for(int i=0;i<zvalues.size();i++){
                        bw.write("S\t"+scan+"\t0000"+(i+1)+"\n");
                      //  System.out.println(zvalues.get(i));
                        bw.write(zvalues.get(i));
                        for(int j=0;j<massList.size();j++){
                            bw.write(massList.get(j));
                        }
                    }
    }*/
    
    public static int getMS2ByMass(double mass){
        int ms=(int)mass;
        int startValue=402;
        if(ms<=startValue){
            return startValue;
        }
        
        double diff=(ms-startValue)/4.0;
        
        if(diff-(int)diff >0.5){
            return (int)Math.ceil(diff)*4+startValue;
        }
        return (int)diff*4+startValue;
        
    }
    
    
      
    
}
