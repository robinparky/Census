       
/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.scripps.dia;

//import com.cluster.DBSCAN;
import java.io.File;
import edu.scripps.pms.census.CensusConstants;
import static edu.scripps.pms.census.ChroGenerator.createIndexedFiles;
import edu.scripps.pms.census.conf.Configuration;
import edu.scripps.pms.census.exception.CensusIndexOutOfBoundException;
import edu.scripps.pms.census.hash.IndexedFile;
//import edu.scripps.pms.census.labelFree.LabelfreeMissingPeptideBuilder;
import edu.scripps.pms.census.labelFree.util.LabelfreeChroUtil;
import edu.scripps.pms.census.model.SpecRange;
import edu.scripps.pms.census.tools.Formatter;
import edu.scripps.pms.census.util.CalcUtil;
import static edu.scripps.pms.census.util.CalcUtil.parseSpectra;
import static edu.scripps.pms.census.util.CalcUtil.readFullSpectrum;
import edu.scripps.pms.census.util.dtaselect.Peptide;
import edu.scripps.pms.census.util.io.MzxmlSpectrumReader;
import gnu.trove.TIntDoubleHashMap;
import gnu.trove.TIntLongHashMap;
import java.io.*;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStream;
import java.io.RandomAccessFile;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.apache.commons.lang3.Range;
import org.apache.commons.math3.analysis.function.Gaussian;
import org.apache.commons.math3.ml.clustering.DoublePoint;
import rpark.statistics.GaussianFitting;
import rpark.statistics.Smooth;
import rpark.statistics.model.GaussianPeakModel;


/**
 *
 * @author Sanghvi
 */

public class DIAAnalysis implements Serializable{
    

    
    public class MS1Peak implements Serializable { 
        private static final long serialVersionUID = 1L;
        double mass;
        List<Double> massLst;
        List<Double> retTime;
        List<Integer> scanNum;
        List<Integer> intensity;
        List<Integer> chargeState;

        public MS1Peak(double mass) {
           // super(new double[]{mass});
            this.mass=mass;
            massLst=new ArrayList<Double>();
            retTime=new ArrayList<Double>();
            scanNum=new ArrayList<Integer>();
            intensity=new ArrayList<Integer>();
            chargeState=new ArrayList<Integer>();
        }
        
        public double getMass(){
            return this.mass; 
        }
        
        public void addMass(double mass){
            this.massLst.add(mass);
        }
        
        public void addretTime(double retTime){
            this.retTime.add(retTime);
        }
        
        public void addscanNum(int scanNum){
            this.scanNum.add(scanNum);
        }
        
        public void addIntensity(int intensity){
            this.intensity.add(intensity);
        }

        public void addChargeState(int cs){
            this.chargeState.add(cs);
        }

        
    }
    
    public static void main(String[] args) throws IOException, FileNotFoundException, ClassNotFoundException{


        if(args.length<1) {
            //usage error
            System.err.append("Not enough Argument");
            return;
        }
        String inputMS1=args[0];

       // String input="/home/rampuria/data/DIA/DIAwithoutPrecur/DIA_400_1200_8mz_3e6_120ms_Rep3_2016_08_08_16_83999/carolina/a1.ms1";
        int tolerance=200;


        DIAAnalysis da=new DIAAnalysis();
        da.createMS1peak(inputMS1,tolerance);
        
    }
    
    public void createMS1peak(String filePath,int tolerance) throws FileNotFoundException, IOException, ClassNotFoundException{
           
            BufferedReader br = new BufferedReader(new FileReader(filePath));
            String eachLine;
            int scan=0;
            double rettime=0.0;
            boolean flag=false;
            
            HashMap<Double,MS1Peak> peak=new HashMap<>();
            
           // List<Double> mass=new ArrayList<Double>();
            
            
            
            while((eachLine=br.readLine()) != null){
                MS1Peak ms1Peak;
                if(eachLine.startsWith("H\t")){
                    continue;
                }
                
                String[] words=eachLine.split("\t");
                
                if(eachLine.startsWith("S\t")){
                    scan=Integer.parseInt(words[1]);
                    //System.out.println("scanNum:"+scan);
                    flag=true;
                }
                
                else if(eachLine.startsWith("I\t") && flag==true){
                    rettime=Double.parseDouble(words[2]);
                   // System.out.println("rettime:"+rettime);
                    flag=false;
                }
                
                else if(eachLine.startsWith("I\t") && flag==false){
                    continue;
                }
                
                 else if(eachLine.startsWith("Z\t")){
                    continue;
                }
                
                
                else{                
                        String[] values=words[0].split(" ");
                        double ms=Double.parseDouble(values[0]);
                        int intensity=(int)Double.parseDouble(values[1]);
                        int cs=Integer.parseInt(values[2]);
                        double key=getKeyOfPeak(peak, ms,scan,tolerance);
                        if(key == -1.0){

                            ms1Peak=new MS1Peak(ms);
                            peak.put(ms, ms1Peak);
                            ms1Peak.addMass(ms);
                            ms1Peak.addIntensity(intensity);
                            ms1Peak.addretTime(rettime);
                            ms1Peak.addscanNum(scan);
                            ms1Peak.addChargeState(cs);
                        }else{

                            ms1Peak=peak.get(key);
                            ms1Peak.addMass(ms);
                            ms1Peak.addIntensity(intensity);
                            ms1Peak.addretTime(rettime);
                            ms1Peak.addscanNum(scan);
                            ms1Peak.addChargeState(cs);
                        }

           
                }
            }
        //printobject(peak);
        writeObjectInFile(peak,filePath);
        //readObjectFromFile();
    }
    
    public static double getKeyOfPeak(HashMap<Double,MS1Peak> map,double key,int scan,int tolerancePPM) throws IOException{
        if(map.isEmpty()){
            return -1.0;
        }
        double tolerance=(key*tolerancePPM)/1000000;
        
        for(double k:map.keySet()){
                if(Math.abs(k-key)<=tolerance) {
                    MS1Peak pk=map.get(k);
                    if(scan-pk.scanNum.get(pk.scanNum.size()-1)>101){
                        return -1;
                    }
                    return k;
                }        
        }
       // writeObjectInFile(map, key);
        //map.remove(key);
        return -1.0;
    }
    
    
    /* public static void writeObjectInFile(HashMap<Double,MS1Peak> map,double k) throws FileNotFoundException, IOException{
        
         //BufferedWriter bw  = new BufferedWriter(new FileWriter("/home/rampuria/data/DIA/ex3.txt",false));
        
        List<MS1Peak> peakList=new ArrayList<>();
       // try (FileOutputStream fout = new FileOutputStream("/home/rampuria/test_data/test_data2/peak.ser", false)) {
         //   ObjectOutputStream oos = new ObjectOutputStream(fout);
           // for(double k:map.keySet()){
                MS1Peak pk=map.get(k);
                // bw.write(pk.getMass()+"\t"+pk.massLst+"\t"+pk.scanNum+"\t"+pk.retTime+"\n");
                
                 int size=pk.massLst.size();
                 if(size>=3){
                 bw.write(pk.getMass()+"\n");
                 for(int i=0;i<size;i++){
                     bw.write(pk.massLst.get(i)+"\t"+pk.intensity.get(i)+"\t"+pk.scanNum.get(i)+"\t"+pk.retTime.get(i)+"\n");
                 }
                 bw.write("\n cluster \n");
                 }
              //  peakList.add(map.get(k));
           // }
           // oos.writeObject(peakList);
            //oos.close();
            //fout.close();
            bw.close();
       // }
    }*/
    
    public static void writeObjectInFile(HashMap<Double,MS1Peak> map,String inputMS1) throws FileNotFoundException, IOException{
        

        BufferedWriter bw  = null;


        String path = inputMS1.substring(0, inputMS1.lastIndexOf(File.separator));

        String filename = inputMS1.substring(inputMS1.lastIndexOf(File.separator)+1);

        String filenameWOExtension = filename.substring(filename.lastIndexOf("."));
        String filenameFromMS = filename.substring(0,filename.indexOf(filenameWOExtension));
        String newFile=path+File.separator+filenameFromMS+"Cluster.txt";

       // bw=new BufferedWriter(new FileWriter("/home/rampuria/data/DIA/DIAwithoutPrecur/DIA_400_1200_8mz_3e6_120ms_Rep3_2016_08_08_16_83999/carolina/exMS_2.txt",false));
        bw=new BufferedWriter(new FileWriter(newFile,false));
        List<MS1Peak> peakList=new ArrayList<>();

            for(double k:map.keySet()){
                MS1Peak pk=map.get(k);
                // bw.write(pk.getMass()+"\t"+pk.massLst+"\t"+pk.scanNum+"\t"+pk.retTime+"\n");
                
                 int size=pk.massLst.size();
                 if(size>=3){
                     System.out.println(size);
                 bw.write("keyMass:"+pk.getMass()+"\n");
                  bw.write("Mass\t"+"SumOfIntensity\t"+"scanNum\t"+"retTime\t\t"+"intensity+chargeState"+"\n");
                     int i=0;
                 while(i<size){
                     int j=i+1;
                     int intensity=pk.intensity.get(i);
                     while(j<size){
                         if(pk.scanNum.get(i).compareTo(pk.scanNum.get(j))==0){

                             intensity=intensity+pk.intensity.get(j);
                             j++;
                         }else{
                             break;
                         }
                     }

                     bw.write(pk.massLst.get(i)+"\t"+intensity+"\t"+pk.scanNum.get(i)+"\t"+pk.retTime.get(i)+"\t");
                     StringBuilder sb=new StringBuilder();
                     for(int l=i;l<j;l++){
                         sb.append(pk.intensity.get(l)+"\t"+pk.chargeState.get(l)+",");
                     }
                     sb.deleteCharAt(sb.length()-1);
                     bw.write(sb.toString());
                     bw.write("\n");
                     i=j;
                 }
                 bw.write("cluster \n");
                 bw.write("\n");
                 }
               // peakList.add(map.get(k));
            }

            bw.close();

    }
    
    public static void readObjectFromFile() throws FileNotFoundException, IOException, ClassNotFoundException{
         HashMap<Double, MS1Peak> map = null;
         MS1Peak peak;
         FileInputStream fis = new FileInputStream("/home/rampuria/test_data/test_data2/peak.ser"); 
         ObjectInputStream ois = new ObjectInputStream(fis);
           // map = (HashMap) ois.readObject();
           // peak = (MS1Peak) ois.readObject();
           // System.out.println(ois.readObject().toString());
            List<MS1Peak> peakList=(List<MS1Peak>) ois.readObject();
            for(MS1Peak ms1peak:peakList){
                System.out.println("mass:"+ms1peak.getMass());
                System.out.println("masses:"+ms1peak.massLst);
                System.out.println("Intensity:"+ms1peak.intensity);
                System.out.println("scanNum:"+ms1peak.scanNum);
                System.out.println("rettime:"+ms1peak.retTime);
                System.out.println();
            }
            ois.close();
            fis.close();
         
        // Set set = map.entrySet();
         //Iterator iterator = set.iterator();
         
         /*while(iterator.hasNext()) {
            Map.Entry mentry = (Map.Entry)iterator.next();
            peak=(MS1Peak) mentry.getValue();
            System.out.println("mass:"+peak.mass);
            System.out.println("Intensity:"+peak.intensity);
            System.out.println("scanNum:"+peak.scanNum);
            System.out.println("rettime:"+peak.retTime);
            System.out.println();
        }*/
    }
    
    /*public static void printobject(HashMap<Double,MS1Peak> map){
        MS1Peak peak;
        for(double k:map.keySet()){
            peak=map.get(k);
            System.out.println("mass:"+peak.getMass());
            System.out.println("Intensity:"+peak.intensity);
            System.out.println("scanNum:"+peak.scanNum);
            System.out.println("rettime:"+peak.retTime);
            System.out.println();
        }
    }*/
}
