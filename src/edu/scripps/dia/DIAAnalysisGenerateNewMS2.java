/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.scripps.dia;

import edu.scripps.pms.census.ChroGenerator;
import edu.scripps.pms.census.exception.CensusGeneralException;
import edu.scripps.pms.census.exception.CensusIndexOutOfBoundException;
import edu.scripps.pms.census.hash.IndexedFile;
import edu.scripps.pms.census.util.CalcUtilGeneric;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.List;

/**
 *
 * @author Sanghvi
 */
public class DIAAnalysisGenerateNewMS2 {

    static double PROTON_MASS=1.00728;

    public static void main(String[] args) throws IOException, CensusGeneralException, CensusIndexOutOfBoundException{

        if(args.length<1) {
            //usage error
            System.err.append("Not enough Argument");
            return;
        }
        String inputMS2=args[0];

        DIAAnalysisGenerateNewMS2.readFromOriginalMS2FileandCreateNewMS(inputMS2);

    }


    public static int getMS1Scan(Hashtable<String, IndexedFile> ht, String filename,int scan) throws IOException, CensusIndexOutOfBoundException {

        int keyIndex=0;
        int[] keys = null;


        IndexedFile file=ht.get(filename);
        keys = file.getKeys();

        keyIndex = Arrays.binarySearch(keys, scan);

        if (keyIndex < 0) //Cannot find index
        {
            keyIndex = -(++keyIndex); //Math.abs(++keyIndex);
        }
        if (keyIndex >= keys.length) {
            keyIndex--;
        }

        if(keyIndex>0) keyIndex--;

        return keys[keyIndex];
    }

    public static double getFinalMass(double precursorMass,int cs){
        return (precursorMass*cs)-((cs-1)*PROTON_MASS);
    }

    public static void readFromOriginalMS2FileandCreateNewMS(String inputMS2) throws IOException, CensusGeneralException, CensusIndexOutOfBoundException{

        BufferedReader brMS2=null;

        brMS2=new BufferedReader(new FileReader(inputMS2));

        String path = inputMS2.substring(0, inputMS2.lastIndexOf(File.separator));

        String filename = inputMS2.substring(inputMS2.lastIndexOf(File.separator)+1);

        String filenameWOExtension = filename.substring(filename.lastIndexOf("."));
        String filenameFromMS = filename.substring(0,filename.indexOf(filenameWOExtension));
        String newMS2=path+File.separator+filenameFromMS+"_mp.ms2";

        String inputMS1=path+File.separator+filenameFromMS+".ms1";

        String pathMS1 = inputMS1.substring(0, inputMS1.lastIndexOf(File.separator))+File.separator;

        String filenameMS1 = inputMS1.substring(inputMS1.lastIndexOf(File.separator)+1);

        Hashtable<String, IndexedFile> ht=ChroGenerator.createIndexedFiles(pathMS1, filenameMS1);

        IndexedFile file=ht.get(pathMS1+filenameMS1);

        BufferedWriter bw=null;

        bw  = new BufferedWriter(new FileWriter(newMS2,false));

        String eachLine;
        int scan = 0,prevscan=0,count=0;
        int prevScanfromMS1=0,countH=0;

        double mass=0;
        List<String> zvalues=null;
        List<String> massList=null;
        StringBuilder sb=new StringBuilder();

        while((eachLine=brMS2.readLine()) != null){


            if(eachLine!=null && eachLine.startsWith("S\t")){
                countH=1;
                if(zvalues!=null && massList!=null && scan!=0 && count!=0){
                    writeIntoFile(massList, zvalues, scan,sb,mass,bw);
                    sb=null;
                }
                count=0;
                sb=new StringBuilder();
                String[] words=eachLine.split("\t");
                scan=Integer.parseInt(words[1]);
                if(prevscan==scan){
                    while((eachLine=brMS2.readLine()) != null){
                        if(eachLine.startsWith("S\t")){
                            String[] words1=eachLine.split("\t");
                            scan=Integer.parseInt(words1[1]);
                            if(scan!=prevscan){
                                break;
                            }
                        }
                    }
                }
                prevscan=scan;

                massList=new ArrayList<String>();
                mass=getMS2ByMass(Double.parseDouble(words[3]));



                sb.append("S\t0000"+scan+"\t0000"+scan+"\t"+Double.parseDouble(words[3])+"\n");
                //  bw.write(eachLine+"\n");

            }

            if(eachLine!=null && eachLine.startsWith("H\t")){
                bw.write(eachLine+"\n");
            }

            else if(eachLine!=null && eachLine.startsWith("I\t")){
                sb.append(eachLine+"\n");
                // bw.write(eachLine+"\n");
            }


            else if(eachLine!=null && eachLine.startsWith("Z\t")){

                prevScanfromMS1=getMS1Scan(ht,pathMS1+filenameMS1,scan);

                zvalues=new ArrayList<String>();

                String str=CalcUtilGeneric.getSpectrumString(file, prevScanfromMS1);

                String[] data=str.split("\n");

                for(String s:data){

                    if(!s.startsWith("S") && !s.startsWith("I")){
                        String[] massData=s.split("\t");


                        double massVal=Double.parseDouble(massData[0].split(" ")[0]);
                        if(Math.abs(mass-massVal)<=2){
                            int chargeState=Integer.parseInt(massData[0].split(" ")[2]);
                            count++;
                            if(chargeState==0){
                                double precursormass1=getFinalMass(massVal, 2);
                                double precursormass2=getFinalMass(massVal, 3);
                                zvalues.add("Z\t"+2+"\t"+precursormass1+"\n"+"Z\t"+3+"\t"+precursormass2+"\n");


                            }else{
                                double precursorFinalmass=getFinalMass(massVal, chargeState);
                                zvalues.add("Z\t"+chargeState+"\t"+precursorFinalmass+"\n");

                            }
                        }
                    }
                }

            }

            else if(eachLine!=null){

                massList.add(eachLine+"\n");

            }

        }
        brMS2.close();
        bw.close();

    }

    public static void writeIntoFile(List<String> massList,List<String> zvalues,int scan,StringBuilder sb,double mass, BufferedWriter bw) throws IOException{

        bw.write(sb.toString());

        for(int i=0;i<zvalues.size();i++){
            // sb.replace(sb.indexOf("\t")+1,sb.indexOf("\n")-1,"0000"+(scan*1000+i+1)+"\t"+"0000"+(scan*1000+i+1)+"\t"+mass);
            bw.write(zvalues.get(i));

            //bw.write(sb.toString());
        }

        for(int j=1;j<massList.size();j++){
            bw.write(massList.get(j));
        }

    }

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
