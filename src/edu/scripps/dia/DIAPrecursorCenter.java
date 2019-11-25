/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.scripps.dia;

import static edu.scripps.dia.DIAAnalysisGetIndexFromMS1.brMS2;
import static edu.scripps.dia.DIAAnalysisGetIndexFromMS1.bw;
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
import java.util.Set;
import java.util.HashSet;

/**
 *
 * @author Sanghvi
 */
public class DIAPrecursorCenter {

    static double PROTON_MASS=1.00728;

    public static void main(String[] args) throws IOException, CensusGeneralException, CensusIndexOutOfBoundException{

        if(args.length<1) {
            //usage error
            System.err.append("Not enough Argument");
            return;
        }
        String inputMS2=args[0];

        DIAPrecursorCenter.readFromOriginalMS2FileandCreateNemMS(inputMS2);

    }



    public static void readFromOriginalMS2FileandCreateNemMS(String inputMS2) throws IOException, CensusGeneralException, CensusIndexOutOfBoundException{

        BufferedReader brMS2=null;

        brMS2=new BufferedReader(new FileReader(inputMS2));

        String path = inputMS2.substring(0, inputMS2.lastIndexOf(File.separator));

        String filename = inputMS2.substring(inputMS2.lastIndexOf(File.separator)+1);

        String filenameWOExtension = filename.substring(filename.lastIndexOf("."));
        String filenameFromMS = filename.substring(0,filename.indexOf(filenameWOExtension));
        String newMS2=path+File.separator+filenameFromMS+"_mp.ms2";


        BufferedWriter bw=null;

        bw  = new BufferedWriter(new FileWriter(newMS2,false));

        String eachLine;
        int scan = 0,prevscan=0,count=0;
        int prevScanfromMS1=0,countH=0;

        double ms=0;
        Set<String> sLine=new HashSet<String>();



        while((eachLine=brMS2.readLine()) != null){


            if(eachLine!=null && eachLine.startsWith("S\t")){

                String[] words=eachLine.split("\t");
                scan=Integer.parseInt(words[1]);



                if(prevscan==scan){
                    while(eachLine!=null && (eachLine=brMS2.readLine()) != null){
                        if(eachLine.startsWith("S\t")){
                            String[] words1=eachLine.split("\t");
                            scan=Integer.parseInt(words1[1]);
                            if(scan!=prevscan){
                                break;
                            }
                        }
                    }
                }

                if(eachLine!=null && !sLine.contains(eachLine)){
                    bw.write(eachLine+"\n");
                    sLine.add(eachLine);
                }

                prevscan=scan;


            }

           else if(eachLine!=null && eachLine.startsWith("H\t")){
                bw.write(eachLine+"\n");
            }




            else if(eachLine!=null && eachLine.startsWith("I\t") && !eachLine.startsWith("I\tFilter")){

                bw.write(eachLine+"\n");
            }

            else if(eachLine!=null && eachLine.startsWith("I\tFilter")){
                String[] mass=eachLine.split("ms2");
                ms=Double.parseDouble(mass[1].split("@")[0]);
                bw.write(eachLine+"\n");
                continue;
            }


            else if(eachLine!=null && eachLine.startsWith("Z\t")){

                double precursormass1=getFinalMass(ms, 2);
                double precursormass2=getFinalMass(ms, 3);

                bw.write("Z\t"+2+"\t"+precursormass1+"\n");
                bw.write("Z\t"+3+"\t"+precursormass2+"\n");
            }

            else if(eachLine!=null){
                bw.write(eachLine+"\n");
            }

        }
        brMS2.close();
        bw.close();

    }

    public static double getFinalMass(double precursorMass,int cs){
        return (precursorMass*cs)-((cs-1)*PROTON_MASS);
    }
}
