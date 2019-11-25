package edu.scripps.dia;

import edu.scripps.pms.census.ChroGenerator;
import edu.scripps.pms.census.exception.CensusGeneralException;
import edu.scripps.pms.census.exception.CensusIndexOutOfBoundException;
import edu.scripps.pms.census.hash.IndexedFile;
import edu.scripps.pms.census.util.CalcUtilGeneric;
import edu.scripps.pms.census.util.Dom4jUtil;

import javax.management.Query;
import java.io.*;
import java.util.HashMap;
import java.util.Hashtable;

/**
 * Created by Sanghvi on 1/3/17.
 */
public class DIAGetMS2FromClusteredMS1 {


    public static void main(String[] args) throws IOException, CensusGeneralException, CensusIndexOutOfBoundException {
        if(args.length<1) {
            //usage error
            System.err.append("Not enough Argument");
            return;
        }
        String inputMS2=args[0];

        String clusteredMS1File=args[1];

        int tolerance=200;

        DIAGetMS2FromClusteredMS1.getMS2Peaks(inputMS2,clusteredMS1File,tolerance);


    }

    public static void getMS2Peaks(String inputMS2,String clusteredMS1File,int tolerance) throws IOException, CensusIndexOutOfBoundException, CensusGeneralException {


        BufferedReader br=new BufferedReader(new FileReader(clusteredMS1File));

        BufferedReader brMS2=new BufferedReader(new FileReader(inputMS2));



        String pathMS2 = inputMS2.substring(0, inputMS2.lastIndexOf(File.separator))+File.separator;

        String filenameMS2 = inputMS2.substring(inputMS2.lastIndexOf(File.separator)+1);

        String output=pathMS2+"clusterfromMS2.txt";

        BufferedWriter bw=null;

        bw  = new BufferedWriter(new FileWriter(output,false));

        Hashtable<String, IndexedFile> ht= ChroGenerator.createIndexedFiles(pathMS2, filenameMS2);

        IndexedFile file=ht.get(pathMS2+filenameMS2);

        StringBuilder sb=null;
        String eachLine;

        double mass=0.0;
        int massFromMS2=0,numberofPeak=0;
        boolean flag=false,flagForMass=false,counter3=false;
        HashMap<String,Integer> masses=null;
        HashMap<String,StringBuilder> multipleMass=null;
        StringBuilder mulMass=null;
        int scanNum,begin=0,end=0,count=0,scan=0;



        while((eachLine=brMS2.readLine())!=null){
            if (eachLine.startsWith("H\t")){
                bw.write(eachLine+"\n");
            }else{
                break;
            }
        }

        while((eachLine=br.readLine()) != null) {

            flag = false;
            if (eachLine.startsWith("keyMass")) {
                massFromMS2 = getMS2ByMass(Double.parseDouble(eachLine.split(":")[1]));
                flagForMass=false;

                if(masses==null){
                    masses=new HashMap<>();
                    multipleMass=new HashMap<>();
                    sb=new StringBuilder();

                }else{

                    for(String k:masses.keySet()){

                        if(masses.get(k)==numberofPeak && masses.get(k)>=3){
                            if(counter3==false){
                                bw.write(sb.toString());
                                counter3=true;
                            }
                            bw.write(multipleMass.get(k).toString());
                        }
                    }
                    sb=new StringBuilder();
                    counter3=false;
                    masses.clear();
                    numberofPeak=0;
                }
               // bw.write("\n");
            } else if (eachLine.startsWith("Mass\t") || eachLine.startsWith("cluster")) {
                continue;
            } else if (eachLine.isEmpty()) {
                continue;
            } else {
                String[] massData = eachLine.split("\t");
                scanNum = Integer.parseInt(massData[2]);
                begin = scanNum + 1;
                end = scanNum + 101;
                numberofPeak++;

                for (int i = begin; i < end; i++) {

                    String str = CalcUtilGeneric.getSpectrumString(file, i);

                    String[] data = str.split("\n");
                    for (String s : data) {

                        if (s.startsWith("S")) {

                                mass = Double.parseDouble(s.split("\t")[3]);
                                scan = Integer.parseInt(s.split("\t")[1]);



                            if(!flag){

                            }

                            if ((int) mass != massFromMS2) {
                                break;
                            } else {
                                if(flagForMass==false) {
                                    sb.append("S\t0" + scan +"\t0"+scan+"\t"+mass+"\n");
                                    //bw.write("S\t0" + scan +"\t0"+scan+"\t"+mass+"\n");
                                }
                                flag = true;
                            }
                        }

                        else if(s.startsWith("I\t") && flagForMass==false){
                            sb.append(s+"\n");
                           // bw.write(s+"\n");

                        }

                        else if(s.startsWith("Z\t") && flagForMass==false){
                            sb.append(s+"\n");
                           // bw.write(s+"\n");

                        }

                        else if(flag && !s.startsWith("I\t") && !s.startsWith("Z\t")){
                            flagForMass=true;
                            String[] massIntensity=s.split("\t");
                            double ms=Double.parseDouble(massIntensity[0].split(" ")[0]);
                            String key=getKeyOfPeak(masses, ms,tolerance);
                            if(key == null){
                                masses.put(s, 1);
                                mulMass=new StringBuilder();
                                mulMass.append(s+"\n");
                                multipleMass.put(s,mulMass);
                            }else{
                                masses.put(key,masses.get(key)+1);
                                multipleMass.put(key,multipleMass.get(key).append(s+"\n"));
                            }

                           // bw.write(s+"\n");
                        }


                    }
                    if (flag) {
                        break;
                    }
                }
            }
        }
        br.close();
        bw.close();

    }

    public static String getKeyOfPeak(HashMap<String,Integer> map,double key,int tolerancePPM) throws IOException{
        if(map.isEmpty()){
            return null;
        }
        double tolerance=(key*tolerancePPM)/1000000;

        for(String str:map.keySet()){
            String[] massIntensity=str.split("\t");
            double k=Double.parseDouble(massIntensity[0].split(" ")[0]);

            if(Math.abs(k-key)<=tolerance) {
                return str;
            }
        }

        return null;
    }

    public static int getMS2ByMass(double mass){
        int ms=(int)mass;
        int startValue=502;
        if(ms<=startValue){
            return 0;
        }

        double diff=(ms-startValue)/4.0;

        if(diff-(int)diff >0.5){
            return (int)Math.ceil(diff)*4+startValue;
        }
        return (int)diff*4+startValue;

    }

}
