/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package com.TmtFilter.DAN;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Harshil
 * 
 * This class will map the Key(Peptide )  to the value(List of Proteins that contains the Peptide.).
 * Its for DAN's Prooject...
 */
public class multihash
{
    Map<String,ArrayList<String>> multiMap = new HashMap();
    String peptide_list[]=new String[100000];
    int count=0;

    public void add_data(String key,String value)
    {
        ArrayList<String> ar= new ArrayList();
        if(multiMap.containsKey(key))
        {
            ar=multiMap.get(key);
            if(!ar.contains(value))
            {
                ar.add(value);
                multiMap.put(key, ar);                    
                
            }
        }
        else
        {
            ar.add(value);
            multiMap.put(key, ar);   
            peptide_list[count]=key;
            count++;

        }
        
    }
    public void display(Map<String,pep_details> peptide_data,ArrayList<String> tests,String directory)
    {
        
        String key=null;
        BufferedWriter bufferedWriter = null;
        pep_details data=new pep_details();
        try
        {
            //              c:/work/output_peptideAnalysis.txt
            if(directory == "")
                bufferedWriter = new BufferedWriter(new FileWriter("C:\\Users\\Harshil\\Documents\\NetBeansProjects\\TmtFilter\\input_data\\output_peptideAnalysis.txt"));
            else
                bufferedWriter = new BufferedWriter(new FileWriter(directory+ File.separator+"output_peptideAnalysis.txt"));//            System.out.println("Peptide \t Protine");
            bufferedWriter.write("H"+"\tPeptide_Name\t");
               for(int k=0;k<tests.size();k++)
                {
                      bufferedWriter.write("Avg_" + tests.get(k)+"\t");
                      bufferedWriter.write("STD_" + tests.get(k) + "\t");
                }
                bufferedWriter.write("SpC"+"\t"+"ScanNum"+"\t"+"Cstate"+"\t"+"Filename"+"\t"+"Protein_Set"+"\t");
                bufferedWriter.write("Ratio_avg1\tRatio_avg2\t"+"Ratio_avg3\tRatio_stdev1\t"+"Ratio_stdev2\tRatio_stdev3\n");
            for(int i=0;i<count;i++)
            {   
                key=peptide_list[i];
                data=peptide_data.get(key);
                
                bufferedWriter.write("S\t"+ peptide_list[i]);
                for(int k=0;k<tests.size();k++)
                {
                      bufferedWriter.write("\t" + data.getIntensity(tests.get(k)) );
                      bufferedWriter.write("\t" + data.getstat(tests.get(k)).getStandardDeviation());
                      
                    //  System.out.println(key +"\t" + tests.get(k)+ "\t" + data.getstat(tests.get(k)));
                }
//                data.printratio();
//                System.out.print("==>" + key);
  //              for(int k=0;k<6;k++)
                       // System.out.println("\t" + data.getallIntensity("norm_m/z_127.1248_int"));

 //               System.out.println("The average :"+data.getratio_avg()[0]+" "+ data.getratio_avg()[1]+" "+data.getratio_avg()[2]);
  //              System.out.println("The stdev: "+data.getratio_stdev()[0]+" "+ data.getratio_stdev()[1]+" "+data.getratio_stdev()[2]);


                bufferedWriter.write("\t"+data.getSpc() +"\t"+ data.getScannum()+"\t"+data.getCstate()+"\t"+data.getFilename());
                bufferedWriter.write("\t"+ multiMap.get(peptide_list[i]));
                bufferedWriter.write("\t"+data.getratio_avg()[0]+"\t"+ data.getratio_avg()[1]+"\t"+data.getratio_avg()[2]);
                bufferedWriter.write("\t" + data.getratio_stdev()[0]+"\t"+ data.getratio_stdev()[1]+"\t"+data.getratio_stdev()[2]);
                
                bufferedWriter.newLine();
            }
            
            //System.out.println("S\t"+ peptide_list[i]+"\t" + data.getIntensity(tests.get(k))+"\t" + data.getstat(tests.get(k)).getStandardDeviation()+"\t"+data.getSpc() +"\t"+ data.getScannum()+"\t"+data.getCstate()+"\t"+data.getFilename());
		System.out.println("Success"+ directory +File.separator +"output_peptideAnalysis.txt" );
            bufferedWriter.close();
        } 
        catch (Exception ex) 
        {
            
            Logger.getLogger(multihash.class.getName()).log(Level.SEVERE, null, ex);
        } 
        
        
        
        
    }
//    public static void main(String args[])
//    {
//        Map<String,ArrayList<String>> multiMap;
//        multiMap = new HashMap<>();
//        String key;
//        ArrayList<String> ar= new ArrayList();
//        ArrayList<String> ar1= new ArrayList();
//        ar.add("asdsa");
//        ar.add("qwewq");
//       // ar=multiMap.get("p1");
//        
//        multiMap.put("p1",ar);
//        ar.add("qwewq");
//        multiMap.put("p2",ar);
//        ar1=multiMap.get("p2");
//        ar1.add("qwewq2222");
//        ar1.add("qwewq55656");
//        multiMap.put("p2",ar1);
//
//        System.out.println(multiMap);
//        
//    }
//        
    }

