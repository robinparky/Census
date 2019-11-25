/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author Rohan
 */
package edu.scripps.pms.census.labelFree;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.paukov.combinatorics.*;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.inference.TestUtils;

public class LabelfreeReportGenerator {
    
    
    public static void main (String [] args ) {
        
        
        
        ArrayList<Hashtable<String, Double>> info = new ArrayList<Hashtable<String, Double>>(); 
        JsonReader p = new JsonReader();
        ArrayList q = p.getFilesBySuffix(args[0], "JSON");
         
        File file = new File(args[1]);
        try{
                    FileWriter wf = new FileWriter(file.getAbsoluteFile(),true);
                    BufferedWriter wb = new BufferedWriter(wf);
                    wb.append("H\tACCESSION\tp-value\tPEPTIDE_NUM\tDESC");
                    wb.close();
        }catch (IOException e) {
			e.printStackTrace();
		}  
        Iterator t = q.iterator();
        while (t.hasNext())
        {
        JSONObject pep = (JSONObject)p.parser(args[0],t.next());
       JSONObject proteinList = (JSONObject)pep.get("protein");
       JSONArray ob = (JSONArray) pep.get("peptideList");
            System.out.println(ob.size());
        ArrayList<Double> pep2= new ArrayList<>();
        Iterator j = ob.iterator();
        while (j.hasNext())
        {
            
            Hashtable<String,Double> myHash = new Hashtable<>();
            JSONArray in = (JSONArray) j.next();
           
            Iterator k = in.iterator();
            while(k.hasNext())
            {
                
                JSONObject inn = (JSONObject) k.next();
                pep2.add((Double) inn.get("intensity"));
                
            }
       		ICombinatoricsVector<Integer> initialVector = Factory.createVector(new Integer[] { 2,3});
		Generator<Integer> generator = Factory.createPermutationGenerator(initialVector);
                Iterator<ICombinatoricsVector<Integer>> itr = generator.createIterator();
                int r=0;
                int o=1;
                while (itr.hasNext()) {
                    ICombinatoricsVector<Integer> permutation = itr.next();
                    List<Integer> myvector= permutation.getVector();
                    Iterator vItr = myvector.iterator(); 
                    while(vItr.hasNext())
                    {   
                        Object obj = vItr.next();
                        double ratio;
                        ratio= pep2.get(r)/pep2.get((int)obj);
                        myHash.put(""+o,ratio);     
                        r++;
                        o++;
                    }            
                    r =0;          
        } 
        o=1;
        info.add(myHash);
        pep2.clear();
          
    }
       int l=1;
       ArrayList peptideList = new ArrayList();
       ArrayList medianList = new ArrayList();
       for(l=1;l<=4;l++)
        {
            
        Iterator h = info.iterator();
        while(h.hasNext())
       {
           
           Hashtable obj = (Hashtable) h.next();
           peptideList.add(obj.get(""+l));
           
       }
                 DescriptiveStatistics stats = new DescriptiveStatistics();

                // Add the data from the array
                for( int z = 0; z < peptideList.size(); z++) {
                        stats.addValue((Double) peptideList.get(z));
               }

                // Compute some statistics
                double median = stats.getPercentile(50);
                medianList.add(median);
        }    
          ArrayList pValues = new ArrayList();
          double[] arrDoubleRatios = new double[medianList.size()];
           for(int b=0; b<medianList.size();b+=2)
           {
                for(int f=b; f<=(b+1);f++){    
                       
                    arrDoubleRatios[f] = (double) medianList.get(f);       
                }
               
                double tmp2;
                tmp2 = TestUtils.tTest(0.0,arrDoubleRatios);
                pValues.add(tmp2);
           }
           
           DescriptiveStatistics stats2 = new DescriptiveStatistics();

                // Add the data from the array
                for( int a = 0; a < pValues.size(); a++) {
                        stats2.addValue((Double) pValues.get(a));
               }

                // Compute some statistics
                double median2 = stats2.getPercentile(50);
                
                
                // Creating text file
                
                 try{  
                    FileWriter fw = new FileWriter(file.getAbsoluteFile(),true);
                    BufferedWriter bw = new BufferedWriter(fw);
			// if file doesnt exists, then create it
			if (!file.exists()) {
				file.createNewFile();
			}
                        bw.append("\nP\t"+proteinList.get("accession")+"\t"+median2+"\t"+ob.size()+"\t"+proteinList.get("desc"));
			bw.close();
                
                    }catch (IOException e) {
			e.printStackTrace();
                    }    
        }    
   }
}
