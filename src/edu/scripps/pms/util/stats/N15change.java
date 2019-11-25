/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.scripps.pms.util.stats;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author rampuria
 */
public class N15change {
    
    
    
    public static void main (String [] args){
        
        String filename = "/home/rampuria/data/census-out_n15_5.txt";
        
        N15change n = new N15change();
        n.findchange(filename);
        
    }
    
    
    public void findchange (String filename){
        
        
        File file = new File(filename);
        List<Double> ratio = new ArrayList<>();
        
        
        try{
            BufferedReader br = new BufferedReader(new FileReader(file));
            String eachLine = null;
            while((eachLine=br.readLine()) != null){
                
                if(eachLine.startsWith("P\t")){                    
                    String [] words = eachLine.split("\t");
                    double compositeVal = Double.parseDouble(words[6]);
                    double fa = 100*(1/(1+compositeVal));
                    System.out.println(""+fa);
                    ratio.add(fa);
                }
            }
            Histogram h = new Histogram(20, 0, 100);
                h.loadData(ratio);
                String result = h.getResult();
                
            
            
        }catch(Exception e){
            e.printStackTrace();
        }
        
    }
    
}
