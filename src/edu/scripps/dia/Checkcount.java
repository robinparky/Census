/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.scripps.dia;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.*;

/**
 *
 * @author rampuria
 */
public class Checkcount {
    
    public static void main(String [] args){
        
        double [] genrated = new double[297782];
        double [] compare = new double[57947];
        try{
            BufferedReader br = new BufferedReader(new FileReader("/home/rampuria/data/carolfc/DIA_400_1200_4mz_1e5_120ms_Rep1_2016_08_08_16_83980/spectra/dd"));
            String eachLine = null;
            int count=0;
            
            while((eachLine=br.readLine()) != null){
                genrated[count] = Double.parseDouble(eachLine);
                count++;
            }
            br.close();
            br = new BufferedReader(new FileReader("/home/rampuria/data/carolfc/cc"));
             eachLine = null;
             count=0;
            
            while((eachLine=br.readLine()) != null){
                compare[count] = Double.parseDouble(eachLine);
                count++;
            }
            
            int trueCount=0;
            for(int i =0;i<genrated.length;i++){
                double ppm = (genrated[i]/1000)*12;
                int key = getClosestK(compare,genrated[i]);
                if(Math.abs(compare[key]-genrated[i]) <= ppm){
                    trueCount++;
                }
            }
            System.out.println("");
            
            
            
            
        }catch(Exception e){
            e.printStackTrace();
        }
        
        
    }
    public static int getClosestK(double[] a, double x) {

    int low = 0;
    int high = a.length - 1;

    if (high < 0)
        throw new IllegalArgumentException("The array cannot be empty");

    while (low < high) {
        int mid = (low + high) / 2;
        assert(mid < high);
        double d1 = Math.abs(a[mid  ] - x);
        double d2 = Math.abs(a[mid+1] - x);
        if (d2 <= d1)
        {
            low = mid+1;
        }
        else
        {
            high = mid;
        }
    }
    return high;
}
    
    
}
