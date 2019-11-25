/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package com.TmtFilter.DAN;

import org.apache.commons.math3.stat.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
/**
 *
 * @author Harshil
 * 
 * This method contains the Each row of the Peptide sequence from the "census-out.txt" file.
 */
public class pep_details 
{
    private Map<String, Float> intensity = new HashMap();
    // The tests that are carried on the Peptide sequence. The key is the peptide sequence.
    private Map<String, ArrayList> all_intensity = new HashMap();// This one will store all the intensities.
    private Map<String, DescriptiveStatistics> stats=new HashMap();// The key is name of the test carried on the peptide. 
    private int spc=0;
    private ArrayList<Integer> scannum = new ArrayList();
    private ArrayList<Integer> cstate = new ArrayList();
    private Set<String> filename = new HashSet();
    private Map<String, Integer> peptide_counter=new HashMap(); // The key is the peptide name. 
    private Map<String, ArrayList> ratio=new HashMap();
    private float ratio_avg[]=new float[3];     // This one will store the average the average of the Ratios 1 2 and 3.
     private float ratio_stdev[]= new float[3]; // This one will store the standadr deviation of the ratios 1 2 3 .
     private DescriptiveStatistics stdev = new DescriptiveStatistics();

    public void setRatioAvg(List total)
    {
        ratio_avg[0] = (float)total.get(0);
        ratio_avg[1] = (float)total.get(1);
        ratio_avg[2] = (float)total.get(2);
    }
    public void setratio_avg()
    {
        ArrayList ar= new ArrayList();
        ar=ratio.get("0");
        float sum=0,number;
        Iterator it = ar.iterator();
        stdev=new DescriptiveStatistics();
        while(it.hasNext())
        {
            number=(float) it.next();
            stdev.addValue(number);
            sum+= number;
        }
        ratio_avg[0]= sum/ar.size();
        ratio_stdev[0]=(float) stdev.getStandardDeviation();
        
        sum=0;stdev=new DescriptiveStatistics();
        ar=ratio.get("1");
        it= ar.iterator();
        while(it.hasNext())
        {
           number=(float) it.next();
            stdev.addValue(number);
            sum+= number;
        }
        ratio_avg[1]= sum/ar.size();
        ratio_stdev[1]=(float) stdev.getStandardDeviation();
        sum=0;stdev=new DescriptiveStatistics();
        ar=ratio.get("2");
        it= ar.iterator();
        while(it.hasNext())
        {
            number=(float) it.next();
            stdev.addValue(number);
            sum+= number;
        }
        ratio_avg[2]= sum/ar.size();
        ratio_stdev[2]=(float) stdev.getStandardDeviation();
    }
    public void setratio_avg(float[] ratio_average)
    {
        ArrayList ar= new ArrayList();
        ar=ratio.get("0");
        float sum=0,number;
        Iterator it = ar.iterator();
        stdev=new DescriptiveStatistics();
        while(it.hasNext())
        {
            number=(float) it.next();
            stdev.addValue(number);
            sum+= number;
        }
        ratio_avg[0]= ratio_average[0];
        ratio_stdev[0]=(float) stdev.getStandardDeviation();
        
        sum=0;stdev=new DescriptiveStatistics();
        ar=ratio.get("1");
        it= ar.iterator();
        while(it.hasNext())
        {
           number=(float) it.next();
            stdev.addValue(number);
            sum+= number;
        }
        ratio_avg[1]= ratio_average[1];
        ratio_stdev[1]=(float) stdev.getStandardDeviation();
        sum=0;stdev=new DescriptiveStatistics();
        ar=ratio.get("2");
        it= ar.iterator();
        while(it.hasNext())
        {
            number=(float) it.next();
            stdev.addValue(number);
            sum+= number;
        }
        ratio_avg[2]= ratio_average[2];
        ratio_stdev[2]=(float) stdev.getStandardDeviation();
    }
    
     public void setratio_avg1(float[] ratio_average)
    {
        ArrayList ar= new ArrayList();
        ar=ratio.get("0");
        float sum=0,number;
        Iterator it = ar.iterator();
        stdev=new DescriptiveStatistics();
        while(it.hasNext())
        {
            number=(float) it.next();
            stdev.addValue(number);
            sum+= number;
        }
        ratio_avg[0]= ratio_average[0];
        ratio_stdev[0]=(float) stdev.getStandardDeviation();
        
        sum=0;stdev=new DescriptiveStatistics();
        ar=ratio.get("1");
        it= ar.iterator();
        while(it.hasNext())
        {
           number=(float) it.next();
            stdev.addValue(number);
            sum+= number;
        }
        ratio_avg[1]= ratio_average[1];
        ratio_stdev[1]=(float) stdev.getStandardDeviation();
        sum=0;stdev=new DescriptiveStatistics();
        ar=ratio.get("2");
        it= ar.iterator();
        while(it.hasNext())
        {
            number=(float) it.next();
            stdev.addValue(number);
            sum+= number;
        }
        ratio_avg[2]= ratio_average[2];
        ratio_stdev[2]=(float) stdev.getStandardDeviation();
    }
    public float[] getratio_avg()
    {
        return ratio_avg;
    }
    public float[] getratio_stdev()
    {
        return ratio_stdev;
    }
    
    public float getIntensity(String key) 
    {
         return intensity.get(key);
    }
    public Map<String,Float> getIntensityMap()
    {
           return intensity;
    }     
    public void setallIntensity(String key, float value)
    {
        ArrayList ar=new ArrayList();
        
        if(all_intensity.containsKey(key))
        {
            ar=all_intensity.get(key);
        }
  
        ar.add(value);
        
        all_intensity.put(key, ar);
    }
    
    public void setratio(List norm_list)
    {
        int currentrow=0;
        ArrayList ar=new ArrayList();
        float sum=0;
      //  System.out.println((all_intensity.get(norm_list.get(0))).get(0));
        int middle = norm_list.size()/2;
        for(int i=0;i<norm_list.size();i++)
        {
            ar= new ArrayList();
         //   sum=((float) all_intensity.get(norm_list.get(i)).get(0) + (float)all_intensity.get(norm_list.get(middle+i)).get(0)) / ( (float)all_intensity.get(norm_list.get(i+1)).get(0)+ (float)all_intensity.get(norm_list.get(4)).get(0));
        }
        
        if(!ratio.containsKey("0"))
        {
            ar=new ArrayList();
            sum=((float) all_intensity.get(norm_list.get(0)).get(0) + (float)all_intensity.get(norm_list.get(3)).get(0)) / ( (float)all_intensity.get(norm_list.get(1)).get(0)+ (float)all_intensity.get(norm_list.get(4)).get(0));
            ar.add(sum);
            ratio.put("0", ar);
            ar=new ArrayList();
            sum=((float) all_intensity.get(norm_list.get(0)).get(0) + (float) all_intensity.get(norm_list.get(3)).get(0)) / ((float)all_intensity.get(norm_list.get(2)).get(0)+ (float)all_intensity.get(norm_list.get(5)).get(0));
            ar.add(sum);
            ratio.put("1",ar);
            ar=new ArrayList();
            sum=((float)all_intensity.get(norm_list.get(1)).get(0) + (float) all_intensity.get(norm_list.get(4)).get(0)) / ((float) all_intensity.get(norm_list.get(2)).get(0)+ (float) all_intensity.get(norm_list.get(5)).get(0));
            ar.add(sum);
            ratio.put("2",ar);
            
        }        
        else
        {
            currentrow=ratio.get("0").size();
            ar=ratio.get("0");
            sum=((float) all_intensity.get(norm_list.get(0)).get(currentrow) + (float)all_intensity.get(norm_list.get(3)).get(currentrow)) / ( (float)all_intensity.get(norm_list.get(1)).get(currentrow)+ (float)all_intensity.get(norm_list.get(4)).get(currentrow));
            ar.add(sum);
            ratio.put("0", ar);
            ar=ratio.get("1");;
            sum=((float) all_intensity.get(norm_list.get(0)).get(currentrow) + (float) all_intensity.get(norm_list.get(3)).get(currentrow)) / ((float)all_intensity.get(norm_list.get(2)).get(currentrow)+ (float)all_intensity.get(norm_list.get(5)).get(currentrow));
            ar.add(sum);
            ratio.put("1",ar);
            ar=ratio.get("2");;
            sum=((float)all_intensity.get(norm_list.get(1)).get(currentrow) + (float) all_intensity.get(norm_list.get(4)).get(currentrow)) / ((float) all_intensity.get(norm_list.get(2)).get(currentrow)+ (float) all_intensity.get(norm_list.get(5)).get(currentrow));
            ar.add(sum);
            ratio.put("2",ar);
            
        }
  
    }
    void printratio()
    {
        System.out.println(ratio);
    }
    public ArrayList getallIntensity(String key)
    {
        return all_intensity.get(key);
    }
     public Map getallIntensity()
    {
        return all_intensity;
    }
    public void setIntensity(String key, Float result) 
    {
         intensity.put(key, result);
    }
//    public Float getaverage(String key)
//    {
//        float avg=0;
//        ArrayList<Float> ar=new ArrayList();
//        ar=intensity.get(key);
//        for(int temp=0;temp<ar.size();temp++)
//        {
//            avg+=ar.get(temp);
//        }
//        return avg;
//    }

    
    public int getSpc() 
    {
        return spc;
    }

    public void setSpc(int spc)
    {
        this.spc += spc;
    }

    public ArrayList<Integer> getScannum() 
    {
        return scannum;
    }


    public void addScannum(int number)
    {
        scannum.add(number);
        
    }

 
    public ArrayList<Integer> getCstate() 
    {
        return cstate;
    }

    public void setCstate(int number)
    {
        cstate.add(number);
    }

    public Set<String> getFilename()
    {
        return filename;
    }


    public void setFilename(String filename)
    {
        this.filename.add(filename);
    }

    public int getPeptide_counter(String key) {
        return peptide_counter.get(key);
    }

 
    public void setPeptide_counter(String key) 
    {
        int c=0;
        if(peptide_counter.containsKey(key))
        {
            c=peptide_counter.get(key);
            
        }
         c++;
        peptide_counter.put(key, c);
    }


    public DescriptiveStatistics getstat(String key) 
    {
        return stats.get(key);
    }
    public void setstat(String key,  DescriptiveStatistics stat) 
    {
        this.stats.put(key, stat);
    }
}
