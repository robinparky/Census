/*
 * DataDependentPeakFinder.java
 *
 * Created on September 30, 2005, 3:34 PM
 */

package edu.scripps.pms.census.util;

import java.io.*;

/**
 *
 * @author rpark
 */
public class DataDependentPeakFinder extends BasePeakFinder {
    
    private int keyIndex;
    private float steepRatioThreshold;
    
    //keyindex is the location that should be included in the peak area
    public DataDependentPeakFinder(int[][] inputArr, int startIndex, int endIndex, float steepRatioThreshold)
    {
        this.keyIndex = inputArr.length/2-startIndex;
        
        this.steepRatioThreshold = steepRatioThreshold;
        size = endIndex-startIndex+1;
       // this.smoothArr = new long[size];
        this.arr = new long[size];
        
        for(int i=0;i<size;i++)
        {
            this.arr[i] = (long)inputArr[startIndex+i][1]; //sample used for peak finding
            //System.out.println(this.arr[i]);
        }   
/*        
        smoothCurve();        
        for(int i=0;i<size;i++)
        {
            System.out.println(smoothArr[i]);
        }
    */
    }
   
    public void calc() throws IOException
    {
        //Find peaks
        double mean;
        
        //if the peak is closer to the mid, it is more likely real peak.  
        //we give weight on the each peak based on its distance from the middle point
        // factor : 1 - abs(mid-i)/length
        
        int steepCount=3;
        
        start = keyIndex;
        end = keyIndex;
        long area1;
        long area2;
        
        for(int i=keyIndex;i>0;i--)
        {
            if(i-3<0)
                break;
            
            area1 = arr[i] + arr[i+1] + arr[i+2];
            area2 = arr[i-1] + arr[i-2] + arr[i-3];
            
            if(area2==0 && area1==0)
                break;

            if((float)area2/area1>steepRatioThreshold)
                break;
            
            start = i;
            
        //    System.out.println("start" + start + " " + arr[i]);
        }
       
        for(int i=keyIndex+1;i<arr.length;i++)
        {
            if(i+3>=arr.length)                
                break;
           
            area1 = arr[i] + arr[i-1] + arr[i-2];
              //  System.out.println("area1 " + arr[i] + " " + arr[i-1] + " " + arr[i-2]);
            area2 = arr[i+1] + arr[i+2] + arr[i+3];
              //  System.out.println("area2 " + arr[i+1] + " " + arr[i+2] + " " + arr[i+3]);
            
            if(area2==0 && area1==0)
                break;

            if((float)area2/area1>steepRatioThreshold)
                break;
            
            end = i;
          //  System.out.println("end" + end + " " + arr[i]);
        }        
                
    }   
/*
    public void calc() throws IOException
    {
        //Find peaks
        double mean;
        double total = 0;
        
        for(int i=0;i<size;i++)
            total += smoothArr[i];
        
        //mean = total / size;
        //double threshold = mean + mean * this.AVG_THRESHOLD;
        
        //if the peak is closer to the mid, it is more likely real peak.  
        //we give weight on the each peak based on its distance from the middle point
        // factor : 1 - abs(mid-i)/length
        
        int steepCount=3;
        
        start = keyIndex;
        end = keyIndex;
        long area1;
        long area2;
        
        for(int i=keyIndex;i>0;i--)
        {
            if(i-3<0)
                break;
            
            area1 = smoothArr[i] + smoothArr[i+1] + smoothArr[i+2];
            area2 = smoothArr[i-1] + smoothArr[i-2] + smoothArr[i-3];
            
            if(area2==0 && area1==0)
                break;

            if((float)area2/area1>steepRatioThreshold)
                break;
            
            start = i;
            
            System.out.println("start" + start + " " + smoothArr[i]);
        }
       
       System.out.println("keyindex + " + keyIndex + " " +  arr.length);
        for(int i=keyIndex+1;i<arr.length;i++)
        {
           System.out.println("fdsafds" + (i+3));
            if(i+3>arr.length)                
                break;
           
           System.out.println("fdsafds");
            area1 = smoothArr[i] + smoothArr[i-1] + smoothArr[i-2];
                System.out.println("area1 " + smoothArr[i] + " " + smoothArr[i-1] + " " + smoothArr[i-2]);
            area2 = smoothArr[i+1] + smoothArr[i+2] + smoothArr[i+3];
                System.out.println("area2 " + smoothArr[i+1] + " " + smoothArr[i+2] + " " + smoothArr[i+3]);
            
            if(area2==0 && area1==0)
                break;

            if((float)area2/area1>steepRatioThreshold)
                break;
            
            end = i;
            System.out.println("end" + end + " " + smoothArr[i]);
        }        
                
    }   
*/    
}
