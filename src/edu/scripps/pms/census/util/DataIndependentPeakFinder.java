/*
 * DataIndependentPeakFinder.java
 *
 * Created on July 23, 2005, 12:39 PM
 */

package edu.scripps.pms.census.util;

import java.io.IOException;

/**
 *
 * @author rpark
 */
public class DataIndependentPeakFinder extends BasePeakFinder {

    private int startIndex;
    private int endIndex;    
    
    /** Creates a new instance of DataIndependentPeakFinder 
        it will find peak only within startIndex and endIndex range.
     */
                
    //find peak in the given range - startIndex and endIndex
    public DataIndependentPeakFinder(long[] arr, int startIndex, int endIndex)
    {
        size = arr.length;
        
        this.smoothArr = new long[size];
        this.arr = arr;
        this.startIndex = startIndex;
        this.endIndex = endIndex;
        
        smoothCurve();
    }

    public void calc() throws IOException
    {
        //Find peaks
        double mean;
        double total = 0;
        
        for(int i=0;i<size;i++)
            total += smoothArr[i];
        
        mean = total / size;
        double threshold = mean + mean * this.AVG_THRESHOLD;
        
        //if the peak is closer to the mid, it is more likely real peak.  
        //we give weight on the each peak based on its distance from the middle point
        // factor : 1 - abs(mid-i)/length
        
        int peak = 0;
        //int mid = size/2;
        double prev = 0;
        double current = 0;
        for (int i=startIndex;i<endIndex-1;i++) 
        {
            if (smoothArr[i]>threshold && smoothArr[i]>smoothArr[i-1] && smoothArr[i]>smoothArr[i+1])
            {
                //current = smoothArr[i];// * (1- (double)Math.abs(mid-i)/size);
                //current = smoothArr[i] * (1- (double)Math.abs(mid-i)/size);
                
//                System.out.println(current + "\t" + (1- (double)Math.abs(mid-i)/size));// + "\t" + mid + "\t" + i);
                if(current > prev)
                {
                    peak = i;
                    prev = current;
                }
                //System.out.println(i); //extremities.add(point);
            }
        }
        
        //peak==0 means (endIndex-startIndex)<=1
        //So, just take middle value.(in this case, it is startIndex)
        if(peak==0)
            peak = (startIndex + endIndex)/2;
        
        //for(int i=0;i<size;i++)
        //    System.out.println("--->>" + peak + "\t" + smoothArr[i] + "\t" + threshold);    
        
        //find the peak area
        for(int i=peak;i>=0;i--)
        {
            if(smoothArr[i]>threshold)
              start = i;
            else
                break;
        }
        
        for(int i=peak;i<size;i++)
        {
            if(smoothArr[i]>threshold)
              end = i;
            else
                break;
        }
    }       
}
