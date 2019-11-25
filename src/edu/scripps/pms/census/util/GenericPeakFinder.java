/*
 * GenericPeakFinder.java
 *
 * Created on June 29, 2005, 10:27 PM
 */

package edu.scripps.pms.census.util;

import java.io.*;

/**
 *
 * @author  Robin Park
 * @version $Id: GenericPeakFinder.java,v 1.2 2007/02/14 01:20:46 rpark Exp $
 */

public class GenericPeakFinder extends BasePeakFinder {
    
    private long[] samArr;
    private long[] refArr;
   
    
    //private long[] smoothRefArr;

    public GenericPeakFinder() {}
    
    public static void main(String args[]) throws Exception
    {

    //        PeakFinder s = new PeakFinder();
    //        s.calc();
    }

    public GenericPeakFinder(long[] samArr, long[] refArr)
    {
        this.samArr = samArr;
        this.refArr = refArr;
        
        size = samArr.length;
        
        this.smoothArr = new long[size];
        //this.smoothRefArr = new long[size];
        this.arr = new long[size];
        
        for(int i=0;i<size;i++)
            this.arr[i] = this.samArr[i] * this.refArr[i];
        
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
        int mid = size/2;
        double prev = 0;
        double current = 0;
        for (int i=1;i<size-1;i++) 
        {
            if (smoothArr[i]>threshold && smoothArr[i]>smoothArr[i-1] && smoothArr[i]>smoothArr[i+1])
            {
                current = smoothArr[i] * (1- (double)Math.abs(mid-i)/size);
                
                if(current > prev)
                {
                    peak = i;
                    prev = current;
                }
            }
        }
        
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
