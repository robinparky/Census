/*
 * BasePeakFinder.java
 *
 * Created on July 24, 2005, 9:27 AM
 */

package edu.scripps.pms.census.util;

import java.io.IOException;
/**
 *
 * @author rpark
 */
public class BasePeakFinder {

    //PEAK FIND VARIABLE
    protected static final double AVG_THRESHOLD = 0.15; //threashold factor for average.  Below the average is not considered as peaks.
    protected static final int NUM_POINT= 7;	
    protected static final int WIDTH = (NUM_POINT-1) / 2;	
    protected static final int DERIVATIVE_ORDER=0; // 0 for smooth 
    protected static final int POLYNOMIAL_ORDER=2; // 2 for quadradic, 3 for cubic
    static protected double[] weightArr = new double[NUM_POINT];

    protected long[] smoothArr;
    protected long[] arr; //combined array    
    protected int size;

    protected int start=0;
    protected int end=0;    
    //static protected double[] weightArr = new double[NUM_POINT];
    
    static {
        for(int i=-WIDTH;i<=WIDTH;i++)
            weightArr[i+WIDTH] = sgWeight(i, 0, WIDTH, POLYNOMIAL_ORDER, DERIVATIVE_ORDER);        
    }
    /** Creates a new instance of BasePeakFinder */
    public BasePeakFinder() {}
    
    public BasePeakFinder(long[] arr) {
        size = arr.length;
        
        this.smoothArr = new long[size];
        this.arr = arr;
        
        smoothCurve();
    }

    public void calc() throws IOException {}
    
    public static int[] smoothCurve(int[] arr)
    {        
        int smoothArr[] = new int[arr.length];
        
        //Smoothing for the points 0 to m
        for(int i=0;i<=WIDTH;i++)	
        {
            int t = -WIDTH;
            int cnt3=0;

            for(int j=-WIDTH;j<=WIDTH;j++)
            {
                smoothArr[i] = (int)(sgWeight(j, -WIDTH, WIDTH, POLYNOMIAL_ORDER, DERIVATIVE_ORDER) * arr[i+cnt3]) + smoothArr[i];
                 cnt3++;
            }
        }

        //Smoothing for the bulk of the chromatogram at t = 0
        for(int i=WIDTH+1;i<arr.length-(WIDTH+1);i++)
        {
            int cnt2 = -WIDTH;

            for(int j=1;j<=NUM_POINT;j++)
            {
                smoothArr[i] = (int)(weightArr[j-1] * arr[i+cnt2]) + smoothArr[i];
                cnt2 = cnt2 + 1;
            }
        }

        //Smoothing for the points Ubound(chro)-m to UBound(chro)
        for(int i=arr.length-WIDTH-1;i<arr.length;i++)
        {
            int t=1;
            int cnt3 = -NUM_POINT;

            for(int j=-WIDTH;j<=WIDTH;j++)
            {
                smoothArr[i] = (int)(sgWeight(j, t, WIDTH, POLYNOMIAL_ORDER, DERIVATIVE_ORDER) * arr[i+cnt3]) + smoothArr[i];
    //		System.out.println(newArr[i][refsam] + " " + j + " " + t + " " + WIDTH + " " + POLYNOMIAL_ORDER + " " + DERIVATIVE_ORDER + " " + i + " " + cnt3);
                cnt3 = cnt3 + 1;

            }
        }         
 
        return smoothArr;
        
    }
    
        //smooth the curve
    public void smoothCurve()
    {
        //for(int i=-WIDTH;i<=WIDTH;i++)
        //        weightArr[i+WIDTH] = sgWeight(i, 0, WIDTH, POLYNOMIAL_ORDER, DERIVATIVE_ORDER);

        //Smoothing for the points 0 to m
        for(int i=0;i<=WIDTH;i++)	
        {
            int t = -WIDTH;
            int cnt3=0;

            for(int j=-WIDTH;j<=WIDTH;j++)
            {
                smoothArr[i] = (long)(sgWeight(j, -WIDTH, WIDTH, POLYNOMIAL_ORDER, DERIVATIVE_ORDER) * arr[i+cnt3]) + smoothArr[i];
                    //System.out.println(i + " " + refsam + " " + newArr[i][refsam] + " " + j + " " + (-WIDTH) + " " + WIDTH + " " + POLYNOMIAL_ORDER + " " + DERIVATIVE_ORDER + " " + (i+cnt3) + " " + refsam + " " + cnt3);
                cnt3++;
            }
        }

        //Smoothing for the bulk of the chromatogram at t = 0
        for(int i=WIDTH+1;i<size-(WIDTH+1);i++)
        {
            int cnt2 = -WIDTH;

            for(int j=1;j<=NUM_POINT;j++)
            {
                smoothArr[i] = (long)(weightArr[j-1] * arr[i+cnt2]) + smoothArr[i];
                cnt2 = cnt2 + 1;
            }
        }

        //Smoothing for the points Ubound(chro)-m to UBound(chro)
        for(int i=size-WIDTH-1;i<size;i++)
        {
            int t=1;
            int cnt3 = -NUM_POINT;

            for(int j=-WIDTH;j<=WIDTH;j++)
            {
                smoothArr[i] = (long)(sgWeight(j, t, WIDTH, POLYNOMIAL_ORDER, DERIVATIVE_ORDER) * arr[i+cnt3]) + smoothArr[i];
    //		System.out.println(newArr[i][refsam] + " " + j + " " + t + " " + WIDTH + " " + POLYNOMIAL_ORDER + " " + DERIVATIVE_ORDER + " " + i + " " + cnt3);
                cnt3 = cnt3 + 1;

            }
        }         
    }
    
    public static double sgWeight(int a, int t, int m, int n, int s)
    {
        double sum=0;

        //Calculates the Savitsky-Golay weight of the i'th data point
        //for the t'th Least-Square point of the s'th derivative
        //over 2m+1 points, order n.

        for(int i=0;i<=n;i++)
        {
            sum = sum + (double)(2 * i + 1) * (sgGenFact(2 * m, i) / (double)sgGenFact(2 * m + i + 1, i + 1)) * sgGramPoly(a, m, i, 0) * sgGramPoly(t, m, i, s);
        }

        return sum;

        //Calculate the Savitsky-Golay Weight

        //System.out.println("==>>");
        //System.out.println(sgGenFact(5, 3));
        //System.out.println(sgGramPoly(2,2,3,2));
        //System.out.println(sgGramPoly(5,7,10,8));

    }	

    public static int sgGenFact(int a, int b)
    {
        int j;
        int factor=1;

        int temp = a-b+1;

        for(int i=temp;i<=a;i++)
        {
                factor *= i;
        }

        return factor;

    }

    public static double sgGramPoly(int i, int m, int k, int s)
    {

        //Calculates the Savitsky-Golay Gram Polynomial (s = 0) or it's s'th
        //derivative evaluated at i, order k, over 2m+1 points.

        if(k>0)
                return (double)(4 * k - 2) / (k * (2 * m - k + 1)) * (i * sgGramPoly(i, m, k - 1, s) + s * sgGramPoly(i, m, k - 1, s - 1)) - (double)((k - 1) * (2 * m + k)) / (k * (2 * m - k + 1)) * sgGramPoly(i, m, k - 2, s);
//(4*k-2)/(k*(2*m-k+1))*(i*sgGramPoly(i,m,k-1,s) + s*sgGramPoly(i,m,k-1,s-1)) - ((k-1) * (2*m+k)) / (k*(2*m-k+1)) * sgGramPoly(i,m,k-2,s);
        else if(k == 0 && s == 0)
                return 1.0;
        else
                return 0.0;

    }
    
    //the index to start peak area
    public int getStart()
    {
        return start;
    }
    
    //the index to end peak area
    public int getEnd()
    {
        return end;
    }     

    public long[] getSmoothArr() {
        return smoothArr;
    }
}
