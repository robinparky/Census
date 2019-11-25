/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rpark.statistics;

import java.util.*;
import org.apache.commons.math3.stat.inference.TestUtils;
/**
 *
 * @author rpark
 */
public class AnovaUtil {


  public static double calculateAnovaPvalueConvertToLog(List classes)throws Exception {

    List newList = new ArrayList();

    for(int i=0 ;i<classes.size(); i++)
    {
        double [] arr = (double[]) classes.get(i);
        double[] logArr = new double[arr.length];
        for(int j=0; j<arr.length; j++)
        {
            logArr[j]= Math.log(arr[j])/Math.log(2);
        }
        newList.add(logArr);
    }

    /*
    double[] logArr = new double[arr.length];
    for(int i=0;i<logArr.length;i++) {
      logArr[i] = Math.log(arr[i])/Math.log(2);
      //    System.out.println(arr[i] + " " + logArr[i]);
    }
*/

    return AnovaUtil.calculateAnovaPvalue(newList);
  }

    public static double calculateAnovaPvalue(List classes)throws Exception {


        if(classes.size()<=1) return 1.0;


// for two group, it is same as two sample ttet with equal variance.  it is good for biological replicates (independent samples).  Not good for technical replicates (paired t-test is good for this case

        //List classes = new ArrayList();
        boolean enoughData = true;
        for(Iterator itr=classes.iterator(); itr.hasNext(); ) {
            double[] arr = (double[])itr.next();
            if(arr.length<=1) enoughData=false;
        }

        if(!enoughData) return 1.0;
      if(classes.size()==2 && (((double[])classes.get(0)).length == ((double[])classes.get(1)).length) ) {
        double[] group1 = (double[])classes.get(0);
        double[] group2 = (double[])classes.get(1);
        // if two group, run paired ttest
       // return TTestUtil.calculatePairedTTest(group1, group2);
        return TTestUtil.calculateUnpairedTTest(group1, group2);
      }

        double pvalue = TestUtils.oneWayAnovaPValue(classes);

        return pvalue;
    }

    public static void main(String[] args) throws Exception {


        double[] group1 = {1.111111111111112,1.111111111111};
        //double[] group1 = {15190, 7200};
        double[] group2 = {1.111111111,1.1111};


        System.out.println( TTestUtil.calculatePairedTTest(group1, group2) );
      System.out.println(  TestUtils.tTest(group1, group2) );



      if(true) return;

        List l = new ArrayList();
        double[] d1 = new double[2];
        d1[0] = 7200;
        d1[1] = 15190;


        double[] d2 = new double[2];
        d2[0] = 0;
        d2[1] = 185780;
        l.add(d1);
        l.add(d2);



        System.out.println( "===" + AnovaUtil.calculateAnovaPvalue(l) );
    }
}
