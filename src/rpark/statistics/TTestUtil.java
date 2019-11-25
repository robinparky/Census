/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rpark.statistics;

import java.util.List;
import org.apache.commons.math3.stat.inference.TestUtils;

/**
 *
 * @author rpark
 */
public class TTestUtil {


    public static double calculatePairedTTest(double[] arr1, double[] arr2)throws Exception {
      return TestUtils.pairedTTest(arr1, arr2);
    }

    public static double calculateUnpairedTTest(double[] arr1, double[] arr2)throws Exception {
      return TestUtils.tTest(arr1, arr2);
    }


  public static double oneSampleTTestBasedOnLog(double[] arr)throws Exception {

            if(arr.length<2) return 1;

            return TestUtils.tTest(0, arr);
        }


        public static double oneSampleTTest(double[] arr)throws Exception {

            return TestUtils.tTest(1, arr);
        }

        public static double oneSampleTTestConvertToLog(double[] arr)throws Exception {

            double[] logArr = new double[arr.length];
            for(int i=0;i<logArr.length;i++) {
                logArr[i] = Math.log(arr[i])/Math.log(2);
            //    System.out.println(arr[i] + " " + logArr[i]);
            }

            return TestUtils.tTest(0, logArr);
        }

        public static void main(String[] args) throws Exception {

        /*
            double[] arr = {
                    0.689102152698532,
                            -1.66361362900047,
                            -4.82334411177401,
                            -4.19139893157585,
                            -5.23811607994774,
                    3.08001337610177,
                            -9.47050362810661,
                    Double.NaN,
                            -1.535012825203,
                            -2.97777863332743,
                            -11.7825013055753,
                            -4.15273549526359,
                            -1.92078633598089,
                            -2.00237213578737,
                            -2.4272784800244,
                            -2.01820164730668,
                            -1.90183035145152,
                            -1.03461967889646

            };
            */

        double[] arr={20,10,4.9};

            //System.out.println( TTestUtil.oneSampleTTest(arr));

            System.out.println( TTestUtil.oneSampleTTestConvertToLog(arr));
            //System.out.println( TTestUtil.oneSampleTTestConvertToLog(arr));


        }
}
