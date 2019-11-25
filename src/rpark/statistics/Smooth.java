/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rpark.statistics;

/**
 *
 * @author rpark
 */
public class Smooth {
    public static long[] smooth(long[] input, int windowSize) {
            //long[] input = {3, 4, 5, 2, 3, 4, 5, 6, 7, 4, };
            
            long[] out = new long[input.length-2];

            for(int i=0;i<input.length-windowSize+1;i++) {
                    long sum = input[i+0];
                    sum += input[i+1];
                    sum += input[i+2];

                    long average = sum / windowSize;

                    out[i] = average;
            }

            return out;

    }

    public static void main(String[] args) {
        /*

        0
773376
1836191
738394
1727735
1243792
1537636
1131183
1658066
1588888
961123
1551891
1419243
2309482
1292776
1212351
1238360
1105376
1046221
3204661
330720
1520253
754539
0





        smoothChromArr = {double[24]@836}
 0 = 521913.0
 1 = 669592.0
 2 = 860464.0
 3 = 741984.0
 4 = 901832.0
 5 = 782522.0
 6 = 865377.0
 7 = 875627.0
 8 = 841615.0
 9 = 820380.0
 10 = 786451.0
 11 = 1056123.0
 12 = 1004300.0
 13 = 962921.0
 14 = 748697.0
 15 = 711217.0
 16 = 677991.0
 17 = 1071251.0
 18 = 916320.0
 19 = 1011126.0
 20 = 0.0
 21 = 0.0
 22 = 754539.0
 23 = 0.0
         */
    }

    public static double[] smoothAsDouble(long[] input, int windowSize) {
            //long[] input = {3, 4, 5, 2, 3, 4, 5, 6, 7, 4, };


            double[] out = new double[input.length];

            for(int i=0;i<input.length-windowSize+1;i++) {
                    long sum = input[i+0];
                    sum += input[i+1];
                    sum += input[i+2];

                    long average = sum / windowSize;

                    out[i] = average;
            }

            //fillout last missing two elements
            if(input.length>0)
                out[input.length-1] = input[input.length-1];

            if(input.length>1)
                out[input.length-2] = input[input.length-2];

            return out;

    }


}
