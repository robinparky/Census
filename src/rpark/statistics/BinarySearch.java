/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rpark.statistics;

import java.util.Arrays;

/**
 *
 * @author rpark
 */
public class BinarySearch {
    
    public static int binarySearch(double[] arr, double value) {
        
        int keyIndex = Arrays.binarySearch(arr, value);        
        if(keyIndex<0) //Cannot find index
            keyIndex=-(++keyIndex); //Math.abs(++keyIndex);

        if(keyIndex>=arr.length)
            keyIndex--;
        
        return keyIndex;
                 
    }
    
}
