/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.scripps.pms.census.util;

import edu.scripps.pms.census.hash.IndexedFile;
import edu.scripps.pms.census.tmtFilter.TMTUtil;
import java.util.Arrays;
import edu.scripps.pms.census.hash.IndexedFile;

/**
 *
 * @author rpark
 */
public class CommonTools {
    
    public static int getKeyindex (IndexedFile iFile, int scanNum) {
            
        int[] keys = iFile.getKeys();
        int keyIndex = Arrays.binarySearch(keys, scanNum);
                    
        if(keyIndex<0) //Cannot find index
            keyIndex=-(++keyIndex); //Math.abs(++keyIndex);                 

        return keyIndex;
    }        
}
