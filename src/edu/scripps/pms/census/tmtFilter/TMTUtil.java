/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.scripps.pms.census.tmtFilter;

import edu.scripps.pms.census.hash.IndexedFile;

import java.util.HashSet;
/**
 *
 * @author rpark
 */
public class TMTUtil {
        
    public static int getCorrespondingMs3Scan(IndexedFile iFile, int scanNum) {
        
        
        return iFile.getPrecursorScan(scanNum);
    }
    
    public static boolean isSamePrecursorWithFilename(IndexedFile iFile, int scanShift, int scanNum, HashSet<String> ht, String fileName) {
                                
        int shiftedScanNum = scanNum+scanShift;
        boolean hasScan = ht.contains(fileName + shiftedScanNum);
        boolean isHcd = "HCD".equals(iFile.getScanType(shiftedScanNum));

        double sPrec = iFile.getPrecursorByScanNum(shiftedScanNum);
        double prec = iFile.getPrecursorByScanNum(scanNum);
        boolean samePrec = (sPrec ==prec);

        return (hasScan && isHcd && samePrec);
    }

    public static boolean isSamePrecursor(IndexedFile iFile, int scanShift, int scanNum, String ms2ScanType) {
                                
        int shiftedScanNum = scanNum+scanShift;
        boolean isHcd = ms2ScanType.equals(iFile.getScanType(shiftedScanNum));

        double sPrec = iFile.getPrecursorByScanNum(shiftedScanNum);
        double prec = iFile.getPrecursorByScanNum(scanNum);
        boolean samePrec = (sPrec ==prec);

        return (isHcd && samePrec);
    }
    
}
