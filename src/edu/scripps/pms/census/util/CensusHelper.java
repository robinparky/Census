/*
 * CensusHelper.java
 *
 * Created on April 25, 2007, 4:22 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package edu.scripps.pms.census.util;

import java.text.DecimalFormat;
import java.text.NumberFormat;

/**
 *
 * @author rpark
 */
public class CensusHelper {

    /** Creates a new instance of CensusHelper */
    public CensusHelper() {
    }


    public static DecimalFormat format = new DecimalFormat("0.00");
    public static DecimalFormat d3format = new DecimalFormat("0.000");
    public static DecimalFormat d4format = new DecimalFormat("0.0000");
    public static DecimalFormat d5format = new DecimalFormat("0.00000");
    //public static DecimalFormat twoDigitFormat = new DecimalFormat("0.00");
    public static NumberFormat scientificFormat = new DecimalFormat("0.###E0");


}
