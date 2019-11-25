/*
 * Formatter.java
 *
 * Created on January 2, 2007, 10:19 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package edu.scripps.pms.census.tools;


import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.text.NumberFormat;

/**
 *
 * @author rpark
 */
public class Formatter {
    
    /** Creates a new instance of Formatter */
    private Formatter() {
    }
    
    private static DecimalFormat format = new DecimalFormat("0.00");    
    private static DecimalFormat formatThree = new DecimalFormat("0.000");    
    private static DecimalFormat formatFour = new DecimalFormat("0.0000");    
    
    public static String formatDecimal(double d)
    {
        return format.format(d);
    }
    
    public static String formatThreeDecimal(double d)
    {
        return formatThree.format(d);
    }
    public static String formatFourDecimal(double d)
    {
        return formatFour.format(d);
    }

    private static DecimalFormat threeDigitFormat =  new DecimalFormat("#.###");


    public static String halfRoundUpThreeDigitFormat(String d) {
        return halfRoundUpThreeDigitFormat(Double.parseDouble(d));

    }

    public static String halfRoundUpThreeDigitFormat(double d) {
        threeDigitFormat.setRoundingMode(RoundingMode.HALF_UP);

        return threeDigitFormat.format(d);
    }

    public static String threeDigitFormat(String s) {

        if(Character.isDigit(s.charAt(0)))

            return threeDigitFormat.format( Double.parseDouble(s) );
        else
            return s;
    }


}
