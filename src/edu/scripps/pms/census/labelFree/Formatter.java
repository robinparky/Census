/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package edu.scripps.pms.census.labelFree;
/*
 * Formatter.java
 *
 * Created on January 2, 2007, 10:19 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

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
    private static DecimalFormat format5 = new DecimalFormat("0.00000");    
    private static DecimalFormat sciformat =  new DecimalFormat("0.###E0");
    
    public static String formatDecimal(double d)
    {
        return format.format(d);
    }
    
    public static String formatDecimal5(double d)
    {
        return format5.format(d);
    }
 
    
		public static double sciRound(double d){      	
        return Double.parseDouble(sciformat.format(d));      
    }
 
    public static double round(double rval, int rpl) {
	double p = (double)Math.pow(10,rpl);
	rval = rval * p;
	double tmp = Math.round(rval);
	return (double)tmp/p;
    }

    public static void main(String[] args) {
        System.out.println(round(321.437732, 3));
    }
}
