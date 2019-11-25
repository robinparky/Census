
/*
* Copyright (c) 2008 The Scripps Research Institute, Yates Lab.  All rights reserved.  
*/

package edu.scripps.pms.util;

/**
 *
 * @author Sung Kyu, Robin, Park
 * @email rpark@scripps.edu
 * Created on Jun 24, 2008 
 * $Revision: 1.1 $
 * $Date: 2008/09/09 22:31:05 $
 */

import java.math.BigDecimal;

public class MathUtil {

    public MathUtil() {

    }
    
        /**
     * Scale decimal number via the rounding mode BigDecimal.ROUND_HALF_UP.
     * 
     * @param value Decimal value.
     * @param scale New scale.
     * @return Scaled number.
     * @since 1.8.3
     */
    static public double getScaled(double value, int scale) {
        double result = value; //default: unscaled
 
        //use BigDecimal String constructor as this is the only exact way for double values
        result = new BigDecimal(""+value).setScale(scale, BigDecimal.ROUND_HALF_UP).doubleValue();
        
        return result;
    }//getScaled()
}
