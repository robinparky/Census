
package edu.scripps.pms.census.plot;

import java.util.Calendar;
import java.io.*;
import java.util.*;

/**
 * @author Robin Park
 * @version $Id: HeatmapUtil.java,v 1.1 2006/10/02 21:59:43 rpark Exp $
 */
public class HeatmapUtil 
{
    final static String GREEN="00FF00";
    final static String RED="FF0000";

    public static String getHeatColor(double min, double max, double input)
    {
	if(input<=min)
	    return GREEN;

	if(input>=max)
	    return RED;

	double mid = (min+max)/2;

	if(input<=mid) //Green color
	{
	    String str = Integer.toHexString( (int)((mid-input)/(max-mid)*255) );
	    str = (str.length()<2? "0" + str: str);
	    return "00" + str + "00";
	}
	else //Red color
	{
	    String str = Integer.toHexString( (int)((input-mid)/(max-mid)*255) );
	    str = (str.length()<2? "0" + str: str);
	    return str + "0000";
	}
    }
}
