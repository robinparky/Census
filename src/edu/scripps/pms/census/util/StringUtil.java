package edu.scripps.pms.census.util;

import java.util.regex.Pattern;
/**
 *
 * @author  Robin Park
 * @version $Id: StringUtil.java,v 1.3 2014/09/03 18:04:13 rpark Exp $
 */

public class StringUtil
{
    public static boolean startsWithDigitOrUpper(String s) {
	return Pattern.compile("^[A-Z0-9]").matcher(s).find();
    }

    public static boolean startsWithDigit(String s) {
	return Pattern.compile("^[0-9]").matcher(s).find();
    }

    public static String removeIsoControlChar(String s) {
	
	char[] arr = s.toCharArray();

	int count=0;
	for(int i=0;i<arr.length;i++)
	{
	    if(Character.isISOControl(arr[i]))
		arr[i] = ' ';
	}

	return new String(arr);
    }
}

