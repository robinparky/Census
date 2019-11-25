package scripts.histogram;
/*
 * CenSusReportReader.java
 *
 * Created on February 23, 2006, 3:03 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */
import java.io.*;
import java.util.*;

/**
 *
 * @author rpark
 * @version $Id: PeptideHistogram.java,v 1.2 2014/06/24 22:33:37 rpark Exp $
 */
public class PeptideHistogram {
    
    public static void main(String args[]) throws Exception
    {
        BufferedReader br = new BufferedReader(new FileReader(args[0]));
        String eachLine;

	Histogram h = new Histogram(80, -2, 2); //number 80 is for 0.05 bin size
       
        while( (eachLine = br.readLine()) != null)
	{
	    if(!eachLine.startsWith("S"))
		continue;

	    String[] arr = eachLine.split("\t");
	    double d = Math.log( Double.parseDouble(arr[3]) )/Math.log(2);
	    h.setData(d);
	}

	h.print();
    }
}
