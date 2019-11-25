/*
 * RelaxMainFrame.java
 *
 * Created on March 18, 2005, 3:52 PM
 */

package edu.scripps.pms.census;

import javax.swing.*;
import java.io.*;
import java.awt.Color;
import java.util.*;
import java.text.DecimalFormat;

import edu.scripps.pms.census.util.SimpleFileNameFilter;

import edu.scripps.pms.census.io.ChroReader;
import edu.scripps.pms.census.io.parse.ChroXMLParser;
        
import edu.scripps.pms.census.util.RelExFileFilter;
import edu.scripps.pms.census.util.PostCalculation;
import edu.scripps.pms.census.util.CalcUtil;

import edu.scripps.pms.census.model.ChroProtein;
import edu.scripps.pms.census.model.ChroPeptide;
import edu.scripps.pms.census.model.ChroData;

import javax.swing.table.DefaultTableModel;

import ptolemy.plot.*;
import ptolemy.plot.plotml.PlotBoxMLParser;
import ptolemy.plot.plotml.PlotMLParser;
import edu.scripps.pms.census.plot.*;
import edu.scripps.pms.census.conf.*;

import edu.scripps.pms.census.util.LinearRegression;
import edu.scripps.pms.census.util.*;
/**
 *
 * @author  Robin Park
 * @version $Id: RelexStand.java,v 1.2 2006/10/20 18:18:44 rpark Exp $
 */

public class RelexStand {
    
    /** Creates new form inRelaxMainFrame */

    private float massTolerance;
    private final String MS1_FILE = "ms1";
    private final String MS2_FILE = "ms2";
    private Configuration conf;
    private String filePath;
    private String elementFile;
    private DecimalFormat twoDigitFormat = new DecimalFormat("0.00");
    
    //private final String PARAM_FILE="relex.param";
    private final String PARAM_FILE="census.param";

    public static void main(String args[]) throws IOException, Exception
    {
//	this.filePath = args[0];

//        long start = System.currentTimeMillis();
//        System.out.println( start );        
        RelexStand relex;

        if(args.length>0)
	    relex = new RelexStand(args[0]);
        else
	    relex = new RelexStand(".");
        
 //       System.out.println( System.currentTimeMillis() -start );
    }

    public RelexStand(String filePath) throws IOException, Exception
    {
	this.filePath = filePath;
        
        init();
        run();
    }

    public void init() throws IOException, Exception
    {
        conf = Configuration.getInstance();
	conf.setStartTime();
        
        conf.readParam(filePath, PARAM_FILE);        
        elementFile = filePath + File.separator + conf.getElementCompFile();
    }
    
    private void run() throws IOException
    {
        // TODO add your handling code here:       

        try {
            JProgressBar bar = new JProgressBar();
            //ChroGenerator chro = new ChroGenerator(bar, null, 50, 50, filePath + File.separator + "N15isotope.ini", massTolerance);
            ChroGenerator chro = new ChroGenerator(bar, null, 50, 50, elementFile, massTolerance);

	    printMode(conf.isDataIndependent(), conf.isHighRes());

            if( conf.isDataIndependent() )
                chro.createMsmsXmlChro(null);
            else
                chro.createFullscanXmlChro();

        }catch(IOException e) {
            System.out.println("Failed to generate a chro file: " + e);
            e.printStackTrace();
        }catch(Exception e) {
            System.out.println("Failed to generate a chro file: " + e);
            e.printStackTrace();
        }

        System.out.print("Running time : ");
        System.out.print( twoDigitFormat.format((float)(System.currentTimeMillis()-conf.getStartTime())/1000/60) );
        System.out.println(" minutes ");
    }

    private void printMode(boolean isDataIndependent, boolean isHighRes)
    {
	System.out.println("\n****************************************************************");
	System.out.println("*                                                              *");
	    System.out.println("*  Running CenSus...                                           *");

	if(isDataIndependent)
	    System.out.println("*  in Data independent mode                                    *");
	else
	    System.out.println("*  in Data dependent mode                                      *");

	if( conf.isHighRes() )
	    System.out.println("*  in High Resolution mode                                     *");
	else
	    System.out.println("*  in Low Resolution mode                                      *");

	System.out.println("*                                                              *");
	System.out.println("****************************************************************");

    }
}
