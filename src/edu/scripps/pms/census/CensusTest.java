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
import java.awt.Toolkit;

import edu.scripps.pms.census.util.SimpleFileNameFilter;

import edu.scripps.pms.census.io.*;
import edu.scripps.pms.census.io.parse.ChroXMLParser;
        
import edu.scripps.pms.census.util.RelExFileFilter;
import edu.scripps.pms.census.util.PostCalculation;
import edu.scripps.pms.census.util.CalcUtil;

import edu.scripps.pms.census.model.ChroProtein;
import edu.scripps.pms.census.model.ChroPeptide;
import edu.scripps.pms.census.model.ChroData;
import edu.scripps.pms.census.dialog.*;

import javax.swing.table.DefaultTableModel;

import ptolemy.plot.*;
import ptolemy.plot.plotml.PlotBoxMLParser;
import ptolemy.plot.plotml.PlotMLParser;
import edu.scripps.pms.census.plot.*;
import edu.scripps.pms.census.conf.*;

import edu.scripps.pms.census.util.LinearRegression;
import edu.scripps.pms.census.util.*;
import edu.scripps.pms.census.exception.*;
import org.jdom.*;
import org.jdom.input.*;

import edu.scripps.pms.census.chroalign.*;
import edu.scripps.pms.census.model.NonLabelMappingModel;

import org.apache.commons.cli.*;

/**
 *
 * @author  Robin Park
 * @version $Id: CensusTest.java,v 1.3 2008/03/10 16:09:27 rpark Exp $
 */

public class CensusTest {
    
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
   
    public static void main(String args[]) throws IOException, CensusGeneralException, Exception
    {
	    UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
	    RelExMainFrame census = new RelExMainFrame();
	    //census.openChroFile("/home/rpark/0000.xml");
	    census.openChroFile("/home/rpark/rpark_on_data/project/census/census_paper/data/input/all_none/GNF_johnv/rep1/census_chro.xml");
	    //census.openChroFile("/home/rpark/rpark_on_data/project/census/census_paper/data/input/all_none/GNF_johnv/rep2/census_chro.xml");
	    //census.openChroFile("/home/rpark/rpark_on_data/project/census/census_paper/data/input/all_none/GNF_johnv/rep1/census_chro.xml");
	   // census.exportReportANFinal(false, true, true, true, 0.5, 0.1, 0, false, true, false, 0.2, 5.0);
	    /*
		    final boolean noFilter, 
		    final boolean detSelect, 
		    final boolean pValueSelect, 
		    final boolean filterFragmentIons, 
		    final double detValue, 
		    final double pValue, 
		    final double correctFactorValue, 
		    final boolean isUniquePeptide, 
		    final boolean removeNegative,
		    final double allNoneLowerBound,
		    final double allNoneUpperBound   
	    */
	   


/*
	    census.setVisible(true);
	    census.validate();
	    census.pack();
	    //            relex.maximize();
	    census.setSize(Toolkit.getDefaultToolkit().getScreenSize());
	*/
    }
}
