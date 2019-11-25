package edu.scripps.pms.census.junit;

import junit.framework.*;

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
import gnu.trove.TDoubleArrayList;
import gnu.trove.TIntDoubleHashMap;
import gnu.trove.TIntLongHashMap;

import edu.scripps.pms.census.hash.*;
import java.util.*;
import java.io.*;

import edu.scripps.pms.census.conf.*;

public class ParamTestCase extends TestCase {
    
    public ParamTestCase(String name)
    {
        super(name);
    }

    protected void setUp() throws Exception
    {
        super.setUp();
        /**@todo verify the constructors*/
    }


    public void paramTest() throws IOException, Exception
    {
        System.out.println("Validating Census param...");
        
        try
        {
            Configuration conf = Configuration.getInstance();
            conf.readParam(TestConstants.TEST_DATA_HOME, TestConstants.CENSUS_PARAM_FILE);

            Assert.assertEquals( "Is data independent check Error ", conf.isDataIndependent(), true );
            Assert.assertEquals( "Element composition file name check Error ", conf.getElementCompFile(), "N15isotope.ini" );
            Assert.assertEquals( "enrichment check Error ", conf.getEnrichment(), 0.5 );
            Assert.assertEquals( "resolution check Error ", conf.getResolution(), 60000 );
            Assert.assertEquals( "is high res check Error ", conf.isHighRes(), true );
            Assert.assertEquals( "Start mass range check Error ", conf.getStartMassRange(), 400.0 );
            Assert.assertEquals( "End mass range check Error ", conf.getEndMassRange(), 2000.0 );
            Assert.assertEquals( "Tolerance check Error ", conf.getMassTolerance(), 20*0.001 );
            Assert.assertEquals( "Max Window check Error ", conf.getMaxWindow(), 50 );            
            Assert.assertEquals( "window margin check Error ", conf.getMargin(),  300);
        }
        catch(Exception e)
        {
            Assert.fail("Failed " + e);
        }
    }

    protected void tearDown() throws Exception
    {
        super.tearDown();
    }
}
