/*
 * ChromatographPlot.java
 *
 * Created on January 30, 2006, 2:54 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package edu.scripps.pms.census.plot;

import ptolemy.plot.*;

import java.awt.*;
import java.util.*;
import java.io.*;

import edu.scripps.pms.census.RelExMainFrame;
import edu.scripps.pms.census.util.BasePeakFinder;
import edu.scripps.pms.census.util.LinearRegression;
import edu.scripps.pms.census.util.DataIndependentPeakFinder;
import gnu.trove.TObjectLongHashMap;

/**
 *
 * @author rpark
 */
public class ChromatographPlot  extends BaseChroPlot {
            
    public ChromatographPlot(RelExMainFrame mainFrame) {
        super(mainFrame);
        
        //this.dtaselectFile = dtaselectFile;
      //  map = new TObjectLongHashMap();
    }
    
    //this is for data dependent...
    protected synchronized void _drawPeak(Graphics graphics)
    {   
        
        drawDTASelect(graphics);
           
    }    
}
