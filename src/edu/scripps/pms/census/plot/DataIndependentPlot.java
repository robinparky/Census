/*
 * DataIndependentPlot.java
 *
 * Created on July 14, 2005, 12:02 AM
 */

package edu.scripps.pms.census.plot;

import ptolemy.plot.*;

import java.awt.*;
import java.util.*;
import java.io.*;

import edu.scripps.pms.census.util.BasePeakFinder;
import edu.scripps.pms.census.util.LinearRegression;
import edu.scripps.pms.census.util.DataIndependentPeakFinder;
import gnu.trove.TObjectLongHashMap;

/**
 *
 * @author  Robin Park
 * @version $Id: DataIndependentPlot.java,v 1.1 2006/10/02 21:59:43 rpark Exp $
 */

public class DataIndependentPlot extends BaseChroPlot {
    
    //private String dtaselectFile;
    //private TObjectLongHashMap map;
    
    private LinearRegression regression=null;
    
    /** Creates a new instance of DataIndependentPlot */
    public DataIndependentPlot() {
        super(null);
        //this.dtaselectFile = dtaselectFile;
      //  map = new TObjectLongHashMap();
    }
    
    //this is for data dependent...
    protected synchronized void _drawPeak(Graphics graphics)
    {   
        //Draw peak area July 31. robin
        /// read DTASelect.txt file        
        //////Peak range calculated from DTASelect.txt file based on XCorr
        
        drawDTASelect(graphics);
                
        /*
        Vector<PlotPoint> samData = (Vector)_points.elementAt(0);
        Vector<PlotPoint> refData = (Vector)_points.elementAt(1);

        long[] samArr = new long[samData.size()];
        long[] refArr = new long[refData.size()];

        //peak range calculated by peak finding algorithm
        double startRange = Double.parseDouble(this.getStartRange());
        double endRange = Double.parseDouble(this.getEndRange());
        
        //int startIndex=0;
        //int endIndex=0;
        

            
            //startRange = samData.elementAt(pFinder.getStart()).x;
            //endRange = samData.elementAt(pFinder.getEnd()).x;
            
            
            
  
            /*
        for (int pointnum = 0; pointnum < samData.size(); pointnum++) {
            samArr[pointnum] = (long)((PlotPoint)samData.get(pointnum)).y;
            refArr[pointnum] = (long)((PlotPoint)refData.get(pointnum)).y;
        }
        

             long xstart = _ulx + (long)((startRange - _xMin) * _xscale);
            long xend = _ulx + (long)((endRange - _xMin) * _xscale);

             graphics.setColor(Color.YELLOW);
            //graphics.fillRect((int)xstart, _uly, (int)(xend-xstart), (int)(_lry-_uly)+50);
            graphics.fillRect((int)xstart, (int)_uly, (int)(xend-xstart), (int)(_lry-_uly));

             xstart = _ulx + (long)((this.getDtaStartRange() - _xMin) * _xscale);
            xend = _ulx + (long)((this.getDtaEndRange() - _xMin) * _xscale);

             graphics.setColor(Color.RED);
            long scaledScanNum = _ulx + (long)((this.getScanNum() - _xMin) * _xscale);
            
            //draw dtaselect range
            graphics.fillRect((int)xstart, (int)(_lry-3), (int)(xend-xstart), 3);

            //draw center scan #
            graphics.setColor(Color.GREEN);            
            graphics.fillRect((int)scaledScanNum-1, (int)_uly, (int)2, (int)(_lry-_uly));

*/         
    }    
}
