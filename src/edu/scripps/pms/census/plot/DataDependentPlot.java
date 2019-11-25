/*
 * DataDependentPlot.java
 *
 * Created on July 14, 2005, 1:53 PM
 */

package edu.scripps.pms.census.plot;

import ptolemy.plot.*;

import java.awt.*;
import java.util.*;
import java.io.*;

import edu.scripps.pms.census.util.GenericPeakFinder;
/**
 *
 * @author  Robin Park
 * @version $Id: DataDependentPlot.java,v 1.1 2006/10/02 21:59:43 rpark Exp $
 */

public class DataDependentPlot extends BaseChroPlot {
    
    /** Creates a new instance of DataDependentPlot */
    public DataDependentPlot() {
        super(null);
    }

    //this is for data dependent...
    protected synchronized void _drawPeak(Graphics graphics)
    {
        //Draw peak area July 31. robin
        //for (int dataset = _points.size() - 1; dataset >= 0 ; dataset--) {


        super.drawDTASelect(graphics);
        
        /*
         *
         
        long xstart = _ulx + (long)((startRange - _xMin) * _xscale);
        long xend = _ulx + (long)((endRange - _xMin) * _xscale);

         graphics.setColor(Color.yellow);
        graphics.fillRect((int)xstart, _uly, (int)(xend-xstart), (int)(_lry-_uly));

        xstart = _ulx + (long)((this.getDtaStartRange() - _xMin) * _xscale);
        xend = _ulx + (long)((this.getDtaEndRange() - _xMin) * _xscale);

        long scaledScanNum = _ulx + (long)((this.getScanNum() - _xMin) * _xscale);
        
        graphics.setColor(Color.RED);
        
        int dtaWidth = (int)(xend-xstart);
        if(dtaWidth == 0)
            dtaWidth = 2;
        
        graphics.fillRect((int)xstart, (int)(_lry-3), dtaWidth, 3);

        graphics.setColor(Color.GREEN);            
        //long scaledScanNum = _ulx + (long)((this.getScanNum() - _xMin) * _xscale);
        graphics.fillRect((int)scaledScanNum-1, (int)_uly, (int)2, (int)(_lry-_uly));
*/                
    }    
}
