/*
 * PeptideDistPlot.java
 *
 * Created on October 19, 2006, 3:43 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package edu.scripps.pms.census.plot;


import java.awt.*;
import java.util.*;
import ptolemy.plot.*;

import edu.scripps.pms.census.RelExMainFrame;
import edu.scripps.pms.census.util.CalcUtil;
import edu.scripps.pms.census.model.WeightedProtein;
/**
 *
 * @author rpark
 */
public class PeptideDistPlot extends Plot {
    
    private RelExMainFrame mainFrame;
    
    private ArrayList arrayList = new ArrayList();
    
    public void clearData()
    {
        this.arrayList.clear();
    }
    
    public void addData(double x, double y)
    {
        this.arrayList.add(new Data(x, y));
    }
    
    public class Data
    {
        private double x;
        private double y;
        
        public Data(double x, double y)
        {
            this.setX(x);
            this.setY(y);
        }

        public double getX() {
            return x;
        }

        public void setX(double x) {
            this.x = x;
        }

        public double getY() {
            return y;
        }

        public void setY(double y) {
            this.y = y;
        }
        
        
        
    }
    /** Creates a new instance of PeptideDistPlot */
    public PeptideDistPlot() {
    }

    public PeptideDistPlot(RelExMainFrame mainFrame) {
        super();
        
        this.mainFrame = mainFrame;
    }
        
    public void drawAdditional(Graphics graphics)
    {
        //center line
        graphics.setColor(Color.GRAY);        
        long center = _ulx + (long)((0 - _xMin) * _xscale);
        //long xend = _ulx + (long)((0.005 - _xMin) * _xscale);
        graphics.fillRect((int)(center-1), (int)_uly, (int)1, (int)(_lry-_uly));

        WeightedProtein.ProteinModel pModel = new WeightedProtein.ProteinModel();        
        
        graphics.setColor(Color.BLUE);
        for(Iterator<Data> itr=this.arrayList.iterator(); itr.hasNext(); )
        {
            Data data = itr.next();
            data.getX();
               
            double invStdev = CalcUtil.getWeightedStdev(data.getY());
            
            pModel.add(invStdev, Math.exp(data.getX()));                                            
            
            long eachY = _lry - (int)((data.getY()-_yMin)*_yscale);
            long eachX = _ulx + (long)((data.getX() - _xMin) * _xscale);
            
            graphics.fillRect((int)(eachX-2), (int)(eachY-2), 4, 4);    
        }
        
        if(this.arrayList.size()>0)
        {
            graphics.setColor(Color.RED);
            long ratio = _ulx + (long)((Math.log( pModel.getStandardWeightedAverage() ) - _xMin) * _xscale);
            graphics.fillRect((int)(ratio-1), (int)_uly, (int)2, (int)(_lry-_uly));
            
        }

        /*
        xstart = _ulx + (long)((this.getDtaStartRange() - _xMin) * _xscale);
        xend = _ulx + (long)((this.getDtaEndRange() - _xMin) * _xscale);
        
        graphics.setColor(Color.GREEN);            
        //long scaledScanNum = _ulx + (long)((this.getScanNum() - _xMin) * _xscale);
        graphics.fillRect((int)scaledScanNum-1, (int)_uly, (int)2, (int)(_lry-_uly));
        
        this.mainFrame.changePeakArea((int)startRange, (int)endRange);
        
        graphics.setColor(Color.RED);
        int dtaWidth = (int)(xend-xstart);
        if(dtaWidth == 0)
            dtaWidth = 3;
        
        */
        
        //long xstart = _ulx + (long)((0.5 - _xMin) * _xscale);
        //long xend = _ulx + (long)((0.7 - _xMin) * _xscale);
    /*    
        //long yPoint = _uly + (long)((3 - _yMin) * _yscale);
        long yPoint = _lry - (int)((3-_yMin)*_yscale);
        //long yvalue = _ul + (long)((0.7 - _xMin) * _xscale);
        
        //graphics.fillRect((int)xstart-3, (int)(_lry-3), dtaWidth, 3);      
        
        
        
        //long scaledScanNum = _ulx + (long)((this.getScanNum() - _xMin) * _xscale);            
                    
        //graphics.fillRect((int)xstart, (int)_uly, (int)(xend-xstart), (int)(_lry-_uly));

        long xstart1 = _ulx + (long)((0.1 - _xMin) * _xscale);
                
        //graphics.fillRect((int)xstart, (int)(_uly+3), 20, 20);
        
        graphics.fillRect((int)xstart1, (int)(yPoint), 10, 10);
  */      
        
        
        //xstart = _ulx + (long)((this.getDtaStartRange() - _xMin) * _xscale);
        //xend = _ulx + (long)((this.getDtaEndRange() - _xMin) * _xscale);
        
        
        
     //   super.drawDTASelect(graphics);
    }
    
}
