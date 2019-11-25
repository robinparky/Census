/*
 * SignalToNoisePanel.java
 *
 * Created on October 17, 2005, 9:11 AM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package edu.scripps.pms.census.plot;

import java.awt.*;
import java.util.*;

import edu.scripps.pms.census.model.*;

/**
 *
 * @author rpark
 */
public class SigNoisePlot extends javax.swing.JPanel {
    
    /** Creates a new instance of SignalToNoisePanel */
    private int width;
    private int height;
    private FragIonList ionList;
    private int thresholdIndex;
    public SigNoisePlot(int width, int height, FragIonList ionList, int thresholdIndex) {
        
        //g.drawRect(1, 1, width-18, height-33);                
        this.width = width-18;
        this.height = height-33;
        
        this.ionList = ionList;
        this.thresholdIndex = thresholdIndex;
        
    }

    
    public void paint(Graphics g) {

        double maxReg = 0;
       
        
        for(Iterator<FragIon> itr=ionList.iterator(); itr.hasNext(); )
        {
            FragIon ion = itr.next();

            if(maxReg<ion.getRegScore())
                maxReg = ion.getRegScore();
        }
        
        double xGridSize=(double)width/ionList.size();
        
        int index=0;
        double scale=0.9;
        
        //System.out.println(width);
        
        g.setColor(Color.GREEN);
        
        for(Iterator<FragIon> itr=ionList.iterator(); itr.hasNext(); )
        {
            FragIon ion = itr.next();
            
            g.drawLine((int)(3+index*xGridSize), (int)(height*scale-ion.getRegScore()*height/maxReg*scale), (int)(3+index*xGridSize), (int)(height*scale));            
            
//            System.out.println("==>>" + (int)(3+index*xGridSize) + " " + (int)(height*scale-ion.getRegScore()*height/maxReg*scale) + " " + (int)(height*scale));
          
            index++;            
            
            //System.out.println("index==>>" + index + " " + this.thresholdIndex + " " + g.getColor().getRGB());
            if(index>this.thresholdIndex)
                g.setColor(Color.GRAY);
                
        }
        
        g.drawRect(1, 1, width, height);                                
    }
    
}
