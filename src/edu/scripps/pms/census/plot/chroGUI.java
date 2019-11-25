/*
 * chroGUI.java
 *
 * Created on March 3, 2006, 1:02 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package edu.scripps.pms.census.plot;

import javax.swing.UIManager;
import java.awt.Toolkit;

/**
 *
 * @author jvenable
 */
public class chroGUI {
    
    /** Creates a new instance of chroGUI */
    public chroGUI() {
    }
    
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        try 
        {
            UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());

            SimplePlot sPlot = new SimplePlot();
            sPlot.Go(args[0]);
        }
        catch (Exception e)
        {
            e.printStackTrace();
        }
    }
    
}
