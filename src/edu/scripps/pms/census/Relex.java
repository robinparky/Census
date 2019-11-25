/*
 * Relax.java
 *
 * Created on March 18, 2005, 3:51 PM
 */

package edu.scripps.pms.census;

import javax.swing.UIManager;
import java.awt.Toolkit;


/**
 *
 * @author  Robin Park
 * @version $Id: Relex.java,v 1.2 2006/10/13 05:50:30 rpark Exp $
 *
 */
public class Relex {
    
    /** Creates a new instance of Relax */
    public Relex() {
    }
    
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
    
        try 
        {

            UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
            RelExMainFrame relex = new RelExMainFrame();
            
            
//            RelExMainFrame.runsimpletest();
//	    //relex.runtest();
            

            relex.setVisible(true);
            relex.validate();
            relex.pack();
//            relex.maximize();
            relex.setSize(Toolkit.getDefaultToolkit().getScreenSize());
//            relex.MAXIMIZED_HORIZ             
            
        //    relex.dummyOpenChroFile();
            
            //MRM /data/1/CF/Quant/Q8WTaspike

//            RelExMainFrame.runNonLabelActionPerformed(null);
        }
        catch (Exception e)
        {
            e.printStackTrace();
        }
    }
    
}
