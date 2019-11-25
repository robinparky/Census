/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package edu.scripps.pms.census;

import javax.swing.*;
import java.awt.*;
import java.io.*;
import java.util.*;
import edu.scripps.pms.census.io.*;
import edu.scripps.pms.census.model.*;
import org.jdom.*;

/**
 *
 * @author rpark
 */
public class ProgressTask extends SwingWorker<Void, Void> {
    /*
     * Main task. Executed in background thread.
     */
    private RelExMainFrame mFrame; 
    private JDialog progressDialog;
    private File chroFile;
    private ArrayList<ChroProtein> proteinList;
    private ChroXmlReader cr;

    public ProgressTask(RelExMainFrame mFrame, File chroFile, JDialog progressDialog)
    {
        this.mFrame = mFrame;
        this.progressDialog = progressDialog;
        this.chroFile = chroFile;
    }

    public void updateProgress(int value) {
        setProgress(value);
    }
    
    public Void doInBackground() {

        try {
            mFrame.preProcessChroOpen();
            
            int initSize=5;
            setProgress(initSize);

	    if(null != chroFile)
		cr = new ChroXmlReader(chroFile);                
	    else 
		cr = new ChroXmlReader(mFrame.getFileUri());                
            //setProgress(15);

            proteinList = cr.getProteinList(this, initSize);
            //list = cr.getProteinList(this);
           // mFrame.setProteinList(proteinList);

        }
        catch(IOException ioe)
        {
            JOptionPane.showMessageDialog(mFrame, "Failed to open a chro file: " + ioe, "Failed to open a chro file", JOptionPane.ERROR_MESSAGE);
            ioe.printStackTrace();

	    progressDialog.setVisible(false);
        }
        catch(JDOMException je)
        {
            JOptionPane.showMessageDialog(mFrame, "Format of chro xml file is incorrect.", "Failed to open a chro file", JOptionPane.ERROR_MESSAGE);
            je.printStackTrace();
	    progressDialog.setVisible(false);
        } 
        catch(Exception e)
        {
            e.printStackTrace();
            JOptionPane.showMessageDialog(mFrame, e.toString(), "Error", JOptionPane.ERROR_MESSAGE);
            JOptionPane.showMessageDialog(mFrame, "System will exit", "Error", JOptionPane.ERROR_MESSAGE);
	    progressDialog.setVisible(false);
            System.exit(0);
        } 

        return null;
    }

    /*
     * Executed in event dispatch thread
     */
    public void done() {
//	    progress.setVisible(false);  
        
        mFrame.postProcessOpenChroFile(cr, proteinList, progressDialog);
        
    }
}
