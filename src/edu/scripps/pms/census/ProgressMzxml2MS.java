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

import edu.scripps.pms.util.Mzxml2Ms;

/**
 *
 * @author rpark
 */
public class ProgressMzxml2MS extends SwingWorker<Void, Void> {
    /*
     * Main task. Executed in background thread.
     */
    private Container container;
    private JDialog progressDialog;
    private String path;

    public ProgressMzxml2MS(Container container, String path, JDialog progressDialog)
    {
        this.container = container;
        this.progressDialog = progressDialog;
        this.path = path;
        
    }

    public void updateProgress(int value) {
        setProgress(value);
    }
    
    public Void doInBackground() {

        try {
            
            System.out.println("Start converting mzXML to MS1 files...");
            //Mzxml2Ms.converMzXML2MS(path, this);            
            System.out.println("Converting mzXML to MS1 files is completed.");
            
        }
        catch(Exception e)
        {
            e.printStackTrace();
            JOptionPane.showMessageDialog(container, "Failed to convert mzXML to MS1", "Error", JOptionPane.ERROR_MESSAGE);
	    progressDialog.setVisible(false);

        } 

        return null;
    }

    /*
     * Executed in event dispatch thread
     */
    public void done() {
	    progressDialog.setVisible(false);  
        
        
    }
}
