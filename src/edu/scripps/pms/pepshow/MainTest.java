/* Top-level window containing a plotter.

 Copyright (c) 1998-2002 The Regents of the University of California.
 All rights reserved.
 Permission is hereby granted, without written agreement and without
 license or royalty fees, to use, copy, modify, and distribute this
 software and its documentation for any purpose, provided that the above
 copyright notice and the following two paragraphs appear in all copies
 of this software.

 IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY
 FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES
 ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
 THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF
 SUCH DAMAGE.

 THE UNIVERSITY OF CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES,
 INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE
 PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
 CALIFORNIA HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES,
 ENHANCEMENTS, OR MODIFICATIONS.

                                        PT_COPYRIGHT_VERSION_2
                                        COPYRIGHTENDKEY
@ProposedRating Yellow (cxh@eecs.berkeley.edu)
@AcceptedRating Yellow (cxh@eecs.berkeley.edu)
*/

package edu.scripps.pms.pepshow;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Toolkit;
import java.awt.Event;
import java.awt.Graphics;
import java.awt.PrintJob;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.print.PrinterJob;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;

import java.net.URL;
import java.util.StringTokenizer;
import java.util.Vector;
import javax.swing.KeyStroke;
import javax.swing.JPanel;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.JMenuBar;
import javax.swing.JOptionPane;
import javax.swing.filechooser.FileFilter;
import javax.swing.JOptionPane;
import javax.swing.UIManager;

import java.net.MalformedURLException;

import ptolemy.plot.*;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.util.StringTokenizer;
import java.net.URL;
import java.net.MalformedURLException;


import ptolemy.plot.plotml.PlotBoxMLParser;
import ptolemy.plot.plotml.PlotMLParser;

import com.microstar.xml.XmlException;
import java.io.BufferedInputStream;
import java.awt.Toolkit;

public class MainTest extends JFrame {
    ///////////////////////////////////////////////////////////////////
    ////                         public variables                  ////
    public PlotBox plot;
    protected File _directory = null;
    protected File _file = null;

    public static void main(String args[]) {
        try {
            UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());

            //PepShowApplication plot = new PepShowApplication(new Plot());
            MainTest plot1 = new MainTest(new Plot());

            //relex.setVisible(true);
            //relex.validate();
            //relex.pack();
            //plot1.maximize();
            plot1.setSize(Toolkit.getDefaultToolkit().getScreenSize());

        } catch (Exception ex) {
            System.err.println(ex.toString());
            ex.printStackTrace();
        }
    }

    /*
    protected void _read(InputStream in) throws IOException {
                System.out.println("lower _Read");
        // Create a buffered input stream so that mark and reset
        // are supported.
        BufferedInputStream bin = new BufferedInputStream(in);
        // Peek at the file...
        bin.mark(9);
        // Read 8 bytes in case 16-bit encoding is being used.
        byte[] peek = new byte[8];
        bin.read(peek);
        bin.reset();

        System.out.println("===========fdsafdsaf=====");

        if ((new String(peek)).startsWith("<?xm")) {
		System.out.println("fdsafds");
            // file is an XML file.
            PlotBoxMLParser parser = new PlotMLParser((Plot)plot);
            //_newParser();
            try {
                parser.parse(null, bin);
            } catch (Exception ex) {
                String msg;
                if (ex instanceof XmlException) {
                    XmlException xmlex = (XmlException)ex;
                    msg =
                        "PepShow: failed to parse PlotML data:\n"
                        + "line: " + xmlex.getLine()
                        + ", column: " + xmlex.getColumn()
                        + "\nIn entity: " + xmlex.getSystemId()
                        + "\n";
                } else {
                    msg = "PepShow: failed to parse PlotML data:\n";
                }
                System.err.println(msg + ex.toString());
                ex.printStackTrace();
            }
        } else {
		System.out.println("11fdsafds");

            _read(bin);
        }
    }
*/
    /*
    protected PlotBoxMLParser _newParser() {
        if (plot instanceof Plot) {
            return new PlotMLParser((Plot)plot);
        } else {
            return new PlotBoxMLParser(plot);
        }
    }
  */
    
    public MainTest(PlotBox plot) throws Exception {

        // invoke the base class constructor and pass in the argument a Plot
        // object. This makes sure that the plot field is an instance of
        // Plot class.

        //super("PepShowApplication");
                
        this.plot = plot;
        // Background color is a light grey.
        plot.setBackground(new Color(0xe5e5e5));        
        getContentPane().add(plot, BorderLayout.CENTER);    
        
                // FIXME: This should not be hardwired in here.
        setSize(700, 300);
        
        
        // Center.
        Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
        Dimension frameSize = getSize();
        int x = (screenSize.width - frameSize.width) / 2;
        int y = (screenSize.height - frameSize.height) / 2;
        setLocation(x, y);
        
        
        // Handle window closing by exiting the application.
        addWindowListener(new WindowAdapter() {
                public void windowClosing(WindowEvent e) {
                    // Strangely, calling _close() here sends javac into
                    // an infinite loop (in jdk 1.1.4).
                    //              _close();
                    System.exit(0);
                }
            });
        //msPlot();
            
            
        FileInputStream fis = new FileInputStream("e:\\test.xml");
        BufferedInputStream bin = new BufferedInputStream(fis);
            // file is an XML file.
            PlotBoxMLParser parser = new PlotMLParser((Plot)plot);
            //_newParser();
            try {
                parser.parse(null, bin);
            } catch (Exception ex) {
                String msg;
                if (ex instanceof XmlException) {
                    XmlException xmlex = (XmlException)ex;
                    msg =
                        "PepShow: failed to parse PlotML data:\n"
                        + "line: " + xmlex.getLine()
                        + ", column: " + xmlex.getColumn()
                        + "\nIn entity: " + xmlex.getSystemId()
                        + "\n";
                } else {
                    msg = "PepShow: failed to parse PlotML data:\n";
                }
                System.err.println(msg + ex.toString());
                ex.printStackTrace();
            }
        //_read(fis);
        
        

        setVisible(true);
    }
/*
    public void msPlot() {

        try {
            FileInputStream fis = new FileInputStream("e:\\test.xml");

            _read(fis);
        } catch (MalformedURLException e) {
            System.err.println(e.toString());
        } catch (FileNotFoundException e) {
            System.err.println("PlotApplet: file not found: " +e);
        } catch (IOException e) {
            System.err.println("PlotApplet: error reading input file: " +e);
        }
        
        plot.drawLegendPanel();
    }
*/
    /** Print the plot.
     */
    /*
    protected void _print() {
        PrinterJob job = PrinterJob.getPrinterJob();
        job.setPrintable(plot);
        if (job.printDialog()) {
            try {
                job.print();
            } catch (Exception ex) {
                JOptionPane.showMessageDialog(this,
                        "Printing failed:\n" + ex.toString(),
                        "Print Error", JOptionPane.WARNING_MESSAGE);
            }
        }
    }
*/
    /** Read the specified stream.  Derived classes may override this
     *  to support other file formats.
     *  @param base The base for relative file references, or null if
     *   there are not relative file references.
     *  @param in The input stream.
     *  @exception IOException If the stream cannot be read.
     */
    /*
    protected void _read(URL base, InputStream in) throws IOException {
	System.out.println("gog11");
        plot.read(in);
	System.out.println("gog11");
    }*/

}
