/**
 * @file Mzxml2Ms2.java
 * This is the source file for edu.scripps.pms.util.io.Mzxml2Ms2
 * @author Tao Xu
 * @author Robin Park 
 * @date $Date
 */

package edu.scripps.pms.util;


import java.io.*;
import java.util.List;
import java.util.Stack;
import java.util.Set;
import java.util.HashMap;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.ListIterator;
import java.util.Collections;
import edu.scripps.pms.util.spectrum.*;
import edu.scripps.pms.census.util.io.MzxmlSpectrumReader;
import edu.scripps.pms.util.MZXmlHandler;
import org.jdom.Document;
import org.jdom.Element;
import org.jdom.JDOMException;
//import org.jdom.Namespace;
import org.jdom.input.SAXBuilder;
import edu.scripps.pms.census.ProgressMzxml2MS;
import edu.scripps.pms.census.*;

public class Mzxml2Ms {

    public static final String USAGE = "USAGE: mzxml2ms2 foldername mzxmlfileextension";
    private static ArrayList<String> getFiles(String dir, String extension) throws IOException {
        
         ArrayList<String> files = new ArrayList<String>();
         File currentDir = new File(dir);
         
         for(String s : currentDir.list()) {
             if(s.endsWith(extension)) {
                 String [] arr = s.split("." + extension);
                 files.add(arr[0]);
             }
         }
         return files;    
    }

    public static void converMzXML2MS(String folder, ChroProgressDialog progressDialog) throws Exception {
       
       
       if(!folder.endsWith(File.separator))
	    folder += File.separator;
	    
        String fileextension = "mzXML";
        ArrayList<String> files = getFiles(folder, fileextension);

	if(files.size()<=0)
	    throw new Exception("Failed to find mzXML files.  Please check if the file extensions are mzXML");
            
        if(null != progressDialog)
            progressDialog.addMessage("Preparing to convert mzXML files in " + folder + " to MS1 files ...");
        
        System.out.println("Preparing to convert mzXML files in " + folder + " to MS1 files ...");

        int fileSize = files.size();
        int count=0;
        
        for(String file : files) {
            count++;
            String mzxmlfile = folder + file + "." + fileextension;

            String ms1file = folder + file + ".ms1";
            MzxmlSpectrumReader ms1msr = new MzxmlSpectrumReader(mzxmlfile);
            System.out.println("Converting " + mzxmlfile + " to " + ms1file + " now. It may take a while...");
            
            if(null != progressDialog)
                progressDialog.addMessage("Converting " + mzxmlfile + " to " + ms1file + " now. It may take a while...");
            PrintStream psms1 = new PrintStream(ms1file);
	    psms1.println("H\t");
            int numSpectra = 0;
       
            int mslevel = 1;
            Iterator<MzxmlPeakList> it = ms1msr.getSpectra(mslevel);
            int i = 0;
            while(it.hasNext()) {
                MzxmlPeakList mpl = it.next();
                psms1.print(mpl.getSpectrumWithoutHlines()); 
            }         

            psms1.close();
            
            if(null != progressDialog)
                progressDialog.setProgress( (int)((float)count/fileSize*100) );
            
           // System.out.println("===" + 100*(int)((float)count/fileSize) );
        }

        
       // progress.updateProgress(aa);
    }

    public static void converMzXML2MS2(String folder, ChroProgressDialog progressDialog) throws Exception {
       
       
       if(!folder.endsWith(File.separator))
	    folder += File.separator;
	    
        String fileextension = "mzXML";
        ArrayList<String> files = getFiles(folder, fileextension);

	if(files.size()<=0)
	    throw new Exception("Failed to find mzXML files.  Please check if the file extensions are mzXML");
            
        if(null != progressDialog)
            progressDialog.addMessage("Preparing to convert mzXML files in " + folder + " to MS2 files ...");
        
        System.out.println("Preparing to convert mzXML files in " + folder + " to MS2 files ...");

        int fileSize = files.size();
        int count=0;
        
        for(String file : files) {
            count++;
            String mzxmlfile = folder + file + "." + fileextension;

            String ms2file = folder + file + ".ms2";
            MzxmlSpectrumReader ms2msr = new MzxmlSpectrumReader(mzxmlfile);
            System.out.println("Converting " + mzxmlfile + " to " + ms2file + " now. It may take a while...");
            
            if(null != progressDialog)
                progressDialog.addMessage("Converting " + mzxmlfile + " to " + ms2file + " now. It may take a while...");
            PrintStream ps = new PrintStream(ms2file);
	    ps.println("H\t");
            int numSpectra = 0;
       
            int mslevel = 2;
            Iterator<MzxmlPeakList> it = ms2msr.getSpectra(mslevel);
            int i = 0;
            while(it.hasNext()) {
                MzxmlPeakList mpl = it.next();
                ps.print(mpl.getSpectrumWithoutHlines()); 
            }         

            ps.close();
            
            if(null != progressDialog)
                progressDialog.setProgress( (int)((float)count/fileSize*100) );
            
        }
    }
    
    public static void converMzXML2MS(String folder) throws Exception {

        converMzXML2MS(folder, null);
        
    }

    public static void converMzXML2MS2(String folder) throws Exception {

        converMzXML2MS2(folder, null);
        
    }

    // for testing 
    public static void main(String args[]) throws Exception {

	long start = System.currentTimeMillis();
	System.out.println( start );

        System.out.println(USAGE);
        String fileextension = args[1];
        String folder = args[0];
        ArrayList<String> files = getFiles(folder, fileextension);
        System.out.println("Preparing to convert mzxml files in " + folder + " to ms2 files ...");
        for(String file : files) {
            String mzxmlfile = file + "." + fileextension;

            String ms1file = file + ".ms1";
            MzxmlSpectrumReader ms1msr = new MzxmlSpectrumReader(mzxmlfile);
            System.out.println("Converting " + mzxmlfile + " to " + ms1file + " now. It may take a while, so please be patient ...");
            PrintStream psms1 = new PrintStream(ms1file);
	    psms1.println("H\t");
            //System.out.println("Number of scans: " + scan2Position.size());
            int numSpectra = 0;
       
            int mslevel = 1;
            Iterator<MzxmlPeakList> it = ms1msr.getSpectra(mslevel);
            int i = 0;
            //for(it = msr.getSpectraWithChildren(); it.hasNext();) {
            while(it.hasNext()) {
                //msr.printMzxmlPeakList(it.next());
                MzxmlPeakList mpl = it.next();
                psms1.print(mpl.getSpectrumWithoutHlines()); 
                //System.out.println("\nscanNumber: " + mpl.getLoscan() + "\tmsLevel: " + mpl.getMsLevel() + "\tnumChildren: " + mpl.getNumChildSpectrum());
                //msr.printMzxmlPeakList(it.next());      
        
            }         
            psms1.close();

	System.out.println( System.currentTimeMillis() -  start);
        }
    }

}
