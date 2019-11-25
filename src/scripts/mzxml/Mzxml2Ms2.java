/**
 * @file Mzxml2Ms2.java
 * This is the source file for edu.scripps.pms.util.io.Mzxml2Ms2
 * @author Tao Xu
 * @author Robin Park 
 * @date $Date
 */



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


public class Mzxml2Ms2 {

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

    // for testing 
    public static void main(String args[]) throws Exception {
        System.out.println(USAGE);
        String fileextension = args[1];
        String folder = args[0];
        ArrayList<String> files = getFiles(folder, fileextension);
        System.out.println("Preparing to convert mzxml files in " + folder + " to ms2 files ...");
        for(String file : files) {
            String mzxmlfile = file + "." + fileextension;
            String ms2file = file + ".ms2";
            MzxmlSpectrumReader msr = new MzxmlSpectrumReader(mzxmlfile);
            System.out.println("Converting " + mzxmlfile + " to " + ms2file + " now. It may take a while, so please be patient ...");
            PrintStream ps = new PrintStream(ms2file);
            //System.out.println("Number of scans: " + scan2Position.size());
            int numSpectra = 0;
       
            int mslevel = 2;
            Iterator<MzxmlPeakList> it = msr.getSpectra(mslevel);
            int i = 0;
            //for(it = msr.getSpectraWithChildren(); it.hasNext();) {
            while(it.hasNext()) {
                //msr.printMzxmlPeakList(it.next());
                MzxmlPeakList mpl = it.next();
                ps.print(mpl.getSpectrumWithoutHlines()); 
                //System.out.println("\nscanNumber: " + mpl.getLoscan() + "\tmsLevel: " + mpl.getMsLevel() + "\tnumChildren: " + mpl.getNumChildSpectrum());
                //msr.printMzxmlPeakList(it.next());      
        
            }         
            ps.close();
        }
    }

}
