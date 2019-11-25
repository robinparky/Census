package scripts.pepXML;

import edu.scripps.pms.census.exception.CensusIndexOutOfBoundException;
import edu.scripps.pms.census.hash.IndexedFile;
import edu.scripps.pms.census.util.CalcUtilGeneric;

import java.io.*;
import java.nio.Buffer;

/**
 * Created by yateslab on 10/12/17.
 */
public class ExtractScansFromMS2 {

    /*
        Target file is a text file, with tab delimited columns. The first column is ms2 file, 2nd and 3rd column are
        low and high scans, 4th is charge state
        i.e.
        example.ms2 110 110 2

     */

    public static void main(String [] args) throws IOException, CensusIndexOutOfBoundException {
        String directory = args[0];
        String targetsFile = args[1];
        ExtractScans(directory,targetsFile);
    }

    public static void ExtractScans(String directory, String targetsFile) throws IOException, CensusIndexOutOfBoundException {
        File dire = new File(directory);
        File targets = new File(targetsFile);
        String outPutPath = directory+File.separatorChar+"extractedScans";
        File outPutDirectory = new File(outPutPath);
        if(!outPutDirectory.exists()) outPutDirectory.mkdir();
        assert(dire.isDirectory());
        assert(targets.exists());
        BufferedReader br = new BufferedReader(new FileReader(targetsFile));
        String line;
        while((line = br.readLine())!=null)
        {
            String [] arr = line.split("\t");
            String ms2File = arr[0];
            int scanNumber = Integer.parseInt(arr[1]);
            String ms2Path = directory+File.separatorChar+ms2File+".ms2";
            String index = ms2Path+".index";
            IndexedFile ifile = new IndexedFile(index,ms2Path);
            String scan = CalcUtilGeneric.getSpectrumString(ifile,scanNumber);

            String outputFile = outPutPath+File.separatorChar+ms2File+"_"+scanNumber+".ms2";
            printScan(outputFile,scan);
        }
    }


    public static void printScan(String filePath, String scan) throws IOException {
        BufferedWriter bw = new BufferedWriter(new FileWriter(filePath));
        bw.write(scan);
        bw.close();
    }
}
