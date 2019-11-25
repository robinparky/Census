
package edu.scripps.ptm;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 *
 * @author Harshil
 */
public class FastaProteinReader
 {
    public static void main(String args[])
    {
        FastaProteinReader fr = new FastaProteinReader();
        System.out.println(fr.getProteinSeq("UniProt_mouse_08-23-2010_reversed.fasta", "Q9CSU0-3|RPR1B_MOUSE"));
       
    }
    public String getProteinSeq(String fileName,String access)
    {
    File fastaFile = new File(fileName);
    File fastaIndexFile = new File(fileName + ".index");
    BufferedReader br =null;
    String proteinSeq = null;
    if(fastaIndexFile.exists()){
        try {
            br = new BufferedReader(new FileReader(fastaIndexFile));
            String eachLine;
            int start = 0;
            int end = 1;
            while (null != (eachLine = br.readLine())) {
                String[] ss = eachLine.split("\t");
                
                if (ss[0].equals(access) || ss[0].contains(access)) {
                    start = Integer.parseInt(ss[ss.length - 1]);
                    eachLine = br.readLine();
                    ss = eachLine.split("\t");
                    //end = Integer.parseInt(ss[1]);
                    end = Integer.parseInt(ss[ss.length - 1]);
                    break;
                }
            }
            br.close();
            
            RandomAccessFile rfile = new RandomAccessFile(fileName, "r");
            rfile.seek(start);
            
            int byteSize = (int) (end - start);
            byte[] sequenceBytes = new byte[byteSize];
            rfile.readFully(sequenceBytes);
            rfile.close();
            proteinSeq = new String(sequenceBytes);
            if (proteinSeq.indexOf("\n") < 0) {
                proteinSeq = "\n";
            }
            proteinSeq = proteinSeq.substring(proteinSeq.indexOf("\n"));
//                                              proteinSeq = proteinSeq.replace("\n","");
            Pattern p = Pattern.compile("\\s*|\r|\n");
            Matcher m = p.matcher(proteinSeq);
            proteinSeq = m.replaceAll("");

        } catch (FileNotFoundException ex) {
            Logger.getLogger(FastaProteinReader.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(FastaProteinReader.class.getName()).log(Level.SEVERE, null, ex);
        }
            }
                return proteinSeq;

    }
}