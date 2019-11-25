package edu.scripps.pms.census.util.io;

import edu.scripps.pms.util.seq.Fasta;
import gnu.trove.TObjectIntHashMap;
import gnu.trove.TObjectLongHashMap;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.*;
import java.io.IOException;
import java.io.FileReader;
import java.io.FileInputStream;
import java.io.File;


/**
 * @author  Tao Xu
 * @author  Robin Park
 * @version $Id: FastaReader.java,v 1.1 2014/09/09 19:29:52 rpark Exp $
 */
public class FastaReader
{

    public static final char FIRSTCHAROFDEFLINE = '>';
    public static final int DEFAULTSEQENCELENGTH = 1000;

    // Becareful, might need lots of memory
    public static List <Fasta> getFastaList(InputStream is) throws IOException {
        List fastaList = new LinkedList();
        for (Iterator <Fasta> fastas = getFastas(is); fastas.hasNext();) {
            fastaList.add(fastas.next());
        }
        return fastaList;
    }
    public static Iterator<Fasta> getFastas(String fastaFileName) throws IOException {
        FileInputStream fis = new FileInputStream(fastaFileName);
        return getFastas(fis);
    } 
    public static Iterator <Fasta> getFastas(final InputStream is) throws IOException {
        return new Iterator() {

            private String lastLine = ""; // remember the last line read
            private BufferedReader br;

            {
                br = new BufferedReader(new InputStreamReader(is));
                // remove the potential empty lines and get the first defline
                while ((lastLine = br.readLine()) != null && lastLine.equals(""));

                if (lastLine.charAt(0) != FIRSTCHAROFDEFLINE) {
                    throw new FileFormatUnknownException();
                }
            }

            public boolean hasNext() {
                return lastLine != null;
            }

            public Object next() {

                Fasta fasta = null;
                try {
                    fasta = getFasta();
                }
                catch (IOException e) {
                    e.printStackTrace();
                }

                return fasta;
            }

            public void remove() {
                throw new UnsupportedOperationException("Not supported");
            }

            private Fasta getFasta() throws IOException {

                StringBuffer sb = new StringBuffer(DEFAULTSEQENCELENGTH);
		String defline = lastLine;

                // if the line read is a empty string, ignore it
                while ( (lastLine = br.readLine()) != null && (lastLine.equals("")
                    || lastLine.charAt(0) != FIRSTCHAROFDEFLINE)) {
                    //System.out.println(lastLine);
                    if (!lastLine.equals("")) {
                        String line = lastLine.trim();
                        sb.append(line);
                    }
                }

                // the lastLine should be the defline
                // and sb.toString should be the sequence
                return new Fasta(defline, sb.toString());
            }

            protected void finalize() throws IOException {
                br.close();
                //System.out.println("Finalized");
            }
        };
    }

    public static HashMap<String, Integer> getProteinLengthMap(String filename) throws IOException {
        HashMap<String, Integer> hashMap = new HashMap<>();
        BufferedReader br = new BufferedReader(new FileReader(filename));
        String eachLine;

        while ( ((eachLine = br.readLine()) != null && !eachLine.equals("")) ) {
            String[] arr = eachLine.split("\t");
            hashMap.put(arr[0], Integer.parseInt(arr[1]));
        }

        br.close();

        return hashMap;
    }

    public static TObjectIntHashMap getProteinLengthTMap(String filename) throws IOException {
        TObjectIntHashMap hashMap = new TObjectIntHashMap();
        BufferedReader br = new BufferedReader(new FileReader(filename));
        String eachLine;

        while ( ((eachLine = br.readLine()) != null && !eachLine.equals("")) ) {
            String[] arr = eachLine.split("\t");
            hashMap.put(arr[0], Integer.parseInt(arr[1]));
        }

        br.close();

        return hashMap;
    }

    public static void main(String args[]) throws IOException
    {

        Fasta fasta;
        String defLine;
        //String fastafile = "/home/rpark/test_data/venom_animals_Bothorops_Lonomia_hybrid_NCBI.fasta";
        //String fastafile = "/home/rpark/test_data/venom_before_fixed_2.fasta";
        String fastafile = args[0];
        //String fastafile = "/home/rpark/test_data/test2.fasta";
//        for (Iterator itr = FastaReader.getFastas(new FileInputStream(args[0])); itr.hasNext(); ) {
        
        java.util.HashSet<String> set = new java.util.HashSet<String>();
        
        for (Iterator itr = FastaReader.getFastas(fastafile); itr.hasNext(); ) {
            fasta = (Fasta) itr.next();
            defLine = fasta.getDefline();

            
            if(set.contains(defLine)) continue;
            
            set.add(defLine);
            //System.out.println("defLine==>>" + fasta.getAccession() + "\t\t" + defLine);
            System.out.println(">" + defLine);            
            System.out.println(fasta.getSequence());
        }


    }
}
