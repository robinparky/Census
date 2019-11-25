package edu.scripps.pms.census.io;

import edu.scripps.pms.census.util.io.FastaReader;
import edu.scripps.pms.util.seq.Fasta;

import java.io.BufferedOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.Iterator;

/**
 * Created by rpark on 12/12/16.
 */
public class FastaFileProteinLengthCalc {
    public static void main(String args[]) throws Exception
    {
        if(args.length<1)
        {
            printError();

            return;
        }

        //  String file = "/data/2/rpark/ip2_data/rpark/database/wormbase_GutBacteria_32_01-01-2015.fasta";
        System.out.println("calculating protein length...");

        FastaFileProteinLengthCalc f = new FastaFileProteinLengthCalc();
        //f.generateProteinLength(file);
        f.generateProteinLength(args[0]);

        System.out.println("done.");


    }


    public static void generateProteinLength(String filename) throws Exception {

        BufferedOutputStream bo = new BufferedOutputStream(new FileOutputStream(filename + ".length"));
        PrintStream p = new PrintStream(bo);

        for (Iterator itr = FastaReader.getFastas(new FileInputStream(filename)); itr.hasNext(); ) {
            Fasta fasta = (Fasta) itr.next();

            p.println(fasta.getAccession() + "\t" + fasta.getSequence().length());
        }

        bo.close();
        p.close();
    }

    private static void printError()
    {
        System.out.println("Usage : java FastaFileProteinLengthCalc filename");
    }

}
