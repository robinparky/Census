package edu.scripps.pms.census.util.io;

import edu.scripps.pms.util.sqt.SQTPeptide;

import java.io.IOException;
import java.util.*;

/**
 * Created by yateslab on 2/21/17.
 */
public class SQTCompare {

    public static void main(String [] args) throws IOException
    {
        String sqt1 = args[0];
        String sqt2 = args[1];

        SQTParser sqtParser1 = new SQTParser(sqt1);
        SQTParser sqtParser2 = new SQTParser(sqt2);

        List<SQTPeptide> sqt1peptide = new ArrayList<>();
        List<SQTPeptide> sqt2peptide = new ArrayList<>();

        List<String> sqt1Strings = new ArrayList<>();
        List<String> sqt2Strings = new ArrayList<>();

        Set<String> sqtSet1 = new HashSet<>();
        Set<String> sqtSet2 = new HashSet<>();

        for(Iterator<SQTPeptide> sqtIt =sqtParser1.getSQTPeptide(); sqtIt.hasNext();)
        {
            StringBuilder sb = new StringBuilder();
            SQTPeptide peptide = sqtIt.next();
            sb.append(peptide.getHiScan());
            sb.append("\t");
            sb.append(peptide.getCalMZ());
            sb.append("\t");
            sb.append(peptide.getNumMlines());
            sqtSet1.add(sb.toString());

        }
        for(Iterator<SQTPeptide> sqtIt =sqtParser2.getSQTPeptide(); sqtIt.hasNext();)
        {
            StringBuilder sb = new StringBuilder();
            SQTPeptide peptide = sqtIt.next();
            sb.append(peptide.getHiScan());
            sb.append("\t");
            sb.append(peptide.getCalMZ());
            sb.append("\t");
            sb.append(peptide.getNumMlines());
            if(sqtSet1.add(sb.toString()))
            {
                System.out.println("difference: "+sb.toString());
            }

        }

    }


}
