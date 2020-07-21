package edu.scripps.pms.census.conf;

import edu.scripps.pms.census.util.TimsTOFIndex;
import edu.scripps.pms.census.util.TimsTOFXICDB;
import edu.scripps.pms.census.util.dtaselect.Peptide;
import edu.scripps.pms.census.util.dtaselect.Protein;
import edu.scripps.pms.census.util.io.DTASelectFilterReader;
import org.apache.commons.lang3.tuple.Pair;


import java.io.*;
import java.sql.SQLException;
import java.text.DecimalFormat;
import java.util.*;

import static edu.scripps.pms.census.util.TimsTOFXICDB.getPeakAreaEstimate;


public class DTASelectLFQUtil {

    public final static DecimalFormat df = new DecimalFormat("#.##");;

    public static Set<String> retrievePSM(String dtaSelectPath) throws IOException {
        DTASelectFilterReader reader = new DTASelectFilterReader(dtaSelectPath);
        Set<String> psmSet = new HashSet<>();
        for(Iterator<Protein> proteinIterator = reader.getProteins(); proteinIterator.hasNext(); )
        {
            Protein protein = proteinIterator.next();
            for(Peptide peptide: protein.getPeptideList())
            {
                String psm = peptide.getFileNameWithScan();
                psmSet.add(psm);
            }
        }
        reader.close();
        return psmSet;
    }

    public static Map<String, Pair<Double,Double>> buildPSMTimstofPeakAreaMap(String sqlitePath, String spectraFolder, Set<String> psmSet) throws SQLException, IOException {
        TimsTOFXICDB timsTOFXICDB = new TimsTOFXICDB(sqlitePath);
        Map<String, TimsTOFIndex> indexMap = new HashMap<>();
        Map<String,  Pair<Double,Double>> psmTimstoFPeakArea = new HashMap<>();
        for(String psmStr: psmSet)
        {
            String [] arr = psmStr.split("\\.");
            String fileName = arr[0] +".ms2";

            int scanNumber = Integer.parseInt(arr[1]);
            TimsTOFIndex index = indexMap.get(fileName);
            if(index == null)
            {
               index= new TimsTOFIndex(spectraFolder + File.separator+ fileName);
            }

            int parentId = index.getParentId(scanNumber);
            TimsTOFXICDB.TimstofQueryResult queryResult = timsTOFXICDB.queryAndSumParentId(parentId);
            double peakArea = queryResult.getGaussianPeakArea();
            double estimatedPeakArea = getPeakAreaEstimate(queryResult.getSummedList());
            psmTimstoFPeakArea.put(psmStr, Pair.of(peakArea, estimatedPeakArea));

            indexMap.put(fileName, index);
        }
        timsTOFXICDB.close();

        return psmTimstoFPeakArea;
    }

    public static void writeToDTASelectOutput(String dtaSelectInput, String dtaSelectOutput,
                                              Map<String,Pair<Double, Double>> psmPeakAreaMap) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(dtaSelectInput));
        BufferedWriter bw = new BufferedWriter(new FileWriter(dtaSelectOutput));
        String line;
        boolean writeXICMode = false;
        while((line = br.readLine())!=null)
        {
            bw.write(line);
            if(line.startsWith("Unique\t"))
            {
                bw.append("\tXIC\tEstimated_XIC");
                writeXICMode = true;
            }
            else if(line.startsWith("\tProteins"))
            {
                writeXICMode = false;
            }
            else if(writeXICMode == true && (line.charAt(0) == '\t' ||  line.charAt(0) == '*'))
            {
                String [] arr = line.split("\t");
                String psm = arr[1];
                double peakArea = psmPeakAreaMap.get(psm).getLeft();
                double estimate = psmPeakAreaMap.get(psm).getRight();
                String peakAreaStr = df.format( peakArea);
                String estimateStr = df.format(estimate);
                bw.append("\t").append(peakAreaStr).append("\t").append(estimateStr);
            }
            bw.newLine();
        }
        br.close();
        bw.close();
    }


    public static void processDTASelect(String inputPath, String spectraFolder, String sqlitePath, String output) throws IOException, SQLException {
        Set<String> psmSet = retrievePSM(inputPath);
        Map<String, Pair<Double, Double>> psmPeakAreaMap = buildPSMTimstofPeakAreaMap(sqlitePath, spectraFolder, psmSet);
        writeToDTASelectOutput(inputPath, output, psmPeakAreaMap);

    }


    public static void main(String [] args) throws IOException, SQLException {
        if(args.length==4)
        {
            String dtaSelectInput = args[0];
            String ms2DirPath = args[1];
            String sqlitePath = args[2];
            String output = args[3];
            processDTASelect(dtaSelectInput, ms2DirPath, sqlitePath, output);
        }
        else
        {
            System.out.println("Improper input; Proper input is dtaselect_path ms2_dir_path sqlite_path output_path ");
        }

    }


}
