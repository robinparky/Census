package scripts;

import edu.scripps.pms.census.io.ChroXmlReader;
import edu.scripps.pms.census.model.ChroPeptide;
import edu.scripps.pms.census.model.ChroProtein;
import org.apache.commons.math3.analysis.integration.IterativeLegendreGaussIntegrator;

import java.io.*;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;

public class TMTLabelfreeCombilePfizer {


    public static void main(String[] args) throws Exception {

        //String inputXmlFile = "/data/2/rpark/ip2_data/carolfc/012017_Proteome1_5ug_Velos2/20180524_S2600_GFP_2018_07_19_18_237783/search/temp//census_chro_temp.xml";
        //String inputCensusOutFile = "/data/2/rpark/ip2_data/carolfc/012017_Proteome1_5ug_Velos2/20180524_S2600_GFP_2018_07_19_18_237783/search/temp/census-out.txt";
        //String inputXmlFile = "/home/rpark/test_data/pfizer/lfree_john/quant/silac/small/census_chro_temp.xml";
        //String inputCensusOutFile = "/home/rpark/test_data/pfizer/lfree_john/quant/silac/small/census-out.txt";

        //String inputXmlFile = "/home/rpark/test_data/pfizer/lfree_john/quant/tmt/small/census_chro_temp.xml";
        //String inputCensusOutFile = "/home/rpark/test_data/pfizer/lfree_john/quant/tmt/small/census-out.txt";
        String inputXmlFile = "/home/rpark/test_data/pfizer/lfree_john/quant/tmt//census_chro_temp.xml";
        String inputCensusOutFile = "/home/rpark/test_data/pfizer/lfree_john/quant/tmt/census-out.txt";

     //   String inputXmlFile = "/home/rpark/test_data/pfizer/lfree_john/quant/tmt//census_chro_temp.xml";
      //  String inputCensusOutFile = "/home/rpark/test_data/pfizer/lfree_john/quant/tmt/a.txt";


        String outputFile = inputCensusOutFile.substring(0, inputCensusOutFile.lastIndexOf("/")+1) + "census_tmt_lfree_merged.txt";
        ChroXmlReader cr = new ChroXmlReader(inputXmlFile);

        //ArrayList<ChroProtein> list = cr.getProteinList();

        //System.out.println(pepList.size());
        Hashtable<String, ChroPeptide> peptideHt = new Hashtable<>();


        for(Iterator<ChroProtein> itr=cr.getProteinList().iterator(); ((Iterator) itr).hasNext(); ) {

            ChroProtein chroProtein = itr.next();
            List pepList = chroProtein.getPeptideList();

            for (int i = 0; i < pepList.size(); i++) {

                ChroPeptide peptide = (ChroPeptide) pepList.get(i);

                peptideHt.put(peptide.getSequence() + "_" +  peptide.getChargeState() + "_" +  peptide.getScanNum() + "_" +  peptide.getFileName(), peptide);

            }
        }
        PrintStream p = new PrintStream( new BufferedOutputStream(new FileOutputStream(outputFile)));
//        p.print(sb.toString());




        BufferedReader br = new BufferedReader(new FileReader(inputCensusOutFile));
        String eachLine;

        while( !(eachLine = br.readLine()).startsWith("H\tSLINE\tUNIQUE") )
            p.println(eachLine);

        p.println(eachLine + "\tLabelfree_PEAK_AREA");

        String[] sarr = eachLine.split("\t");

        int seqIndex=-1;
        int scanIndex=-1;
        int csIndex = -1;
        int fnameIndex = -1;

        for(int i=0;i<sarr.length;i++) {

            if(sarr[i].equals("SEQUENCE"))
                seqIndex = i-1;
            else if(sarr[i].equals("ScanNum"))
                scanIndex = i-1;
            else if(sarr[i].equals("CState"))
                csIndex = i-1;
            else if(sarr[i].equals("Filename"))
                fnameIndex = i-1;
        }

        while( null != (eachLine=br.readLine()) ) {

            if(!eachLine.startsWith("S\t")) {
                p.println(eachLine);
                continue;
            }


            String[] arr = eachLine.split("\t");

            //peptideHt.put(peptide.getSequence() + peptide.getChargeState() + peptide.getScanNum() + peptide.getFileName(), peptide);

            String key = arr[seqIndex] + "_" +  arr[csIndex] + "_" +  arr[scanIndex] + "_" +  arr[fnameIndex];

            ChroPeptide chroPeptide = peptideHt.get(key);
            if(null == chroPeptide)
                p.println(eachLine + "\tNA");
            else
                p.println(eachLine + "\t" + chroPeptide.getPeakArea());

        }


        br.close();
        if(null != p)
            p.close();


        //System.out.println(peptideHt);
        /*
        ArrayList sampleList = cr.getSampleList();
        for (Iterator<String> itr = sampleList.iterator(); itr.hasNext(); ) {
            String each = itr.next();
        }*/
    }



    //for silac
    public static void runSilac() throws Exception {

        //String inputXmlFile = "/data/2/rpark/ip2_data/carolfc/012017_Proteome1_5ug_Velos2/20180524_S2600_GFP_2018_07_19_18_237783/search/temp//census_chro_temp.xml";
        //String inputCensusOutFile = "/data/2/rpark/ip2_data/carolfc/012017_Proteome1_5ug_Velos2/20180524_S2600_GFP_2018_07_19_18_237783/search/temp/census-out.txt";
        //String inputXmlFile = "/home/rpark/test_data/pfizer/lfree_john/quant/silac/small/census_chro_temp.xml";
        //String inputCensusOutFile = "/home/rpark/test_data/pfizer/lfree_john/quant/silac/small/census-out.txt";

        String inputXmlFile = "/home/rpark/test_data/pfizer/lfree_john/quant/tmt/small/census_chro_temp.xml";
        String inputCensusOutFile = "/home/rpark/test_data/pfizer/lfree_john/quant/tmt/small/census-out.txt";
        //String inputXmlFile = "/home/rpark/test_data/pfizer/lfree_john/quant/tmt//census_chro_temp.xml";
        //String inputCensusOutFile = "/home/rpark/test_data/pfizer/lfree_john/quant/tmt/census-out.txt";


        String outputFile = inputCensusOutFile.substring(0, inputCensusOutFile.lastIndexOf("/")+1) + "census_tmt_lfree_merged.txt";
        ChroXmlReader cr = new ChroXmlReader(inputXmlFile);

        ArrayList<ChroProtein> list = cr.getProteinList();
        List pepList = list.get(0).getPeptideList();

        ChroPeptide peptide;
        //System.out.println(pepList.size());
        Hashtable<String, ChroPeptide> peptideHt = new Hashtable<>();
        for (int i = 0; i < pepList.size(); i++) {

            peptide = (ChroPeptide) pepList.get(i);
            // if (42945 == peptide.getScanNum()) {

            //  String[] chroDataArr = peptide.getChroData().split(";");

            /*

            //  peptide.get
            for(int j=1;j<chroDataArr.length;j++) {
                String[] peakArr = chroDataArr[j].split(" ");

                System.out.println(peakArr[1] + "\t" + peakArr[2]);
            }
            */

            peptideHt.put(peptide.getSequence() + "_" +  peptide.getChargeState() + "_" +  peptide.getScanNum() + "_" +  peptide.getFileName(), peptide);

            //System.out.println("================================"  +  " " + peptide.getPeakArea());
            //System.out.println("================================"  + peptide.getFileName() + " " + );

            //  String[] peakStr = peptide.getGaussianPeakString().split(";");

            //for(String each:peakStr) {
            //   String[] arr = each.split(" ");
            // System.out.println(arr[0] + "\t" + arr[1]);
            //  }

        }

        PrintStream p = new PrintStream( new BufferedOutputStream(new FileOutputStream(outputFile)));
//        p.print(sb.toString());




        BufferedReader br = new BufferedReader(new FileReader(inputCensusOutFile));
        String eachLine;

        while( !(eachLine = br.readLine()).startsWith("H\tSLINE\tUNIQUE") )
            p.println(eachLine);

        p.println(eachLine + "\tLabelfree_PEAK_AREA");

        String[] sarr = eachLine.split("\t");

        int seqIndex=-1;
        int scanIndex=-1;
        int csIndex = -1;
        int fnameIndex = -1;

        for(int i=0;i<sarr.length;i++) {

            if(sarr[i].equals("SEQUENCE"))
                seqIndex = i-1;
            else if(sarr[i].equals("SCAN"))
                scanIndex = i-1;
            else if(sarr[i].equals("CS"))
                csIndex = i-1;
            else if(sarr[i].equals("FILE_NAME"))
                fnameIndex = i-1;
        }

        while( null != (eachLine=br.readLine()) ) {

            if(!eachLine.startsWith("S\t")) {
                p.println(eachLine);
                continue;
            }


            String[] arr = eachLine.split("\t");

            //peptideHt.put(peptide.getSequence() + peptide.getChargeState() + peptide.getScanNum() + peptide.getFileName(), peptide);

            String key = arr[seqIndex] + "_" +  arr[csIndex] + "_" +  arr[scanIndex] + "_" +  arr[fnameIndex];

            ChroPeptide chroPeptide = peptideHt.get(key);
            if(null == chroPeptide)
                p.println(eachLine + "\tNA");
            else
                p.println(eachLine + "\t" + chroPeptide.getPeakArea());

        }


        br.close();
        if(null != p)
            p.close();


        //System.out.println(peptideHt);
        /*
        ArrayList sampleList = cr.getSampleList();
        for (Iterator<String> itr = sampleList.iterator(); itr.hasNext(); ) {
            String each = itr.next();
        }*/
    }

}
