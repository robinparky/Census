package scripts;


import edu.scripps.pms.census.io.ChroXmlReader;
import edu.scripps.pms.census.model.ChroData;
import edu.scripps.pms.census.model.ChroPeptide;
import edu.scripps.pms.census.model.ChroProtein;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * Created by rpark on 5/2/16.
 */
public class CensusChroParser {
    public static void main(String[] args) throws Exception {
        //ChroXmlReader cr = new ChroXmlReader(args[0]);
        ChroXmlReader cr =
                new ChroXmlReader(
      //            "/home/rpark/census_chro_temp.xml");
      "/data/2/rpark/ip2_data/carolfc/012017_Proteome1_5ug_Velos2/20180524_S2600_GFP_2018_07_19_18_237783/search/temp//census_chro_temp.xml");
                        //"/data/2/rpark/ip2_data//xmeng/histone_test_4_11_2017/DDA_2017_04_12_15_228053/search/temp/census_chro_temp.xml");

        ArrayList<ChroProtein> list = cr.getProteinList();

   //     System.out.println(list.size());
   //     System.out.println(list.get(0).getProteinLine());

        List pepList = list.get(0).getPeptideList();

        ChroPeptide peptide;
        System.out.println(pepList.size());
        for (int i = 0; i < pepList.size(); i++) {

            peptide = (ChroPeptide) pepList.get(i);
           // if (42945 == peptide.getScanNum()) {

          String[] chroDataArr = peptide.getChroData().split(";");


        //  peptide.get
          for(int j=1;j<chroDataArr.length;j++) {
            String[] peakArr = chroDataArr[j].split(" ");

            System.out.println(peakArr[1] + "\t" + peakArr[2]);
          }

          System.out.println("================================");
          String[] peakStr = peptide.getGaussianPeakString().split(";");

          for(String each:peakStr) {
            String[] arr = each.split(" ");
           // System.out.println(arr[0] + "\t" + arr[1]);
          }


          System.out.println("================================");

          /*
                System.out.println(peptide.getDataList().get(0));
                List<ChroData> l = peptide.getDataList();

                System.out.println(peptide.getStartRange() + " " + peptide.getEndRange());
                int start = (int)Double.parseDouble( peptide.getStartRange() );
                int end = (int)Double.parseDouble( peptide.getEndRange() );


                for(ChroData d:l) {

                    if(d.getScanNum()>=start && d.getScanNum()<=end)
                        System.out.println(d.getSampleIntensity() +"\t" + d.getRefIntensity());
                }
                */
              //  System.exit(0);
           // }


        }

        /*
        ArrayList sampleList = cr.getSampleList();
        for (Iterator<String> itr = sampleList.iterator(); itr.hasNext(); ) {
            String each = itr.next();
        }*/
    }
}
