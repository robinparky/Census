package scripts.labelfree;

import com.mongodb.util.Hash;
import edu.scripps.pms.census.labelFree.LabelfreeMissingPeptideBuilderSplit;
import edu.scripps.pms.census.labelFree.ProteinModel;
import edu.scripps.pms.census.labelFree.model.LabelfreePeptide;
import edu.scripps.pms.census.model.ChroPeptide;
import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.*;

/**
 * Created by rpark on 4/15/16.
 */
public class LabelfreeFilledParser {
    public static void main(String[] args) throws Exception {


        String fname ="/home/rpark/test_data/berkeley_lori/census_labelfree_out_162_filled.txt";

        edu.scripps.pms.census.labelFree.LabelfreeFilledParser l = new edu.scripps.pms.census.labelFree.LabelfreeFilledParser(fname);

        List<ProteinModel> list = l.readWholeFile();
        int count=0;
        Set<String> set = new HashSet<>();

        for(ProteinModel pm:list) {



            List<Double> ilist = pm.getNormIntensityList();

            double avg1 = (ilist.get(0) + ilist.get(1))/2;
            double avg2 = (ilist.get(2) + ilist.get(3))/2;
            if(avg1<=0|| avg2<=0)
                continue;

            System.out.println(Math.log10( avg1) + "\t" + Math.log10( avg2 ));

            /*
            for(LabelfreePeptide labelfreePeptide :pm.getPeptideList()) {

              //  System.out.println(labelfreePeptide);

                String key = labelfreePeptide.getSequence() + labelfreePeptide.getChargeState();

                if(set.contains(key)) continue;

                set.add(key);

                double intensity1 = ( labelfreePeptide.getPeptideList().get(0).getAverageIntensity() +
                        labelfreePeptide.getPeptideList().get(1).getAverageIntensity() )/2;
                if(intensity1>0)
                    intensity1 = Math.log10(intensity1);

                double intensity2 = ( labelfreePeptide.getPeptideList().get(2).getAverageIntensity() +
                        labelfreePeptide.getPeptideList().get(3).getAverageIntensity() )/2;
                if(intensity2>0)
                    intensity2 = Math.log10(intensity2);

                double normIntensity1 = ( labelfreePeptide.getPeptideList().get(0).getNormIntensity() +
                        labelfreePeptide.getPeptideList().get(1).getNormIntensity() )/2;

                if(normIntensity1>0)
                    normIntensity1 = Math.log10(normIntensity1);

                double normIntensity2 = ( labelfreePeptide.getPeptideList().get(2).getNormIntensity() +
                        labelfreePeptide.getPeptideList().get(3).getNormIntensity() )/2;
                if(normIntensity2>0)
                    normIntensity2 = Math.log10(normIntensity2);



                System.out.println(intensity1 + "\t" + intensity2 + "\t" + normIntensity1 + "\t" + normIntensity2);
           //    System.out.println(Math.log10(normInt1) + "\t" + Math.log10(normInt2));

            }

            */

           // if(count>2) break;
           // count++;
        }


        if(true) return;

//        String fname = "/data/2/rpark/ip2_data/xudong/TMTMS3/labelfree_quant/labelfree_12334/census_labelfree_out_12334_filled.txt.small";
        //String fname = "/home/rpark/test_data/amt_analysis/small.txt";

    //    String fname = "/home/rpark/test_data/amt_analysis/nofilled_results/census_labelfree_out_100_filled.txt";
//S       [1]     R.DQFTTTEVNMAR.V        20141026_HeLa_100ng_BEH60_140min_35ms_IC_DE5_5e3_1      75441   2       4903750.0       19251   0.0     1412.6476
// 1412.6475       92436.9 4.2538  0.6422  0.1     1
//       10.820911       2       75394.0 0.0     92.79   4.57    [2]     R.DQFTTTEVNMAR.V        20141026_HeLa_100ng_BEH60_140min_35ms_IC_DE5_5e3_2      76757
// 2       5294877.0       22264   0.0     1412.6473
//       1412.6475       180788.8        4.3165  0.6214  -0.1    1       10.925089       2       76723.0 0.0     92.979  2.75
// [3]     R.DQFTTTEVNMAR.V        20141026_HeLa_100ng_BEH60_140min_35ms_IC_DE5_5e3_3
//      78739   2       5233913.0       23526   0.0     1413.6494       1412.6475       94122.4 4.5748  0.5994  -1.0    1       10.694483       3       78698.0 0.0     93.691  3.89    4935329.771482516       5338781.336675488       5158483.70026052        1.0

        BufferedReader br = new BufferedReader(new FileReader(fname));
        String eachLine;
        Set<String> keyset;
        keyset = new HashSet<>();
        /*int count111=0;
        int count011=0;
        int count101=0;
        int count110=0;
        int count001=0;
        int count010=0;
        int count100=0;

        */
        Hashtable<String, Integer> ht = new Hashtable<>();

        // H	PLINE	LOCUS	AVERAGE_RATIO	AVERAGE_RATIO_REV	STANDARD_DEVIATION	STANDARD_DEVIATION_REV	COMPOSITE_RATIO	COMPOSITE_RATIO_STANDARD_DEVIATION	WEIGHTED_AVERAGE	LOG_INV_AVERAGE	LOG_INV_AVERAGE_REV	PEPTIDE_NUM	SPEC_COUNT	LSPEC_COUNT	HSPEC_COUNT	AREA_RATIO	DESCRIPTION
        while (null != (eachLine=br.readLine())) {

            if(!eachLine.startsWith("S\t"))
                continue;

            String[] arr = eachLine.split("\t");


            String key = arr[2] + "_" + arr[5];
            if(keyset.contains(key))
                continue;

            ((HashSet) keyset).add(key);

            StringBuffer vennInfo = new StringBuffer();
            if(Double.parseDouble(arr[6])>0)
                vennInfo.append("1_");
            else
                vennInfo.append("0_");

            if(Double.parseDouble(arr[6+21])>0)
                vennInfo.append("1_");
            else
                vennInfo.append("0_");

            if(Double.parseDouble(arr[6+21*2])>0)
                vennInfo.append("1_");
            else
                vennInfo.append("0_");

            String vennStr = vennInfo.toString();
            Integer integer = ht.get(vennStr);
            if(null == integer)
                ht.put(vennStr, 1);
            else
                ht.put(vennStr, integer+1);



            System.out.println(key + "\t" + arr[6] + "\t" + arr[6+21] + "\t" + arr[6+21*2] + "\t" + vennStr);
           // System.out.println(eachLine);
          //  break;
        }

   //     List<ProteinModel>  l = LabelfreeFilledParser.read(fname);

        //System.out.println(keyset);
        System.out.println(ht);
    }


    /*
    public static List<ProteinModel> read(String filename) {
        List<ProteinModel> qproList = new ArrayList<>();
        try {
            BufferedReader br = new BufferedReader(new FileReader(filename));
            String eachLine;
            String[] tmpArr;
            // H	PLINE	LOCUS	AVERAGE_RATIO	AVERAGE_RATIO_REV	STANDARD_DEVIATION	STANDARD_DEVIATION_REV	COMPOSITE_RATIO	COMPOSITE_RATIO_STANDARD_DEVIATION	WEIGHTED_AVERAGE	LOG_INV_AVERAGE	LOG_INV_AVERAGE_REV	PEPTIDE_NUM	SPEC_COUNT	LSPEC_COUNT	HSPEC_COUNT	AREA_RATIO	DESCRIPTION
            int locusIndex = -1;
            int ratioIndex = -1;
            int revRatioIndex = -1;
            int stdIndex = -1;
            int revStdIndex = -1;
            int weightRatioIndex = -1;
            int pepNumIndex = -1;
            int pepNumTotalIndex = -1;
            int specCountIndex = -1;
            int descIndex = -1;
            int compositeRatioIndex = -1;
            int compositeRatioStdevIndex = -1;
            List<ProteinQuant> proList = new ArrayList<>();
            eachLine = br.readLine();
            while (eachLine != null) {
                if (eachLine.startsWith("PLINE\tACCESSION")) {
                    tmpArr = eachLine.split("\t");
                    for (int i = 0; i < tmpArr.length; i++) {
                        if (tmpArr[i].equals("ACCESSION")) {
                            locusIndex = i - 1;
                        } /*else if (tmpArr[i].equals("AVERAGE_RATIO_REV")) {
                            revRatioIndex = i - 1;
                        } else if (tmpArr[i].equals("STANDARD_DEVIATION")) {
                            stdIndex = i - 1;
                        } else if (tmpArr[i].equals("STANDARD_DEVIATION_REV")) {
                            revStdIndex = i - 1;
                        } else if (tmpArr[i].equals("WEIGHTED_AVERAGE")) {
                            weightRatioIndex = i - 1;
                        } else if (tmpArr[i].equals("PEPTIDE_NUM")) {
                            pepNumIndex = i - 1;
                        } else if (tmpArr[i].equals("TOTAL_PEPTIDE_NUM")) {
                            pepNumTotalIndex = i - 1;
                        } else if (tmpArr[i].equals("SPEC_COUNT")) {
                            specCountIndex = i - 1;
                        } else if (tmpArr[i].equals("DESCRIPTION")) {
                            descIndex = i - 1;
                        } else if (tmpArr[i].equals("COMPOSITE_RATIO")) {
                            compositeRatioIndex = i - 1;
                        } else if (tmpArr[i].equals("COMPOSITE_RATIO_STANDARD_DEVIATION")) {
                            compositeRatioStdevIndex = i - 1;
                        }
                    }
                }
                //     H	SLINE	UNIQUE	SEQUENCE	RATIO	REV_SLOPE_RATIO	REGRESSION_FACTOR	DETERMINANT_FACTOR	PROBABILITY_SCORE	XCorr	deltaCN	SAM_INT	REF_INT	PEAK_INT	AREA_RATIO	PROFILE_SCORE	FILE_NAME	SCAN	CS	ENRICHMENT	Mod Sequence	Localization Score	Debunker Score

                if (eachLine.startsWith("H\tSLINE\tUNIQUE")) {
                    CensusPeptide.setFeatureIndices(eachLine);
                }
                if (eachLine.startsWith("H\t&SLINE\tUNIQUE")) {
                    CensusPeptide.setSingletonFeatureIndices(eachLine);
                }
                if (eachLine.startsWith("P\t")) {
                    QuantProteinModel qProtein = new QuantProteinModel();
                    while (eachLine.startsWith("P\t")){
                        //    System.out.println("aaaaaaaabbbb " + eachLine);
                        String[] words;
                        words = eachLine.split("\t");
                        ProteinQuant pro = new ProteinQuant();
                        pro.setLocus(words[locusIndex]);

                        qProtein.addQuantProtein(pro);
                        eachLine = br.readLine();
                    }
                    if (eachLine.startsWith("S\t") || eachLine.startsWith("&S")){
                        // read peptide line
                        com.ipa.ip2.util.dtaselect.CensusPeptide pep = null;
                        while(eachLine != null){

                            if (eachLine.startsWith("S\t")) {
                                String[] arr = eachLine.split("\t");
                                pep = new com.ipa.ip2.util.dtaselect.CensusPeptide(arr);
                                qProtein.addQuantPeptide(pep);
                                eachLine = br.readLine();
                            }
                            else if(eachLine.startsWith("&S")) {
                                String[] arr = eachLine.split("\t");
                                pep = new com.ipa.ip2.util.dtaselect.CensusPeptide(arr,"singleton");
                                qProtein.addQuantPeptide(pep);
                                eachLine = br.readLine();
                            }
                            else{
                                break;
                            }

                        }

                        if (eachLine == null ||(eachLine.startsWith("P\t"))){
                            qproList.add(qProtein);
                            proList.clear();
                        }
                    }
                }
                else {
                    eachLine = br.readLine();
                }
            }

        } catch (IOException | NumberFormatException e) {
            e.printStackTrace();
        }
        return qproList;
    }

                        */

/*
    public static double getRatioMedianValue(List<Double> list) {

        DescriptiveStatistics stats = new DescriptiveStatistics();

        for (Double d : list) {
            // if(d<=0) d = nonZeroValue;

            stats.addValue(Math.log(d) / Math.log(2));
        }

        double median = stats.getPercentile(50);

        return Math.pow(2, median);

    }

    */

}

