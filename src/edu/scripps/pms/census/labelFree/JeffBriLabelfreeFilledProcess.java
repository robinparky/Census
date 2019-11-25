package edu.scripps.pms.census.labelFree;

import edu.scripps.pms.census.conf.Configuration;
import edu.scripps.pms.census.labelFree.model.*;
import edu.scripps.pms.census.model.ChroPeptide;
import edu.scripps.pms.census.labelFree.model.LabelfreePeptide;
import edu.scripps.pms.census.model.ChroProtein;
import edu.scripps.pms.stats.TTest;
import rpark.statistics.AnovaUtil;
import rpark.statistics.BHCorrection;
import rpark.statistics.Outlier;
import rpark.statistics.TTestUtil;

import java.io.*;
import java.util.*;

/**
 * Created by Titus Jung on 2/13/17.
 */
public class JeffBriLabelfreeFilledProcess {

    public static void main(String [] args) throws Exception
    {

        String fname = args[0];
        String filename = args[1];



        LabelfreeFilledParser filledParser = new LabelfreeFilledParser(fname);
        List<ProteinModel> proteinList = filledParser.readWholeFile();
        List<double []> ratioList = new ArrayList<>();
        //List<List<double[]>> ratioList = new ArrayList<>();
        Map<String,String> locusRatioMap = new HashMap<>();
        for(ProteinModel protein: proteinList)
        {
            List<LabelfreePeptide> peptideList = protein.getPeptideList();
            //List<double[]> ratios = new ArrayList<>();
            List<Double> ratios = new ArrayList<>();
            for(LabelfreePeptide peptide: peptideList)
            {
                List<ChroPeptide> chroPeptideList = peptide.getPeptideList();
                int i=0;
                double [] intensityArr = new double[8];
                Double [] ratioArr = new Double[4];
                for(ChroPeptide chroPeptide: chroPeptideList)
                {

                    intensityArr[i++]=chroPeptide.getAverageIntensity();
                }
                ratioArr[0] = intensityArr[0]/intensityArr[4];
                ratioArr[1] = intensityArr[1]/intensityArr[5];
                ratioArr[2] = intensityArr[2]/intensityArr[6];
                ratioArr[3] = intensityArr[3]/intensityArr[7];
                ratios.addAll(Arrays.asList(ratioArr));
            }
            //List<double []> outlierList = Outlier.removeOutlierDoubleArr(ratios,0.01);
            List<Double> outlierList = Outlier.removeOutlierZero(ratios,0.01);

            double [] input = new double[outlierList.size()];
            for(int i=0; i<input.length;i++)
            {
                input[i] = outlierList.get(i);
            }
            ratioList.add(input);



            //ratioList.add(outlierList);
            for(int i=0; i<protein.getRedundnatProteinList().size();i++)
            {
                locusRatioMap.put(protein.getRedundnatProteinList().get(i).getLocus(),createRatioString(outlierList));
            }
        }

        SortedMap<String,double[]> ratioMap = new TreeMap<>();
        SortedMap<String,Double> ratioPValueMap = new TreeMap<>();
        SortedMap<String,Double> ratioQValueMap = new TreeMap<>();

        List<Double> pvalueList = new ArrayList<>();
        for(double[] dlist:ratioList)
        {
            ratioMap.put(createRatioString(dlist),dlist);
        }
        for(String key: ratioMap.keySet())
        {
            if(ratioMap.get(key).length<2)
            {
            }
            else
            {
                double pvalue = TTestUtil.oneSampleTTestConvertToLog(ratioMap.get(key));
                //pvalueList.add(pvalue);
                ratioPValueMap.put(key,pvalue);
            }
        }
        for(double d:ratioPValueMap.values())
        {
            pvalueList.add(d);
        }
        /*
        for(double [] arr: ratioMap.values())
        {
            double pvalue = TTestUtil.oneSampleTTestConvertToLog(arr);
            pvalueList.add(pvalue);
        }*/
        int i=0;



        List<Double> qvalueList = BHCorrection.runBhCorrection(pvalueList);

        for(String key: ratioPValueMap.keySet())
        {
            ratioQValueMap.put(key,qvalueList.get(i++));
        }
        BufferedWriter bw = new BufferedWriter(new FileWriter(filename));

        bw.append("PLINE\tACCESSION\tDESCRIPTION\tRATIO_STRING\tP_VALUE\tQ_VALUE");
        bw.newLine();
        for (Iterator<ProteinModel> itr = proteinList.iterator(); itr.hasNext();) {
            ProteinModel protein = itr.next();
            for (ChroProtein cprotein : protein.getRedundnatProteinList()) {
                bw.append("P\t");
                bw.append(cprotein.getLocus());
                bw.append("\t");
                bw.append(cprotein.getDescription());
                bw.append("\t");
                String ratio = locusRatioMap.get(cprotein.getLocus());
                bw.append(ratio);
                bw.append("\t");
                if(ratioPValueMap.get(ratio)!=null)
                {
                    bw.append(ratioPValueMap.get(ratio).toString());
                }
                else
                {
                    bw.append("Na");
                }
                bw.append("\t");
                if(ratioQValueMap.get(ratio)!=null)
                {
                    bw.append(ratioQValueMap.get(ratio).toString());
                }
                else
                {
                    bw.append("Na");
                }
            }
            bw.flush();
            bw.newLine();
        }
        bw.close();

        System.out.println("output in :"+filename);
    }
    public static String createRatioString(List<Double> list)
    {
        StringBuilder sb = new StringBuilder();
        for(double d: list)
        {
            sb.append(Double.toString(d)+",");
        }
        return sb.toString();
    }

    public static String createRatioString(double[] list)
    {
        StringBuilder sb = new StringBuilder();
        for(double d: list)
        {
            sb.append(Double.toString(d)+",");
        }
        return sb.toString();
    }

    public static String createRatioStringArr(List<double[]> list)
    {
        StringBuilder sb = new StringBuilder();
        for(double [] darr: list)
        {
            for(double d: darr)
                sb.append(Double.toString(d)+",");
        }
        return sb.toString();
    }
}
