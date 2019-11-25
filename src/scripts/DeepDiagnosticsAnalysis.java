package scripts;

import rpark.statistics.TTestUtil;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;

public class DeepDiagnosticsAnalysis {

    public static void main(String[] args) throws Exception {


        //String inputFile = "/home/rpark/test_data/deep_diagonastic_solutions/labelfree_18910/pep_AI.txt";
        String inputFile = "/home/rpark/test_data/deep_diagonastic_solutions/labelfree_18910/pep_AI_norm.txt";
        BufferedReader br = new BufferedReader(new FileReader(inputFile));
        String eachLine;

        while( null != (eachLine=br.readLine()) ) {

            String[] arr = eachLine.split("\t");

          //  System.out.println(arr.length);

            double[] group1 = new double[18];
            double[] group2 = new double[18];

            List<Double> list1 = new ArrayList<>();
            List<Double> list2 = new ArrayList<>();

            double sum = 0;

            for(int i=3;i<21;i++) {
                list1.add(Math.log10(Double.parseDouble(arr[i])));
                sum += Math.log10(Double.parseDouble(arr[i]));
            }

            for(int i=21;i<arr.length;i++) {
                list2.add(Math.log10(Double.parseDouble(arr[i])));
                sum += Math.log10(Double.parseDouble(arr[i]));
            }



            for(int i=0;i<group1.length;i++)
                group1[i] = list1.get(i);

            for(int i=0;i<group2.length;i++)
                group2[i] = list2.get(i);

       //     System.out.println(group1.length + " " + group2.length);
         //   System.out.println(list1.size() +  " " + list2.size());



            double pvalue = TTestUtil.calculateUnpairedTTest(group1, group2);

            System.out.print(pvalue + "\t");
            for(int i=3;i<21;i++) {
                System.out.print( Math.log10(Double.parseDouble(arr[i])) / sum + "\t");
            }

            for(int i=21;i<arr.length;i++) {
                System.out.print( Math.log10(Double.parseDouble(arr[i])) / sum + "\t");
            }

            System.out.println("");

            //System.out.println(pvalue + "\t" + eachLine);
        }



    }
}
