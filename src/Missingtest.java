
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.distribution.NormalDistribution;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author rampuria
 */
public class Missingtest {
    
    
    public static void main(String [] args) throws FileNotFoundException, IOException{
         BufferedReader br = new BufferedReader(new FileReader("/home/rpark/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                census_labelfree_out_12575_filled.txt-"));
         String eachLine = null;
         List<Integer> index = new ArrayList<>();
         List<Double> sumInt = new ArrayList<>();
         int count=0;
         while((eachLine=br.readLine()) != null){
             if(eachLine.startsWith("H\t")){
                 continue;
             }
             else if(eachLine.startsWith("SLINE\t")){
                 String [] words = eachLine.split("\t");
                 
                 System.out.println("");
                 for(int i =0;i<words.length;i++){
                     if(words[i].startsWith("INTENSITY")){
                         index.add(i);
                     }
                 }
             }
             else if(eachLine.startsWith("S\t")){
                 count++;
                 String [] words = eachLine.split("\t");
                 double sum =0;
                 for(int i:index){
                     sum = sum+Double.parseDouble(words[i]);
                 }
                 sumInt.add(sum);
             }
         }
         BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream("/home/rpark/slinecount.txt"));
        PrintStream ps = new PrintStream(out);


         
         DescriptiveStatistics stat = new DescriptiveStatistics();
         double [] d = new double[sumInt.size()];
         for(int i=0;i<sumInt.size();i++){
             stat.addValue(Math.log10(sumInt.get(i)));
             d[i] =Math.log10(sumInt.get(i));
         }
         double avg = stat.getMean();
         double stdev = stat.getStandardDeviation();
         NormalDistribution ndist = new NormalDistribution(avg, stdev);
         Arrays.sort(d);
         for(double l :d){
             double pval= ndist.cumulativeProbability(l); //1-ndist.cumulativeProbability(l);
             ps.println(l+"\t"+pval);
         }
         
         ps.close();
         out.close();
    }
    
}
