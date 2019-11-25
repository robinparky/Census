/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.scripps.pms.util.stats;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.apache.commons.math3.stat.inference.TestUtils;

/**
 * To 
 *
 * @author Harshil
 */
public class OneTTest {

    public static void main(String args[]) throws Exception {
        
        StringBuffer sb = new StringBuffer();
        sb.append("How to run: java OneTTest inputfile").append("\n\n");
        sb.append("Note. input file format must be as below").append("\n");
        sb.append("1. 1st line should be header, which is not used for ttest").append("\n");
        sb.append("2. 1st column should be keys(e.g. protein accession, peptide sequence, etc.").append("\n");
        sb.append("3. from 2nd column to the 2nd last column should be numbers to be used for one sample ttest").append("\n");
        
        if(args.length<1) {
            System.out.println("Usage: java OneTTest inputfile ");
            System.out.println(sb.toString());
            System.exit(0);
        }
            

        System.out.println(sb.toString());
        
        System.out.println("running now...");
        
            
        /*
        if(args.length<=0) {
            System.out.println("Usage: java OneTTest [filePath] [list of Columns]");
            System.exit(0);
        }*/
     
      
       // args = new String[1];
       // args[0] = "/home/rpark/test_data/oneSamplet_test_EFG1-22.txt";
        
        double mu = 0d;
        String outputFile = args[0]+ ".out";
        PrintStream p = new PrintStream( new BufferedOutputStream(new FileOutputStream(outputFile)));

        
//        BufferedReader br = new BufferedReader(new FileReader(new File(args[0])));
        
        BufferedReader br = new BufferedReader(new FileReader(new File(args[0])));
       // outSb.append(br.readLine()).append("\n");
        p.println(br.readLine()); //read out header
        

        String eachLine;
        while( null !=(eachLine=br.readLine()) ) {
        
           // System.out.println(eachLine);
            String[] arr = eachLine.split("\t");
            p.print(arr[0]); p.print("\t");
        //    outSb.append(arr[0]).append("\t");
            
            int size = arr.length;
            double[] group = new double[size-2];
            
            for (int i = 1; i < size-1; i++) {
                    //currentGroup[i] = Double.parseDouble(words[Integer.parseInt(args[i + 1])-1]);
                //}
                
                p.print (arr[i] + "\t" );
          //      outSb.append(arr[i]).append("\t");
                
                group[i-1] = Double.parseDouble(arr[i]);
                
            }

            p.print(TestUtils.tTest(mu, group));p.print("\t");
            p.println(arr[arr.length-1]);
            
            //outSb.append(TestUtils.tTest(mu, group)).append("\t");
            //outSb.append(arr[arr.length-1]).append("\n");
            
            
        }

        br.close();
        
        System.out.println("See output file " +outputFile);
        
        p.close();
        //System.out.println(outSb.toString());

    }
    
    /*
        public static void main(String args[]) {

        
        args = new String[4];
        args[0] = "/cheezer_share/rpark/oneSamplet_test_EFG1-2.txt";
        args[1] = "a";
        args[1] = "2";
        args[1] = "3";
        
        BufferedReader br = null;
        try {
            double mu = 0d;
            int groupSize = args.length - 1;

            double[] currentGroup = new double[groupSize];
            br = new BufferedReader(new FileReader(new File(args[0])));
            br.readLine();
            
            String currentLine = br.readLine();
            while (currentLine != null) 
            {
                String words[] = currentLine.split("\t");
                //if (currentLine.startsWith("S"))
                {
                    for (int i = 0; i < groupSize; i++) {
                        currentGroup[i] = Double.parseDouble(words[Integer.parseInt(args[i + 1])-1]);
                    }
                    System.out.println(currentLine + "\t" +TestUtils.tTest(mu, currentGroup)) ;
                }

                currentLine = br.readLine();
            }
        } catch (Exception e) {
            StringBuffer sb = new StringBuffer();
            sb.append("Format ")
            Logger.getLogger(OneTTest.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(OneTTest.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            try {
                br.close();
            } catch (IOException ex) {
                Logger.getLogger(OneTTest.class.getName()).log(Level.SEVERE, null, ex);
            }
        }


    }
        */
        
        
}
