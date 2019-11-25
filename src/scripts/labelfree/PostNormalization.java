/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package scripts.labelfree;

/**
 *
 * @author rpark
 */
import java.io.*;

public class PostNormalization {

    public static void main(String[] args) throws Exception {

//        args[0] = "/data/2/rpark/ip2_data/pankows/CFTR/labelfree_quant/census_labelfree_out_1123.txt";
//        args[1] = "IPI00302383.3";
/*
        String args1 = "/data/2/rpark/ip2_data/xmhan/FAd_labelfree/labelfree_quant/census_labelfree_out_1596.txt";
        String args2 = "IPI00220030.1";
        PostNormalization pn = new PostNormalization();
        pn.normalizeWithProtein(args1, args2);
       */ 
        String args1 = "/data/2/rpark/ip2_data/xmhan/FAd_labelfree/labelfree_quant/census_labelfree_out_1596_stat.txt";
        String args2 = "IPI00220030.1";
        PostNormalization pn = new PostNormalization();
        pn.normalizeWithProteinSpecCount(args1, args2);
    }

    public void normalizeWithProteinSpecCount(String filename, String locus) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(filename));

        String eachLine; // = br.readLine();
        String searchKey = "P\t" + locus;

        while( (eachLine = br.readLine())!=null && !eachLine.startsWith("H\tNormalization") );

        int expSize = eachLine.split("\t").length -2;

        if(expSize<=0) {
            System.out.println("Error: experiment size is zero");
            throw new IOException();
        }

        while( (eachLine = br.readLine())!=null && !eachLine.startsWith("PLINE\t") );
	String[] arr = eachLine.split("\t");
	int[] intensityIndexArr = new int[expSize];

	int tmpIndex=0;
	for(int i=0;i<arr.length;i++) {
		String s = arr[i];
		//if(!s.startsWith("NORM_INTENSITY_") && s.startsWith("NORM_INTENSITY_CORRECT_"))
		if(s.startsWith("NORM_INTENSITY_CORRECT_"))
			intensityIndexArr[tmpIndex++] = i;
	}

//	for(int i:intensityIndexArr) 
//		System.out.println(arr[i]);

	while( (eachLine = br.readLine())!=null )
	{
            if(!eachLine.startsWith("P\t"))
                continue;

            if(eachLine.startsWith(searchKey)) {
                //System.out.println(eachLine);
                break;
            }
        }

	double[] normIntArr = new double[expSize];
	String[] targetProteinArr = eachLine.split("\t");
	double intSum = 0;
	for(int i:intensityIndexArr) {
		double normInt = Double.parseDouble(targetProteinArr[i]);
		intSum += normInt;

	}
	double avgInt = intSum/expSize;
		System.out.println(avgInt + " " + intSum);
    }


    public void normalizeWithProtein(String filename, String locus) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(filename));

        String eachLine; // = br.readLine();
        String searchKey = "P\t" + locus;

        while( (eachLine = br.readLine())!=null && !eachLine.startsWith("H\tNormalization") );

        int expSize = eachLine.split("\t").length -2;

        if(expSize<=0) {
            System.out.println("Error: experiment size is zero");
            throw new IOException();
        }

        while( (eachLine = br.readLine())!=null && !eachLine.startsWith("PLINE\t") );
	String[] arr = eachLine.split("\t");
	int[] intensityIndexArr = new int[expSize];

	int tmpIndex=0;
	for(int i=0;i<arr.length;i++) {
		String s = arr[i];
		//if(!s.startsWith("NORM_INTENSITY_") && s.startsWith("NORM_INTENSITY_CORRECT_"))
		if(s.startsWith("NORM_INTENSITY_CORRECT_"))
			intensityIndexArr[tmpIndex++] = i;
	}

//	for(int i:intensityIndexArr) 
//		System.out.println(arr[i]);

	while( (eachLine = br.readLine())!=null )
	{
            if(!eachLine.startsWith("P\t"))
                continue;

            if(eachLine.startsWith(searchKey)) {
                //System.out.println(eachLine);
                break;
            }
        }

	double[] normIntArr = new double[expSize];
	String[] targetProteinArr = eachLine.split("\t");
	double intSum = 0;
	for(int i:intensityIndexArr) {
		double normInt = Double.parseDouble(targetProteinArr[i]);
		intSum += normInt;

	}
	double avgInt = intSum/expSize;
		System.out.println(avgInt + " " + intSum);
    }

}
