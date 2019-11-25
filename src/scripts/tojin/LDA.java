package scripts.tojin;

import java.util.*;
import java.io.*;


import java.io.BufferedReader;
import java.util.*;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.FileReader;
import java.io.RandomAccessFile;

import edu.scripps.pms.census.util.dtaselect.Protein;
import edu.scripps.pms.census.util.dtaselect.Peptide;

import edu.scripps.pms.census.model.*;

import org.jdom.*;
import org.jdom.output.*;
import org.jdom.input.*;

public class LDA 
{
    public static void run(double[] linear_coefficients, double con) throws Exception {
	String fname = "/home/rpark/rpark_on_data/project/tojin/LDA/0421_jin_peptide_gray.txt";

	BufferedReader br = new BufferedReader(new FileReader(fname));
	String each;

	while( null != (each = br.readLine()) )
	{
	    String[] arr = each.split("\t");

	    double scan = Double.parseDouble(arr[2]);
	    double kd = Double.parseDouble(arr[3]);

	    double output = scan*linear_coefficients[0] + kd*linear_coefficients[1] + con;

	    System.out.println(each + "\t" + output);

	}
    }


    public static void main(String args[]) throws Exception
    {
	///home/rpark/rpark_on_data/project/tojin/LDA/LDA_data.txt

//	String fname = "/home/rpark/rpark_on_data/project/tojin/LDA/LDA_data.txt";
	//String fname = "/home/rpark/rpark_on_data/project/tojin/LDA/new_training_data.txt";
	String fname = "/home/rpark/rpark_on_data/project/tojin/LDA/more_data/train_New_sol_5pm_top_peptide_0424.txt";

	BufferedReader br = new BufferedReader(new FileReader(fname));

	int count=0;
	String each;
	while( null != (each = br.readLine()) )
	    count++;

	double[][] data_points = new double[2][count];
	double projections[] = new double[count];
	double linear_coefficients[] = new double[2];
	int N_dim = 2;
	int N_points = count;
	byte[] positive = new byte[count];
	br = new BufferedReader(new FileReader(fname));
	int index=0;
	while( null != (each = br.readLine()) )
	{
	    String[] arr = each.split("\t");
	    data_points[0][index] = Double.parseDouble(arr[0]);
	    data_points[1][index] = Double.parseDouble(arr[1]);
	    positive[index] = Byte.parseByte(arr[2]); 

	    index++;

	}

	double dout = edu.scripps.pms.stats.NRUtils.LDA (data_points, projections, linear_coefficients, N_dim, N_points, positive);

//	for(double d:projections)
//	    System.out.println("project:" + d);

	for(double d:linear_coefficients)
	    System.out.println("linear corr:" + d);

	double con = Double.parseDouble(args[0]);

	br = new BufferedReader(new FileReader(fname));

	int posCount=0;
	int negCount=0;
	while( null != (each = br.readLine()) )
	{
	    String[] arr = each.split("\t");

	    double x = Double.parseDouble(arr[0]);
	    double y = Double.parseDouble(arr[1]);

	    double output = x*linear_coefficients[0] + y*linear_coefficients[1] + con;

	    double posneg = Double.parseDouble(arr[2]);
//	    System.out.println(x + "\t" + y + " " + output);

	    if(output>0) {
		if(posneg==1)
		    posCount++;
		else if(posneg==0)
		    negCount++;
	    }
//	    else if(output<0 && posneg==0)

	}

	System.out.println("pos: " + posCount);
	System.out.println("neg: " + negCount);
	System.out.println("FP: " + ((double)negCount/posCount));

//	double test = Double.parseDouble(args[1])*linear_coefficients[0] + Double.parseDouble(args[2])*linear_coefficients[1] + con;

//	System.out.println("test : " + test);
	System.out.println("eq: y=-" + (linear_coefficients[0]/linear_coefficients[1]) + " x - " + (con/linear_coefficients[1]));
	System.out.println("eq: " + (linear_coefficients[0]) + " x + " + (linear_coefficients[1]) + " y + " + con);

	//run(linear_coefficients, con);
    }

}

