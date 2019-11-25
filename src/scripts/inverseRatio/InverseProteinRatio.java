import java.util.*;
import java.io.*;
import java.text.*;

public class InverseProteinRatio {
    
    public static void main(String[] args) throws Exception {
	DecimalFormat formatter = new DecimalFormat("0.00");

	BufferedReader br = new BufferedReader(new FileReader(args[0]));

	String eachLine;

	while( (eachLine = br.readLine()) != null ) {

	    if( !eachLine.startsWith("P\t") ) {
	    
		System.out.println(eachLine);
		continue;

	    }

	    String[] arr = eachLine.split("\t");

	    try {
		double d = Double.parseDouble(arr[2]);

		for(int i=0;i<2;i++) {
		    System.out.print(arr[i]);
		    System.out.print("\t");
		}

		String inv = formatter.format(1/d);
		System.out.print(inv);
		System.out.print("\t");


		for(int i=3;i<arr.length;i++) {
		    System.out.print(arr[i]);
		    System.out.print("\t");
		}

		System.out.print("\n");
	    } catch (Exception e) {

		System.out.print(eachLine);

	    }

	}


    }
}
