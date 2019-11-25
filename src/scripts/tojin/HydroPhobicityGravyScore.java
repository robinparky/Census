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

public class HydroPhobicityGravyScore 
{
    public static void main(String args[]) throws Exception
    {
	//http://gcat.davidson.edu/DGPB/kd/aminoacidscores.htm
	BufferedReader br = new BufferedReader(new FileReader(args[0]));
	String each;

	while( null != (each = br.readLine()) )
	{
	    String[] arr = each.split("\t");
	    System.out.print(each);
	    System.out.print("\t");

	    double score = getGravy(arr[4]);
	    System.out.print(score);
	    System.out.print("\t");
	    System.out.print(arr[4].toCharArray().length);
	    System.out.print("\t");
	    System.out.println(score/arr[4].toCharArray().length);
	}
    }

    public static double getGravy(String seq) {
	char[] ch = seq.toCharArray();

	double score=0;

	for(char c:ch) {
	    switch(c) {
		case 'I': score += 4.5; break;
		case 'V': score += 4.2; break;
		case 'L': score += 3.8; break;
		case 'F': score += 2.8; break;
		case 'C': score += 2.5; break;
		case 'M': score += 1.9; break;
		case 'A': score += 1.8; break;
		case 'G': score += -0.4; break;
		case 'T': score += -0.7; break;
		case 'W': score += -0.9; break;
		case 'S': score += -0.8; break;
		case 'Y': score += -1.3; break;
		case 'P': score += -1.6; break;
		case 'H': score += -3.2; break;
		case 'E': score += -3.5; break;
		case 'Q': score += -3.5; break;
		case 'D': score += -3.5; break;
		case 'N': score += -3.5; break;
		case 'K': score += -3.9; break;
		case 'R': score += -4.5; break;

		default :
		    System.out.println("invalid residue " + c);
		    System.exit(0);
		  
		  /*
		Valine 	V 	4.2
		Leucine 	L 	3.8
		Phenylalanine 	F 	2.8
		Cysteine 	C 	2.5
		Methionine 	M 	1.9
		Alanine 	A 	1.8
		Glycine 	G 	-0.4
		Threonine 	T 	-0.7
		Tryptophan 	W 	-0.9
		Serine 	S 	-0.8
		Tyrosine 	Y 	-1.3
		Proline 	P 	-1.6
		Histidine 	H 	-3.2
		Glutamic acid 	E 	-3.5
		Glutamine 	Q 	-3.5
		Aspartic acid 	D 	-3.5
		Asparagine 	N 	-3.5
		Lysine 	K 	-3.9
		Arginine 	R 	-4.5
		*/
		    
		break;


	    }
	}


	return score;
    
    }
}

