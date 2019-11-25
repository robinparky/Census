import java.io.*;
import edu.scripps.pms.census.hash.IndexedFile;
import edu.scripps.pms.census.util.*;

public class SRMTest {
	public static void main(String[] args) throws Exception {

		File indexFile = new File("kiran/2009-0722-31-51-AB-TMT-O2-s01.ms2.index");
		IndexedFile iFile = new IndexedFile(indexFile, "kiran/2009-0722-31-51-AB-TMT-O2-s01.ms2");
		double[][] arr = CalcUtilGeneric.getSpectrumArr(iFile, 3510);

		for(int i=0;i<arr.length;i++) 
			for(int j=0;j<arr[i].length;j++)
				System.out.println(arr[i][j]);

	}

}
