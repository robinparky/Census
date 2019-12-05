package scripts;

import edu.scripps.pms.util.io.CensusOutputReader;

import org.apache.commons.math3.stat.regression.SimpleRegression;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.nio.Buffer;
import java.util.HashMap;

/**
 * Created by rpark on 7/19/18.
 */
public class DIALibCompare {



  public static void main(String[] args) throws Exception {



    //String file1 = "/ip2_data4/rpark/libsearch/quant/DB_based_quant/DIA1/census-out.txt";
    //String file2 = "/ip2_data4/rpark/libsearch/quant/DB_based_quant/DIA2/census-out.txt";


    String file1 = args[0];
    String file2 = args[1];

    DIALibCompare compare = new DIALibCompare();
    HashMap<String, Double> map1 = compare.getPeptideAbundance(file1);
    HashMap<String, Double> map2 = compare.getPeptideAbundance(file2);

//    HashMap<String, Double> map1 = compare.getProteinAbundance(file1);
 //   HashMap<String, Double> map2 = compare.getProteinAbundance(file2);


 //   System.out.println(map1.size() + " " + map2.size());

    SimpleRegression regression = new SimpleRegression();
    int count=0;
    for(String key:map1.keySet()) {
      Double v1 = map1.get(key);
      Double v2 = map2.get(key);

      if(v1 != null && v2 != null) {
        double d1 = v1.doubleValue();
        double d2 = v2.doubleValue();
        if(d1 <=0 || d2<=0) continue;

        regression.addData(v1.doubleValue(), v2.doubleValue());
        count++;
       System.out.println(key + "\t" + v1 + "\t" + v2);
      }

    }

    System.out.println("data size:\t" + count);
    System.out.println("Regression:\t" + regression.getR());


  }

  public HashMap getPeptideAbundance(String file) throws Exception {
    BufferedReader br = new BufferedReader(new FileReader(file));

    String eachLine;
    HashMap<String, Double> hashMap = new HashMap<>();
    while (null != (eachLine = br.readLine())) {

      if (!eachLine.startsWith("S\t")) continue;

      String[] arr = eachLine.split("\t");
      hashMap.put(arr[2] + arr[arr.length-1], Double.parseDouble(arr[3]));


    }


    br.close();

    return hashMap;
  }


  public HashMap getProteinAbundance(String file) throws Exception {
    BufferedReader br = new BufferedReader(new FileReader(file));

    String eachLine;
    HashMap<String, Double> hashMap = new HashMap<>();
    while( null != (eachLine=br.readLine()) ) {

      if(!eachLine.startsWith("P\t")) continue;

      String[] arr = eachLine.split("\t");
      hashMap.put(arr[1], Double.parseDouble(arr[4]));


    }

   // System.out.println(hashMap);

    br.close();

    return hashMap;
  }

}
