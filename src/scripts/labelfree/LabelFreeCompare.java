package scripts.labelfree;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * Created by rpark on 8/10/18.
 */
public class LabelFreeCompare {

  public static void main(String[] args) throws Exception {


   // String filename = "/ip2_data4/rpark/libsearch/quant/Lib_based_quant/lfree_DDA/re-run/census_labelfree_out_filled.txt";
    String filename = args[0];

   // System.out.println("===" + filename);

    //labelfreeFilledFileComparePeptide(filename);
    labelfreeFilledFileCompareProtein(filename);

  }

  //compare in peptide level
  public static void labelfreeFilledFileCompareProtein(String filename) throws Exception {

    BufferedReader br = new BufferedReader(new FileReader(filename));

    String eachLine;
    while (null != (eachLine = br.readLine()) && !eachLine.startsWith("PLINE"));

    String[] arr = eachLine.split("\t");
    List<Integer> intensityIndexList = new ArrayList<>();
    List<Integer> normIntensityIndexList = new ArrayList<>();

    for(int i=0;i<arr.length;i++) {
      if(arr[i].startsWith("INTENSITY_"))
        intensityIndexList.add(i);
      else if(arr[i].startsWith("NORM_INTENSITY_"))
        normIntensityIndexList.add(i);

    }


    //System.out.println(intensityIndexList);
    //System.out.println(normIntensityIndexList);



    HashMap<String, Double> hashMap = new HashMap<>();
    while (null != (eachLine = br.readLine())) {

      if (!eachLine.startsWith("P\t")) continue;

      String[] tmparr = eachLine.split("\t");
      String accession = tmparr[1];
      String desc = tmparr[tmparr.length-1];


      System.out.print(accession);

      for(int each:intensityIndexList)
        System.out.print("\t" + tmparr[each]);

      for(int each:normIntensityIndexList)
        System.out.print("\t" + tmparr[each]);

      System.out.println();

      //hashMap.put(arr[2] + arr[arr.length-1], Double.parseDouble(arr[3]));
    }

    br.close();

    // return hashMap;


  }

  //compare in peptide level
  public static void labelfreeFilledFileComparePeptide(String filename) throws Exception {

    BufferedReader br = new BufferedReader(new FileReader(filename));

    String eachLine;
    while (null != (eachLine = br.readLine()) && !eachLine.startsWith("SLINE"));

    String[] arr = eachLine.split("\t");
    List<Integer> intensityIndexList = new ArrayList<>();
    List<Integer> normIntensityIndexList = new ArrayList<>();

    for(int i=0;i<arr.length;i++) {
      if(arr[i].equals("INTENSITY"))
        intensityIndexList.add(i);
      else if(arr[i].startsWith("NORM_INTENSITY"))
        normIntensityIndexList.add(i);

    }


    //System.out.println(intensityIndexList);
    //System.out.println(normIntensityIndexList);



    HashMap<String, Double> hashMap = new HashMap<>();
    while (null != (eachLine = br.readLine())) {

      if (!eachLine.startsWith("S\t")) continue;

      String[] tmparr = eachLine.split("\t");
      String sequence = tmparr[2];
      String cs = tmparr[5];

      boolean tooBigSmall=false;
      for(int each:intensityIndexList) {
        double d = Double.parseDouble(tmparr[each]);
/*
        //if(d<=10000 || d>=1000000000) {
        if(d>=1000000000) {
          tooBigSmall = true;
          System.out.print("x");
         // System.out.println("\t" + tmparr[each]);
        } else {
          System.out.print("0");

        }
        */
        // System.out.print("..\t" + tmparr[each]);


      }


      //  if(tooBigSmall) continue;

//if(true) continue;

      System.out.print(sequence + cs);

      for(int each:intensityIndexList)
        System.out.print("\t" + tmparr[each]);

      for(int each:normIntensityIndexList)
        System.out.print("\t" + tmparr[each]);

      System.out.println();

      //hashMap.put(arr[2] + arr[arr.length-1], Double.parseDouble(arr[3]));
    }

    br.close();

    // return hashMap;


  }

}
