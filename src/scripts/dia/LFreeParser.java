package scripts.dia;

import edu.scripps.pms.util.StringUtil;
import rpark.statistics.CommonStat;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by rpark on 6/12/18.
 */
public class LFreeParser {

  public static void main(String[] args) throws Exception {


    labelfreeIntensityParser("");

    //labelfreeSpecCountParser("");

  }

  public static void labelfreeSpecCountParser(String file) throws Exception {
    if ("".equals(file))
      file = "/ip2_data4/rpark/libsearch/quant/Lib_based_quant/ident_compare_scdda.txt";

    //  file = "/ip2_data4/rpark/libsearch/quant/DB_based_quant/lfree_DIA_spec_count/ident_compare_sc9635.txt";
      //file = "/ip2_data4/rpark/libsearch/quant/DB_based_quant/lfree_DDA_spec_count/ident_compare_sc9634.txt";

    String line;

    BufferedReader br = new BufferedReader(new FileReader(file));

    while (null != (line = br.readLine())) {

      if(!line.startsWith("P\t")) continue;

      String[] arr = line.split("\t");
      String desc = arr[8];

      String gene = StringUtil.getGene(desc);

      List<Double> scList = new ArrayList<>();
      boolean zero=false;

      for(int i=0;i<6;i+=2) {

        double value =Double.parseDouble(arr[2+i]);
        if(value<=0) zero=true;
        scList.add(value);
      }

      if(zero) continue;


      System.out.print(gene + "\t" + arr[2] + "\t" + arr[4] + "\t" + arr[6] + "\t");
      System.out.println( CommonStat.getCoefficientOfVariance(scList) + "\t");



    }

    br.close();

  }

  public static void labelfreeIntensityParser(String file) throws Exception {

    if ("".equals(file))
      file = "/ip2_data4/rpark/libsearch/quant/Lib_based_quant/lfree_DDA/census_labelfree_out_100_stat.txt";

      //file = "/ip2_data4/rpark/libsearch/quant/DB_based_quant/lfree_DDA/census_labelfree_out_100_stat.txt";


    String line;

    BufferedReader br = new BufferedReader(new FileReader(file));

    int intCol = 14;
    int specCol = 11;
    while (null != (line = br.readLine())) {

      if (!line.startsWith("P\t"))
        continue;

      String[] arr = line.split("\t");

      List<Double> intList = new ArrayList<>();

      boolean specCountTwo=false;
      for (int i = 0; i < 3; i++) {
        int intValue = Integer.parseInt(arr[specCol + i]);

        //for db search //if(intValue<2) specCountTwo=true;
        //for lib search
        if(intValue<30) specCountTwo=true;

      }
      if(specCountTwo) continue;


      boolean zeroInt = false;


      for (int i = 0; i < 3; i++) {
        double intValue = Double.parseDouble(arr[intCol + i]);
        if (intValue <= 0) zeroInt = true;

        intList.add(intValue);

      }

      if (zeroInt) continue;

      //System.out.println(arr[col]);
      String desc = arr[arr.length - 1];
      String gene = StringUtil.getGene(desc);

      boolean extreamRatio=false;

      for (double d : intList) {
        double ratio = Math.log(d)/Math.log(2);
        if(ratio>=40) extreamRatio=true;

      }

      if(extreamRatio) continue;


      System.out.print(gene + "\t"); // +  + "\t");


      for (double d : intList) {
        System.out.print(d + "\t");
      }

      for (double d : intList) {
        double ratio = Math.log(d)/Math.log(2);
        System.out.print(ratio + "\t");
      }



      System.out.println(CommonStat.getCoefficientOfVariance(intList));


    }


    br.close();
  }
}
