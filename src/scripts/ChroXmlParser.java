package scripts;

import org.jdom.input.SAXBuilder;
import org.jdom.*;
import rpark.statistics.CommonStat;

import java.io.File;
import java.util.*;

public class ChroXmlParser
{


    public static void main(String args[]) throws Exception
    {
      //db search
      /*
        String file1 = "/ip2_data4/rpark/libsearch/quant/DB_based_quant/DIA1/census_chro.xml";
        String file2 = "/ip2_data4/rpark/libsearch/quant/DB_based_quant/DIA2/census_chro.xml";
        String file3 = "/ip2_data4/rpark/libsearch/quant/DB_based_quant/DIA3/census_chro.xml";
*/

       //lib search

        String file1 = "/ip2_data4/rpark/libsearch/quant/Lib_based_quant/DIA1/census_chro.xml";
        String file2 = "/ip2_data4/rpark/libsearch/quant/Lib_based_quant/DIA2/census_chro.xml";
        String file3 = "/ip2_data4/rpark/libsearch/quant/Lib_based_quant/DIA3/census_chro.xml";



        Map<String, Double> geneMap1 = getGeneIntMap(file1);
        Map<String, Double> geneMap2 = getGeneIntMap(file2);
        Map<String, Double> geneMap3 = getGeneIntMap(file3);



      for(Iterator<String> itr=geneMap1.keySet().iterator(); itr.hasNext(); ) {
        String key = itr.next();
        Double intensity1 = geneMap1.get(key);
        Double intensity2 = geneMap2.get(key);
        Double intensity3 = geneMap3.get(key);


        if(null != intensity2 && null != intensity3) {
          List<Double> intList = new ArrayList<>();
          intList.add(intensity1);
          intList.add(intensity2);
          intList.add(intensity3);

          double rsd = CommonStat.getCoefficientOfVariance(intList);
          //double rsd = CommonStat.getRelativeStandardDeviation(intList);

          System.out.println(key + "\t" + intensity1 + "\t" + intensity2 + "\t" + intensity3 + "\t" +
            Math.log(intensity1)/Math.log(2) + "\t" + Math.log(intensity2)/Math.log(2) + "\t" + Math.log(intensity3)/Math.log(2) + "\t" +
            rsd); // + "\t" + intList);

        //  return;
        }


      }



         //System.out.println(geneMap1.size() + " " + geneMap2.size() + " " + geneMap3.size());
    }


  public static Map<String, Double> getGeneIntMap(String file) throws Exception {


    SAXBuilder sb = new SAXBuilder();
    Document doc = sb.build(new File(file));

    Element rootEle = doc.getRootElement();
    List<Element> proList = rootEle.getChildren("protein");
    Map<String, Double> proteinMap = new HashMap<>();

    for(Iterator<Element> itr=proList.iterator(); itr.hasNext(); ) {
      Element proEle = itr.next();
      List<Element> pepList = proEle.getChildren("peptide");

      List<Double> pepIntList = new ArrayList<>();

      for(Iterator<Element> pepItr=pepList.iterator(); pepItr.hasNext(); ) {

        Element pepEle = pepItr.next();
        String chro = pepEle.getChildText("chro");

        String[] chroArr = chro.split(";");

        //  System.out.println(chroArr[0]);

        String[] arr = chroArr[0].split(" ");
        int startScan = Integer.parseInt(arr[1]);
        int endScan = Integer.parseInt(arr[2]);

        //if(true) return;

        double intSum = 0;
        for(int i=1;i<chroArr.length;i++) {

          String[] eachArr = chroArr[i].split(" ");

          int curScan = Integer.parseInt(eachArr[0]);

          if(curScan<startScan || curScan>endScan) continue;


          //  System.out.println(eachArr[0] + "\t" + eachArr[1] + "\t" + eachArr[2]);
          intSum += Double.parseDouble(eachArr[1]);


        }

        pepIntList.add(intSum);
        //System.out.println(intSum);
      }

      //  System.out.println(pepIntList);
      //  System.out.println(pepIntList.size());
      double mInt = CommonStat.getMedianValue(pepIntList);

      String locus = proEle.getAttributeValue("locus");
      String desc = proEle.getAttributeValue("desc");
      String gene = desc.substring(desc.indexOf("GN=")+3);
      gene = gene.substring(0, gene.indexOf(" "));

      proteinMap.put(gene, mInt);
      // System.out.println(locus + "\t>" + gene + "<\t" + mInt);

      //        P 42487.0 42907.0 66.288 66.788;42508 66.314 0 -1.0 0.0;42529 66.338 0 -1.0 0.0;42550 66.364 316742 -1.0 0.0;42571 66.39 633485 -1.0 0.0;42592 66.413 950228

    }


    return proteinMap;


  }
}

