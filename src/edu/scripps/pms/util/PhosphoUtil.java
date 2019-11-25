package edu.scripps.pms.util;

/**
 * Created by rpark on 3/20/18.
 */
public class PhosphoUtil {

  public static void main(String[] args) {


    String seq = "(12.3)ABCD(67.6)XYZ(23.44)FDSA(22.3)";
    double ptmMass = getTotalPTMMass(seq);

    System.out.println(cleanSequence(seq));
    System.out.println(ptmMass);
  }


  public static String cleanSequence(String seq) {
    seq = removeNumerics(seq);
    seq = seq.replaceAll("\\(\\)", "");
    seq = seq.replaceAll("\\(-\\)", "");
    //seq = seq.replaceAll("\\(\\)", "\\*");
    seq = seq.replaceAll("\\*", "");
    seq = seq.replaceAll("#", "");
    seq = seq.replaceAll("@", "");

    return seq;
  }
  public static String removeNumerics(String str) {
    if (str == null) {
      return null;
    }else {
      return str.replaceAll("[0-9|.]", "");
    }
  }

  public static double getTotalPTMMass(String str) {

    if (str == null) {
      return 0;
    }

    StringBuffer strBuff = new StringBuffer();
    char c;
    double totalPtmMass=0;

    for (int i = 0; i < str.length() ; i++) {
      c = str.charAt(i);

      if(c == '(') {

        String num ="";

        c = str.charAt(++i);
        while(c != ')') {

          num += c;
          c = str.charAt(++i);
        }

        totalPtmMass += Double.parseDouble(num);
        //System.out.println(num);

      }
    }

    return totalPtmMass;
  }

}
