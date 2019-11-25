package scripts;

import edu.scripps.pms.census.labelFree.LabelfreeFilledParserTemp;
import edu.scripps.pms.census.labelFree.LabelfreeMissingPeptideBuilderSplit;
import edu.scripps.pms.census.labelFree.ProteinModel;
import edu.scripps.pms.census.util.CensusHelper;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.inference.TestUtils;
import org.jdom.input.SAXBuilder;
import org.jdom.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.*;
import edu.scripps.pms.census.model.*;
import edu.scripps.pms.util.*;
import rpark.statistics.BHCorrection;
import rpark.statistics.TTestUtil;

import edu.scripps.pms.census.util.dtaselect.Peptide;
import java.util.Arrays;
import org.apache.commons.math3.fitting.*;
import org.apache.commons.math3.fitting.WeightedObservedPoints;
public class Test
{
    public static void main(String args[]) throws Exception
    {

        /*
        0 = 19.011958814763226
1 = 19.042690715207804

0 = 19.66099517321764
1 = 19.164023036671388

         */


        DescriptiveStatistics list = new DescriptiveStatistics();
        list.addValue(6.071794734163179);
        list.addValue(5.004431069214665);
        list.addValue(4.456425311423719);
        list.addValue(5.343509483777055);
        list.addValue(5.762108965249347);
        list.addValue(6.64385619);
        list.addValue(7.644880455939909);
        list.addValue(5.76669104972917);
        list.addValue(4.947944999448616);
        list.addValue(5.1440947777568695);
        list.addValue(10.580840440009345);
        list.addValue(11.013358634276747);
        list.addValue(9.506495176168114);
        list.addValue(9.057530369395943);
        list.addValue(8.688017148420359);
        list.addValue(6.64385619);
        list.addValue(11.010070799618942);
        list.addValue(2.000355655690796);
        list.addValue(10.955989074296165);
        list.addValue(7.456423120166089);
        list.addValue(4.186268426927853);
        list.addValue(6.42454740411225);
        list.addValue(7.172152926996229);


        double dd = Math.pow(2, list.getPercentile(50));

        System.out.println(dd);


        if(true) return;





        // double[] arr1 = {19.25563603663021,  19.25563603663021, };
        double[] arr1 = {0, 0, };
        double[] arr2 = {0, 100004, };
        System.out.println("===" + TestUtils.tTest(arr1, arr2) );



        if(true) return;


        WeightedObservedPoints obs = new WeightedObservedPoints();

/*
25 31 26	==	90.054	45.0
25 31 27	==	90.098	45.0
25 31 28	==	90.141	96.0
25 31 29	==	90.184	50.0
25 31 30	==	90.227	50.0
25 31 31	==	90.27	0.0
*/
            obs.add(90.00, 0.0);
            obs.add(90.014, 0.0);
            obs.add(90.054, 45.0);
            obs.add(90.098, 45.0);
            obs.add(90.141, 96.0);
            obs.add(90.184, 50.0);
            obs.add(90.227, 50.0);
            obs.add(90.27, 0.0);

            //System.out.println(xArr[i] + "\t" + yArr[i]);
            //        }
            //
            //
            //
            //              //  long startTime = System.currentTimeMillis(); //fetch starting time
                                  double[] parameters = null;
            //
                                                      parameters = GaussianCurveFitter.create().withMaxIterations(1000).fit(obs.toList());
                                                                      System.out.println("111\t" + parameters[0]);
                                                                      System.out.println("111\t" + parameters[1]);
                                                                      System.out.println("111\t" + parameters[2]);
      System.out.println("111\t" + parameters[0]*parameters[2]*106.3472);



            //
            //
            //
            //


if(true) return;

      String configFile = "/data/2/rpark/ip2_data/carolfc/Marta_UCSD_2018_Secretome/labelfree_quant/labelfree_16169/temp/census_config_labelfree_16169.xml";
      String tmpFile = "/data/2/rpark/ip2_data/carolfc/Marta_UCSD_2018_Secretome/labelfree_quant/labelfree_16169/temp/census_labelfree_out_16169.txttmp"; //http://192.168.2.9/ip2/viewLabelfree.html?pid=127&projectName=JonB_Cox
    //  String jsonFile = "/data/2/rpark/ip2_data/carolfc/Marta_UCSD_2018_Secretome/labelfree_quant/labelfree_16169/temp/labelfree_dt_pep_16169.JSON";

      String mainFileName = tmpFile.substring(0, tmpFile.length() - 7);
      String filledFile = mainFileName + "_filled.txt";
      //fill IIT intensities for unidentified peptides




      LabelfreeFilledParserTemp l = new LabelfreeFilledParserTemp(filledFile);
      List<ProteinModel> proteinList = l.readWholeFile(configFile);



      LabelfreeMissingPeptideBuilderSplit.generateLabelfreeOutputFile(proteinList, mainFileName + "_stat.txt", configFile);


      if(true) return;

      float precMass = 100;
      float tolerance = 50;

      float res = tolerance*precMass/1000000f;

      System.out.println(res);

      /*
      double ppm = tolerance/1000;
      ppm += 1;
      ppm = precMass - precMass/ppm;

      System.out.println(ppm/1000);
      */



    if(true) return;



      //BufferedReader br = new BufferedReader(new FileReader("/home/rpark/new.tt"));
      BufferedReader br = new BufferedReader(new FileReader("/home/rpark/new2.txt"));

      String eachLine="";
      List<Double> pValues = new ArrayList<>();
      while( null != (eachLine=br.readLine()) ) {

        /*
        String[] arr = eachLine.split("\t");

        List<Double> list = new ArrayList<>();
        for(String each:arr) {
          if("X".equals(each)) continue;;

          list.add(Double.parseDouble(each));
        }

        if(list.size()<=1) continue;

        double[] darr = new double[list.size()];
        for(int i=0;i<darr.length;i++) {
          darr[i] = list.get(i);
        }

      pValues.add(TTestUtil.oneSampleTTest(darr));
*/
        if("NA".equals(eachLine)) continue;

        pValues.add(Double.parseDouble(eachLine));

      }
      List<Double> qvalues = BHCorrection.runBhCorrection(pValues);
      for(Double d:qvalues)
        System.out.println(d);

      /*
      try {
        long int1 = (long)Double.parseDouble("02510457600");
        System.out.println(int1);
      } catch(Exception ee) {
        ee.printStackTrace();
      }


	/*
	Object[][] arr = new Object[50][200];
	int hw = 5;

        // i, j for matrix array coordinate
        // x, y for chromatograph cooridate


	hw = Integer.parseInt(args[0]);
	int x = Integer.parseInt(args[1]);
	int y = Integer.parseInt(args[2]);


	rotateOrigin(hw, x,y);
	*/

	//System.out.println("== " + MathUtil.getScaled(1.5995, 4) );

    }



	public static void test()
    {
        /* test
        int x=0,y=0;
        int hw=500;

        Object[][] matrix = nwe Object[500][10000];

        for(int i=0;i<1000;i++)
        {
            //increase x
            //increase y
            if( AlignNode.checkWithinBound(hw, x, y) )
            {
                AlignNode node = new AlignNode(x, y, hw);

                matrix[node.getI(), node.getJ()];

            }

        }
        */


    }

    //convert from matrix array coordinate to chromatograph coordinate
    public static void rotateNew(int hw, int i, int j)
    {

	int x;
	int y;

	if(hw>=j)
	{
	    x=i;
	    y=(i+j)-hw;
	}
	else
	{
	    //X=2*hw-1-j;
	    x=i-j+hw;
	    y=i;
	}

	System.out.println(x + " " + y);
    }

    //convert from chromatograph coordinate to matrix array coordinate
    public static void rotateOrigin(int hw, int x, int y)
    {
	int i,j;

	if(y>=x)
	{
	    i=y;
	    j=y-x+hw;
	}
	else
	{
	    i=x;
	    j=y-x+hw;
	}

	System.out.println(i + " " + j);
    }

    //convert from matrix array coordinate to chromatograph coordinate
    public static boolean checkWithinBound(int hw, int x, int y)
    {

	int i,j;

	if(y>=x)
	{
	    i=y;
	    j=y-x+hw;
	}
	else
	{
	    i=x;
	    j=y-x+hw;
	}

	if(j<0||j>=2*hw)
	    return false;
	else
	    return true;
/*
	if(hw>=j)
	{
	    x=i;
	    y=(i+j)-hw;
	}
	else
	{
	    //X=2*hw-1-j;
	    x=i-j+hw;
	    y=i;
	}

	if(y<0 || y>=hw*2)
	    return false;
	else
	    return true;

*/
    }
}

