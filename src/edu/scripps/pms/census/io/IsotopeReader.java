package edu.scripps.pms.census.io;

import java.io.IOException;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.File;
import java.util.*;

import edu.scripps.pms.census.model.IsotopeTable;
import edu.scripps.pms.census.*;
import edu.scripps.pms.census.exception.InvalidAAException;
import edu.scripps.pms.census.util.IsotopeDist;

import org.jdom.Element;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: Yates Lab</p>
 *
 * @author Robin Park
 * @version 1.0
 */
public class IsotopeReader
{
    private BufferedReader br;
    //private static final int ISOTOPE_SIZE = 9;
    private static final int ISOTOPE_SIZE = 10;
    private IsotopeTable<String, int[]> isotope;
    private Hashtable<String, Double> aaHt = new Hashtable<String, Double>();

    public double getAAMonoMass(char residue)
    {
        return getAAMonoMass(String.valueOf(residue));
    }
    
    public double getAAMonoMass(String residue)
    {
/*
        if(null == aaHt.get(residue))
        {
            System.out.println(aaHt);
            
            System.out.println("===>>" + " " + residue);
        }
*/
        
        return aaHt.get(residue);
    }
    
    private Element rootEle;
    
    /*
         array[0] : Carbon
         array[1] : Hydrogen
         array[2] : Oxygen
         array[3] : Nitrogen
         array[4] : Sulfur
         array[5] : Phosphorous
         array[6] : 15N
         array[7] : 2H
         array[8] : 13C
     */

    public IsotopeReader(String fileName) throws IOException
    {
        br = new BufferedReader(new FileReader(fileName));
        init();
        
    }
    
    public IsotopeReader(Element rootEle) throws IOException, InvalidAAException
    {
        this.rootEle = rootEle;
        
        initXml();                
    }
       

    public IsotopeReader(File file) throws IOException
    {
        br = new BufferedReader(new FileReader(file));
        init();
    }
    
    
    //private final String[] PROTEIN_COLUMNS = {"C", "H", "O", "N", "S", "P", "15N", "2H", "13C", };
    private final static String[] AA = {"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "NTERM", "CTERM",};

    private final static int[][] DEFAULT_ELE_COMP = {                
        {3, 5, 1, 1, 0, 0, 0, 0, 0, },
        {5, 8, 2, 2, 1, 0, 0, 0, 0, },
        {4, 5, 3, 1, 0, 0, 0, 0, 0, },
        {5, 7, 3, 1, 0, 0, 0, 0, 0, },
        {9, 9, 1, 1, 0, 0, 0, 0, 0, },
        {2, 3, 1, 1, 0, 0, 0, 0, 0, },
        {6, 7, 1, 3, 0, 0, 0, 0, 0, },
        {6, 11, 1, 1, 0, 0, 0, 0, 0, },
        {6, 12, 1, 2, 0, 0, 0, 0, 0, },
        {6, 11, 1, 1, 0, 0, 0, 0, 0, },
        {5, 9, 1, 1, 1, 0, 0, 0, 0, },
        {4, 6, 2, 2, 0, 0, 0, 0, 0, },
        {5, 7, 1, 1, 0, 0, 0, 0, 0, },
        {5, 8, 2, 2, 0, 0, 0, 0, 0, },
        {6, 12, 1, 4, 0, 0, 0, 0, 0, },
        {3, 5, 2, 1, 0, 0, 0, 0, 0, },
        {4, 7, 2, 1, 0, 0, 0, 0, 0, },
        {5, 9, 1, 1, 0, 0, 0, 0, 0, },
        {11, 10, 1, 2, 0, 0, 0, 0, 0, },
        {9, 9, 2, 1, 0, 0, 0, 0, 0, },
        {0, 1, 0, 0, 0, 0, 0, 0, 0, },
        {0, 1, 1, 0, 0, 0, 0, 0, 0, },
    };
    
    //this is standard isotope table
    public static IsotopeTable<String, int[]> getStandardIsotopeTable()
    {
        IsotopeTable<String, int[]> isotope = new IsotopeTable<String, int[]>();
        
        //System.out.println(AA.length + " " + DEFAULT_CELL_DATA.length);
    //        System.out.println("----*********---->>");
        for(int i=0;i<AA.length;i++)
        {
            isotope.put(AA[i], DEFAULT_ELE_COMP[i]);
      //      System.out.println("-------->>" + AA[i] + " " + DEFAULT_ELE_COMP[i]);
            //DEFAULT_ELE_COMP
        }
        
        return isotope;
    }
    
    private void initXml() throws IOException, InvalidAAException
    {        
        isotope = new IsotopeTable<String, int[]>();
        String eachLine;
        String sampleName;
        String refName;
        String[] arr;
        //Use integer array due to intensive calculation later
        int[] sampleArr = null;
        int[] refArr = null;
        int[] heavyArr = null;
        
        List<Element> l = this.rootEle.getChild("element_comp").getChildren("each_sample");

        boolean isLabeling = false;
	Element labelTypeEle = this.rootEle.getChild("label_type");
	if(null != labelTypeEle)
		isLabeling = labelTypeEle.getAttributeValue("labeling").equals("false")?false:true;
        
        Element sampleEle = l.get(0);        
        
        Element refEle = null;
        Element heavyEle = null;
        List<Element> refResidueList = null;
        List<Element> heavyResidueList = null;

	Element expTypeEle = this.rootEle.getChild("experiment_type");
	String expType = "";
	if(null != expTypeEle)
	    expType = expTypeEle.getText();

	if(expType.equals("13") || 
		expType.equals("14") || 
		expType.equals("201") )
	    isLabeling = false;

        if(isLabeling && l.size()>1)
        {
            refEle = l.get(1);
            refResidueList = refEle.getChildren("residue");
        }
        if(l.size()>2)
        {
            heavyEle = l.get(2);
            heavyResidueList = heavyEle.getChildren("residue");
        }
        
        List<Element> samResidueList = sampleEle.getChildren("residue");
        
        for(int i=0;i<samResidueList.size();i++)
        {
            Element eachSamEle = samResidueList.get(i);            
            sampleArr = new int[ISOTOPE_SIZE];
            
            Element eachRefEle = null;
            Element eachHeavyEle = null;
            
            if(isLabeling && null != refResidueList)
            {
                eachRefEle = refResidueList.get(i);
                refArr = new int[ISOTOPE_SIZE];
            }
	    if(l.size()>2)
	    {
                eachHeavyEle = heavyResidueList.get(i);
                heavyArr = new int[ISOTOPE_SIZE];
	    }

            for(int j=0;j<chemArr.length;j++)             {
		String tmpEle = eachSamEle.getChildText("ele_" + chemArr[j]);
                
                try {
		if(null != tmpEle)
		    sampleArr[j] = Integer.parseInt( tmpEle); //eachSamEle.getChildText("ele_" + chemArr[j]) );
                
                }catch (Exception e) {
                    System.out.println(tmpEle);
                }

	    }
            
            String eachName = eachSamEle.getAttributeValue("name");
            
            isotope.put("sample" + eachName, sampleArr);
            
            //ElementComposition element = new ElementComposition(eachName, null); //isoReader.getIsotope());

            IsotopeDist iDist = new IsotopeDist(sampleArr, 0, false);
        
            getAaHt().put(eachName, iDist.getStartMass());                   
        
            if(isLabeling && null != eachRefEle)
            {
                for(int j=0;j<chemArr.length;j++) {
		    String tmpEle = eachRefEle.getChildText("ele_" + chemArr[j]); 
		    if(null != tmpEle)
			refArr[j] = Integer.parseInt( tmpEle); //eachRefEle.getChildText("ele_" + chemArr[j]) );
		}
		
                isotope.put("ref" + eachRefEle.getAttributeValue("name"), refArr);

//		System.out.println("===" + eachRefEle.getAttributeValue("name"));
            }            
            if(l.size()>2)
            {
                for(int j=0;j<chemArr.length;j++) {
		    String tmpEle = eachHeavyEle.getChildText("ele_" + chemArr[j]); 
		    if(null != tmpEle)
			heavyArr[j] = Integer.parseInt( tmpEle); 
		}
		
                isotope.put("heavy" + eachHeavyEle.getAttributeValue("name"), heavyArr);

//		System.out.println("===" + eachRefEle.getAttributeValue("name"));
            }            
        }   



/*
        int[] tempArr = isotope.get("sampleK");
for(int iii:tempArr)
    System.out.println("==========" + iii);
*/
    }
    
    private void init() throws IOException
    {        
        isotope = new IsotopeTable<String, int[]>();
        String eachLine;
        String sampleName;
        String refName;
        String[] arr;
        //Use integer array due to intensive calculation later
        int[] sampleArr;
        int[] refArr;

        while ( (eachLine = br.readLine()) != null)// && eachLine.startsWith("<"))
        {
            if(eachLine.startsWith("<"))
            {
                sampleArr = new int[ISOTOPE_SIZE];
                refArr = new int[ISOTOPE_SIZE];

                String str = eachLine.substring(1, eachLine.indexOf(">"));

                //Read carbon
                eachLine = br.readLine();
                arr = eachLine.split("\t");
                sampleArr[0] = Integer.parseInt(arr[1]);
                refArr[0] = Integer.parseInt(arr[2]);

                //Read Hydrogen
                eachLine = br.readLine();
                arr = eachLine.split("\t");
                sampleArr[1] = Integer.parseInt(arr[1]);
                refArr[1] = Integer.parseInt(arr[2]);

                //Read Oxygen
                eachLine = br.readLine();
                arr = eachLine.split("\t");
                sampleArr[2] = Integer.parseInt(arr[1]);
                refArr[2] = Integer.parseInt(arr[2]);

                //Read Nitrogen
                eachLine = br.readLine();
                arr = eachLine.split("\t");
                sampleArr[3] = Integer.parseInt(arr[1]);
                refArr[3] = Integer.parseInt(arr[2]);

                //Read Sulfur
                eachLine = br.readLine();
                arr = eachLine.split("\t");
                sampleArr[4] = Integer.parseInt(arr[1]);
                refArr[4] = Integer.parseInt(arr[2]);

                //Read Phosphorous
                eachLine = br.readLine();
                arr = eachLine.split("\t");
                sampleArr[5] = Integer.parseInt(arr[1]);
                refArr[5] = Integer.parseInt(arr[2]);

                //Read 15N
                eachLine = br.readLine();
                arr = eachLine.split("\t");
                sampleArr[6] = Integer.parseInt(arr[1]);
                refArr[6] = Integer.parseInt(arr[2]);

                //Read 2H
                eachLine = br.readLine();
                arr = eachLine.split("\t");
                sampleArr[7] = Integer.parseInt(arr[1]);
                refArr[7] = Integer.parseInt(arr[2]);

                //Read 13C
                eachLine = br.readLine();
                arr = eachLine.split("\t");
                sampleArr[8] = Integer.parseInt(arr[1]);
                refArr[8] = Integer.parseInt(arr[2]);

                isotope.put("sample" + str, sampleArr);
                isotope.put("ref" + str, refArr);

            }

            //System.out.println("==>" + eachLine + "<--");
        }
    }

    String[] chemArr = {"C", "H", "O", "N", "S", "P", "15N", "2H", "13C", "18O",};

    public static void test(IsotopeTable<String, int[]> ht) throws Exception
    {

String[][] arr = {
{"A", },
{"C", },
{"D", },
{"E", },
{"F", },
{"G", },
{"H", },
{"I", },
{"K", },
{"L", },
{"M", },
{"N", },
{"P", },
{"Q", },
{"R", },
{"S", },
{"T", },
{"V", },
{"W", },
{"Y", },
{"NTERM", },
{"CTERM", },
{"*", },
{"#", },
{"@", },
};

for(int j=0;j<arr.length;j++)
{

        ElementComposition element = new ElementComposition(arr[j][0], ht);
        int[] tempA = element.getElementSampleArr();


	System.out.print("{\"" + arr[j][0] + "\", ");
        for(int i=0;i<tempA.length;i++)
            System.out.print("\"" + ht.get("sample" + arr[j][0])[i] + "\", ");
	System.out.println("},");
}

	System.out.println("");
for(int j=0;j<arr.length;j++)
{

        ElementComposition element = new ElementComposition(arr[j][0], ht);
        int[] tempA = element.getElementSampleArr();


	System.out.print("{\"" + arr[j][0] + "\", ");
        for(int i=0;i<tempA.length;i++)
            System.out.print("\"" + ht.get("ref" + arr[j][0])[i] + "\", ");
	System.out.println("},");
}
    }

    public static void main(String args[]) throws Exception
    {
        IsotopeReader reader = new IsotopeReader(args[0]);
        //IsotopeReader reader = new IsotopeReader("/data/1/rpark/relex_run/RelEx_new_data/data-independent/N15isotope.ini");

        IsotopeTable<String, int[]> ht = reader.getIsotope();

System.out.println("start");

test(ht);

/*
        ElementComposition element = new ElementComposition("C", ht);
        //ElementComposition element = new ElementComposition("AAADLMAYCEAHAKEDPLLTPVPASENPFR", ht);

        int[] tempA = element.getElementSampleArr();

        for(int i=0;i<tempA.length;i++)
            System.out.print("\"" + tempA[i] + "\", ");

        /*

        for(java.util.Enumeration enu = ht.keys(); enu.hasMoreElements(); )
        {
            Object obj = enu.nextElement();
            System.out.println (obj + "====>>" +  ht.get(obj)[0] );
        }

        System.out.println(ht.get("sampleA")[0]);
        System.out.println(ht.get("sampleR")[0]);
        System.out.println(ht.get("sampleN")[0]);
        System.out.println(ht.get("sampleD")[0]);
        System.out.println(ht.get("sampleNTERM")[0]);
        System.out.println(ht.get("refNTERM")[0]);
        System.out.println(ht.get("ref*")[0]);
        System.out.println(ht.get("ref@")[0]);
        System.out.println(ht.get("ref#")[0]);
*/
    }

    public void setIsotope(IsotopeTable isotope)
    {
        this.isotope = isotope;
    }

    public IsotopeTable<String, int[]> getIsotope()
    {
        return isotope;
    }

    public Hashtable<String, Double> getAaHt() {
        return aaHt;
    }

    public void setAaHt(Hashtable<String, Double> aaHt) {
        this.aaHt = aaHt;
    }

}
