package edu.scripps.pms.census;

import edu.scripps.pms.census.model.IsotopeTable;
import edu.scripps.pms.census.conf.Configuration;

import edu.scripps.pms.census.exception.InvalidAAException;


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
public class ElementComposition
{
    //private String peptide;
    public static final String[] AMINO_ACIDS =
    { "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "NTERM", "CTERM", "*", "#", "@", };

    private static final int ISOTOPE_SIZE = 10;
    private int[] sampleArr = new int[ISOTOPE_SIZE];
    private int[] refArr = new int[ISOTOPE_SIZE];
    private int[] heavyArr = new int[ISOTOPE_SIZE];
    private IsotopeTable<String, int[]> isoTable=null;
    private char[] peptide;
    private int start;
    private int length;
    private boolean quantifiable=true;
    private Configuration conf;
    private double modShift=0;

    public ElementComposition(char residue, IsotopeTable<String, int[]> isoTable) throws InvalidAAException
    {
        this(String.valueOf(residue), isoTable);
    }

    public ElementComposition(String peptideStr, IsotopeTable<String, int[]> isoTable) throws InvalidAAException
    {
        this(peptideStr.toCharArray(), 0, peptideStr.length(), isoTable);
    }

    /*
     * length : length of char array to calculate
     */
    public ElementComposition(char[] peptide, int start, int length, IsotopeTable<String, int[]> isoTable) throws InvalidAAException
    {
        this.peptide = peptide;
        this.start = start;
        this.length = length;
        this.isoTable = isoTable;
        this.conf = Configuration.getInstance();
    }

    public int[] getElementSampleArr()
    {
        return sampleArr;
    }

    public int[] getElementRefArr()
    {
        return refArr;
    }

    private void sumSampleElement(int[] arr)
    {
        for(int i=0; i<this.ISOTOPE_SIZE; i++)
        {
            sampleArr[i] += arr[i];
        }
    }

    private void sumHeavyElement(int[] arr)
    {
        if(conf.getSampleNum()<3)
            return;

        for (int i = 0; i < this.ISOTOPE_SIZE; i++)
        {
            heavyArr[i] += arr[i];
        }

    }

    private void sumRefElement(int[] arr)
    {
        if(!conf.isLabeling() || arr == null)
            return;

        for (int i = 0; i < this.ISOTOPE_SIZE; i++)
        {
            refArr[i] += arr[i];
        }

    }


    private void minusSampleElement(int[] arr)
    {
        if(!conf.isLabeling())
            return;

        for (int i = 0; i < this.ISOTOPE_SIZE; i++)
        {
            sampleArr[i] -= arr[i];
        }

    }

    private void minusRefElement(int[] arr)
    {
        if(!conf.isLabeling())
            return;

        for (int i = 0; i < this.ISOTOPE_SIZE; i++)
        {
            refArr[i] -= arr[i];
        }

    }

    public void calculate() throws InvalidAAException
    {
        char c;
        int j;
        int[] arr;

        boolean isMod=false;
        StringBuffer modValue=null;


        try
        {
            for(int i=start; i<length; i++)
            {
                if( peptide[i]=='(' || peptide[i]=='[')
                {
                    isMod = true;
                    modValue = new StringBuffer();
                    continue;
                }

                if( peptide[i]==')' || peptide[i]==']')
                {
                    isMod = false;
                    //System.out.println(isMod + "\t" + peptide[i] + "\t" + modValue.toString());

                    this.modShift += Double.parseDouble(modValue.toString());
                    continue;
                }

                if(isMod)
                {
//                    modValue.append(peptide[i]);
		    modValue.append(peptide[i]);
                    continue;
                }


                if( peptide[i]=='*' || peptide[i]=='@' || peptide[i]=='#' || 'X' == peptide[i])
                    continue;

               // System.out.println("======" + "sample" + peptide[i]);
                sumSampleElement(isoTable.get("sample" + peptide[i]));
                sumRefElement(isoTable.get("ref" + peptide[i]));

		if(conf.getSampleNum()>2)
			sumHeavyElement(isoTable.get("heavy" + peptide[i]));

                j = i+1;

                if(j<peptide.length && (peptide[j]=='*' || peptide[j]=='@' || peptide[j]=='#') )
                {
                    sumSampleElement(isoTable.get("sample" + peptide[j]));
                    sumRefElement(isoTable.get("ref" + peptide[j]));
		    if(conf.getSampleNum()>2)
			    sumHeavyElement(isoTable.get("heavy" + peptide[j]));
                }

            }

            sumSampleElement(isoTable.get("sampleNTERM"));
            sumSampleElement(isoTable.get("sampleCTERM"));

            sumRefElement(isoTable.get("refNTERM"));
            sumRefElement(isoTable.get("refCTERM"));

	    if(conf.getSampleNum()>2) {
		    sumHeavyElement(isoTable.get("heavyNTERM"));
		    sumHeavyElement(isoTable.get("heavyCTERM"));
	    }




        } catch(Exception e)
        {
            quantifiable=false;

         //   e.printStackTrace();
            throw new InvalidAAException();
        }
    }

    public void calculateSample() throws InvalidAAException
    {
        char c;
        int j;
        int[] arr;

        boolean isMod=false;
        StringBuffer modValue=null;


        try
        {
            for(int i=start; i<length; i++)
            {
                if( peptide[i]=='(' || peptide[i]=='[')
                {
                    isMod = true;
                    modValue = new StringBuffer();
                    continue;
                }

                if( peptide[i]==')' || peptide[i]==']')
                {
                    isMod = false;
                    //System.out.println(isMod + "\t" + peptide[i] + "\t" + modValue.toString());

                    this.modShift += Double.parseDouble(modValue.toString());
                    continue;
                }

                if(isMod)
                {
//                    modValue.append(peptide[i]);
		    modValue.append(peptide[i]);
                    continue;
                }


                if( peptide[i]=='*' || peptide[i]=='@' || peptide[i]=='#' || 'X' == peptide[i])
                    continue;

               // System.out.println("======" + "sample" + peptide[i]);
                sumSampleElement(isoTable.get("sample" + peptide[i]));
//                sumRefElement(isoTable.get("ref" + peptide[i]));

		if(conf.getSampleNum()>2)
			sumHeavyElement(isoTable.get("heavy" + peptide[i]));

                j = i+1;

                if(j<peptide.length && (peptide[j]=='*' || peptide[j]=='@' || peptide[j]=='#') )
                {
                    sumSampleElement(isoTable.get("sample" + peptide[j]));
//                    sumRefElement(isoTable.get("ref" + peptide[j]));
		    if(conf.getSampleNum()>2)
			    sumHeavyElement(isoTable.get("heavy" + peptide[j]));
                }

            }

            sumSampleElement(isoTable.get("sampleNTERM"));
            sumSampleElement(isoTable.get("sampleCTERM"));

//            sumRefElement(isoTable.get("refNTERM"));
//            sumRefElement(isoTable.get("refCTERM"));

	    if(conf.getSampleNum()>2) {
		    sumHeavyElement(isoTable.get("refNTERM"));
		    sumHeavyElement(isoTable.get("refCTERM"));
	    }




        } catch(Exception e)
        {
            quantifiable=false;

            e.printStackTrace();
            throw new InvalidAAException();
        }
    }

    public void lightCalculate() throws InvalidAAException
    {
        char c;
        int j;
        int[] arr;

        boolean isMod=false;
        StringBuffer modValue=null;

        try
        {
            for(int i=start; i<length; i++)
            {
                if( peptide[i]=='(' || peptide[i]=='[')
                {
                    isMod = true;
                    modValue = new StringBuffer();
                    continue;
                }

                if( peptide[i]==')' || peptide[i]==']')
                {
                    isMod = false;
                    //System.out.println(isMod + "\t" + peptide[i] + "\t" + modValue.toString());

                    this.modShift += Double.parseDouble(modValue.toString());
                    continue;
                }

                if(isMod)
                {
//                    modValue.append(peptide[i]);
		    modValue.append(peptide[i]);
                    continue;
                }

                if( peptide[i]=='*' || peptide[i]=='@' || peptide[i]=='#' || 'X' == peptide[i])
                    continue;

                sumSampleElement(isoTable.get("sample" + peptide[i]));

                j = i+1;

                if(j<peptide.length && (peptide[j]=='*' || peptide[j]=='@' || peptide[j]=='#') )
                    sumSampleElement(isoTable.get("sample" + peptide[j]));

            }

            sumSampleElement(isoTable.get("sampleNTERM"));
            sumSampleElement(isoTable.get("sampleCTERM"));


        } catch(Exception e)
        {
            quantifiable=false;

            //e.printStackTrace();
            throw new InvalidAAException();
        }
    }

    public void calculateYion()
    {

        char c;
        int j;
        int[] arr;

        boolean isMod=false;
        StringBuffer modValue=null;

        try
        {
            for(int i=start; i<length; i++)
            {
                if( peptide[i]=='(' || peptide[i]=='[')
                {
                    isMod = true;
                    modValue = new StringBuffer();
                    continue;
                }

                if( peptide[i]==')' || peptide[i]==']')
                {
                    isMod = false;
                    //System.out.println(isMod + "\t" + peptide[i] + "\t" + modValue.toString());

                    this.modShift += Double.parseDouble(modValue.toString());
                    continue;
                }

                if(isMod)
                {
//                    modValue.append(peptide[i]);
		    modValue.append(peptide[i]);
                    continue;
                }

                if( peptide[i]=='*' || peptide[i]=='@' || peptide[i]=='#' || 'X' == peptide[i])
                    continue;

                sumSampleElement(isoTable.get("sample" + peptide[i]));
                sumRefElement(isoTable.get("ref" + peptide[i]));

                j = i+1;

                if(j<peptide.length && (peptide[j]=='*' || peptide[j]=='@' || peptide[j]=='#') ) {
                    sumSampleElement(isoTable.get("sample" + peptide[j]));
                    sumRefElement(isoTable.get("ref" + peptide[j]));
                }

            }

            int[] add = {0,2,0,0,0,0,0,0,0,0,};
            sumSampleElement(add);
            sumSampleElement(isoTable.get("sampleCTERM"));
            sumRefElement(add);
            sumRefElement(isoTable.get("refCTERM"));



        } catch(Exception e)
        {
            quantifiable=false;

            //e.printStackTrace();
            //throw new InvalidAAException();
        }
    }

    public void calculateBion()
    {
        // calculate b ion from y ion
        sampleArr[1] -= 1;  //Hydrogen
        sampleArr[2] -= 1;  //Oxygen

  //      for(int i:refArr)
    //        System.out.println("b\t" + i);

        refArr[1] -= 1;
        refArr[2] -= 1;

        //minusSampleElement(isoTable.get("sampleNTERM"));
        //minusRefElement(isoTable.get("refNTERM"));
       // for(int i:refArr)
      //      System.out.println("a\t" + i);

        heavyArr[1] -= 1;
        heavyArr[2] -= 1;

    }

    public void printComposition()
    {
        System.out.println("light");
        for(int i=0; i<this.ISOTOPE_SIZE; i++)
		System.out.print(" " + sampleArr[i]);

		System.out.println("");

		System.out.println("ref");
		for(int i=0; i<this.ISOTOPE_SIZE; i++)
		    System.out.print(" " + refArr[i]);

		System.out.println("heavy");
		for(int i=0; i<this.ISOTOPE_SIZE; i++)
		    System.out.print(" " + heavyArr[i]);

        System.out.println("");
    }

    public boolean isQuantifiable() {
        return quantifiable;
    }

    public double getModShift() {
        return modShift;
    }

    public void setModShift(double modShift) {
        this.modShift = modShift;
    }

    public int[] getHeavyArr() {
        return heavyArr;
    }

    public void setHeavyArr(int[] heavyArr) {
        this.heavyArr = heavyArr;
    }

}
