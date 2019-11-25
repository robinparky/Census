package edu.scripps.pms.census.util;

import edu.scripps.pms.census.conf.*;
import edu.scripps.pms.census.util.*;
import edu.scripps.pms.util.*;
import edu.scripps.pms.census.hash.*;
import java.util.*;

import edu.scripps.pms.census.util.LinearRegression;
import edu.scripps.pms.census.util.LinearRegressionDouble;

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

public class N15EnrichmentCalc extends IsotopeDist
{
    private final double PROTON_MASS = 1.00728; 
    private static final double abund[][] = new double[9][4];
    private static final int DIST_SIZE = 200;
    private static int npeak[] = new int[9];
    private double startMass;
    private double endMass;
    private double avgMass;
    private double[] masslist;
    //private double[] relabun = new double[100];
    private double[] relabun = new double[DIST_SIZE];
    private double modShift;

    private int[] element;
    private static double c13=1.003354;
    private static double N15=0.997034;
    private static double H2=1.00627675;
    //private static double H2=1.005175; 
    //private static double H2=2.013;  //find accurate number later

    private final double ISODIST_THRESHOLD=0.05; // remove masses below 10%
    
//    private static double enrichment = 0.98;
    private double enrichment;
    private double[] bestEnrichMassArr = null;
    private double bestEnrichRatio=-1;
    private double bestEnrichDelCN=-1;
    private double fixedEnrichRatio=-1;
    private double corrOnePlus=0;
    private double corrOneMinus=0;


    public void init()
    {
//	enrichment = 0.94;

        // Carbon
        abund[0][0] = 100.0;
        abund[0][1] = 1.0958793;
        npeak[0] = 2;

        // Hydrogen
        abund[1][0] = 100.0;
        abund[1][1] = 0.014502102;
        npeak[1] = 2;

        // Oxygen
        abund[2][0] = 100.0;
        abund[2][1] = 0.03799194;
        abund[2][2] = 0.20499609;
        npeak[2] = 3;

        // Nitrogen
        abund[3][0] = 100.0;
        abund[3][1] = 0.368351851;
        npeak[3] = 2;

        // Sulfur
        abund[4][0] = 100.0;
        abund[4][1] = 0.784;
        abund[4][2] = 4.442;
        abund[4][3] = 0.014;
        npeak[4] = 4;

        // Phosphorous
        abund[5][0] = 100.0;
        npeak[5] = 4;

        // 15N Nitrogens
        abund[6][0] = 100 - enrichment * 100;
        abund[6][1] = (enrichment * 100) + (abund[3][0] * abund[3][1] / 100);
        npeak[6] = 2;

        // 2H Deuterium
        abund[7][0] = 100 - 98.27;
        //abund[7][1] = 98.27 + (abund[6][0] * abund[1][1] / 100);
        //abund[7][1] = (enrichment*100) + (abund[6][0] * abund[1][1] / 100);
        abund[7][1] = (enrichment*100) + (abund[1][0] * abund[1][1] / 100);
        npeak[7] = 2;

        // 13C Carbons
        abund[8][0] = 100 - enrichment * 100;
        abund[8][1] = (enrichment * 100) + (abund[8][0] * abund[0][1] / 100);
        npeak[8] = 2;
    }

    public N15EnrichmentCalc(
        int[] element, 
	double modShift, 
	double enrichmentMedian, 
	double enrichmentStart, 
	double enrichmentEnd, 
	double enrichmentMaxDeviation, 
	IndexedFile iFile, 
	int scanNum,
	int chargeState) 
	    throws java.io.IOException, edu.scripps.pms.census.exception.CensusIndexOutOfBoundException
    {
        this.element = element;
	this.modShift = modShift;

//System.out.print(scanNum + "\t" + chargeState + "\t");
	for(double d=enrichmentStart; d<=enrichmentEnd; d+=0.01)
        //double d=0.98;
	{
	    d = MathUtil.getScaled(d, 2);	
	    this.enrichment = d;
	    init();
	    calculateHeavy();
            
	    double corr = getEnrichCorrelation(iFile, scanNum, chargeState);

//if(corr>0)
//	System.out.print(MathUtil.getScaled(corr,2) +"\t");
//else
//	System.out.print("0\t");
//	    double[] intList = getHighIntList();
//	    System.exit(0);
            
	    if(enrichCorr<corr)
	    {
                if(enrichCorr<0)
                    this.bestEnrichDelCN = corr;
                else
                    this.bestEnrichDelCN = corr-enrichCorr;
		enrichCorr=corr;
		bestEnrichRatio=d;
                
		bestEnrichMassArr = tempMassArr;
	    }
	
            //System.out.println("enrichment==\t" + corr + "\t" + d + "\t" + bestEnrichRatio + "\t" +enrichCorr + "\t" + this.bestEnrichDelCN);
	}
        
        if(bestEnrichRatio>0) {
            this.enrichment+=0.01;
            init();
            calculateHeavy();
            corrOnePlus = getEnrichCorrelation(iFile, scanNum, chargeState);
        }

        if(bestEnrichRatio<1) {
            this.enrichment-=0.01;
            init();
            calculateHeavy();
            corrOneMinus = getEnrichCorrelation(iFile, scanNum, chargeState);
        }

        
//System.out.println("\n\nenrichment==\t" + bestEnrichRatio + "\t" +enrichCorr + "\t" + this.bestEnrichDelCN);
//System.out.println("");
//if(true) return;

//	System.out.println(enrichmentMaxDeviation);

	if(enrichmentMedian>0) {
	    double upperBound = enrichmentMedian + enrichmentMaxDeviation;
	    double lowerBound = enrichmentMedian - enrichmentMaxDeviation;


	    if(lowerBound>bestEnrichRatio)
		fixedEnrichRatio = lowerBound;
	    else if(upperBound<bestEnrichRatio)
		fixedEnrichRatio = upperBound;
	    else
		fixedEnrichRatio = bestEnrichRatio;

	    this.enrichment = fixedEnrichRatio;
	    init();
	    calculateHeavy();

/*
	    System.out.println("==" + upperBound);
	    System.out.println("==" + lowerBound);
	    System.out.println("==" + enrichmentMedian);
	    System.out.println("==" + enrichmentMaxDeviation);
	    System.out.println("==" + bestEnrichRatio);
	    System.out.println("==" + fixedEnrichRatio);
	    */
	}
	    
//	System.out.println(bestEnrichRatio);

    }

    public double getEnrichCorrelation(IndexedFile iFile, int scanNum, int chargeState) throws java.io.IOException, edu.scripps.pms.census.exception.CensusIndexOutOfBoundException {

	//needed
	//1. spectrum
	double[][] spArr = CalcUtilGeneric.getSpectrumArr(iFile, scanNum);

	double[] refIsoArr = getHighMassList();
	double[] refIntArr = getHighIntList();

	for(int i=0;i<refIsoArr.length;i++) {
	    refIsoArr[i] = (refIsoArr[i]+chargeState*PROTON_MASS)/chargeState;

	}

	setTempMassArr(refIsoArr);
        
	/*
        for(double d : spArr[1])
		System.out.println(d);                
	    
	System.out.println("=========================");                
        for(double d : refIntArr)
		System.out.println(d);                
	System.out.println("=========================");                
	*/
	
	    return calculateCorrelation(spArr[0], spArr[1], refIsoArr, refIntArr);

	}

	private double[] tempMassArr;

	public static double calculateCorrelation(double[] massArr, double[] intArr, double[] isoArr, double[] relIntArr)
	{
	    int totalPeakFound=0;
	    double toleranceSum=0;
	    
	    Configuration conf = Configuration.getInstance();
	    double massTolerance = conf.getMassTolerance();
	    double[] expIntArr = new double[isoArr.length];

	    int plusCount=0;

	    for(int i=0;i<isoArr.length;i++)
	    {
		double tempTolerance = isoArr[i]/1000*massTolerance;

		double sumIntensity=0;
		
		int start = Arrays.binarySearch(massArr, isoArr[i]-tempTolerance);
		if(start<0)
		    start = -start -1;

		int j=0;
		double small=100;
		double massSmall=-1;
		double highinten=0;
		
		boolean isFound = false;
		
		while(true)
		{
		    if(start>=massArr.length)
			break;

		    double temp = isoArr[i]-massArr[start];
		    if(temp<0)
			temp = -temp;

		    if(temp<tempTolerance)
		    {
			sumIntensity+=intArr[start];
			totalPeakFound++;
			toleranceSum += temp;
			
			isFound = true;
			
			double diff = (massArr[start] - isoArr[i]);
			//if(diff<0)
			    //diff = -diff;
			if(Math.abs(diff)<Math.abs(small))
			{
			    small = massArr[start] - isoArr[i];
			    massSmall = massArr[start];
			    highinten = intArr[start];
			}

		    }

		    if(massArr[start]>isoArr[i])
			break;

		    start++;
		}
		
		//if(isFound)
		  //  foundIsoNum++;

		//System.out.println(isoArr[i] + "\t" + sumIntensity + "\t" + relIntArr[i]);
		expIntArr[i] = sumIntensity;
		if(expIntArr[i]>0)
		    plusCount++;
	    }

	    //if experimental peaks are less or equal to 3, make it invalid
	    if(plusCount<=3)
		return -1;

	    LinearRegressionDouble reg = new LinearRegressionDouble(expIntArr, relIntArr, 0, (relIntArr.length-1), 0);

/*
        for(int i=0;i<relIntArr.length;i++) {
            System.out.println("sp==" + isoArr[i] + "\t" + expIntArr[i] + "\t" + relIntArr[i]);
        }
        
        /*
        System.out.println("----------"); 
        for(int i=0;i<massArr.length;i++) {
            System.out.println("e\t" + massArr[i]);
        }

        System.out.println("----------"); 
        for(int i=0;i<expIntArr.length;i++) {
            System.out.println("e\t" + expIntArr[i]);
        }
        System.out.println("----------");                
        for(int i=0;i<relIntArr.length;i++) {
            System.out.println("e\t" + relIntArr[i]);
        }
        System.out.println("----------");          
            System.out.println("e\t" + reg.getCorr());
                    
*/
	return reg.getCorrWithNeg();
    }
    public static void main(String args[])
    {
        //int[] element = {77, 129, 27, 23, 1, 0, 0, 0, 0, };
        //int[] element = {77, 129, 27, 0, 1, 0, 23, 0, 0, };
        //int[] element = {4, 8, 2, 1, 0, 5, 0, 0, 0, };
        //int[] element = {4, 9, 5, 1, 0, 1, 0, 0, 0, };

	/*
        //int[] element = {57, 99, 22, 17, 0, 0, 0, 0, 0, };
        int[] element = {53,111, 19,11, 1, 0, 6, 0,12, };
        //int[] element = {57, 99, 22, 0, 0, 0, 17, 0, 0, };
//C, H, O, N, S, Ph, 15N, 2H, 13C
        i15EnrichmentCalc iso = new N15EnrichmentCalc(element, 0, false, 0.4);
	iso.setEnrichment(0.98);


	double[] arr = iso.getHighMassList();

	for(int i=0;i<arr.length;i++)
	   System.out.println(arr[i] + " " + (arr[i]+1.00728*2)/2); 

	   System.out.println(""); 
	arr = iso.getMasslist();
*/

        //System.out.println(iso.getEachAbund(0, 0));

        //long start = System.currentTimeMillis();

        //System.out.println( iso.getStartMass() + " " + iso.getEndMass());
    }

    private void calculateLight()
    {
        int i0, i1, j, k, p, q, l, ii, cnt, i;
        double[] CPATT = new double[DIST_SIZE];
        double[] D = new double[DIST_SIZE];
        double mass = 0.0;

        double max, prec;
        double sum, maxmass, dif;
        double[] masses = new double[6];
        double[] fracabun = new double[DIST_SIZE];
        masslist = new double[DIST_SIZE];

        p = 0;
        q = 0;
        prec = 0.00001;
        CPATT[0] = 1.0;
        for (j = 0; j <= 8; ++j) {

        int count=1;
            if (element[j] > 0) {

                for (ii = 0; ii < element[j]; ++ii) {

                    for (i1 = 0; i1 <DIST_SIZE; ++i1) D[i1] = 0.0; /* Erase D[] */
                    //D = new double[200]; //this slows down

                    /* Calculate Isotope Distribution */
                    for (k = p; k <= q; ++k) {
                        count++;

                        for (l = 0; l < npeak[j]; ++l) {
                            cnt = k + l;
                            D[cnt] = D[cnt] + CPATT[k] * abund[j][l];
                        }
                    }

                    /* Normalize Intensities */
                    q = q + npeak[j] - 1;

                    max = 0.0;
                    for (k = p; k <= q; ++k) {
                        if (D[k] > max) {
                            max = D[k];
                        }
                    }

                    for (k = p; k <= q; ++k) {
                        D[k] /= max;
                    }

                    /* Eliminate small peaks to the left */
                    for (k = p; k <= q; ++k) {
                        if (D[k] > prec) {
                            p = k;
                            //k = q;
                            break;
                        }
                    }
                    /* Eliminate small peaks to the right */
                    for (k = q; D[k] < prec; --k) {
                        q = k;
                        k = k - 1;
                    }

                    /* Create new isotope pattern */
                    for (i0 = 0; i0 < DIST_SIZE; ++i0) CPATT[i0] = 0.0;
                        /* Clear CPATT[] */

                    for (k = p; k <= q; ++k) {
                        CPATT[k] = D[k];
                    }
                }

            }

        }

        /* Calculate mono-isotopic mass */
        mass = (element[0] * 12) + (element[1] * 1.007825) +
            (element[2] * 15.9949146)
            + (element[3] * 14.003074) + (element[4] * 31.9720718) +
            (element[5] * 30.9737620) +
            (element[6] * 14.003074 + (element[7] * 1.007825) +
             (element[8] * 12.000000000));

        masses[0] = mass;

	/*
	for(int iii=0;iii<element.length;iii++)
	    System.out.println("light" + element[iii]);

	System.out.println(element[6] + "\t" + N15 + "\t" + element[7] + "\t" + H2 + " " + element[8] + " " + c13);
	System.out.println(mass +"\t"  );
*/

        /* Calculate average and fractional masses*/
        sum = 0;
        maxmass = 0;
        dif = 0;
        max = 0;
        for (k = p; k <= q; k++) {
            sum = sum + CPATT[k];
            if (CPATT[k] > max) {
                max = CPATT[k];
                dif = k;
            }
        }

        i = 0;
        avgMass = 0;

        for (k = p; k <= q; k++) {
            relabun[i] = CPATT[k]; //calculates relative abundances
            fracabun[i] = CPATT[k] / sum; //calculates fractional abundances
            masslist[i] = mass + c13*k; //masses of isotopes
            //masslist[i] = mass + k; //masses of isotopes
          //  System.out.println("light------------->" + masslist[i] + " " + fracabun[k-p] + " " + relabun[i] + " " + mass  + " " + k);
            avgMass = fracabun[i] * masslist[i] + avgMass; //determines avg. mass
            i = i + 1;
        }
        //masses[1] = avgMass; //avg. mass
        //masses[2] = mass + dif; //calculate mass at maximum intensity
        //masses[3] = dif - p + 1; //calculate which ion has the max intensity

        //calculates beginning mass for integration
        ii = 0;

        for (ii = 0; ii <= q; ii++) {
            if (relabun[ii] > ISODIST_THRESHOLD) {
                //beginmass = masslist[ii] - 0.5;
                this.startMass = masslist[ii];
                break;
            }
        }
//        masses[4] = beginmass;
        //this.startMass = beginmass;


        for (ii = 0; ii <= q; ii++) {
            
            if (relabun[ii] > ISODIST_THRESHOLD) {                
                this.endMass = masslist[ii];
            }
        }

    }

    private void calculateHeavy()
    {
//        this.enrichment = conf.getEnrichment();
        int i0, i1, j, k, p, q, l, ii, cnt, i;
        double[] CPATT = new double[200];
        double[] D = new double[DIST_SIZE];
        double mass = 0.0;
//abund[6][0]
        double max, prec;
        double sum, maxmass, dif;
        double[] masses = new double[6];
        double[] fracabun = new double[DIST_SIZE];
        masslist = new double[DIST_SIZE];

        p = 0;
        q = 0;
        prec = 0.00001;
        CPATT[0] = 1.0;
        for (j = 0; j <= 8; ++j) {


        int count=1;
            if (element[j] > 0) {


                for (ii = 0; ii < element[j]; ++ii) {

                    for (i1 = 0; i1 <DIST_SIZE; ++i1) D[i1] = 0.0; /* Erase D[] */
                    //D = new double[200]; //this slows down

                    /* Calculate Isotope Distribution */
                    for (k = p; k <= q; ++k) {
                        count++;

                        for (l = 0; l < npeak[j]; ++l) {
                            cnt = k + l;

			    if(k>=(DIST_SIZE-1))
				break;

                            D[cnt] = D[cnt] + CPATT[k] * abund[j][l];

                        }
                    }



                    /* Normalize Intensities */
                    q = q + npeak[j] - 1;

                    max = 0.0;

		    if(k>=(DIST_SIZE-1))
			break;

                    for (k = p; k <= q; ++k) {
                        if (D[k] > max) {
                            max = D[k];
                        }
                    }

                    for (k = p; k <= q; ++k) {
                        D[k] /= max;
                    }

                    /* Eliminate small peaks to the left */
                    for (k = p; k <= q; ++k) {
                        if (D[k] > prec) {
                            p = k;
                            //k = q;
                            break;
                        }
                    }
                    /* Eliminate small peaks to the right */
                    for (k = q; D[k] < prec; --k) {
                        q = k;
                        k = k - 1;
                    }

                    /* Create new isotope pattern */
                    for (i0 = 0; i0 < DIST_SIZE; ++i0) CPATT[i0] = 0.0;
                        /* Clear CPATT[] */

                    for (k = p; k <= q; ++k) {
                        CPATT[k] = D[k];
                    }
                }

            }

        }


/*
        for(int cpi=0; cpi<D.length; cpi++)
        {
        //	System.out.println(D[cpi]);
        }
 */
        /* Calculate mono-isotopic mass */
        mass = (element[0] * 12) + (element[1] * 1.007825) +
            (element[2] * 15.9949146)
            + (element[3] * 14.003074) + (element[4] * 31.9720718) +
            (element[5] * 30.9737620) +
            (element[6] * 14.003074 + (element[7] * 1.007825) + //true monoisotopic
            //(element[6] * 15.00010897 + (element[7] * 1.007825) + //100% enriched base peak
             (element[8] * 12.000000000));

        masses[0] = mass;

        /* Calculate average and fractional masses*/
        sum = 0;
        maxmass = 0;
        dif = 0;
        max = 0;
        for (k = p; k <= q; k++) {
            sum = sum + CPATT[k];
            if (CPATT[k] > max) {
                max = CPATT[k];
                dif = k;
            }
        }

        i = 0;
        

/*
        for (k = p; k <= q; k++) {
            relabun[i] = CPATT[k]; //calculates relative abundances
            fracabun[i] = CPATT[k] / sum; //calculates fractional abundances
//            masslist[i] = mass + c13*k; //masses of isotopes
//            masslist[i] = mass + N15*k; //masses of isotopes
            masslist[i] = mass + k; //masses of isotopes
System.out.println(masslist[i] + " " + relabun[i] + " " + mass  + " " + k);
            avgMass = fracabun[i] * masslist[i] + avgMass; //determines avg. mass
            i = i + 1;
        }
*/
	double massshift = element[6]*N15 + element[7]*H2 + element[8]*c13;

	int shiftNum = element[6] + element[7] + element[8];

	if(shiftNum>DIST_SIZE-1)
	    shiftNum=DIST_SIZE-1;

	if((shiftNum-p)<0 || relabun.length<=(shiftNum-p)) return;
        //System.out.println(shiftNum + " " + CPATT.length + " " + relabun.length + " " + (shiftNum-p));        
	relabun[shiftNum-p] = CPATT[shiftNum];
	fracabun[shiftNum-p] = CPATT[shiftNum] / sum; //calculates fractional abundances
	double monoMass = mass + massshift;
	masslist[shiftNum-p] = monoMass;

	i=1;
        //avgMass = 0;
        //System.out.println("ssss");

	for(k=shiftNum-1;k>=p;k--)
	{
            relabun[k-p] = CPATT[k]; //calculates relative abundances
            fracabun[k-p] = CPATT[k] / sum; //calculates fractional abundances
            masslist[k-p] = monoMass - c13*i;
            //avgMass = fracabun[k-p] * masslist[k-p] + avgMass; //determines avg. mass

	    i++;
	}
        

	i=1;

        for (k=shiftNum+1; k <= q; k++) {
            relabun[k-p] = CPATT[k]; //calculates relative abundances
            fracabun[k-p] = CPATT[k] / sum; //calculates fractional abundances
            masslist[k-p] = monoMass + c13*i;
            //avgMass = fracabun[i] * masslist[i] + avgMass; //determines avg. mass
	    i++;
        }
        
        
        avgMass = 0;
        for(i=0;i<fracabun.length;i++)
        {
            if(fracabun[i]<=0)
                break;
            
            avgMass += fracabun[i] * masslist[i];
            
        }
        
            
/*
for(i=0;i<40;i++)
    System.out.println(i + "\t" + masslist[i] + " " + relabun[i]);
    System.out.println("----------------");
*/
        //calculates beginning mass for integration
        //ii = 0;

	//double[] tempArr = new double[100];
	//int tempIndex=0;
	/*for(ii=element[6];ii>0;ii--)
	{
	    //tempArr[tempIndex++] = masslist[0]-ii*N15;
	    tempArr[tempIndex++] = masslist[0]-(ii)*N15;
	}

	for(ii=0;ii<masslist.length;ii++)
	{
	    if(masslist[ii]<=0)
		break;
	}
	     
	for(ii=0;ii<masslist.length;ii++)
	{
	    if(masslist[ii]<=0 || relabun[element[6] + ii]<=0)
		break;

	    tempArr[element[6]+ii] = masslist[ii];
	}

	masslist = tempArr;
        */

        for (ii = 0; ii <= q; ii++) {
            if (relabun[ii] > ISODIST_THRESHOLD) {
                startMass = masslist[ii]; // - 0.5; //tolerance
                break;
            }
        }
        
//        masses[4] = beginmass;
        //this.startMass = beginmass;

	if(q>99) //distribution elements are not more than 100
	    q = 99;

        for (ii = 0; ii <= q; ii++) {
            if (relabun[ii] > ISODIST_THRESHOLD) {
                endMass = masslist[ii]; // + 0.5; //tolerance
            }
        }

/*

	for(int u=0;u<masslist.length;u++)
	{
	    System.out.println(u + "mass==>>\t" + masslist[u] + "\t" + relabun[u]);
	    if(masslist[u]==0)
		break;
	}
	for(int u=0;u<relabun.length;u++)
	{
	    System.out.println(u + "abun==>>\t" + relabun[u]);
	    if(masslist[u]==0)
		break;
	}
	
//        masses[5] = endmass;
        this.endMass = endmass;
//        System.out.println(endmass);

        /*
        for (int ii0 = 0; ii0 < masslist.length; ii0++) {
            System.out.println(masslist[ii0]);

        }


        for (int ii0 = 0; ii0 < masses.length; ii0++) {
            System.out.println("mass " + ii0 + "\t" + masses[ii0]);

        }
*/
    }

    public double getEachAbund(int i, int j)
    {
        return abund[i][j];
    }

    public double getStartMass()
    {
        return startMass+modShift;
    }

    public double getEndMass()
    {
        return endMass+modShift;
    }

    public double[] getMasslist()
    {
	for(int i=0;i<masslist.length;i++)
	    masslist[i] += modShift;

        return masslist;
    }

    public double[] getHighMassList()
    {
	int start=-1, end=-1;

	for(int i=0;i<relabun.length;i++)
	{
	    if(relabun[i]>=ISODIST_THRESHOLD)
	    {
		start = i;
		break;
	    }
	}

	for(int i=start;i<relabun.length;i++)
	{
	    if(relabun[i]<=ISODIST_THRESHOLD)
	    {
		end = i;		
		break;
	    }
	}

	double[] arr = new double[end-start];

	int index=0;
	for(int i=start;i<end;i++)
	{
	    arr[index] = masslist[i];
	    arr[index] += modShift;
	    index++;
	}

	return arr;
    }
   
    //return intensity list higher than threshold
    public double[] getHighIntList()
    {
	int start=-1, end=-1;

	for(int i=0;i<relabun.length;i++)
	{
	    if(relabun[i]>=ISODIST_THRESHOLD)
	    {
		start = i;
		break;
	    }
	}

	for(int i=start;i<relabun.length;i++)
	{
	    if(relabun[i]<=ISODIST_THRESHOLD)
	    {
		end = i;		
		break;
	    }
	}

	double[] arr = new double[end-start];

	int index=0;
	for(int i=start;i<end;i++)
	{
	    arr[index] = relabun[i];
	    index++;
	}

	return arr;
    }
    public double getAvgMass() {
        return avgMass + modShift;
    }

    public void setEnrichment(double enrichment)
    {
	this.enrichment = enrichment;
    }

    public double getBestEnrichRatio() {
        return bestEnrichRatio;
    }

    public void setBestEnrichRatio(double bestEnrichRatio) {
        this.bestEnrichRatio = bestEnrichRatio;
    }

    public double getFixedEnrichRatio() {
        return fixedEnrichRatio;
    }

    public void setFixedEnrichRatio(double fixedEnrichRatio) {
        this.fixedEnrichRatio = fixedEnrichRatio;
    }

/*
    public double getEnrichCorr() {
        return enrichCorr;
    }

    public void setEnrichCorr(double enrichCorr) {
        this.enrichCorr = enrichCorr;
    }
*/
    public double[] getTempMassArr() {
        return tempMassArr;
    }

    public void setTempMassArr(double[] tempMassArr) {
        this.tempMassArr = tempMassArr;
    }

    public double[] getBestEnrichMassArr() {
        return bestEnrichMassArr;
    }

    public void setBestEnrichMassArr(double[] bestEnrichMassArr) {
        this.bestEnrichMassArr = bestEnrichMassArr;
    }

    public double getBestEnrichDelCN() {
        return bestEnrichDelCN;
    }

    public void setBestEnrichDelCN(double bestEnrichDelCN) {
        this.bestEnrichDelCN = bestEnrichDelCN;
    }

    public double getCorrOnePlus() {
        return corrOnePlus;
    }

    public void setCorrOnePlus(double corrOnePlus) {
        this.corrOnePlus = corrOnePlus;
    }

    public double getCorrOneMinus() {
        return corrOneMinus;
    }

    public void setCorrOneMinus(double corrOneMinus) {
        this.corrOneMinus = corrOneMinus;
    }
    
    
}
