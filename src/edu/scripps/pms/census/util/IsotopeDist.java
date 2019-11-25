package edu.scripps.pms.census.util;

import edu.scripps.pms.census.conf.*;
import edu.scripps.pms.census.hash.*;

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

public class IsotopeDist
{
    private static final double abund[][] = new double[10][4];
    private static final int DIST_SIZE = 200;
    private static int npeak[] = new int[10];
    private double startMass;
    private double endMass;
    private double avgMass;
    private double[] masslist;
    private double[] relabun = new double[DIST_SIZE];
    private double modShift;

    private int[] element;
    private static double c13=1.003354;
    private static double N15=0.997034;
    //private static double H2=1.005175;
    private static double H2=1.00627675;
    private static double O18=2.004244;  //mass shift from mono
    //private static double H2=2.013;  //find accurate number later

    //private final double isodistThreshold=0.05; // remove masses below 10%
    private static double isodistThreshold; // remove masses below 10%
    private static int isotopePeak = 0;

//    private static double enrichment = 0.98;
    private static double enrichment;
    protected double enrichCorr = -1;
    private long [] chromPeakArr;


    public long[] getChromPeakArr() {
        return chromPeakArr;
    }

    public void setChromPeakArr(long[] chromPeakArr) {
        this.chromPeakArr = chromPeakArr;
    }

    //= 0.98;
    static
    {
	Configuration conf = Configuration.getInstance();
	enrichment = conf.getEnrichment();
	isodistThreshold = conf.getIsodistThreshold();
	isotopePeak = conf.getIsotopePeak();

//System.out.println("=====" + isotopePeak);
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

	// 18O
        abund[9][0] = 100 - enrichment * 100;
        abund[9][1] = (enrichment * 100) + (abund[9][0] * abund[2][1] / 100);
//        abund[9][2] = (enrichment * 100) + (abund[9][0] * abund[2][2] / 100);
        npeak[9] = 2;
	//O(16)     15.994915   99.76    O(17)     16.999131   0.038    O(18)     17.999159    0.20
    }

    public IsotopeDist(int[] element, double modShift, boolean isLight)
    {


        this.element = element;
	this.modShift = modShift;

	if(isLight)
	    calculateLight();
	else
	    calculateHeavy();
    }

    public IsotopeDist() {}

    public static void main(String args[])
    {
        //int[] element = {77, 129, 27, 23, 1, 0, 0, 0, 0, };
        //int[] element = {77, 129, 27, 0, 1, 0, 23, 0, 0, };
        //int[] element = {4, 8, 2, 1, 0, 5, 0, 0, 0, };
        //int[] element = {4, 9, 5, 1, 0, 1, 0, 0, 0, };
        //int[] element = {57, 99, 22, 17, 0, 0, 0, 0, 0, };
        //int[] element = {53,111, 19,11, 1, 0, 6, 0,12, };
	//C, H, O, N, S, Ph, 15N, 2H, 13C
        //int[] element = {71, 108, 28, 19, 0, 1, 0, 0, 0, };
        int[] element = {108, 144, 40, 33, 0, 1, 0, 0, 0, };
        //int[] element = {57, 99, 22, 0, 0, 0, 17, 0, 0, };
        IsotopeDist iso = new IsotopeDist(element, 0, false);
	iso.setEnrichment(0.98);


	double[] arr = iso.getHighMassList();

	for(int i=0;i<arr.length;i++)
	   System.out.println(arr[i] + " " + (arr[i]+1.00728*2)/2);

	   System.out.println("");
	arr = iso.getMasslist();



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
        for (j = 0; j <= 9; ++j) {

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
        mass = 	(element[0] * 12) +
		(element[1] * 1.007825) +
            	(element[2] * 15.9949146) +
            	(element[3] * 14.003074) +
		(element[4] * 31.9720718) +
            	(element[5] * 30.9737620) +
            	(element[6] * 15.00010898) +
		(element[7] * 2.014101787) + //true monoisotopic
             	(element[8] * 13.00335484) +
              	(element[9] * 15.994915);

        masses[0] = mass;

/*
System.out.println("light calc");
System.out.println(element[0] + "\t" + (element[0] * 12));
System.out.println(element[1] + "\t" + (element[1] * 1.007825));
System.out.println(element[2] + "\t" + (element[2] * 15.9949146));
System.out.println(element[3] + "\t" + (element[3] * 14.003074));
System.out.println(element[4] + "\t" + (element[4] * 31.9720718));
System.out.println(element[5] + "\t" + (element[5] * 30.9737620));
System.out.println(element[6] + "\t" + (element[6] * 15.00010898));
System.out.println(element[7] + "\t" + (element[7] * 2.014101787));
System.out.println(element[8] + "\t" + (element[8] * 13.00335484));
System.out.println(element[9] + "\t" + (element[9] * 15.994915));
/*
int sumaa=0;
for(int iii:element) {
System.out.println("===>>" + iii);
sumaa+= iii;
}
System.out.println("l===========" + mass);
System.out.println("sum===========" + sumaa);
System.out.println("masslist 0===========" + masslist[0]);

*/

//try { throw new Exception("00"); } catch (Exception e ) { e.printStackTrace(); }

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


	//double massshift = element[6]*N15 + element[7]*H2 + element[8]*c13 + element[9]*O18;
        double massshift = element[6]*c13 + element[7]*c13 + element[8]*c13 + element[9]*c13;
	int shiftNum = element[6] + element[7] + element[8] + element[9];

	if(shiftNum>DIST_SIZE-1)
	    shiftNum=DIST_SIZE-1;

	relabun[shiftNum-p] = CPATT[shiftNum];
	fracabun[shiftNum-p] = CPATT[shiftNum] / sum; //calculates fractional abundances
	double monoMass = mass + massshift;
	masslist[shiftNum-p] = monoMass;

	i=1;

	for(k=shiftNum-1;k>=p;k--)
	{
            relabun[k-p] = CPATT[k]; //calculates relative abundances
            fracabun[k-p] = CPATT[k] / sum; //calculates fractional abundances
            masslist[k-p] = monoMass - c13*i;
	    i++;
	}


	i=1;

        for (k=shiftNum+1; k <= q; k++) {
            relabun[k-p] = CPATT[k]; //calculates relative abundances
            fracabun[k-p] = CPATT[k] / sum; //calculates fractional abundances
            masslist[k-p] = monoMass + c13*i;
	    i++;
        }

        avgMass = 0;
        for(i=0;i<fracabun.length;i++)
        {
            if(fracabun[i]<=0)
                break;

            avgMass += fracabun[i] * masslist[i];

        }













        //calculates beginning mass for integration
        ii = 0;

        for (ii = 0; ii <= q; ii++) {
//            if (relabun[ii] > ISODIST_THRESHOLD) {
               //beginmass = masslist[ii] - 0.5;

            if (relabun[ii] > isodistThreshold) {
                //beginmass = masslist[ii] - 0.5;
                this.startMass = masslist[ii];
                break;
            }
        }
//        masses[4] = beginmass;
        //this.startMass = beginmass;


        for (ii = 0; ii <= q; ii++) {
            if (relabun[ii] > isodistThreshold) {
                this.endMass = masslist[ii];
                break;
            }
        }


	for(int u=0;u<masslist.length;u++)
	{
//	    System.out.println(u + "mass==>>\t" + masslist[u] + "\t" + relabun[u]);
	    if(masslist[u]==0)
		break;
	}

	for(int u=0;u<relabun.length;u++)
	{
	    if(masslist[u]==0)
		break;
	}
/*
//        masses[5] = endmass;
        this.endMass = endmass;
//        System.out.println(endmass);

        /*
System.out.println("masslist");
        for (int ii0 = 0; ii0 < masslist.length; ii0++) {
            System.out.println(masslist[ii0]);

        }


        for (int ii0 = 0; ii0 < masses.length; ii0++) {
            System.out.println("mass " + ii0 + "\t" + masses[ii0]);

        }
*/
    }

	//cbamberg
    public void calculateLightSimple()
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
        for (j = 0; j <= 9; ++j) {

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

        mass = 	(element[0] * 12) +
		(element[1] * 1.007825) +
            	(element[2] * 15.9949146) +
            	(element[3] * 14.003074) +
		(element[4] * 31.9720718) +
            	(element[5] * 30.9737620) +
            	(element[6] * 15.00010898) +
		(element[7] * 2.014101787) + //true monoisotopic
             	(element[8] * 13.00335484) +
              	(element[9] * 15.994915);

        masses[0] = mass;

/*
System.out.println("light calc");
System.out.println(element[0] + "\t" + (element[0] * 12));
System.out.println(element[1] + "\t" + (element[1] * 1.007825));
System.out.println(element[2] + "\t" + (element[2] * 15.9949146));
System.out.println(element[3] + "\t" + (element[3] * 14.003074));
System.out.println(element[4] + "\t" + (element[4] * 31.9720718));
System.out.println(element[5] + "\t" + (element[5] * 30.9737620));
System.out.println(element[6] + "\t" + (element[6] * 15.00010898));
System.out.println(element[7] + "\t" + (element[7] * 2.014101787));
System.out.println(element[8] + "\t" + (element[8] * 13.00335484));
System.out.println(element[9] + "\t" + (element[9] * 15.994915));

int sumaa=0;
for(int iii:element) {
System.out.println("===>>" + iii);
sumaa+= iii;
}
System.out.println("l===========" + mass);
System.out.println("sum===========" + sumaa);
System.out.println("masslist 0===========" + masslist[0]);



//try { throw new Exception("00"); } catch (Exception e ) { e.printStackTrace(); }

	/*
	for(int iii=0;iii<element.length;iii++)
	    System.out.println("light" + element[iii]);

	System.out.println(element[6] + "\t" + N15 + "\t" + element[7] + "\t" + H2 + " " + element[8] + " " + c13);
	System.out.println(mass +"\t"  );
*/

        /* Calculate average and fractional masses*/
	mass -= 0.00054857; //remove electron
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


	//double massshift = element[6]*N15 + element[7]*H2 + element[8]*c13 + element[9]*O18;
        double massshift = element[6]*c13 + element[7]*c13 + element[8]*c13 + element[9]*c13;
	int shiftNum = element[6] + element[7] + element[8] + element[9];

	if(shiftNum>DIST_SIZE-1)
	    shiftNum=DIST_SIZE-1;

	relabun[shiftNum-p] = CPATT[shiftNum];
	fracabun[shiftNum-p] = CPATT[shiftNum] / sum; //calculates fractional abundances
	double monoMass = mass + massshift;
	masslist[shiftNum-p] = monoMass;

	i=1;

	for(k=shiftNum-1;k>=p;k--)
	{
            relabun[k-p] = CPATT[k]; //calculates relative abundances
            fracabun[k-p] = CPATT[k] / sum; //calculates fractional abundances
            masslist[k-p] = monoMass - c13*i;
	    i++;
	}


	i=1;

        for (k=shiftNum+1; k <= q; k++) {
            relabun[k-p] = CPATT[k]; //calculates relative abundances
            fracabun[k-p] = CPATT[k] / sum; //calculates fractional abundances
            masslist[k-p] = monoMass + c13*i;
	    i++;
        }

        avgMass = 0;
        for(i=0;i<fracabun.length;i++)
        {
            if(fracabun[i]<=0)
                break;

            avgMass += fracabun[i] * masslist[i];

        }


        //calculates beginning mass for integration
        ii = 0;

        for (ii = 0; ii <= q; ii++) {
//            if (relabun[ii] > ISODIST_THRESHOLD) {
               //beginmass = masslist[ii] - 0.5;

            if (relabun[ii] > isodistThreshold) {
                //beginmass = masslist[ii] - 0.5;
                this.startMass = masslist[ii];
                break;
            }
        }
//        masses[4] = beginmass;
        //this.startMass = beginmass;


        for (ii = 0; ii <= q; ii++) {
            if (relabun[ii] > isodistThreshold) {
                this.endMass = masslist[ii];
                break;
            }
        }


	for(int u=0;u<masslist.length;u++)
	{
//	    System.out.println(u + "mass==>>\t" + masslist[u] + "\t" + relabun[u]);
	    if(masslist[u]==0)
		break;
	}

	for(int u=0;u<relabun.length;u++)
	{
	    if(masslist[u]==0)
		break;
	}

	masslist[0] = mass;
    }

//cbamberg
    public void calculateHeavySimple()
    {
//        this.enrichment = conf.getEnrichment();
        int i0, i1, j, k, p, q, l, ii, cnt, i;
        double[] CPATT = new double[200];
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
        for (j = 0; j <= 9; ++j) {

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
			if(k<=0) break;
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
System.out.println("start=====================");
        for(int cpi=0; cpi<D.length; cpi++)
        {
        //	System.out.println(D[cpi]);
        }
    System.out.println("end=====================");
 */
        /* Calculate mono-isotopic mass */
        mass = 	(element[0] * 12) +
		(element[1] * 1.007825) +
            	(element[2] * 15.9949146) +
            	(element[3] * 14.003074) +
		(element[4] * 31.9720718) +
            	(element[5] * 30.9737620) +
            	(element[6] * 15.00010898) +
		(element[7] * 2.014101787) + //true monoisotopic
             	(element[8] * 13.00335484) +
              	(element[9] * 15.994915);

        masses[0] = mass;



/*
System.out.println("heavy calc");
System.out.println(element[0] + "\t" + (element[0] * 12));
System.out.println(element[1] + "\t" + (element[1] * 1.007825));
System.out.println(element[2] + "\t" + (element[2] * 15.9949146));
System.out.println(element[3] + "\t" + (element[3] * 14.003074));
System.out.println(element[4] + "\t" + (element[4] * 31.9720718));
System.out.println(element[5] + "\t" + (element[5] * 30.9737620));
System.out.println(element[6] + "\t" + (element[6] * 15.00010898));
System.out.println(element[7] + "\t" + (element[7] * 2.014101787));
System.out.println(element[8] + "\t" + (element[8] * 13.00335484));
System.out.println(element[9] + "\t" + (element[9] * 15.994915));

int sumaa=0;
for(int iii:element) {
System.out.println("===>>" + iii);
sumaa+= iii;
}
System.out.println("l===========" + mass);
System.out.println("sum===========" + sumaa);
System.out.println("masslist 0===========" + masslist[0]);
*/

//        for (int iii:element)
//            System.out.println("element===" + iii);


//`System.out.println("ref mass======" + mass);

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


	//double massshift = element[6]*N15 + element[7]*H2 + element[8]*c13;
	//double massshift = element[6]*N15 + element[7]*H2 + element[8]*c13 + element[9]*O18;
	double massshift = 0;
//element[6]*N15 + element[7]*H2 + element[8]*c13 + element[9]*O18;
	int shiftNum = element[6] + element[7] + element[8] + element[9];

//	System.out.println("==========" + element[6] + " " + element[7] + " " + element[8] + " " + element[9]);
//	System.out.println("==========" + N15 +  " " + H2  +  " " + c13 + " " + element[9]);
	//int shiftNum = element[6] + element[7] + element[8];

	if(shiftNum>DIST_SIZE-1)
	    shiftNum=DIST_SIZE-1;

//System.out.println("=========" + shiftNum + " " + CPATT.length + " " + relabun.length);
	relabun[shiftNum-p] = CPATT[shiftNum];
	fracabun[shiftNum-p] = CPATT[shiftNum] / sum; //calculates fractional abundances
	double monoMass = mass + massshift;
	masslist[shiftNum-p] = monoMass;

	monoMass -= 0.00054857; //remove electron
	i=1;
        //avgMass = 0;
        //System.out.println("ssss");

	for(k=shiftNum-1;k>=p;k--)
	{
            relabun[k-p] = CPATT[k]; //calculates relative abundances
            fracabun[k-p] = CPATT[k] / sum; //calculates fractional abundances
            masslist[k-p] = monoMass - c13*i;
//System.out.println("------------->" + (k-p) + " " + masslist[k-p] + " " + fracabun[k-p] + " " + relabun[i] + " " + mass  + " " + k);
            //avgMass = fracabun[k-p] * masslist[k-p] + avgMass; //determines avg. mass

	    i++;
	}

        //System.out.println("avg" + avgMass);

	i=1;

        for (k=shiftNum+1; k <= q; k++) {
            relabun[k-p] = CPATT[k]; //calculates relative abundances
            fracabun[k-p] = CPATT[k] / sum; //calculates fractional abundances
            masslist[k-p] = monoMass + c13*i;
//System.out.println("22------------->" + i + " " + masslist[i] + " " + fracabun[i] + " " + relabun[i] + " " + mass  + " " + k);
            //avgMass = fracabun[i] * masslist[i] + avgMass; //determines avg. mass
	    i++;
        }

        //System.out.println("avg" + avgMass);

        avgMass = 0;
        for(i=0;i<fracabun.length;i++)
        {
            if(fracabun[i]<=0)
                break;

            avgMass += fracabun[i] * masslist[i];

        }

        //System.out.println("new avg" + avgMass);

        for (ii = 0; ii <= q; ii++) {
            if (relabun[ii] > isodistThreshold) {
                startMass = masslist[ii]; // - 0.5; //tolerance
                break;
            }
        }

	if(q>99) //distribution elements are not more than 100
	    q = 99;

        for (ii = 0; ii <= q; ii++) {
            if (relabun[ii] > isodistThreshold) {
                endMass = masslist[ii]; // + 0.5; //tolerance
            }
        }


	for(int u=0;u<masslist.length;u++)
	{
	    if(masslist[u]==0)
		break;
	}
	masslist[0] = monoMass;
    }

    private void calculateHeavy()
    {
//        this.enrichment = conf.getEnrichment();
        int i0, i1, j, k, p, q, l, ii, cnt, i;
        double[] CPATT = new double[200];
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
        for (j = 0; j <= 9; ++j) {

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
			if(k<=0) break;
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
System.out.println("start=====================");
        for(int cpi=0; cpi<D.length; cpi++)
        {
        //	System.out.println(D[cpi]);
        }
    System.out.println("end=====================");
 */
        /* Calculate mono-isotopic mass */
        mass = 	(element[0] * 12) +
		(element[1] * 1.007825) +
            	(element[2] * 15.9949146) +
            	(element[3] * 14.003074) +
		(element[4] * 31.9720718) +
            	(element[5] * 30.9737620) +
            	(element[6] * 15.00010898) +
		(element[7] * 2.014101787) + //true monoisotopic
             	(element[8] * 13.00335484) +
              	(element[9] * 15.994915);

        masses[0] = mass;

/*
System.out.println("heavy calc");
System.out.println(element[0] + "\t" + (element[0] * 12));
System.out.println(element[1] + "\t" + (element[1] * 1.007825));
System.out.println(element[2] + "\t" + (element[2] * 15.9949146));
System.out.println(element[3] + "\t" + (element[3] * 14.003074));
System.out.println(element[4] + "\t" + (element[4] * 31.9720718));
System.out.println(element[5] + "\t" + (element[5] * 30.9737620));
System.out.println(element[6] + "\t" + (element[6] * 15.00010898));
System.out.println(element[7] + "\t" + (element[7] * 2.014101787));
System.out.println(element[8] + "\t" + (element[8] * 13.00335484));
System.out.println(element[9] + "\t" + (element[9] * 15.994915));

/*
int sumaa=0;
for(int iii:element) {
System.out.println("===>>" + iii);
sumaa+= iii;
}
System.out.println("l===========" + mass);
System.out.println("sum===========" + sumaa);
System.out.println("masslist 0===========" + masslist[0]);
*/

//        for (int iii:element)
//            System.out.println("element===" + iii);


//`System.out.println("ref mass======" + mass);

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


	double massshift = element[6]*c13 + element[7]*c13 + element[8]*c13;
	//double massshift = element[6]*N15 + element[7]*H2 + element[8]*c13 + element[9]*O18;
	int shiftNum = element[6] + element[7] + element[8];// + element[9];

//	System.out.println("==========" + element[6] + " " + element[7] + " " + element[8] + " " + element[9]);
//	System.out.println("==========" + N15 +  " " + H2  +  " " + c13 + " " + element[9]);
	//int shiftNum = element[6] + element[7] + element[8];

	if(shiftNum>DIST_SIZE-1)
	    shiftNum=DIST_SIZE-1;

//System.out.println("=========" + shiftNum + " " + CPATT.length + " " + relabun.length);
	relabun[shiftNum-p] = CPATT[shiftNum];
	fracabun[shiftNum-p] = CPATT[shiftNum] / sum; //calculates fractional abundances
	double monoMass = mass + massshift;
	masslist[shiftNum-p] = monoMass;

//System.out.println("ref mass======" + monoMass);
	i=1;
        //avgMass = 0;
        //System.out.println("ssss");

	for(k=shiftNum-1;k>=p;k--)
	{
            relabun[k-p] = CPATT[k]; //calculates relative abundances
            fracabun[k-p] = CPATT[k] / sum; //calculates fractional abundances
            masslist[k-p] = monoMass - c13*i;
//System.out.println("------------->" + (k-p) + " " + masslist[k-p] + " " + fracabun[k-p] + " " + relabun[i] + " " + mass  + " " + k);
            //avgMass = fracabun[k-p] * masslist[k-p] + avgMass; //determines avg. mass

	    i++;
	}

        //System.out.println("avg" + avgMass);

	i=1;

        for (k=shiftNum+1; k <= q; k++) {
            relabun[k-p] = CPATT[k]; //calculates relative abundances
            fracabun[k-p] = CPATT[k] / sum; //calculates fractional abundances
            masslist[k-p] = monoMass + c13*i;
//System.out.println("22------------->" + i + " " + masslist[i] + " " + fracabun[i] + " " + relabun[i] + " " + mass  + " " + k);
            //avgMass = fracabun[i] * masslist[i] + avgMass; //determines avg. mass
	    i++;
        }

        //System.out.println("avg" + avgMass);

        avgMass = 0;
        for(i=0;i<fracabun.length;i++)
        {
            if(fracabun[i]<=0)
                break;

            avgMass += fracabun[i] * masslist[i];

        }

        //System.out.println("new avg" + avgMass);

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
            if (relabun[ii] > isodistThreshold) {
                startMass = masslist[ii]; // - 0.5; //tolerance
                break;
            }
        }

//        masses[4] = beginmass;
        //this.startMass = beginmass;

	if(q>99) //distribution elements are not more than 100
	    q = 99;

        for (ii = 0; ii <= q; ii++) {
            if (relabun[ii] > isodistThreshold) {
                endMass = masslist[ii]; // + 0.5; //tolerance
            }
        }


	for(int u=0;u<masslist.length;u++)
	{
	    if(masslist[u]==0)
		break;
	}

        /*
	for(int u=0;u<relabun.length;u++)
	{
	    System.out.println(u + "abun==>>\t" + masslist[u] + "\t" + relabun[u]);
	    if(masslist[u]==0)
		break;
	}
	/*

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

    public double[] getRelabun() {
	return relabun;
    }

    public double[] getRelabun(int size) {


            int start=-1, end=-1;


	for(int i=0;i<relabun.length;i++)
	{
	    if(relabun[i]>=isodistThreshold)
	    {
		start = i;
		break;
	    }
	}

	for(int i=start;i<relabun.length;i++)
	{
	    if(relabun[i]<=isodistThreshold)
	    {
		end = i;
		break;
	    }
	}

	double[] arr = new double[end-start];

	int index=0;


	if(isotopePeak<=0)
		for(int i=start;i<end;i++)
		{
			arr[index] = masslist[i];
			arr[index] += modShift;
			//	    System.out.println(arr[index] + "==\t" + modShift);
			//	    System.out.println(arr[index] + "\t" + modShift);
			index++;
		}
	else {
		int newEnd = start + isotopePeak;
		if(newEnd<=arr.length)
		end = newEnd;
//		end = isotopePeak;

		for(int i=start;i<end;i++)
		{
			arr[index] = masslist[i];
			arr[index] += modShift;
			index++;
		}
	}

	return arr;
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
	    if(relabun[i]>=isodistThreshold)
	    {
		start = i;
		break;
	    }
	}

	for(int i=start;i<relabun.length;i++)
	{
	    if(relabun[i]<=isodistThreshold)
	    {
		end = i;
		break;
	    }
	}

	double[] arr = null;

	int index=0;


	if(isotopePeak<=0) {
    arr = new double[end-start];
    for (int i = start; i < end; i++) {
      arr[index] = masslist[i];
      arr[index] += modShift;
      //	    System.out.println(arr[index] + "==\t" + modShift);
      //	    System.out.println(arr[index] + "\t" + modShift);
      index++;
    }
  } else {
		int newEnd = start + isotopePeak;
		//if(newEnd<=arr.length)
		end = newEnd;
//		end = isotopePeak;

    arr = new double[end-start];
		for(int i=start;i<end;i++)
		{
			arr[index] = masslist[i];
			arr[index] += modShift;
			index++;
		}
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

    public double[] getBestEnrichMassArr()
    {
        return null;
    }

    public double getFixedEnrichRatio() {
        return -1;
    }

    public double getBestEnrichRatio() { return -1; }
    public double getBestEnrichDelCN() { return -1; }

    public double getCorrOnePlus() { return 0; }

    public double getCorrOneMinus() { return 0; }

    public double getEnrichCorr() {
        return enrichCorr;
    }

    public void setEnrichCorr(double enrichCorr) {
        this.enrichCorr = enrichCorr;
    }

 	public static double[] getHeavySilacDist(double[] sampleIsoArr, int[] elementArr) {

		double massshift = elementArr[6]*N15 + elementArr[7]*H2 + elementArr[8]*c13 + elementArr[9]*O18;

		double[] arr = new double[sampleIsoArr.length];
		for(int i=0;i<sampleIsoArr.length;i++)
			arr[i] = sampleIsoArr[i] + massshift;

		return arr;


	}


public void setElement(int[] element) {
	this.element = element;
}
public void setModShift(double  modShift) {
	this.modShift = modShift;
}

}
