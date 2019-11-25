package edu.scripps.pms.averagine;

/**
 * Created by rpark on 9/15/16.
 */
public class AssignMass {
    /*
    public static final int SIZE= 256;
    public static final float DIFFMASSC12C13 = 1.003354826f;

    public static final float MADD_DIFF_C12C13 = 1.003354826f;
    public static final int MADD_DIFF_C12C13_PPM = 1003;
    private static final float FLOAT_ZERO = 0.000001f;


    private static final float [] aaMassAvg = new float[SIZE];
    private static final float [] aaMassMono = new float[SIZE];
    public static final float NH3 = 17.02647f;
    public static final float H2O = 18.01051f;
    public static final float CO = 27.99491f;
    public static final float NH3_CS2 = 8.513235f;
    public static final float H2O_CS2 = 9.005255f;
    public static final float CO_CS2 = 13.997455f;

    public static final int NH3_INT = 17026;
    public static final int H2O_INT = 18011;
    public static final int CO_INT = 27995;
    public static final int NH3_CS2_INT = 8513;
    public static final int H2O_CS2_INT = 9005;
    public static final int CO_CS2_INT = 13998;


    private static float[] aaMasses;
    private static AssignMass assignMass;
    private static float cTerm;
    private static float nTerm;
    private static float yionfragment;
    private static float bionfragment;

    public static final int[][] intensePeaks = new int[20][];
    public static final int[] mostIntensePeaks = new int[20];

    static {
        aaMassAvg['G'] =  57.05192f;   aaMassMono['G'] =  57.0214636f;
        aaMassAvg['A'] =  71.07880f;   aaMassMono['A'] =  71.0371136f;
        aaMassAvg['S'] =  87.07820f;   aaMassMono['S'] =  87.0320282f;
        aaMassAvg['P'] =  97.11668f;   aaMassMono['P'] =  97.0527636f;
        aaMassAvg['V'] =  99.13256f;   aaMassMono['V'] =  99.0684136f;
        aaMassAvg['T'] = 101.10508f;   aaMassMono['T'] = 101.0476782f;
        aaMassAvg['C'] = 103.13880f;   aaMassMono['C'] = 103.0091854f;
        aaMassAvg['L'] = 113.15944f;   aaMassMono['L'] = 113.0840636f;
        aaMassAvg['I'] = 113.15944f;   aaMassMono['I'] = 113.0840636f;
        aaMassAvg['X'] = 113.15944f;   aaMassMono['X'] = 113.0840636f;
        aaMassAvg['N'] = 114.10384f;   aaMassMono['N'] = 114.0429272f;
        aaMassAvg['O'] = 114.14720f;   aaMassMono['O'] = 114.0793126f;
        aaMassAvg['B'] = 114.59622f;   aaMassMono['B'] = 114.5349350f;
        aaMassAvg['D'] = 115.08860f;   aaMassMono['D'] = 115.0269428f;
        aaMassAvg['Q'] = 128.13072f;   aaMassMono['Q'] = 128.0585772f;
        aaMassAvg['K'] = 128.17408f;   aaMassMono['K'] = 128.0949626f;
        aaMassAvg['Z'] = 128.62310f;   aaMassMono['Z'] = 128.5505850f;
        aaMassAvg['E'] = 129.11548f;   aaMassMono['E'] = 129.0425928f;
        aaMassAvg['M'] = 131.19256f;   aaMassMono['M'] = 131.0404854f;
        aaMassAvg['H'] = 137.14108f;   aaMassMono['H'] = 137.0589116f;
        aaMassAvg['F'] = 147.17656f;   aaMassMono['F'] = 147.0684136f;
        aaMassAvg['R'] = 156.18748f;   aaMassMono['R'] = 156.1011106f;
        aaMassAvg['Y'] = 163.17596f;   aaMassMono['Y'] = 163.0633282f;
        aaMassAvg['W'] = 186.21320f;   aaMassMono['W'] = 186.0793126f;
    }
    */

//    private static final double [][] isotopicDistribution = new double[];
    private static final double[] isotopicDistribution0 = new double[3];
    private static final double[] isotopicDistribution1 = new double[4];
    private static final double[] isotopicDistribution2 = new double[5];
    private static final double[] isotopicDistribution3 = new double[5];
    private static final double[] isotopicDistribution4 = new double[6];
    private static final double[] isotopicDistribution5 = new double[6];
    private static final double[] isotopicDistribution6 = new double[7];
    private static final double[] isotopicDistribution7 = new double[8];
    private static final double[] isotopicDistribution8 = new double[8];
    private static final double[] isotopicDistribution9 = new double[9];
    private static final double[] isotopicDistribution10 = new double[9];
    private static final double[] isotopicDistribution11 = new double[10];
    private static final double[] isotopicDistribution12 = new double[10];
    private static final double[] isotopicDistribution13 = new double[11];
    private static final double[] isotopicDistribution14 = new double[12];
    private static final double[] isotopicDistribution15 = new double[12];
    private static final double[] isotopicDistribution16 = new double[12];
    private static final double[] isotopicDistribution17 = new double[13];
    private static final double[] isotopicDistribution18 = new double[13];
    private static final double[] isotopicDistribution19 = new double[14];

    static {
        // from averagine
        
        isotopicDistribution0[0] = 1.000000;
        isotopicDistribution0[1] = 0.276950;
        isotopicDistribution0[2] = 0.051217;
        //isotopicDistribution0[3] = 0.007109;

        isotopicDistribution1[0] = 1.000000;
        isotopicDistribution1[1] = 0.555776;
        isotopicDistribution1[2] = 0.178136;
        isotopicDistribution1[3] = 0.041790;
        //isotopicDistribution1[4] = 0.007894;
        //isotopicDistribution1[5] = 0.001263;

        isotopicDistribution2[0] = 1.000000;
        isotopicDistribution2[1] = 0.845212;
        isotopicDistribution2[2] = 0.437964;
        isotopicDistribution2[3] = 0.168570;
        isotopicDistribution2[4] = 0.052286;
        //isotopicDistribution2[5] = 0.013663;
        //isotopicDistribution2[6] = 0.003097;

        isotopicDistribution3[0] = 0.891137;
        isotopicDistribution3[1] = 1.000000;
        isotopicDistribution3[2] = 0.644526;
        isotopicDistribution3[3] = 0.303221;
        isotopicDistribution3[4] = 0.114261;
       // isotopicDistribution3[5] = 0.036227;
        //isotopicDistribution3[6] = 0.009967;
        //isotopicDistribution3[7] = 0.002431;

        isotopicDistribution4[0] = 0.712904;
        isotopicDistribution4[1] = 1.000000;
        isotopicDistribution4[2] = 0.775825;
        isotopicDistribution4[3] = 0.432162;
        isotopicDistribution4[4] = 0.191232;
        isotopicDistribution4[5] = 0.070896;
        //isotopicDistribution4[6] = 0.022756;
        //isotopicDistribution4[7] = 0.006467;
        //isotopicDistribution4[8] = 0.001654;

        isotopicDistribution5[0] = 0.595410;
        isotopicDistribution5[1] = 1.000000;
        isotopicDistribution5[2] = 0.909613;
        isotopicDistribution5[3] = 0.587263;
        isotopicDistribution5[4] = 0.299194;
        isotopicDistribution5[5] = 0.127213;
        //isotopicDistribution5[6] = 0.046724;
        //isotopicDistribution5[7] = 0.015174;
        //isotopicDistribution5[8] = 0.004431;
        //isotopicDistribution5[9] = 0.001178;

        isotopicDistribution6[0] = 0.484866;
        isotopicDistribution6[1] = 0.953177;
        isotopicDistribution6[2] = 1.000000;
        isotopicDistribution6[3] = 0.737825;
        isotopicDistribution6[4] = 0.427036;
        isotopicDistribution6[5] = 0.205462;
        isotopicDistribution6[6] = 0.085169;
        //isotopicDistribution6[7] = 0.031161;
        //isotopicDistribution6[8] = 0.010239;
        //isotopicDistribution6[9] = 0.003062;

        isotopicDistribution7[0] = 0.368757;
        isotopicDistribution7[1] = 0.828938;
        isotopicDistribution7[2] = 1.000000;
        isotopicDistribution7[3] = 0.851429;
        isotopicDistribution7[4] = 0.570087;
        isotopicDistribution7[5] = 0.317911;
        isotopicDistribution7[6] = 0.152962;
        isotopicDistribution7[7] = 0.065031;
        //isotopicDistribution7[8] = 0.024851;
        //isotopicDistribution7[9] = 0.008647;
        //isotopicDistribution7[10] = 0.002767;

        isotopicDistribution8[0] = 0.295369;
        isotopicDistribution8[1] = 0.745770;
        isotopicDistribution8[2] = 1.000000;
        isotopicDistribution8[3] = 0.939923;
        isotopicDistribution8[4] = 0.691490;
        isotopicDistribution8[5] = 0.422289;
        isotopicDistribution8[6] = 0.221980;
        isotopicDistribution8[7] = 0.102927;
//        isotopicDistribution8[8] = 0.042845;
 //       isotopicDistribution8[9] = 0.016224;
        //isotopicDistribution8[10] = 0.005647;
        //isotopicDistribution8[11] = 0.001821;

        isotopicDistribution9[0] = 0.234948;
        isotopicDistribution9[1] = 0.658283;
        isotopicDistribution9[2] = 0.971762;
        isotopicDistribution9[3] = 1.000000;
        isotopicDistribution9[4] = 0.802247;
        isotopicDistribution9[5] = 0.532683;
        isotopicDistribution9[6] = 0.303778;
        isotopicDistribution9[7] = 0.152559;
        isotopicDistribution9[8] = 0.068695;
  //      isotopicDistribution9[9] = 0.028112;
     //   isotopicDistribution9[10] = 0.010566;
     //   isotopicDistribution9[11] = 0.003678;
     //   isotopicDistribution9[12] = 0.001194;

        isotopicDistribution10[0] = 0.178508;
        isotopicDistribution10[1] = 0.551598;
        isotopicDistribution10[2] = 0.891811;
        isotopicDistribution10[3] = 1.000000;
        isotopicDistribution10[4] = 0.870814;
        isotopicDistribution10[5] = 0.625779;
        isotopicDistribution10[6] = 0.385337;
        isotopicDistribution10[7] = 0.208576;
        isotopicDistribution10[8] = 0.101081;
//        isotopicDistribution10[9] = 0.044468;
//        isotopicDistribution10[10] = 0.017951;
//        isotopicDistribution10[11] = 0.006706;
        //isotopicDistribution10[12] = 0.002335;

        isotopicDistribution11[0] = 0.139841;
        isotopicDistribution11[1] = 0.470844;
        isotopicDistribution11[2] = 0.825470;
        isotopicDistribution11[3] = 1.000000;
        isotopicDistribution11[4] = 0.938111;
        isotopicDistribution11[5] = 0.724608;
        isotopicDistribution11[6] = 0.478745;
        isotopicDistribution11[7] = 0.277647;
        isotopicDistribution11[8] = 0.144001;
        isotopicDistribution11[9] = 0.067736;
//        isotopicDistribution11[10] = 0.029214;
  //      isotopicDistribution11[11] = 0.011655;
    //    isotopicDistribution11[12] = 0.004331;
      //  isotopicDistribution11[13] = 0.001508;

        isotopicDistribution12[0] = 0.110880;
        isotopicDistribution12[1] = 0.404250;
        isotopicDistribution12[2] = 0.764124;
        isotopicDistribution12[3] = 0.994691;
        isotopicDistribution12[4] = 1.000000;
        isotopicDistribution12[5] = 0.825970;
        isotopicDistribution12[6] = 0.582522;
        isotopicDistribution12[7] = 0.360096;
        isotopicDistribution12[8] = 0.198832;
        isotopicDistribution12[9] = 0.099472;
//        isotopicDistribution12[10] = 0.045591;
  //      isotopicDistribution12[11] = 0.019314;
    //    isotopicDistribution12[12] = 0.007618;
      //  isotopicDistribution12[13] = 0.002814;

        isotopicDistribution13[0] = 0.081380;
        isotopicDistribution13[1] = 0.319627;
        isotopicDistribution13[2] = 0.652303;
        isotopicDistribution13[3] = 0.918383;
        isotopicDistribution13[4] = 1.000000;
        isotopicDistribution13[5] = 0.895645;
        isotopicDistribution13[6] = 0.685621;
        isotopicDistribution13[7] = 0.460421;
        isotopicDistribution13[8] = 0.276378;
        isotopicDistribution13[9] = 0.150407;
        isotopicDistribution13[10] = 0.075030;
//        isotopicDistribution13[11] = 0.034611;
  //      isotopicDistribution13[12] = 0.014871;
    //    isotopicDistribution13[13] = 0.005987;
      //  isotopicDistribution13[14] = 0.002270;

        isotopicDistribution14[0] = 0.062581;
        isotopicDistribution14[1] = 0.263702;
        isotopicDistribution14[2] = 0.575328;
        isotopicDistribution14[3] = 0.863482;
        isotopicDistribution14[4] = 1.000000;
        isotopicDistribution14[5] = 0.950817;
        isotopicDistribution14[6] = 0.771506;
        isotopicDistribution14[7] = 0.548469;
        isotopicDistribution14[8] = 0.348161;
        isotopicDistribution14[9] = 0.200189;
        isotopicDistribution14[10] = 0.105433;
        isotopicDistribution14[11] = 0.051316;
//        isotopicDistribution14[12] = 0.023251;
  //      isotopicDistribution14[13] = 0.009867;
    //    isotopicDistribution14[14] = 0.003941;
      //  isotopicDistribution14[15] = 0.001489;

        isotopicDistribution15[0] = 0.049015;
        isotopicDistribution15[1] = 0.220215;
        isotopicDistribution15[2] = 0.510673;
        isotopicDistribution15[3] = 0.812617;
        isotopicDistribution15[4] = 0.995728;
        isotopicDistribution15[5] = 1.000000;
        isotopicDistribution15[6] = 0.855816;
        isotopicDistribution15[7] = 0.640920;
        isotopicDistribution15[8] = 0.428150;
        isotopicDistribution15[9] = 0.258844;
        isotopicDistribution15[10] = 0.143228;
        isotopicDistribution15[11] = 0.073195;
//        isotopicDistribution15[12] = 0.034802;
  //      isotopicDistribution15[13] = 0.015490;
    //    isotopicDistribution15[14] = 0.006487;
      //  isotopicDistribution15[15] = 0.002568;

        isotopicDistribution16[0] = 0.037099;
        isotopicDistribution16[1] = 0.176952;
        isotopicDistribution16[2] = 0.434584;
        isotopicDistribution16[3] = 0.730909;
        isotopicDistribution16[4] = 0.945007;
        isotopicDistribution16[5] = 1.000000;
        isotopicDistribution16[6] = 0.900675;
        isotopicDistribution16[7] = 0.709153;
        isotopicDistribution16[8] = 0.497628;
        isotopicDistribution16[9] = 0.315790;
        isotopicDistribution16[10] = 0.183301;
        isotopicDistribution16[11] = 0.098210;
//        isotopicDistribution16[12] = 0.048934;
  ///      isotopicDistribution16[13] = 0.022815;
     //   isotopicDistribution16[14] = 0.010006;
       // isotopicDistribution16[15] = 0.004146;
        //isotopicDistribution16[16] = 0.001630;

        isotopicDistribution17[0] = 0.028489;
        isotopicDistribution17[1] = 0.143772;
        isotopicDistribution17[2] = 0.372798;
        isotopicDistribution17[3] = 0.660813;
        isotopicDistribution17[4] = 0.899122;
        isotopicDistribution17[5] = 1.000000;
        isotopicDistribution17[6] = 0.945614;
        isotopicDistribution17[7] = 0.780955;
        isotopicDistribution17[8] = 0.574359;
        isotopicDistribution17[9] = 0.381740;
        isotopicDistribution17[10] = 0.231932;
        isotopicDistribution17[11] = 0.130003;
        isotopicDistribution17[12] = 0.067735;
//        isotopicDistribution17[13] = 0.033010;
  //      isotopicDistribution17[14] = 0.015127;
    //    isotopicDistribution17[15] = 0.006548;
      //  isotopicDistribution17[16] = 0.002688;
      //  isotopicDistribution17[17] = 0.001050;
      //  isotopicDistribution17[18] = 0.021473;

        isotopicDistribution18[0] = 0.114423;
        isotopicDistribution18[1] = 0.313596;
        isotopicDistribution18[2] = 0.588027;
        isotopicDistribution18[3] = 0.846974;
        isotopicDistribution18[4] = 0.997820;
        isotopicDistribution18[5] = 1.000000;
        isotopicDistribution18[6] = 0.875691;
        isotopicDistribution18[7] = 0.683169;
        isotopicDistribution18[8] = 0.481829;
        isotopicDistribution18[9] = 0.310749;
        isotopicDistribution18[10] = 0.184951;
        isotopicDistribution18[11] = 0.102348;
        isotopicDistribution18[12] = 0.052989;
//        isotopicDistribution18[13] = 0.025802;
//        isotopicDistribution18[14] = 0.011870;
//        isotopicDistribution18[15] = 0.005179;
//        isotopicDistribution18[16] = 0.002150;

        isotopicDistribution19[0] = 0.016073;
        isotopicDistribution19[1] = 0.090246;
        isotopicDistribution19[2] = 0.260110;
        isotopicDistribution19[3] = 0.512071;
        isotopicDistribution19[4] = 0.773268;
        isotopicDistribution19[5] = 0.953917;
        isotopicDistribution19[6] = 1.000000;
        isotopicDistribution19[7] = 0.915161;
        isotopicDistribution19[8] = 0.745551;
        isotopicDistribution19[9] = 0.548712;
        isotopicDistribution19[10] = 0.369064;
        isotopicDistribution19[11] = 0.228958;
        isotopicDistribution19[12] = 0.132003;
        isotopicDistribution19[13] = 0.071173;
//        isotopicDistribution19[14] = 0.036079;
//        isotopicDistribution19[15] = 0.017273;
//        isotopicDistribution19[16] = 0.007841;
//        isotopicDistribution19[17] = 0.003386;
//        isotopicDistribution19[18] = 0.001395;
    }

    public static double[] getIsotopeModel(double mass) {

        int averagineIndex = (int)(mass/500);
        double[] arr = null;

        switch(averagineIndex) {
            case 0:
                arr = isotopicDistribution0;
                break;
            case 1:
                arr = isotopicDistribution1;
                break;
            case 2:
                arr = isotopicDistribution2;
                break;
            case 3:
                arr = isotopicDistribution3;
                break;
            case 4:
                arr = isotopicDistribution4;
                break;
            case 5:
                arr = isotopicDistribution5;
                break;
            case 6:
                arr = isotopicDistribution6;
                break;
            case 7:
                arr = isotopicDistribution7;
                break;
            case 8:
                arr = isotopicDistribution8;
                break;
            case 9:
                arr = isotopicDistribution9;
                break;
            case 10:
                arr = isotopicDistribution10;
                break;
            case 11:
                arr = isotopicDistribution11;
                break;
            case 12:
                arr = isotopicDistribution12;
                break;
            case 13:
                arr = isotopicDistribution13;
                break;
            case 14:
                arr = isotopicDistribution14;
                break;
            case 15:
                arr = isotopicDistribution15;
                break;
            case 16:
                arr = isotopicDistribution16;
                break;
            case 17:
                arr = isotopicDistribution17;
                break;
            case 18:
                arr = isotopicDistribution18;
                break;

        }

        return arr;
    }
}
