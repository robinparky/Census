package edu.scripps.pms.census.util;

//import com.sun.jndi.ldap.PersistentSearchControl;

import edu.scripps.pms.census.labelFree.LabelfreeMissingPeptideBuilderSplit;
import edu.scripps.pms.census.tandem.MZValues;
import edu.scripps.pms.census.util.io.MzxmlSpectrumReader;

import java.io.*;

import edu.scripps.pms.util.PmsUtil;
import edu.scripps.pms.util.sqlite.spectra.SpectraDB;
import gnu.trove.TDoubleArrayList;

import java.io.RandomAccessFile;

import edu.scripps.pms.census.hash.IndexedFile;
import edu.scripps.pms.census.CensusConstants;
import edu.scripps.pms.census.conf.*;
import edu.scripps.pms.census.model.*;
import edu.scripps.pms.census.model.mrm.*;
import gnu.trove.TIntLongHashMap;
import gnu.trove.TIntDoubleHashMap;
import edu.scripps.pms.census.exception.PrecursorNotFoundException;
import edu.scripps.pms.census.exception.CensusIndexOutOfBoundException;
//import edu.scripps.pms.census.labelFree.LabelfreeMissingPeptideBuilder;
import edu.scripps.pms.census.tmtFilter.TMTUtil;

import java.sql.SQLException;
import java.util.*;
import java.text.*;

import edu.scripps.pms.util.spectrum.*;
import edu.scripps.pms.util.MZXmlHandler;
import edu.scripps.pms.util.stats.GaussianFitting;
//import rpark.statistics.GaussianFitting;
import gnu.trove.TDoubleDoubleHashMap;
import rpark.statistics.Smooth;

/**
 * <p>Title: </p>
 * <p>
 * <p>Description: </p>
 * <p>
 * <p>Copyright: Copyright (c) 2004</p>
 * <p>
 * <p>Company: Yates Lab</p>
 *
 * @author Robin Park
 * @version 1.0
 */
public class CalcUtil {
  public static final char SPACE = ' ';
  public static int ION_START_INDEX = 3;
  public static int ION_START_DIA_LFREE = 2;

  public static final char CARRIAGE_RETURN = '\n';
  public static final char WINDOW_CR = '\r';
  public static final char DOT = '.';
  private static TIntDoubleHashMap precursorMap;
  private static double samplePrecursor;
  private static double refPrecursor;
  private static int chargeState;
  public static double SignalToNoise;
  public static double ReportIonSignalToNoise;
  public static List<MZValues> MZValuesList;
  public static double NoiseLevel;
  protected static double[] massArr = null;
  protected static double[] intArr = null;
  private static boolean hasChargeCol;
  public static double massTolerance = -1;

  public static float steepRatioThreshold;
  private static int[][][] pathArray;
  private static Hashtable<String, IndexedFile> ht;

  private static String pathFileName;
  private static String revSearchFile;
  private static int prevScanNumber;
  private static int scanToSearch;
  public static double[] xyMassArr;
  public static double[] xyIntArr;
  private static int lowbound;
  private static int highbound;


  private static boolean align;
  protected static Configuration conf;

  private static DecimalFormat formatter = new DecimalFormat("0.0000");

  public static void findMRMPeakScan(
    List<MRMPeptideModel> pepList,
    String[] flist,
    Hashtable<String, IndexedFile> ht
  )
    throws PrecursorNotFoundException, IOException, CensusIndexOutOfBoundException, Exception {
    for (Iterator<MRMPeptideModel> pitr = pepList.iterator(); pitr.hasNext(); ) {
      MRMPeptideModel model = pitr.next();
      IndexedFile iFile = ht.get(conf.getFilePath() + model.getFileName());
      double[] rtArr = iFile.getRtArr();
      double mass = model.getParentMass();

      RandomAccessFile rfile = iFile.getFile();
      TIntLongHashMap indexMap = iFile.getMsIndex();
      TDoubleDoubleHashMap rtPrecursorMap = iFile.getRtPrecursorMap();
      //long startPos = indexMap.get( Integer.parseInt(model.getScanNum()) );

      int keyIndex = Arrays.binarySearch(rtArr, Double.parseDouble(model.getRt()));
      int startKeyIndex = Arrays.binarySearch(rtArr, Double.parseDouble(model.getRt()) - model.getRtTolerance());
      int endKeyIndex = Arrays.binarySearch(rtArr, Double.parseDouble(model.getRt()) + model.getRtTolerance());
      keyIndex = fixIndex(keyIndex, rtArr.length);
      model.setKeyIndex(keyIndex);
      startKeyIndex = fixIndex(startKeyIndex, rtArr.length);
      endKeyIndex = fixIndex(endKeyIndex, rtArr.length);

      long maxInt = 0;
      double maxRt = 0;


      double[][] bs = model.getBionArr();
      double[][] ys = model.getYionArr();


      //double[][] br = model.getBionRef();
      //double[][] yr = model.getYionRef();

            /*
            for(int l=0;l<bs.length;l++)
            {
                for(int m=0;m<bs[l].length;m++)
            }
            */

      for (int j = startKeyIndex; j <= endKeyIndex; j++) {

        double pcursor = rtPrecursorMap.get(rtArr[j]);

        if (pcursor != model.getParentMass())
          continue;
        double[] result = readSpectrumByRt(rtArr, j, iFile, model); //.getBionSample(), model.getBionRef(), model.getYionSample(), model.getYionRef());
        if (maxInt < result[1]) {
          maxInt = (long) result[1];
          maxRt = result[0];
        }
      }

      keyIndex = Arrays.binarySearch(rtArr, maxRt);
      keyIndex = fixIndex(keyIndex, rtArr.length);
      model.setKeyIndex(keyIndex);
      model.setStartRt(maxRt);

    }
  }

  public static int fixIndex(int i, int length) {
    if (i < 0) //Cannot find index
      i = -(++i); //Math.abs(++keyIndex);

    if (i >= length)
      i--;

    return i;
  }

  public static String calculateMRMNoId(
    List<MRMPeptideModel> pepList,
    String[] flist,
    //double sPrecursor, double rPrecursor,
    Configuration con,
    Hashtable<String, IndexedFile> ht
  ) throws PrecursorNotFoundException, IOException, CensusIndexOutOfBoundException, Exception

  {
    conf = con;
    int maxWindow = conf.getMaxWindow();
    int margin = conf.getMargin();

    // double[][] result = new double[maxWindow*2+1 + margin*2][4*bioSample.length+3]; //scan #, sample intensity, ref intensity

    findMRMPeakScan(pepList, flist, ht);

    StringBuffer sb = new StringBuffer();

    for (Iterator<MRMPeptideModel> pitr = pepList.iterator(); pitr.hasNext(); ) {
      MRMPeptideModel model = pitr.next();
      IndexedFile iFile = ht.get(conf.getFilePath() + model.getFileName());
      sb.append(peakFindingMRM(model, conf, iFile));
      sb.append("&");
    }



        /*
            precursorMap = iFile.getPrecursorMap();

            int margin = conf.getMargin();
            int steepArea = conf.getSteepArea();
            steepRatioThreshold = conf.getSteepRatioThreshold();
            int maxWindow = conf.getMaxWindow();

            double[][] result = new double[maxWindow*2+1 + margin*2][4*bioSample.length+3]; //scan #, sample intensity, ref intensity

         *
         *
    */

        /*
            IndexedFile iFile = ht.get(conf.getFilePath() + model.getFileName());

            precursorMap = iFile.getPrecursorMap();

            int[] keys = iFile.getKeys();



            SpectrumModel sModel = new SpectrumModel();
            sModel.setHighRes(false);
            sModel.setKeys(keys);
            sModel.setIndex(index);
            sModel.setIFile(iFile);
            sModel.setFile(file);
            sModel.setBioSample(bioSample);
            sModel.setYioSample(yioSample);
            sModel.setBioRef(bioRef);
            sModel.setYioRef(yioRef);


            return peakFindingMSMS(sModel, range, conf, keyIndex);

            */

    return sb.toString();
  }


  //MSMS Labeled cbamberg
  public static String calculateSingleMS2Mass(
    IndexedFile iFile,
    int keyIndex,
    ///int diff,
    double[][] massList,
    Configuration conf,
    int cState,
    double[][] refMassList,
    int scanNum
  ) throws PrecursorNotFoundException, IOException, CensusIndexOutOfBoundException, Exception {
    TIntLongHashMap index = iFile.getMsIndex();
    precursorMap = iFile.getPrecursorMap();

    int[] keys = iFile.getKeys();
    RandomAccessFile file = iFile.getFile();

    double wideTolerance = conf.getMs2WideTolerance();
    double narrowTolerance = conf.getMs2NarrowTolerance();
    double intensityThreshold = conf.getIntensityThreshold();

//        wideTolerance = 0.005;
//        narrowTolerance =0.001;


    if ((conf.getMs2ScanTypeValue().equals(iFile.getScanType(scanNum)) && conf.getScanShift() != 0)) {

      int newScanNum = scanNum + conf.getScanShift();
      String sType = iFile.getScanType(newScanNum);


      boolean isSamePrecursor = TMTUtil.isSamePrecursor(iFile, conf.getScanShift(), scanNum, conf.getMs2ScanTypeValue());
      if (!isSamePrecursor) return null;

      //check if it is cid or hcd
      //if cid, return null;
      //if hcd, check precursor
      //if different, return null;
      //else if same, run quant
      keyIndex = CommonTools.getKeyindex(iFile, newScanNum);
    } else if (conf.isMs2ScanTypeValueUnassigned() && conf.getScanShift() != 0) {
      int newScanNum = scanNum + conf.getScanShift();
      String sType = iFile.getScanType(newScanNum);
      keyIndex = CommonTools.getKeyindex(iFile, newScanNum);

    }


    double[][] arr = readSpectrumOnly(keyIndex, iFile, conf);
    if (null == arr) return null;

    double[] massArr = arr[0];
    double[] intArr = arr[1];

    double[] standardArr = conf.getMs2Standards();
    boolean isFound = false;
    if (null != standardArr) {


      for (double each : standardArr) {
        double tempTolerance = each / 1000000 * narrowTolerance;
        double startMass = each - tempTolerance;
        double endMass = each + tempTolerance;


        int start = 0;
        while (true) {
          if (start >= massArr.length)
            break;

          double temp = each - massArr[start];
          if (temp < 0)
            temp = -temp;

          if (temp < tempTolerance) {

            isFound = true;
          }

          if (massArr[start] > endMass)
            break;

          start++;
        }
      }
    }

    //    System.out.println("============standard peak(s) are not found for Scan number\t" + isFound);

    if (null != standardArr && !isFound) {
      System.out.println("standard peak(s) are not found for Scan number\t" + scanNum + "\tat " + iFile.getFileName());
      return null;
    }


    StringBuffer sb = new StringBuffer();
    for (int i = 0; i < massList.length; i++) {
      for (int j = 0; j < massList[i].length; j++) {
        //double sum = intensitySumForSinglePeak(massArr, intArr, massList[i][j], tolerance);
        //public static double intensitySumForSinglePeak(double[] massArr, double[] intArr, double mass, double tolerance)
        double sumIntensity = 0;
        double smallCount = 0;
        double mass = massList[i][j];

//if(203 != (int)mass) continue;
        double refMass = refMassList[i][j];
        double massDiff = refMass - mass;
        double startMass = mass - mass / 1000000 * wideTolerance;
        double endMass = mass + mass / 1000000 * wideTolerance;

        int start = Arrays.binarySearch(massArr, startMass);
        int end = Arrays.binarySearch(massArr, endMass);

        start = ~start;
        end = ~end;
        end--;
        //modify the code
        //if(start<0)
        //    start = -start -1;

        //if(end<0)
        //   end = -end -2;


        // System.out.println("mass start wide " + mass + " " + startMass  + " "  +endMass  + " " + start + " " + end + " " + massArr[start] + " " + massArr[end]);

        int pairCount = 0;
        int samPeakCount = 0;
        for (int k = start; k <= end; k++) {
          //System.out.println("data\t" + mass + "\t" + massArr[k] + "\t" + intArr[k] + "\t" + mass + " " + massDiff);

          double eachRefMass = massArr[k] + massDiff;
          //double driftedMass = massArr[k]-wideTolerance;
          //       System.out.println(intArr[k]);

          int subRefStart = Arrays.binarySearch(massArr, eachRefMass - eachRefMass / 1000000 * narrowTolerance);
          int subRefEnd = Arrays.binarySearch(massArr, eachRefMass + eachRefMass / 1000000 * narrowTolerance);

          subRefStart = ~subRefStart;
          subRefEnd = ~subRefEnd;

          subRefEnd--;
                      /*
                    if(subRefStart<0) //Cannot find index
                        subRefStart=-(++subRefStart); //Math.abs(++keyIndex);

                    if(subRefStart>=keys.length)
                        subRefStart--;
                    */
          //if(labeledMass<massArr[subStart])
          //    subStart--;


                    /*
                    if(subEnd<0) //Cannot find index
                        subEnd=-(++subEnd); //Math.abs(++keyIndex);

		    subEnd--;  //subtract one more

                    if(subEnd>=keys.length)
                        subEnd--;
                    */

          int subSamStart = Arrays.binarySearch(massArr, massArr[k] - massArr[k] / 1000000 * narrowTolerance);
          ;
          int subSamEnd = Arrays.binarySearch(massArr, massArr[k] + massArr[k] / 1000000 * narrowTolerance);

          subSamStart = ~subSamStart;
          subSamEnd = ~subSamEnd;
          subSamEnd--;
  /*
                    if(subSamStart<0) //Cannot find index
                        subSamStart=-(++subSamStart); //Math.abs(++keyIndex);

                    if(subSamStart>=keys.length)
                        subSamStart--;

                    if(eachRefMass<massArr[subSamStart])
                        subSamStart--;

                    if(subSamEnd<0) //Cannot find index
                        subSamEnd=-(++subSamEnd); //Math.abs(++keyIndex);

		    subSamEnd--;  //subtract one more

                    if(subSamEnd>=keys.length)
                        subSamEnd--;
*/
          double tempSum = 0;

          //        System.out.println("mass==>" + mass + " " + (massArr[subSamStart]-mass)  + " " +  (massArr[subSamEnd]-mass)  + " " +  subSamStart + " " + subSamEnd + " " + massArr[subSamStart]  + " " + massArr[subSamEnd]);
          //        System.out.println("labeled mass ======" + eachRefMass + " " +  (massArr[subRefStart]-eachRefMass)  + " " + (massArr[subRefEnd]-eachRefMass) + " " );

          //  System.out.println("start");
          for (int l = subRefStart; l <= subRefEnd; l++) {
            if (intArr[l] < intensityThreshold) continue;
            tempSum += intArr[l];
          }

          if (tempSum > 0) {
            pairCount++;
          }

          //  System.out.println( "=------ " + pairCount);
          int tempCount = 0;
          for (int l = subSamStart; l <= subSamEnd; l++) {
            // System.out.println(" " + narrowTolerance);
            if (intArr[l] < intensityThreshold) continue;

            sumIntensity += intArr[l];
            samPeakCount++;
            tempCount++;
          }

          if (tempCount > 1)
            smallCount = 2;

//                       System.out.println("==============" + massArr[k] +  " "+ massDiff + " " + sumIntensity);
        }


//                    System.out.println("=====" + smallCount + " " + samPeakCount  + " " + pairCount);
        if (smallCount > 1)
          sb.append(massList[i][j]).append(" ").append("MSC").append(" ");
        else if (samPeakCount > 1)
          sb.append(massList[i][j]).append(" ").append("MP").append(" ");
        else if (samPeakCount == 1 || pairCount == 1) //singleton
          sb.append(massList[i][j]).append(" ").append(sumIntensity).append(" ");
        else if (pairCount < 1)
          sb.append(massList[i][j]).append(" ").append("NP").append(" ");
        else
          sb.append(massList[i][j]).append(" ").append("MP" + pairCount).append(" ");

        //if(sumIntensity>0)        System.exit(0);


      }
    }

    return sb.toString();

  }

  //MSMS Labeled Low res
  public static String calculateMS2Labelfree(
    IndexedFile iFile,
    SpecRange range,
    int keyIndex,
    ///int diff,
    double[][] bioSample,
    double[][] yioSample,
    Configuration conf, int cState
  ) throws PrecursorNotFoundException, IOException, CensusIndexOutOfBoundException, Exception {
//System.out.println("==-");

    TIntLongHashMap index = iFile.getMsIndex();
    precursorMap = iFile.getPrecursorMap();

    int[] keys = iFile.getKeys();
    RandomAccessFile file = iFile.getFile();

    //this condition is not default
    if (conf.getMsmsFragType() < 0 || conf.getMsmsFragType() == conf.AUTOMATIC_FRAGMENT_ION) //data independent
    {
      //hasChargeCol = conf.isChargeColumn();
      //massTolerance = conf.getMassTolerance();
      //int numIsoWindows = conf.getNumOfIsolationWindow();
      //samplePrecursor = sPrecursor;
      //refPrecursor = rPrecursor;
      //chargeState = cState;
      precursorMap = iFile.getPrecursorMap();
      TIntDoubleHashMap scanToRetMap = iFile.getScanRtMap();

      int margin = conf.getMargin();
      int steepArea = conf.getSteepArea();
      steepRatioThreshold = conf.getSteepRatioThreshold();
      int maxWindow = conf.getMaxWindow();

      double[][] result = new double[maxWindow * 2 + 1 + margin * 2][4 * bioSample.length + 3]; //scan #, sample intensity, ref intensity


      SpectrumModel sModel = new SpectrumModel();
      sModel.setHighRes(false);
      sModel.setKeys(keys);
      sModel.setIndex(index);
      sModel.setIFile(iFile);
      sModel.setFile(file);
      sModel.setBioSample(bioSample);
      sModel.setYioSample(yioSample);
      //sModel.setDiff(diff);





      return peakFindingMSMSLabelfree(sModel, range, conf, keyIndex, scanToRetMap);
    }
    //this option is default
    else if (conf.getMsmsFragType() == conf.SPECIFIC_FRAGMENT_ION) {
      //there are two different types.  1. single spectrum 2. multiple spectra

      List<ReportIon> massList = conf.getReportIonList();

      if (conf.getMsmsSpectrumNum() == conf.MSMS_SINGLE_SPECTRUM)  //single spectrum.  No peak finding
      {
        String result = readSpectrumSpecificPeaksWithPurityCorrection(
          keyIndex,
          iFile,
          massList,
          conf
        );

        if (null == result) return null;

        StringBuffer sb = new StringBuffer();


        for (Iterator<ReportIon> mItr = massList.iterator(); mItr.hasNext(); ) {
          ReportIon reportIon = mItr.next();
          double mass = reportIon.getMass();

          sb.append(mass).append(" ");
        }

        sb.setCharAt(sb.length() - 1, ',');

        sb.append(result);
                /*
                for(double d : result)
                {
                    sb.append((long)d).append(" ");
                }*/

        return sb.toString();

      }

    }

    return null;

  }

  //MSMS Labeled Low res, or TMT isobaric data analysis
  public static String calculateMS2Mass(
    IndexedFile iFile,
    SpecRange range,
    int keyIndex,
    ///int diff,
    double[][] bioSample, double[][] bioRef,
    double[][] yioSample, double[][] yioRef,
    //double sPrecursor, double rPrecursor,
    Configuration conf, int cState
  ) throws PrecursorNotFoundException, IOException, CensusIndexOutOfBoundException, Exception {
//System.out.println("==-");

    TIntLongHashMap index = iFile.getMsIndex();
    precursorMap = iFile.getPrecursorMap();

    int[] keys = iFile.getKeys();
    RandomAccessFile file = iFile.getFile();

    //this condition is not default
    if (conf.getMsmsFragType() < 0 || conf.getMsmsFragType() == conf.AUTOMATIC_FRAGMENT_ION) //data independent
    {
      //hasChargeCol = conf.isChargeColumn();
      //massTolerance = conf.getMassTolerance();
      //int numIsoWindows = conf.getNumOfIsolationWindow();
      //samplePrecursor = sPrecursor;
      //refPrecursor = rPrecursor;
      //chargeState = cState;
      precursorMap = iFile.getPrecursorMap();

      int margin = conf.getMargin();
      int steepArea = conf.getSteepArea();
      steepRatioThreshold = conf.getSteepRatioThreshold();
      int maxWindow = conf.getMaxWindow();

      double[][] result = new double[maxWindow * 2 + 1 + margin * 2][4 * bioSample.length + 3]; //scan #, sample intensity, ref intensity

      SpectrumModel sModel = new SpectrumModel();
      sModel.setHighRes(false);
      sModel.setKeys(keys);
      sModel.setIndex(index);
      sModel.setIFile(iFile);
      sModel.setFile(file);
      sModel.setBioSample(bioSample);
      sModel.setYioSample(yioSample);
      sModel.setBioRef(bioRef);
      sModel.setYioRef(yioRef);
      //sModel.setDiff(diff);

      return peakFindingMSMS(sModel, range, conf, keyIndex);
    }
    //this option is default
    else if (conf.getMsmsFragType() == conf.SPECIFIC_FRAGMENT_ION) {
      //there are two different types.  1. single spectrum 2. multiple spectra

      List<ReportIon> massList = conf.getReportIonList();

      if (conf.getMsmsSpectrumNum() == conf.MSMS_SINGLE_SPECTRUM)  //single spectrum.  No peak finding
      {
        String result = readSpectrumSpecificPeaksWithPurityCorrection(
          keyIndex,
          iFile,
          massList,
          conf
        );

        if (null == result) return null;

        StringBuffer sb = new StringBuffer();


        for (Iterator<ReportIon> mItr = massList.iterator(); mItr.hasNext(); ) {
          ReportIon reportIon = mItr.next();
          double mass = reportIon.getMass();

          sb.append(mass).append(" ");
        }

        sb.setCharAt(sb.length() - 1, ',');

        sb.append(result);
                /*
                for(double d : result)
                {
                    sb.append((long)d).append(" ");
                }*/

        return sb.toString();

      } /*else if(conf.getMsmsSpectrumNum() == conf.MSMS_MULTIPLE_SPECTRA) //multiple spectra.  Peak finding
            {
                int margin = conf.getMargin();
                int steepArea = conf.getSteepArea();
                steepRatioThreshold = conf.getSteepRatioThreshold();
                int maxWindow = conf.getMaxWindow();

                double[][] result = new double[maxWindow*2+1 + margin*2][massList.size()+3]; //scan #, sample intensity, ref intensity

                SpectrumModel sModel = new SpectrumModel();
                sModel.setHighRes(false);
                sModel.setKeys(keys);
                sModel.setIndex(index);
                sModel.setIFile(iFile);
                sModel.setFile(file);
                sModel.setMsmsSpecificMassList(massList);

                return peakFindingMSMSMultipleSpecific(sModel, range, conf, keyIndex);
            }*/

    }

    return null;

  }

  protected static String buildDDResult(int peakStart, int peakEnd, int leftIndex, int rightIndex, double[][] result) {
    StringBuffer sb = new StringBuffer();
    sb.append("P ").append((long) result[peakStart][0]).append(" ").append((long) result[peakEnd][0]).append(";");

    for (int i = leftIndex + 1; i < rightIndex - 1; i++) {
      for (int j = 0; j < result[i].length - 2; j++) {
        sb.append((long) result[i][j]).append(" ");
      }

      sb.append((double) result[i][result[i].length - 2]).append(" ");
      sb.append((double) result[i][result[i].length - 1]).append(";");

      //sb.append((int)result[i][0]).append(" ").append((long)result[i][1]).append(" ").append((long)result[i][2]).append(";");
    }

    return sb.toString();
  }

  protected static String buildDDResultPeakFinding(int peakStart, int peakEnd, int leftIndex, int rightIndex, double[][] result, int identifiedIndex, DisplayData.DisplayChroData chroData) {
    if (conf.isSmooth()) {
      long[] testArr = new long[result.length];
      long[] testArr2 = new long[result.length];
      for (int i = 0; i < testArr.length; i++) {
        testArr[i] = (long) result[i][1];
      }
      for (int i = 0; i < testArr2.length; i++) {
        testArr2[i] = (long) result[i][2];
      }

      //for(int i=0;i<testArr2.length;i++){
      //    System.out.println(testArr[i]+"\t"+testArr2[i]);
      //}
      double[] smoothArr = Smooth.smoothAsDouble(testArr, LabelfreeMissingPeptideBuilderSplit.SMOOTH_WINDOW_SIZE);
      double[] smoothArr2 = Smooth.smoothAsDouble(testArr2, LabelfreeMissingPeptideBuilderSplit.SMOOTH_WINDOW_SIZE);
      // for(int i=0;i<testArr.length;i++){
      //    System.out.println(smoothArr[i]+"\t"+smoothArr2[i]);
      // }


      for (int i = 0; i < testArr.length; i++) {
        result[i][1] = smoothArr[i];
      }
      for (int i = 0; i < testArr2.length; i++) {
        result[i][2] = smoothArr2[i];
      }
    }

    Range range = GaussianFitting.getGaussianPeakRangeIndex(result, peakStart, peakEnd, identifiedIndex);


    StringBuffer sb = new StringBuffer();
    //  sb.append("P ").append((long)result[peakStart][0]).append(" ").append((long)result[peakEnd][0]).append(";");
    sb.append("P ").append((long) range.getLowBound()).append(" ").append((long) range.getHighBound()).append(";");
    setHighbound((int) range.getHighBound());
    setLowbound((int) range.getLowBound());
    long maxIntensity = 0;

    for (int i = leftIndex + 1; i < rightIndex - 1; i++) {
      for (int j = 0; j < result[i].length - 2; j++) {
        if (chroData != null) {
          double d = result[i][j];
          if (d > maxIntensity) maxIntensity = (long) d;
          if (j == 1) {
            chroData.getData1().add(new DisplayData.DisplayChroDataXY((int) result[i][0], (long) result[i][1]));
          } else if (j == 2) {
            chroData.getData2().add(new DisplayData.DisplayChroDataXY((int) result[i][0], (long) result[i][2]));
          }
        }

        sb.append((long) result[i][j]).append(" ");
      }

      sb.append((double) result[i][result[i].length - 2]).append(" ");
      sb.append((double) result[i][result[i].length - 1]).append(";");

      //sb.append((int)result[i][0]).append(" ").append((long)result[i][1]).append(" ").append((long)result[i][2]).append(";");
    }
    if (chroData != null) {
      chroData.setMaxIntensity(maxIntensity);
      chroData.smooth();
    }
    return sb.toString();
  }


  //build non labeling result
  //expSize : number of experiments
  private static String buildNLResult(int peakStart, int peakEnd, int leftIndex, int rightIndex, double[][] result, int expSize) {
    return buildNLResult(peakStart, peakEnd, leftIndex, rightIndex, result, -1, expSize);
  }

  private static String buildNLResult(int peakStart, int peakEnd, int leftIndex, int rightIndex, double[][] result, int ionLength, int expSize) {
    StringBuffer resultBuffer = new StringBuffer();
    resultBuffer.append("P ").append((int) result[peakStart][0]).append(" ").append((int) result[peakEnd][0]).append(";");

    if (conf.getQuantLevel() == 1) {
      if (leftIndex < 0) leftIndex = 0;
      if (rightIndex >= result.length) rightIndex = result.length - 1;

      for (int i = leftIndex + 1; i < rightIndex - 1; i++) {
        for (int j = 0; j < result[i].length; j++) {
          resultBuffer.append((int) result[i][j]).append(" ");

          //System.out.print(result[i][j] + "\t");
        }

        //System.out.println(result[231].length + "\t" + result[i].length + " \t" + i);
        //          resultBuffer.append((int)result[i][0]).append(" "); //spec num
        //        resultBuffer.append((long)result[i][1]).append(" "); //sample intensity
        //      resultBuffer.append((long)result[i][2]).append(","); //ref intensity

        resultBuffer.append(";");
      }


    } else if (conf.getQuantLevel() == 2) {
      int ionMaxIndex = ionLength - 1;

      if (leftIndex < 0) leftIndex = 0;
      if (rightIndex >= result.length) rightIndex = result.length - 1;

      for (int i = leftIndex + 1; i < rightIndex - 1; i++) {
        for (int j = 0; j < expSize - 1; j++) {
          resultBuffer.append((long) result[i][j]).append(" ");
//		    System.out.println(j);
        }

        resultBuffer.append((long) result[i][expSize - 1]).append(",");
//		    System.out.println(expSize-1);

        for (int j = expSize; j < 2 * expSize - 1; j++) {
          resultBuffer.append((long) result[i][j]).append(" ");
//		    System.out.println(j);
        }

        resultBuffer.append((long) result[i][2 * expSize - 1]).append(",");
//		System.out.println(2*expSize-1);

        for (int j = 0; j < expSize * 2; j++) {
          for (int k = 6; k < ionLength + 6; k++) {
            resultBuffer.append(result[i][k + j * (ionMaxIndex + 1)]).append(" ");
//			System.out.println(k+j*(ionMaxIndex+1));
          }

          resultBuffer.append(result[i][ionLength + 6]).append(",");

        }

        resultBuffer.setCharAt(resultBuffer.length() - 1, ';');

//	       System.out.println(resultBuffer.toString());

      }

    }

    return resultBuffer.toString();

  }

  private static String buildDIResultLabelfree(int peakStart, int peakEnd, int leftIndex, int rightIndex, double[][] result, int ionLength, DisplayData.DisplayChroData chrodata) {
    StringBuffer resultBuffer = new StringBuffer();
    resultBuffer.append("P ").append((int) result[peakStart][0]).append(" ").append((int) result[peakEnd][0]).append(";");

/*
    for(int i=0;i<result.length;i++) {
      for(int j=0;j<result[i].length;j++) {
        System.out.print(result[i][j] + " ");
      }

      System.out.println("");
    } */

    int ionMaxIndex = ionLength - 1;

    if (leftIndex < 0) leftIndex = 0;
    if (rightIndex >= result.length) rightIndex = result.length - 1;

    for (int i = leftIndex + 1; i < rightIndex - 1; i++) {

      //System.out.println("==\t" + result[i][0]);


      for(int j=0;j<result[i].length;j++) {
        resultBuffer.append(result[i][j]).append(" "); //spec num
        //resultBuffer.append((long) result[i][j]).append(" "); //spec num

 //       System.out.println(i + " " + j  + "==\t" + result[i][j]);
      }

      resultBuffer.append(";");



    }

    return resultBuffer.toString();
  }

  public static String buildDIResult(int peakStart, int peakEnd, int leftIndex, int rightIndex, double[][] result, int ionLength, DisplayData.DisplayChroData chrodata) {
    StringBuffer resultBuffer = new StringBuffer();
    resultBuffer.append("P ").append((int) result[peakStart][0]).append(" ").append((int) result[peakEnd][0]).append(";");

    int ionMaxIndex = ionLength - 1;

    if (leftIndex < 0) leftIndex = 0;
    if (rightIndex >= result.length) rightIndex = result.length - 1;

    for (int i = leftIndex + 1; i < rightIndex - 1; i++) {
      resultBuffer.append((int) result[i][0]).append(" "); //spec num
      resultBuffer.append((long) result[i][1]).append(" "); //sample intensity
      resultBuffer.append((long) result[i][2]).append(","); //ref intensity

      for (int j = 0; j < 3; j++) //3 ions b_sample, b_ref, y_sample
      {
        for (int k = 3; k < ionMaxIndex + 3; k++) {
          resultBuffer.append(result[i][k + j * (ionMaxIndex + 1)]).append(" ");
        }

        resultBuffer.append(result[i][3 + ionMaxIndex + j * (ionMaxIndex + 1)]).append(",");
      }

      //y_ref ion. do this separately due to delimitation
      for (int k = 3; k < ionMaxIndex + 3; k++) {
        resultBuffer.append(result[i][k + ionMaxIndex * 3 + 3]).append(" ");
      }

      resultBuffer.append(result[i][6 + ionMaxIndex + ionMaxIndex * 3]).append(";");
    }

    return resultBuffer.toString();
  }

  private static String buildMRMResult(int peakStart, int peakEnd, int leftIndex, int rightIndex, double[][] result, int ionLength) {
    StringBuffer resultBuffer = new StringBuffer();
    resultBuffer.append("P ").append(result[peakStart][0]).append(" ").append(result[peakEnd][0]).append(";");

    int ionMaxIndex = ionLength - 1;

    if (leftIndex < 0) leftIndex = 0;
    if (rightIndex >= result.length) rightIndex = result.length - 1;

    for (int i = leftIndex + 1; i < rightIndex - 1; i++) {
      resultBuffer.append((double) result[i][0]).append(" "); //spec num
      resultBuffer.append((long) result[i][1]).append(" "); //sample intensity
      resultBuffer.append((long) result[i][2]).append(","); //ref intensity

      for (int j = 0; j < 3; j++) //3 ions b_sample, b_ref, y_sample
      {
        for (int k = 3; k < ionMaxIndex + 3; k++) {
          //	    System.out.println(i + " " + result.length + " " + (k+j*(ionMaxIndex+1)) + " " + result[i].length);
          resultBuffer.append(result[i][k + j * (ionMaxIndex + 1)]).append(" ");
        }

        resultBuffer.append(result[i][3 + ionMaxIndex + j * (ionMaxIndex + 1)]).append(",");
      }

      //y_ref ion. do this separately due to delimitation
      for (int k = 3; k < ionMaxIndex + 3; k++) {
        resultBuffer.append(result[i][k + ionMaxIndex * 3 + 3]).append(" ");
      }

      resultBuffer.append(result[i][6 + ionMaxIndex + ionMaxIndex * 3]).append(";");
    }

    return resultBuffer.toString();
  }

  public static IrisDataModel readSpectrumPeaks(IndexedFile iFile, int keyIndex, int[] keys) throws IOException {

    //double[] result = new double[2];
    RandomAccessFile rfile = iFile.getFile();
    TIntLongHashMap indexMap = iFile.getMsIndex();

    long startPos = indexMap.get(keys[keyIndex]);

    int byteSize = -1;
    byte[] bytes;

    char ch;
    int pos;

    rfile.seek(startPos);

    //if( (keyIndex+1)>=keys.length )
    if ((keyIndex + 1) >= keys.length)
      byteSize = (int) (rfile.length() - startPos);

    if (indexMap.get(keys[keyIndex + 1]) > 0)
      byteSize = (int) (indexMap.get(keys[keyIndex + 1]) - startPos);

    bytes = new byte[byteSize];

    rfile.readFully(bytes);

    pos = 0;

    ch = (char) bytes[pos];

    //        pos++;
    //Remove Z, S, I, D lines
    //System.out.println("sample");
    while ((ch = (char) bytes[pos]) == 'Z' || ch == 'S' || ch == 'I' || ch == 'D') {

      while (ch != CARRIAGE_RETURN) {
        pos++;
        ch = (char) bytes[pos];
      }

      pos++;
    }

    int arrSize = 0;
    for (int j = pos; j < byteSize; j++) {
      if (CARRIAGE_RETURN == (char) bytes[j])
        arrSize++;

    }

    double[][] result = parseSpectra(bytes);
    double[] massArr = result[0];
    double[] intArr = result[1];

    //pos++;
    StringBuilder mass = new StringBuilder(10);
    StringBuilder intensity = new StringBuilder(10);
    intensity.append('0');


    return new IrisDataModel(massArr, intArr);
  }


  public static String getPrecursorMass(TIntLongHashMap index, int i, long sampleIndex, RandomAccessFile file, int[] keys)
    throws IOException {

    file.seek(sampleIndex);

    int increIndex = 1;
    int byteSize;
    byte[] bytes;

    while (true) {
      if ((i + increIndex) >= keys.length) {
        byteSize = (int) (file.length() - sampleIndex);

        break;
      }

      if (index.get(keys[i + increIndex]) > 0) {
//System.out.println("11--------------" + i + " " + keys[i] + " " + (int)(index.get(keys[i])) + " " + sampleIndex);
        byteSize = (int) (index.get(keys[i]) - sampleIndex);
//System.out.println("22-------------->>" + byteSize + " " + i + " " + keys[i+increIndex]  + " " + index.get(keys[i]));
        break;
      }

      increIndex++;
    }


    //byteSize = (int)(index.get(keys[i+1]) - sampleIndex);
    bytes = new byte[byteSize];

    file.readFully(bytes);
    int pos = 0;

    char ch = (char) bytes[pos];

    //        pos++;
    //Remove Z, S, I, D lines
    //System.out.println("sample");
    while ((ch = (char) bytes[pos]) == 'S') {
      StringBuffer tmpsb = new StringBuffer();
      while (ch != CARRIAGE_RETURN) {
        pos++;
        ch = (char) bytes[pos];
        tmpsb.append(ch);

      }

      String tmpstr = tmpsb.toString();
      String[] arr = tmpstr.split("\t");
      if (null != arr[3]) return arr[3];

      pos++;
    }

    return null;

  }


  public static double[] readSpectrumLabelfree (
    int[] keys,
    int i,
    TIntLongHashMap index,
    Object ofile,
    double[][] bioSample,
    double[][] yioSample,
    Configuration conf,
    TIntDoubleHashMap scanToRetMap
  ) throws IOException {
    int ionLength = bioSample.length;



    int ionSize = bioSample[0].length;

    double[] result = new double[bioSample.length * ionSize * 2 + ION_START_DIA_LFREE+1];

    //if( (i+diff)>=keys.length )
    //    return result;

    RandomAccessFile file = (RandomAccessFile) ofile;

    int byteSize;
    byte[] bytes;
    char ch;
    int pos;

    if (i < 0 || i >= keys.length) {
      return result;
    }

    long sampleIndex = index.get(keys[i]);

    //System.exit(1);

    StringBuilder mass, intensity;

    if (sampleIndex < 0) {
      return result;
    }

    result[0] = keys[i]; //scan
    result[1] = scanToRetMap.get(keys[i]); //ret time


    /**********************************
     * Go to Spectrum for sample, and calculate intensity
     **********************************/
    file.seek(sampleIndex);

    int increIndex = 1;

    while (true) {
      if ((i + increIndex) >= keys.length) {
        byteSize = (int) (file.length() - sampleIndex);

        break;
      }

      if (index.get(keys[i + increIndex]) > 0) {
        byteSize = (int) (index.get(keys[i + increIndex]) - sampleIndex);
        break;
      }

      increIndex++;
    }

    //byteSize = (int)(index.get(keys[i+1]) - sampleIndex);
    bytes = new byte[byteSize];

    file.readFully(bytes);
    pos = 0;

    ch = (char) bytes[pos];

    //        pos++;
    //Remove Z, S, I, D lines
    //System.out.println("sample");
    while ((ch = (char) bytes[pos]) == 'Z' || ch == 'S' || ch == 'I' || ch == 'D') {
      while (ch != CARRIAGE_RETURN) {
        pos++;
        ch = (char) bytes[pos];
      }

      pos++;
    }

    int arrSize = 0;
    for (int j = pos; j < byteSize; j++) {
      if (CARRIAGE_RETURN == (char) bytes[j])
        arrSize++;

    }

    double[][] specResult = parseSpectra(bytes);
    double[] massArr = specResult[0];
    double[] intArr = specResult[1];

    int sampleSum = 0;

//        int startMassIndex = 3*chargeState-2;
//        int endMassIndex = startMassIndex+1;


    //boolean tempBool=false;


    int indexCount=3;
    for (int k = 0; k < bioSample.length; k++) {
      for (double d : bioSample[k]) {
        double tolerance = conf.getMassToleranceInPPM() * d / 1000000;

        result[indexCount] = intensitySum(massArr, intArr, (d - tolerance), (d + tolerance));
        sampleSum += result[indexCount];

        indexCount++;
      }


    }

    for (int k = 0; k < yioSample.length; k++) {
      for (double d : yioSample[k]) {
        double tolerance = conf.getMassToleranceInPPM() * d / 1000000;

        result[indexCount] = intensitySum(massArr, intArr, d - tolerance, d + tolerance);
        sampleSum += result[indexCount];

        indexCount++;
      }
    }




    result[2] = sampleSum;



    return result;
  }


  public static double[] readSpectrum(
    int[] keys,
    int i,
    TIntLongHashMap index,
    int refIndexKey,
    Object ofile,
    double[][] bioSample, double[][] bioRef,
    double[][] yioSample, double[][] yioRef
  ) throws IOException {
    //int[] result = new int[3];
    //ION_START_INDEX == 3
    int ionLength = bioSample.length;
    double[] result = new double[ionLength * 4 + ION_START_INDEX];

    //if( (i+diff)>=keys.length )
    //    return result;

    RandomAccessFile file = (RandomAccessFile) ofile;

    int byteSize;
    byte[] bytes;
    char ch;
    int pos;

    if (i < 0 || i >= keys.length) {
      return result;
    }

    long sampleIndex = index.get(keys[i]);

    //System.exit(1);

    StringBuilder mass, intensity;

    if (sampleIndex < 0) {
      return result;
    }

    result[0] = keys[i];


    /**********************************
     * Go to Spectrum for sample, and calculate intensity
     **********************************/
    file.seek(sampleIndex);

    int increIndex = 1;

    while (true) {
      if ((i + increIndex) >= keys.length) {
        byteSize = (int) (file.length() - sampleIndex);

        break;
      }

      if (index.get(keys[i + increIndex]) > 0) {
        byteSize = (int) (index.get(keys[i + increIndex]) - sampleIndex);
        break;
      }

      increIndex++;
    }

    //byteSize = (int)(index.get(keys[i+1]) - sampleIndex);
    bytes = new byte[byteSize];

    file.readFully(bytes);
    pos = 0;

    ch = (char) bytes[pos];

    //        pos++;
    //Remove Z, S, I, D lines
    //System.out.println("sample");
    while ((ch = (char) bytes[pos]) == 'Z' || ch == 'S' || ch == 'I' || ch == 'D') {
      while (ch != CARRIAGE_RETURN) {
        pos++;
        ch = (char) bytes[pos];
      }

      pos++;
    }

    int arrSize = 0;
    for (int j = pos; j < byteSize; j++) {
      if (CARRIAGE_RETURN == (char) bytes[j])
        arrSize++;

    }

    double[][] specResult = parseSpectra(bytes);
    double[] massArr = specResult[0];
    double[] intArr = specResult[1];

    int sampleSum = 0;

//        int startMassIndex = 3*chargeState-2;
//        int endMassIndex = startMassIndex+1;

    //boolean tempBool=false;
    for (int k = 0; k < bioSample.length; k++) {

      result[ION_START_INDEX + k] = intensitySum(massArr, intArr, bioSample[k][1], bioSample[k][2]);
      result[ION_START_INDEX + ionLength + k] = intensitySum(massArr, intArr, yioSample[k][1], yioSample[k][2]);

      sampleSum += result[ION_START_INDEX + k];
      sampleSum += result[ION_START_INDEX + bioSample.length + k];

    }

    result[1] = sampleSum;

    if (refIndexKey < 0)
      return result;


    long refIndex = index.get(keys[refIndexKey]);

    if (refIndex < 0)
      return result;


    file.seek(refIndex);

    increIndex = 1;

    while (true) {
      //if( (i+increIndex+diff)>=keys.length )
      if ((increIndex + refIndexKey) >= keys.length) {
        byteSize = (int) (file.length() - refIndex);
        break;
      }

      //if(index.get(keys[i+increIndex+diff])>0)
      if (index.get(keys[refIndexKey + increIndex]) > 0) {
        byteSize = (int) (index.get(keys[refIndexKey + increIndex]) - refIndex);
        break;
      }

      increIndex++;
    }

    bytes = new byte[byteSize];

    file.readFully(bytes);


    specResult = parseSpectra(bytes);
    massArr = specResult[0];
    intArr = specResult[1];

    pos = 0;
    ch = (char) bytes[pos];

        /*
        //Remove S line
        while( ch != CARRIAGE_RETURN )
        {
        pos++;
        ch = (char)bytes[pos];
        }
         */
    //        pos++;
    //Remove Z lines
    //System.out.println("reference");
    while ((ch = (char) bytes[pos]) == 'Z' || ch == 'S' || ch == 'I' || ch == 'D') {
      while (ch != CARRIAGE_RETURN) {
        pos++;
        ch = (char) bytes[pos];
      }

      pos++;
    }

//        System.out.println("end print ch");


    //            refPeakList.setMassList(massArr);

    int refSum = 0;

    for (int k = 0; k < bioSample.length; k++) {
      result[ION_START_INDEX + ionLength * 2 + k] = intensitySum(massArr, intArr, bioRef[k][1], bioRef[k][2]);
      result[ION_START_INDEX + ionLength * 3 + k] = intensitySum(massArr, intArr, yioRef[k][1], yioRef[k][2]);

      refSum += result[ION_START_INDEX + ionLength * 2 + k];
      refSum += result[ION_START_INDEX + ionLength * 3 + k];

    }

    result[2] = refSum;

    return result;
  }

  public static String readSpectrumSpecificPeaksWithPurityCorrection(
    int i,
    IndexedFile iFile,
    List massList,
    Configuration conf
  ) throws IOException {


    RandomAccessFile file = iFile.getFile();
    TIntLongHashMap index = iFile.getMsIndex();

    int[] keys = iFile.getKeys();
    StringBuffer result = new StringBuffer();

    int byteSize;
    byte[] bytes;
    char ch;
    int pos;

    if (i < 0)
      return null;

    long sampleIndex = index.get(keys[i]);

    StringBuilder mass, intensity;


    if (sampleIndex < 0) {
      return null;
    }

    //System.out.println("============" + );
    if (conf.getScanShift() != 0) {


      int tmpIndex = i + conf.getScanShift();
      if (tmpIndex < 0 || tmpIndex >= keys.length) return null;

//System.out.println("============" + index.get(keys[i+conf.getScanShift()]));
//System.out.println("===" + i  + " " + keys[i] + " "  + keys[i+1] + " "  + index.get(keys[i+1]) + " " + conf.getScanShift());
//System.out.println("===" + i + " " + (i+conf.getScanShift()));
      //	String scanShiftPrecursor = getPrecursorMass( index, i+1, index.get(keys[i]), file, keys );
      //	String currentPrecursor = getPrecursorMass( index, i+1-conf.getScanShift(), index.get(keys[i-conf.getScanShift()]), file, keys );
    }

    result.append(keys[i]).append(" ");

    file.seek(sampleIndex);

    int increIndex = 1;

    while (true) {
      if ((i + increIndex) >= keys.length) {
        byteSize = (int) (file.length() - sampleIndex);

        break;
      }

      if (index.get(keys[i + increIndex]) > 0) {
        byteSize = (int) (index.get(keys[i + increIndex]) - sampleIndex);
        break;
      }

      increIndex++;
    }

    //byteSize = (int)(index.get(keys[i+1]) - sampleIndex);
    bytes = new byte[byteSize];

    file.readFully(bytes);

//for(byte b:bytes)
//System.out.print((char)b);

    pos = 0;

    ch = (char) bytes[pos];
    //Remove Z, S, I, D lines
    while ((ch = (char) bytes[pos]) == 'Z' || ch == 'S' || ch == 'I' || ch == 'D') {
      while (ch != CARRIAGE_RETURN) {
        pos++;
        ch = (char) bytes[pos];
      }

      pos++;
    }

    double[][] specResult = parseSpectra(bytes);
    double[] massArr = specResult[0];
    double[] intArr = specResult[1];


    //pos++;
    mass = new StringBuilder(10);
    intensity = new StringBuilder(10);
    intensity.append('0');

    int tmpIndex = 0;

    double[] correctIntArr = new double[massList.size()];
    double[] intSumArr = new double[massList.size()];

    //  List<MZValues> mzList = new ArrayList<>();
    for (Iterator<ReportIon> itr = massList.iterator(); itr.hasNext(); ) {
      //  MZValues mzValues = new MZValues();

      ReportIon ri = itr.next();
      double d = ri.getMass();
      //  mzValues.setMzHeader("m/z_"+Double.toString(d)+"_int");
      //System.out.println("===========" + ri.getcPlusList());
      //System.out.println("===========" + ri.getcMinusList());
      //    mzList.add(mzValues);
      double tolerance = conf.getMsmsSpecificTolerance();
      if (conf.isMsmsTolerancePpm())
        tolerance = PmsUtil.getMZFromPpm((float) d, (float) tolerance);

      intSumArr[tmpIndex] = intensitySumForSinglePeak(massArr, intArr, d, tolerance);

      //for(double each:correctIntArr)
      //    System.out.print(each +"\t");
      //  System.out.println("");

      if (tmpIndex > 0 && ri.getcMinusList().size() > 0) {
        double beforeInt = 0;
        beforeInt = intSumArr[tmpIndex] * 0.01 * ri.getcMinusList().get(0);

        correctIntArr[tmpIndex - 1] -= beforeInt;
        //  System.out.println("b\t" + beforeInt + " " + sum + " " + ri.getcMinusList().get(0));

      }
      if (tmpIndex < correctIntArr.length - 1 && ri.getcPlusList().size() > 0) {
        double afterInt = 0;
        afterInt = intSumArr[tmpIndex] * 0.01 * ri.getcPlusList().get(0);

        correctIntArr[tmpIndex + 1] -= afterInt;
        //System.out.println("a\t" + afterInt + " " + sum + " " + ri.getcPlusList().get(0));
      }

      //for(double each:correctIntArr)
      //    System.out.print(each +"\t");
      //  System.out.println("");

      //result.append((long)sum).append(" ");

      tmpIndex++;
    }

    double[] correctedIntSumArr = new double[correctIntArr.length];

    for (int ii = 0; ii < correctIntArr.length; ii++) {
      correctedIntSumArr[ii] = intSumArr[ii] + correctIntArr[ii];

      if (correctedIntSumArr[ii] <= 0) {
        result.append("0 ");
        //     mzList.get(ii).setMzValue(Double.toString(0));
      } else {
        result.append((long) correctedIntSumArr[ii]).append(" ");
        //       mzList.get(ii).setMzValue(Integer.toString((int)correctedIntSumArr[ii]));

      }

      //   System.out.println("=================" + (correctedIntSumArr[ii]<=0) + "\t" + (long)correctedIntSumArr[ii] +"\t"+ correctedIntSumArr[ii] + "\t" + intSumArr[ii] +"\t"+ correctIntArr[ii]);
      //   result.append((long)correctedIntSumArr[ii]).append(" ");
    }

    result.append(",");
    Arrays.sort(intArr);
    int bottom10 = (int) (intArr.length * 0.1);
    double sum = 0;
    for (int ii = 0; ii < bottom10; ii++) {
      sum += intArr[ii];
    }
    double avg = sum / bottom10;
    if (sum == 0) avg = 1;
    NoiseLevel = avg;

    double maxIntensitySum = Double.MIN_VALUE;
        /*if(Double.isInfinite(avg) || Double.isNaN(avg))
        {
            for(byte b:bytes)
            {
                System.out.print((char)b);
            }
            System.out.println(avg);
        }*/

    for (double each : intSumArr) {
      if (each > maxIntensitySum) maxIntensitySum = each;
      result.append((long) each).append(" ");
    }

    //      MZValuesList = mzList;
    SignalToNoise = maxIntensitySum / avg;
    // SignalToNoise =Math.log10(maxIntensitySum/avg);
    SignalToNoise = Double.isInfinite(SignalToNoise) ? -10 : SignalToNoise;
   /* if(Double.compare(SignalToNoise,0)==0)
    {
      System.out.println(">>>>DEBUG: first peak: "+massArr[0]+" Scan: "+keys[i]);
    }*/
    double noise = maxIntensitySum / SignalToNoise;
    List<Double> reportIonSNList = new ArrayList<>();
    double sumSN = 0;
    for (double each : intSumArr) {
      sumSN += each / noise;
      reportIonSNList.add(each / noise);
    }
    double reportIonSN = sumSN / intSumArr.length;
    ReportIonSignalToNoise = reportIonSN;
    return result.toString();

  }

/*
    public static String readSpectrumSpecificPeaksWithPurityCorrection(
            int i,
            IndexedFile iFile,
            List massList,
	    Configuration conf
            ) throws IOException
    {

	double tolerance = conf.getMsmsSpecificTolerance();
	RandomAccessFile file = iFile.getFile();
        TIntLongHashMap index = iFile.getMsIndex();

        int[] keys = iFile.getKeys();
        StringBuffer result = new StringBuffer();

        int byteSize;
        byte[] bytes;
        char ch;
        int pos;

        if(i<0)
            return null;

        long sampleIndex = index.get(keys[i]);

        StringBuilder mass, intensity;


        if(sampleIndex<0)
        {
            return null;
        }

        //System.out.println("============" + );
	if(conf.getScanShift()!=0) {


                int tmpIndex = i+conf.getScanShift();
		if(tmpIndex<0 || tmpIndex>=keys.length) return null;

	}

        result.append(keys[i]).append(" ");

        file.seek(sampleIndex);

        int increIndex=1;

        while(true)
        {
            if( (i+increIndex)>=keys.length )
            {
                byteSize = (int)(file.length() - sampleIndex);

                break;
            }

            if(index.get(keys[i+increIndex])>0)
            {
                byteSize = (int)(index.get(keys[i+increIndex]) - sampleIndex);
                break;
            }

            increIndex++;
        }

        //byteSize = (int)(index.get(keys[i+1]) - sampleIndex);
        bytes = new byte[byteSize];

        file.readFully(bytes);

//for(byte b:bytes)
//System.out.print((char)b);

        pos=0;

        ch = (char)bytes[pos];

        //        pos++;
        //Remove Z, S, I, D lines
        //System.out.println("sample");
	boolean isWrongScan=false;
        while( (ch=(char)bytes[pos]) == 'Z' || ch == 'S' || ch == 'I' || ch == 'D')
	{
		StringBuffer tmpsb = new StringBuffer();
		while( ch != CARRIAGE_RETURN )
		{
			pos++;
			ch = (char)bytes[pos];
			tmpsb.append(ch);

		}

		String tmpstr = tmpsb.toString();
		if(tmpstr.startsWith("\tActiv")) {

			if(!tmpstr.trim().endsWith(conf.getScanShiftType())) {
				isWrongScan=true;
				break;
			}

		}

            pos++;
        }

	if(isWrongScan)
		return null;

	    int arrSize=0;
	    for(int j=pos;j<byteSize;j++)
	    {
		if( CARRIAGE_RETURN == (char)bytes[j] )
		    arrSize++;

	    }
	    double[] massArr = new double[arrSize];
	    double[] intArr = new double[arrSize];

	    //pos++;
	    mass = new StringBuilder(10);
	    intensity = new StringBuilder(10);
	    intensity.append('0');

	    boolean isMass=true;
	    boolean isInt=true;

	    int massIndex=0;
	    //double tempMass;
	    int spaceCount=0;
	    for(int j=pos;j<byteSize;j++)
	    {
		ch = (char)bytes[j];
                //System.out.print(ch);

		//fixed on new ms2 files
		switch(ch)
		{
		    case WINDOW_CR:
			break;

		    case SPACE:
			spaceCount++;
			isMass=false;
			if(spaceCount==1) isInt=true;

			break;

		    case CARRIAGE_RETURN:
			isMass=true;
			isInt=false;
			spaceCount=0;

                        //try {
			intArr[massIndex] = Long.parseLong(intensity.toString());
			massArr[massIndex++] = Double.parseDouble(mass.toString());
			mass.delete(0, mass.length());  //This is faster than creating new StringBuilder object
			intensity.delete(0, intensity.length()).append('0');

			break;

		    case DOT:
			isInt=false;

		    default:
			if(isMass)
			    mass.append(ch);
			else if(isInt) //remove decimal value of intensity
			    intensity.append(ch);

                    break;
            }

        }


        int tmpIndex=0;

        double[] correctIntArr = new double[massList.size()];
        double[] intSumArr = new double[massList.size()];

    //    for(int i1=0;i1<massArr.length;i1++)
     //       System.out.println(massArr[i1] +"\t" + intArr[i1]);


        for(Iterator<ReportIon> itr=massList.iterator(); itr.hasNext(); )
        {
		    ReportIon ri = itr.next();
            double d = ri.getMass();

            //System.out.println("===========" + ri.getcPlusList());
            //System.out.println("===========" + ri.getcMinusList());

            intSumArr[tmpIndex] = intensitySumForSinglePeak(massArr, intArr, d, tolerance);

            //for(double each:correctIntArr)
            //    System.out.print(each +"\t");
            //  System.out.println("");

            if(tmpIndex>0 && ri.getcMinusList().size()>0) {
                double beforeInt = 0;
                beforeInt = intSumArr[tmpIndex] * 0.01 * ri.getcMinusList().get(0);

                correctIntArr[tmpIndex-1] -= beforeInt;
              //  System.out.println("b\t" + beforeInt + " " + sum + " " + ri.getcMinusList().get(0));

            }
            if(tmpIndex<correctIntArr.length-1 && ri.getcPlusList().size()>0) {
                double afterInt = 0;
                afterInt = intSumArr[tmpIndex] * 0.01 * ri.getcPlusList().get(0);

                correctIntArr[tmpIndex+1] -= afterInt;
                //System.out.println("a\t" + afterInt + " " + sum + " " + ri.getcPlusList().get(0));
            }

            //for(double each:correctIntArr)
            //    System.out.print(each +"\t");
            //  System.out.println("");

            //result.append((long)sum).append(" ");

            tmpIndex++;
        }

        double[] correctedIntSumArr = new double[correctIntArr.length];

        for(int ii=0;ii<correctIntArr.length;ii++) {
            correctedIntSumArr[ii] = intSumArr[ii] + correctIntArr[ii];

            if(correctedIntSumArr[ii]<=0)
                result.append("0 ");
            else
                result.append((long)correctedIntSumArr[ii]).append(" ");


         //   System.out.println("=================" + (correctedIntSumArr[ii]<=0) + "\t" + (long)correctedIntSumArr[ii] +"\t"+ correctedIntSumArr[ii] + "\t" + intSumArr[ii] +"\t"+ correctIntArr[ii]);
         //   result.append((long)correctedIntSumArr[ii]).append(" ");
        }

        result.append(",");

        for(double each:intSumArr) {
            result.append((long)each).append(" ");
        }

        return result.toString();

    }
*/

  public static double[][] readSpectrumOnly(
    int i,
    IndexedFile iFile,
    Configuration conf
  ) throws IOException {


    RandomAccessFile file = iFile.getFile();
    TIntLongHashMap index = iFile.getMsIndex();

    int[] keys = iFile.getKeys();

    int byteSize;
    byte[] bytes;
    char ch;
    int pos;
    if (i < 0) {
      return null;
    }

    long sampleIndex = index.get(keys[i]);

    StringBuilder mass, intensity;

    if (sampleIndex < 0) {
      return null;
    }

    if (conf.getScanShift() != 0) {

      int tmpIndex = i + conf.getScanShift();

      if (tmpIndex < 0 || tmpIndex >= keys.length) return null;

//                String scanShiftPrecursor = getPrecursorMass( index, i+1, index.get(keys[i]), file, keys );
//                String currentPrecursor = getPrecursorMass( index, i+1-conf.getScanShift(), index.get(keys[i-conf.getScanShift()]), file, keys );
//compare thise
    }


    file.seek(sampleIndex);

    int increIndex = 1;

    while (true) {
      if ((i + increIndex) >= keys.length) {
        byteSize = (int) (file.length() - sampleIndex);

        break;
      }

      if (index.get(keys[i + increIndex]) > 0) {
        byteSize = (int) (index.get(keys[i + increIndex]) - sampleIndex);
        break;
      }

      increIndex++;
    }

    //byteSize = (int)(index.get(keys[i+1]) - sampleIndex);
    bytes = new byte[byteSize];

    file.readFully(bytes);


    pos = 0;

    ch = (char) bytes[pos];

    //        pos++;
    //Remove Z, S, I, D lines
    //System.out.println("sample");
    boolean isWrongScan = false;
    while ((ch = (char) bytes[pos]) == 'Z' || ch == 'S' || ch == 'I' || ch == 'D') {
      StringBuffer tmpsb = new StringBuffer();
      while (ch != CARRIAGE_RETURN) {
        pos++;
        ch = (char) bytes[pos];
        tmpsb.append(ch);
      }

      String tmpstr = tmpsb.toString();
      if (tmpstr.startsWith("\tActiv")) {
        if (!tmpstr.trim().endsWith(conf.getScanShiftType()) && !conf.isMs2ScanTypeValueUnassigned()) {
          isWrongScan = true;
          break;
        }
      }

      pos++;
    }

    if (isWrongScan)
      return null;

        /*
      int arrSize=0;
	    for(int j=pos;j<byteSize;j++)
	    {
		if( CARRIAGE_RETURN == (char)bytes[j] )
		    arrSize++;

	    }*/


    double[][] result = parseSpectra(bytes);

    return result;

  }


  public static double[][] parseSpectra(byte[] bytes) {


    String str = new String(bytes);
    String[] lines = str.split(System.getProperty("line.separator"));


    //    for(String each:lines)
    //      System.out.println(each);


    int count = 0;

    for (String each : lines) {

      if (Character.isDigit(each.charAt(0)))
        break;

      count++;
    }

    double[] massArr = new double[lines.length - count];
    double[] intArr = new double[lines.length - count];

    int massIndex = 0;
    for (int i = count; i < lines.length; i++) {

      String[] arr = lines[i].split(" ");
      //       System.out.println(lines[i]);

      massArr[massIndex] = Double.parseDouble(arr[0]);
      intArr[massIndex] = Double.parseDouble(arr[1]);
      massIndex++;
    }

    double[][] result = new double[2][2];
    result[0] = massArr;
    result[1] = intArr;

    return result;
  }

  public static double[] readSpectrumSpecificPeaks(
    int i,
    IndexedFile iFile,
    List massList,
    Configuration conf
  ) throws IOException {


    RandomAccessFile file = iFile.getFile();
    TIntLongHashMap index = iFile.getMsIndex();

    int[] keys = iFile.getKeys();
    double[] result = new double[massList.size() + 1];

    int byteSize;
    byte[] bytes;
    char ch;
    int pos;

    if (i < 0) {
      return result;
    }

    long sampleIndex = index.get(keys[i]);

    StringBuilder mass, intensity;

    if (sampleIndex < 0) {
      return result;
    }

    if (conf.getScanShift() != 0) {

      int tmpIndex = i + conf.getScanShift();

      if (tmpIndex < 0 || tmpIndex >= keys.length) return result;

      //String scanShiftPrecursor = getPrecursorMass( index, tmpIndex, index.get(keys[i+conf.getScanShift()]), file, keys );
      //String currentPrecursor = getPrecursorMass( index, i, index.get(keys[i]), file, keys );
      String scanShiftPrecursor = getPrecursorMass(index, i + 1, index.get(keys[i]), file, keys);
      String currentPrecursor = getPrecursorMass(index, i + 1 - conf.getScanShift(), index.get(keys[i - conf.getScanShift()]), file, keys);
//compare thise
    }

    result[0] = keys[i];

    file.seek(sampleIndex);

    int increIndex = 1;

    while (true) {
      if ((i + increIndex) >= keys.length) {
        byteSize = (int) (file.length() - sampleIndex);

        break;
      }

      if (index.get(keys[i + increIndex]) > 0) {
        byteSize = (int) (index.get(keys[i + increIndex]) - sampleIndex);
        break;
      }

      increIndex++;
    }

    //byteSize = (int)(index.get(keys[i+1]) - sampleIndex);
    bytes = new byte[byteSize];

    file.readFully(bytes);
    pos = 0;

    ch = (char) bytes[pos];

    //        pos++;
    //Remove Z, S, I, D lines
    //System.out.println("sample");
    boolean isWrongScan = false;
    while ((ch = (char) bytes[pos]) == 'Z' || ch == 'S' || ch == 'I' || ch == 'D') {
      StringBuffer tmpsb = new StringBuffer();
      while (ch != CARRIAGE_RETURN) {
        pos++;
        ch = (char) bytes[pos];
        tmpsb.append(ch);

      }

      String tmpstr = tmpsb.toString();
      if (tmpstr.startsWith("\tActiv")) {
        if (!tmpstr.trim().endsWith(conf.getScanShiftType())) {
          isWrongScan = true;
          break;
        }

      }

      pos++;
    }

    if (isWrongScan)
      return null;


    int tmpIndex = 1;


    double[][] tmpArr = parseSpectra(bytes);
    double[] massArr = tmpArr[0];
    double[] intArr = tmpArr[1];


    for (Iterator<ReportIon> itr = massList.iterator(); itr.hasNext(); ) {
      ReportIon ri = itr.next();
      double d = ri.getMass();
      double tolerance = conf.getMsmsSpecificTolerance();
      if (conf.isMsmsTolerancePpm())
        tolerance = PmsUtil.getMZFromPpm((float) d, (float) tolerance);

      double sum = intensitySumForSinglePeak(massArr, intArr, d, tolerance);

      result[tmpIndex++] = sum;
    }


    return result;

  }


  public static double intensitySumForSinglePeakNewMassTolerance(double[] massArr, double[] intArr, double mass, double tolerance) {
    double sumIntensity = 0;

    double mzTolerance = PmsUtil.getMZFromPpm(mass, tolerance);

    double startMass = mass - mzTolerance;
    double endMass = mass + mzTolerance;

    //System.out.println(startMass + " " + endMass);
    int start = Arrays.binarySearch(massArr, startMass);
    int end = Arrays.binarySearch(massArr, endMass);

    //System.out.println(startMass + " " + endMass);
    if (start < 0)
      start = -start - 1;

    if (end < 0)
      end = -end - 2;

    //double inttemp=0;
//	double masstemp=0;


    for (int i = start; i <= end; i++) {


//           System.out.println("data\t" + mass + "\t" + massArr[i] + "\t" + intArr[i]);
      //  if(inttemp< intArr[i]) {
//		masstemp=massArr[i];
      //	inttemp = intArr[i];
      //  }


      sumIntensity += intArr[i];
    }

//	if(masstemp!=0)
    //     System.out.println("mass============\t" + mass + "\t" + sumIntensity); //

    return sumIntensity;
  }
  private static int count =0;
  public static double intensitySumForSinglePeak(double[] massArr, double[] intArr, double mass, double tolerance) {
    double sumIntensity = 0;

    double startMass = mass - tolerance;
    double endMass = mass + tolerance;

    //System.out.println(startMass + " " + endMass);
    int start = Arrays.binarySearch(massArr, startMass);
    int end = Arrays.binarySearch(massArr, endMass);

    start = ~start;
    end = ~end;
    end--;

    /*
    //System.out.println(startMass + " " + endMass);
    if (start < 0)
      start = -start - 1;

    if (end < 0)
      end = -end - 2;
    if(start>end)
    {
      //count++;
     // System.out.println("DEBUG>>>>>>"+count);
    }
    */

    //double inttemp=0;
//	double masstemp=0;
    for (int i = start; i <= end; i++) {
//           System.out.println("data\t" + mass + "\t" + massArr[i] + "\t" + intArr[i]);
      //  if(inttemp< intArr[i]) {
//		masstemp=massArr[i];
      //	inttemp = intArr[i];
      //  }


      sumIntensity += intArr[i];
    }

//	if(masstemp!=0)
    //     System.out.println("mass============\t" + mass + "\t" + sumIntensity); //

    return sumIntensity;
  }


  //read frag ions only for light one.  For MRM CRV
  public static double[] readSimpleSpectrum(
    MRMPeptideModel model,
    //.getBionArr(), model.getYionArr()
    ListIterator<edu.scripps.pms.util.spectrum.Peak> peakList
  ) throws IOException {
    double[][] bionArr = model.getBionSample();
    double[][] yionArr = model.getYionSample();

    TDoubleArrayList massList = new TDoubleArrayList();
    TDoubleArrayList intList = new TDoubleArrayList();

    for (Iterator<edu.scripps.pms.util.spectrum.Peak> itr = peakList; itr.hasNext(); ) {
      edu.scripps.pms.util.spectrum.Peak p = itr.next();

      massList.add(p.getM2z());
      intList.add(p.getIntensity());
    }

    double[] massArr = massList.toNativeArray();
    double[] intArr = intList.toNativeArray();


    int bionSum = 0;
    int yionSum = 0;

    //boolean tempBool=false;
    double mass = model.getParentMass();

    double[] ionArr = new double[bionArr.length];


    for (int k = 0; k < bionArr.length; k++) {
      if (bionArr[k][1] > mass) {
        double bion = intensitySum(massArr, intArr, bionArr[k][1], bionArr[k][2]);
        bionSum += bion;

        //    System.out.print(bionArr[k][1] + "\t" + bion  + "\t");
      }

      if (yionArr[k][1] > mass) {
        double yion = intensitySum(massArr, intArr, yionArr[k][1], yionArr[k][2]);
        yionSum += yion;

        //      System.out.print(yionArr[k][1]  + "\t" + yion  + "\t");
      }

//            sampleSum += result[ION_START_INDEX+k];
      //          sampleSum += result[ION_START_INDEX+bioSample.length+k];

    }

    return null;
  }

  //read frag ions only for light one.  For MRM CRV
  public static double[] readSpectrum(
    int[] keys,
    int i,
    IndexedFile iFile,
    //TIntLongHashMap index,
    //int diff,
    double[][] bioSample,
    double[][] yioSample
  ) throws IOException {
    RandomAccessFile file = iFile.getFile();
    TIntLongHashMap index = iFile.getMsIndex();

    //int[] result = new int[3];
    ION_START_INDEX = 2; //no labeling
    int ionLength = bioSample.length;
    double[] result = new double[ionLength * 2 + ION_START_INDEX];

    //if( (i+diff)>=keys.length )
    //    return result;

    int byteSize;
    byte[] bytes;
    char ch;
    int pos;

    if (i < 0) {
      return result;
    }

    long sampleIndex = index.get(keys[i]);

    StringBuilder mass, intensity;

    if (sampleIndex < 0) {
      return result;
    }


    result[0] = keys[i];

    /**********************************
     * Go to Spectrum for sample, and calculate intensity
     **********************************/
    file.seek(sampleIndex);

    int increIndex = 1;

    while (true) {
      if ((i + increIndex) >= keys.length) {
        byteSize = (int) (file.length() - sampleIndex);

        break;
      }

      if (index.get(keys[i + increIndex]) > 0) {
        byteSize = (int) (index.get(keys[i + increIndex]) - sampleIndex);
        break;
      }

      increIndex++;
    }

    //byteSize = (int)(index.get(keys[i+1]) - sampleIndex);
    bytes = new byte[byteSize];

    file.readFully(bytes);
    pos = 0;

    ch = (char) bytes[pos];

    //        pos++;
    //Remove Z, S, I, D lines
    //System.out.println("sample");
    while ((ch = (char) bytes[pos]) == 'Z' || ch == 'S' || ch == 'I' || ch == 'D') {
      while (ch != CARRIAGE_RETURN) {
        pos++;
        ch = (char) bytes[pos];
      }

      pos++;
    }

    double[][] specArr = parseSpectra(bytes);


    int sampleSum = 0;


    double[][] specResult = parseSpectra(bytes);
    double[] massArr = specResult[0];
    double[] intArr = specResult[1];


    //boolean tempBool=false;
    for (int k = 0; k < bioSample.length; k++) {
      //result[ION_START_INDEX+k] = intensitySum(massArr, intArr, bioSample[k][startMassIndex], bioSample[k][endMassIndex]);
      //result[ION_START_INDEX+ionLength+k] = intensitySum(massArr, intArr, yioSample[k][startMassIndex], yioSample[k][endMassIndex]);


      result[ION_START_INDEX + k] = intensitySum(massArr, intArr, bioSample[k][1], bioSample[k][2]);
      result[ION_START_INDEX + ionLength + k] = intensitySum(massArr, intArr, yioSample[k][1], yioSample[k][2]);

      sampleSum += result[ION_START_INDEX + k];
      sampleSum += result[ION_START_INDEX + bioSample.length + k];


      //sampleSum += intensitySum(massArr, intArr, bioSample[k][1], bioSample[k][2]);
      //sampleSum += intensitySum(massArr, intArr, yioSample[k][1], yioSample[k][2]);
    }
    result[1] = sampleSum;

    return result;
  }


  public static double[] readSpectrumByRt(
    double[] keys,
    int i,
    IndexedFile iFile,
//            Object ofile,
    MRMPeptideModel model
  ) throws IOException {
    double[][] bionArr = model.getBionArr();
    double[][] yionArr = model.getYionArr();

    //System.out.println(bionArr.length);

    int ionLength = bionArr.length;
    double[] result = new double[ionLength * 4 + ION_START_INDEX];

    RandomAccessFile file = iFile.getFile();

    int byteSize;
    byte[] bytes;
    char ch;
    int pos;

    if (i < 0 || i >= keys.length) {
      return result;
    }

    long index = iFile.getPosByRt(keys[i]);

    StringBuilder mass, intensity;

    if (index < 0) {
      return result;
    }

    result[0] = keys[i];

    /**********************************
     * Go to Spectrum for sample, and calculate intensity
     *********************************/
    file.seek(index);

    int increIndex = 1;

    while (true) {
      if ((i + increIndex) >= keys.length) {
        byteSize = (int) (file.length() - index);

        break;
      }

      if (iFile.getPosByRt(keys[i + increIndex]) > 0) {
        byteSize = (int) (iFile.getPosByRt(keys[i + increIndex]) - index);
        break;
      }

      increIndex++;
    }

    bytes = new byte[byteSize];

    file.readFully(bytes);
    pos = 0;

    ch = (char) bytes[pos];


    while ((ch = (char) bytes[pos]) == 'Z' || ch == 'S' || ch == 'I' || ch == 'D') {
      while (ch != CARRIAGE_RETURN) {
        pos++;
        ch = (char) bytes[pos];
      }

      pos++;
    }


    int arrSize = 0;
    for (int j = pos; j < byteSize; j++) {
      if (CARRIAGE_RETURN == (char) bytes[j])
        arrSize++;

    }


    double[][] specResult = parseSpectra(bytes);
    double[] massArr = specResult[0];
    double[] intArr = specResult[1];

    int intSum = 0;

    //boolean tempBool=false;
    for (int k = 0; k < bionArr.length; k++) {

      result[ION_START_INDEX + k] = intensitySum(massArr, intArr, bionArr[k][1], bionArr[k][2]);
      result[ION_START_INDEX + ionLength + k] = intensitySum(massArr, intArr, yionArr[k][1], yionArr[k][2]);

      intSum += result[ION_START_INDEX + k];
      intSum += result[ION_START_INDEX + bionArr.length + k];

    }

    result[1] = intSum;


//        for(int ii=0;ii<result.length;ii++)


    return result;

  }

  //mzxml
  public static String calculateFullMS(
    int scanNum,
    double samStartMass, double samEndMass,
    double refStartMass, double refEndMass,
    Configuration conf,
    SpecRange range)
    throws IOException, CensusIndexOutOfBoundException {

        /*
        TIntLongHashMap index = iFile.getMsIndex();
        RandomAccessFile file = iFile.getFile();
        int[] keys = iFile.getKeys();

        SpectrumModel sModel = new SpectrumModel();
        sModel.setHighRes(false);
        sModel.setKeys(keys);
        sModel.setIndex(index);
        sModel.setIFile(iFile);
        sModel.setFile(file);
        sModel.setSamStartMass(samStartMass);
        sModel.setSamEndMass(samEndMass);
        sModel.setRefStartMass(refStartMass);
        sModel.setRefEndMass(refEndMass);

        return peakFinding(sModel, range, conf, keyIndex);
*/
    return null;
    //return sModel.getResult();
  }

  //low resolution
  public static String calculateFullMS(
    int keyIndex,
    IndexedFile iFile,
    double samStartMass, double samEndMass,
    double refStartMass, double refEndMass,
    SpecRange range)
    throws IOException, CensusIndexOutOfBoundException, Exception {
    conf = Configuration.getInstance();
    TIntLongHashMap index = iFile.getMsIndex();
    RandomAccessFile file = iFile.getFile();
    MzxmlSpectrumReader mzReader = iFile.getMzreader();
    int[] keys = iFile.getKeys();

    SpectrumModel sModel = new SpectrumModel();
    sModel.setHighRes(false);
    sModel.setKeys(keys);
    sModel.setIndex(index);
    sModel.setIFile(iFile);
    sModel.setFile(file);
    sModel.setMzxmlreader(mzReader);
    sModel.setSamStartMass(samStartMass);
    sModel.setSamEndMass(samEndMass);
    sModel.setRefStartMass(refStartMass);
    sModel.setRefEndMass(refEndMass);

    return peakFinding(sModel, range, keyIndex, null);

    //return sModel.getResult();
  }

  //high resolution
  public static class ResultList extends ArrayList<double[][]> {
    public ResultList(int size, int maxWindow, int margin) {
      super();

      for (int i = 0; i < size; i++) {
        add(new double[maxWindow * 2 + 1 + margin * 2][2]);
      }
    }

    public void addElement(int chromIndex, int scanIndex, double scanNumber, double intensity) {
      double[][] arr = get(chromIndex);
      arr[0][scanIndex] = scanNumber;
      arr[1][scanIndex] = intensity;
    }
  }

  //for non-labeling peak finding labelfree
  public static String peakFinding(
    SpectrumModel sModel, SpecRange range, Configuration conf, int keyIndex, double[] samIsoArr)
    throws IOException, Exception {
    steepRatioThreshold = conf.getSteepRatioThreshold();

    int maxWindow = conf.getMaxWindow();
    int margin = conf.getMargin();

    NonLabelMappingModel mapModel = conf.getMapModel();
    //Hashtable<String, Hashtable> ms2ms1Ht = mapModel.getMs2ms1Ht();

    // System.out.println(mapModel.getMs2ms1Ht());
    // System.exit(0);

    Vector<String> pathFileNameList = mapModel.getPathFileNameList();

    //ResultList result = new ResultList(pathArray.length+1, maxWindow, margin);
    int resultSize = pathFileNameList.size() * 2;
    int resultIntensityStart = pathFileNameList.size();

    double[][] result = null;
    if (conf.getQuantLevel() == 1)
      result = new double[maxWindow * 2 + 1 + margin * 2][resultSize];
    else if (conf.getQuantLevel() == 2)
      result = new double[maxWindow * 2 + 1 + margin * 2][pathFileNameList.size() * 2 * sModel.getBioSample().length + resultSize]; //scan #, sample intensity, ref intensity

    int maxScanIndex = mapModel.getMaxScanIndex(pathFileName);

    if (maxScanIndex < 0) {
      System.out.println("Error : Spectral file name in the config file is invalid : " + pathFileName);
      throw new Exception();
    }


    //int maxScanIndex = result.length;

    int leftIndex = maxWindow + margin;
    int rightIndex = maxWindow + margin + 1;

    double totalIntensity = 0;

    int steepArea = conf.getSteepArea();

    int moveLeftKeyIndex = keyIndex;
    int moveRightKeyIndex = -1;
    //int[] keys = sModel.getKeys();


    if (conf.getQuantLevel() == 1)
      moveRightKeyIndex = keyIndex + 1;
    else if (conf.getQuantLevel() == 2)
      moveRightKeyIndex = findLightKey(sModel.getIFile().getKeys(), keyIndex, MOVE_RIGHT, sModel.getIFile(), conf);

    int initWin = 2;
    Hashtable<Integer, double[]> resultHt = new Hashtable<Integer, double[]>();

    for (int i = 0; i < initWin; i++) {
      if (moveLeftKeyIndex <= 0 || leftIndex <= 0) {
        if (conf.getQuantLevel() == 1)
          moveLeftKeyIndex++;
        else
          moveLeftKeyIndex = findLightKey(sModel.getIFile().getKeys(), moveLeftKeyIndex, MOVE_RIGHT, sModel.getIFile(), conf);

        leftIndex++;

        break;
      }


      if (conf.isHighRes())
        result[leftIndex] = readFullSpectrum(moveLeftKeyIndex, null, null, samIsoArr, mapModel, null, resultHt);//
      else
        result[leftIndex] = readFullSpectrum(moveLeftKeyIndex, null, null, null, mapModel, sModel, resultHt);

      // System.exit(0);

      //result = new double[maxWindow*2+1+margin*2][4*sModel.getBioSample().length+resultSize]; //scan #, sample intensity, ref intensity

      for (int j = resultIntensityStart; j < result[leftIndex].length; j++)
        totalIntensity += result[leftIndex][j];

      if (conf.getQuantLevel() == 1)
        moveLeftKeyIndex--;
      else
        moveLeftKeyIndex = findLightKey(sModel.getIFile().getKeys(), moveLeftKeyIndex, MOVE_LEFT, sModel.getIFile(), conf);

      //if it is null, continue without changing leftIndex value;
      if (0 == result[leftIndex][0]) {
        i--;
        continue;
      }

      leftIndex--;
    }

    for (int i = 0; i < initWin; i++) //leftIndex=maxWindow;keys[leftPeakIndex]>=leftMinValue; )
    {
      if (moveRightKeyIndex > maxScanIndex || moveRightKeyIndex < 0) {
        if (conf.getQuantLevel() == 1)
          moveRightKeyIndex--;
        else
          moveRightKeyIndex = findLightKey(sModel.getIFile().getKeys(), moveRightKeyIndex, MOVE_LEFT, sModel.getIFile(), conf);

        rightIndex--;

        break;
      }

      if (conf.isHighRes())
        result[rightIndex] = readFullSpectrum(moveRightKeyIndex, null, null, samIsoArr, mapModel, null, resultHt);
      else
        result[rightIndex] = readFullSpectrum(moveRightKeyIndex, null, null, null, mapModel, sModel, resultHt);

      for (int j = resultIntensityStart; j < result[leftIndex].length; j++)
        totalIntensity += result[leftIndex][j];

      if (conf.getQuantLevel() == 1)
        moveRightKeyIndex++;
      else
        moveRightKeyIndex = findLightKey(sModel.getIFile().getKeys(), moveRightKeyIndex, MOVE_RIGHT, sModel.getIFile(), conf);

      //if it is null, continue without changing leftIndex value;
      if (0 == result[rightIndex][0]) {
        i--;
        continue;
      }

      rightIndex++;
    }

    if (moveRightKeyIndex < 0 && moveLeftKeyIndex < 0)
      return null;

    boolean isGoingUp = true;
    boolean isHighIntensity = true;  //if the intensity is lower than one third of average intensity, this becomes false
    double[][] arr = new double[steepArea][resultSize];
    double[][] prevArr = new double[steepArea][resultSize];
    Hashtable<Integer, double[]> tmpKeyHt = new Hashtable<Integer, double[]>();

    if (leftIndex + steepArea >= 0) {
      //for(int i=0;i<300;i++)
      while (true) {
        double area1 = 0;
        double area2 = 0;

        if (moveLeftKeyIndex - steepArea < 0)
          break;

        int steepCount = 0;
        int tempKeyIndex = moveLeftKeyIndex;
        while (true) {
          if (tempKeyIndex < 0 || steepCount < 0 || steepCount >= steepArea - 1)
            break;

          if (prevArr[0][0] != 0) {
            for (int l = 0; l < 2; l++)
              for (int m = 0; m < resultSize; m++)
                arr[l][m] = prevArr[l + 1][m];
          } else {

            if (conf.getQuantLevel() == 2) {
              if (conf.isHighRes()) {
                System.out.println("Error : we do not support high resolution analysis for MSMS yet.  Email to rpark@scripps.edu");
                System.exit(0);
              } else {
                arr[steepCount] = readSpectrum(sModel, tempKeyIndex);
              }
            } else if (conf.getQuantLevel() == 1) {
              double[] tmpKeyArr = tmpKeyHt.get(tempKeyIndex);
              if (null != tmpKeyArr)
                arr[steepCount] = tmpKeyArr;
              else {
                if (conf.isHighRes())
                  arr[steepCount] = readFullSpectrum(tempKeyIndex, null, null, samIsoArr, mapModel, null, resultHt);
                else
                  arr[steepCount] = readFullSpectrum(tempKeyIndex, null, null, null, mapModel, sModel, resultHt);

                tmpKeyHt.put(tempKeyIndex, arr[steepCount]);
              }
            }
          }

          steepCount++;
          if (conf.getQuantLevel() == 2)
            tempKeyIndex = findLightKey(sModel.getIFile().getKeys(), tempKeyIndex, MOVE_LEFT, sModel.getIFile(), conf);
          else
            tempKeyIndex--;
        }

        if (tempKeyIndex < 0)
          break;

        if (conf.getQuantLevel() == 2) {
          if (conf.isHighRes()) {
            System.out.println("Error : we do not support high resolution analysis for MSMS yet.  Email to rpark@scripps.edu");
            System.exit(0);
          } else {
            arr[steepCount] = readSpectrum(sModel, tempKeyIndex);
          }
        } else {
          double[] tmpKeyArr = tmpKeyHt.get(tempKeyIndex);
          if (null != tmpKeyArr)
            arr[steepCount] = tmpKeyArr;
          else {
            if (conf.isHighRes())
              arr[steepCount] = readFullSpectrum(tempKeyIndex, null, null, samIsoArr, mapModel, null, resultHt);
            else
              arr[steepCount] = readFullSpectrum(tempKeyIndex, null, null, null, mapModel, sModel, resultHt);

            tmpKeyHt.put(tempKeyIndex, arr[steepCount]);
          }
        }

        for (int k = 0; k < 3; k++) {
          for (int l = 0; l < resultSize; l++)
            prevArr[k][l] = arr[k][l];
        }

        if (conf.getQuantLevel() == 2) {

          if (conf.isHighRes()) {
            System.out.println("Error : we do not support high resolution analysis for MSMS yet.  Email to rpark@scripps.edu");
            System.exit(0);
          } else {
            for (int j = leftIndex + 1; j <= leftIndex + steepCount + 1; j++) {
              for (int k = 1; k <= resultIntensityStart; k++) {
                area1 += result[j][k];
              }
            }

            for (int j = 0; j <= steepCount; j++) {
              for (int k = 1; k <= resultIntensityStart; k++) {
                area2 += arr[j][k];
              }
            }
          }
        } else {

          for (int j = leftIndex + 1; j <= leftIndex + steepCount + 1; j++) {
            for (int k = resultIntensityStart; k < result[j].length; k++) {
              area1 += result[j][k];
            }
          }

          for (int j = 0; j <= steepCount; j++) {
            for (int k = resultIntensityStart; k < result[j].length; k++) {
              area2 += arr[j][k];
            }
          }

        }

        if (area2 == 0 && area1 == 0)
          break;

        double areaRatio = (double) area2 / area1;

        /*****************************************
         * AREA1 is right side area
         * AREA2 is left side area
         *****************************************/
        if (isHighIntensity && ((double) totalIntensity / (rightIndex - leftIndex - 1) / 1.5 > area2 / 3))
          isHighIntensity = false;

        if ((isGoingUp && isHighIntensity) || areaRatio < steepRatioThreshold) {
          if (leftIndex <= 0) {
            isGoingUp = false;
            break;
          }

          result[leftIndex] = arr[0];

          if (conf.getQuantLevel() == 2)
            moveLeftKeyIndex = findLightKey(sModel.getIFile().getKeys(), moveLeftKeyIndex, MOVE_LEFT, sModel.getIFile(), conf);
          else
            moveLeftKeyIndex--; // -= numIsoWindow;

          if (0 == result[leftIndex][0])
            continue;

          leftIndex--;

          for (int j = 0; j < steepCount; j++)
            arr[j] = arr[j + 1];

          for (int j = resultIntensityStart; j < arr.length; j++)
            totalIntensity += arr[0][j];

          continue;
        } else
          isGoingUp = false;


        if (areaRatio > steepRatioThreshold) {
          for (int k = 0; k < 3; k++) {
            if (result[leftIndex + 1][resultIntensityStart] < arr[k][resultIntensityStart] || arr[k][0] == 0)
              break;

            if (leftIndex < 0)
              break;

            result[leftIndex] = arr[k];

            if (conf.getQuantLevel() == 2)
              moveLeftKeyIndex = findLightKey(sModel.getIFile().getKeys(), moveLeftKeyIndex, MOVE_LEFT, sModel.getIFile(), conf);
            else
              moveLeftKeyIndex--; // -= numIsoWindow;

            leftIndex--;
          }

          break;
        }

        result[leftIndex--] = arr[0];
        if (conf.getQuantLevel() == 2)
          moveLeftKeyIndex = findLightKey(sModel.getIFile().getKeys(), moveLeftKeyIndex, MOVE_LEFT, sModel.getIFile(), conf);
        else
          moveLeftKeyIndex--; // -= numIsoWindow;

        for (int j = 0; j < steepCount; j++)
          arr[j] = arr[j + 1];
      }
    }

    prevArr = new double[steepArea][resultSize]; //clean up the array

    isGoingUp = true;
    isHighIntensity = true;

    if (rightIndex + steepArea < result.length) {
      while (true) {
        double area1 = 0;
        double area2 = 0;

        if (moveRightKeyIndex + steepArea >= maxScanIndex)
          break;

        int steepCount = 0;

        int tempKeyIndex = moveRightKeyIndex;
        while (true) {
          if (tempKeyIndex >= maxScanIndex || steepCount >= steepArea - 1 || steepCount >= maxScanIndex)
            break;

          if (prevArr[0][0] != 0) {
            for (int l = 0; l < 2; l++)
              for (int m = resultIntensityStart; m < resultSize; m++)
                arr[l][m] = prevArr[l + 1][m];
          } else {
            if (conf.getQuantLevel() == 2) {
              if (conf.isHighRes()) {
                System.out.println("Error : we do not support high resolution analysis for MSMS yet.  Email to rpark@scripps.edu");
                System.exit(0);
              } else {
                arr[steepCount] = readSpectrum(sModel, tempKeyIndex);
              }
            } else if (conf.getQuantLevel() == 1) {
              double[] tmpKeyArr = tmpKeyHt.get(tempKeyIndex);
              if (null != tmpKeyArr)
                arr[steepCount] = tmpKeyArr;
              else {
                if (conf.isHighRes())
                  arr[steepCount] = readFullSpectrum(tempKeyIndex, null, null, samIsoArr, mapModel, null, resultHt);
                else
                  arr[steepCount] = readFullSpectrum(tempKeyIndex, null, null, null, mapModel, sModel, resultHt);

                tmpKeyHt.put(tempKeyIndex, arr[steepCount]);
              }
            }
          }

          steepCount++;

          if (conf.getQuantLevel() == 2)
            tempKeyIndex = findLightKey(sModel.getIFile().getKeys(), tempKeyIndex, MOVE_RIGHT, sModel.getIFile(), conf);
          else
            tempKeyIndex++; // = tempKeyIndex + 1*numIsoWindow;
        }

        if (tempKeyIndex >= maxScanIndex)
          break;

        if (conf.getQuantLevel() == 2) {
          if (conf.isHighRes()) {
            System.out.println("Error : we do not support high resolution analysis for MSMS yet.  Email to rpark@scripps.edu");
            System.exit(0);
          } else {
            arr[steepCount] = readSpectrum(sModel, tempKeyIndex);
          }
        } else {
          double[] tmpKeyArr = tmpKeyHt.get(tempKeyIndex);
          if (null != tmpKeyArr)
            arr[steepCount] = tmpKeyArr;
          else {
            if (conf.isHighRes())
              arr[steepCount] = readFullSpectrum(tempKeyIndex, null, null, samIsoArr, mapModel, null, resultHt);
            else
              arr[steepCount] = readFullSpectrum(tempKeyIndex, null, null, null, mapModel, sModel, resultHt);

            tmpKeyHt.put(tempKeyIndex, arr[steepCount]);
          }
        }

        if (null == arr[steepCount])
          break;

        for (int k = 0; k < 3; k++) {
          for (int l = resultIntensityStart; l < resultSize; l++)
            prevArr[k][l] = arr[k][l];
        }


        if (conf.getQuantLevel() == 2) {
          if (conf.isHighRes()) {
            System.out.println("Error : we do not support high resolution analysis for MSMS yet.  Email to rpark@scripps.edu");
            System.exit(0);
          } else {
            for (int j = rightIndex - 1; j > rightIndex - steepCount - 2; j--) {
              for (int k = 1; k <= resultIntensityStart; k++) {
                area1 += result[j][k];
              }
            }

            for (int j = 0; j <= steepCount; j++) {
              //for(int k=resultIntensityStart;k<result[j].length;k++)
              for (int k = 1; k <= resultIntensityStart; k++) {
                area2 += arr[j][k];
              }
            }

          }
        } else {
          for (int j = rightIndex - 1; j > rightIndex - steepCount - 2; j--) {
            for (int k = resultIntensityStart; k < result[j].length; k++) {
              area1 += result[j][k];
            }
          }

          for (int j = 0; j <= steepCount; j++) {
            for (int k = resultIntensityStart; k < result[j].length; k++) {
              area2 += arr[j][k];
            }
          }

        }

        if (area2 == 0 && area1 == 0)
          break;

        double areaRatio = (double) area2 / area1;

        /*****************************************
         * AREA1 is left side area
         * AREA2 is right side area
         *****************************************/

        if (isHighIntensity && ((double) totalIntensity / (rightIndex - leftIndex - 1) / 1.5 > area2 / 3))
          isHighIntensity = false;

        //if( ((isGoingUp && (double)totalIntensity/(rightIndex-leftIndex-1)/3<area2/3)) || areaRatio<steepRatioThreshold )
        if ((isGoingUp && isHighIntensity) || areaRatio < steepRatioThreshold) {
          if ((rightIndex + 1) >= result.length) {
            isGoingUp = false;
            break;
          }

          result[rightIndex] = arr[0];
          if (conf.getQuantLevel() == 2)
            moveRightKeyIndex = findLightKey(sModel.getIFile().getKeys(), moveRightKeyIndex, MOVE_RIGHT, sModel.getIFile(), conf);
          else
            moveRightKeyIndex++;

          if (0 == result[rightIndex][0])
            continue;

          rightIndex++;

          for (int j = 0; j < steepCount; j++)
            arr[j] = arr[j + 1];

          //rightTotalIntensity += arr[0][1];


          for (int j = resultIntensityStart; j < resultSize; j++) {
            totalIntensity += arr[0][j];
          }

          //totalIntensity += arr[0][1];
          //totalIntensity += arr[0][2];

          continue;
        } else
          isGoingUp = false;

        if (area2 / area1 > steepRatioThreshold) {
          for (int k = 0; k < 3; k++) {
            if (result[rightIndex - 1][resultIntensityStart] < arr[k][resultIntensityStart])
              break;

            if (result.length <= rightIndex)
              break;

            result[rightIndex] = arr[k];


            if (conf.getQuantLevel() == 2)
              moveRightKeyIndex = findLightKey(sModel.getIFile().getKeys(), moveRightKeyIndex, MOVE_RIGHT, sModel.getIFile(), conf);
            else
              moveRightKeyIndex++;

            rightIndex++;
          }
          break;
        }

        result[rightIndex++] = arr[0];

        if (conf.getQuantLevel() == 2)
          moveRightKeyIndex = findLightKey(sModel.getIFile().getKeys(), moveRightKeyIndex, MOVE_RIGHT, sModel.getIFile(), conf);
        else
          moveRightKeyIndex++;

        for (int j = 0; j < steepCount; j++)
          arr[j] = arr[j + 1];
      }
    }

/*
        System.out.println("===============");
        for(int iii=0;iii<result.length;iii++)
        {
            //if(result[iii][0]==0)
              //  continue;

            for(int jjj=0;jjj<result[iii].length;jjj++)
                System.out.print(result[iii][jjj] + " ");

            System.out.println("");
        }
        System.out.println("===============");

         System.exit(0);
  */

    int peakStart = leftIndex + 1;
    int peakEnd = rightIndex - 1;

    //NonLabelMappingModel conf.getMapModel();
    int dtaStart = range.getMin();
    int dtaEnd = range.getMax();

    //for(int i=0;i<result.length;i++)
    //   for(int j=0;j<result[i].length;j++)

    int specSpace = 1; //(int)(result[peakStart+1][0]-result[peakStart][0]);

    int diff = (int) (result[peakStart][0] - dtaStart) / specSpace;

    int leftMargin = margin;

    int rightMargin = margin;
    if (diff > 0)
      leftMargin += diff;

    diff = (int) (dtaEnd - result[rightIndex - 1][0]) / specSpace;
    if (diff > 0)
      rightMargin += diff;


    for (int i = 0; i < leftMargin; i++) {
      if (leftIndex < 0 || moveLeftKeyIndex <= 0)
        break;

      if (conf.getQuantLevel() == 2) {

        if (conf.isHighRes())
          ;
        else
          result[leftIndex] = readSpectrum(sModel, moveLeftKeyIndex); //, null, null);
      } else {
        if (conf.isHighRes())
          result[leftIndex] = readFullSpectrum(moveLeftKeyIndex, null, null, samIsoArr, mapModel, null, resultHt);
        else
          result[leftIndex] = readFullSpectrum(moveLeftKeyIndex, null, null, null, mapModel, sModel, resultHt);

      }


      if (0 == result[leftIndex][0]) {
        if (conf.getQuantLevel() == 2)
          moveLeftKeyIndex = findLightKey(sModel.getIFile().getKeys(), moveLeftKeyIndex, MOVE_LEFT, sModel.getIFile(), conf);
        else
          moveLeftKeyIndex--;

        continue;
      }

      leftIndex--;

      if (conf.getQuantLevel() == 2)
        moveLeftKeyIndex = findLightKey(sModel.getIFile().getKeys(), moveLeftKeyIndex, MOVE_LEFT, sModel.getIFile(), conf);
      else
        moveLeftKeyIndex--;
    }


    for (int i = 0; i < rightMargin; i++) {
//        System.out.println(moveRightKeyIndex);
      if (rightIndex >= result.length || moveRightKeyIndex >= maxScanIndex || moveRightKeyIndex < 0)
        break;

      if (conf.getQuantLevel() == 2) {

        if (conf.isHighRes())
          ;
        else
          result[rightIndex] = readSpectrum(sModel, moveRightKeyIndex);
      } else {
        if (conf.isHighRes())
          result[rightIndex] = readFullSpectrum(moveRightKeyIndex, null, null, samIsoArr, mapModel, null, resultHt);
        else
          result[rightIndex] = readFullSpectrum(moveRightKeyIndex, null, null, null, mapModel, sModel, resultHt);
      }

      if (0 == result[rightIndex][0]) {
        if (conf.getQuantLevel() == 2)
          moveRightKeyIndex = findLightKey(sModel.getIFile().getKeys(), moveRightKeyIndex, MOVE_RIGHT, sModel.getIFile(), conf);
        else
          moveRightKeyIndex++;

        continue;
      }

      rightIndex++;

      if (conf.getQuantLevel() == 2)
        moveRightKeyIndex = findLightKey(sModel.getIFile().getKeys(), moveRightKeyIndex, MOVE_RIGHT, sModel.getIFile(), conf);
      else
        moveRightKeyIndex++;
    }


        /*
        System.out.println("===============");
        for(int iii=0;iii<result.length;iii++)
        {
            //if(result[iii][0]==0)
              //  continue;

            for(int jjj=0;jjj<result[iii].length;jjj++)
                System.out.print(result[iii][jjj] + " ");

            System.out.println("");
        }
        */


    if (conf.getQuantLevel() == 1)
      return buildNLResult(peakStart, peakEnd, leftIndex, rightIndex, result, pathFileNameList.size());
    else //if(conf.getQuantLevel()==2)
      return buildNLResult(peakStart, peakEnd, leftIndex, rightIndex, result, sModel.getBioSample().length, pathFileNameList.size());
  }


  public final static int MOVE_LEFT = 1;
  public final static int MOVE_RIGHT = 2;

  public static int findNextKeyIndex(int keys[], int keyIndex, int move, IndexedFile iFile, Configuration conf, double avgMass, boolean findHeavy) //find heavy ion or just move to next window
  {
    TIntDoubleHashMap precurMap = iFile.getPrecursorMap();
    TIntLongHashMap byteMap = iFile.getMsIndex();


  ///  System.out.println("-----------------\t" + conf.getCalcSamAvgMass());

    //  System.out.println("keyindex=\t" + keyIndex);

    double tolerance = conf.getIsolationWindow() / 2;


    //System.out.println("===\t" + avgMass);

    double startMass = avgMass - tolerance;
    double endMass = avgMass + tolerance;

    int startLeftIndex = -1;
    int startRightIndex = -1;

    switch (conf.getExpType()) {
      //    case CensusConstants.MRM_EXPERIMENT :
      //System.out.println(finHeavy);
      case CensusConstants.MSMS_DATA_INDEPENDENT:
      case CensusConstants.MSMS_DATA_INDEPENDENT_LFREE:
        if (findHeavy) {
          startLeftIndex = keyIndex;
          startRightIndex = keyIndex;
        } else {
          startLeftIndex = keyIndex - 1;
          startRightIndex = keyIndex + 1;

        }
        break;

      default:
        startLeftIndex = keyIndex - 1;
        startRightIndex = keyIndex + 1;

        break;
    }

    if (move == MOVE_LEFT) {
      for (int i = startLeftIndex; i >= 0; i--) {
        double prec = precurMap.get(keys[i]);  //retrieve precursor

    //    System.out.println("LL=====\t" + prec + " " + keys[i] + "\t" + startMass + " " + endMass);

        if (startMass <= prec && endMass >= prec) {
      //    System.out.println("L=====\t" + prec + " " + keys[i] + "\t" + startMass + " " + endMass);
          return i;
        }
      }

    } else //MOVE_RIGHT
    {
      for (int i = startRightIndex; i < keys.length; i++) {
        double prec = precurMap.get(keys[i]);  //retrieve precursor

     //   System.out.println("RR=====\t" + prec + " " + keys[i] + "\t" + startMass + " " + endMass);

        if (startMass <= prec && endMass >= prec) {
       //   System.out.println("R=====\t" + prec + " " + keys[i] + "\t" + startMass + " " + endMass);
          return i;
        }
      }
    }


    return -1;
  }

  public static int findHeavyKey(int keys[], int keyIndex, IndexedFile iFile, Configuration conf, double avgMass) {
    return findNextKeyIndex(keys, keyIndex, MOVE_RIGHT, iFile, conf, avgMass, true);
  }

  public static int findLightKey(int keys[], int keyIndex, int move, IndexedFile iFile, Configuration conf) //throws PrecursorNotFoundException
  {
    if (conf.getQuantLevel() == 1)
      return keyIndex + 1;
    //findNextKey(keys, keyIndex, MOVE_RIGHT, sModel.getIFile(), conf, conf.getCalcSamAvgMass());


 //   System.out.println("-------->>>>>---------\t" + conf.getCalcSamAvgMass());
    int key = findNextKeyIndex(keys, keyIndex, move, iFile, conf, conf.getCalcSamAvgMass(), false);

    return key;

  }

  //for mrm data including iTRAQ multiple spectra
  public static int findMRMNextKey(int keys[], int keyIndex, int move, IndexedFile iFile, Configuration conf) //throws PrecursorNotFoundException
  {

    int key = findMRMNextKeyIndex(keys, keyIndex, move, iFile, conf);

    return key;
  }


  public static int findMRMNextKeyIndex(int keys[], int keyIndex, int move, IndexedFile iFile, Configuration conf) {
    if (keyIndex < 0 || keyIndex >= keys.length)
      return -1;

    TIntDoubleHashMap precurMap = iFile.getPrecursorMap();
    TIntLongHashMap byteMap = iFile.getMsIndex();

    double tolerance = conf.MRM_PRECURSOR_TOLERANCE;

    //double startMass = conf.getCalcSamAvgMass() - tolerance;
    //double endMass = conf.getCalcSamAvgMass() + tolerance;

    double precursorMass = precurMap.get(keys[keyIndex]);
    double startMass = precursorMass - tolerance;
    double endMass = precursorMass + tolerance;

    if (move == MOVE_LEFT) {
      for (int i = keyIndex - 1; i >= 0; i--) {
        double prec = precurMap.get(keys[i]);  //retrieve precursor

        if (startMass <= prec && endMass >= prec)
          return i;
      }

    } else //MOVE_RIGHT
    {
      for (int i = keyIndex + 1; i < keys.length; i++) {
        double prec = precurMap.get(keys[i]);  //retrieve precursor


        if (startMass <= prec && endMass >= prec)
          return i;
      }

    }

    return -1;
  }

  public static int findMRMNextKeyIndexByPrecursor(double rtArr[], int keyIndex, int move, IndexedFile iFile, Configuration conf, double pCursor) {
    if (keyIndex < 0 || keyIndex >= rtArr.length)
      if (MOVE_RIGHT == move)
        return rtArr.length;
      else
        return -1;

    TDoubleDoubleHashMap precurMap = iFile.getRtPrecursorMap();

    //double startMass = conf.getCalcSamAvgMass() - tolerance;
    //double endMass = conf.getCalcSamAvgMass() + tolerance;

    double precursorMass = precurMap.get(rtArr[keyIndex]);
    //double startMass = precursorMass - tolerance;

    if (move == MOVE_LEFT) {
      for (int i = keyIndex - 1; i >= 0; i--) {
        double prec = precurMap.get(rtArr[i]);  //retrieve precursor

        if (pCursor == prec) //startMass<=prec && endMass>=prec)
          return i;
      }
    } else //MOVE_RIGHT
    {
      for (int i = keyIndex + 1; i < rtArr.length; i++) {
        double prec = precurMap.get(rtArr[i]);  //retrieve precursor

        if (pCursor == prec) //if(startMass<=prec && endMass>=prec)
          return i;
      }
    }

    if (MOVE_RIGHT == move)
      return rtArr.length;

    return -1;
  }


  //MSMS peak finding for iTRAQ
  public static String peakFindingMSMSMultipleSpecific(SpectrumModel sModel, SpecRange range, Configuration conf, int keyIndex)
    throws PrecursorNotFoundException, IOException, CensusIndexOutOfBoundException {

    TIntLongHashMap index = sModel.getIndex();

    RandomAccessFile file = sModel.getFile();
    steepRatioThreshold = conf.getSteepRatioThreshold();

//	int numIsoWindow = conf.getNumOfIsolationWindow();
    int maxWindow = conf.getMaxWindow();
    int margin = conf.getMargin();
    List massList = sModel.getMsmsSpecificMassList();

    double[][] result = new double[maxWindow * 2 + 1 + margin * 2][massList.size() + 1];

    int leftIndex = maxWindow + margin;
    int rightIndex = maxWindow + margin + 1;


    double totalIntensity = 0;
    int[] keys = sModel.getKeys();
    int initWin = 2;

    int steepArea = conf.getSteepArea();

    int moveLeftKeyIndex = keyIndex;
    int moveRightKeyIndex = findMRMNextKey(keys, keyIndex, MOVE_RIGHT, sModel.getIFile(), conf);

    IndexedFile iFile = sModel.getIFile();

    for (int i = 0; i < initWin; i++) ///Index=maxWindow;keys[leftPeakIndex]>=leftMinValue; )
    {
      if (moveLeftKeyIndex <= 0 || leftIndex <= 0) {
        //move back to right because it cannot move to left any more
        moveLeftKeyIndex = findMRMNextKey(keys, moveLeftKeyIndex, MOVE_RIGHT, sModel.getIFile(), conf);
        leftIndex++;

        break;
      }

      result[leftIndex] = readSpectrumSpecificPeaks(moveLeftKeyIndex, sModel.getIFile(), massList, conf);

      for (int j = 0; j < result[leftIndex].length; j++)
        totalIntensity += result[leftIndex][j];

      moveLeftKeyIndex = findMRMNextKey(keys, moveLeftKeyIndex, MOVE_LEFT, sModel.getIFile(), conf);

      leftIndex--;
    }

    for (int i = 0; i < initWin; i++) //leftIndex=maxWindow;keys[leftPeakIndex]>=leftMinValue; )
    {
      //System.out.println(moveRightKeyIndex + " " + keys.length + " ---" + i);
      if (moveRightKeyIndex >= keys.length || moveRightKeyIndex < 0) {
        //      System.out.println(moveRightKeyIndex + " " + keys.length + " " + i);
        //    System.out.println(moveRightKeyIndex + " " + keys.length + " " + i);
        //  System.out.println(moveRightKeyIndex + " " + keys.length + " " + i);
        //  System.out.println(moveRightKeyIndex + " " + keys.length + " " + i);


        moveRightKeyIndex = findMRMNextKey(keys, moveRightKeyIndex, MOVE_LEFT, sModel.getIFile(), conf);

        rightIndex--;

        break;
      }

      result[rightIndex] = readSpectrumSpecificPeaks(moveRightKeyIndex, sModel.getIFile(), massList, conf);

      for (int j = 0; j < result[rightIndex].length; j++)
        totalIntensity += result[rightIndex][j];


      //moveRightKeyIndex += numIsoWindow;
      moveRightKeyIndex = findMRMNextKey(keys, moveRightKeyIndex, MOVE_RIGHT, sModel.getIFile(), conf);
      rightIndex++;
    }

    if (moveRightKeyIndex < 0 && moveLeftKeyIndex < 0)
      return null;

    boolean isGoingUp = true;
    boolean isHighIntensity = true;  //if the intensity is lower than one third of average intensity, this becomes false


    double[][] arr = new double[steepArea][massList.size() + 1];
    double[][] prevArr = new double[steepArea][massList.size() + 1];

    //steepRatioThreshold = 0.1f;  //give very strigent threshold for iTRAQ
    if (leftIndex + steepArea >= 0) {
      //for(int i=0;i<300;i++)
      while (true) {
        double area1 = 0;
        double area2 = 0;

        if (moveLeftKeyIndex - steepArea < 0) {
          break;
        }

        int steepCount = 0;
        int tempKeyIndex = moveLeftKeyIndex;
        while (true) {
          if (tempKeyIndex < 0 || steepCount < 0 || steepCount >= steepArea - 1)
            break;

          if (prevArr[0][0] != 0) {
            for (int l = 0; l < 2; l++) {
              //                  arr[l] = prevArr[l+1];
              for (int m = 0; m < 3; m++) {
                arr[l][m] = prevArr[l + 1][m];

              }
            }
          } else {
            arr[steepCount] = readSpectrumSpecificPeaks(tempKeyIndex, sModel.getIFile(), massList, conf);
            //arr[steepCount] = readSpectrum(keys, tempKeyIndex, index, findHeavyKey(keys, tempKeyIndex, iFile, conf, conf.getCalcRefAvgMass()), file, sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());
          }

          steepCount++;
          tempKeyIndex = findMRMNextKey(keys, tempKeyIndex, MOVE_LEFT, sModel.getIFile(), conf);
        }

        if (tempKeyIndex < 0)
          break;

        arr[steepCount] = readSpectrumSpecificPeaks(tempKeyIndex, sModel.getIFile(), massList, conf);
        //readSpectrum(keys, tempKeyIndex, index, findHeavyKey(keys, tempKeyIndex, iFile, conf, conf.getCalcRefAvgMass()), file, sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());

        for (int k = 0; k < 3; k++) {
          for (int l = 0; l < massList.size() + 1; l++)
            prevArr[k][l] = arr[k][l];
        }

        for (int j = leftIndex + 1; j <= leftIndex + steepCount + 1; j++) {
          for (int k = 1; k < result[j].length; k++) {
            area1 += result[j][k];
          }
        }

        for (int j = 0; j <= steepCount; j++) {
          for (int k = 1; k < result[j].length; k++) {
            area2 += arr[j][k];
          }
        }
        //            System.out.println("");
        //          System.out.println(area2 + " " + area1 + " " + (double)area2/area1);
        //        System.out.println("");

        if (area2 == 0 && area1 == 0)
          break;

        double areaRatio = (double) area2 / area1;

        /*****************************************
         * AREA1 is right side area
         * AREA2 is left side area
         *****************************************/
        double intHeight = (double) totalIntensity / (rightIndex - leftIndex - 1) / 1.5;

        if (isHighIntensity && (intHeight > area2 / 3))
          isHighIntensity = false;


        if (areaRatio < 0.01)
          isGoingUp = false;

        //This is comparison between most intensive peak vs. new peak
        double heightRatio = area2 / steepArea / intHeight;

        //if( (isGoingUp && isHighIntensity) || areaRatio<steepRatioThreshold )
        if (heightRatio > 0.01 && ((isGoingUp || isHighIntensity) || areaRatio < steepRatioThreshold)) {
          //                System.out.println("2222");
          if (leftIndex <= 0) {
            isGoingUp = false;
            break;
          }

          result[leftIndex] = arr[0];
          moveLeftKeyIndex = findMRMNextKey(keys, moveLeftKeyIndex, MOVE_LEFT, sModel.getIFile(), conf);

          if (0 == result[leftIndex][0])
            continue;

          leftIndex--;

          for (int j = 0; j < steepCount; j++)
            arr[j] = arr[j + 1];

          for (int j = 1; j < arr[0].length; j++) {
            totalIntensity += arr[0][j];
          }
          //            System.out.print("\n");
          continue;
        } else
          isGoingUp = false;

        //      System.out.println("00000000000 " + areaRatio + " " + steepRatioThreshold + " " + isGoingUp);
        if (areaRatio > steepRatioThreshold && !isGoingUp) {
          for (int k = 0; k < 3; k++) {
            if (result[leftIndex + 1][1] < arr[k][1] || arr[k][0] == 0)
              break;

            if (leftIndex < 0)
              break;

            result[leftIndex] = arr[k];

            //moveLeftKeyIndex -= numIsoWindow;
            moveLeftKeyIndex = findMRMNextKey(keys, moveLeftKeyIndex, MOVE_LEFT, sModel.getIFile(), conf);

            leftIndex--;
          }

          break;
        }

        //System.out.println("000000000000");

        result[leftIndex--] = arr[0];
        moveLeftKeyIndex = findMRMNextKey(keys, moveLeftKeyIndex, MOVE_LEFT, sModel.getIFile(), conf);

        for (int j = 0; j < steepCount; j++)
          arr[j] = arr[j + 1];

        //System.out.println("000000000000");
      }
    }

    isGoingUp = true;
    isHighIntensity = true;

    prevArr = new double[steepArea][massList.size() + 1]; //clean up the array

    if (rightIndex + steepArea < result.length) {
      //for(int i=0;i<500;i++)
      while (true) {
        double area1 = 0;
        double area2 = 0;

        if (moveRightKeyIndex + steepArea >= keys.length)
          break;

        int steepCount = 0;

        int tempKeyIndex = moveRightKeyIndex;
        while (true) {
          if (tempKeyIndex >= keys.length || steepCount >= steepArea - 1 || steepCount >= keys.length || tempKeyIndex < 0)
            break;

          if (prevArr[0][0] != 0) {
            for (int l = 0; l < 2; l++)
              for (int m = 0; m < 3; m++)
                arr[l][m] = prevArr[l + 1][m];
          } else {
            //arr[steepCount] = readSpectrum(keys, tempKeyIndex, index, findHeavyKey(keys, tempKeyIndex, iFile, conf, conf.getCalcRefAvgMass()), file, sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());

            arr[steepCount] = readSpectrumSpecificPeaks(tempKeyIndex, sModel.getIFile(), massList, conf);
          }

          steepCount++;
          //tempKeyIndex = tempKeyIndex + 1*numIsoWindow;

          tempKeyIndex = findMRMNextKey(keys, tempKeyIndex, MOVE_RIGHT, sModel.getIFile(), conf);
          //tempKeyIndex++;
        }

        if (tempKeyIndex >= keys.length)
          break;

        arr[steepCount] = readSpectrumSpecificPeaks(tempKeyIndex, sModel.getIFile(), massList, conf);

        for (int k = 0; k < 3; k++) {
          for (int l = 0; l < massList.size() + 1; l++)
            prevArr[k][l] = arr[k][l];
        }


        for (int j = rightIndex - 1; j > rightIndex - steepCount - 2; j--) {
          for (int k = 1; k < result[j].length; k++) {
            area1 += result[j][k];
          }
        }

        for (int j = 0; j <= steepCount; j++) {
          for (int k = 1; k < result[j].length; k++) {
            area2 += arr[j][k];
          }
        }

        if (area2 == 0 && area1 == 0)
          break;

        double areaRatio = (double) area2 / area1;

        /*****************************************
         * AREA1 is left side area
         * AREA2 is right side area
         *****************************************/

        double intHeight = (double) totalIntensity / (rightIndex - leftIndex - 1) / 1.5;
        if (isHighIntensity && (intHeight > area2 / 3))
          isHighIntensity = false;

        //if( (isGoingUp && isHighIntensity) || areaRatio<steepRatioThreshold )


        if (areaRatio < 0.01)
          isGoingUp = false;

        //This is comparison between most intensive peak vs. new peak
        double heightRatio = area2 / steepArea / intHeight;
        if (heightRatio > 0.01 && ((isGoingUp || isHighIntensity) || areaRatio < steepRatioThreshold)) {
          if ((rightIndex + 1) >= result.length || moveRightKeyIndex < 0) {
            isGoingUp = false;
            break;
          }

          result[rightIndex] = arr[0];
          //System.out.println("000>>" + arr[0][0]);

//		    moveRightKeyIndex += numIsoWindow;
          moveRightKeyIndex = findMRMNextKey(keys, moveRightKeyIndex, MOVE_RIGHT, sModel.getIFile(), conf);


          if (0 == result[rightIndex][0])
            continue;

          rightIndex++;

          for (int j = 0; j < steepCount; j++)
            arr[j] = arr[j + 1];

          for (int j = 0; j < arr[0].length; j++)
            totalIntensity += arr[0][j];


          continue;
        } else
          isGoingUp = false;

        if (area2 / area1 > steepRatioThreshold) {
          for (int k = 0; k < 3; k++) {
            if (result[rightIndex - 1][1] < arr[k][1])
              break;

            if (result.length <= rightIndex)
              break;

            result[rightIndex] = arr[k];

            //System.out.println(arr[0][1] + " " + arr[1][1] + " " + arr[2][1] + " " + result[rightIndex-1][1]);

            //moveRightKeyIndex += numIsoWindow;
            moveRightKeyIndex = findMRMNextKey(keys, moveRightKeyIndex, MOVE_RIGHT, sModel.getIFile(), conf);

            rightIndex++;
          }
          break;
        }

        result[rightIndex++] = arr[0];

        //moveRightKeyIndex += numIsoWindow;
        moveRightKeyIndex = findMRMNextKey(keys, moveRightKeyIndex, MOVE_RIGHT, sModel.getIFile(), conf);

        for (int j = 0; j < steepCount; j++)
          arr[j] = arr[j + 1];
      }
    }


    int peakStart = leftIndex + 1;
    int peakEnd = rightIndex - 1;
    int dtaStart = range.getMin();
    int dtaEnd = range.getMax();


    int specSpace = (int) (result[peakStart + 1][0] - result[peakStart][0]);

    if (specSpace == 0)
      specSpace = (int) (result[peakStart][0] - result[peakStart - 1][0]);

    int diff = (int) (result[peakStart][0] - dtaStart) / specSpace;
    int leftMargin = margin;

    int rightMargin = margin;
    if (diff > 0)
      leftMargin += diff;

    diff = (int) (dtaEnd - result[rightIndex - 1][0]) / specSpace;

    if (diff > 0)
      rightMargin += diff;


    for (int i = 0; i < leftMargin; i++) {
      if (leftIndex < 0 || moveLeftKeyIndex < 0)
        break;

      //if(conf.isDataIndependent()) //data independent
      result[leftIndex] = readSpectrumSpecificPeaks(moveLeftKeyIndex, sModel.getIFile(), massList, conf);

      if (0 == result[leftIndex][0]) {
        //moveLeftKeyIndex -= numIsoWindow;
        moveLeftKeyIndex = findMRMNextKey(keys, moveLeftKeyIndex, MOVE_LEFT, sModel.getIFile(), conf);

        continue;
      }


      //moveLeftKeyIndex -= numIsoWindow;
      moveLeftKeyIndex = findMRMNextKey(keys, moveLeftKeyIndex, MOVE_LEFT, sModel.getIFile(), conf);

      if (moveLeftKeyIndex < 0)
        break;

      leftIndex--;
    }
    for (int i = 0; i < rightMargin; i++) {
      if (rightIndex >= result.length || moveRightKeyIndex >= keys.length || moveRightKeyIndex < 0)
        break;

      result[rightIndex] = readSpectrumSpecificPeaks(moveRightKeyIndex, sModel.getIFile(), massList, conf);

      if (0 == result[rightIndex][0]) {
        //moveRightKeyIndex += numIsoWindow;
        moveRightKeyIndex = findMRMNextKey(keys, moveRightKeyIndex, MOVE_RIGHT, sModel.getIFile(), conf);

        continue;
      }

      moveRightKeyIndex = findMRMNextKey(keys, moveRightKeyIndex, MOVE_RIGHT, sModel.getIFile(), conf);

      if (moveRightKeyIndex < 0)
        break;

      rightIndex++;
    }


        /*
        for(int i=0;i<result.length;i++)
        {
            for(int j=0;j<result[i].length;j++)
            {
                System.out.print(result[i][j] + " ");

            }

            System.out.println(" ");
        }


*/

    return buildMsmsSpecificResult(peakStart, peakEnd, leftIndex, rightIndex, result, massList);
  }

  //e.g. iTRAQ throughout the time
  private static String buildMsmsSpecificResult(int peakStart, int peakEnd, int leftIndex, int rightIndex, double[][] result, List massList) {
    StringBuffer sb = new StringBuffer();
    sb.append("P ").append((long) result[peakStart][0]).append(" ").append((long) result[peakEnd][0]).append(";");

    for (Iterator itr = massList.iterator(); itr.hasNext(); ) {
      sb.append(itr.next()).append(" ");
    }

    sb.setCharAt(sb.length() - 1, ';');

    for (int i = leftIndex + 1; i < rightIndex - 1; i++) {
      //spectrum number
      sb.append((int) result[i][0]).append(" ");

      for (int j = 1; j < result[i].length - 1; j++) {
        sb.append((long) result[i][j]).append(" ");
      }

      sb.append((long) result[i][result[i].length - 1]).append(";");
    }

    return sb.toString();
  }

  //MSMS peak finding for MRM
  public static String peakFindingMRM(MRMPeptideModel model, Configuration conf, IndexedFile iFile)
    throws PrecursorNotFoundException, IOException, CensusIndexOutOfBoundException, Exception {
    RandomAccessFile file = iFile.getFile();
    steepRatioThreshold = conf.getSteepRatioThreshold();
    int maxWindow = conf.getMaxWindow();
    int margin = conf.getMargin();

    double[][] result = new double[maxWindow * 2 + 1 + margin * 2][4 * model.getBionArr().length + 3]; //scan #, sample intensity, ref intensity
    int leftIndex = maxWindow + margin;
    int rightIndex = maxWindow + margin + 1;

    double totalIntensity = 0;
    int[] keys = iFile.getKeys(); //sModel.getKeys();
    int initWin = 2;

    int steepArea = conf.getSteepArea();

    int keyIndex = model.getKeyIndex();

    //iFile.getPrecursorMap()
    int moveLeftKeyIndex = keyIndex;
    double[] rtArr = iFile.getRtArr();

    for (int i = 0; i < initWin; i++) ///Index=maxWindow;keys[leftPeakIndex]>=leftMinValue; )
    {
      if (moveLeftKeyIndex <= 0 || leftIndex <= 0) {
        //move back to right because it cannot move to left any more
        //= findLightKey(keys, moveLeftKeyIndex, MOVE_RIGHT, sModel.getIFile(), conf);
        moveLeftKeyIndex = findMRMNextKeyIndexByPrecursor(rtArr, keyIndex, MOVE_RIGHT, iFile, conf, model.getParentMass());

        leftIndex++;

        break;
      }

      result[leftIndex] = readSpectrumByRt(rtArr, moveLeftKeyIndex, iFile, model);
      //readSpectrum(keys, moveLeftKeyIndex, index, findHeavyKey(keys, moveLeftKeyIndex, iFile, conf, conf.getCalcRefAvgMass()), file, sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());
      totalIntensity += result[leftIndex][1];
      System.out.println(iFile.getScanNumByRtIndex(moveLeftKeyIndex) + " " + totalIntensity + " " + result[leftIndex][1]);

      moveLeftKeyIndex = findMRMNextKeyIndexByPrecursor(rtArr, moveLeftKeyIndex, MOVE_LEFT, iFile, conf, model.getParentMass());

      leftIndex--;
    }

    int moveRightKeyIndex = findMRMNextKeyIndexByPrecursor(rtArr, keyIndex, MOVE_RIGHT, iFile, conf, model.getParentMass());

    for (int i = 0; i < initWin; i++) //leftIndex=maxWindow;keys[leftPeakIndex]>=leftMinValue; )
    {
      if (moveRightKeyIndex >= rtArr.length || moveRightKeyIndex < 0) {
        moveRightKeyIndex = findMRMNextKeyIndexByPrecursor(rtArr, moveRightKeyIndex, MOVE_LEFT, iFile, conf, model.getParentMass());
        rightIndex--;
        break;
      }

      result[rightIndex] = readSpectrumByRt(rtArr, moveRightKeyIndex, iFile, model);

      totalIntensity += result[rightIndex][1];
      //totalIntensity += result[rightIndex][2];

//System.out.println(iFile.getScanNumByRtIndex(moveRightKeyIndex) + " " + totalIntensity + " " + result[rightIndex][1]);
      moveRightKeyIndex = findMRMNextKeyIndexByPrecursor(rtArr, moveRightKeyIndex, MOVE_RIGHT, iFile, conf, model.getParentMass());
      rightIndex++;
    }

    if (moveRightKeyIndex < 0 && moveLeftKeyIndex < 0)
      return null;

    boolean isGoingUp = true;
    boolean isHighIntensity = true;  //if the intensity is lower than one third of average intensity, this becomes false
    double[][] arr = new double[steepArea][3];
    double[][] prevArr = new double[steepArea][3];

    if (leftIndex + steepArea >= 0) {
      //for(int i=0;i<300;i++)
      while (true) {
        double area1 = 0;
        double area2 = 0;

        if (moveLeftKeyIndex - steepArea < 0) {
          break;
        }

        int steepCount = 0;
        int tempKeyIndex = moveLeftKeyIndex;
        while (true) {
          if (tempKeyIndex < 0 || steepCount < 0 || steepCount >= steepArea - 1)
            break;

          if (prevArr[0][0] != 0) {
            for (int l = 0; l < 2; l++)
              for (int m = 0; m < 3; m++)
                arr[l][m] = prevArr[l + 1][m];
          } else {
            arr[steepCount] = readSpectrumByRt(rtArr, tempKeyIndex, iFile, model);
          }

          steepCount++;

          tempKeyIndex = findMRMNextKeyIndexByPrecursor(rtArr, tempKeyIndex, MOVE_LEFT, iFile, conf, model.getParentMass());
        }

        if (tempKeyIndex < 0)
          break;

        arr[steepCount] = readSpectrumByRt(rtArr, tempKeyIndex, iFile, model);

        for (int k = 0; k < 3; k++) {
          for (int l = 0; l < 3; l++)
            prevArr[k][l] = arr[k][l];
        }

        for (int j = leftIndex + 1; j <= leftIndex + steepCount + 1; j++) {
          area1 += result[j][1];
          area1 += result[j][2];
        }

        for (int j = 0; j <= steepCount; j++) {
          area2 += arr[j][1];
          area2 += arr[j][2];
        }

        if (area2 == 0 && area1 == 0)
          break;

        double areaRatio = (double) area2 / area1;

        /*****************************************
         * AREA1 is right side area
         * AREA2 is left side area
         *****************************************/
        if (isHighIntensity && ((double) totalIntensity / (rightIndex - leftIndex - 1) / 1.5 > area2 / 3))
          isHighIntensity = false;

        if ((isGoingUp && isHighIntensity) || areaRatio < steepRatioThreshold) {
          if (leftIndex <= 0) {
            isGoingUp = false;
            break;
          }

          result[leftIndex] = arr[0];


          moveLeftKeyIndex = findMRMNextKeyIndexByPrecursor(rtArr, moveLeftKeyIndex, MOVE_LEFT, iFile, conf, model.getParentMass());

          if (0 == result[leftIndex][0])
            continue;

          leftIndex--;

          for (int j = 0; j < steepCount; j++)
            arr[j] = arr[j + 1];

          totalIntensity += arr[0][1];
          totalIntensity += arr[0][2];

          continue;
        } else
          isGoingUp = false;


        if (areaRatio > steepRatioThreshold) {
          for (int k = 0; k < 3; k++) {
            if (result[leftIndex + 1][1] < arr[k][1] || arr[k][0] == 0)
              break;

            if (leftIndex < 0)
              break;

            result[leftIndex] = arr[k];

            moveLeftKeyIndex = findMRMNextKeyIndexByPrecursor(rtArr, moveLeftKeyIndex, MOVE_LEFT, iFile, conf, model.getParentMass());

            leftIndex--;
          }

          break;
        }

        result[leftIndex--] = arr[0];

        moveLeftKeyIndex = findMRMNextKeyIndexByPrecursor(rtArr, moveLeftKeyIndex, MOVE_LEFT, iFile, conf, model.getParentMass());

        for (int j = 0; j < steepCount; j++)
          arr[j] = arr[j + 1];
      }
    }

    isGoingUp = true;
    isHighIntensity = true;

    prevArr = new double[steepArea][3]; //clean up the array

    if (rightIndex + steepArea < result.length) {
      //for(int i=0;i<500;i++)
      while (true) {

        double area1 = 0;
        double area2 = 0;

        if (moveRightKeyIndex + steepArea >= rtArr.length)
          break;

        int steepCount = 0;

        int tempKeyIndex = moveRightKeyIndex;
        while (true) {
          if (tempKeyIndex >= rtArr.length || steepCount >= steepArea - 1 || steepCount >= rtArr.length)
            break;

          if (prevArr[0][0] != 0) {
            for (int l = 0; l < 2; l++)
              for (int m = 0; m < 3; m++)
                arr[l][m] = prevArr[l + 1][m];
          } else {
            arr[steepCount] = readSpectrumByRt(rtArr, tempKeyIndex, iFile, model);
          }

          steepCount++;
          //tempKeyIndex = tempKeyIndex + 1*numIsoWindow;
          tempKeyIndex = findMRMNextKeyIndexByPrecursor(rtArr, tempKeyIndex, MOVE_RIGHT, iFile, conf, model.getParentMass());
          //tempKeyIndex++;
        }

        if (tempKeyIndex >= rtArr.length)
          break;

        arr[steepCount] = readSpectrumByRt(rtArr, tempKeyIndex, iFile, model);

        for (int k = 0; k < 3; k++) {
          for (int l = 0; l < 3; l++)
            prevArr[k][l] = arr[k][l];
        }

        for (int j = rightIndex - 1; j > rightIndex - steepCount - 2; j--) {
          area1 += result[j][1];
          area1 += result[j][2];
        }

        for (int j = 0; j <= steepCount; j++) {
          area2 += arr[j][1];
          area2 += arr[j][2];
        }

        if (area2 == 0 && area1 == 0)
          break;

        double areaRatio = (double) area2 / area1;

        /*****************************************
         * AREA1 is left side area
         * AREA2 is right side area
         *****************************************/

        if (isHighIntensity && ((double) totalIntensity / (rightIndex - leftIndex - 1) / 1.5 > area2 / 3))
          isHighIntensity = false;

        //if( ((isGoingUp && (double)totalIntensity/(rightIndex-leftIndex-1)/3<area2/3)) || areaRatio<steepRatioThreshold )
        if ((isGoingUp && isHighIntensity) || areaRatio < steepRatioThreshold) {
          if ((rightIndex + 1) >= result.length) {
            isGoingUp = false;
            break;
          }

          result[rightIndex] = arr[0];
          moveRightKeyIndex = findMRMNextKeyIndexByPrecursor(rtArr, moveRightKeyIndex, MOVE_RIGHT, iFile, conf, model.getParentMass());

          if (0 == result[rightIndex][0])
            continue;

          rightIndex++;

          for (int j = 0; j < steepCount; j++)
            arr[j] = arr[j + 1];

          //rightTotalIntensity += arr[0][1];
          totalIntensity += arr[0][1];
          totalIntensity += arr[0][2];

          continue;
        } else
          isGoingUp = false;

        if (area2 / area1 > steepRatioThreshold) {
          for (int k = 0; k < 3; k++) {
            if (result[rightIndex - 1][1] < arr[k][1])
              break;

            if (result.length <= rightIndex)
              break;

            result[rightIndex] = arr[k];

            //System.out.println(arr[0][1] + " " + arr[1][1] + " " + arr[2][1] + " " + result[rightIndex-1][1]);

            //moveRightKeyIndex += numIsoWindow;
            moveRightKeyIndex = findMRMNextKeyIndexByPrecursor(rtArr, moveRightKeyIndex, MOVE_RIGHT, iFile, conf, model.getParentMass());

            rightIndex++;
          }
          break;
        }

        result[rightIndex++] = arr[0];

        //moveRightKeyIndex += numIsoWindow;
        moveRightKeyIndex = findMRMNextKeyIndexByPrecursor(rtArr, moveRightKeyIndex, MOVE_RIGHT, iFile, conf, model.getParentMass());

        for (int j = 0; j < steepCount; j++)
          arr[j] = arr[j + 1];

        //System.exit(0);

      }
    }
    int peakStart = leftIndex + 1;
    int peakEnd = rightIndex - 1;

    int leftMargin = margin;
    int rightMargin = margin;


    for (int i = 0; i < leftMargin; i++) {
      if (leftIndex < 0 || moveLeftKeyIndex < 0)
        break;

      result[leftIndex] = readSpectrumByRt(rtArr, moveLeftKeyIndex, iFile, model);

      if (0 == result[leftIndex][0]) {
        //moveLeftKeyIndex -= numIsoWindow;
        moveLeftKeyIndex = findMRMNextKeyIndexByPrecursor(rtArr, moveLeftKeyIndex, MOVE_LEFT, iFile, conf, model.getParentMass());
        continue;
      }

      //moveLeftKeyIndex -= numIsoWindow;
      moveLeftKeyIndex = findMRMNextKeyIndexByPrecursor(rtArr, moveLeftKeyIndex, MOVE_LEFT, iFile, conf, model.getParentMass());

      if (moveLeftKeyIndex < 0)
        break;

      leftIndex--;
    }
    for (int i = 0; i < rightMargin; i++) {
      if (rightIndex >= result.length || moveRightKeyIndex >= rtArr.length)
        break;

      result[rightIndex] = readSpectrumByRt(rtArr, moveRightKeyIndex, iFile, model);

      if (0 == result[rightIndex][0]) {
        //moveRightKeyIndex += numIsoWindow;
        moveRightKeyIndex = findMRMNextKeyIndexByPrecursor(rtArr, moveRightKeyIndex, MOVE_RIGHT, iFile, conf, model.getParentMass());
        continue;
      }

      //moveRightKeyIndex += numIsoWindow;
      moveRightKeyIndex = findMRMNextKeyIndexByPrecursor(rtArr, moveRightKeyIndex, MOVE_RIGHT, iFile, conf, model.getParentMass());

      if (moveRightKeyIndex < 0)
        break;

      rightIndex++;
    }

  /*
        for(int i=leftIndex;i<=rightIndex;i++)
        {
            for(int j=0;j<5;j++)
                System.out.print(result[i][j] + "\t");

            System.out.println("");
        }

            System.out.println("");
        for(int i=peakStart;i<=peakEnd;i++)
        {
            for(int j=0;j<5;j++)
                System.out.print(result[i][j] + "\t");

            System.out.println("");
        }
*/

    return buildMRMResult(peakStart, peakEnd, leftIndex, rightIndex, result, model.getBionArr().length);
  }


  //MSMS peak finding.  labelfree. data independent
  public static String peakFindingMSMSLabelfree(SpectrumModel sModel, SpecRange range, Configuration conf,
                                                int keyIndex, TIntDoubleHashMap scanToRetMap)
    throws PrecursorNotFoundException, IOException, CensusIndexOutOfBoundException, Exception {
    TIntLongHashMap index = sModel.getIndex();
    RandomAccessFile file = sModel.getFile();
    steepRatioThreshold = conf.getSteepRatioThreshold();

//	int numIsoWindow = conf.getNumOfIsolationWindow();
    int maxWindow = conf.getMaxWindow();
    int margin = conf.getMargin();

    double[][] result = null;
    //conf.getQuantLevel()



   // if (conf.getQuantLevel() == 2)
      //result = new double[maxWindow * 2 + 1 + margin * 2][4 * sModel.getBioSample().length + 3]; //scan #, sample intensity, ref intensity
      result = new double[maxWindow * 2 + 1 + margin * 2][4 * sModel.getBioSample().length + 3 + 1]; //scan #, ret time, sample intensity, ref intensity
   // else
   //   result = new double[maxWindow * 2 + 1 + margin * 2][3];

    int leftIndex = maxWindow + margin;
    int rightIndex = maxWindow + margin + 1;

    //double leftTotalIntensity=0;
    //double rightTotalIntensity=0;
    double totalIntensity = 0;
    int[] keys = sModel.getKeys();
    int initWin = 2;

    int steepArea = conf.getSteepArea();

    int moveLeftKeyIndex = keyIndex;
    int moveRightKeyIndex = findLightKey(keys, keyIndex, MOVE_RIGHT, sModel.getIFile(), conf);

    IndexedFile iFile = sModel.getIFile();

    for (int i = 0; i < initWin; i++) ///Index=maxWindow;keys[leftPeakIndex]>=leftMinValue; )
    {
      if (moveLeftKeyIndex <= 0 || leftIndex <= 0) {
        //move back to right because it cannot move to left any more
        moveLeftKeyIndex = findLightKey(keys, moveLeftKeyIndex, MOVE_RIGHT, sModel.getIFile(), conf);
        leftIndex++;

        break;
      }

        //result[leftIndex] = readSpectrum(keys, moveLeftKeyIndex, index, sModel.getDiff(), file, sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());

        result[leftIndex] = readSpectrumLabelfree(keys, moveLeftKeyIndex, index, file, sModel.getBioSample(), sModel.getYioSample(), conf, scanToRetMap);

      totalIntensity += result[leftIndex][1]; //bion sum
   //   totalIntensity += result[leftIndex][2];  //yion sum


    //  System.out.println("===========" + result[leftIndex][1] + "\t" + result[leftIndex][2]);

      moveLeftKeyIndex = findLightKey(keys, moveLeftKeyIndex, MOVE_LEFT, sModel.getIFile(), conf);

      leftIndex--;
    }

    for (int i = 0; i < initWin; i++) //leftIndex=maxWindow;keys[leftPeakIndex]>=leftMinValue; )
    {
      if (moveRightKeyIndex >= keys.length || moveRightKeyIndex < 0) {
        //	moveRightKeyIndex -= 1*numIsoWindow;
        moveRightKeyIndex = findLightKey(keys, moveRightKeyIndex, MOVE_LEFT, sModel.getIFile(), conf);

        rightIndex--;

        break;
      }

      //if(conf.isDataIndependent()) //data independent

      result[rightIndex] = readSpectrumLabelfree(keys, moveRightKeyIndex, index, file, sModel.getBioSample(), sModel.getYioSample(), conf, scanToRetMap);

      //readSpectrumLabelfree(keys, moveLeftKeyIndex, index, file, sModel.getBioSample(), sModel.getYioSample());

      totalIntensity += result[rightIndex][1];
    //  totalIntensity += result[rightIndex][2];

      //moveRightKeyIndex += numIsoWindow;
      moveRightKeyIndex = findLightKey(keys, moveRightKeyIndex, MOVE_RIGHT, sModel.getIFile(), conf);
      rightIndex++;
    }

    if (moveRightKeyIndex < 0 && moveLeftKeyIndex < 0)
      return null;

    boolean isGoingUp = true;
    boolean isHighIntensity = true;  //if the intensity is lower than one third of average intensity, this becomes false
    double[][] arr = new double[steepArea][3];
    double[][] prevArr = new double[steepArea][3];

    if (leftIndex + steepArea >= 0) {
      //for(int i=0;i<300;i++)
      while (true) {
        double area1 = 0;
        double area2 = 0;

        if (moveLeftKeyIndex - steepArea < 0) {
          break;
        }

        int steepCount = 0;
        int tempKeyIndex = moveLeftKeyIndex;
        while (true) {
          if (tempKeyIndex < 0 || steepCount < 0 || steepCount >= steepArea - 1)
            break;

          if (prevArr[0][0] != 0) {
            for (int l = 0; l < 2; l++)
              for (int m = 0; m < 3; m++)
                arr[l][m] = prevArr[l + 1][m];
          } else {

              arr[steepCount] = readSpectrumLabelfree(keys, tempKeyIndex, index, file, sModel.getBioSample(), sModel.getYioSample(), conf, scanToRetMap);

          }

          steepCount++;
          //tempKeyIndex--;
//		    tempKeyIndex = tempKeyIndex - 1*numIsoWindow;
          tempKeyIndex = findLightKey(keys, tempKeyIndex, MOVE_LEFT, sModel.getIFile(), conf);
        }

        if (tempKeyIndex < 0)
          break;

        //	if(conf.isDataIndependent()) //data independent

          arr[steepCount] = readSpectrumLabelfree(keys, tempKeyIndex, index, file, sModel.getBioSample(), sModel.getYioSample(), conf, scanToRetMap);


        for (int k = 0; k < 3; k++) {
          for (int l = 0; l < 3; l++)
            prevArr[k][l] = arr[k][l];
        }


        for (int j = leftIndex + 1; j <= leftIndex + steepCount + 1; j++) {
          area1 += result[j][1];
          area1 += result[j][2];
        }

        for (int j = 0; j <= steepCount; j++) {
          area2 += arr[j][1];
          area2 += arr[j][2];
        }

        if (area2 == 0 && area1 == 0)
          break;

        double areaRatio = (double) area2 / area1;

        /*****************************************
         * AREA1 is right side area
         * AREA2 is left side area
         *****************************************/
        if (isHighIntensity && ((double) totalIntensity / (rightIndex - leftIndex - 1) / 1.5 > area2 / 3))
          isHighIntensity = false;

        //if( (isGoingUp && (float)leftTotalIntensity/(rightIndex-leftIndex-1)/3<area2/3) || areaRatio<steepRatioThreshold )
        //if( (isGoingUp && (float)totalIntensity/(rightIndex-leftIndex-1)/3<area2/3) || areaRatio<steepRatioThreshold )
        if ((isGoingUp && isHighIntensity) || areaRatio < steepRatioThreshold) {
          if (leftIndex <= 0) {
            isGoingUp = false;
            break;
          }

          result[leftIndex] = arr[0];
          moveLeftKeyIndex = findLightKey(keys, moveLeftKeyIndex, MOVE_LEFT, sModel.getIFile(), conf);

          if (0 == result[leftIndex][0])
            continue;

          leftIndex--;

          for (int j = 0; j < steepCount; j++)
            arr[j] = arr[j + 1];

          totalIntensity += arr[0][1];
          totalIntensity += arr[0][2];

          continue;
        } else
          isGoingUp = false;


        if (areaRatio > steepRatioThreshold) {
          for (int k = 0; k < 3; k++) {
            if (result[leftIndex + 1][1] < arr[k][1] || arr[k][0] == 0)
              break;

            if (leftIndex < 0)
              break;

            result[leftIndex] = arr[k];

            //moveLeftKeyIndex -= numIsoWindow;
            moveLeftKeyIndex = findLightKey(keys, moveLeftKeyIndex, MOVE_LEFT, sModel.getIFile(), conf);
            leftIndex--;
          }

          break;
        }

        result[leftIndex--] = arr[0];
        moveLeftKeyIndex = findLightKey(keys, moveLeftKeyIndex, MOVE_LEFT, sModel.getIFile(), conf);

        for (int j = 0; j < steepCount; j++)
          arr[j] = arr[j + 1];
      }
    }

    isGoingUp = true;
    isHighIntensity = true;

    prevArr = new double[steepArea][3]; //clean up the array

    if (rightIndex + steepArea < result.length) {
      //for(int i=0;i<500;i++)
      while (true) {
        double area1 = 0;
        double area2 = 0;

        if (moveRightKeyIndex + steepArea >= keys.length)
          break;

        int steepCount = 0;

        int tempKeyIndex = moveRightKeyIndex;
        while (true) {
          if (tempKeyIndex < 0 || tempKeyIndex >= keys.length || steepCount >= steepArea - 1 || steepCount >= keys.length)
            break;

          if (prevArr[0][0] != 0) {
            for (int l = 0; l < 2; l++)
              for (int m = 0; m < 3; m++)
                arr[l][m] = prevArr[l + 1][m];
          } else {

            arr[steepCount] = readSpectrumLabelfree(keys, tempKeyIndex, index, file, sModel.getBioSample(), sModel.getYioSample(), conf, scanToRetMap);

          }
          steepCount++;
          //tempKeyIndex = tempKeyIndex + 1*numIsoWindow;
          tempKeyIndex = findLightKey(keys, tempKeyIndex, MOVE_RIGHT, sModel.getIFile(), conf);
          //tempKeyIndex++;
        }

        if (tempKeyIndex >= keys.length || tempKeyIndex < 0)
          break;

        //if(conf.isDataIndependent()) //data independent

          arr[steepCount] = readSpectrumLabelfree(keys, tempKeyIndex, index, file, sModel.getBioSample(), sModel.getYioSample(), conf, scanToRetMap);
//        else if (sModel.isHighRes())
  //        arr[steepCount] = readFullSpectrum(keys, tempKeyIndex, index, file, sModel.getSamIsoArr(), sModel.getRefIsoArr());

        for (int k = 0; k < 3; k++) {
          for (int l = 0; l < 3; l++)
            prevArr[k][l] = arr[k][l];
        }

        for (int j = rightIndex - 1; j > rightIndex - steepCount - 2; j--) {
          area1 += result[j][1];
          area1 += result[j][2];
        }

        for (int j = 0; j <= steepCount; j++) {
          area2 += arr[j][1];
          area2 += arr[j][2];
        }

        if (area2 == 0 && area1 == 0)
          break;

        double areaRatio = (double) area2 / area1;

        /*****************************************
         * AREA1 is left side area
         * AREA2 is right side area
         *****************************************/

        if (isHighIntensity && ((double) totalIntensity / (rightIndex - leftIndex - 1) / 1.5 > area2 / 3))
          isHighIntensity = false;

        //if( ((isGoingUp && (double)totalIntensity/(rightIndex-leftIndex-1)/3<area2/3)) || areaRatio<steepRatioThreshold )
        if ((isGoingUp && isHighIntensity) || areaRatio < steepRatioThreshold) {
          if ((rightIndex + 1) >= result.length) {
            isGoingUp = false;
            break;
          }

          result[rightIndex] = arr[0];
//		    moveRightKeyIndex += numIsoWindow;
          moveRightKeyIndex = findLightKey(keys, moveRightKeyIndex, MOVE_RIGHT, sModel.getIFile(), conf);

          if (0 == result[rightIndex][0])
            continue;

          rightIndex++;

          for (int j = 0; j < steepCount; j++)
            arr[j] = arr[j + 1];

          //rightTotalIntensity += arr[0][1];
          totalIntensity += arr[0][1];
          totalIntensity += arr[0][2];

          continue;
        } else
          isGoingUp = false;

        if (area2 / area1 > steepRatioThreshold) {
          for (int k = 0; k < 3; k++) {
            if (result[rightIndex - 1][1] < arr[k][1])
              break;

            if (result.length <= rightIndex)
              break;

            result[rightIndex] = arr[k];

            //System.out.println(arr[0][1] + " " + arr[1][1] + " " + arr[2][1] + " " + result[rightIndex-1][1]);

            //moveRightKeyIndex += numIsoWindow;
            moveRightKeyIndex = findLightKey(keys, moveRightKeyIndex, MOVE_RIGHT, sModel.getIFile(), conf);
            rightIndex++;
          }
          break;
        }

        //System.out.println("2-->>peak gogo" + " " + rightIndex + " " + moveRightKeyIndex + " " + index + " " + sModel.getDiff());
        result[rightIndex++] = arr[0];

        //moveRightKeyIndex += numIsoWindow;
        moveRightKeyIndex = findLightKey(keys, moveRightKeyIndex, MOVE_RIGHT, sModel.getIFile(), conf);

        for (int j = 0; j < steepCount; j++)
          arr[j] = arr[j + 1];
      }
    }
    int peakStart = leftIndex + 1;
    int peakEnd = rightIndex - 1;

    int dtaStart = range.getMin();
    int dtaEnd = range.getMax();


    int specSpace = (int) (result[peakStart + 1][0] - result[peakStart][0]);

    if (specSpace == 0)
      specSpace = (int) (result[peakStart][0] - result[peakStart - 1][0]);

    int diff = (int) (result[peakStart][0] - dtaStart) / specSpace;
    int leftMargin = margin;

    int rightMargin = margin;
    if (diff > 0)
      leftMargin += diff;

    diff = (int) (dtaEnd - result[rightIndex - 1][0]) / specSpace;

    if (diff > 0)
      rightMargin += diff;


    for (int i = 0; i < leftMargin; i++) {
      if (leftIndex < 0 || moveLeftKeyIndex < 0)
        break;

      //if(conf.isDataIndependent()) //data independent

        result[leftIndex] = readSpectrumLabelfree(keys, moveLeftKeyIndex, index, file, sModel.getBioSample(), sModel.getYioSample(), conf, scanToRetMap);

      if (0 == result[leftIndex][0]) {
        //moveLeftKeyIndex -= numIsoWindow;
        moveLeftKeyIndex = findLightKey(keys, moveLeftKeyIndex, MOVE_LEFT, sModel.getIFile(), conf);
        continue;
      }


      //moveLeftKeyIndex -= numIsoWindow;
      moveLeftKeyIndex = findLightKey(keys, moveLeftKeyIndex, MOVE_LEFT, sModel.getIFile(), conf);

      if (moveLeftKeyIndex < 0)
        break;

      leftIndex--;
    }
    for (int i = 0; i < rightMargin; i++) {
      if (rightIndex >= result.length || moveRightKeyIndex >= keys.length || moveRightKeyIndex < 0)
        break;

      //if(conf.isDataIndependent()) //data independent
        result[rightIndex] = readSpectrumLabelfree(keys, moveRightKeyIndex, index, file, sModel.getBioSample(), sModel.getYioSample(), conf, scanToRetMap);


      if (0 == result[rightIndex][0]) {
        //moveRightKeyIndex += numIsoWindow;
        moveRightKeyIndex = findLightKey(keys, moveRightKeyIndex, MOVE_RIGHT, sModel.getIFile(), conf);
        continue;
      }


      //moveRightKeyIndex += numIsoWindow;
      moveRightKeyIndex = findLightKey(keys, moveRightKeyIndex, MOVE_RIGHT, sModel.getIFile(), conf);
      if (moveRightKeyIndex < 0)
        break;

      rightIndex++;
    }


    //if(conf.isDataIndependent())
 //   if (conf.getQuantLevel() == 2)
      return buildDIResultLabelfree(peakStart, peakEnd, leftIndex, rightIndex, result, sModel.getBioSample().length, null);
   // else
   //   return buildDDResult(peakStart, peakEnd, leftIndex, rightIndex, result);

    //return buildDDResult(peakStart, peakEnd, leftIndex, rightIndex, result);

  }

  //MSMS peak finding.  labeled. data independent
  public static String peakFindingMSMS(SpectrumModel sModel, SpecRange range, Configuration conf, int keyIndex)
    throws PrecursorNotFoundException, IOException, CensusIndexOutOfBoundException, Exception {
    TIntLongHashMap index = sModel.getIndex();
    RandomAccessFile file = sModel.getFile();
    steepRatioThreshold = conf.getSteepRatioThreshold();

//	int numIsoWindow = conf.getNumOfIsolationWindow();
    int maxWindow = conf.getMaxWindow();
    int margin = conf.getMargin();

    double[][] result = null;
    //conf.getQuantLevel()

    if (conf.getQuantLevel() == 2)
      result = new double[maxWindow * 2 + 1 + margin * 2][4 * sModel.getBioSample().length + 3]; //scan #, sample intensity, ref intensity
    else
      result = new double[maxWindow * 2 + 1 + margin * 2][3];

    int leftIndex = maxWindow + margin;
    int rightIndex = maxWindow + margin + 1;

    //double leftTotalIntensity=0;
    //double rightTotalIntensity=0;
    double totalIntensity = 0;
    int[] keys = sModel.getKeys();
    int initWin = 2;

    int steepArea = conf.getSteepArea();

    int moveLeftKeyIndex = keyIndex;
    int moveRightKeyIndex = findLightKey(keys, keyIndex, MOVE_RIGHT, sModel.getIFile(), conf);

    IndexedFile iFile = sModel.getIFile();

    for (int i = 0; i < initWin; i++) ///Index=maxWindow;keys[leftPeakIndex]>=leftMinValue; )
    {
      if (moveLeftKeyIndex <= 0 || leftIndex <= 0) {
        //move back to right because it cannot move to left any more
        moveLeftKeyIndex = findLightKey(keys, moveLeftKeyIndex, MOVE_RIGHT, sModel.getIFile(), conf);
        leftIndex++;

        break;
      }
/////
      if (conf.getQuantLevel() == 2) {

        //result[leftIndex] = readSpectrum(keys, moveLeftKeyIndex, index, sModel.getDiff(), file, sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());
        result[leftIndex] = readSpectrum(keys, moveLeftKeyIndex, index, findHeavyKey(keys, moveLeftKeyIndex, iFile, conf, conf.getCalcRefAvgMass()), file, sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());
      } else //data dependent
      {
        if (sModel.isHighRes())
          result[leftIndex] = readFullSpectrum(keys, moveLeftKeyIndex, index, file, sModel.getSamIsoArr(), sModel.getRefIsoArr());
        else
          result[leftIndex] = readFullSpectrum(keys, moveLeftKeyIndex, index, file, sModel.getSamStartMass(), sModel.getSamEndMass(), sModel.getRefStartMass(), sModel.getRefEndMass());

      }
      totalIntensity += result[leftIndex][1];
      totalIntensity += result[leftIndex][2];

      moveLeftKeyIndex = findLightKey(keys, moveLeftKeyIndex, MOVE_LEFT, sModel.getIFile(), conf);

      leftIndex--;
    }

    for (int i = 0; i < initWin; i++) //leftIndex=maxWindow;keys[leftPeakIndex]>=leftMinValue; )
    {
      if (moveRightKeyIndex >= keys.length || moveRightKeyIndex < 0) {
        //	moveRightKeyIndex -= 1*numIsoWindow;
        moveRightKeyIndex = findLightKey(keys, moveRightKeyIndex, MOVE_LEFT, sModel.getIFile(), conf);

        rightIndex--;

        break;
      }

      //if(conf.isDataIndependent()) //data independent
      if (conf.getQuantLevel() == 2) {
        result[rightIndex] = readSpectrum(keys, moveRightKeyIndex, index, findHeavyKey(keys, moveRightKeyIndex, iFile, conf, conf.getCalcRefAvgMass()), file, sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());
      } else //data dependent
      {
        if (sModel.isHighRes())
          result[rightIndex] = readFullSpectrum(keys, moveRightKeyIndex, index, file, sModel.getSamIsoArr(), sModel.getRefIsoArr());
        else
          result[rightIndex] = readFullSpectrum(keys, moveRightKeyIndex, index, file, sModel.getSamStartMass(), sModel.getSamEndMass(), sModel.getRefStartMass(), sModel.getRefEndMass());
      }

      totalIntensity += result[rightIndex][1];
      totalIntensity += result[rightIndex][2];

      //moveRightKeyIndex += numIsoWindow;
      moveRightKeyIndex = findLightKey(keys, moveRightKeyIndex, MOVE_RIGHT, sModel.getIFile(), conf);
      rightIndex++;
    }

    if (moveRightKeyIndex < 0 && moveLeftKeyIndex < 0)
      return null;

    boolean isGoingUp = true;
    boolean isHighIntensity = true;  //if the intensity is lower than one third of average intensity, this becomes false
    double[][] arr = new double[steepArea][3];
    double[][] prevArr = new double[steepArea][3];

    if (leftIndex + steepArea >= 0) {
      //for(int i=0;i<300;i++)
      while (true) {
        double area1 = 0;
        double area2 = 0;

        if (moveLeftKeyIndex - steepArea < 0) {
          break;
        }

        int steepCount = 0;
        int tempKeyIndex = moveLeftKeyIndex;
        while (true) {
          if (tempKeyIndex < 0 || steepCount < 0 || steepCount >= steepArea - 1)
            break;

          if (prevArr[0][0] != 0) {
            for (int l = 0; l < 2; l++)
              for (int m = 0; m < 3; m++)
                arr[l][m] = prevArr[l + 1][m];
          } else {
            if (conf.getQuantLevel() == 2)
              arr[steepCount] = readSpectrum(keys, tempKeyIndex, index, findHeavyKey(keys, tempKeyIndex, iFile, conf, conf.getCalcRefAvgMass()), file, sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());
            else {
              if (sModel.isHighRes())
                arr[steepCount] = readFullSpectrum(keys, tempKeyIndex, index, file, sModel.getSamIsoArr(), sModel.getRefIsoArr());
              else
                arr[steepCount] = readFullSpectrum(keys, tempKeyIndex, index, file, sModel.getSamStartMass(), sModel.getSamEndMass(), sModel.getRefStartMass(), sModel.getRefEndMass());
            }
          }

          steepCount++;
          //tempKeyIndex--;
//		    tempKeyIndex = tempKeyIndex - 1*numIsoWindow;
          tempKeyIndex = findLightKey(keys, tempKeyIndex, MOVE_LEFT, sModel.getIFile(), conf);
        }

        if (tempKeyIndex < 0)
          break;

        //	if(conf.isDataIndependent()) //data independent
        if (conf.getQuantLevel() == 2)
          arr[steepCount] = readSpectrum(keys, tempKeyIndex, index, findHeavyKey(keys, tempKeyIndex, iFile, conf, conf.getCalcRefAvgMass()), file, sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());
        else {
          if (sModel.isHighRes())
            arr[steepCount] = readFullSpectrum(keys, tempKeyIndex, index, file, sModel.getSamIsoArr(), sModel.getRefIsoArr());
          else
            arr[steepCount] = readFullSpectrum(keys, tempKeyIndex, index, file, sModel.getSamStartMass(), sModel.getSamEndMass(), sModel.getRefStartMass(), sModel.getRefEndMass());
        }

        for (int k = 0; k < 3; k++) {
          for (int l = 0; l < 3; l++)
            prevArr[k][l] = arr[k][l];
        }


        for (int j = leftIndex + 1; j <= leftIndex + steepCount + 1; j++) {
          area1 += result[j][1];
          area1 += result[j][2];
        }

        for (int j = 0; j <= steepCount; j++) {
          area2 += arr[j][1];
          area2 += arr[j][2];
        }

        if (area2 == 0 && area1 == 0)
          break;

        double areaRatio = (double) area2 / area1;

        /*****************************************
         * AREA1 is right side area
         * AREA2 is left side area
         *****************************************/
        if (isHighIntensity && ((double) totalIntensity / (rightIndex - leftIndex - 1) / 1.5 > area2 / 3))
          isHighIntensity = false;

        //if( (isGoingUp && (float)leftTotalIntensity/(rightIndex-leftIndex-1)/3<area2/3) || areaRatio<steepRatioThreshold )
        //if( (isGoingUp && (float)totalIntensity/(rightIndex-leftIndex-1)/3<area2/3) || areaRatio<steepRatioThreshold )
        if ((isGoingUp && isHighIntensity) || areaRatio < steepRatioThreshold) {
          if (leftIndex <= 0) {
            isGoingUp = false;
            break;
          }

          result[leftIndex] = arr[0];
          moveLeftKeyIndex = findLightKey(keys, moveLeftKeyIndex, MOVE_LEFT, sModel.getIFile(), conf);

          if (0 == result[leftIndex][0])
            continue;

          leftIndex--;

          for (int j = 0; j < steepCount; j++)
            arr[j] = arr[j + 1];

          totalIntensity += arr[0][1];
          totalIntensity += arr[0][2];

          continue;
        } else
          isGoingUp = false;


        if (areaRatio > steepRatioThreshold) {
          for (int k = 0; k < 3; k++) {
            if (result[leftIndex + 1][1] < arr[k][1] || arr[k][0] == 0)
              break;

            if (leftIndex < 0)
              break;

            result[leftIndex] = arr[k];

            //moveLeftKeyIndex -= numIsoWindow;
            moveLeftKeyIndex = findLightKey(keys, moveLeftKeyIndex, MOVE_LEFT, sModel.getIFile(), conf);
            leftIndex--;
          }

          break;
        }

        result[leftIndex--] = arr[0];
        moveLeftKeyIndex = findLightKey(keys, moveLeftKeyIndex, MOVE_LEFT, sModel.getIFile(), conf);

        for (int j = 0; j < steepCount; j++)
          arr[j] = arr[j + 1];
      }
    }

    isGoingUp = true;
    isHighIntensity = true;

    prevArr = new double[steepArea][3]; //clean up the array

    if (rightIndex + steepArea < result.length) {
      //for(int i=0;i<500;i++)
      while (true) {
        double area1 = 0;
        double area2 = 0;

        if (moveRightKeyIndex + steepArea >= keys.length)
          break;

        int steepCount = 0;

        int tempKeyIndex = moveRightKeyIndex;
        while (true) {
          if (tempKeyIndex < 0 || tempKeyIndex >= keys.length || steepCount >= steepArea - 1 || steepCount >= keys.length)
            break;

          if (prevArr[0][0] != 0) {
            for (int l = 0; l < 2; l++)
              for (int m = 0; m < 3; m++)
                arr[l][m] = prevArr[l + 1][m];
          } else {
            if (conf.getQuantLevel() == 2)
              //if(conf.isDataIndependent()) //data independent
              arr[steepCount] = readSpectrum(keys, tempKeyIndex, index, findHeavyKey(keys, tempKeyIndex, iFile, conf, conf.getCalcRefAvgMass()), file, sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());
            else if (sModel.isHighRes())
              arr[steepCount] = readFullSpectrum(keys, tempKeyIndex, index, file, sModel.getSamIsoArr(), sModel.getRefIsoArr());
            else
              arr[steepCount] = readFullSpectrum(keys, tempKeyIndex, index, file, sModel.getSamStartMass(), sModel.getSamEndMass(), sModel.getRefStartMass(), sModel.getRefEndMass());
          }

          steepCount++;
          //tempKeyIndex = tempKeyIndex + 1*numIsoWindow;
          tempKeyIndex = findLightKey(keys, tempKeyIndex, MOVE_RIGHT, sModel.getIFile(), conf);
          //tempKeyIndex++;
        }

        if (tempKeyIndex >= keys.length || tempKeyIndex < 0)
          break;

        //if(conf.isDataIndependent()) //data independent
        if (conf.getQuantLevel() == 2)
          arr[steepCount] = readSpectrum(keys, tempKeyIndex, index, findHeavyKey(keys, tempKeyIndex, iFile, conf, conf.getCalcRefAvgMass()), file, sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());
        else if (sModel.isHighRes())
          arr[steepCount] = readFullSpectrum(keys, tempKeyIndex, index, file, sModel.getSamIsoArr(), sModel.getRefIsoArr());
        else
          arr[steepCount] = readFullSpectrum(keys, tempKeyIndex, index, file, sModel.getSamStartMass(), sModel.getSamEndMass(), sModel.getRefStartMass(), sModel.getRefEndMass());

        for (int k = 0; k < 3; k++) {
          for (int l = 0; l < 3; l++)
            prevArr[k][l] = arr[k][l];
        }

        for (int j = rightIndex - 1; j > rightIndex - steepCount - 2; j--) {
          area1 += result[j][1];
          area1 += result[j][2];
        }

        for (int j = 0; j <= steepCount; j++) {
          area2 += arr[j][1];
          area2 += arr[j][2];
        }

        if (area2 == 0 && area1 == 0)
          break;

        double areaRatio = (double) area2 / area1;

        /*****************************************
         * AREA1 is left side area
         * AREA2 is right side area
         *****************************************/

        if (isHighIntensity && ((double) totalIntensity / (rightIndex - leftIndex - 1) / 1.5 > area2 / 3))
          isHighIntensity = false;

        //if( ((isGoingUp && (double)totalIntensity/(rightIndex-leftIndex-1)/3<area2/3)) || areaRatio<steepRatioThreshold )
        if ((isGoingUp && isHighIntensity) || areaRatio < steepRatioThreshold) {
          if ((rightIndex + 1) >= result.length) {
            isGoingUp = false;
            break;
          }

          result[rightIndex] = arr[0];
//		    moveRightKeyIndex += numIsoWindow;
          moveRightKeyIndex = findLightKey(keys, moveRightKeyIndex, MOVE_RIGHT, sModel.getIFile(), conf);

          if (0 == result[rightIndex][0])
            continue;

          rightIndex++;

          for (int j = 0; j < steepCount; j++)
            arr[j] = arr[j + 1];

          //rightTotalIntensity += arr[0][1];
          totalIntensity += arr[0][1];
          totalIntensity += arr[0][2];

          continue;
        } else
          isGoingUp = false;

        if (area2 / area1 > steepRatioThreshold) {
          for (int k = 0; k < 3; k++) {
            if (result[rightIndex - 1][1] < arr[k][1])
              break;

            if (result.length <= rightIndex)
              break;

            result[rightIndex] = arr[k];

            //System.out.println(arr[0][1] + " " + arr[1][1] + " " + arr[2][1] + " " + result[rightIndex-1][1]);

            //moveRightKeyIndex += numIsoWindow;
            moveRightKeyIndex = findLightKey(keys, moveRightKeyIndex, MOVE_RIGHT, sModel.getIFile(), conf);
            rightIndex++;
          }
          break;
        }

        //System.out.println("2-->>peak gogo" + " " + rightIndex + " " + moveRightKeyIndex + " " + index + " " + sModel.getDiff());
        result[rightIndex++] = arr[0];

        //moveRightKeyIndex += numIsoWindow;
        moveRightKeyIndex = findLightKey(keys, moveRightKeyIndex, MOVE_RIGHT, sModel.getIFile(), conf);

        for (int j = 0; j < steepCount; j++)
          arr[j] = arr[j + 1];
      }
    }
    int peakStart = leftIndex + 1;
    int peakEnd = rightIndex - 1;

    int dtaStart = range.getMin();
    int dtaEnd = range.getMax();


    int specSpace = (int) (result[peakStart + 1][0] - result[peakStart][0]);

    if (specSpace == 0)
      specSpace = (int) (result[peakStart][0] - result[peakStart - 1][0]);

    int diff = (int) (result[peakStart][0] - dtaStart) / specSpace;
    int leftMargin = margin;

    int rightMargin = margin;
    if (diff > 0)
      leftMargin += diff;

    diff = (int) (dtaEnd - result[rightIndex - 1][0]) / specSpace;

    if (diff > 0)
      rightMargin += diff;


    for (int i = 0; i < leftMargin; i++) {
      if (leftIndex < 0 || moveLeftKeyIndex < 0)
        break;

      //if(conf.isDataIndependent()) //data independent
      if (conf.getQuantLevel() == 2)
        result[leftIndex] = readSpectrum(keys, moveLeftKeyIndex, index, findHeavyKey(keys, moveLeftKeyIndex, iFile, conf, conf.getCalcRefAvgMass()), file, sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());
      else if (sModel.isHighRes())
        result[leftIndex] = readFullSpectrum(keys, moveLeftKeyIndex, index, file, sModel.getSamIsoArr(), sModel.getRefIsoArr());
      else
        result[leftIndex] = readFullSpectrum(keys, moveLeftKeyIndex, index, file, sModel.getSamStartMass(), sModel.getSamEndMass(), sModel.getRefStartMass(), sModel.getRefEndMass());

      if (0 == result[leftIndex][0]) {
        //moveLeftKeyIndex -= numIsoWindow;
        moveLeftKeyIndex = findLightKey(keys, moveLeftKeyIndex, MOVE_LEFT, sModel.getIFile(), conf);
        continue;
      }


      //moveLeftKeyIndex -= numIsoWindow;
      moveLeftKeyIndex = findLightKey(keys, moveLeftKeyIndex, MOVE_LEFT, sModel.getIFile(), conf);

      if (moveLeftKeyIndex < 0)
        break;

      leftIndex--;
    }
    for (int i = 0; i < rightMargin; i++) {
      if (rightIndex >= result.length || moveRightKeyIndex >= keys.length || moveRightKeyIndex < 0)
        break;

      //if(conf.isDataIndependent()) //data independent
      if (conf.getQuantLevel() == 2)
        result[rightIndex] = readSpectrum(keys, moveRightKeyIndex, index, findHeavyKey(keys, moveRightKeyIndex, iFile, conf, conf.getCalcRefAvgMass()), file, sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());
      else if (sModel.isHighRes())
        result[rightIndex] = readFullSpectrum(keys, moveRightKeyIndex, index, file, sModel.getSamIsoArr(), sModel.getRefIsoArr());
      else
        result[rightIndex] = readFullSpectrum(keys, moveRightKeyIndex, index, file, sModel.getSamStartMass(), sModel.getSamEndMass(), sModel.getRefStartMass(), sModel.getRefEndMass());

      if (0 == result[rightIndex][0]) {
        //moveRightKeyIndex += numIsoWindow;
        moveRightKeyIndex = findLightKey(keys, moveRightKeyIndex, MOVE_RIGHT, sModel.getIFile(), conf);
        continue;
      }


      //moveRightKeyIndex += numIsoWindow;
      moveRightKeyIndex = findLightKey(keys, moveRightKeyIndex, MOVE_RIGHT, sModel.getIFile(), conf);
      if (moveRightKeyIndex < 0)
        break;

      rightIndex++;
    }


    //if(conf.isDataIndependent())
    if (conf.getQuantLevel() == 2)
      return buildDIResult(peakStart, peakEnd, leftIndex, rightIndex, result, sModel.getBioSample().length, null);
    else
      return buildDDResult(peakStart, peakEnd, leftIndex, rightIndex, result);

  }

  //peak finding for labeling + high + low res + labelfree
  public static String peakFinding(SpectrumModel sModel, SpecRange range, int keyIndex, DisplayData.DisplayChroData chroData)
    throws IOException, CensusIndexOutOfBoundException, Exception {

    TIntLongHashMap index = sModel.getIndex();
    Object file = sModel.getGenericIndexFile(conf);

    steepRatioThreshold = conf.getSteepRatioThreshold();

    int numIsoWindow = conf.getNumOfIsolationWindow();
    int maxWindow = conf.getMaxWindow();
    int margin = conf.getMargin();

    double[][] result = null;
    //conf.getQuantLevel()

    if (conf.isDataIndependent())
      result = new double[maxWindow * 2 + 1 + margin * 2][4 * sModel.getBioSample().length + 3]; //scan #, sample intensity, ref intensity
    else
      result = new double[maxWindow * 2 + 1 + margin * 2][3];

    int leftIndex = maxWindow + margin;
    int rightIndex = maxWindow + margin + 1;


    //double leftTotalIntensity=0;
    //double rightTotalIntensity=0;
    double totalIntensity = 0;

    int steepArea = conf.getSteepArea();

    int moveLeftKeyIndex = keyIndex;

    int moveRightKeyIndex = keyIndex + 1 * numIsoWindow;

    int[] keys = sModel.getKeys();
    int initWin = 2;

    for (int i = 0; i < initWin; i++) ///Index=maxWindow;keys[leftPeakIndex]>=leftMinValue; )
    {
      if (moveLeftKeyIndex <= 0 || leftIndex <= 0) {
        moveLeftKeyIndex += 1 * numIsoWindow;
        leftIndex++;

        break;
      }

      if (conf.isDataIndependent()) //data independent
      {
        result[leftIndex] = readSpectrum(keys, moveLeftKeyIndex, index, sModel.getDiff(), file, sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());

      } else //data dependent
      {
        if (sModel.isHighRes())
          result[leftIndex] = readFullSpectrum(keys, moveLeftKeyIndex, index, file, sModel.getSamIsoArr(), sModel.getRefIsoArr());
        else
          result[leftIndex] = readFullSpectrum(keys, moveLeftKeyIndex, index, file, sModel.getSamStartMass(), sModel.getSamEndMass(), sModel.getRefStartMass(), sModel.getRefEndMass());

      }
      totalIntensity += result[leftIndex][1];
      totalIntensity += result[leftIndex][2];

      // if(moveLeftKeyIndex<=0 || leftIndex<=0)
      //	break;

      moveLeftKeyIndex -= numIsoWindow;
      leftIndex--;
    }

    for (int i = 0; i < initWin; i++) //leftIndex=maxWindow;keys[leftPeakIndex]>=leftMinValue; )
    {
      if (moveRightKeyIndex >= keys.length) {
        moveRightKeyIndex -= 1 * numIsoWindow;
        rightIndex--;

        break;
      }

      if (conf.isDataIndependent()) //data independent
      {
        result[rightIndex] = readSpectrum(keys, moveRightKeyIndex, index, sModel.getDiff(), file, sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());
      } else //data dependent
      {
        if (sModel.isHighRes())
          result[rightIndex] = readFullSpectrum(keys, moveRightKeyIndex, index, file, sModel.getSamIsoArr(), sModel.getRefIsoArr());
        else
          result[rightIndex] = readFullSpectrum(keys, moveRightKeyIndex, index, file, sModel.getSamStartMass(), sModel.getSamEndMass(), sModel.getRefStartMass(), sModel.getRefEndMass());
      }

      totalIntensity += result[rightIndex][1];
      totalIntensity += result[rightIndex][2];

      moveRightKeyIndex += numIsoWindow;
      rightIndex++;
    }

    boolean isGoingUp = true;
    boolean isHighIntensity = true;  //if the intensity is lower than one third of average intensity, this becomes false
    double[][] arr = new double[steepArea][3];
    double[][] prevArr = new double[steepArea][3];

    if (leftIndex + steepArea >= 0) {
      //for(int i=0;i<300;i++)
      while (true) {
        double area1 = 0;
        double area2 = 0;

        if (moveLeftKeyIndex - steepArea < 0) {
          break;
        }

        int steepCount = 0;
        int tempKeyIndex = moveLeftKeyIndex;
        while (true) {
          if (tempKeyIndex < 0 || steepCount < 0 || steepCount >= steepArea - 1)
            break;

          if (prevArr[0][0] != 0) {
            for (int l = 0; l < 2; l++)
              for (int m = 0; m < 3; m++)
                arr[l][m] = prevArr[l + 1][m];
          } else {
            if (conf.isDataIndependent()) //data independent
              arr[steepCount] = readSpectrum(keys, tempKeyIndex, index, sModel.getDiff(), file, sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());
            else {
              if (sModel.isHighRes())
                arr[steepCount] = readFullSpectrum(keys, tempKeyIndex, index, file, sModel.getSamIsoArr(), sModel.getRefIsoArr());
              else
                arr[steepCount] = readFullSpectrum(keys, tempKeyIndex, index, file, sModel.getSamStartMass(), sModel.getSamEndMass(), sModel.getRefStartMass(), sModel.getRefEndMass());
            }
          }

          steepCount++;
          //tempKeyIndex--;
          tempKeyIndex = tempKeyIndex - 1 * numIsoWindow;
        }

        if (tempKeyIndex < 0)
          break;

        if (conf.isDataIndependent()) //data independent
          arr[steepCount] = readSpectrum(keys, tempKeyIndex, index, sModel.getDiff(), file, sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());
        else {
          if (sModel.isHighRes())
            arr[steepCount] = readFullSpectrum(keys, tempKeyIndex, index, file, sModel.getSamIsoArr(), sModel.getRefIsoArr());
          else
            arr[steepCount] = readFullSpectrum(keys, tempKeyIndex, index, file, sModel.getSamStartMass(), sModel.getSamEndMass(), sModel.getRefStartMass(), sModel.getRefEndMass());
        }

        for (int k = 0; k < 3; k++) {
          for (int l = 0; l < 3; l++)
            prevArr[k][l] = arr[k][l];
        }


        for (int j = leftIndex + 1; j <= leftIndex + steepCount + 1; j++) {
          area1 += result[j][1];
          area1 += result[j][2];
        }

        for (int j = 0; j <= steepCount; j++) {
          area2 += arr[j][1];
          area2 += arr[j][2];
        }

        if (area2 == 0 && area1 == 0)
          break;

        double areaRatio = (double) area2 / area1;

        /*****************************************
         * AREA1 is right side area
         * AREA2 is left side area
         *****************************************/
        if (isHighIntensity && ((double) totalIntensity / (rightIndex - leftIndex - 1) / 1.5 > area2 / 3))
          isHighIntensity = false;

        //if( (isGoingUp && (float)leftTotalIntensity/(rightIndex-leftIndex-1)/3<area2/3) || areaRatio<steepRatioThreshold )
        //if( (isGoingUp && (float)totalIntensity/(rightIndex-leftIndex-1)/3<area2/3) || areaRatio<steepRatioThreshold )
        if ((isGoingUp && isHighIntensity) || areaRatio < steepRatioThreshold) {
          if (leftIndex <= 0) {
            isGoingUp = false;
            break;
          }

          //System.out.println("peak gogo" + " " + leftIndex + " " + moveLeftKeyIndex + " " + index + " " + sModel.getDiff());
          result[leftIndex] = arr[0];
          moveLeftKeyIndex -= numIsoWindow;

          if (0 == result[leftIndex][0])
            continue;

          leftIndex--;

          for (int j = 0; j < steepCount; j++)
            arr[j] = arr[j + 1];

          totalIntensity += arr[0][1];
          totalIntensity += arr[0][2];

          continue;
        } else
          isGoingUp = false;


        if (areaRatio > steepRatioThreshold) {
          for (int k = 0; k < 3; k++) {
            if (result[leftIndex + 1][1] < arr[k][1] || arr[k][0] == 0)
              break;

            if (leftIndex < 0)
              break;

            result[leftIndex] = arr[k];

            moveLeftKeyIndex -= numIsoWindow;
            leftIndex--;
          }

          break;
        }

        result[leftIndex--] = arr[0];

        moveLeftKeyIndex -= numIsoWindow;

        for (int j = 0; j < steepCount; j++)
          arr[j] = arr[j + 1];
      }
    }

    prevArr = new double[steepArea][3]; //clean up the array

    isGoingUp = true;
    isHighIntensity = true;
    if (rightIndex + steepArea < result.length) {
      //for(int i=0;i<500;i++)
      while (true) {
        double area1 = 0;
        double area2 = 0;

        if (moveRightKeyIndex + steepArea >= keys.length)
          break;

        int steepCount = 0;

        int tempKeyIndex = moveRightKeyIndex;
        while (true) {
          if (tempKeyIndex >= keys.length || steepCount >= steepArea - 1 || steepCount >= keys.length)
            break;

          if (prevArr[0][0] != 0) {
            for (int l = 0; l < 2; l++)
              for (int m = 0; m < 3; m++)
                arr[l][m] = prevArr[l + 1][m];
          } else {
            if (conf.isDataIndependent()) //data independent
              arr[steepCount] = readSpectrum(keys, tempKeyIndex, index, sModel.getDiff(), file, sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());
            else if (sModel.isHighRes())
              arr[steepCount] = readFullSpectrum(keys, tempKeyIndex, index, file, sModel.getSamIsoArr(), sModel.getRefIsoArr());
            else
              arr[steepCount] = readFullSpectrum(keys, tempKeyIndex, index, file, sModel.getSamStartMass(), sModel.getSamEndMass(), sModel.getRefStartMass(), sModel.getRefEndMass());
          }

          steepCount++;
          tempKeyIndex = tempKeyIndex + 1 * numIsoWindow;
          //tempKeyIndex++;
        }

        if (tempKeyIndex >= keys.length)
          break;

        if (conf.isDataIndependent()) //data independent
          arr[steepCount] = readSpectrum(keys, tempKeyIndex, index, sModel.getDiff(), file, sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());
        else if (sModel.isHighRes())
          arr[steepCount] = readFullSpectrum(keys, tempKeyIndex, index, file, sModel.getSamIsoArr(), sModel.getRefIsoArr());
        else
          arr[steepCount] = readFullSpectrum(keys, tempKeyIndex, index, file, sModel.getSamStartMass(), sModel.getSamEndMass(), sModel.getRefStartMass(), sModel.getRefEndMass());

        for (int k = 0; k < 3; k++) {
          for (int l = 0; l < 3; l++)
            prevArr[k][l] = arr[k][l];
        }

        for (int j = rightIndex - 1; j > rightIndex - steepCount - 2; j--) {
          area1 += result[j][1];
          area1 += result[j][2];
        }

        for (int j = 0; j <= steepCount; j++) {
          area2 += arr[j][1];
          area2 += arr[j][2];
        }

        if (area2 == 0 && area1 == 0)
          break;

        double areaRatio = (double) area2 / area1;

        /*****************************************
         * AREA1 is left side area
         * AREA2 is right side area
         *****************************************/

        if (isHighIntensity && ((double) totalIntensity / (rightIndex - leftIndex - 1) / 1.5 > area2 / 3))
          isHighIntensity = false;


        //if( ((isGoingUp && (double)totalIntensity/(rightIndex-leftIndex-1)/3<area2/3)) || areaRatio<steepRatioThreshold )
        if ((isGoingUp && isHighIntensity) || areaRatio < steepRatioThreshold) {
          if ((rightIndex + 1) >= result.length) {
            isGoingUp = false;
            break;
          }

          result[rightIndex] = arr[0];
          moveRightKeyIndex += numIsoWindow;

          if (0 == result[rightIndex][0])
            continue;

          rightIndex++;

          for (int j = 0; j < steepCount; j++)
            arr[j] = arr[j + 1];

          //rightTotalIntensity += arr[0][1];
          totalIntensity += arr[0][1];
          totalIntensity += arr[0][2];

          continue;
        } else
          isGoingUp = false;

        if (area2 / area1 > steepRatioThreshold) {
          for (int k = 0; k < 3; k++) {
            if (result[rightIndex - 1][1] < arr[k][1])
              break;

            if (result.length <= rightIndex)
              break;

            result[rightIndex] = arr[k];

            //System.out.println(arr[0][1] + " " + arr[1][1] + " " + arr[2][1] + " " + result[rightIndex-1][1]);

            moveRightKeyIndex += numIsoWindow;
            rightIndex++;
          }
          break;
        }

        result[rightIndex++] = arr[0];

        moveRightKeyIndex += numIsoWindow;

        for (int j = 0; j < steepCount; j++)
          arr[j] = arr[j + 1];
      }
    }
    int peakStart = leftIndex + 1;
    int peakEnd = rightIndex - 1;

    int dtaStart = (range == null ? 0 : range.getMin());
    int dtaEnd = (range == null ? 0 : range.getMax());
    int specSpace = (int) (result[peakStart + 1][0] - result[peakStart][0]);


//System.out.println( result[peakStart+1][0] + " " + result[peakStart][0] );
    int diff = 0;
    if (specSpace > 0)
      diff = (int) (result[peakStart][0] - dtaStart) / specSpace;

    int leftMargin = margin;

    int rightMargin = margin;
    if (diff > 0)
      leftMargin += diff;

    if (specSpace > 0)
      diff = (int) (dtaEnd - result[rightIndex - 1][0]) / specSpace;

    if (diff > 0)
      rightMargin += diff;

    for (int i = 0; i < leftMargin; i++) {
      if (leftIndex < 0 || moveLeftKeyIndex < 0)
        break;

      if (conf.isDataIndependent()) //data independent
        result[leftIndex] = readSpectrum(keys, moveLeftKeyIndex, index, sModel.getDiff(), file, sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());
      else if (sModel.isHighRes())
        result[leftIndex] = readFullSpectrum(keys, moveLeftKeyIndex, index, file, sModel.getSamIsoArr(), sModel.getRefIsoArr());
      else
        result[leftIndex] = readFullSpectrum(keys, moveLeftKeyIndex, index, file, sModel.getSamStartMass(), sModel.getSamEndMass(), sModel.getRefStartMass(), sModel.getRefEndMass());

      if (0 == result[leftIndex][0]) {
        moveLeftKeyIndex -= numIsoWindow;
        continue;
      }

      leftIndex--;
      moveLeftKeyIndex -= numIsoWindow;
    }
    for (int i = 0; i < rightMargin; i++) {
      if (rightIndex >= result.length || moveRightKeyIndex >= keys.length)
        break;

      if (conf.isDataIndependent()) //data independent
        result[rightIndex] = readSpectrum(keys, moveRightKeyIndex, index, sModel.getDiff(), file, sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());
      else if (sModel.isHighRes())
        result[rightIndex] = readFullSpectrum(keys, moveRightKeyIndex, index, file, sModel.getSamIsoArr(), sModel.getRefIsoArr());
      else
        result[rightIndex] = readFullSpectrum(keys, moveRightKeyIndex, index, file, sModel.getSamStartMass(), sModel.getSamEndMass(), sModel.getRefStartMass(), sModel.getRefEndMass());

      if (0 == result[rightIndex][0]) {
        moveRightKeyIndex += numIsoWindow;
        continue;
      }

      rightIndex++;
      moveRightKeyIndex += numIsoWindow;
    }
    if (conf.isDataIndependent())
      return buildDIResult(peakStart, peakEnd, leftIndex, rightIndex, result, sModel.getBioSample().length, chroData);
    else
      //return buildDDResult(peakStart, peakEnd, leftIndex, rightIndex, result);
      return buildDDResultPeakFinding(peakStart, peakEnd, leftIndex, rightIndex, result, keys[keyIndex], chroData);

  }

  public static String peakFinding3Plex(SpectrumModel sModel, int keyIndex)
    throws IOException, CensusIndexOutOfBoundException, Exception {

    TIntLongHashMap index = sModel.getIndex();
    Object file = sModel.getGenericIndexFile(conf);

    steepRatioThreshold = conf.getSteepRatioThreshold();

    int numIsoWindow = conf.getNumOfIsolationWindow();
    int maxWindow = conf.getMaxWindow();
    int margin = conf.getMargin();

    //double[][] result = null;
    // result = new double[maxWindow*2+1+margin*2][3];
    //conf.getQuantLevel()

    LabelingResultModel[] resultArr = new LabelingResultModel[maxWindow * 2 + 1 + margin * 2];

    int leftIndex = maxWindow + margin;
    int rightIndex = maxWindow + margin + 1;
    double totalIntensity = 0;

    int steepArea = conf.getSteepArea();

    int moveLeftKeyIndex = keyIndex;

    int moveRightKeyIndex = keyIndex + 1 * numIsoWindow;

    int[] keys = sModel.getKeys();

    int initWin = 2;

    //LabelingResultModel resultModel=null;
    for (int i = 0; i < initWin; i++) ///Index=maxWindow;keys[leftPeakIndex]>=leftMinValue; )
    {
      if (moveLeftKeyIndex <= 0 || leftIndex <= 0) {
        moveLeftKeyIndex += 1 * numIsoWindow;
        leftIndex++;

        break;
      }

      resultArr[leftIndex] = readFullSpectrum3Plex(keys, moveLeftKeyIndex, index, file, sModel);

      // result[leftIndex] =

      totalIntensity += resultArr[leftIndex].getTotalIntensity();
      //totalIntensity += result[leftIndex][1];
      //totalIntensity += result[leftIndex][2];
      //totalIntensity += result[leftIndex][3];

      // if(moveLeftKeyIndex<=0 || leftIndex<=0)
      //	break;

      moveLeftKeyIndex -= numIsoWindow;
      leftIndex--;
    }

    for (int i = 0; i < initWin; i++) //leftIndex=maxWindow;keys[leftPeakIndex]>=leftMinValue; )
    {
      if (moveRightKeyIndex >= keys.length) {
        moveRightKeyIndex -= 1 * numIsoWindow;
        rightIndex--;

        break;
      }

      resultArr[rightIndex] = readFullSpectrum3Plex(keys, moveRightKeyIndex, index, file, sModel);
      totalIntensity += resultArr[rightIndex].getTotalIntensity();

      moveRightKeyIndex += numIsoWindow;
      rightIndex++;
    }

    boolean isGoingUp = true;
    boolean isHighIntensity = true;  //if the intensity is lower than one third of average intensity, this becomes false

    LabelingResultModel[] arrModel = new LabelingResultModel[steepArea];
    LabelingResultModel[] prevArrModel = new LabelingResultModel[steepArea];
    //double[][] arr = new double[steepArea][3];
    //double[][] prevArr = new double[steepArea][3];

    if (leftIndex + steepArea >= 0) {
      //for(int i=0;i<300;i++)
      while (true) {
        double area1 = 0;
        double area2 = 0;

        if (moveLeftKeyIndex - steepArea < 0) {
          break;
        }

        int steepCount = 0;
        int tempKeyIndex = moveLeftKeyIndex;
        while (true) {
          if (tempKeyIndex < 0 || steepCount < 0 || steepCount >= steepArea - 1)
            break;

          if (null != prevArrModel[0] && prevArrModel[0].getTotalIntensity() != 0) {
            for (int l = 0; l < 2; l++)
              arrModel[l] = prevArrModel[l + 1];
          } else {
            arrModel[steepCount] = readFullSpectrum3Plex(keys, tempKeyIndex, index, file, sModel);
          }

          steepCount++;
          //tempKeyIndex--;
          tempKeyIndex = tempKeyIndex - 1 * numIsoWindow;
        }

        if (tempKeyIndex < 0)
          break;

        arrModel[steepCount] = readFullSpectrum3Plex(keys, tempKeyIndex, index, file, sModel);


        for (int k = 0; k < 3; k++) {
          prevArrModel[k] = arrModel[k];
        }


        for (int j = leftIndex + 1; j <= leftIndex + steepCount + 1; j++) {
          if (null == resultArr[j]) break;

          area1 += resultArr[j].getTotalIntensity();
        }

        for (int j = 0; j <= steepCount; j++) {
          area2 += arrModel[j].getTotalIntensity();
        }

        if (area2 == 0 && area1 == 0)
          break;

        double areaRatio = (double) area2 / area1;

        /*****************************************
         * AREA1 is right side area
         * AREA2 is left side area
         *****************************************/
        if (isHighIntensity && ((double) totalIntensity / (rightIndex - leftIndex - 1) / 1.5 > area2 / 3))
          isHighIntensity = false;

        if ((isGoingUp && isHighIntensity) || areaRatio < steepRatioThreshold) {
          if (leftIndex <= 0) {
            isGoingUp = false;
            break;
          }

          resultArr[leftIndex] = arrModel[0];
          moveLeftKeyIndex -= numIsoWindow;

          if (0 == resultArr[leftIndex].getTotalIntensity())
            continue;

          leftIndex--;

          for (int j = 0; j < steepCount; j++)
            arrModel[j] = arrModel[j + 1];

          totalIntensity += arrModel[0].getTotalIntensity();

          continue;
        } else
          isGoingUp = false;


        if (areaRatio > steepRatioThreshold) {
          for (int k = 0; k < 3; k++) {
            if (resultArr[leftIndex + 1].getTotalIntensity() < arrModel[k].getTotalIntensity() || arrModel[k].getTotalIntensity() == 0)
              break;

            if (leftIndex < 0)
              break;

            resultArr[leftIndex] = arrModel[k];

            moveLeftKeyIndex -= numIsoWindow;
            leftIndex--;
          }

          break;
        }

        resultArr[leftIndex--] = arrModel[0];

        moveLeftKeyIndex -= numIsoWindow;

        for (int j = 0; j < steepCount; j++)
          arrModel[j] = arrModel[j + 1];
      }
    }

    prevArrModel = new LabelingResultModel[steepArea]; //clean up the array

    isGoingUp = true;
    isHighIntensity = true;
    if (rightIndex + steepArea < resultArr.length) {
      //for(int i=0;i<500;i++)
      while (true) {
        double area1 = 0;
        double area2 = 0;

        if (moveRightKeyIndex + steepArea >= keys.length)
          break;

        int steepCount = 0;

        int tempKeyIndex = moveRightKeyIndex;
        while (true) {
          if (tempKeyIndex >= keys.length || steepCount >= steepArea - 1 || steepCount >= keys.length)
            break;

          if (prevArrModel[0] != null && prevArrModel[0].getTotalIntensity() != 0) {
            for (int l = 0; l < 2; l++)
              arrModel[l] = prevArrModel[l + 1];
          } else {
            arrModel[steepCount] = readFullSpectrum3Plex(keys, tempKeyIndex, index, file, sModel);
          }

          steepCount++;
          tempKeyIndex = tempKeyIndex + 1 * numIsoWindow;
          //tempKeyIndex++;
        }

        if (tempKeyIndex >= keys.length)
          break;

        arrModel[steepCount] = readFullSpectrum3Plex(keys, tempKeyIndex, index, file, sModel);

        for (int k = 0; k < 3; k++) {
          prevArrModel[k] = arrModel[k];
        }

        for (int j = rightIndex - 1; j > rightIndex - steepCount - 2; j--) {
          area1 += resultArr[j].getTotalIntensity();
        }

        for (int j = 0; j <= steepCount; j++) {
          area2 += arrModel[j].getTotalIntensity();
        }

        if (area2 == 0 && area1 == 0)
          break;

        double areaRatio = (double) area2 / area1;

        /*****************************************
         * AREA1 is left side area
         * AREA2 is right side area
         *****************************************/

        if (isHighIntensity && ((double) totalIntensity / (rightIndex - leftIndex - 1) / 1.5 > area2 / 3))
          isHighIntensity = false;


        //if( ((isGoingUp && (double)totalIntensity/(rightIndex-leftIndex-1)/3<area2/3)) || areaRatio<steepRatioThreshold )
        if ((isGoingUp && isHighIntensity) || areaRatio < steepRatioThreshold) {
          if ((rightIndex + 1) >= resultArr.length) {
            isGoingUp = false;
            break;
          }

          resultArr[rightIndex] = arrModel[0];
          moveRightKeyIndex += numIsoWindow;

          if (0 == resultArr[rightIndex].getTotalIntensity())
            continue;

          rightIndex++;

          for (int j = 0; j < steepCount; j++)
            arrModel[j] = arrModel[j + 1];

          totalIntensity += arrModel[0].getTotalIntensity();

          continue;
        } else
          isGoingUp = false;

        if (area2 / area1 > steepRatioThreshold) {
          for (int k = 0; k < 3; k++) {
            if (resultArr[rightIndex - 1].getTotalIntensity() < arrModel[k].getTotalIntensity())
              break;

            if (resultArr.length <= rightIndex)
              break;

            resultArr[rightIndex] = arrModel[k];

            //System.out.println(arr[0][1] + " " + arr[1][1] + " " + arr[2][1] + " " + result[rightIndex-1][1]);

            moveRightKeyIndex += numIsoWindow;
            rightIndex++;
          }
          break;
        }

        resultArr[rightIndex++] = arrModel[0];

        moveRightKeyIndex += numIsoWindow;

        for (int j = 0; j < steepCount; j++)
          arrModel[j] = arrModel[j + 1];
      }
    }
    int peakStart = leftIndex + 1;
    int peakEnd = rightIndex - 1;

    // int specSpace = (int)(result[peakStart+1][0]-result[peakStart][0]);

    int leftMargin = margin;

    int rightMargin = margin;

    for (int i = 0; i < leftMargin; i++) {
      if (leftIndex < 0 || moveLeftKeyIndex < 0)
        break;


      resultArr[leftIndex] = readFullSpectrum3Plex(keys, moveLeftKeyIndex, index, file, sModel);

      if (0 == resultArr[leftIndex].getTotalIntensity()) {
        moveLeftKeyIndex -= numIsoWindow;
        continue;
      }

      leftIndex--;
      moveLeftKeyIndex -= numIsoWindow;
    }
    for (int i = 0; i < rightMargin; i++) {
      if (rightIndex >= resultArr.length || moveRightKeyIndex >= keys.length)
        break;


      resultArr[rightIndex] = readFullSpectrum3Plex(keys, moveRightKeyIndex, index, file, sModel);

      if (0 == resultArr[rightIndex].getTotalIntensity()) {
        moveRightKeyIndex += numIsoWindow;
        continue;
      }

      rightIndex++;
      moveRightKeyIndex += numIsoWindow;
    }

        /*LabelingResult result = new LabelingResult();
        result.setLeftIndex(leftIndex);;
        result.setPeakEnd(peakEnd);
        result.setPeakStart(peakStart);
        result.setRightIndex(rightIndex);*/

    //peak finding...
    //Range range = GaussianFitting.getGaussianPeakRangeIndex(resultArr, peakStart, peakEnd, keyIndex);

    // System.out.println(peakEnd);

        /*
        for(int i=peakStart;i<=peakEnd;i++) { //LabelingResultModel resultModel : resultArr) {
            System.out.println("------------->>\t" + resultArr[i].getScanNum() + "\t" + resultArr[i].getRetTime() +
                    resultArr[i].getResultList().get(1).getSumIntensity());
        }*/

    Range range = GaussianFitting.getGaussianPeakRangeIndex(resultArr, peakStart, peakEnd, keys[keyIndex]);


    StringBuffer sb = new StringBuffer();
//        sb.append("P ").append(resultArr[peakStart].getScanNum()).append(" ").append(resultArr[peakEnd].getScanNum()).append(" ")
//                       .append(resultArr[peakStart].getRetTime()).append(" ").append(resultArr[peakEnd].getRetTime()).append(";");

    //due to smoothing, we need to adjust peak region
    int upperBound = (int) range.getHighBound();

    if ((upperBound + 3) < resultArr.length && null != resultArr[upperBound + 3])
      upperBound += 3;

    //System.out.println(upperBound + " ========== " + resultArr.length);
    //System.out.println(resultArr[upperBound] + " ");

    sb.append("P ").append(resultArr[(int) range.getLowBound()].getScanNum()).append(" ").append(resultArr[upperBound].getScanNum()).append(" ")
      .append(resultArr[(int) range.getLowBound()].getRetTime()).append(" ").append(resultArr[upperBound].getRetTime()).append(";");


    //System.out.println(resultArr[peakStart].getScanNum() + "\t" + resultArr[peakEnd].getScanNum());

    if (leftIndex + 1 >= rightIndex - 1)
      return null;

    for (int i = leftIndex + 1; i < rightIndex - 1; i++) {
      sb.append(resultArr[i].getContent()).append(";");

      //  System.out.println(resultArr[i].getIntensityArr());
      //sb.append((double)result[i][result[i].length-2]).append(" ");
      //sb.append((double)result[i][result[i].length-1]).append(";");

      //sb.append((int)result[i][0]).append(" ").append((long)result[i][1]).append(" ").append((long)result[i][2]).append(";");
    }

    return sb.toString();

    //return buildDDResultModel(peakStart, peakEnd, leftIndex, rightIndex, resultArr);
  }


  /*
    public static String peakFinding3Plex(SpectrumModel sModel, SpecRange range, int keyIndex)
        throws IOException, CensusIndexOutOfBoundException, Exception
    {

	TIntLongHashMap index = sModel.getIndex();
        Object file = sModel.getGenericIndexFile(conf);

        steepRatioThreshold = conf.getSteepRatioThreshold();

	int numIsoWindow = conf.getNumOfIsolationWindow();
        int maxWindow=conf.getMaxWindow();
        int margin = conf.getMargin();

	double[][] result = null;
        //conf.getQuantLevel()

	if(conf.isDataIndependent())
	    result = new double[maxWindow*2+1+margin*2][4*sModel.getBioSample().length+3]; //scan #, sample intensity, ref intensity
	else
	    result = new double[maxWindow*2+1+margin*2][3];

        int leftIndex=maxWindow+margin;
        int rightIndex=maxWindow+margin+1;


        //double leftTotalIntensity=0;
        //double rightTotalIntensity=0;
        double totalIntensity=0;

        int steepArea = conf.getSteepArea();

        int moveLeftKeyIndex = keyIndex;

        int moveRightKeyIndex = keyIndex+1*numIsoWindow;

        int[] keys = sModel.getKeys();
        int initWin=2;

        for(int i=0;i<initWin;i++) ///Index=maxWindow;keys[leftPeakIndex]>=leftMinValue; )
        {
	    if(moveLeftKeyIndex<=0 || leftIndex<=0)
	    {
		moveLeftKeyIndex += 1*numIsoWindow;
		leftIndex++;

		break;
	    }

	    if(conf.isDataIndependent()) //data independent
	    {
		result[leftIndex] = readSpectrum(keys, moveLeftKeyIndex, index, sModel.getDiff(), file, sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());

	    }
	    else //data dependent
	    {
		if(sModel.isHighRes())
		    result[leftIndex] = readFullSpectrum3Plex(keys, moveLeftKeyIndex, index, file, sModel.getSamIsoArr(), sModel.getRefIsoArr(), sModel.getRef2IsoArr());
		else
		    result[leftIndex] = readFullSpectrum(keys, moveLeftKeyIndex, index, file, sModel.getSamStartMass(), sModel.getSamEndMass(), sModel.getRefStartMass(), sModel.getRefEndMass());

	    }
            totalIntensity += result[leftIndex][1];
            totalIntensity += result[leftIndex][2];
            totalIntensity += result[leftIndex][3];

	   // if(moveLeftKeyIndex<=0 || leftIndex<=0)
	//	break;

	    moveLeftKeyIndex -= numIsoWindow;
            leftIndex--;
        }

        for(int i=0;i<initWin;i++) //leftIndex=maxWindow;keys[leftPeakIndex]>=leftMinValue; )
        {
            if(moveRightKeyIndex>=keys.length)
            {
		moveRightKeyIndex -= 1*numIsoWindow;
                rightIndex--;

                break;
            }

	    if(conf.isDataIndependent()) //data independent
	    {
		result[rightIndex] = readSpectrum(keys, moveRightKeyIndex, index, sModel.getDiff(), file, sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());
	    }
	    else //data dependent
	    {
		if(sModel.isHighRes())
		    result[rightIndex] = readFullSpectrum3Plex(keys, moveRightKeyIndex, index, file, sModel.getSamIsoArr(), sModel.getRefIsoArr(), sModel.getRef2IsoArr());
		else
		    result[rightIndex] = readFullSpectrum(keys, moveRightKeyIndex, index, file, sModel.getSamStartMass(), sModel.getSamEndMass(), sModel.getRefStartMass(), sModel.getRefEndMass());
	    }

            totalIntensity += result[rightIndex][1];
            totalIntensity += result[rightIndex][2];
            totalIntensity += result[rightIndex][3];

            moveRightKeyIndex += numIsoWindow;
            rightIndex++;
        }

        boolean isGoingUp=true;
        boolean isHighIntensity=true;  //if the intensity is lower than one third of average intensity, this becomes false
        double[][] arr = new double[steepArea][3];
        double[][] prevArr = new double[steepArea][3];

        if(leftIndex+steepArea>=0)
        {
            //for(int i=0;i<300;i++)
	    while(true)
            {
                double area1=0;
                double area2=0;

                if(moveLeftKeyIndex-steepArea<0)
                {
                    break;
                }

                int steepCount=0;
                int tempKeyIndex = moveLeftKeyIndex;
                while(true)
                {
		    if(tempKeyIndex<0 || steepCount<0 || steepCount>=steepArea-1)
			break;

		    if(prevArr[0][0]!=0)
		    {
			for(int l=0;l<2;l++)
			    for(int m=0;m<3;m++)
				arr[l][m] = prevArr[l+1][m];
		    }
		    else
		    {
			if(conf.isDataIndependent()) //data independent
			    arr[steepCount] = readSpectrum(keys, tempKeyIndex, index, sModel.getDiff(), file, sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());
			else
			{
			    if(sModel.isHighRes())
				arr[steepCount] = readFullSpectrum3Plex(keys, tempKeyIndex, index, file, sModel.getSamIsoArr(), sModel.getRefIsoArr(), sModel.getRef2IsoArr());
			    else
				arr[steepCount] = readFullSpectrum(keys, tempKeyIndex, index, file, sModel.getSamStartMass(), sModel.getSamEndMass(), sModel.getRefStartMass(), sModel.getRefEndMass());
			}
		    }

                    steepCount++;
                    //tempKeyIndex--;
		    tempKeyIndex = tempKeyIndex - 1*numIsoWindow;
                }

		if(tempKeyIndex<0)
		    break;

		if(conf.isDataIndependent()) //data independent
		    arr[steepCount] = readSpectrum(keys, tempKeyIndex, index, sModel.getDiff(), file, sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());
		else
		{
		    if(sModel.isHighRes())
			arr[steepCount] = readFullSpectrum3Plex(keys, tempKeyIndex, index, file, sModel.getSamIsoArr(), sModel.getRefIsoArr(), sModel.getRef2IsoArr());
		    else
			arr[steepCount] = readFullSpectrum(keys, tempKeyIndex, index, file, sModel.getSamStartMass(), sModel.getSamEndMass(), sModel.getRefStartMass(), sModel.getRefEndMass());
		}

		for(int k=0;k<3;k++)
		{
		    for(int l=0;l<3;l++)
			prevArr[k][l] = arr[k][l];
		}


                for(int j=leftIndex+1;j<=leftIndex+steepCount+1;j++)
                {
                    area1 += result[j][1];
                    area1 += result[j][2];
                    area1 += result[j][3];
                }

                for(int j=0;j<=steepCount;j++)
                {
                    area2 += arr[j][1];
                    area2 += arr[j][2];
                    area2 += arr[j][3];
                }

                if(area2==0 && area1==0)
                    break;

                double areaRatio = (double)area2/area1;

                ****************************************
                 * AREA1 is right side area
                 * AREA2 is left side area
                 ****************************************
		if(isHighIntensity && ((double)totalIntensity/(rightIndex-leftIndex-1)/1.5>area2/3))
		    isHighIntensity = false;

                if( (isGoingUp && isHighIntensity) || areaRatio<steepRatioThreshold )
                {
		    if(leftIndex<=0)
		    {
                        isGoingUp=false;
                        break;
		    }

		    //System.out.println("peak gogo" + " " + leftIndex + " " + moveLeftKeyIndex + " " + index + " " + sModel.getDiff());
                    result[leftIndex] = arr[0];
		    moveLeftKeyIndex -= numIsoWindow;

		    if(0 == result[leftIndex][0])
			continue;

		    leftIndex--;

                    for(int j=0;j<steepCount;j++)
                        arr[j] = arr[j+1];

                    totalIntensity += arr[0][1];
                    totalIntensity += arr[0][2];

                    continue;
                }
                else
                    isGoingUp=false;


                if(areaRatio>steepRatioThreshold)
                {
                    for(int k=0;k<3;k++)
                    {
                        if(result[leftIndex+1][1]<arr[k][1] || arr[k][0]==0)
                            break;

			if(leftIndex<0)
			    break;

                        result[leftIndex] = arr[k];

			moveLeftKeyIndex -= numIsoWindow;
                        leftIndex--;
                    }

                    break;
                }

                result[leftIndex--] = arr[0];

		moveLeftKeyIndex -= numIsoWindow;

                for(int j=0;j<steepCount;j++)
                    arr[j] = arr[j+1];
            }
        }

	prevArr = new double[steepArea][3]; //clean up the array

        isGoingUp = true;
	isHighIntensity = true;
        if(rightIndex+steepArea<result.length)
        {
            //for(int i=0;i<500;i++)
            while(true)
            {
                double area1=0;
                double area2=0;

                if(moveRightKeyIndex+steepArea>=keys.length)
                    break;

                int steepCount=0;

                int tempKeyIndex = moveRightKeyIndex;
                while(true)
                {
		    if(tempKeyIndex>=keys.length || steepCount>= steepArea-1 || steepCount>=keys.length)
			break;

		    if(prevArr[0][0]!=0)
		    {
			for(int l=0;l<2;l++)
			    for(int m=0;m<3;m++)
				arr[l][m] = prevArr[l+1][m];
		    }
		    else
		    {
			if(conf.isDataIndependent()) //data independent
			    arr[steepCount] = readSpectrum(keys, tempKeyIndex, index, sModel.getDiff(), file, sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());
			else if(sModel.isHighRes())
			    arr[steepCount] = readFullSpectrum3Plex(keys, tempKeyIndex, index, file, sModel.getSamIsoArr(), sModel.getRefIsoArr(), sModel.getRef2IsoArr());
			else
			    arr[steepCount] = readFullSpectrum(keys, tempKeyIndex, index, file, sModel.getSamStartMass(), sModel.getSamEndMass(), sModel.getRefStartMass(), sModel.getRefEndMass());
		    }

                    steepCount++;
		    tempKeyIndex = tempKeyIndex + 1*numIsoWindow;
                    //tempKeyIndex++;
                }

		if(tempKeyIndex>=keys.length)
		    break;

		if(conf.isDataIndependent()) //data independent
		    arr[steepCount] = readSpectrum(keys, tempKeyIndex, index, sModel.getDiff(), file, sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());
		else if(sModel.isHighRes())
		    arr[steepCount] = readFullSpectrum3Plex(keys, tempKeyIndex, index, file, sModel.getSamIsoArr(), sModel.getRefIsoArr(), sModel.getRef2IsoArr());
		else
		    arr[steepCount] = readFullSpectrum(keys, tempKeyIndex, index, file, sModel.getSamStartMass(), sModel.getSamEndMass(), sModel.getRefStartMass(), sModel.getRefEndMass());

		for(int k=0;k<3;k++)
		{
		    for(int l=0;l<3;l++)
			prevArr[k][l] = arr[k][l];
		}

                for(int j=rightIndex-1;j>rightIndex-steepCount-2;j--)
                {
                    area1 += result[j][1];
                    area1 += result[j][2];
                    area1 += result[j][3];
                }

                for(int j=0;j<=steepCount;j++)
                {
                    area2 += arr[j][1];
                    area2 += arr[j][2];
                }

                if(area2==0 && area1==0)
                    break;

                double areaRatio = (double)area2/area1;

                ****************************************
                 * AREA1 is left side area
                 * AREA2 is right side area
                 ****************************************

		if(isHighIntensity && ((double)totalIntensity/(rightIndex-leftIndex-1)/1.5>area2/3))
		    isHighIntensity = false;


                //if( ((isGoingUp && (double)totalIntensity/(rightIndex-leftIndex-1)/3<area2/3)) || areaRatio<steepRatioThreshold )
                if( (isGoingUp && isHighIntensity) || areaRatio<steepRatioThreshold )
                {
                    if((rightIndex+1)>=result.length)
                    {
                        isGoingUp=false;
                        break;
                    }

                    result[rightIndex] = arr[0];
		    moveRightKeyIndex += numIsoWindow;

		    if(0 == result[rightIndex][0])
			continue;

		    rightIndex++;

                    for(int j=0;j<steepCount;j++)
                        arr[j] = arr[j+1];

                    //rightTotalIntensity += arr[0][1];
                    totalIntensity += arr[0][1];
                    totalIntensity += arr[0][2];
                    totalIntensity += arr[0][3];

                    continue;
                }
                else
                    isGoingUp=false;

                if(area2/area1>steepRatioThreshold)
                {
                    for(int k=0;k<3;k++)
                    {
                        if(result[rightIndex-1][1]<arr[k][1])
                            break;

			if(result.length<=rightIndex)
			    break;

                        result[rightIndex] = arr[k];

                        //System.out.println(arr[0][1] + " " + arr[1][1] + " " + arr[2][1] + " " + result[rightIndex-1][1]);

			moveRightKeyIndex += numIsoWindow;
                        rightIndex++;
                    }
                     break;
                }

                result[rightIndex++] = arr[0];

		moveRightKeyIndex += numIsoWindow;

                for(int j=0;j<steepCount;j++)
                    arr[j] = arr[j+1];
            }
        }
	int peakStart = leftIndex+1;
	int peakEnd = rightIndex-1;

        int dtaStart = (range==null?0:range.getMin());
        int dtaEnd = (range==null?0:range.getMax());
        int specSpace = (int)(result[peakStart+1][0]-result[peakStart][0]);

//System.out.println( result[peakStart+1][0] + " " + result[peakStart][0] );
        int diff = (int)(result[peakStart][0]-dtaStart)/specSpace;
        int leftMargin = margin;

        int rightMargin = margin;
        if(diff>0)
            leftMargin += diff;

        diff = (int)(dtaEnd-result[rightIndex-1][0])/specSpace;

        if(diff>0)
            rightMargin += diff;

	for(int i=0;i<leftMargin;i++)
        {
            if(leftIndex<0 || moveLeftKeyIndex<0)
                break;

	    if(conf.isDataIndependent()) //data independent
		result[leftIndex] = readSpectrum(keys, moveLeftKeyIndex, index, sModel.getDiff(), file, sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());
	    else if(sModel.isHighRes())
		result[leftIndex] = readFullSpectrum3Plex(keys, moveLeftKeyIndex, index, file, sModel.getSamIsoArr(), sModel.getRefIsoArr(), sModel.getRef2IsoArr());
	    else
		result[leftIndex] = readFullSpectrum(keys, moveLeftKeyIndex, index, file, sModel.getSamStartMass(), sModel.getSamEndMass(), sModel.getRefStartMass(), sModel.getRefEndMass());

	    if(0 == result[leftIndex][0])
	    {
		moveLeftKeyIndex -= numIsoWindow;
		continue;
	    }

            leftIndex--;
	    moveLeftKeyIndex -= numIsoWindow;
        }
        for(int i=0;i<rightMargin;i++)
        {
            if(rightIndex>=result.length || moveRightKeyIndex>=keys.length)
                break;

	    if(conf.isDataIndependent()) //data independent
		result[rightIndex] = readSpectrum(keys, moveRightKeyIndex, index, sModel.getDiff(), file, sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());
	    else if(sModel.isHighRes())
		result[rightIndex] = readFullSpectrum3Plex(keys, moveRightKeyIndex, index, file, sModel.getSamIsoArr(), sModel.getRefIsoArr(), sModel.getRef2IsoArr());
	    else
		result[rightIndex] = readFullSpectrum(keys, moveRightKeyIndex, index, file, sModel.getSamStartMass(), sModel.getSamEndMass(), sModel.getRefStartMass(), sModel.getRefEndMass());

	    if(0 == result[rightIndex][0])
	    {
		moveRightKeyIndex += numIsoWindow;
		continue;
	    }

            rightIndex++;
	    moveRightKeyIndex += numIsoWindow;
        }
	if(conf.isDataIndependent())
	    return buildDIResult(peakStart, peakEnd, leftIndex, rightIndex, result, sModel.getBioSample().length);
	else
	    return buildDDResult(peakStart, peakEnd, leftIndex, rightIndex, result);
    }
   */
  //non-labeling for low resolution and full ms
  public static String calculateNonlabelMS(
    int keyIndex,
    String fileName,
    Hashtable<String, IndexedFile> htInput,
    int[][][] pathArrayInput,
    double startMass,
    double endMass,
    SpecRange range,
    String pathInput)
    throws IOException, Exception {
    ht = htInput;

    conf = Configuration.getInstance();
    align = conf.isAlign();

    massTolerance = conf.getMassTolerance();
    pathArray = pathArrayInput;
    pathFileName = pathInput + fileName;

    SpectrumModel sModel = new SpectrumModel();
    sModel.setHighRes(false);
    sModel.setSamStartMass(startMass);
    sModel.setSamEndMass(endMass);

    return peakFinding(sModel, range, conf, keyIndex, null);
  }

  //non-labeling for full ms
  public static String calculateNonlabelMS(
    int keyIndex,
    String fileName,
    Hashtable<String, IndexedFile> htInput,
    int[][][] pathArrayInput,
    double[] samIsoArr,
    SpecRange range,
    String pathInput)
    throws IOException, Exception {
    ht = htInput;

    conf = Configuration.getInstance();
    align = conf.isAlign();


    massTolerance = conf.getMassTolerance();
    pathArray = pathArrayInput;
    pathFileName = pathInput + fileName;

    return peakFinding(null, range, conf, keyIndex, samIsoArr);
  }

  //non-labeling for tandem scans
  public static String calculateNonlabelMS(
    int keyIndex,
    String fileName,
    Hashtable<String, IndexedFile> htInput,
    int[][][] pathArrayInput,
    double[] samIsoArr,
    SpecRange range,
    String pathInput,
    double[][] bionSample,
//            double[][] bionRef,
    double[][] yionSample,
//            double[][] yionRef,
    IndexedFile iFile)
    throws IOException, Exception {
    ht = htInput;

    conf = Configuration.getInstance();
    align = conf.isAlign();


    massTolerance = conf.getMassTolerance();
    pathArray = pathArrayInput;
    pathFileName = pathInput + fileName;

    SpectrumModel sModel = new SpectrumModel();
    sModel.setHighRes(false);

    //sModel.setKeys(keys);
    //sModel.setIndex(index);
    sModel.setIFile(iFile);
    //sModel.setFile(file);
  /*
      for(int i=0;i<bionSample.length;i++)
      {
	  for(int j=0;j<bionSample[i].length;j++)
	      System.out.print(" " + bionSample[i][j]);

	  System.out.println("");
      }
*/
    sModel.setBioSample(bionSample);
    sModel.setYioSample(yionSample);
    //       sModel.setBioRef(bionRef);
    //       sModel.setYioRef(yionRef);

    //sModel.setDiff(diff);

    return peakFinding(sModel, range, conf, keyIndex, samIsoArr);
    //return null;


  }

  //MultipleMs1Labeling
  public static String calculateFullMS3Plex(
    int keyIndex,
    IndexedFile iFile,
    double[] samIsoArr,
    double[] refIsoArr,
    double[] ref2IsoArr
  )

    throws IOException, CensusIndexOutOfBoundException, Exception {
    TIntLongHashMap index = iFile.getMsIndex();
    RandomAccessFile file = iFile.getFile();
    MzxmlSpectrumReader mzReader = iFile.getMzreader();
    int[] keys = iFile.getKeys();
    double[] retArr = iFile.getRtArr();

    //TDoubleIntHashMap rtScanMap = iFile.getRetentonToScanMap();
    //TIntDoubleHashMap scanRtMap = iFile.getScanRtMap();
    //TDoubleIntHashMap rtScanMap2 = iFile.getRtScanNumMap();

    conf = Configuration.getInstance();
    massTolerance = conf.getMassTolerance();

    SpectrumModel sModel = new SpectrumModel();
    sModel.setHighRes(true);
    sModel.setKeys(keys);
    sModel.setIndex(index);
    sModel.setIFile(iFile);
    sModel.setFile(file);
    sModel.setMzxmlreader(mzReader);
    sModel.setSamIsoArr(samIsoArr);
    sModel.setRefIsoArr(refIsoArr);
    sModel.setRef2IsoArr(ref2IsoArr);

    //sModel.setRetArr(retArr);
    //   sModel.setRtScanMap(rtScanMap);


    return peakFinding3Plex(sModel, keyIndex);
  }

  //labelfree
  public static String calculateLabelfreeFullMS(
    int keyIndex,
    IndexedFile iFile,
    double[] samIsoArr,
    SpecRange range)

    throws IOException, CensusIndexOutOfBoundException, Exception {
    TIntLongHashMap index = iFile.getMsIndex();
    RandomAccessFile file = iFile.getFile();
    MzxmlSpectrumReader mzReader = iFile.getMzreader();
    int[] keys = iFile.getKeys();
    conf = Configuration.getInstance();
    massTolerance = conf.getMassTolerance();

    SpectrumModel sModel = new SpectrumModel();
    sModel.setHighRes(true);
    sModel.setKeys(keys);
    sModel.setIndex(index);
    sModel.setIFile(iFile);
    sModel.setFile(file);
    sModel.setMzxmlreader(mzReader);
    sModel.setSamIsoArr(samIsoArr);

    return peakFinding(sModel, range, keyIndex, null);
  }

  public static String calculateFullMS(
    int keyIndex,
    IndexedFile iFile,
    double[] samIsoArr,
    double[] refIsoArr,
    SpecRange range, DisplayData.DisplayChroData chroData)

    throws IOException, CensusIndexOutOfBoundException, Exception {
    TIntLongHashMap index = iFile.getMsIndex();
    RandomAccessFile file = iFile.getFile();
    MzxmlSpectrumReader mzReader = iFile.getMzreader();
    int[] keys = iFile.getKeys();
    conf = Configuration.getInstance();
    massTolerance = conf.getMassTolerance();

    SpectrumModel sModel = new SpectrumModel();
    sModel.setHighRes(true);
    sModel.setKeys(keys);
    sModel.setIndex(index);
    sModel.setIFile(iFile);
    sModel.setFile(file);
    sModel.setMzxmlreader(mzReader);
    sModel.setSamIsoArr(samIsoArr);
    sModel.setRefIsoArr(refIsoArr);

    return peakFinding(sModel, range, keyIndex, chroData);
  }

  public static void getSpectrumArr(int[] keys, int curIndex, TIntLongHashMap index, Object ofile) throws IOException, CensusIndexOutOfBoundException {
    RandomAccessFile file = (RandomAccessFile) ofile;

    if (keys.length <= curIndex)
      throw new CensusIndexOutOfBoundException();

    long startPos = index.get(keys[curIndex]);
    long endPos;

    file.seek(startPos);

    if ((curIndex + 1) >= keys.length)
      endPos = file.length();
    else
      endPos = index.get(keys[curIndex + 1]);

    int byteSize = (int) (endPos - startPos);

//System.out.println("size...  \t" + endPos + "\t" + startPos + "\t" + byteSize + "\t" +  index.get(keys[curIndex])  + "\t" + keys[curIndex]);

    byte[] bytes = new byte[byteSize];
    file.readFully(bytes);
    // System.out.println(keys[curIndex] + " " + startPos + " " + endPos + " " + byteSize);

    char ch;
    int pos = 0;

    ch = (char) bytes[pos];

    try {
      //Remove Z, S, I, D lines
      while ((ch = (char) bytes[pos]) == 'Z' || ch == 'S' || ch == 'I' || ch == 'D') {
        while (ch != CARRIAGE_RETURN) {
          pos++;
          ch = (char) bytes[pos];
        }

        pos++;
      }


      double[][] specResult = parseSpectra(bytes);
      massArr = specResult[0];
      intArr = specResult[1];
      if (keys[curIndex] == CalcUtil.getScanToSearch()) {
        CalcUtil.setXyMassArr(massArr);
        CalcUtil.setXyIntArr(intArr);
      }
    } catch (Exception e) {
      System.out.println("no peaks");
    }

  }


  //Low resolution for label free only for now
  //for both full and tandem
  //
  public static double[] readSpectrum(
    SpectrumModel sModel,
    //int[] keys,
    int curIndex
    //TIntLongHashMap index
    //int diff,
    //int refIndexKey,
    //RandomAccessFile file

  ) throws IOException, CensusIndexOutOfBoundException {
    if (curIndex < 0)
      return null;

    double[][] bioSample = sModel.getBioSample();
    double[][] yioSample = sModel.getYioSample();

    NonLabelMappingModel mapModel = conf.getMapModel();
    Vector<String> arr = mapModel.getPathFileNameList();

    int ionLength = bioSample.length;
    double[] tmpArr = null;

    Hashtable<String, Set> childHt = null;
    double precursor = -1;
    IndexedFile iFile = null;

    //ION_START_INDEX = arr.size();


    if (conf.getQuantLevel() == 2) {
      childHt = mapModel.getChildHashtableByMS2(pathFileName, curIndex);
      ION_START_INDEX = arr.size();


      tmpArr = new double[ionLength * 2 * arr.size() + ION_START_INDEX * 2]; //scan arr, intsum arr, sam1 b y arr, sam2 by arr, ..; repeat

      iFile = sModel.getIFile();
      precursor = iFile.getPrecursorMap().get(iFile.getKeys()[curIndex]);
    } else {
      childHt = mapModel.getChildHashtable(pathFileName, curIndex);

      tmpArr = new double[arr.size() * 2];
    }

    if (null == childHt)
      return tmpArr;

    int refIndex = mapModel.getRefIndex();

    // for(int i=0;i<samIsoArr.length;i++)
    //     System.out.println("-->> " + samIsoArr[i]);

    int tmpIndex = 0;
    int targetIndex = 1;

    for (Iterator<String> itr = arr.iterator(); itr.hasNext(); ) {
      String eachName = itr.next();

      //    if(conf.getQuantLevel()==2)
      //	eachName = eachName.replace(".ms2", ".ms1");

      Set<Integer> tmpSet = childHt.get(eachName);


      IndexedFile tmpIndexFile = ht.get(eachName);
      int[] keys = tmpIndexFile.getKeys();
      RandomAccessFile file = tmpIndexFile.getFile();
      TIntLongHashMap index = tmpIndexFile.getMsIndex();
      TIntDoubleHashMap p = tmpIndexFile.getPrecursorMap();

      int count = 0;
      double intensity = 0;

      int firstIndex = 0;
      int firstScanIndex = 0;

      int ms2MapIndex = -1; // for tandem quant


      for (Iterator<Integer> indItr = tmpSet.iterator(); indItr.hasNext(); ) {
        Integer eachIndex = indItr.next();

        //if(tandem)
        if (conf.getQuantLevel() == 2) {
          if (firstIndex == 0) {
            firstScanIndex = eachIndex.intValue() - 1;
            if (eachIndex <= 0)
              firstScanIndex = 0;

            firstScanIndex = mapModel.getMs2Index(eachName, firstScanIndex, precursor, conf.getMassTolerance(), ht.get(eachName)); //index
            if (firstScanIndex < 0)
              return tmpArr;

            firstIndex = ht.get(eachName).getKeys()[firstScanIndex]; //scan number

            //System.out.println(eachIndex.intValue()  + " " + firstScanIndex);
          }

          ms2MapIndex = mapModel.getMs2Index(eachName, eachIndex.intValue(), precursor, conf.getMassTolerance(), ht.get(eachName));

          tmpArr[tmpIndex] = keys[ms2MapIndex]; //scan number
          //System.out.println(keys[ms2MapIndex]);

          //update massarr
          getSpectrumArr(keys, ms2MapIndex, index, file);
        } else {
          if (firstIndex == 0) {
            firstScanIndex = eachIndex.intValue() - 1;

            firstIndex = keys[firstScanIndex];
          }

          //update massarr
          getSpectrumArr(keys, eachIndex.intValue(), index, file);
        }

        //System.out.println("each name--" + mapModel.getMs2Index(eachName, eachIndex.intValue(),  );

        //                intensity += intensitySum(massArr, intArr, samIsoArr);
        //System.out.println(intensity);

        count++;
      }


      count++;

      long sampleSum = 0;
      if (tmpIndex == refIndex) {
        int tmpScanIndex;
        if (conf.getQuantLevel() == 2)
          tmpScanIndex = ms2MapIndex;
        else
          tmpScanIndex = tmpSet.iterator().next().intValue() - 1;

        if (conf.getQuantLevel() == 2) {
          if (conf.isAlign()) {
            //tmpArr[0] = keys[tmpScanIndex];


            //        System.out.println(tmpScanIndex + " ----"  + keys[tmpScanIndex]);

            for (int k = 0; k < bioSample.length; k++) {
              //tmpArr[k] = keys[tmpScanIndex];  //scan number

              int bIndex = ION_START_INDEX * 2 + k;
              int yIndex = ION_START_INDEX * 2 + ionLength + k;

              tmpArr[bIndex] = intensitySum(massArr, intArr, bioSample[k][1], bioSample[k][2]);
              sampleSum += tmpArr[bIndex];
              tmpArr[yIndex] = intensitySum(massArr, intArr, yioSample[k][1], yioSample[k][2]);
              sampleSum += tmpArr[yIndex];


//System.out.println("==11== " + bIndex + " " + yIndex);
              //	sampleSum += tmpArr[ION_START_INDEX+bioSample.length+k];
            }

            tmpArr[ION_START_INDEX] = sampleSum;
            //                      System.out.println("sample1--" + ION_START_INDEX + " " + sampleSum);
          } else {
            System.out.println("non alignment for tandem spectra is not supported yet. contact to rpark@scripps.edu");
          }

        } else {
          if (conf.isAlign()) {
            tmpArr[0] = keys[tmpScanIndex];
          } else {
            tmpArr[0] = conf.getRetArr()[tmpIndex][tmpScanIndex] * 100;

          }

        }
      } else {
        if (conf.getQuantLevel() == 2) {
          for (int k = 0; k < bioSample.length; k++) {

            int bIndex = ION_START_INDEX * 2 + k + ionLength * (targetIndex * 2);
            int yIndex = ION_START_INDEX * 2 + k + ionLength * (targetIndex * 2 + 1);

            tmpArr[bIndex] = intensitySum(massArr, intArr, bioSample[k][1], bioSample[k][2]);
            sampleSum += tmpArr[bIndex];
            tmpArr[yIndex] = intensitySum(massArr, intArr, yioSample[k][1], yioSample[k][2]);
            sampleSum += tmpArr[yIndex];
//System.out.println("==22== " + bIndex + " " + yIndex + " " + targetIndex);

          }


          tmpArr[ION_START_INDEX + targetIndex] = sampleSum;
          //System.out.println("sample--" + (ION_START_INDEX+ targetIndex) + " " + sampleSum);
        } else {
          if (conf.isAlign()) {
            tmpArr[targetIndex] = firstIndex;

            System.out.println(firstIndex);

          } else {
            tmpArr[targetIndex] = conf.getRetArr()[tmpIndex][firstScanIndex] * 100;
          }

        }

        targetIndex++;
      }

      tmpArr[tmpArr.length / 2 + tmpIndex] = intensity / count;
      tmpIndex++;
    }


    return tmpArr;
  }

  //High resolution for non-labeling only
  //Low resolution, too now.
  //There is no reference
  public static double[] readFullSpectrum(
    int curIndex,
    TIntLongHashMap nouseindex,
    RandomAccessFile nousefile,
    double[] samIsoArr,
    NonLabelMappingModel mapModel,
    SpectrumModel sModel,
    Hashtable<Integer, double[]> resultHt
  ) throws IOException, CensusIndexOutOfBoundException {
    //populate massarr field

//	getSpectrumArr(keys, curIndex, index, file);
    Vector<String> arr = mapModel.getPathFileNameList();

    double[] tmpArr = new double[arr.size() * 2];

    Hashtable<String, Set> childHt = mapModel.getChildHashtable(pathFileName, curIndex);

    if (null == childHt)
      return tmpArr;

    int refIndex = mapModel.getRefIndex();

//	System.out.println( "11--->>" + arr.get(refIndex) );
//	System.out.println( "22--->>" + childHt.get(arr.get(refIndex)) );

    int hashCode = childHt.get(arr.get(refIndex)).hashCode();
    double[] setArr = resultHt.get(hashCode);
    if (null != setArr)
      return setArr;

    // for(int i=0;i<samIsoArr.length;i++)
    //     System.out.println("-->> " + samIsoArr[i]);

    int tmpIndex = 0;
    int targetIndex = 1;

    for (Iterator<String> itr = arr.iterator(); itr.hasNext(); ) {
      String eachName = itr.next();
      Set<Integer> tmpSet = childHt.get(eachName);


      IndexedFile tmpIndexFile = ht.get(eachName);

      int[] keys = tmpIndexFile.getKeys();

      RandomAccessFile file = tmpIndexFile.getFile();
      TIntLongHashMap index = tmpIndexFile.getMsIndex();

      int count = 0;
      double intensity = 0;

      int firstIndex = 0;
      int firstScanIndex = 0;
      for (Iterator<Integer> indItr = tmpSet.iterator(); indItr.hasNext(); ) {
        Integer eachIndex = indItr.next();

        if (firstIndex == 0) {
          firstScanIndex = eachIndex.intValue() - 1;

          if (firstScanIndex >= keys.length)
            firstIndex = keys[keys.length - 1];
          else
            firstIndex = keys[firstScanIndex];
        }


        if (eachIndex >= keys.length)
          continue;

        getSpectrumArr(keys, eachIndex.intValue(), index, file);


        //double[] tempArr = null;

        if (conf.isHighRes()) {
          double[] tempArr = intensitySum(massArr, intArr, samIsoArr);
          intensity += tempArr[0];
        } else
          intensity += intensitySum(massArr, intArr, sModel.getSamStartMass(), sModel.getSamEndMass());

        //
        //System.out.println(intensity);

        count++;
      }

      if (tmpIndex == refIndex) {
        int tmpScanIndex = tmpSet.iterator().next().intValue() - 1;

        if (tmpScanIndex >= keys.length)
          continue;

        if (conf.isAlign())
          tmpArr[0] = keys[tmpScanIndex];
        else
          tmpArr[0] = conf.getRetArr()[tmpIndex][tmpScanIndex] * 100;
      } else {
        if (conf.isAlign()) {
          tmpArr[targetIndex] = firstIndex;

        } else {
          tmpArr[targetIndex] = conf.getRetArr()[tmpIndex][firstScanIndex] * 100;
        }

        targetIndex++;

      }

      tmpArr[tmpArr.length / 2 + tmpIndex] = intensity / count;
      tmpIndex++;
      //if(true)
      //	throw new CensusIndexOutOfBoundException();

    }

    resultHt.put(hashCode, tmpArr);

    return tmpArr;
  }


  //High resolution //lalebling data
  public static double[] readFullSpectrum(
    int[] keys,
    int curIndex,
    TIntLongHashMap index,
    //RandomAccessFile file,
    Object reader,
    double[] samIsoArr,
    double[] refIsoArr
  ) throws IOException, CensusIndexOutOfBoundException, Exception {

    double result[] = new double[9];
    if (curIndex < 0)
      return result;

    if (reader instanceof RandomAccessFile) {
      getSpectrumArr(keys, curIndex, index, reader);
    } else if (reader instanceof MzxmlSpectrumReader) {
      MzxmlSpectrumReader mzReader = (MzxmlSpectrumReader) reader;
      readFullSpectrumMzXml(keys, curIndex, mzReader);
    }


    result[0] = keys[curIndex];
    double[] tempArr = intensitySum(massArr, intArr, samIsoArr);
    result[1] = tempArr[0];     //intensity
    result[3] = samIsoArr.length;
    result[5] = tempArr[1];      //Found iso num
    result[7] = tempArr[2];       //tolerance

    if (null != refIsoArr) {
      tempArr = intensitySum(massArr, intArr, refIsoArr);
      result[2] = tempArr[0];
      result[4] = refIsoArr.length;
      result[6] = tempArr[1];
      result[8] = tempArr[2];
    }

      /*
    System.out.println("============");
    System.out.println("============");
    System.out.println("============");
    System.out.println("============");
    System.out.println("============");
    System.out.println("============");
    System.out.println("============" + keys[curIndex]);

    if(keys[curIndex]==4294) {
    for(double d : result)
	System.out.println(d);

    System.exit(0);
}
*/

    return result;
  }

  /*
    //High resolution //lalebling data
    public static double[] readFullSpectrum3Plex(
            int[] keys,
            int curIndex,
            TIntLongHashMap index,
            //RandomAccessFile file,
            Object reader,
            double[] samIsoArr,
	    double[] refIsoArr,
	    double[] ref2IsoArr
            ) throws IOException, CensusIndexOutOfBoundException, Exception
    {

	double result[] = new double[13];
	if(curIndex<0)
	    return result;

	if(reader instanceof RandomAccessFile)
	{
	    getSpectrumArr(keys, curIndex, index, reader);
	}
	else if(reader instanceof MzxmlSpectrumReader)
	{
	    MzxmlSpectrumReader mzReader = (MzxmlSpectrumReader)reader;
	    readFullSpectrumMzXml(keys, curIndex, mzReader);
	}

	result[0] = keys[curIndex];
	double[] tempArr = intensitySum(massArr, intArr, samIsoArr);
        result[1] = tempArr[0];     //intensity
        result[4] = samIsoArr.length;
        result[7] = tempArr[1];      //Found iso num
        result[10] = tempArr[2];	     //tolerance

	if(null != refIsoArr) {
	    tempArr = intensitySum(massArr, intArr, refIsoArr);
	    result[2] = tempArr[0];
	    result[5] = refIsoArr.length;
	    result[8] = tempArr[1];
	    result[11] = tempArr[2];
	}

	if(null != ref2IsoArr) {
	    tempArr = intensitySum(massArr, intArr, ref2IsoArr);
	    result[3] = tempArr[0];
	    result[6] = ref2IsoArr.length;
	    result[9] = tempArr[1];
	    result[12] = tempArr[2];
	}

	return result;
    }
*/
  //High resolution //lalebling data
  public static LabelingResultModel readFullSpectrum3Plex(
    int[] keys,
    int curIndex,
    TIntLongHashMap index,
    //RandomAccessFile file,
    Object reader,
    SpectrumModel sModel
  ) throws IOException, CensusIndexOutOfBoundException, Exception {

    //double result[] = new double[13];
    if (curIndex < 0)
      return null;

    getSpectrumArr(keys, curIndex, index, reader);

    LabelingResultModel result = new LabelingResultModel();
    result.setScanNum(keys[curIndex]);
    result.setRetTime(sModel.getIFile().getScanRtMap().get(keys[curIndex]));


    //  System.out.println("=======================" + sModel.getIFile().getFileName() + " " + result.getScanNum() + " " + result.getRetTime());
    //result[0] = keys[curIndex];
    IsotopeModel imodel = intensitySumModel(massArr, intArr, sModel.getSamIsoArr());
    result.addResult(imodel);
    //result[1] = tempArr[0];     //intensity
    //result[4] = samIsoArr.length;
    //result[7] = tempArr[1];      //Found iso num
    //result[10] = tempArr[2];	     //tolerance

    //   System.out.print("---\t" + keys[curIndex] + "\t");
    imodel = intensitySumModel(massArr, intArr, sModel.getRefIsoArr());
    result.addResult(imodel);
    imodel = intensitySumModel(massArr, intArr, sModel.getRef2IsoArr());
    result.addResult(imodel);
        /*
	if(null != refIsoArr) {
	    tempArr = intensitySum(massArr, intArr, refIsoArr);
	    result[2] = tempArr[0];
	    result[5] = refIsoArr.length;
	    result[8] = tempArr[1];
	    result[11] = tempArr[2];
	}

	if(null != ref2IsoArr) {
	    tempArr = intensitySum(massArr, intArr, ref2IsoArr);
	    result[3] = tempArr[0];
	    result[6] = ref2IsoArr.length;
	    result[9] = tempArr[1];
	    result[12] = tempArr[2];
	}*/

    return result;
  }

  //high resolution
  public static void readFullSpectrumMzXml(
    int[] keys,
    int curIndex,
    MzxmlSpectrumReader mzReader
  ) throws IOException, Exception {
    double result[] = new double[3];

    MzxmlPeakList pList = mzReader.scanNum2PeakList(keys[curIndex], false);

    MZXmlHandler.NativePeakList nPeakList = MZXmlHandler.decode32ToArr(pList.getEncodedM2zAndIntensities());

    intArr = nPeakList.getIntArr();
    massArr = nPeakList.getMassArr();

  }

  //Low resolution
  public static double[] readFullSpectrum(
    int[] keys,
    int curIndex,
    TIntLongHashMap index,
//           RandomAccessFile file,
    Object reader,
    double samStartMass, double samEndMass,
    double refStartMass, double refEndMass
  ) throws IOException, Exception {
    if (reader instanceof RandomAccessFile) {
      RandomAccessFile rFile = (RandomAccessFile) reader;
      return readFullSpectrumMS(keys, curIndex, index, rFile, samStartMass, samEndMass, refStartMass, refEndMass);
    } else if (reader instanceof MzxmlSpectrumReader) {
      MzxmlSpectrumReader mzReader = (MzxmlSpectrumReader) reader;

      return readFullSpectrumMzXml(keys, curIndex, mzReader, samStartMass, samEndMass, refStartMass, refEndMass);
    }

    return null;
  }

  //Low resolution
  public static double[] readFullSpectrumMzXml(
    int[] keys,
    int curIndex,
    MzxmlSpectrumReader mzReader,
    double samStartMass, double samEndMass,
    double refStartMass, double refEndMass
  ) throws IOException, Exception {

    double result[] = new double[3];

    if (curIndex < 0)
      return result;

    MzxmlPeakList pList = mzReader.scanNum2PeakList(keys[curIndex], false);

    MZXmlHandler.NativePeakList nPeakList = MZXmlHandler.decode32ToArr(pList.getEncodedM2zAndIntensities());

    double[] intArr = nPeakList.getIntArr();
    double[] massArr = nPeakList.getMassArr();

    result[0] = keys[curIndex];
    result[1] = intensitySum(massArr, intArr, samStartMass, samEndMass);
    result[2] = intensitySum(massArr, intArr, refStartMass, refEndMass);

    return result;
  }

  //Low resolution
  public static double[] readFullSpectrumMS(
    int[] keys,
    int curIndex,
    TIntLongHashMap index,
    RandomAccessFile file,
    double samStartMass, double samEndMass,
    double refStartMass, double refEndMass
  ) throws IOException {

//	RandomAccessFile file = (RandomAccessFile)ofile;

    double result[] = new double[3];

    if (curIndex < 0)
      return result;

    long startPos = index.get(keys[curIndex]);
    long endPos;

    file.seek(startPos);

    if ((curIndex + 1) >= keys.length)
      endPos = file.length();
    else
      endPos = index.get(keys[curIndex + 1]);

    int byteSize = (int) (endPos - startPos);
    byte[] bytes = new byte[byteSize];

    file.readFully(bytes);

    char ch;
    int pos = 0;

    ch = (char) bytes[pos];

    //Remove Z, S, I, D lines
    while ((ch = (char) bytes[pos]) == 'Z' || ch == 'S' || ch == 'I' || ch == 'D') {
      while (ch != CARRIAGE_RETURN) {
        pos++;
        ch = (char) bytes[pos];
      }

      pos++;
    }

    int arrSize = 0;
    for (int j = pos; j < byteSize; j++) {
      if (CARRIAGE_RETURN == (char) bytes[j])
        arrSize++;
    }


    buildSpec(arrSize, byteSize, bytes, pos);


    result[0] = keys[curIndex];
    result[1] = intensitySum(massArr, intArr, samStartMass, samEndMass);
    result[2] = intensitySum(massArr, intArr, refStartMass, refEndMass);

    if (conf.isUseProline()) {
      double[] addArr = conf.getAddtionalRefMassArr();
      if (null != addArr) {
        for (int i = 0; i < addArr.length; i += 2) {
          result[2] += intensitySum(massArr, intArr, addArr[i], addArr[i + 1]);

          //System.out.println(i + " " + (i+1) + " " + addArr[i] +  " " + addArr[i+1]);
        }
      }
    }

    return result;

  }

  private static void buildSpec(int arrSize, int byteSize, byte[] bytes, int pos) {

    massArr = new double[arrSize];
    intArr = new double[arrSize];
    int massIndex = 0;
    char ch;

    StringBuilder mass = new StringBuilder(10);
    StringBuilder intensity = new StringBuilder(10);
    intensity.append('0');

    boolean isMass = true;
    boolean isInt = true;

    int spaceCount = 0;

    if (hasChargeCol) {
      for (int i = pos; i < byteSize; i++) {
        ch = (char) bytes[i];

        //fixed on new ms2 files
        switch (ch) {
          case WINDOW_CR:
            break;

          case SPACE:
            spaceCount++;
            isMass = false;
            if (spaceCount == 1) isInt = true;

            break;

          case CARRIAGE_RETURN:
            isMass = true;
            isInt = false;
            spaceCount = 0;

            //try {
            intArr[massIndex] = Long.parseLong(intensity.toString());
            massArr[massIndex++] = Double.parseDouble(mass.toString());
            mass.delete(0, mass.length());  //This is faster than creating new StringBuilder object
            intensity.delete(0, intensity.length()).append('0');
            //} catch(Exception e) { System.out.println("-->>" + mass.toString() + " " + intensity.toString()); System.exit(0); };

            break;

          case DOT:
            isInt = false;

          default:
            if (isMass)
              mass.append(ch);
            else if (isInt) //remove decimal value of intensity
              intensity.append(ch);

            break;
        }
/*
                switch(ch)
                {
                    case WINDOW_CR:
                        break;

                    case SPACE:
                        spaceCount++;
                        isMass=false;
                        isInt=true;
                        break;

                    case CARRIAGE_RETURN:
                        spaceCount=0;
                        isMass=true;
                        isInt=true;

                        intArr[massIndex] = Long.parseLong(intensity.toString());
                        massArr[massIndex++] = Double.parseDouble(mass.toString());

                        mass.delete(0, mass.length());  //This is faster than creating new StringBuilder object
                        intensity.delete(0, intensity.length()).append('0');

                        break;

                    case DOT:
                        isInt=false;

                    default:
                        if(spaceCount>=2)
                            break;

                        if(isMass)
                            mass.append(ch);
                        else if(isInt) //remove decimal value of intensity
                            intensity.append(ch);

                        break;
                }
*/
      }
    } else {
      for (int i = pos; i < byteSize; i++) {
        ch = (char) bytes[i];

        //fixed on new ms2 files
        switch (ch) {
          case WINDOW_CR:
            break;

          case SPACE:
            spaceCount++;
            isMass = false;
            if (spaceCount == 1) isInt = true;

            break;

          case CARRIAGE_RETURN:
            isMass = true;
            isInt = false;
            spaceCount = 0;

            //try {
            intArr[massIndex] = Long.parseLong(intensity.toString());
            massArr[massIndex++] = Double.parseDouble(mass.toString());
            mass.delete(0, mass.length());  //This is faster than creating new StringBuilder object
            intensity.delete(0, intensity.length()).append('0');
            //} catch(Exception e) { System.out.println("-->>" + mass.toString() + " " + intensity.toString()); System.exit(0); };

            break;

          case DOT:
            isInt = false;

          default:
            if (isMass)
              mass.append(ch);
            else if (isInt) //remove decimal value of intensity
              intensity.append(ch);

            break;
        }
/*
                switch(ch)
                {
                    case WINDOW_CR:
                        break;

                    case SPACE:
                        isMass=false;
                        isInt=true;
                        break;

                    case CARRIAGE_RETURN:
                        isMass=true;
                        isInt=true;

                        intArr[massIndex] = Long.parseLong(intensity.toString());
                        massArr[massIndex++] = Double.parseDouble(mass.toString());

                        mass.delete(0, mass.length());  //This is faster than creating new StringBuilder object
                        intensity.delete(0, intensity.length()).append('0');

                        break;

                    case DOT:
                        isInt=false;

                    default:
                        if(isMass)
                            mass.append(ch);
                        else if(isInt) //remove decimal value of intensity
                            intensity.append(ch);

                        break;
                }
*/
      }
    }
  }

  public static double intensitySum(double[] massArr, double[] intArr, double startMass, double endMass) {
    double sumIntensity = 0;


    //System.out.println("b-->>" + bioSample[k][1] + " " + bioSample[k][2]);
    //System.out.println("y-->>" + yioSample[k][1] + " " + yioSample[k][2]);


    //Sum sample intensity
    int start = Arrays.binarySearch(massArr, startMass);
    int end = Arrays.binarySearch(massArr, endMass);

    //System.out.println(startMass + " " + endMass);
    if (start < 0)
      start = -start - 1;

    if (end < 0)
      end = -end - 2;

    for (int i = start; i <= end; i++) {
      sumIntensity += intArr[i];
      //  System.out.println(intArr[i]);
    }

    //System.out.println("sum---" + sumIntensity);

    return sumIntensity;
  }

  //return arr of intensity from spectrum
  public static double[] intensityArr(double[] massArr, double[] intArr, double[] isoArr, double massTolerance) {
    double[] resultArr = new double[isoArr.length];

    int foundIsoNum = 0;
    int totalPeakFound = 0;
    double toleranceSum = 0;

    int isoIndex = 0;

    for (int i = 0; i < isoArr.length; i++) {
      double tempTolerance = isoArr[i] / 1000 * massTolerance * 0.001;
//            double tempTolerance = massTolerance;

      int start = Arrays.binarySearch(massArr, isoArr[i] - tempTolerance);
      if (start < 0)
        start = -start - 1;

      int j = 0;
      double small = 100;
      double massSmall = -1;
      double highinten = 0;

      boolean isFound = false;

      while (true) {
        if (start >= massArr.length)
          break;

        double temp = isoArr[i] - massArr[start];
        if (temp < 0)
          temp = -temp;

        if (temp < tempTolerance) {
          totalPeakFound++;
          toleranceSum += temp;
          //	    System.out.println("==" + resultArr[i] +"\t" + isoArr[i] + "\t" + massArr[start] + "\t" + intArr[start] +"\t" + tempTolerance);

          resultArr[i] += intArr[start];
          //	    System.out.println("--" + resultArr[i] +"\t" + isoArr[i] + "\t" + massArr[start] + "\t" + intArr[start]);

          isFound = true;

          double diff = (massArr[start] - isoArr[i]);
          //if(diff<0)
          //diff = -diff;
          if (Math.abs(diff) < Math.abs(small)) {
            small = massArr[start] - isoArr[i];
            massSmall = massArr[start];
            highinten = intArr[start];
          }

        }

        if (massArr[start] > isoArr[i])
          break;

        start++;
      }

      if (isFound)
        foundIsoNum++;
    }

/*
	for(int i1=0;i1<isoArr.length;i1++)
	    System.out.println("==" + isoArr[i1] + "\t" + massTolerance);
	for(int i1=0;i1<resultArr.length;i1++)
	    System.out.println("==" + resultArr[i1]);
	for(int i1=0;i1<massArr.length;i1++)
	    System.out.println(massArr[i1] + "\t" + intArr[i1]);

System.exit(0);
*/
    return resultArr;
  }

  public static double[] intensitySum(double[] massArr, double[] intArr, double[] isoArr) {

    if (null == massArr || massArr.length <= 0) return null;

    double[] resultArr = new double[3];

    double sumIntensity = 0;
    int foundIsoNum = 0;
    int totalPeakFound = 0;
    double toleranceSum = 0;

    for (int i = 0; i < isoArr.length; i++) {
      //double tempTolerance = isoArr[i]/1000*massTolerance;
      //double tempTolerance = isoArr[i]/1000*massTolerance*0.001;
      double tempTolerance = isoArr[i] / 1000 * massTolerance;
      int start = Arrays.binarySearch(massArr, isoArr[i] - tempTolerance);
      if (start < 0)
        start = -start - 1;

      int j = 0;
      double small = 100;
      double massSmall = -1;
      double highinten = 0;

      boolean isFound = false;

      while (true) {
        if (start >= massArr.length)
          break;

        double temp = isoArr[i] - massArr[start];
        if (temp < 0)
          temp = -temp;

        if (temp < tempTolerance) {
          sumIntensity += intArr[start];
          totalPeakFound++;
          toleranceSum += temp;

          isFound = true;

          double diff = (massArr[start] - isoArr[i]);
          //if(diff<0)
          //diff = -diff;
          if (Math.abs(diff) < Math.abs(small)) {
            small = massArr[start] - isoArr[i];
            massSmall = massArr[start];
            highinten = intArr[start];
          }

        }

        if (massArr[start] > isoArr[i])
          break;

        start++;
      }

      if (isFound)
        foundIsoNum++;
    }

    resultArr[0] = sumIntensity;
    resultArr[1] = foundIsoNum;
    resultArr[2] = (totalPeakFound > 0) ? (toleranceSum / totalPeakFound) : (-1);

    return resultArr;
  }

  public static IsotopeModel intensitySumModel(double[] massArr, double[] intArr, double[] isoArr) {

    //double[] resultArr = new double[3];

    double sumIntensity = 0;
    int foundIsoNum = 0;
    int totalPeakFound = 0;
    double toleranceSum = 0;
    double[] isoFoundArr = new double[isoArr.length];

    for (int i = 0; i < isoArr.length; i++) {
      //double tempTolerance = isoArr[i]/1000*massTolerance;
      //double tempTolerance = isoArr[i]/1000*massTolerance*0.001;
      double tempTolerance = isoArr[i] / 1000 * massTolerance;
      int start = Arrays.binarySearch(massArr, isoArr[i] - tempTolerance);
      if (start < 0)
        start = -start - 1;

      int j = 0;
      double small = 100;
      double massSmall = -1;
      double highinten = 0;

      boolean isFound = false;

      while (true) {
        if (start >= massArr.length)
          break;

        double temp = isoArr[i] - massArr[start];
        if (temp < 0)
          temp = -temp;

        if (temp < tempTolerance) {
          sumIntensity += intArr[start];
          isoFoundArr[i] += intArr[start];

          totalPeakFound++;
          toleranceSum += temp;

          isFound = true;

          double diff = (massArr[start] - isoArr[i]);
          //if(diff<0)
          //diff = -diff;
          if (Math.abs(diff) < Math.abs(small)) {
            small = massArr[start] - isoArr[i];
            massSmall = massArr[start];
            highinten = intArr[start];
          }

        }

        if (massArr[start] > isoArr[i])
          break;

        start++;
      }

      if (isFound)
        foundIsoNum++;
    }


    IsotopeModel imodel = new IsotopeModel();
    imodel.setSumIntensity(sumIntensity);
    imodel.setFoundIsoNum(foundIsoNum);
    imodel.setIsoArr(isoArr);
    imodel.setIsoFoundArr(isoFoundArr);
    //resultArr[2] = (totalPeakFound>0)?(toleranceSum/totalPeakFound):(-1);


    return imodel;
  }


  public static void getPeaks(ChroPeptide peptide, int startWindow, int endWindow, int threshold) {
    //Point point= null;
    double prev = peptide.getData(startWindow).getSampleIntensity();
    double current;

    double peakValue = peptide.getData(startWindow).getSampleIntensity();
    int peakIndex = startWindow;
    double prevPeakValue;
    int prevPeakIndex;
    int chroCenter = peptide.getChroCenter();

    for (int i = startWindow; i <= endWindow; i++) {
      ChroData data = peptide.getData(i);

      if (data.getSampleIntensity() < threshold)
        continue;

      current = data.getSampleIntensity();

      if (current > peakValue && (Math.abs(peakIndex - chroCenter) > Math.abs(i - chroCenter))) {
        peakValue = current;
        peakIndex = i;
      }
    }
  }

  public static FragIonList getBestFragIons(long[][] bsTempArr, long[][] ysTempArr, long[][] brTempArr, long[][] yrTempArr, int startPeakIndex, int endPeakIndex, int maxShift) {
    PostOptions options = PostOptions.getInstance();

    double[] bsRegArr = getRegressArr(bsTempArr, startPeakIndex, endPeakIndex, maxShift);
    double[] ysRegArr = getRegressArr(ysTempArr, startPeakIndex, endPeakIndex, maxShift);
    double[] brRegArr = getRegressArr(brTempArr, startPeakIndex, endPeakIndex, maxShift);
    double[] yrRegArr = getRegressArr(yrTempArr, startPeakIndex, endPeakIndex, maxShift);

    List<FragIon> iList = new Vector<FragIon>();

    int index = 0;
    boolean isBion = true;

    for (int i = 0; i < bsRegArr.length; i++) {
      iList.add(new FragIon(i + 1, bsTempArr[i], brTempArr[i], bsRegArr[i] + brRegArr[i], true)); //bion sample + ref
      iList.add(new FragIon(bsRegArr.length - i, ysTempArr[i], yrTempArr[i], ysRegArr[i] + yrRegArr[i], false)); //yion sample + ref
    }


    Object[] fArr = iList.toArray();
    Arrays.sort(fArr);

    iList.clear();

    iList.add((FragIon) fArr[fArr.length - 1]);

    long[] samArr = new long[iList.get(0).getSArr().length];
    long[] refArr = new long[samArr.length];


    int bestIndex = 1;
    double bestReg = 0;
    int tempIndex = 0;
    FragIonList fList = new FragIonList();

    for (int i = fArr.length - 1; i >= 0; i--) {
      FragIon eachIon = ((FragIon) fArr[i]);
      fList.add(eachIon);

      long[] tempSArr = eachIon.getSArr();
      long[] tempRArr = eachIon.getRArr();

      for (int j = 0; j < tempSArr.length; j++) {
        samArr[j] += tempSArr[j];
        refArr[j] += tempRArr[j];
      }

      LinearRegression reg = new LinearRegression(samArr, refArr, startPeakIndex, (endPeakIndex != 0) ? endPeakIndex : samArr.length, maxShift);

      if (bestReg < reg.getCorr()) {
        bestReg = reg.getCorr();
        bestIndex = tempIndex;
      }

      tempIndex++;
    }

    fList.setBestIndex(bestIndex);

    return fList;

  }

  private static double[] getRegressArr(long[][] arr, int startPeakIndex, int endPeakIndex, int maxShift) {
    double[] regArr = new double[arr.length];

    for (int i = 0; i < arr.length; i++) {
      double sum = 0;
      for (int j = 0; j < arr.length; j++) {
        if (j == i)
          continue;

        LinearRegression reg = new LinearRegression(arr[i], arr[j], startPeakIndex, endPeakIndex, maxShift);

        if (reg.getCorr() > 0)
          sum += reg.getCorr();
      }

      regArr[i] = sum;
    }

    return regArr;
  }


  public static double calculateSNRatio(long[] samArr, long[] refArr, int startIndex, int endIndex) {
        /*
         long[] snData = new long[samArr.length * 2]; // data for signal to noise calculation

        for(int i=0;i<samArr.length;i++)
        {
            snData[i] = samArr[i];
            snData[samArr.length + i] = refArr[i];
        }

        Arrays.sort(snData);
         **/


    startIndex = (startIndex >= 0) ? startIndex : 0;
    endIndex = (endIndex < samArr.length) ? endIndex : samArr.length - 1;

    long noise = 0;
    long signal = 0;

    for (int i = startIndex; i <= endIndex; i++) {
      signal += samArr[i];
      signal += refArr[i];
    }

    int peakWidth = (endIndex - startIndex + 1);
    signal = signal / peakWidth;

    //Arrays.sort(samArr);
    //Arrays.sort(refArr);

    //int sampleSize = (int)(snData.length*0.05);


    int noiseStart = startIndex - peakWidth;
    noiseStart = (noiseStart >= 0) ? noiseStart : 0;

    int noiseEnd = endIndex + peakWidth;
    noiseEnd = (noiseEnd >= samArr.length) ? (samArr.length - 1) : noiseEnd;

    long[] tempSamData = new long[noiseEnd - noiseStart + 1];
    long[] tempRefData = new long[tempSamData.length];

    int index = 0;
    for (int i = noiseStart; i <= noiseEnd; i++) {
      tempSamData[index] = samArr[i];
      tempRefData[index] = refArr[i];
      index++;
    }

    Arrays.sort(tempSamData);
    Arrays.sort(tempRefData);

    int sampleSize = (int) (tempSamData.length * 0.3);

    for (int i = 0; i < sampleSize; i++) {
      //noise += snData[i];
      //signal += snData[snData.length-i-1];

      noise += tempSamData[i];
      noise += tempRefData[i];

      //    signal += samArr[samArr.length-i-1];
      //    signal += refArr[refArr.length-i-1];

    }


    noise = noise / (sampleSize * 2);
    //signal = signal / (sampleSize*2);

    double snRatio;
    if (noise == 0)
      snRatio = Math.log10(signal);
    else
      snRatio = Math.log10(signal / noise);

    return snRatio;

  }

  public static class SpectrumModel {
    private boolean highRes;

    private IndexedFile iFile;
    private double samStartMass;
    private double samEndMass;
    private double refStartMass;
    private double refEndMass;

    private double[] samIsoArr;
    private double[] refIsoArr;
    private double[] ref2IsoArr;

    private double[][] bioSample;
    private double[][] yioSample;
    private double[][] bioRef;
    private double[][] yioRef;

    private int[] keys;
    private int diff;

    private TIntLongHashMap index;
    private RandomAccessFile file;

    private int steepArea;
    private String result;
    private List msmsSpecificMassList;

    private MzxmlSpectrumReader mzxmlreader;

    public int getDiff() {
      return this.diff;
    }

    public void setDiff(int diff) {
      this.diff = diff;
    }

    public Object getGenericIndexFile(Configuration conf) {
      if (conf.getSpectrumFormat() == conf.MS_FILE_FORMAT)
        return file;
      else if (conf.getSpectrumFormat() == conf.MZXML_FILE_FORMAT)
        return mzxmlreader;
      else
        return null;

    }

    public IndexedFile getIFile() {
      return iFile;
    }

    public void setIFile(IndexedFile iFile) {
      this.iFile = iFile;
    }

    public int[] getKeys() {
      return keys;
    }

    public void setKeys(int[] keys) {
      this.keys = keys;
    }

    public double getSamStartMass() {
      return samStartMass;
    }

    public void setSamStartMass(double samStartMass) {
      this.samStartMass = samStartMass;
    }

    public double getSamEndMass() {
      return samEndMass;
    }

    public void setSamEndMass(double samEndMass) {
      this.samEndMass = samEndMass;
    }

    public double getRefStartMass() {
      return refStartMass;
    }

    public void setRefStartMass(double refStartMass) {
      this.refStartMass = refStartMass;
    }

    public double getRefEndMass() {
      return refEndMass;
    }

    public void setRefEndMass(double refEndMass) {
      this.refEndMass = refEndMass;
    }

    public boolean isHighRes() {
      return highRes;
    }

    public void setHighRes(boolean highRes) {
      this.highRes = highRes;
    }

    public TIntLongHashMap getIndex() {
      return index;
    }

    public void setIndex(TIntLongHashMap index) {
      this.index = index;
    }

    public RandomAccessFile getFile() {
      return file;
    }

    public void setFile(RandomAccessFile file) {
      this.file = file;
    }

    public int getSteepArea() {
      return steepArea;
    }

    public void setSteepArea(int steepArea) {
      this.steepArea = steepArea;
    }

    public double[] getSamIsoArr() {
      return samIsoArr;
    }

    public void setSamIsoArr(double[] samIsoArr) {
      this.samIsoArr = samIsoArr;
    }

    public double[] getRefIsoArr() {
      return refIsoArr;
    }

    public void setRefIsoArr(double[] refIsoArr) {
      this.refIsoArr = refIsoArr;
    }

    public double[] getRef2IsoArr() {
      return ref2IsoArr;
    }

    public void setRef2IsoArr(double[] ref2IsoArr) {
      this.ref2IsoArr = ref2IsoArr;
    }

    public String getResult() {
      return this.result;
    }

    public void setResult(String result) {
      this.result = result;
    }

    public double[][] getBioSample() {
      return bioSample;
    }

    public void setBioSample(double[][] bioSample) {
      this.bioSample = bioSample;
    }

    public double[][] getYioSample() {
      return yioSample;
    }

    public void setYioSample(double[][] yioSample) {
      this.yioSample = yioSample;
    }

    public double[][] getBioRef() {
      return bioRef;
    }

    public void setBioRef(double[][] bioRef) {
      this.bioRef = bioRef;
    }

    public double[][] getYioRef() {
      return yioRef;
    }

    public void setYioRef(double[][] yioRef) {
      this.yioRef = yioRef;
    }

    public List getMsmsSpecificMassList() {
      return msmsSpecificMassList;
    }

    public void setMsmsSpecificMassList(List msmsSpecificMassList) {
      this.msmsSpecificMassList = msmsSpecificMassList;
    }

    public MzxmlSpectrumReader getMzxmlreader() {
      return mzxmlreader;
    }

    public void setMzxmlreader(MzxmlSpectrumReader mzxmlreader) {
      this.mzxmlreader = mzxmlreader;
    }

  }

  //calculate weighted standard deviation for the labeled data
  public static double getWeightedStdev(double corr) {
    double rsqrtLog = Math.log(corr * corr);
    double stdevLog = -0.84 * rsqrtLog + 0.43;
    //double invStdev = Math.exp(stdevLog);

    return Math.exp(stdevLog);

  }


  public static void main(String[] args) throws Exception {
    //CalcUtil.calculateMass(new RandomAccessFile(args[0], "r"), Long.parseLong(args[1]), Long.parseLong(args[2]),
    //	Double.parseDouble(args[3]), Double.parseDouble(args[4]));
  }


  public static double[][] readSpectrumOnly(
    int i,
    IndexedFile iFile
  ) throws IOException {
    RandomAccessFile file = iFile.getFile();
    TIntLongHashMap index = iFile.getMsIndex();

    int[] keys = iFile.getKeys();

    int byteSize;
    byte[] bytes;
    char ch;
    int pos;

    if (i < 0)
      return null;

    long sampleIndex = index.get(keys[i]);

    StringBuilder mass, intensity;

    if (sampleIndex < 0)
      return null;


    file.seek(sampleIndex);

    int increIndex = 1;

    while (true) {
      if ((i + increIndex) >= keys.length) {
        byteSize = (int) (file.length() - sampleIndex);

        break;
      }

      if (index.get(keys[i + increIndex]) > 0) {
        byteSize = (int) (index.get(keys[i + increIndex]) - sampleIndex);
        break;
      }

      increIndex++;
    }

    //byteSize = (int)(index.get(keys[i+1]) - sampleIndex);
    bytes = new byte[byteSize];

    file.readFully(bytes);
    pos = 0;

    ch = (char) bytes[pos];

    //        pos++;
    //Remove Z, S, I, D lines
    //System.out.println("sample");
    while ((ch = (char) bytes[pos]) == 'Z' || ch == 'S' || ch == 'I' || ch == 'D') {
      while (ch != CARRIAGE_RETURN) {
        pos++;
        ch = (char) bytes[pos];
      }

      pos++;
    }


    int arrSize = 0;
    for (int j = pos; j < byteSize; j++) {
      if (CARRIAGE_RETURN == (char) bytes[j])
        arrSize++;

    }
    double[] massArr = new double[arrSize];
    double[] intArr = new double[arrSize];

    //pos++;
    mass = new StringBuilder(10);
    intensity = new StringBuilder(10);
    intensity.append('0');

    boolean isMass = true;
    boolean isInt = true;

    int massIndex = 0;
    //double tempMass;
    int spaceCount = 0;

    for (int j = pos; j < byteSize; j++) {
      ch = (char) bytes[j];
      //System.out.print(ch);
      //fixed on new ms2 files
      switch (ch) {
        case WINDOW_CR:
          break;

        case SPACE:
          spaceCount++;
          isMass = false;
          if (spaceCount == 1) isInt = true;

          break;

        case CARRIAGE_RETURN:
          isMass = true;
          isInt = false;
          spaceCount = 0;

          //try {
          intArr[massIndex] = Long.parseLong(intensity.toString());
          massArr[massIndex++] = Double.parseDouble(mass.toString());
          mass.delete(0, mass.length());  //This is faster than creating new StringBuilder object
          intensity.delete(0, intensity.length()).append('0');
          //} catch(Exception e) { System.out.println("-->>" + mass.toString() + " " + intensity.toString()); System.exit(0); };

          break;

        case DOT:
          isInt = false;

        default:
          if (isMass)
            mass.append(ch);
          else if (isInt) //remove decimal value of intensity
            intensity.append(ch);

          break;
      }

/*
		switch(ch)
		{
		    case WINDOW_CR:
			break;

		    case SPACE:
			isMass=false;
			isInt=true;
			break;

		    case CARRIAGE_RETURN:
			isMass=true;
			isInt=false;

                        //try {
			intArr[massIndex] = Long.parseLong(intensity.toString());
			massArr[massIndex++] = Double.parseDouble(mass.toString());
			mass.delete(0, mass.length());  //This is faster than creating new StringBuilder object
			intensity.delete(0, intensity.length()).append('0');
                        //} catch(Exception e) { System.out.println("-->>" + mass.toString() + " " + intensity.toString()); System.exit(0); };

			break;

		    case DOT:
			isInt=false;

		    default:
			if(isMass)
			    mass.append(ch);
			else if(isInt) //remove decimal value of intensity
			    intensity.append(ch);

                    break;
            }
*/
    }

    double[][] massIntArr = new double[2][];
    massIntArr[0] = massArr;
    massIntArr[1] = intArr;


//	for(int i1=0;i1<massArr.length;i1++)
//	    System.out.println(massArr[i1] + "\t" + intArr[i1]);

//	System.exit(0);

    return massIntArr;
  }

  public static double[] getMassArr() {
    return massArr;
  }

  public static double[] getIntArr() {
    return intArr;
  }

    /*
	public static long[] smooth(long[] input) {
		//long[] input = {3, 4, 5, 2, 3, 4, 5, 6, 7, 4, };
		int size=3;
		long[] out = new long[input.length-2];

		for(int i=0;i<input.length-size+1;i++) {
			long sum = input[i+0];
			sum += input[i+1];
			sum += input[i+2];

			long average = sum / size;

			out[i] = average;
		}

		return out;

	}*/


  public static int getScanToSearch() {
    return scanToSearch;
  }

  public static void setScanToSearch(int scanToSearch) {
    CalcUtil.scanToSearch = scanToSearch;
  }

  public static double[] getXyMassArr() {
    return xyMassArr;
  }

  public static void setXyMassArr(double[] xyMassArr) {
    CalcUtil.xyMassArr = Arrays.copyOf(xyMassArr, xyMassArr.length);
  }

  public static double[] getXyIntArr() {
    return xyIntArr;
  }

  public static void setXyIntArr(double[] xyIntArr) {
    CalcUtil.xyIntArr = Arrays.copyOf(xyIntArr, xyIntArr.length);
  }

  public static int getLowbound() {
    return lowbound;
  }

  public static void setLowbound(int lowbound) {
    CalcUtil.lowbound = lowbound;
  }

  public static int getHighbound() {
    return highbound;
  }

  public static void setHighbound(int highbound) {
    CalcUtil.highbound = highbound;
  }
}
