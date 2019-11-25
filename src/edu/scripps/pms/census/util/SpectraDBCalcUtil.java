package edu.scripps.pms.census.util;

import edu.scripps.pms.census.conf.Configuration;
import edu.scripps.pms.census.exception.CensusGeneralException;
import edu.scripps.pms.census.exception.CensusIndexOutOfBoundException;
import edu.scripps.pms.census.exception.PrecursorNotFoundException;
import edu.scripps.pms.census.hash.IndexedFile;
import edu.scripps.pms.census.model.ReportIon;
import edu.scripps.pms.census.model.SpecRange;

import edu.scripps.pms.util.PmsUtil;
import edu.scripps.pms.util.sqlite.spectra.SpectraDB;
import gnu.trove.TDoubleArrayList;
import gnu.trove.TIntDoubleHashMap;
import gnu.trove.TIntLongHashMap;

import java.io.*;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

import static edu.scripps.pms.census.util.CalcUtil.*;


public class SpectraDBCalcUtil {

    private static TIntDoubleHashMap precursorMap;


    public static double[] readSpectrum(
            int[] keys,
            int i,
            int refIndexKey,
            SpectraDB spectraDB,
            double[][] bioSample, double[][] bioRef,
            double[][] yioSample, double[][] yioRef
    ) throws IOException, SQLException {

        int ionLength = bioSample.length;
        double[] result = new double[ionLength * 4 + ION_START_INDEX];


        result[0] = keys[i];

        /**********************************
         * Go to Spectrum for sample, and calculate intensity
         **********************************/
        SpectraDB.Spectrum spectrum = spectraDB.getSpectrumFromDB(keys[0]);


        double[] massArr = spectrum.getMzList().toNativeArray();
        double[] intArr = spectrum.getIntensityList().toNativeArray();

        int sampleSum = 0;


        for (int k = 0; k < bioSample.length; k++) {

            result[ION_START_INDEX + k] = intensitySum(massArr, intArr, bioSample[k][1], bioSample[k][2]);
            result[ION_START_INDEX + ionLength + k] = intensitySum(massArr, intArr, yioSample[k][1], yioSample[k][2]);

            sampleSum += result[ION_START_INDEX + k];
            sampleSum += result[ION_START_INDEX + bioSample.length + k];

        }

        result[1] = sampleSum;

        if (refIndexKey < 0)
            return result;

        int scan= keys[refIndexKey];


        SpectraDB.Spectrum spectrum1 = spectraDB.getSpectrumFromDB(scan);
        massArr = spectrum1.getMzList().toNativeArray();
        intArr = spectrum1.getIntensityList().toNativeArray();

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


    public static String calculateMS2Mass(
            IndexedFile iFile,
            SpecRange range,
            int keyIndex,
            ///int diff,
            double[][] bioSample, double[][] bioRef,
            double[][] yioSample, double[][] yioRef,
            //double sPrecursor, double rPrecursor,
            Configuration conf, int cState, SpectraDB spectraDB
    ) throws PrecursorNotFoundException, IOException, CensusIndexOutOfBoundException, Exception {
//System.out.println("==-");


        precursorMap = iFile.getPrecursorMap();

        int[] keys = iFile.getKeys();

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

            CalcUtil.SpectrumModel sModel = new CalcUtil.SpectrumModel();
            sModel.setHighRes(false);
            sModel.setKeys(keys);
            sModel.setIndex(null);
            sModel.setIFile(iFile);
            sModel.setFile(null);
            sModel.setBioSample(bioSample);
            sModel.setYioSample(yioSample);
            sModel.setBioRef(bioRef);
            sModel.setYioRef(yioRef);
            //sModel.setDiff(diff);

            return peakFindingMSMS(sModel, range, conf, keyIndex,spectraDB);
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
                        spectraDB,
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


    public static String peakFindingMSMS(CalcUtil.SpectrumModel sModel, SpecRange range, Configuration conf,
                                         int keyIndex, SpectraDB spectraDB)
            throws IOException, SQLException {

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
                int heavyKey = findHeavyKey(keys, moveLeftKeyIndex, iFile, conf, conf.getCalcRefAvgMass());
                //result[leftIndex] = readSpectrum(keys, moveLeftKeyIndex, index, sModel.getDiff(), file, sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());
                result[leftIndex] = readSpectrum(keys, moveLeftKeyIndex,  heavyKey, spectraDB, sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());
            } else //data dependent
            {
                if (sModel.isHighRes())
                    result[leftIndex] = readFullSpectrum(keys, moveLeftKeyIndex, spectraDB, sModel.getSamIsoArr(), sModel.getRefIsoArr());
                else
                    result[leftIndex] = readFullSpectrum(keys, moveLeftKeyIndex, spectraDB, sModel.getSamStartMass(), sModel.getSamEndMass(), sModel.getRefStartMass(), sModel.getRefEndMass());

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
                int heavyKey = findHeavyKey(keys, moveRightKeyIndex, iFile, conf, conf.getCalcRefAvgMass());
                result[rightIndex] = readSpectrum(keys, moveRightKeyIndex, heavyKey , spectraDB,
                        sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());
            } else //data dependent
            {
                if (sModel.isHighRes())
                    result[rightIndex] = readFullSpectrum(keys, moveRightKeyIndex, spectraDB, sModel.getSamIsoArr(), sModel.getRefIsoArr());
                else
                    result[rightIndex] = readFullSpectrum(keys, moveRightKeyIndex, spectraDB, sModel.getSamStartMass(), sModel.getSamEndMass(), sModel.getRefStartMass(), sModel.getRefEndMass());
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
                        {
                            int heavyKey = findHeavyKey(keys, tempKeyIndex, iFile, conf, conf.getCalcRefAvgMass());
                            arr[steepCount] = readSpectrum(keys, tempKeyIndex, heavyKey, spectraDB,
                                    sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());
                        }
                        else {
                            if (sModel.isHighRes())
                                arr[steepCount] = readFullSpectrum(keys, tempKeyIndex,spectraDB, sModel.getSamIsoArr(), sModel.getRefIsoArr());
                            else
                                arr[steepCount] = readFullSpectrum(keys, tempKeyIndex, spectraDB, sModel.getSamStartMass(), sModel.getSamEndMass(), sModel.getRefStartMass(), sModel.getRefEndMass());
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
                {
                    int heavyKey = findHeavyKey(keys, tempKeyIndex, iFile, conf, conf.getCalcRefAvgMass());
                    arr[steepCount] = readSpectrum(keys, tempKeyIndex, heavyKey, spectraDB, sModel.getBioSample(),
                            sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());

                }
                else {
                    if (sModel.isHighRes())
                        arr[steepCount] = readFullSpectrum(keys, tempKeyIndex, spectraDB, sModel.getSamIsoArr(), sModel.getRefIsoArr());
                    else
                        arr[steepCount] = readFullSpectrum(keys, tempKeyIndex, spectraDB, sModel.getSamStartMass(), sModel.getSamEndMass(), sModel.getRefStartMass(), sModel.getRefEndMass());
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
                        {
                            int heavyKey = findHeavyKey(keys, tempKeyIndex, iFile, conf, conf.getCalcRefAvgMass());
                            arr[steepCount] = readSpectrum(keys, tempKeyIndex, heavyKey, spectraDB,
                                    sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());
                        }
                        else if (sModel.isHighRes())
                            arr[steepCount] = readFullSpectrum(keys, tempKeyIndex, spectraDB, sModel.getSamIsoArr(), sModel.getRefIsoArr());
                        else
                            arr[steepCount] = readFullSpectrum(keys, tempKeyIndex, spectraDB, sModel.getSamStartMass(), sModel.getSamEndMass(), sModel.getRefStartMass(), sModel.getRefEndMass());
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
                {
                    int heavyKey = findHeavyKey(keys, tempKeyIndex, iFile, conf, conf.getCalcRefAvgMass());
                    arr[steepCount] = readSpectrum(keys, tempKeyIndex, heavyKey, spectraDB, sModel.getBioSample(),
                            sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());
                }
                else if (sModel.isHighRes())
                    arr[steepCount] = readFullSpectrum(keys, tempKeyIndex, spectraDB, sModel.getSamIsoArr(), sModel.getRefIsoArr());
                else
                    arr[steepCount] = readFullSpectrum(keys, tempKeyIndex, spectraDB, sModel.getSamStartMass(), sModel.getSamEndMass(), sModel.getRefStartMass(), sModel.getRefEndMass());

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
            {
                int heavyKey = findHeavyKey(keys, moveLeftKeyIndex, iFile, conf, conf.getCalcRefAvgMass());
                result[leftIndex] = readSpectrum(keys, moveLeftKeyIndex, heavyKey, spectraDB, sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());

            }
            else if (sModel.isHighRes())
                result[leftIndex] = readFullSpectrum(keys, moveLeftKeyIndex, spectraDB, sModel.getSamIsoArr(), sModel.getRefIsoArr());
            else
                result[leftIndex] = readFullSpectrum(keys, moveLeftKeyIndex, spectraDB, sModel.getSamStartMass(), sModel.getSamEndMass(), sModel.getRefStartMass(), sModel.getRefEndMass());

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
            {
                int heavyKey = findHeavyKey(keys, moveRightKeyIndex, iFile, conf, conf.getCalcRefAvgMass());
                result[rightIndex] = readSpectrum(keys, moveRightKeyIndex, heavyKey, spectraDB, sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());
            }
            else if (sModel.isHighRes())
                result[rightIndex] = readFullSpectrum(keys, moveRightKeyIndex,spectraDB, sModel.getSamIsoArr(), sModel.getRefIsoArr());
            else
                result[rightIndex] = readFullSpectrum(keys, moveRightKeyIndex, spectraDB, sModel.getSamStartMass(), sModel.getSamEndMass(), sModel.getRefStartMass(), sModel.getRefEndMass());

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


    public static String readSpectrumSpecificPeaksWithPurityCorrection(
            int i,
            IndexedFile iFile,
            SpectraDB spectraDB,
            List massList,
            Configuration conf
    ) throws SQLException {



        int[] keys = iFile.getKeys();
        StringBuffer result = new StringBuffer();


        if (i < 0)
            return null;



        SpectraDB.Spectrum spectrum = spectraDB.getSpectrumFromDB(keys[i]);

        result.append(keys[i]).append(" ");


        double[] massArr = spectrum.getMzList().toNativeArray();
        double[] intArr = spectrum.getIntensityList().toNativeArray();



        int tmpIndex = 0;

        double[] correctIntArr = new double[massList.size()];
        double[] intSumArr = new double[massList.size()];

        //  List<MZValues> mzList = new ArrayList<>();
        for (Iterator<ReportIon> itr = massList.iterator(); itr.hasNext(); ) {
            //  MZValues mzValues = new MZValues();

            ReportIon ri = itr.next();
            double d = ri.getMass();

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

    public static double[] readFullSpectrum(
            int[] keys,
            int curIndex,
            SpectraDB spectraDB,
            double[] samIsoArr,
            double[] refIsoArr
    ) throws SQLException {

        double result[] = new double[9];
        if (curIndex < 0)
            return result;
        SpectraDB.Spectrum spectrum = spectraDB.getSpectrumFromDB(keys[curIndex]);
        massArr = spectrum.getMzList().toNativeArray();
        intArr = spectrum.getIntensityList().toNativeArray();

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

        return result;
    }

    public static double[] readFullSpectrum(
            int[] keys,
            int curIndex,
            SpectraDB spectraDB,
            double samStartMass, double samEndMass,
            double refStartMass, double refEndMass
    ) throws IOException, SQLException {
        return readFullSpectrumMS(keys, curIndex, spectraDB, samStartMass, samEndMass, refStartMass, refEndMass);

    }



    //Low resolution
    public static double[] readFullSpectrumMS(
            int[] keys,
            int curIndex,
            SpectraDB spectraDB,
            double samStartMass, double samEndMass,
            double refStartMass, double refEndMass
    ) throws IOException, SQLException {

//	RandomAccessFile file = (RandomAccessFile)ofile;

        double result[] = new double[3];

        if (curIndex < 0)
            return result;

        SpectraDB.Spectrum spectrum = spectraDB.getSpectrumFromDB(keys[curIndex]);
        massArr = spectrum.getMzList().toNativeArray();
        intArr = spectrum.getIntensityList().toNativeArray();


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

    public static void getSpectrumMS1Arr(SpectraDB spectraDB, int scanNum, double start_pos, double end_pos,
                                         int prevScan, TDoubleArrayList massArr, TDoubleArrayList intArr)
            throws SQLException {

        int  prevScanfromMS1=prevScan;
        SpectraDB.Spectrum spectrum = spectraDB.getSpectrumFromDB(prevScan);
        if(spectrum == null)
        {
            massArr = null;
            intArr = null;
            return;
        }
        TDoubleArrayList mzList = spectrum.getMzList();
        TDoubleArrayList intensityList = spectrum.getIntensityList();
        for(int i=0; i<mzList.size(); i++)
        {
            double mz = mzList.get(i);
            double intensity = intensityList.get(i);
            if(mz>start_pos && mz < end_pos)
            {
                massArr.add(mz);
                intArr.add(intensity);
            }
        }

        //return new double[][]{massArr,intArr};
        //return massInt;
    }
    //Low resolution



}
