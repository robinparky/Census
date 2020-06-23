/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.scripps.pms.census.labelFree.util;

import edu.scripps.pms.census.ChroGenerator;
import edu.scripps.pms.census.ElementComposition;
import edu.scripps.pms.census.conf.Configuration;
import edu.scripps.pms.census.exception.CensusIndexOutOfBoundException;
import edu.scripps.pms.census.exception.InvalidAAException;
import edu.scripps.pms.census.hash.IndexedFile;
import edu.scripps.pms.census.io.IsotopeReader;
//import edu.scripps.pms.census.labelFree.LabelfreeMissingPeptideBuilder;
import edu.scripps.pms.census.labelFree.LabelfreeMissingPeptideBuilderSplit;
import edu.scripps.pms.census.model.SpecRange;
import edu.scripps.pms.census.util.CalcUtil;

import static edu.scripps.pms.census.util.CalcUtil.intensitySum;
import static edu.scripps.pms.census.util.CalcUtil.peakFinding;
import static edu.scripps.pms.census.util.CalcUtil.readFullSpectrum;
import static edu.scripps.pms.census.util.CalcUtil.readFullSpectrumMzXml;
import static edu.scripps.pms.census.util.CalcUtil.readSpectrum;
import static edu.scripps.pms.census.util.TimsTOFXICDB.getPeakAreaEstimate;
import static rpark.statistics.GaussianFitting.getGaussianPeakRangeIndex;
import static rpark.statistics.model.GaussianPeakModel.getGaussianPeakArea;

import edu.scripps.pms.census.util.CalcUtilGeneric;
import edu.scripps.pms.census.util.IsotopeDist;
import edu.scripps.pms.census.util.TimsTOFXICDB;
import edu.scripps.pms.census.util.dtaselect.Peptide;
import edu.scripps.pms.census.util.io.MzxmlSpectrumReader;
import edu.scripps.pms.util.sqlite.spectra.SpectraDB;
import gnu.trove.TDoubleArrayList;
import gnu.trove.TIntDoubleHashMap;
import gnu.trove.TIntLongHashMap;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.sql.SQLException;
import java.util.*;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.analysis.function.Gaussian;
import org.jdom.Element;
import rpark.statistics.GaussianFitting;
import rpark.statistics.Smooth;
import rpark.statistics.model.GaussianPeakModel;

/**
 *
 * @author rpark
 * @author rohan
 */
public class LabelfreeChroUtil {


    private final static int MOVE_LEFT=1;
    private final static int MOVE_RIGHT=2;

    public static Element getPeptideDomElement(Peptide peptide, IsotopeReader isoReader,
                                               String eachPath,
                                               Hashtable<String, IndexedFile> origMs1FileHt,
                                               Hashtable<String, IndexedFile> splitMs1FileHt,
                                               HashMap<String, String> splitSpectraMap,
                                               HashMap<String, HashMap<Integer, Integer>> ms2ToMs1Map) throws IOException, Exception {
        String pepSequence = peptide.getSequence();
        String fileName=peptide.getFileName();
        fileName = ChroGenerator.cleanFileName(fileName);

        //trim additional characters from peptide sequence at both ends

        char[] ch = pepSequence.substring(2, peptide.getSequence().length()-2).toCharArray();


        ElementComposition element = new ElementComposition(ch, 0, ch.length, isoReader.getIsotope());

        try {
            element.calculate();
        }
        catch (InvalidAAException invE)
        {
            System.out.println("Not Quantifiable peptide : " + pepSequence);
            return null;
        }
        catch(Exception e) {
            System.out.println("Not Quantifiable peptide : " + pepSequence);
            return null;
        }

        if(!element.isQuantifiable())
        {
            System.out.print("\nError : ");
            System.out.println(pepSequence + " is not quantifiable.");
            return null;
        }

        Configuration conf = Configuration.getInstance();

        IsotopeDist sampleDist = null;

        sampleDist = new IsotopeDist(element.getElementSampleArr(), element.getModShift(), true);

        //String origFileKey = eachPath + fileName+".ms1";
        String origFile = eachPath + fileName+".ms1"; //splitSpectraMap.get(eachPath + fileName+".ms1");
        int ms2ScanNumber = peptide.getScanNumber();

      //  splitSpectraMap.get(fileKey)

        HashMap<Integer, Integer> ms2Ms1ScanMap = ms2ToMs1Map.get(origFile);

        int ms1ScanNum = ms2Ms1ScanMap.get(ms2ScanNumber);
        String fileKey = ms1ScanNum + "\t" + origFile;
        String splitMs1File = splitSpectraMap.get(fileKey);


        IndexedFile splitIFile = splitMs1FileHt.get(splitMs1File);
        IndexedFile origIFile = origMs1FileHt.get(origFile);


        if (splitIFile == null){
            System.out.println("The file is missing=== "+eachPath+fileName+".ms1");
            System.exit(0);
        }
        TIntDoubleHashMap retentionTimeMap = origIFile.getRetentionTimeMap();
        TIntDoubleHashMap ionInjectionMap = origIFile.getIonInjectionMap();


        double retTime = retentionTimeMap.get( ms2ScanNumber );
        double ionInjectionTime = ionInjectionMap.get(ms2ScanNumber);

        int tmpScanNumber = ms2ScanNumber;
        int tmpCount=0;
        while(retTime<=0) {
            retTime = retentionTimeMap.get(--tmpScanNumber);

	        if(tmpScanNumber<=0) break;

            tmpCount++;
            if(tmpCount>20000) { //check if ms1.index file contains retention time
                System.out.println(">>>" +ms2ScanNumber);
                System.out.println("==retention time is required in ms1.index file.");
                System.exit(0);
            }
        }

        ionInjectionTime = ionInjectionMap.get(tmpScanNumber);


       // System.out.println( retentionTimeMap.get( peptide.getScanNumber()) );
      //  System.out.println( retTime + " == " + tmpScanNumber );

        Element peptideEle=LabelfreeChroUtil.createXmlChroPeptideTitle(peptide); //true is for full scan
        peptideEle.setAttribute("rt", String.valueOf(retTime));
        peptideEle.setAttribute("iit", String.valueOf(ionInjectionTime));
        conf.setCalcSamAvgMass(sampleDist.getAvgMass());
        peptideEle.setAttribute("lightStartMass", String.valueOf(sampleDist.getStartMass()));
        peptideEle.setAttribute("lightAvgMass", String.valueOf(conf.getCalcSamAvgMass()));

        //String fileName = peptide.getFileName().substring(0, peptide.getFileName().indexOf(".")+1) + "ms1";

      //  SpecRange range = null;

        int tmpScanNum = Integer.parseInt(peptide.getScanNum());
    //    range = new SpecRange( tmpScanNum, tmpScanNum );

        int scanNum = Integer.parseInt(peptide.getScanNum());

        switch(conf.getSpectrumFormat())
        {
            case Configuration.MS_FILE_FORMAT:
                fileName += ".ms1";
                break;

            case Configuration.MZXML_FILE_FORMAT:
                fileName += ".mzXML";
                break;

            default:
                break;
        }



        //iFile = ht.get(fileName);
        if(null == splitIFile)
        {
            splitIFile = splitMs1FileHt.get(eachPath + fileName.substring(1, fileName.length()));
            if(null == splitIFile) {
                System.out.println("Error : cannot find the file " + eachPath + fileName.substring(1, fileName.length()));
                throw new IOException("Error : cannot find the file " + eachPath + fileName);
            }
            //System.exit(0);
        }

        int[] keys = origIFile.getKeys();
        int keyIndex=-1;
        keyIndex = Arrays.binarySearch(keys, scanNum);

        if(keyIndex<0) //Cannot find index
            keyIndex=-(++keyIndex); //Math.abs(++keyIndex);

        if(keyIndex>=keys.length)
            keyIndex--;

        if(keyIndex>0) keyIndex--;  //adjust

        Element chro = new Element("chro");

        int chargeState = Integer.parseInt(peptide.getChargeState());

        String chroText="";


        double[] samIsoArr = sampleDist.getHighMassList();

        for(int i=0;i<samIsoArr.length;i++)
        {
            samIsoArr[i] = (samIsoArr[i]+chargeState*ChroGenerator.PROTON_MASS)/chargeState;
        }

        GaussianPeakModel gModel = calculateLabelfreeFullMS(
                peptide,
                keyIndex,
                origIFile,
                samIsoArr,
                //range,
                origFile,
                splitMs1FileHt,
                splitSpectraMap
                );

        chroText = gModel.getChroData();

        //peak area
        double peakAreaTemp = gModel.getPeakArea();
        if(peakAreaTemp<=0) {




        int[] tempKeys = origIFile.getKeys();
        TIntLongHashMap indexTemp = origIFile.getMsIndex();
        RandomAccessFile randomAccessFile = origIFile.getFile();

        if(tempKeys.length<=keyIndex) {
            System.out.println("out of bound exception for split file");
            throw new CensusIndexOutOfBoundException();
        }

        long startPos = indexTemp.get(tempKeys[keyIndex]);
        long endPos;

        randomAccessFile.seek(startPos);

        if( (keyIndex+1)>=tempKeys.length )
            endPos = randomAccessFile.length();
        else
            endPos = indexTemp.get(tempKeys[keyIndex+1]);

        int diff = (int)(endPos-startPos);
            long backGroundNoise = CalcUtilGeneric.getBackGroundNoise(startPos, diff, origIFile);
            gModel.setPeakArea(backGroundNoise);
            gModel.setPeakDetected(false);

        }

        chro.setText( chroText );

        String[] tmpStrArr = chroText.substring(0, chroText.indexOf(";")).split(" ");
        peptideEle.setAttribute("peak_sigma", String.valueOf(gModel.getSigma()));
        peptideEle.setAttribute("peak_x", String.valueOf(gModel.getX()));
        peptideEle.setAttribute("peak_y", String.valueOf(gModel.getY()));
        peptideEle.setAttribute("start_scan", String.valueOf((int)Double.parseDouble(tmpStrArr[1])));
        peptideEle.setAttribute("end_scan", String.valueOf((int)Double.parseDouble(tmpStrArr[2])));
        peptideEle.setAttribute("start_rt", String.valueOf(Double.parseDouble(tmpStrArr[3])));
        peptideEle.setAttribute("end_rt", String.valueOf(Double.parseDouble(tmpStrArr[4])));
        peptideEle.setAttribute("peak_area", String.valueOf(gModel.getPeakArea()));

        peptideEle.setAttribute("peak_sigmaIonInjectionCorrection", String.valueOf(gModel.getSigmaIonInjectCorrection()));
        peptideEle.setAttribute("peak_xIonInjectionCorrection",String.valueOf(gModel.getxIonInjectionCorrection()));
        peptideEle.setAttribute("peak_yIonInjectionCorrection", String.valueOf(gModel.getyIonInjectionCorrection()));
        peptideEle.setAttribute("peak_areaIonInjectionCorrection",String.valueOf(gModel.getPeakAreaIonInjectionCorrection()));
        peptideEle.setAttribute("peak_detected",String.valueOf(gModel.isPeakDetected()));



        Element gaussianPeaksEle = new Element("peaks");
        gaussianPeaksEle.setText(gModel.getGaussianPeakString());
        peptideEle.addContent(chro);
        peptideEle.addContent(gaussianPeaksEle);

        return peptideEle;
    }


    public static Element getPeptideDomElement(Peptide peptide, IsotopeReader isoReader,
                                               String eachPath,
                                               Hashtable<String, IndexedFile> origMs1FileHt,
                                               HashMap<String, HashMap<Integer, Integer>> ms2ToMs1Map)
            throws IOException, SQLException, CensusIndexOutOfBoundException, InvalidAAException {
        String pepSequence = peptide.getSequence();
        String fileName=peptide.getFileName();
        fileName = ChroGenerator.cleanFileName(fileName);

        //trim additional characters from peptide sequence at both ends

        char[] ch = pepSequence.substring(2, peptide.getSequence().length()-2).toCharArray();


        ElementComposition element = new ElementComposition(ch, 0, ch.length, isoReader.getIsotope());

        try {
            element.calculate();
        }
        catch (InvalidAAException invE)
        {
            System.out.println("Not Quantifiable peptide : " + pepSequence);
            return null;
        }
        catch(Exception e) {
            System.out.println("Not Quantifiable peptide : " + pepSequence);
            return null;
        }

        if(!element.isQuantifiable())
        {
            System.out.print("\nError : ");
            System.out.println(pepSequence + " is not quantifiable.");
            return null;
        }

        Configuration conf = Configuration.getInstance();

        IsotopeDist sampleDist = null;

        sampleDist = new IsotopeDist(element.getElementSampleArr(), element.getModShift(), true);

        //String origFileKey = eachPath + fileName+".ms1";
        String origFile = eachPath + fileName+".ms1"; //splitSpectraMap.get(eachPath + fileName+".ms1");
        int ms2ScanNumber = peptide.getScanNumber();

        //  splitSpectraMap.get(fileKey)

        HashMap<Integer, Integer> ms2Ms1ScanMap = ms2ToMs1Map.get(origFile);

        int ms1ScanNum = ms2Ms1ScanMap.get(ms2ScanNumber);
        String fileKey = ms1ScanNum + "\t" + origFile;


        IndexedFile origIFile = origMs1FileHt.get(origFile);



        TIntDoubleHashMap retentionTimeMap = origIFile.getRetentionTimeMap();
        TIntDoubleHashMap ionInjectionMap = origIFile.getIonInjectionMap();


        double retTime = retentionTimeMap.get( ms2ScanNumber );
        double ionInjectionTime = ionInjectionMap.get(ms2ScanNumber);

        int tmpScanNumber = ms2ScanNumber;
        int tmpCount=0;
        while(retTime<=0) {
            retTime = retentionTimeMap.get(--tmpScanNumber);

            if(tmpScanNumber<=0) break;

            tmpCount++;
            if(tmpCount>20000) { //check if ms1.index file contains retention time
                System.out.println(">>>" +ms2ScanNumber);
                System.out.println("==retention time is required in ms1.index file.");
                System.exit(0);
            }
        }

        ionInjectionTime = ionInjectionMap.get(tmpScanNumber);


        // System.out.println( retentionTimeMap.get( peptide.getScanNumber()) );
        //  System.out.println( retTime + " == " + tmpScanNumber );

        Element peptideEle=LabelfreeChroUtil.createXmlChroPeptideTitle(peptide); //true is for full scan
        peptideEle.setAttribute("rt", String.valueOf(retTime));
        peptideEle.setAttribute("iit", String.valueOf(ionInjectionTime));
        conf.setCalcSamAvgMass(sampleDist.getAvgMass());
        peptideEle.setAttribute("lightStartMass", String.valueOf(sampleDist.getStartMass()));
        peptideEle.setAttribute("lightAvgMass", String.valueOf(conf.getCalcSamAvgMass()));

        //String fileName = peptide.getFileName().substring(0, peptide.getFileName().indexOf(".")+1) + "ms1";

        //  SpecRange range = null;

        int tmpScanNum = Integer.parseInt(peptide.getScanNum());
        //    range = new SpecRange( tmpScanNum, tmpScanNum );

        int scanNum = Integer.parseInt(peptide.getScanNum());

        switch(conf.getSpectrumFormat())
        {
            case Configuration.MS_FILE_FORMAT:
                fileName += ".ms1";
                break;

            case Configuration.MZXML_FILE_FORMAT:
                fileName += ".mzXML";
                break;

            default:
                break;
        }



        //iFile = ht.get(fileName);
       /* if(null == splitIFile)
        {
            splitIFile = splitMs1FileHt.get(eachPath + fileName.substring(1, fileName.length()));
            if(null == splitIFile) {
                System.out.println("Error : cannot find the file " + eachPath + fileName.substring(1, fileName.length()));
                throw new IOException("Error : cannot find the file " + eachPath + fileName);
            }
            //System.exit(0);
        }*/

        int[] keys = origIFile.getKeys();
        int keyIndex=-1;
        keyIndex = Arrays.binarySearch(keys, scanNum);

        if(keyIndex<0) //Cannot find index
            keyIndex=-(++keyIndex); //Math.abs(++keyIndex);

        if(keyIndex>=keys.length)
            keyIndex--;

        if(keyIndex>0) keyIndex--;  //adjust

        Element chro = new Element("chro");

        int chargeState = Integer.parseInt(peptide.getChargeState());

        String chroText="";


        double[] samIsoArr = sampleDist.getHighMassList();

        for(int i=0;i<samIsoArr.length;i++)
        {
            samIsoArr[i] = (samIsoArr[i]+chargeState*ChroGenerator.PROTON_MASS)/chargeState;
        }

        GaussianPeakModel gModel = calculateLabelfreeFullMS(
                peptide,
                keyIndex,
                origIFile,
                samIsoArr
        );

        chroText = gModel.getChroData();

        //peak area
        double peakAreaTemp = gModel.getPeakArea();
        if(peakAreaTemp<=0) {



            int[] tempKeys = origIFile.getKeys();
          //  TIntLongHashMap indexTemp = origIFile.getMsIndex();
           // RandomAccessFile randomAccessFile = origIFile.getFile();

            if(tempKeys.length<=keyIndex) {
                System.out.println("out of bound exception for split file");
                throw new CensusIndexOutOfBoundException();
            }

            SpectraDB spectraDB = origIFile.getSpectraDB();

            long backGroundNoise = CalcUtilGeneric.getBackGroundNoise(keys[keyIndex], spectraDB);
            gModel.setPeakArea(backGroundNoise);
            gModel.setPeakDetected(false);
          /*  if(pepSequence.equals("R.SEHEVSEIIDGLSEQENLEK.Q"))
            {
                System.out.println("<<>>< "+backGroundNoise);
            }*/


        }

        chro.setText( chroText );

        String[] tmpStrArr = chroText.substring(0, chroText.indexOf(";")).split(" ");
        peptideEle.setAttribute("peak_sigma", String.valueOf(gModel.getSigma()));
        peptideEle.setAttribute("peak_x", String.valueOf(gModel.getX()));
        peptideEle.setAttribute("peak_y", String.valueOf(gModel.getY()));
        peptideEle.setAttribute("start_scan", String.valueOf((int)Double.parseDouble(tmpStrArr[1])));
        peptideEle.setAttribute("end_scan", String.valueOf((int)Double.parseDouble(tmpStrArr[2])));
        peptideEle.setAttribute("start_rt", String.valueOf(Double.parseDouble(tmpStrArr[3])));
        peptideEle.setAttribute("end_rt", String.valueOf(Double.parseDouble(tmpStrArr[4])));
        peptideEle.setAttribute("peak_area", String.valueOf(gModel.getPeakArea()));

        peptideEle.setAttribute("peak_sigmaIonInjectionCorrection", String.valueOf(gModel.getSigmaIonInjectCorrection()));
        peptideEle.setAttribute("peak_xIonInjectionCorrection",String.valueOf(gModel.getxIonInjectionCorrection()));
        peptideEle.setAttribute("peak_yIonInjectionCorrection", String.valueOf(gModel.getyIonInjectionCorrection()));
        peptideEle.setAttribute("peak_areaIonInjectionCorrection",String.valueOf(gModel.getPeakAreaIonInjectionCorrection()));
        peptideEle.setAttribute("peak_detected",String.valueOf(gModel.isPeakDetected()));



        Element gaussianPeaksEle = new Element("peaks");
        gaussianPeaksEle.setText(gModel.getGaussianPeakString());
        peptideEle.addContent(chro);
        peptideEle.addContent(gaussianPeaksEle);

        return peptideEle;
    }

    public static Element getPeptideDomElement(Peptide peptide, IsotopeReader isoReader,
                                               String eachPath, Map<String, TimsTOFXICDB> timsTOFXICDBMap)
            throws IOException, SQLException, CensusIndexOutOfBoundException, InvalidAAException {
        String pepSequence = peptide.getSequence();
        String fileName=peptide.getFileName();
        fileName = ChroGenerator.cleanFileName(fileName);

        //trim additional characters from peptide sequence at both ends

        char[] ch = pepSequence.substring(2, peptide.getSequence().length()-2).toCharArray();


        ElementComposition element = new ElementComposition(ch, 0, ch.length, isoReader.getIsotope());

        try {
            element.calculate();
        }
        catch (InvalidAAException invE)
        {
            System.out.println("Not Quantifiable peptide : " + pepSequence);
            return null;
        }
        catch(Exception e) {
            System.out.println("Not Quantifiable peptide : " + pepSequence);
            return null;
        }

        if(!element.isQuantifiable())
        {
            System.out.print("\nError : ");
            System.out.println(pepSequence + " is not quantifiable.");
            return null;
        }

        Configuration conf = Configuration.getInstance();

        IsotopeDist sampleDist = null;

        sampleDist = new IsotopeDist(element.getElementSampleArr(), element.getModShift(), true);

        //String origFileKey = eachPath + fileName+".ms1";
        String origFile =  fileName+".ms2"; //splitSpectraMap.get(eachPath + fileName+".ms1");
        int ms2ScanNumber = peptide.getScanNumber();
        TimsTOFXICDB timsTOFXICDB = timsTOFXICDBMap.get(origFile);

        //  splitSpectraMap.get(fileKey)
        TimsTOFXICDB.TimstofQueryResult queryResult = timsTOFXICDB.queryAndSumMS2(ms2ScanNumber);
        List<Pair<Double,Double>> peakList = queryResult.getSummedList();
        TDoubleArrayList xarrayList = new TDoubleArrayList();
        TDoubleArrayList yarrayList = new TDoubleArrayList();
        double max = Double.MIN_VALUE;
        double sum = 0;
        for(Pair< Double,Double> r: peakList)
        {
            xarrayList.add(r.getLeft());
            yarrayList.add(r.getRight());
            if(r.getRight()> max)
                max = r.getRight();
            sum+=r.getRight();
            //System.out.println(r.getLeft()+"\t"+r.getRight());
        }
        double [] xarr = xarrayList.toNativeArray();
        double [] yarr = yarrayList.toNativeArray();
        GaussianPeakModel gModel =  getGaussianPeakRangeIndex(xarr, yarr, -1, peakList.size()-1);
        gModel.setMaxIntensity(max);
        double peakHeight = gModel.getY();
        double sigma = gModel.getSigma();
        double area = getGaussianPeakArea(peakHeight, sigma);


        gModel.setPeakArea(area);

        // System.out.println( retentionTimeMap.get( peptide.getScanNumber()) );
        //  System.out.println( retTime + " == " + tmpScanNumber );

        Element peptideEle=LabelfreeChroUtil.createXmlChroPeptideTitle(peptide); //true is for full scan
       peptideEle.setAttribute("rt", String.valueOf(queryResult.retTime));
        conf.setCalcSamAvgMass(sampleDist.getAvgMass());
        peptideEle.setAttribute("lightStartMass", String.valueOf(sampleDist.getStartMass()));
        peptideEle.setAttribute("lightAvgMass", String.valueOf(conf.getCalcSamAvgMass()));

        //String fileName = peptide.getFileName().substring(0, peptide.getFileName().indexOf(".")+1) + "ms1";

        //  SpecRange range = null;

;

        peptideEle.setAttribute("peak_sigma", String.valueOf(gModel.getSigma()));
        peptideEle.setAttribute("peak_x", String.valueOf(gModel.getX()));
        peptideEle.setAttribute("peak_y", String.valueOf(gModel.getY()));
    /*    peptideEle.setAttribute("start_scan", String.valueOf((int)Double.parseDouble(tmpStrArr[1])));
        peptideEle.setAttribute("end_scan", String.valueOf((int)Double.parseDouble(tmpStrArr[2])));
        peptideEle.setAttribute("start_rt", String.valueOf(Double.parseDouble(tmpStrArr[3])));
        peptideEle.setAttribute("end_rt", String.valueOf(Double.parseDouble(tmpStrArr[4])));*/
        peptideEle.setAttribute("peak_area", String.valueOf(gModel.getPeakArea()));
/*
        peptideEle.setAttribute("peak_sigmaIonInjectionCorrection", String.valueOf(gModel.getSigmaIonInjectCorrection()));
        peptideEle.setAttribute("peak_xIonInjectionCorrection",String.valueOf(gModel.getxIonInjectionCorrection()));
        peptideEle.setAttribute("peak_yIonInjectionCorrection", String.valueOf(gModel.getyIonInjectionCorrection()));
        peptideEle.setAttribute("peak_areaIonInjectionCorrection",String.valueOf(gModel.getPeakAreaIonInjectionCorrection()));
        peptideEle.setAttribute("peak_detected",String.valueOf(gModel.isPeakDetected()));

*/

        Element gaussianPeaksEle = new Element("peaks");
        gaussianPeaksEle.setText(gModel.getGaussianPeakString());
        //peptideEle.addContent(chro);
        peptideEle.addContent(gaussianPeaksEle);

        return peptideEle;
    }



    public static double[] readSplitFastFullSpectrum(
        //    IndexedFile iFile,
            int curIndex,
            double[] isoArr,
            int[] origKeys,
            String origFile,
            Hashtable<String, IndexedFile> splitMs1FileHt,
            HashMap<String, String> splitSpectraMap,
            double pepMass
    ) throws IOException, CensusIndexOutOfBoundException, Exception
    {

        //temp......................................
        int tmpMs1Scan = origKeys[curIndex];

        String fileKey = tmpMs1Scan + "\t" + origFile;
        String spltiMs1File = splitSpectraMap.get(fileKey);
        IndexedFile splitIFile = splitMs1FileHt.get(spltiMs1File);
        int splitCurIndex = splitIFile.getIndexByScan(tmpMs1Scan);

        int[] splitKeys = splitIFile.getKeys();
        TIntLongHashMap index = splitIFile.getMsIndex();
        RandomAccessFile randomAccessFile = splitIFile.getFile();

        if(splitKeys.length<=splitCurIndex) {
            System.out.println("out of bound exception for split file");
            throw new CensusIndexOutOfBoundException();
        }

        long startPos = index.get(splitKeys[splitCurIndex]);
        long endPos;

        randomAccessFile.seek(startPos);

        if( (splitCurIndex+1)>=splitKeys.length )
            endPos = randomAccessFile.length();
        else
            endPos = index.get(splitKeys[splitCurIndex+1]);

        int byteSize = (int)(endPos-startPos);

        byte[] bytes = new byte[byteSize];
        randomAccessFile.readFully(bytes);

        char ch;
        int pos=0;

        ch = (char)bytes[pos];

        double[] massArr = null;
        double[] intArr = null;
        try {
            //Remove Z, S, I, D lines
            while ((ch = (char) bytes[pos]) == 'Z' || ch == 'S' || ch == 'I' || ch == 'D') {
                while (ch != CalcUtil.CARRIAGE_RETURN) {
                    pos++;
                    ch = (char) bytes[pos];
                }

                pos++;
            }


            double[][] specResult = CalcUtil.parseSpectra(bytes);
            massArr = specResult[0];
            intArr = specResult[1];
        } catch (Exception e) {
            System.out.println("no peaks");
        }



        double result[] = new double[9];

        result[0] = splitKeys[splitCurIndex];

        Configuration conf = Configuration.getInstance();


        double[] tempArr = CalcUtilGeneric.intensitySumWithIsotopeModeling(massArr, intArr, isoArr, conf.getMassTolerance(), pepMass);
       // double[] tempArr1 = intensitySum(massArr, intArr, isoArr);

        if(null != tempArr) {
            result[1] = tempArr[0];     //intensity
            result[3] = isoArr.length;
            result[5] = tempArr[1];      //Found iso num
            result[7] = tempArr[2];         //tolerance
        }

        if(null != isoArr) {
            tempArr = intensitySum(massArr, intArr, isoArr);
            if(null != tempArr) {
                result[2] = tempArr[0];
                result[4] = isoArr.length;
                result[6] = tempArr[1];
                result[8] = tempArr[2];
            }
        }

        randomAccessFile.close();

        return result;
    }


    public static double[] readSplitFastFullSpectrum(
            int curIndex,
            double[] isoArr,
            int[] origKeys,
            SpectraDB spectraDB,
            double pepMass
    ) throws SQLException {

        //temp......................................
        int tmpMs1Scan = origKeys[curIndex];
        SpectraDB.Spectrum spectrum = spectraDB.getSpectrumFromDB(tmpMs1Scan);
        double [] massArr = spectrum.getMzList().toNativeArray();
        double [] intArr = spectrum.getIntensityList().toNativeArray();


        double result[] = new double[9];

        result[0] = tmpMs1Scan;

        Configuration conf = Configuration.getInstance();


        double[] tempArr = CalcUtilGeneric.intensitySumWithIsotopeModeling(massArr, intArr, isoArr, conf.getMassTolerance(), pepMass);
        // double[] tempArr1 = intensitySum(massArr, intArr, isoArr);

        if(null != tempArr) {
            result[1] = tempArr[0];     //intensity
            result[3] = isoArr.length;
            result[5] = tempArr[1];      //Found iso num
            result[7] = tempArr[2];         //tolerance
        }

        if(null != isoArr) {
            tempArr = intensitySum(massArr, intArr, isoArr);
            if(null != tempArr) {
                result[2] = tempArr[0];
                result[4] = isoArr.length;
                result[6] = tempArr[1];
                result[8] = tempArr[2];
            }
        }


        return result;
    }



    public static GaussianPeakModel calculateLabelfreeFullMS(Peptide peptide,
                                                         int keyIndex,
                                                         IndexedFile origIFile,
                                                         double[] samIsoArr,
                                                         //SpecRange range,
                                                         String origFile,
                                                         Hashtable<String, IndexedFile> splitMs1FileHt,
                                                         HashMap<String, String> splitSpectraMap

    )
            throws IOException, CensusIndexOutOfBoundException, Exception
    {

        //peptide.getChargeState();
        int[] origIFileKeys = origIFile.getKeys();
        Configuration conf = Configuration.getInstance();
        double pepMass = Double.parseDouble(peptide.getCalcMHplus());
      //  CalcUtil.massTolerance = conf.getMassTolerance();

        /*CalcUtil.SpectrumModel sModel = new CalcUtil.SpectrumModel();
        sModel.setHighRes(true);
        sModel.setKeys(keys);
        sModel.setIndex(index);
        sModel.setIFile(iFile);
        sModel.setFile(file);
        sModel.setMzxmlreader(mzReader);
        sModel.setSamIsoArr(samIsoArr);*/

        double steepRatioThreshold = conf.getSteepRatioThreshold();

        int numIsoWindow = conf.getNumOfIsolationWindow();
        int maxWindow=conf.getMaxWindow();
        int margin = conf.getMargin();

        double[][] result = null;
        //conf.getQuantLevel()

        maxWindow = maxWindow * 3;
        result = new double[maxWindow*2+1+margin*2][3];

        int leftIndex=maxWindow+margin;
        int rightIndex=maxWindow+margin+1;


        //double leftTotalIntensity=0;
        //double rightTotalIntensity=0;
        double totalIntensity=0;

        int steepArea = conf.getSteepArea();

        int moveLeftKeyIndex = keyIndex;

        int moveRightKeyIndex = keyIndex+1*numIsoWindow;
        int identifiedScan = origIFileKeys[keyIndex];
/*
System.out.println("==============" + keyIndex + " " + numIsoWindow);
System.out.println("==============" + keyIndex + " " + numIsoWindow);
System.out.println("==============" + keyIndex + " " + numIsoWindow);
System.out.println("==============" + origIFileKeys[keyIndex]);
System.out.println("==============" + numIsoWindow);
*/
        int initWin=2;

        for(int i=0;i<initWin;i++) ///Index=maxWindow;keys[leftPeakIndex]>=leftMinValue; )
        {
            if(moveLeftKeyIndex<=0 || leftIndex<=0)
            {
                moveLeftKeyIndex += 1*numIsoWindow;
                leftIndex++;

                break;
            }

            result[leftIndex] = readSplitFastFullSpectrum(moveLeftKeyIndex,samIsoArr, origIFileKeys, origFile, splitMs1FileHt, splitSpectraMap, pepMass);

//            result[leftIndex] = readFullSpectrum(keys, moveLeftKeyIndex, index, file, samIsoArr, null);
            totalIntensity += result[leftIndex][1];

            // if(moveLeftKeyIndex<=0 || leftIndex<=0)
            //	break;

            moveLeftKeyIndex -= numIsoWindow;
            leftIndex--;
        }

        for(int i=0;i<initWin;i++) //leftIndex=maxWindow;keys[leftPeakIndex]>=leftMinValue; )
        {
            if(moveRightKeyIndex>=origIFileKeys.length)
            {
                moveRightKeyIndex -= 1*numIsoWindow;
                rightIndex--;

                break;
            }

            result[rightIndex] = readSplitFastFullSpectrum(moveRightKeyIndex, samIsoArr, origIFileKeys, origFile, splitMs1FileHt, splitSpectraMap, pepMass);

            totalIntensity += result[rightIndex][1];

            moveRightKeyIndex += numIsoWindow;
            rightIndex++;

		if(totalIntensity<=0) break;

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
                        arr[steepCount] = result[rightIndex] = readSplitFastFullSpectrum(tempKeyIndex, samIsoArr, origIFileKeys, origFile, splitMs1FileHt, splitSpectraMap, pepMass);
                    }

                    steepCount++;
if(arr[steepCount-1][1]<=0) break;
                    //tempKeyIndex--;
                    tempKeyIndex = tempKeyIndex - 1*numIsoWindow;
                }

                if(tempKeyIndex<0 || arr[steepCount][1]<=0  ) {
                    break;

}

                arr[steepCount] = readSplitFastFullSpectrum(tempKeyIndex, samIsoArr, origIFileKeys, origFile, splitMs1FileHt, splitSpectraMap, pepMass);

if(arr[steepCount][1]<=0) break;

                for(int k=0;k<3;k++)
                {
                    for(int l=0;l<3;l++)
                        prevArr[k][l] = arr[k][l];
                }


                for(int j=leftIndex+1;j<=leftIndex+steepCount+1;j++)
                {
                    area1 += result[j][1];
                    area1 += result[j][2];
                }

                for(int j=0;j<=steepCount;j++)
                {
                    area2 += arr[j][1];
                    area2 += arr[j][2];
                }

                if(area2==0 && area1==0)
                    break;

                double areaRatio = (double)area2/area1;

                /*****************************************
                 * AREA1 is right side area
                 * AREA2 is left side area
                 *****************************************/
                if(isHighIntensity && ((double)totalIntensity/(rightIndex-leftIndex-1)/1.5>area2/3))
                    isHighIntensity = false;

                //if( (isGoingUp && (float)leftTotalIntensity/(rightIndex-leftIndex-1)/3<area2/3) || areaRatio<steepRatioThreshold )
                //if( (isGoingUp && (float)totalIntensity/(rightIndex-leftIndex-1)/3<area2/3) || areaRatio<steepRatioThreshold )
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

                if(moveRightKeyIndex+steepArea>=origIFileKeys.length)
                    break;

                int steepCount=0;

                int tempKeyIndex = moveRightKeyIndex;
                while(true)
                {
                    if(tempKeyIndex>=origIFileKeys.length || steepCount>= steepArea-1 || steepCount>=origIFileKeys.length)
                        break;

                    if(prevArr[0][0]!=0)
                    {
                        for(int l=0;l<2;l++)
                            for(int m=0;m<3;m++)
                                arr[l][m] = prevArr[l+1][m];
                    }
                    else
                    {
                        arr[steepCount] = readSplitFastFullSpectrum(tempKeyIndex, samIsoArr, origIFileKeys, origFile, splitMs1FileHt, splitSpectraMap, pepMass);

                    }

                    steepCount++;
                    tempKeyIndex = tempKeyIndex + 1*numIsoWindow;
                    //tempKeyIndex++;
                }

		if(arr[steepCount-1][1]<=0) break;
                if(tempKeyIndex>=origIFileKeys.length)
                    break;

                arr[steepCount] = readSplitFastFullSpectrum(tempKeyIndex, samIsoArr, origIFileKeys, origFile, splitMs1FileHt, splitSpectraMap, pepMass);
		if(arr[steepCount][1]<=0) break;

                for(int k=0;k<3;k++)
                {
                    for(int l=0;l<3;l++)
                        prevArr[k][l] = arr[k][l];
                }

                for(int j=rightIndex-1;j>rightIndex-steepCount-2;j--)
                {
                    area1 += result[j][1];
                    area1 += result[j][2];
                }

                for(int j=0;j<=steepCount;j++)
                {
                    area2 += arr[j][1];
                    area2 += arr[j][2];
                }

                if(area2==0 && area1==0)
                    break;

                double areaRatio = (double)area2/area1;

                /*****************************************
                 * AREA1 is left side area
                 * AREA2 is right side area
                 *****************************************/

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

        if(leftIndex<0) leftIndex=0;
        if(rightIndex>=result.length) rightIndex=result.length-1;


//        int peakStart = leftIndex+1;
//        int peakEnd = rightIndex-1;
      int tmpLeftIndex = leftIndex;
      int tmpRightIndex = rightIndex;
      int tmpMoveLeftKeyIndex = moveLeftKeyIndex;
      int tmpMoveRightKeyIndex = moveRightKeyIndex;

      for(int i=0;i<100;i++) ///Index=maxWindow;keys[leftPeakIndex]>=leftMinValue; )
      {
        if(tmpMoveLeftKeyIndex<0 || tmpLeftIndex<0)
        {
          tmpMoveLeftKeyIndex += 1*numIsoWindow;
          tmpLeftIndex++;

          break;
        }

        result[tmpLeftIndex] = readSplitFastFullSpectrum(tmpMoveLeftKeyIndex,samIsoArr, origIFileKeys, origFile, splitMs1FileHt, splitSpectraMap, pepMass);

        tmpMoveLeftKeyIndex -= numIsoWindow;
        tmpLeftIndex--;

if(result[tmpLeftIndex+1][1]<=0) {

break;
}
      }

      for(int i=0;i<100;i++) //leftIndex=maxWindow;keys[leftPeakIndex]>=leftMinValue; )
      {
        if(tmpMoveRightKeyIndex>=origIFileKeys.length || tmpRightIndex>=result.length)
        {
          tmpMoveRightKeyIndex -= 1*numIsoWindow;
          tmpRightIndex--;

          break;
        }

        result[tmpRightIndex] = readSplitFastFullSpectrum(tmpMoveRightKeyIndex, samIsoArr, origIFileKeys, origFile, splitMs1FileHt, splitSpectraMap, pepMass);

        tmpMoveRightKeyIndex += numIsoWindow;
        tmpRightIndex++;



if(result[tmpRightIndex-1][1]<=0) break;

      }




        if(moveRightKeyIndex>=origIFileKeys.length)
        {
            moveRightKeyIndex -= 1*numIsoWindow;
            rightIndex--;
        }

        result[rightIndex] = readSplitFastFullSpectrum(moveRightKeyIndex, samIsoArr, origIFileKeys, origFile, splitMs1FileHt, splitSpectraMap, pepMass);


        result[leftIndex] = readSplitFastFullSpectrum(moveLeftKeyIndex, samIsoArr, origIFileKeys, origFile, splitMs1FileHt, splitSpectraMap, pepMass);


        //readFullSpectrum(keys, moveLeftKeyIndex, index, file, samIsoArr, null);
/*
      if(moveRightKeyIndex<5) {

        result[--leftIndex] = readSplitFastFullSpectrum(--moveLeftKeyIndex, samIsoArr, origIFileKeys, origFile, splitMs1FileHt, splitSpectraMap, pepMass);
        result[++rightIndex] = readSplitFastFullSpectrum(++moveRightKeyIndex, samIsoArr, origIFileKeys, origFile, splitMs1FileHt, splitSpectraMap, pepMass);
      } */

        //long[] chromPeakArr = new long[rightIndex-leftIndex+1];
      long[] chromPeakArr = new long[tmpRightIndex-tmpLeftIndex+1];


        int count=0;

        TIntDoubleHashMap scanRtMap = origIFile.getScanRtMap();
        TIntDoubleHashMap scanIonMap = origIFile.getIonInjectionMap();

        double[] retArr = new double[chromPeakArr.length];
	double identifiedRt = scanRtMap.get(identifiedScan);
//System.out.println("ret--------------" + identifiedScan);
//System.out.println("ret--------------" + scanRtMap.get(identifiedScan));
        double[] ionArr = new double[chromPeakArr.length];


        //for(int i=leftIndex;i<=rightIndex;i++) {
      for(int i=tmpLeftIndex;i<=tmpRightIndex;i++) {
            try {
                chromPeakArr[count] = (long)result[i][1];
                retArr[count] = scanRtMap.get((int)result[i][0]);
                count++;
            } catch(Exception e) {

                e.printStackTrace();
            }
        }
        count = 0;
        //for(int i=leftIndex;i<=rightIndex;i++) {
          for(int i=tmpLeftIndex;i<=tmpRightIndex;i++) {
            try {
                ionArr[count++] = scanIonMap.get((int)result[i][0]);
                //count++;
            } catch(Exception e) {
              System.out.println(ionArr.length + "\t" + count + "\t" +result.length + "\t" + i);

                e.printStackTrace();
            }
        }

        long[] chroIonArr = new long[chromPeakArr.length];
        for(int i =0;i<chromPeakArr.length;i++){
            if(i==0){
                chroIonArr[i] =  (long) ((chromPeakArr[i]/ionArr[i])*(retArr[i+1]-retArr[i]));
            }else if(i==chromPeakArr.length-1){
                chroIonArr[i] =  (long) ((chromPeakArr[i]/ionArr[i])*(retArr[i]-retArr[i-1]));
            }
            else{
                chroIonArr[i] =  (long) ((chromPeakArr[i]/ionArr[i])*(((retArr[i]-retArr[i-1])/2)+((retArr[i+1]-retArr[i])/2)));
            }
        }

        double[] smoothchroIonArr = Smooth.smoothAsDouble(chroIonArr, LabelfreeMissingPeptideBuilderSplit.SMOOTH_WINDOW_SIZE);
        double[] smoothChromArr = Smooth.smoothAsDouble(chromPeakArr, LabelfreeMissingPeptideBuilderSplit.SMOOTH_WINDOW_SIZE);

      //for(double d:chromPeakArr)
      //  System.out.println("===\t" + d);

        double ionbasePeak=0;
        double basePeak=0;
        int basePeakIndex=0;
        /*
        for(int i=0;i<smoothChromArr.length;i++) {

       //     System.out.println("==\t" + retArr[i] + "\t" + smoothChromArr[i] + "\t" + chromPeakArr[i]);
            if(basePeak<smoothChromArr[i]) {
                basePeak = smoothChromArr[i];
                basePeakIndex=i;
            }

           // System.out.println(smoothChromArr[i]);
        }

        int ionbasePeakIndex=0;

        for(int i=0;i<smoothchroIonArr.length;i++) {

       //     System.out.println("==\t" + retArr[i] + "\t" + smoothChromArr[i] + "\t" + chromPeakArr[i]);
            if(ionbasePeak<smoothchroIonArr[i]) {
                ionbasePeak = smoothchroIonArr[i];
                ionbasePeakIndex=i;
            }

           // System.out.println(smoothChromArr[i]);
        }
        */
        int identifiedIndex = maxWindow+margin-leftIndex;
      if(identifiedIndex<0) identifiedIndex=0;
        basePeak = smoothChromArr[identifiedIndex];
        ionbasePeak = smoothchroIonArr[identifiedIndex];
        /***********************************
         * Find simple/rough peak range (1/3 of base peak) for Gaussian input
         ***********************************/
        int[] indexResult = LabelfreeChroUtil.getPeakRange(identifiedIndex, basePeak, smoothChromArr);
        int peakStartIndex = indexResult[0];
        int peakEndIndex = indexResult[1];

        int[] indexIonResult = LabelfreeChroUtil.getPeakRange(identifiedIndex, ionbasePeak, smoothchroIonArr);
        int peakStartIonIndex = indexIonResult[0];
        int peakEndIonIndex = indexIonResult[1];


        /***********************************
         * Find peak range based on Gaussian fitting
         ***********************************/
//System.out.println("==================");
        GaussianPeakModel gModel = GaussianFitting.getGaussianPeakRangeIndex2(
          peptide, retArr, smoothChromArr, peakStartIndex, peakEndIndex, ionArr,smoothchroIonArr,peakStartIonIndex,peakEndIonIndex);

      double maxPeakIntensity = 0;
      for(double d:smoothChromArr) {
        if(d>maxPeakIntensity)
          maxPeakIntensity = d;
      }

      gModel.setMaxIntensity(maxPeakIntensity);

        //backGroundNoise = CalcUtilGeneric.getBackGroundNoise(isoArr, currentPos, diff, splitIFile);


        int finalPeakStartIndex=gModel.getPeakStartIndex()+leftIndex;
        //int finalPeakStartIndex=gModel.getPeakStartIndex();
        if(finalPeakStartIndex<0) finalPeakStartIndex=0;
        //int finalPeakEndIndex=gModel.getPeakEndIndex();
        int finalPeakEndIndex=gModel.getPeakEndIndex()+leftIndex;
        if(finalPeakEndIndex>=result.length)
            finalPeakEndIndex=result.length-1;

        double peakArea = gModel.getPeakArea();

       // System.out.println("peak aaaaaaaaaaaaa" + peakArea);
        double peakAreaIon = gModel.getPeakAreaIonInjectionCorrection();

        /***********************************
         * build chro data
         ***********************************/
        StringBuffer sb = new StringBuffer();

        //  try {
        //  sb.append("P ").append((long)result[peakStart][0]).append(" ").append((long)result[peakEnd][0]).append(";");
        sb.append("P ").append(result[finalPeakStartIndex][0]).append(" ").append(result[finalPeakEndIndex][0]).
                append(" ").append(retArr[gModel.getPeakStartIndex()]).append(" ").append(retArr[gModel.getPeakEndIndex()]).append(";");
        //.append( ).append(";");
        /*    } catch(Exception e) {

            }*/


        Gaussian g = null;

        if(gModel.getX()>0)
            g = new Gaussian(gModel.getY(), gModel.getX(), gModel.getSigma());

        double gstart = -4 * gModel.getSigma() + gModel.getX();
        double gend = 4 * gModel.getSigma() + gModel.getX();

        List<Double> gxList = new ArrayList<>();
        List<Double> gyList = new ArrayList<>();

      count=1; //it starts with 1. why??


 //     for(int i=leftIndex+1;i<rightIndex-1;i++)
      for(int i=0;i<smoothChromArr.length;i++)
        {
          if(retArr[i]<=0) continue;
            sb.append((int)result[i+tmpLeftIndex][0]).append(" "); //scan
            sb.append(retArr[i]).append(" "); //ret time
            sb.append((int)smoothChromArr[i]).append(" "); //intensity

            count++;
            sb.append(result[i][result[i].length-2]).append(" ");
            sb.append(result[i][result[i].length-1]).append(";");

        }




        // if(null != g && retArr[count]>=gstart && retArr[count]<=gend) {
        if(null != g) { // && retArr[count]>=gstart && retArr[count]<=gend) {


          ///////////


          double[] allRtArr = origIFile.getRtArr();
          int rtStartkeyIndex = Arrays.binarySearch(allRtArr, gstart);
          int rtEndkeyIndex = Arrays.binarySearch(allRtArr, gend);

          if (rtStartkeyIndex < 0) //Cannot find index
          {
            rtStartkeyIndex = -(++rtStartkeyIndex); //Math.abs(++keyIndex);
          }
          if (rtStartkeyIndex >= allRtArr.length) {
            rtStartkeyIndex--;
          }

          if (rtEndkeyIndex < 0) //Cannot find index
          {
            rtEndkeyIndex = -(++rtEndkeyIndex); //Math.abs(++keyIndex);
          }
          if (rtEndkeyIndex >= allRtArr.length) {
            rtEndkeyIndex--;
          }

          for(int i=rtStartkeyIndex;i<=rtEndkeyIndex;i++) {
            gxList.add(allRtArr[i]);
            gyList.add(g.value(allRtArr[i]));
          }

          /*

          for(double d:retArr) {


                gxList.add(d);
                gyList.add(g.value(d));
            }
            */
        }



        if(null != g) {
            double[] gxArr = new double[gxList.size()];
            double[] gyArr = new double[gxList.size()];
            for(int j=0;j<gxList.size();j++) {
             // try {
                gxArr[j] = gxList.get(j).doubleValue();
                gyArr[j] = gyList.get(j).doubleValue();
            //  } catch(Exception e) {
             //   e.printStackTrace();
             // }

            }

            gModel.setGaussianXArr(gxArr);
            gModel.setGaussianYArr(gyArr);
        }

        gModel.setChroData(sb.toString());

//        System.out.println("peak area =====" + gModel.getPeakArea());
        return gModel;

    }

    public static GaussianPeakModel calculateLabelfreeFullMSOrig(Peptide peptide,
            int keyIndex,
            IndexedFile iFile,
            double[] samIsoArr,
            SpecRange range)

        throws IOException, CensusIndexOutOfBoundException, Exception
    {

        peptide.getChargeState();
        TIntLongHashMap index = iFile.getMsIndex();
        RandomAccessFile file = iFile.getFile();
        MzxmlSpectrumReader mzReader = iFile.getMzreader();
        int[] keys = iFile.getKeys();
        Configuration conf = Configuration.getInstance();
        CalcUtil.massTolerance = conf.getMassTolerance();

        /*CalcUtil.SpectrumModel sModel = new CalcUtil.SpectrumModel();
        sModel.setHighRes(true);
        sModel.setKeys(keys);
        sModel.setIndex(index);
        sModel.setIFile(iFile);
        sModel.setFile(file);
        sModel.setMzxmlreader(mzReader);
        sModel.setSamIsoArr(samIsoArr);*/

        double steepRatioThreshold = conf.getSteepRatioThreshold();

	int numIsoWindow = conf.getNumOfIsolationWindow();
        int maxWindow=conf.getMaxWindow();
        int margin = conf.getMargin();

	double[][] result = null;
        //conf.getQuantLevel()

	result = new double[maxWindow*2+1+margin*2][3];

        int leftIndex=maxWindow+margin;
        int rightIndex=maxWindow+margin+1;


        //double leftTotalIntensity=0;
        //double rightTotalIntensity=0;
        double totalIntensity=0;

        int steepArea = conf.getSteepArea();

        int moveLeftKeyIndex = keyIndex;

        int moveRightKeyIndex = keyIndex+1*numIsoWindow;

        int initWin=2;

        for(int i=0;i<initWin;i++) ///Index=maxWindow;keys[leftPeakIndex]>=leftMinValue; )
        {
	    if(moveLeftKeyIndex<=0 || leftIndex<=0)
	    {
		moveLeftKeyIndex += 1*numIsoWindow;
		leftIndex++;

		break;
	    }

	    result[leftIndex] = readFullSpectrum(keys, moveLeftKeyIndex, index, file, samIsoArr, null);
            totalIntensity += result[leftIndex][1];


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


		    result[rightIndex] = readFullSpectrum(keys, moveRightKeyIndex, index, file, samIsoArr, null);

            totalIntensity += result[rightIndex][1];

            moveRightKeyIndex += numIsoWindow;
            rightIndex++;
        }

        boolean isGoingUp=true;
        boolean isHighIntensity=true;  //if the intensity is lower than one third of average intensity, this becomes false
        double[][] arr = new double[steepArea][3];
        double[][] prevArr = new double[steepArea][3];

        if(leftIndex+steepArea>=0)
        {

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
			    arr[steepCount] = readFullSpectrum(keys, tempKeyIndex, index, file, samIsoArr, null);
		    }

                    steepCount++;
                    //tempKeyIndex--;
		    tempKeyIndex = tempKeyIndex - 1*numIsoWindow;
                }

		if(tempKeyIndex<0)
		    break;

		arr[steepCount] = readFullSpectrum(keys, tempKeyIndex, index, file, samIsoArr, null);

		for(int k=0;k<3;k++)
		{
		    for(int l=0;l<3;l++)
			prevArr[k][l] = arr[k][l];
		}


                for(int j=leftIndex+1;j<=leftIndex+steepCount+1;j++)
                {
                    area1 += result[j][1];
                    area1 += result[j][2];
                }

                for(int j=0;j<=steepCount;j++)
                {
                    area2 += arr[j][1];
                    area2 += arr[j][2];
                }

                if(area2==0 && area1==0)
                    break;

                double areaRatio = (double)area2/area1;

                /*****************************************
                 * AREA1 is right side area
                 * AREA2 is left side area
                 *****************************************/
		if(isHighIntensity && ((double)totalIntensity/(rightIndex-leftIndex-1)/1.5>area2/3))
		    isHighIntensity = false;

                //if( (isGoingUp && (float)leftTotalIntensity/(rightIndex-leftIndex-1)/3<area2/3) || areaRatio<steepRatioThreshold )
                //if( (isGoingUp && (float)totalIntensity/(rightIndex-leftIndex-1)/3<area2/3) || areaRatio<steepRatioThreshold )
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
			    arr[steepCount] = readFullSpectrum(keys, tempKeyIndex, index, file, samIsoArr, null);
		    }

                    steepCount++;
		    tempKeyIndex = tempKeyIndex + 1*numIsoWindow;
                    //tempKeyIndex++;
                }

		if(tempKeyIndex>=keys.length)
		    break;

		arr[steepCount] = readFullSpectrum(keys, tempKeyIndex, index, file, samIsoArr, null);

		for(int k=0;k<3;k++)
		{
		    for(int l=0;l<3;l++)
			prevArr[k][l] = arr[k][l];
		}

                for(int j=rightIndex-1;j>rightIndex-steepCount-2;j--)
                {
                    area1 += result[j][1];
                    area1 += result[j][2];
                }

                for(int j=0;j<=steepCount;j++)
                {
                    area2 += arr[j][1];
                    area2 += arr[j][2];
                }

                if(area2==0 && area1==0)
                    break;

                double areaRatio = (double)area2/area1;

                /*****************************************
                 * AREA1 is left side area
                 * AREA2 is right side area
                 *****************************************/

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

        if(leftIndex<0) leftIndex=0;
        if(rightIndex>=result.length) rightIndex=result.length-1;


	int peakStart = leftIndex+1;
	int peakEnd = rightIndex-1;



	/*
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

	    result[leftIndex] = readFullSpectrum(keys, moveLeftKeyIndex, index, file, samIsoArr, null);

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

	    result[rightIndex] = readFullSpectrum(keys, moveRightKeyIndex, index, file, samIsoArr, null);

	    if(0 == result[rightIndex][0])
	    {
		moveRightKeyIndex += numIsoWindow;
		continue;
	    }

            rightIndex++;
	    moveRightKeyIndex += numIsoWindow;
        }
*/
       //
       if(moveRightKeyIndex>=keys.length)
            {
		moveRightKeyIndex -= 1*numIsoWindow;
                rightIndex--;
            }

        result[rightIndex] = readFullSpectrum(keys, moveRightKeyIndex, index, file, samIsoArr, null);


        //moveLeftKeyIndex=leftIndex;
       // moveLeftKeyIndex -= numIsoWindow;
        result[leftIndex] = readFullSpectrum(keys, moveLeftKeyIndex, index, file, samIsoArr, null);


        long[] chromPeakArr = new long[rightIndex-leftIndex+1];
        int count=0;

        TIntDoubleHashMap scanRtMap = iFile.getScanRtMap();
        TIntDoubleHashMap scanIonMap = iFile.getIonInjectionMap();

        double[] retArr = new double[chromPeakArr.length];
        double[] ionArr = new double[chromPeakArr.length];

        for(int i=leftIndex;i<=rightIndex;i++) {
            try {
            chromPeakArr[count] = (long)result[i][1];
            retArr[count] = scanRtMap.get((int)result[i][0]);
            count++;
            } catch(Exception e) {

                e.printStackTrace();
            }
        }
        count = 0;
        for(int i=leftIndex;i<=rightIndex;i++) {
            try {
            ionArr[count] = scanIonMap.get((int)result[i][0]);
            count++;
            } catch(Exception e) {

                e.printStackTrace();
            }
        }

        long[] chroIonArr = new long[chromPeakArr.length];
        for(int i =0;i<chromPeakArr.length;i++){
            if(i==0){
                chroIonArr[i] =  (long) ((chromPeakArr[i]/ionArr[i])*(retArr[i+1]-retArr[i]));
            }else if(i==chromPeakArr.length-1){
                chroIonArr[i] =  (long) ((chromPeakArr[i]/ionArr[i])*(retArr[i]-retArr[i-1]));
            }
            else{
              chroIonArr[i] =  (long) ((chromPeakArr[i]/ionArr[i])*(((retArr[i]-retArr[i-1])/2)+((retArr[i+1]-retArr[i])/2)));
            }
        }




        double[] smoothchroIonArr = Smooth.smoothAsDouble(chroIonArr, LabelfreeMissingPeptideBuilderSplit.SMOOTH_WINDOW_SIZE);
        double[] smoothChromArr = Smooth.smoothAsDouble(chromPeakArr, LabelfreeMissingPeptideBuilderSplit.SMOOTH_WINDOW_SIZE);

        double ionbasePeak=0;
        double basePeak=0;
        int basePeakIndex=0;
        /*
        for(int i=0;i<smoothChromArr.length;i++) {

       //     System.out.println("==\t" + retArr[i] + "\t" + smoothChromArr[i] + "\t" + chromPeakArr[i]);
            if(basePeak<smoothChromArr[i]) {
                basePeak = smoothChromArr[i];
                basePeakIndex=i;
            }

           // System.out.println(smoothChromArr[i]);
        }

        int ionbasePeakIndex=0;

        for(int i=0;i<smoothchroIonArr.length;i++) {

       //     System.out.println("==\t" + retArr[i] + "\t" + smoothChromArr[i] + "\t" + chromPeakArr[i]);
            if(ionbasePeak<smoothchroIonArr[i]) {
                ionbasePeak = smoothchroIonArr[i];
                ionbasePeakIndex=i;
            }

           // System.out.println(smoothChromArr[i]);
        }
        */
        int identifiedIndex = maxWindow+margin-leftIndex;
        basePeak = smoothChromArr[identifiedIndex];
        ionbasePeak = smoothchroIonArr[identifiedIndex];
        /***********************************
         * Find simple/rough peak range (1/3 of base peak) for Gaussian input
         ***********************************/
        int[] indexResult = LabelfreeChroUtil.getPeakRange(identifiedIndex, basePeak, smoothChromArr);
        int peakStartIndex = indexResult[0];
        int peakEndIndex = indexResult[1];

        int[] indexIonResult = LabelfreeChroUtil.getPeakRange(identifiedIndex, ionbasePeak, smoothchroIonArr);
        int peakStartIonIndex = indexIonResult[0];
        int peakEndIonIndex = indexIonResult[1];


        /***********************************
         * Find peak range based on Gaussian fitting
        ***********************************/


        GaussianPeakModel gModel = GaussianFitting.getGaussianPeakRangeIndex2(peptide, retArr, smoothChromArr, peakStartIndex, peakEndIndex, ionArr,smoothchroIonArr,peakStartIonIndex,peakEndIonIndex);
        int finalPeakStartIndex=gModel.getPeakStartIndex()+leftIndex;
        //int finalPeakStartIndex=gModel.getPeakStartIndex();
        if(finalPeakStartIndex<0) finalPeakStartIndex=0;
        //int finalPeakEndIndex=gModel.getPeakEndIndex();
        int finalPeakEndIndex=gModel.getPeakEndIndex()+leftIndex;
        if(finalPeakEndIndex>=result.length)
            finalPeakEndIndex=result.length-1;

        double peakArea = gModel.getPeakArea();
        double peakAreaIon = gModel.getPeakAreaIonInjectionCorrection();

        /***********************************
         * build chro data
        ***********************************/
        StringBuffer sb = new StringBuffer();

      //  try {
      //  sb.append("P ").append((long)result[peakStart][0]).append(" ").append((long)result[peakEnd][0]).append(";");
        sb.append("P ").append(result[finalPeakStartIndex][0]).append(" ").append(result[finalPeakEndIndex][0]).
		append(" ").append(retArr[gModel.getPeakStartIndex()]).append(" ").append(retArr[gModel.getPeakEndIndex()]).append(";");
                //.append( ).append(";");
        /*    } catch(Exception e) {

            }*/

        count=1; //it starts with 1. why??

        Gaussian g = null;

        if(gModel.getX()>0)
            g = new Gaussian(gModel.getY(), gModel.getX(), gModel.getSigma());

        double gstart = -4 * gModel.getSigma() + gModel.getX();
        double gend = 4 * gModel.getSigma() + gModel.getX();

        List<Double> gxList = new ArrayList<>();
        List<Double> gyList = new ArrayList<>();

        for(int i=leftIndex+1;i<rightIndex-1;i++)
	{
            sb.append((int)result[i][0]).append(" "); //scan
            sb.append(retArr[count]).append(" "); //intensity
            sb.append((int)smoothChromArr[count]).append(" "); //intensity

      //  System.out.println("====\t" + retArr[count] + "\t" + smoothChromArr[count]);

            count++;
              /*
            for(int j=0;j<result[i].length-2;j++)
            {
                sb.append((long)result[i][j]).append(" ");
            }
            */
            sb.append((double)result[i][result[i].length-2]).append(" ");
            sb.append((double)result[i][result[i].length-1]).append(";");


           // if(null != g && retArr[count]>=gstart && retArr[count]<=gend) {
          /*  if(null != g) { // && retArr[count]>=gstart && retArr[count]<=gend) {
                gxList.add(retArr[count]);
                gyList.add(g.value(retArr[count]));
            }
        */

            //sb.append((int)result[i][0]).append(" ").append((long)result[i][1]).append(" ").append((long)result[i][2]).append(";");
	}




            // if(null != g && retArr[count]>=gstart && retArr[count]<=gend) {
            if(null != g) { // && retArr[count]>=gstart && retArr[count]<=gend) {
                for(double d:retArr) {
                    gxList.add(d);
                    gyList.add(g.value(d));
                }
            }



        if(null != g) {
            double[] gxArr = new double[gxList.size()];
            double[] gyArr = new double[gxList.size()];
            for(int j=0;j<gxList.size();j++) {
                gxArr[j] =gxList.get(j).doubleValue();
                gyArr[j] =gyList.get(j).doubleValue();

            }

            gModel.setGaussianXArr(gxArr);
            gModel.setGaussianYArr(gyArr);
        }

        gModel.setChroData(sb.toString());

        return gModel;

    }


    /*public static GaussianPeakModel buildGmodel(){

        GaussianPeakModel gModel = GaussianFitting.getGaussianPeakRangeIndex(retArr, smoothChromArr, peakStartIndex, peakEndIndex);
        int finalPeakStartIndex=gModel.getPeakStartIndex()+leftIndex;
        //int finalPeakStartIndex=gModel.getPeakStartIndex();
        if(finalPeakStartIndex<0) finalPeakStartIndex=0;
        //int finalPeakEndIndex=gModel.getPeakEndIndex();
        int finalPeakEndIndex=gModel.getPeakEndIndex()+leftIndex;
        if(finalPeakEndIndex>=result.length)
            finalPeakEndIndex=result.length-1;

        double peakArea = gModel.getPeakArea();


        StringBuffer sb = new StringBuffer();

      //  try {
      //  sb.append("P ").append((long)result[peakStart][0]).append(" ").append((long)result[peakEnd][0]).append(";");
        sb.append("P ").append(result[finalPeakStartIndex][0]).append(" ").append(result[finalPeakEndIndex][0]).
		append(" ").append(retArr[gModel.getPeakStartIndex()]).append(" ").append(retArr[gModel.getPeakEndIndex()]).append(";");
                //.append( ).append(";");
        /*    } catch(Exception e) {

                gModel = GaussianFitting.getGaussianPeakRangeIndex(retArr, smoothChromArr, peakStartIndex, peakEndIndex);
                e.printStackTrace();
            }

        count=1; //it starts with 1. why??

        Gaussian g = null;

        if(gModel.getX()>0)
            g = new Gaussian(gModel.getY(), gModel.getX(), gModel.getSigma());

        double gstart = -4 * gModel.getSigma() + gModel.getX();
        double gend = 4 * gModel.getSigma() + gModel.getX();

        List<Double> gxList = new ArrayList<>();
        List<Double> gyList = new ArrayList<>();

        for(int i=leftIndex+1;i<rightIndex-1;i++)
	{
            sb.append((int)result[i][0]).append(" "); //scan
            sb.append(retArr[count]).append(" "); //intensity
            sb.append((int)smoothChromArr[count]).append(" "); //intensity

      //  System.out.println("====\t" + retArr[count] + "\t" + smoothChromArr[count]);

            count++;
              /*
            for(int j=0;j<result[i].length-2;j++)
            {
                sb.append((long)result[i][j]).append(" ");
            }

            sb.append((double)result[i][result[i].length-2]).append(" ");
            sb.append((double)result[i][result[i].length-1]).append(";");


            if(null != g && retArr[count]>=gstart && retArr[count]<=gend) {
                gxList.add(retArr[count]);
                gyList.add(g.value(retArr[count]));
            }


            //sb.append((int)result[i][0]).append(" ").append((long)result[i][1]).append(" ").append((long)result[i][2]).append(";");
	}

        if(null != g) {
            double[] gxArr = new double[gxList.size()];
            double[] gyArr = new double[gxList.size()];
            for(int j=0;j<gxList.size();j++) {
                gxArr[j] =gxList.get(j).doubleValue();
                gyArr[j] =gyList.get(j).doubleValue();

            }

            gModel.setGaussianXArr(gxArr);
            gModel.setGaussianYArr(gyArr);
        }

        gModel.setChroData(sb.toString());

        return gModel;
    }*/

    //find peak resion based on 1/3 of base peak intensity
    //double[0] : start peak index
    //double[1]: end peak index
    public static int[] getPeakRange(int basePeakIndex, double basePeak, double[] smoothChromArr) {
                double peakThreshold = basePeak/3;

        int peakStartIndex=-1;
        int peakEndIndex=-1;

        for(int i=basePeakIndex;i<smoothChromArr.length;i++) {
            if(smoothChromArr[i]<peakThreshold) {
                peakEndIndex = i;
                break;
            }
        }

        if(peakEndIndex<0)
            peakEndIndex = smoothChromArr.length-1;

        for(int i=basePeakIndex;i>=0;i--) {
            if(smoothChromArr[i]<peakThreshold) {
                peakStartIndex = i;
                break;
            }
        }

        if(peakStartIndex<0)
            peakStartIndex = 0;

        int[] result = new int[2];
        result[0] = peakStartIndex;
        result[1] = peakEndIndex;

/*
int temp=0;
for(double d:smoothChromArr)
System.out.println("======\t" + temp++ + "\t" + d + " " + basePeakIndex);

System.out.println("==>>====\t" + result[0] + " " + result[1]);
System.exit(0);
      */
        return result;

    }

    public static Element createXmlChroPeptideTitle(Peptide peptide)
    {
        Element peptideEle = new Element("peptide");
        peptideEle.setAttribute("unique", peptide.isUnique()?"*":"");
        peptideEle.setAttribute("file", peptide.getFileName());
        peptideEle.setAttribute("scan", peptide.getScanNum());
        peptideEle.setAttribute("seq", peptide.getSequence());
        peptideEle.setAttribute("xcorr", peptide.getXCorr());
        peptideEle.setAttribute("calcMHplus", peptide.getCalcMHplus());
        peptideEle.setAttribute("MHplus", peptide.getMhPlus());
        peptideEle.setAttribute("totalIntensity", peptide.getTotalIntensity());
        peptideEle.setAttribute("spRank", peptide.getSpRank());
        peptideEle.setAttribute("spScore", peptide.getSpScore());
        peptideEle.setAttribute("redundancy", peptide.getRedundancy());
        peptideEle.setAttribute("deltaCN", peptide.getDeltCN());
        peptideEle.setAttribute("deltaMass", String.valueOf(peptide.getDeltaMass()));
	if(null == peptide.getDeltCN())
	    peptideEle.setAttribute("deltaCN", "");

        peptideEle.setAttribute("charge", peptide.getChargeState());
        peptideEle.setAttribute("spC", peptide.getRedundancy());


	Hashtable<String, String> scoreHt = peptide.getScoreHt();
	for(Iterator<String> itrScr=scoreHt.keySet().iterator(); itrScr.hasNext(); )
	{
	    String score = itrScr.next();
	    String value = scoreHt.get(score);

	    Element scoreEle = new Element("search_score");
	    scoreEle.setAttribute("name", score);
	    scoreEle.setAttribute("value", value);
	    peptideEle.addContent(scoreEle);
	}


	return peptideEle;
    }
    public static GaussianPeakModel calculateLabelfreeFullMS(Peptide peptide, int keyIndex, IndexedFile origIFile,
                                                             double[] samIsoArr)
            throws IOException, SQLException {
        SpectraDB spectraDB = origIFile.getSpectraDB();
        //peptide.getChargeState();
        int[] origIFileKeys = origIFile.getKeys();
        Configuration conf = Configuration.getInstance();
        double pepMass = Double.parseDouble(peptide.getCalcMHplus());


        double steepRatioThreshold = conf.getSteepRatioThreshold();

        int numIsoWindow = conf.getNumOfIsolationWindow();
        int maxWindow=conf.getMaxWindow();
        int margin = conf.getMargin();

        double[][] result = null;
        //conf.getQuantLevel()

        maxWindow = maxWindow * 3;
        result = new double[maxWindow*2+1+margin*2][3];

        int leftIndex=maxWindow+margin;
        int rightIndex=maxWindow+margin+1;


        //double leftTotalIntensity=0;
        //double rightTotalIntensity=0;
        double totalIntensity=0;

        int steepArea = conf.getSteepArea();

        int moveLeftKeyIndex = keyIndex;

        int moveRightKeyIndex = keyIndex+1*numIsoWindow;
        int identifiedScan = origIFileKeys[keyIndex];

        int initWin=2;

        for(int i=0;i<initWin;i++) ///Index=maxWindow;keys[leftPeakIndex]>=leftMinValue; )
        {
            if(moveLeftKeyIndex<=0 || leftIndex<=0)
            {
                moveLeftKeyIndex += 1*numIsoWindow;
                leftIndex++;

                break;
            }

            result[leftIndex] = readSplitFastFullSpectrum(moveLeftKeyIndex,samIsoArr, origIFileKeys, spectraDB, pepMass);

//            result[leftIndex] = readFullSpectrum(keys, moveLeftKeyIndex, index, file, samIsoArr, null);
            totalIntensity += result[leftIndex][1];

            // if(moveLeftKeyIndex<=0 || leftIndex<=0)
            //	break;

            moveLeftKeyIndex -= numIsoWindow;
            leftIndex--;
        }

        for(int i=0;i<initWin;i++) //leftIndex=maxWindow;keys[leftPeakIndex]>=leftMinValue; )
        {
            if(moveRightKeyIndex>=origIFileKeys.length)
            {
                moveRightKeyIndex -= 1*numIsoWindow;
                rightIndex--;

                break;
            }

            result[rightIndex] = readSplitFastFullSpectrum(moveRightKeyIndex, samIsoArr, origIFileKeys, spectraDB, pepMass);

            totalIntensity += result[rightIndex][1];

            moveRightKeyIndex += numIsoWindow;
            rightIndex++;

            if(totalIntensity<=0) break;

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
                        arr[steepCount] = result[rightIndex] = readSplitFastFullSpectrum(tempKeyIndex, samIsoArr, origIFileKeys, spectraDB, pepMass);
                    }

                    steepCount++;
                    if(arr[steepCount-1][1]<=0) break;
                    //tempKeyIndex--;
                    tempKeyIndex = tempKeyIndex - 1*numIsoWindow;
                }

                if(tempKeyIndex<0 || arr[steepCount][1]<=0  ) {
                    break;

                }

                arr[steepCount] = readSplitFastFullSpectrum(tempKeyIndex, samIsoArr, origIFileKeys, spectraDB, pepMass);

                if(arr[steepCount][1]<=0) break;

                for(int k=0;k<3;k++)
                {
                    for(int l=0;l<3;l++)
                        prevArr[k][l] = arr[k][l];
                }


                for(int j=leftIndex+1;j<=leftIndex+steepCount+1;j++)
                {
                    area1 += result[j][1];
                    area1 += result[j][2];
                }

                for(int j=0;j<=steepCount;j++)
                {
                    area2 += arr[j][1];
                    area2 += arr[j][2];
                }

                if(area2==0 && area1==0)
                    break;

                double areaRatio = (double)area2/area1;

                /*****************************************
                 * AREA1 is right side area
                 * AREA2 is left side area
                 *****************************************/
                if(isHighIntensity && ((double)totalIntensity/(rightIndex-leftIndex-1)/1.5>area2/3))
                    isHighIntensity = false;

                //if( (isGoingUp && (float)leftTotalIntensity/(rightIndex-leftIndex-1)/3<area2/3) || areaRatio<steepRatioThreshold )
                //if( (isGoingUp && (float)totalIntensity/(rightIndex-leftIndex-1)/3<area2/3) || areaRatio<steepRatioThreshold )
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

                if(moveRightKeyIndex+steepArea>=origIFileKeys.length)
                    break;

                int steepCount=0;

                int tempKeyIndex = moveRightKeyIndex;
                while(true)
                {
                    if(tempKeyIndex>=origIFileKeys.length || steepCount>= steepArea-1 || steepCount>=origIFileKeys.length)
                        break;

                    if(prevArr[0][0]!=0)
                    {
                        for(int l=0;l<2;l++)
                            for(int m=0;m<3;m++)
                                arr[l][m] = prevArr[l+1][m];
                    }
                    else
                    {
                        arr[steepCount] = readSplitFastFullSpectrum(tempKeyIndex, samIsoArr, origIFileKeys, spectraDB, pepMass);

                    }

                    steepCount++;
                    tempKeyIndex = tempKeyIndex + 1*numIsoWindow;
                    //tempKeyIndex++;
                }

                if(arr[steepCount-1][1]<=0) break;
                if(tempKeyIndex>=origIFileKeys.length)
                    break;

                arr[steepCount] = readSplitFastFullSpectrum(tempKeyIndex, samIsoArr, origIFileKeys, spectraDB, pepMass);
                if(arr[steepCount][1]<=0) break;

                for(int k=0;k<3;k++)
                {
                    for(int l=0;l<3;l++)
                        prevArr[k][l] = arr[k][l];
                }

                for(int j=rightIndex-1;j>rightIndex-steepCount-2;j--)
                {
                    area1 += result[j][1];
                    area1 += result[j][2];
                }

                for(int j=0;j<=steepCount;j++)
                {
                    area2 += arr[j][1];
                    area2 += arr[j][2];
                }

                if(area2==0 && area1==0)
                    break;

                double areaRatio = (double)area2/area1;

                /*****************************************
                 * AREA1 is left side area
                 * AREA2 is right side area
                 *****************************************/

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

        if(leftIndex<0) leftIndex=0;
        if(rightIndex>=result.length) rightIndex=result.length-1;

        int tmpLeftIndex = leftIndex;
        int tmpRightIndex = rightIndex;
        int tmpMoveLeftKeyIndex = moveLeftKeyIndex;
        int tmpMoveRightKeyIndex = moveRightKeyIndex;

        for(int i=0;i<100;i++) ///Index=maxWindow;keys[leftPeakIndex]>=leftMinValue; )
        {
            if(tmpMoveLeftKeyIndex<0 || tmpLeftIndex<0)
            {
                tmpMoveLeftKeyIndex += 1*numIsoWindow;
                tmpLeftIndex++;

                break;
            }

            result[tmpLeftIndex] = readSplitFastFullSpectrum(tmpMoveLeftKeyIndex,samIsoArr, origIFileKeys, spectraDB, pepMass);

            tmpMoveLeftKeyIndex -= numIsoWindow;
            tmpLeftIndex--;

            if(result[tmpLeftIndex+1][1]<=0) {

                break;
            }
        }

        for(int i=0;i<100;i++) //leftIndex=maxWindow;keys[leftPeakIndex]>=leftMinValue; )
        {
            if(tmpMoveRightKeyIndex>=origIFileKeys.length || tmpRightIndex>=result.length)
            {
                tmpMoveRightKeyIndex -= 1*numIsoWindow;
                tmpRightIndex--;

                break;
            }

            result[tmpRightIndex] = readSplitFastFullSpectrum(tmpMoveRightKeyIndex, samIsoArr, origIFileKeys, spectraDB, pepMass);

            tmpMoveRightKeyIndex += numIsoWindow;
            tmpRightIndex++;



            if(result[tmpRightIndex-1][1]<=0) break;

        }




        if(moveRightKeyIndex>=origIFileKeys.length)
        {
            moveRightKeyIndex -= 1*numIsoWindow;
            rightIndex--;
        }

        result[rightIndex] = readSplitFastFullSpectrum(moveRightKeyIndex, samIsoArr, origIFileKeys, spectraDB, pepMass);


        result[leftIndex] = readSplitFastFullSpectrum(moveLeftKeyIndex, samIsoArr, origIFileKeys, spectraDB, pepMass);


        long[] chromPeakArr = new long[tmpRightIndex-tmpLeftIndex+1];


        int count=0;

        TIntDoubleHashMap scanRtMap = origIFile.getScanRtMap();
        TIntDoubleHashMap scanIonMap = origIFile.getIonInjectionMap();

        double[] retArr = new double[chromPeakArr.length];
        double identifiedRt = scanRtMap.get(identifiedScan);
        double[] ionArr = new double[chromPeakArr.length];


        //for(int i=leftIndex;i<=rightIndex;i++) {
        for(int i=tmpLeftIndex;i<=tmpRightIndex;i++) {
            try {
                chromPeakArr[count] = (long)result[i][1];
                retArr[count] = scanRtMap.get((int)result[i][0]);
                count++;
            } catch(Exception e) {

                e.printStackTrace();
            }
        }
        count = 0;
        //for(int i=leftIndex;i<=rightIndex;i++) {
        for(int i=tmpLeftIndex;i<=tmpRightIndex;i++) {
            try {
                ionArr[count++] = scanIonMap.get((int)result[i][0]);
                //count++;
            } catch(Exception e) {
                System.out.println(ionArr.length + "\t" + count + "\t" +result.length + "\t" + i);

                e.printStackTrace();
            }
        }

        long[] chroIonArr = new long[chromPeakArr.length];
        for(int i =0;i<chromPeakArr.length;i++){
            if(i==0){
                chroIonArr[i] =  (long) ((chromPeakArr[i]/ionArr[i])*(retArr[i+1]-retArr[i]));
            }else if(i==chromPeakArr.length-1){
                chroIonArr[i] =  (long) ((chromPeakArr[i]/ionArr[i])*(retArr[i]-retArr[i-1]));
            }
            else{
                chroIonArr[i] =  (long) ((chromPeakArr[i]/ionArr[i])*(((retArr[i]-retArr[i-1])/2)+((retArr[i+1]-retArr[i])/2)));
            }
        }

        double[] smoothchroIonArr = Smooth.smoothAsDouble(chroIonArr, LabelfreeMissingPeptideBuilderSplit.SMOOTH_WINDOW_SIZE);
        double[] smoothChromArr = Smooth.smoothAsDouble(chromPeakArr, LabelfreeMissingPeptideBuilderSplit.SMOOTH_WINDOW_SIZE);

        //for(double d:chromPeakArr)
        //  System.out.println("===\t" + d);

        double ionbasePeak=0;
        double basePeak=0;
        int basePeakIndex=0;

        int identifiedIndex = maxWindow+margin-leftIndex;
        if(identifiedIndex<0) identifiedIndex=0;
        basePeak = smoothChromArr[identifiedIndex];
        ionbasePeak = smoothchroIonArr[identifiedIndex];
        /***********************************
         * Find simple/rough peak range (1/3 of base peak) for Gaussian input
         ***********************************/
        int[] indexResult = LabelfreeChroUtil.getPeakRange(identifiedIndex, basePeak, smoothChromArr);
        int peakStartIndex = indexResult[0];
        int peakEndIndex = indexResult[1];

        int[] indexIonResult = LabelfreeChroUtil.getPeakRange(identifiedIndex, ionbasePeak, smoothchroIonArr);
        int peakStartIonIndex = indexIonResult[0];
        int peakEndIonIndex = indexIonResult[1];


        /***********************************
         * Find peak range based on Gaussian fitting
         ***********************************/
//System.out.println("==================");
        GaussianPeakModel gModel = GaussianFitting.getGaussianPeakRangeIndex2(
                peptide, retArr, smoothChromArr, peakStartIndex, peakEndIndex, ionArr,smoothchroIonArr,peakStartIonIndex,peakEndIonIndex);

        double maxPeakIntensity = 0;
        for(double d:smoothChromArr) {
            if(d>maxPeakIntensity)
                maxPeakIntensity = d;
        }

        gModel.setMaxIntensity(maxPeakIntensity);

        //backGroundNoise = CalcUtilGeneric.getBackGroundNoise(isoArr, currentPos, diff, splitIFile);


        int finalPeakStartIndex=gModel.getPeakStartIndex()+leftIndex;
        //int finalPeakStartIndex=gModel.getPeakStartIndex();
        if(finalPeakStartIndex<0) finalPeakStartIndex=0;
        //int finalPeakEndIndex=gModel.getPeakEndIndex();
        int finalPeakEndIndex=gModel.getPeakEndIndex()+leftIndex;
        if(finalPeakEndIndex>=result.length)
            finalPeakEndIndex=result.length-1;

        double peakArea = gModel.getPeakArea();

        double peakAreaIon = gModel.getPeakAreaIonInjectionCorrection();

        /***********************************
         * build chro data
         ***********************************/
        StringBuffer sb = new StringBuffer();

        sb.append("P ").append(result[finalPeakStartIndex][0]).append(" ").append(result[finalPeakEndIndex][0]).
                append(" ").append(retArr[gModel.getPeakStartIndex()]).append(" ").append(retArr[gModel.getPeakEndIndex()]).append(";");



        Gaussian g = null;

        if(gModel.getX()>0)
            g = new Gaussian(gModel.getY(), gModel.getX(), gModel.getSigma());

        double gstart = -4 * gModel.getSigma() + gModel.getX();
        double gend = 4 * gModel.getSigma() + gModel.getX();

        List<Double> gxList = new ArrayList<>();
        List<Double> gyList = new ArrayList<>();

        count=1; //it starts with 1. why??

        for(int i=0;i<smoothChromArr.length;i++)
        {
            if(retArr[i]<=0) continue;
            sb.append((int)result[i+tmpLeftIndex][0]).append(" "); //scan
            sb.append(retArr[i]).append(" "); //ret time
            sb.append((int)smoothChromArr[i]).append(" "); //intensity

            count++;
            sb.append(result[i][result[i].length-2]).append(" ");
            sb.append(result[i][result[i].length-1]).append(";");

        }

        if(null != g) {


            double[] allRtArr = origIFile.getRtArr();
            int rtStartkeyIndex = Arrays.binarySearch(allRtArr, gstart);
            int rtEndkeyIndex = Arrays.binarySearch(allRtArr, gend);

            if (rtStartkeyIndex < 0) //Cannot find index
            {
                rtStartkeyIndex = -(++rtStartkeyIndex); //Math.abs(++keyIndex);
            }
            if (rtStartkeyIndex >= allRtArr.length) {
                rtStartkeyIndex--;
            }

            if (rtEndkeyIndex < 0) //Cannot find index
            {
                rtEndkeyIndex = -(++rtEndkeyIndex); //Math.abs(++keyIndex);
            }
            if (rtEndkeyIndex >= allRtArr.length) {
                rtEndkeyIndex--;
            }

            for(int i=rtStartkeyIndex;i<=rtEndkeyIndex;i++) {
                gxList.add(allRtArr[i]);
                gyList.add(g.value(allRtArr[i]));
            }

        }



        if(null != g) {
            double[] gxArr = new double[gxList.size()];
            double[] gyArr = new double[gxList.size()];
            for(int j=0;j<gxList.size();j++) {
                gxArr[j] = gxList.get(j).doubleValue();
                gyArr[j] = gyList.get(j).doubleValue();
            }

            gModel.setGaussianXArr(gxArr);
            gModel.setGaussianYArr(gyArr);
        }

        gModel.setChroData(sb.toString());

        return gModel;
    }


}
