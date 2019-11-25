/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.scripps.pms.census.labelFree;

import edu.scripps.pms.census.conf.Configuration;
import edu.scripps.pms.census.exception.CensusIndexOutOfBoundException;
import edu.scripps.pms.census.hash.IndexedFile;
import edu.scripps.pms.census.model.NonLabelMappingModel;
import edu.scripps.pms.census.tools.IsotopeModel;
import edu.scripps.pms.census.util.CalcUtil;
import edu.scripps.pms.census.util.io.MzxmlSpectrumReader;
import gnu.trove.TIntLongHashMap;
import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author harshil
 */
public class LabelFreeCalcUtil extends CalcUtil{

    private static final char CARRIAGE_RETURN = '\n';


    public LabelFreeCalcUtil() {
        
    }

    
    public static void main(String args[])
    {
        
        
        IndexedFile iFile = null;
        double[] mass = {805,806,807};
        File file =null;
        try {
            file =  new File("/home/harshil/data/jolene/projects2014_03_03_23_58406/20140227_MP23_HEK_PBS_1ug_BEH_50cm_140min_35ms_CID_1.ms2.index");
            iFile = new IndexedFile(file,"/home/harshil/data/jolene/projects2014_03_03_23_58406/20140227_MP23_HEK_PBS_1ug_BEH_50cm_140min_35ms_CID_1.ms2");
            System.out.println(LabelFreeCalcUtil.peakFindingIso(iFile,mass, 450));
            System.out.println("");
        } catch (IOException ex) {
            Logger.getLogger(LabelFreeCalcUtil.class.getName()).log(Level.SEVERE, null, ex);
        } catch (Exception ex) {
            Logger.getLogger(LabelFreeCalcUtil.class.getName()).log(Level.SEVERE, null, ex);
        }
        
    }
    
    public static PeakRange peakFindingIso(IndexedFile indexFile, double[] samIsoArr, int keyIndex)
            throws IOException, CensusIndexOutOfBoundException, Exception {

//	TIntLongHashMap index = sModel.getIndex();
        TIntLongHashMap index = indexFile.getMsIndex();
//        Object file = sModel.getGenericIndexFile(conf);
        Object file = indexFile.getFile();
        if(conf == null)
            conf = Configuration.getInstance();
        steepRatioThreshold = conf.getSteepRatioThreshold();

        int numIsoWindow = conf.getNumOfIsolationWindow();
        int maxWindow = conf.getMaxWindow();
        int margin = conf.getMargin();

        double[][] result = null;
        //conf.getQuantLevel()

//	if(conf.isDataIndependent())
//	    result = new double[maxWindow*2+1+margin*2][4*sModel.getBioSample().length+3]; //scan #, sample intensity, ref intensity
//	else
        result = new double[maxWindow * 2 + 1 + margin * 2][3];

        int leftIndex = maxWindow + margin;
        int rightIndex = maxWindow + margin + 1;

        //double leftTotalIntensity=0;
        //double rightTotalIntensity=0;
        double totalIntensity = 0;

        int steepArea = conf.getSteepArea();

        int moveLeftKeyIndex = keyIndex;

        int moveRightKeyIndex = keyIndex + 1 * numIsoWindow;

//        int[] keys = sModel.getKeys();
        int[] keys = indexFile.getKeys();

        int initWin = 2;

        for (int i = 0; i < initWin; i++) ///Index=maxWindow;keys[leftPeakIndex]>=leftMinValue; )
        {
            if (moveLeftKeyIndex <= 0 || leftIndex <= 0) {
                moveLeftKeyIndex += 1 * numIsoWindow;
                leftIndex++;

                break;
            }

//	    if(conf.isDataIndependent()) //data independent
//	    {
//		result[leftIndex] = readSpectrum(keys, moveLeftKeyIndex, index, sModel.getDiff(), file, sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());
//
//	    }	
//	    else //data dependent
//	    {
//		if(sModel.isHighRes())
            result[leftIndex] = readFullSpectrum(keys, moveLeftKeyIndex, index, file, samIsoArr);
//		else
//		    result[leftIndex] = readFullSpectrum(keys, moveLeftKeyIndex, index, file, sModel.getSamStartMass(), sModel.getSamEndMass(), sModel.getRefStartMass(), sModel.getRefEndMass());

//	    }
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

//	    if(conf.isDataIndependent()) //data independent
//	    {
//		result[rightIndex] = readSpectrum(keys, moveRightKeyIndex, index, sModel.getDiff(), file, sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());
//	    }
//	    else //data dependent
//	    {
//		if(sModel.isHighRes())
            result[rightIndex] = readFullSpectrum(keys, moveRightKeyIndex, index, file, samIsoArr);
//		else
//		    result[rightIndex] = readFullSpectrum(keys, moveRightKeyIndex, index, file, sModel.getSamStartMass(), sModel.getSamEndMass(), sModel.getRefStartMass(), sModel.getRefEndMass());
//	    }

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
                    if (tempKeyIndex < 0 || steepCount < 0 || steepCount >= steepArea - 1) {
                        break;
                    }

                    if (prevArr[0][0] != 0) {
                        for (int l = 0; l < 2; l++) {
                            for (int m = 0; m < 3; m++) {
                                arr[l][m] = prevArr[l + 1][m];
                            }
                        }
                    } else {
//			if(conf.isDataIndependent()) //data independent
//			    arr[steepCount] = readSpectrum(keys, tempKeyIndex, index, sModel.getDiff(), file, sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());
//			else
//			{
//			    if(sModel.isHighRes())
                        arr[steepCount] = readFullSpectrum(keys, tempKeyIndex, index, file, samIsoArr);
//			    else
//				arr[steepCount] = readFullSpectrum(keys, tempKeyIndex, index, file, sModel.getSamStartMass(), sModel.getSamEndMass(), sModel.getRefStartMass(), sModel.getRefEndMass());
//			}            
                    }

                    steepCount++;
                    //tempKeyIndex--;
                    tempKeyIndex = tempKeyIndex - 1 * numIsoWindow;
                }

                if (tempKeyIndex < 0) {
                    break;
                }

//		if(conf.isDataIndependent()) //data independent
//		    arr[steepCount] = readSpectrum(keys, tempKeyIndex, index, sModel.getDiff(), file, sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());
//		else
//		{
//		    if(sModel.isHighRes())
                arr[steepCount] = readFullSpectrum(keys, tempKeyIndex, index, file, samIsoArr);
//		    else
//			arr[steepCount] = readFullSpectrum(keys, tempKeyIndex, index, file, sModel.getSamStartMass(), sModel.getSamEndMass(), sModel.getRefStartMass(), sModel.getRefEndMass());
//		}

                for (int k = 0; k < 3; k++) {
                    for (int l = 0; l < 3; l++) {
                        prevArr[k][l] = arr[k][l];
                    }
                }

                for (int j = leftIndex + 1; j <= leftIndex + steepCount + 1; j++) {
                    area1 += result[j][1];
                    area1 += result[j][2];
                }

                for (int j = 0; j <= steepCount; j++) {
                    area2 += arr[j][1];
                    area2 += arr[j][2];
                }

                if (area2 == 0 && area1 == 0) {
                    break;
                }

                double areaRatio = (double) area2 / area1;

                /**
                 * ***************************************
                 * AREA1 is right side area AREA2 is left side area
                 ****************************************
                 */
                if (isHighIntensity && ((double) totalIntensity / (rightIndex - leftIndex - 1) / 1.5 > area2 / 3)) {
                    isHighIntensity = false;
                }

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

                    if (0 == result[leftIndex][0]) {
                        continue;
                    }

                    leftIndex--;

                    for (int j = 0; j < steepCount; j++) {
                        arr[j] = arr[j + 1];
                    }

                    totalIntensity += arr[0][1];
                    totalIntensity += arr[0][2];

                    continue;
                } else {
                    isGoingUp = false;
                }

                if (areaRatio > steepRatioThreshold) {
                    for (int k = 0; k < 3; k++) {
                        if (result[leftIndex + 1][1] < arr[k][1] || arr[k][0] == 0) {
                            break;
                        }

                        if (leftIndex < 0) {
                            break;
                        }

                        result[leftIndex] = arr[k];

                        moveLeftKeyIndex -= numIsoWindow;
                        leftIndex--;
                    }

                    break;
                }

                result[leftIndex--] = arr[0];

                moveLeftKeyIndex -= numIsoWindow;

                for (int j = 0; j < steepCount; j++) {
                    arr[j] = arr[j + 1];
                }
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

                if (moveRightKeyIndex + steepArea >= keys.length) {
                    break;
                }

                int steepCount = 0;

                int tempKeyIndex = moveRightKeyIndex;
                while (true) {
                    if (tempKeyIndex >= keys.length || steepCount >= steepArea - 1 || steepCount >= keys.length) {
                        break;
                    }

                    if (prevArr[0][0] != 0) {
                        for (int l = 0; l < 2; l++) {
                            for (int m = 0; m < 3; m++) {
                                arr[l][m] = prevArr[l + 1][m];
                            }
                        }
                    } else {
//			if(conf.isDataIndependent()) //data independent
//			    arr[steepCount] = readSpectrum(keys, tempKeyIndex, index, sModel.getDiff(), file, sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());
//			else if(sModel.isHighRes())
                        arr[steepCount] = readFullSpectrum(keys, tempKeyIndex, index, file, samIsoArr);
//			else
//			    arr[steepCount] = readFullSpectrum(keys, tempKeyIndex, index, file, sModel.getSamStartMass(), sModel.getSamEndMass(), sModel.getRefStartMass(), sModel.getRefEndMass());
                    }

                    steepCount++;
                    tempKeyIndex = tempKeyIndex + 1 * numIsoWindow;
                    //tempKeyIndex++;
                }

                if (tempKeyIndex >= keys.length) {
                    break;
                }

//		if(conf.isDataIndependent()) //data independent
//		    arr[steepCount] = readSpectrum(keys, tempKeyIndex, index, sModel.getDiff(), file, sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());
//		else if(sModel.isHighRes())
                arr[steepCount] = readFullSpectrum(keys, tempKeyIndex, index, file, samIsoArr);
//		else
//		    arr[steepCount] = readFullSpectrum(keys, tempKeyIndex, index, file, sModel.getSamStartMass(), sModel.getSamEndMass(), sModel.getRefStartMass(), sModel.getRefEndMass());

                for (int k = 0; k < 3; k++) {
                    for (int l = 0; l < 3; l++) {
                        prevArr[k][l] = arr[k][l];
                    }
                }

                for (int j = rightIndex - 1; j > rightIndex - steepCount - 2; j--) {
                    area1 += result[j][1];
                    area1 += result[j][2];
                }

                for (int j = 0; j <= steepCount; j++) {
                    area2 += arr[j][1];
                    area2 += arr[j][2];
                }

                if (area2 == 0 && area1 == 0) {
                    break;
                }

                double areaRatio = (double) area2 / area1;

                /**
                 * ***************************************
                 * AREA1 is left side area AREA2 is right side area
                 ****************************************
                 */
                if (isHighIntensity && ((double) totalIntensity / (rightIndex - leftIndex - 1) / 1.5 > area2 / 3)) {
                    isHighIntensity = false;
                }

                //if( ((isGoingUp && (double)totalIntensity/(rightIndex-leftIndex-1)/3<area2/3)) || areaRatio<steepRatioThreshold )                 
                if ((isGoingUp && isHighIntensity) || areaRatio < steepRatioThreshold) {
                    if ((rightIndex + 1) >= result.length) {
                        isGoingUp = false;
                        break;
                    }

                    result[rightIndex] = arr[0];
                    moveRightKeyIndex += numIsoWindow;

                    if (0 == result[rightIndex][0]) {
                        continue;
                    }

                    rightIndex++;

                    for (int j = 0; j < steepCount; j++) {
                        arr[j] = arr[j + 1];
                    }

                    //rightTotalIntensity += arr[0][1];
                    totalIntensity += arr[0][1];
                    totalIntensity += arr[0][2];

                    continue;
                } else {
                    isGoingUp = false;
                }

                if (area2 / area1 > steepRatioThreshold) {
                    for (int k = 0; k < 3; k++) {
                        if (result[rightIndex - 1][1] < arr[k][1]) {
                            break;
                        }

                        if (result.length <= rightIndex) {
                            break;
                        }

                        result[rightIndex] = arr[k];

                        //System.out.println(arr[0][1] + " " + arr[1][1] + " " + arr[2][1] + " " + result[rightIndex-1][1]);
                        moveRightKeyIndex += numIsoWindow;
                        rightIndex++;
                    }
                    break;
                }

                result[rightIndex++] = arr[0];

                moveRightKeyIndex += numIsoWindow;

                for (int j = 0; j < steepCount; j++) {
                    arr[j] = arr[j + 1];
                }
            }
        }
        int peakStart = leftIndex + 1;
        int peakEnd = rightIndex - 1;

//        int dtaStart = (range == null ? 0 : range.getMin());
//        int dtaEnd = (range == null ? 0 : range.getMax());
        int specSpace = (int) (result[peakStart + 1][0] - result[peakStart][0]);

//System.out.println( result[peakStart+1][0] + " " + result[peakStart][0] );
        //int diff = (int) (result[peakStart][0] - dtaStart) / specSpace;
        int leftMargin = margin;

        int rightMargin = margin;
//        if (diff > 0) {
//            leftMargin += diff;
//        }

//        diff = (int) (dtaEnd - result[rightIndex - 1][0]) / specSpace;
//
//        if (diff > 0) {
//            rightMargin += diff;
//        }

        for (int i = 0; i < leftMargin; i++) {
            if (leftIndex < 0 || moveLeftKeyIndex < 0) {
                break;
            }

//	    if(conf.isDataIndependent()) //data independent
//		result[leftIndex] = readSpectrum(keys, moveLeftKeyIndex, index, sModel.getDiff(), file, sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());
//	    else if(sModel.isHighRes())
            result[leftIndex] = readFullSpectrum(keys, moveLeftKeyIndex, index, file, samIsoArr);
//	    else
//		result[leftIndex] = readFullSpectrum(keys, moveLeftKeyIndex, index, file, sModel.getSamStartMass(), sModel.getSamEndMass(), sModel.getRefStartMass(), sModel.getRefEndMass());

            if (0 == result[leftIndex][0]) {
                moveLeftKeyIndex -= numIsoWindow;
                continue;
            }

            leftIndex--;
            moveLeftKeyIndex -= numIsoWindow;
        }
        for (int i = 0; i < rightMargin; i++) {
            if (rightIndex >= result.length || moveRightKeyIndex >= keys.length) {
                break;
            }

//	    if(conf.isDataIndependent()) //data independent
//		result[rightIndex] = readSpectrum(keys, moveRightKeyIndex, index, sModel.getDiff(), file, sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());
//	    else if(sModel.isHighRes())
            result[rightIndex] = readFullSpectrum(keys, moveRightKeyIndex, index, file, samIsoArr);
//	    else
//		result[rightIndex] = readFullSpectrum(keys, moveRightKeyIndex, index, file, sModel.getSamStartMass(), sModel.getSamEndMass(), sModel.getRefStartMass(), sModel.getRefEndMass());

            if (0 == result[rightIndex][0]) {
                moveRightKeyIndex += numIsoWindow;
                continue;
            }

            rightIndex++;
            moveRightKeyIndex += numIsoWindow;
        }
//	if(conf.isDataIndependent())
//	    return buildDIResult(peakStart, peakEnd, leftIndex, rightIndex, result, sModel.getBioSample().length);
//	else
        return new PeakRange((int)result[peakStart][0],(int)result[peakEnd][0]);
//        return buildDDResult(peakStart, peakEnd, leftIndex, rightIndex, result);
    }

    
    /**
     * 
     * High resolution AND lalebling data 
     * @param keys
     * @param curIndex
     * @param index
     * @param reader
     * @param samIsoArr
     * @return
     * @throws IOException
     * @throws CensusIndexOutOfBoundException
     * @throws Exception 
     */
    public static double[] readFullSpectrum(
            int[] keys, 
            int curIndex, 
            TIntLongHashMap index, 
            //RandomAccessFile file,
            Object reader,
            double[] samIsoArr	    
            ) throws IOException, CensusIndexOutOfBoundException, Exception
    {                       

	double result[] = new double[5];
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
	double[] tempArr = LabelFreeCalcUtil.intensitySum(massArr, intArr, samIsoArr);        
        result[1] = tempArr[0];     //intensity
        result[2] = samIsoArr.length;
        result[3] = tempArr[1];      //Found iso num 
        result[4] = tempArr[2];	     //tolerance
               
	
      
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
    
    
    public static edu.scripps.pms.census.labelFree.SpectrumModel readSpectrumPeak(IndexedFile iFile, int keyIndex, int[] keys) throws IOException
    {

        //double[] result = new double[2];
	RandomAccessFile rfile = iFile.getFile();
	TIntLongHashMap indexMap = iFile.getMsIndex();
	long startPos = indexMap.get(keys[keyIndex]);

        int byteSize=-1;
        byte[] bytes;

        char ch;
        int pos;

        rfile.seek(startPos);

	//if( (keyIndex+1)>=keys.length )
	if( (keyIndex)>=keys.length )
	    byteSize = (int)(rfile.length() - startPos); 

	if(indexMap.get(keys[keyIndex+1])>0)
    //    else
	    byteSize = (int)(indexMap.get(keys[keyIndex+1]) - startPos);    

        bytes = new byte[byteSize];

        rfile.readFully(bytes);

        pos=0;

        ch = (char)bytes[pos];
        int scanNumner = keyIndex;
        double retentionTime = 0.00;
        double ionInjectionTime = 0.00;
        //        pos++;
        //Remove Z, S, I, D lines
        //System.out.println("sample");
        while( (ch=(char)bytes[pos]) == 'Z' || ch == 'S' || ch == 'I' || ch == 'D')
        {
            if((ch == 'S' || ch == 'I' ) && bytes[pos+1] == '\t' )
            {
                StringBuffer sb = new StringBuffer();
                for(int i = pos; bytes[i] != CARRIAGE_RETURN ;i++)
                {
                    sb.append((char) bytes[i]);
                }
                String[] arr = sb.toString().split("\t");
                if(arr[1].equals("RetTime"))
                    retentionTime = Double.parseDouble(arr[2]);
                if(arr[1].equals("IonInjectionTime"))
                    ionInjectionTime = Double.parseDouble(arr[2]);
                if(arr[0].equalsIgnoreCase("S"))
                    scanNumner = Integer.parseInt(arr[1]);
            }
            while( ch != CARRIAGE_RETURN )
            {
                pos++;
                ch = (char)bytes[pos];
            }       

            pos++;
        }

        int arrSize=0;
        for(int j=pos;j<byteSize;j++)
        {
            if( CARRIAGE_RETURN == (char)bytes[j] )
                arrSize++;

        }
        
        double[][] result = parseSpectra(bytes);
        double[] massArr = result[0];
        double[] intArr = result[1];

        //pos++;
//        StringBuilder mass = new StringBuilder(10);
//        StringBuilder intensity = new StringBuilder(10);
//	intensity.append('0');

	return new edu.scripps.pms.census.labelFree.SpectrumModel(massArr, intArr,retentionTime,ionInjectionTime,scanNumner);
    }
    
    
    /**
     * pass the ScanNumebr which is a PEAK ScanNumber ..
     * 
     *      
     * @param iFile
     * @param massArr
     * @param scanNumber (peak) 
     * @param isoTopeModelList
     * @return  PeakRange object thats is start and end range.
     * @throws CensusIndexOutOfBoundException
     * @throws Exception 
     */
    public static PeakRange generateRange(IndexedFile iFile,double[] massArr, int scanNumber,List<IsotopeModel> isoTopeModelList) throws CensusIndexOutOfBoundException, Exception
    {
        
        int[] keys = iFile.getKeys();
        int keyIndex = Arrays.binarySearch(keys, scanNumber);
        if (keyIndex < 0) //Cannot find index
        {
            keyIndex = -(++keyIndex); //Math.abs(++keyIndex);
        }
        if (keyIndex >= keys.length) {
            keyIndex--;
        }
        if (isoTopeModelList.size() <= 0) {
            return null;
        }
        double[] massList = isoTopeModelList.get(0).getIsoArr();
        
        
//        return LabelFreeCalcUtil.peakFindingIso(iFile, massArr,keyIndex );
          return LabelFreeCalcUtil.peakFindingHarshil(iFile, massArr,keyIndex );
        
    }
    
    
    public static PeakRange peakFinding(
            IndexedFile indexFile, double[] samIsoArr, int keyIndex)
//            CalcUtil.SpectrumModel sModel, SpecRange range, Configuration conf, int keyIndex, double[] samIsoArr)            
        throws IOException, Exception
    {
        if(conf == null)
            conf = Configuration.getInstance();
        steepRatioThreshold = conf.getSteepRatioThreshold();
    
        int maxWindow=conf.getMaxWindow();
        int margin = conf.getMargin();
        
        NonLabelMappingModel mapModel = conf.getMapModel();        
        //Hashtable<String, Hashtable> ms2ms1Ht = mapModel.getMs2ms1Ht();
                
       // System.out.println(mapModel.getMs2ms1Ht());
       // System.exit(0);
        
//        Vector<String> pathFileNameList = mapModel.getPathFileNameList();

        //ResultList result = new ResultList(pathArray.length+1, maxWindow, margin);
//        int resultSize = pathFileNameList.size()*2;
//        int resultIntensityStart = pathFileNameList.size();
        int resultSize = 3;
        int resultIntensityStart = 1;
        
        double[][] result = null;
//        result = new double[maxWindow * 2 + 1 + margin * 2][3];
        result = new double[maxWindow*2+1+margin*2][resultSize];
        
//        int maxScanIndex = mapModel.getMaxScanIndex(pathFileName);
        int maxScanIndex = indexFile.getKeys()[indexFile.getKeys().length-1];

	if(maxScanIndex<0)
	{
	    System.out.println("Error : Spectral file name in the config file is invalid : " );
	    throw new Exception();
	}

	
        //int maxScanIndex = result.length;
        
        int leftIndex=maxWindow+margin;
        int rightIndex=maxWindow+margin+1;

        double totalIntensity=0;

        int steepArea = conf.getSteepArea(); 
	
        int moveLeftKeyIndex = keyIndex;
        int moveRightKeyIndex = -1;
        //int[] keys = sModel.getKeys();
        

        moveRightKeyIndex = keyIndex+1;

        int initWin=2;
	Hashtable<Integer, double[]> resultHt = new Hashtable<Integer, double[]>();
        
        for(int i=0;i<initWin;i++)
        {
	    if(moveLeftKeyIndex<=0 || leftIndex<=0)
	    {
		    moveLeftKeyIndex++;                

		leftIndex++;
		
		break;
	    }

            
                result[leftIndex] = readFullSpectrum(moveLeftKeyIndex, null, null, samIsoArr, mapModel, null, resultHt);//
            
           // System.exit(0);
                
            //result = new double[maxWindow*2+1+margin*2][4*sModel.getBioSample().length+resultSize]; //scan #, sample intensity, ref intensity

            for(int j=resultIntensityStart;j<result[leftIndex].length;j++)
                totalIntensity += result[leftIndex][j];
            
		moveLeftKeyIndex--;                

            //if it is null, continue without changing leftIndex value;
            if(0 == result[leftIndex][0])
            {
                i--;
                continue;
            }
	    
            leftIndex--;
        }

        for(int i=0;i<initWin;i++) //leftIndex=maxWindow;keys[leftPeakIndex]>=leftMinValue; )
        {
            if(moveRightKeyIndex>maxScanIndex || moveRightKeyIndex<0)
            {
		    moveRightKeyIndex--;                

                rightIndex--;

                break;
            }
            
                result[rightIndex] = readFullSpectrum(moveRightKeyIndex, null, null, samIsoArr, mapModel, null, resultHt);

            for(int j=resultIntensityStart;j<result[leftIndex].length;j++)
                totalIntensity += result[leftIndex][j];

		moveRightKeyIndex++;                
            
            //if it is null, continue without changing leftIndex value;
            if(0 == result[rightIndex][0])
            {
                i--;
                continue;
            }
            
            rightIndex++;
        }

	if(moveRightKeyIndex<0 && moveLeftKeyIndex<0)
	    return null;

        boolean isGoingUp=true;
        boolean isHighIntensity=true;  //if the intensity is lower than one third of average intensity, this becomes false
        double[][] arr = new double[steepArea][resultSize];
        double[][] prevArr = new double[steepArea][resultSize];
	Hashtable<Integer, double[]> tmpKeyHt = new Hashtable<Integer, double[]>();

        if(leftIndex+steepArea>=0)
        {
            //for(int i=0;i<300;i++)
	    while(true)
            {            
                double area1=0;
                double area2=0;

                if(moveLeftKeyIndex-steepArea<0)                    
                    break;                

                int steepCount=0;
                int tempKeyIndex = moveLeftKeyIndex;
                while(true)
                {
		    if(tempKeyIndex<0 || steepCount<0 || steepCount>=steepArea-1)
			break;
                    
		    if(prevArr[0][0]!=0)
		    {
			for(int l=0;l<2;l++)
			    for(int m=0;m<resultSize;m++)
				arr[l][m] = prevArr[l+1][m];
		    }
		    else
		    {
			
			
			    double[] tmpKeyArr = tmpKeyHt.get(tempKeyIndex);
			    if(null != tmpKeyArr)
				arr[steepCount] = tmpKeyArr;
			    else
			    {
				    arr[steepCount] = readFullSpectrum(tempKeyIndex, null, null, samIsoArr, mapModel, null, resultHt);
			    
				tmpKeyHt.put(tempKeyIndex, arr[steepCount]);
			    }
		    }
		    
                    steepCount++;
			tempKeyIndex--;
                }

		if(tempKeyIndex<0)
		    break;

		
		    double[] tmpKeyArr = tmpKeyHt.get(tempKeyIndex);
		    if(null != tmpKeyArr)
			arr[steepCount] = tmpKeyArr;
		    else
		    {
			    arr[steepCount] = readFullSpectrum(tempKeyIndex, null, null, samIsoArr, mapModel, null, resultHt);
			
			tmpKeyHt.put(tempKeyIndex, arr[steepCount]);
		    }

		for(int k=0;k<3;k++)
		{
		    for(int l=0;l<resultSize;l++)
			prevArr[k][l] = arr[k][l];
		}

		

		    for(int j=leftIndex+1;j<=leftIndex+steepCount+1;j++)
		    {
			for(int k=resultIntensityStart;k<result[j].length;k++)                    
			{                        
			    area1 += result[j][k];
			}
		    }

		    for(int j=0;j<=steepCount;j++)
		    {
			for(int k=resultIntensityStart;k<result[j].length;k++)
			{
			    area2 += arr[j][k];
			}
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
        
                if( (isGoingUp && isHighIntensity) || areaRatio<steepRatioThreshold ) 
                {
		    if(leftIndex<=0)
		    {
                        isGoingUp=false;
                        break;
		    }

                    result[leftIndex] = arr[0];
		    
			moveLeftKeyIndex--; // -= numIsoWindow;

		    if(0 == result[leftIndex][0])
			continue;
		    
		    leftIndex--;

                    for(int j=0;j<steepCount;j++)
                        arr[j] = arr[j+1];
                    
                    for(int j=resultIntensityStart;j<arr.length;j++)
                        totalIntensity += arr[0][j];
                    
                    continue;                    
                }
                else
                    isGoingUp=false;

                                
                if(areaRatio>steepRatioThreshold)
                {
		    for(int k=0;k<3;k++)
		    {   
			if(result[leftIndex+1][resultIntensityStart]<arr[k][resultIntensityStart] || arr[k][0]==0)
			    break;

			if(leftIndex<0)
			    break;

			result[leftIndex] = arr[k];

			    moveLeftKeyIndex--; // -= numIsoWindow;

			leftIndex--;
		    }
                
                    break;
                }

                result[leftIndex--] = arr[0];
		    moveLeftKeyIndex--; // -= numIsoWindow;

                for(int j=0;j<steepCount;j++)
                    arr[j] = arr[j+1];
            }
        }

	prevArr = new double[steepArea][resultSize]; //clean up the array
         
        isGoingUp = true;
	isHighIntensity = true;

        if(rightIndex+steepArea<result.length)
        {
            while(true)
            {            
                double area1=0;
                double area2=0;
                
                if(moveRightKeyIndex+steepArea>=maxScanIndex)                    
                    break;                

                int steepCount=0;

                int tempKeyIndex = moveRightKeyIndex;
                while(true)
                {
		    if(tempKeyIndex>=maxScanIndex || steepCount>= steepArea-1 || steepCount>=maxScanIndex)
			break;

		    if(prevArr[0][0]!=0)
		    {
			for(int l=0;l<2;l++)
			    for(int m=resultIntensityStart;m<resultSize;m++)
				arr[l][m] = prevArr[l+1][m];
		    }
		    else
		    {
			    double[] tmpKeyArr = tmpKeyHt.get(tempKeyIndex);
			    if(null != tmpKeyArr)
				arr[steepCount] = tmpKeyArr;
			    else
			    {
				    arr[steepCount] = readFullSpectrum(tempKeyIndex, null, null, samIsoArr, mapModel, null, resultHt);			    

				tmpKeyHt.put(tempKeyIndex, arr[steepCount]);
			    }
		    }

                    steepCount++;

			tempKeyIndex++; // = tempKeyIndex + 1*numIsoWindow;
                }

		if(tempKeyIndex>=maxScanIndex)
		    break;

		
		    double[] tmpKeyArr = tmpKeyHt.get(tempKeyIndex);
		    if(null != tmpKeyArr)
			arr[steepCount] = tmpKeyArr;
		    else
		    {
			    arr[steepCount] = readFullSpectrum(tempKeyIndex, null, null, samIsoArr, mapModel, null, resultHt);

			tmpKeyHt.put(tempKeyIndex, arr[steepCount]);
		    }

		if(null == arr[steepCount])
		    break;

		for(int k=0;k<3;k++)
		{
		    for(int l=resultIntensityStart;l<resultSize;l++)
			prevArr[k][l] = arr[k][l];
		}


		
		    for(int j=rightIndex-1;j>rightIndex-steepCount-2;j--)
		    {
			for(int k=resultIntensityStart;k<result[j].length;k++)
			{                        
			    area1 += result[j][k];
			}
		    }

		    for(int j=0;j<=steepCount;j++)
		    {
			for(int k=resultIntensityStart;k<result[j].length;k++)
			{                        
			    area2 += arr[j][k];
			}
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
			moveRightKeyIndex++; 

		    if(0 == result[rightIndex][0])
			continue;
		    
		    rightIndex++;
                    
                    for(int j=0;j<steepCount;j++)
                        arr[j] = arr[j+1];

                    //rightTotalIntensity += arr[0][1];
                    
                    
                    for(int j=resultIntensityStart;j<resultSize;j++)
                    {                        
                        totalIntensity += arr[0][j];
                    }
                    
                    //totalIntensity += arr[0][1];
                    //totalIntensity += arr[0][2];
                    
                    continue;                
                }
                else
                    isGoingUp=false;
                   
                if(area2/area1>steepRatioThreshold)
                {                    
                    for(int k=0;k<3;k++)
                    {          
                        if(result[rightIndex-1][resultIntensityStart]<arr[k][resultIntensityStart])
                            break;

			if(result.length<=rightIndex)
			    break;
                       
                        result[rightIndex] = arr[k];
                            

			    moveRightKeyIndex++;
                        
                        rightIndex++;
                    }
                     break;
                }
                    
                result[rightIndex++] = arr[0];
                                
		    moveRightKeyIndex++; 

                for(int j=0;j<steepCount;j++)
                    arr[j] = arr[j+1];
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
      
	int peakStart = leftIndex+1;
	int peakEnd = rightIndex-1;

        //NonLabelMappingModel conf.getMapModel();
//        int dtaStart = range.getMin();
//        int dtaEnd = range.getMax();

	//for(int i=0;i<result.length;i++)
	 //   for(int j=0;j<result[i].length;j++)

        int specSpace = 1; //(int)(result[peakStart+1][0]-result[peakStart][0]);

//        int diff = (int)(result[peakStart][0]-dtaStart)/specSpace;
        
        int leftMargin = margin;

        int rightMargin = margin;
//        if(diff>0)
//            leftMargin += diff;            
//     
//        diff = (int)(dtaEnd-result[rightIndex-1][0])/specSpace;
//        if(diff>0)
//            rightMargin += diff;
        
        
	for(int i=0;i<leftMargin;i++)
        {   
            if(leftIndex<0 || moveLeftKeyIndex<=0)
                break;
  
            
		    result[leftIndex] = readFullSpectrum(moveLeftKeyIndex, null, null, samIsoArr, mapModel, null, resultHt);
                

            
	  
            
	    if(0 == result[leftIndex][0])
	    {
		
		    moveLeftKeyIndex--; 

		continue;
	    }

            leftIndex--;

	    
		moveLeftKeyIndex--;            
        }
        

        for(int i=0;i<rightMargin;i++)
        {
//        System.out.println(moveRightKeyIndex);
            if(rightIndex>=result.length || moveRightKeyIndex>=maxScanIndex || moveRightKeyIndex<0)
                break;

                    result[rightIndex] = readFullSpectrum(moveRightKeyIndex, null, null, samIsoArr, mapModel, null, resultHt);

            if(0 == result[rightIndex][0])
            {
                    moveRightKeyIndex++;

                continue;
            }

            rightIndex++;

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
        
         return new PeakRange((int)result[peakStart][0],(int)result[peakEnd][0]);
                    
    }
    
    /**
     * Pass the index of the PeakScanNumber to generate the start and end range
     * @param iFile
     * @param samIsoArr
     * @param keyIndex
     * @return PeakRange object with start and end scanNUmber 
     * @throws IOException
     * @throws CensusIndexOutOfBoundException
     * @throws Exception 
     */
  //  public static String peakFindingHarshil(SpectrumModel sModel, double[] samIsoArr, SpecRange range, int keyIndex, IndexedFile iFile)
      public static PeakRange peakFindingHarshil( IndexedFile iFile, double[] samIsoArr, int keyIndex)
        throws IOException, CensusIndexOutOfBoundException, Exception
    {
               
//	TIntLongHashMap index = sModel.getIndex();
//        Object file = sModel.getGenericIndexFile(conf);
        TIntLongHashMap index = iFile.getMsIndex();
//        Object file = sModel.getGenericIndexFile(conf);
        Object file = iFile.getFile();
        conf = Configuration.getInstance();
        steepRatioThreshold = conf.getSteepRatioThreshold();

	int numIsoWindow = conf.getNumOfIsolationWindow();
        int maxWindow=conf.getMaxWindow();
        int margin = conf.getMargin();

	double[][] result = null;
        //conf.getQuantLevel()

//	if(conf.isDataIndependent())
//	    result = new double[maxWindow*2+1+margin*2][4*sModel.getBioSample().length+3]; //scan #, sample intensity, ref intensity
//	else
	    result = new double[maxWindow*2+1+margin*2][3];

        int leftIndex=maxWindow+margin;
        int rightIndex=maxWindow+margin+1;
        
        
        //double leftTotalIntensity=0;
        //double rightTotalIntensity=0;
        double totalIntensity=0;

        int steepArea = conf.getSteepArea(); 
	
        int moveLeftKeyIndex = keyIndex;
               
        int moveRightKeyIndex = keyIndex+1*numIsoWindow;
        
       // int[] keys = sModel.getKeys();
         int[] keys = iFile.getKeys();
        int initWin=2;
       
        for(int i=0;i<initWin;i++) ///Index=maxWindow;keys[leftPeakIndex]>=leftMinValue; )
        {
	    if(moveLeftKeyIndex<=0 || leftIndex<=0)
	    {
		moveLeftKeyIndex += 1*numIsoWindow;
		leftIndex++;
		
		break;
	    }	
            
//	    if(conf.isDataIndependent()) //data independent
//	    {
//		result[leftIndex] = readSpectrum(keys, moveLeftKeyIndex, index, sModel.getDiff(), file, sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());
//
//	    }	
//	    else //data dependent
//	    {
//		if(sModel.isHighRes())
		    result[leftIndex] = readFullSpectrum(keys, moveLeftKeyIndex, index, file, samIsoArr);
//		else
//		    result[leftIndex] = readFullSpectrum(keys, moveLeftKeyIndex, index, file, sModel.getSamStartMass(), sModel.getSamEndMass(), sModel.getRefStartMass(), sModel.getRefEndMass());
//
//	    }
            totalIntensity += result[leftIndex][1];
            totalIntensity += result[leftIndex][2];
            
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

//	    if(conf.isDataIndependent()) //data independent
//	    {
//		result[rightIndex] = readSpectrum(keys, moveRightKeyIndex, index, sModel.getDiff(), file, sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());
//	    }
//	    else //data dependent
//	    {
//		if(sModel.isHighRes())
		    result[rightIndex] = readFullSpectrum(keys, moveRightKeyIndex, index, file, samIsoArr);
//		else
//		    result[rightIndex] = readFullSpectrum(keys, moveRightKeyIndex, index, file, sModel.getSamStartMass(), sModel.getSamEndMass(), sModel.getRefStartMass(), sModel.getRefEndMass());
//	    }

            totalIntensity += result[rightIndex][1];
            totalIntensity += result[rightIndex][2];

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
//			if(conf.isDataIndependent()) //data independent
//			    arr[steepCount] = readSpectrum(keys, tempKeyIndex, index, sModel.getDiff(), file, sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());
//			else
//			{
//			    if(sModel.isHighRes())
				arr[steepCount] = readFullSpectrum(keys, tempKeyIndex, index, file, samIsoArr);
//			    else
//				arr[steepCount] = readFullSpectrum(keys, tempKeyIndex, index, file, sModel.getSamStartMass(), sModel.getSamEndMass(), sModel.getRefStartMass(), sModel.getRefEndMass());
//			}            
		    }
		    
                    steepCount++;
                    //tempKeyIndex--;
		    tempKeyIndex = tempKeyIndex - 1*numIsoWindow;
                }

		if(tempKeyIndex<0)
		    break;

//		if(conf.isDataIndependent()) //data independent
//		    arr[steepCount] = readSpectrum(keys, tempKeyIndex, index, sModel.getDiff(), file, sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());
//		else
//		{
//		    if(sModel.isHighRes())
			arr[steepCount] = readFullSpectrum(keys, tempKeyIndex, index, file, samIsoArr);
//		    else
//			arr[steepCount] = readFullSpectrum(keys, tempKeyIndex, index, file, sModel.getSamStartMass(), sModel.getSamEndMass(), sModel.getRefStartMass(), sModel.getRefEndMass());
//		}

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
//			if(conf.isDataIndependent()) //data independent
//			    arr[steepCount] = readSpectrum(keys, tempKeyIndex, index, sModel.getDiff(), file, sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());
//			else if(sModel.isHighRes())
			    arr[steepCount] = readFullSpectrum(keys, tempKeyIndex, index, file, samIsoArr);
//			else
//			    arr[steepCount] = readFullSpectrum(keys, tempKeyIndex, index, file, sModel.getSamStartMass(), sModel.getSamEndMass(), sModel.getRefStartMass(), sModel.getRefEndMass());
		    }

                    steepCount++;
		    tempKeyIndex = tempKeyIndex + 1*numIsoWindow;
                    //tempKeyIndex++;
                }

		if(tempKeyIndex>=keys.length)
		    break;
//
//		if(conf.isDataIndependent()) //data independent
//		    arr[steepCount] = readSpectrum(keys, tempKeyIndex, index, sModel.getDiff(), file, sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());
//		else if(sModel.isHighRes())
                    arr[steepCount] = readFullSpectrum(keys, tempKeyIndex, index, file, samIsoArr);
//		else
//		    arr[steepCount] = readFullSpectrum(keys, tempKeyIndex, index, file, sModel.getSamStartMass(), sModel.getSamEndMass(), sModel.getRefStartMass(), sModel.getRefEndMass());

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
	int peakStart = leftIndex+1;
	int peakEnd = rightIndex-1;
//
//        int dtaStart = (range==null?0:range.getMin());
//        int dtaEnd = (range==null?0:range.getMax());
        int specSpace = (int)(result[peakStart+1][0]-result[peakStart][0]);

//System.out.println( result[peakStart+1][0] + " " + result[peakStart][0] );
//        int diff = (int)(result[peakStart][0]-dtaStart)/specSpace;
        int leftMargin = margin;

        int rightMargin = margin;
//        if(diff>0)
//            leftMargin += diff;            
//     
//        diff = (int)(dtaEnd-result[rightIndex-1][0])/specSpace;
//
//        if(diff>0)
//            rightMargin += diff;
        
	for(int i=0;i<leftMargin;i++)
        {   
            if(leftIndex<0 || moveLeftKeyIndex<0)
                break;

//	    if(conf.isDataIndependent()) //data independent
//		result[leftIndex] = readSpectrum(keys, moveLeftKeyIndex, index, sModel.getDiff(), file, sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());
//	    else if(sModel.isHighRes())
		result[leftIndex] = readFullSpectrum(keys, moveLeftKeyIndex, index, file, samIsoArr);
//	    else
//		result[leftIndex] = readFullSpectrum(keys, moveLeftKeyIndex, index, file, sModel.getSamStartMass(), sModel.getSamEndMass(), sModel.getRefStartMass(), sModel.getRefEndMass());

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
            
//	    if(conf.isDataIndependent()) //data independent
//		result[rightIndex] = readSpectrum(keys, moveRightKeyIndex, index, sModel.getDiff(), file, sModel.getBioSample(), sModel.getBioRef(), sModel.getYioSample(), sModel.getYioRef());
//	    else if(sModel.isHighRes())
		result[rightIndex] = readFullSpectrum(keys, moveRightKeyIndex, index, file, samIsoArr); 
//	    else
//		result[rightIndex] = readFullSpectrum(keys, moveRightKeyIndex, index, file, sModel.getSamStartMass(), sModel.getSamEndMass(), sModel.getRefStartMass(), sModel.getRefEndMass());
//            
	    if(0 == result[rightIndex][0])
	    {
		moveRightKeyIndex += numIsoWindow;
		continue;
	    }

            rightIndex++;
	    moveRightKeyIndex += numIsoWindow;
        }
//	if(conf.isDataIndependent())
//	    return buildDIResult(peakStart, peakEnd, leftIndex, rightIndex, result, sModel.getBioSample().length);
//	else
//	    return buildDDResult(peakStart, peakEnd, leftIndex, rightIndex, result);
            return new PeakRange((int)result[peakStart][0],(int)result[peakEnd][0]);
    }
      
    public static double[] intensitySum(double[] massArr, double[] intArr, double[] isoArr)
    {

        double[] resultArr = new double[3];
        
        double sumIntensity=0;
        int foundIsoNum=0; 
        int totalPeakFound=0;
        double toleranceSum=0;
        double massTolerance = Configuration.getInstance().getMassTolerance();
	for(int i=0;i<isoArr.length;i++)
	{
            double tempTolerance = isoArr[i]/1000*massTolerance;
            //double tempTolerance = isoArr[i]/1000*massTolerance*0.001;
//            double tempTolerance = isoArr[i]*massTolerance;
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
            
            if(isFound)
                foundIsoNum++;
	}

        resultArr[0] = sumIntensity;
        resultArr[1] = foundIsoNum;
        resultArr[2] = (totalPeakFound>0)?(toleranceSum/totalPeakFound):(-1);
        
	return resultArr;
    }
}
