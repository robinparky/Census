package edu.scripps.pms.census.util;

//import com.sun.jndi.ldap.PersistentSearchControl;
import java.io.*;

import edu.scripps.pms.averagine.AssignMass;
import edu.scripps.pms.census.labelFree.SpectrumModel;
import edu.scripps.pms.census.model.Spectrum;
import edu.scripps.pms.util.spectrum.PeakList;
import edu.scripps.pms.util.sqlite.spectra.SpectraDB;
import gnu.trove.*;

import java.io.RandomAccessFile;
import edu.scripps.pms.census.hash.IndexedFile;
import edu.scripps.pms.census.CensusConstants;

import static edu.scripps.pms.census.ChroGenerator.createIndexedFiles;

import edu.scripps.pms.census.conf.*;
import edu.scripps.pms.census.exception.CensusIndexOutOfBoundException;

import java.sql.SQLException;
import java.util.*;
import java.text.*;

import org.apache.commons.math3.stat.regression.SimpleRegression;

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
public class CalcUtilGeneric
{
    private static final char SPACE = ' ';
    private static int ION_START_INDEX = 3;
    private static final char CARRIAGE_RETURN = '\n';
    private static final char WINDOW_CR = '\r';
    private static final char DOT = '.';
    public static final String ILINERETSTRING = "I\tRetTime";
    public static final char TAB = '\t';
    private static DecimalFormat formatter = new DecimalFormat("0.0000");

    /*
    public static void getSpectrumArr(IndexedFile iFile, int scanNum, double tolerance) throws IOException, CensusIndexOutOfBoundException
    {
        //double[][] arr = getSpectrumArray(iFile, scanNum);

        //for(int i=0;i<arr.length;i++) {

        //}

    }*/


    public static boolean isIsolateContaminanted(IndexedFile iFile, int scanNum, double tolerance, double precursor, int chargeState,String sequence) throws IOException, CensusIndexOutOfBoundException
    {
	RandomAccessFile file = iFile.getFile();
        TIntLongHashMap index = iFile.getMsIndex();
	int[] keys = iFile.getKeys();
	int curIndex = Arrays.binarySearch(keys, scanNum);
	if(curIndex<0) //Cannot find index
	    curIndex=-(++curIndex); //Math.abs(++keyIndex);

	if(curIndex>=keys.length)
	    curIndex--;

	curIndex--;
	if(keys.length<=curIndex)
	    throw new CensusIndexOutOfBoundException();

        long startPos = index.get(keys[curIndex]);
        long endPos;

        double startRange = precursor-tolerance;
        double endRange = precursor+tolerance;

        file.seek(startPos);

        if( (curIndex+1)>=keys.length )
            endPos = file.length();
        else
            endPos = index.get(keys[curIndex+1]);

        int byteSize = (int)(endPos-startPos);

        byte[] bytes = new byte[byteSize];
        file.readFully(bytes);

        char ch;
        int pos=0;

        ch = (char)bytes[pos];

        //Remove Z, S, I, D lines
        while( (ch=(char)bytes[pos]) == 'Z' || ch == 'S' || ch == 'I' || ch == 'D')
        {
            while( ch != CARRIAGE_RETURN )
            {
                pos++;
                ch = (char)bytes[pos];
            }

            pos++;
        }

        StringBuilder mass = new StringBuilder(10);
        StringBuilder intensity = new StringBuilder(10);
	intensity.append('0');

        boolean isMass=true;
        boolean isInt=true;

        int arrSize=0;
        for(int j=pos;j<byteSize;j++)
        {
            if( CARRIAGE_RETURN == (char)bytes[j] )
                arrSize++;
        }

	double[] massArr = new double[arrSize];
	double[] intArr = new double[arrSize];
	double[][] resultArr = new double[2][arrSize];

        List<Double> isoMassList = new ArrayList<Double>();
        List<Double> isoIntList = new ArrayList<Double>();

	resultArr[0] = massArr;
	resultArr[1] = intArr;

        int massIndex=0;

	int spaceCount=0;

        double maxIsoPeak = 0.0;
        double maxIsoMass = 0.0;

        for(int i=pos;i<byteSize;i++)
        {
            ch = (char)bytes[i];
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
                    massArr[massIndex] = Double.parseDouble(mass.toString());

                    if(massArr[massIndex]>=startRange && massArr[massIndex]<=endRange) {
                        isoMassList.add(massArr[massIndex]);
                        isoIntList.add(intArr[massIndex]);
                        if(maxIsoPeak<intArr[massIndex]) {
                            maxIsoPeak = intArr[massIndex];
			    maxIsoMass = massArr[massIndex];
			}
                    }

                    massIndex++;

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
        }


        double diff=10000000;

        List<Peak> peakList=new ArrayList();
        double intThreshold = maxIsoPeak*5.0/100.00;

	for(int i=0;i<isoIntList.size();i++) {

            double d = isoIntList.get(i);
            //CHeck if the peack is within the 5% of the Max Peak...

            if(d<=intThreshold) {
                continue;
            }


            peakList.add(new Peak(isoMassList.get(i), isoIntList.get(i)));
        }

if(peakList.size()<=0) return false;

	List matchCounts = new ArrayList();
	double dist = 1.003354826/chargeState;

	int maxPeakIndex=0;
	double maxInt=0;
	int maxIntIndex=0;
	int count=0;
        for(int i=0;i<peakList.size();i++)
        {
           // double d = peakList.get(i).getIntensity();

           int matchCounter = 0;
           double totalIntensity = 0;
           // int missedCounter = 0;

//	    double massDiff = Math.abs(maxIsoMass - peakList.get(i).getMass());



//	    double minDiff = 1000000;

	    Peak p = peakList.get(i);
	    for(int j=-5;j<=5;j++) {
//		    double tmpDiff = Math.abs(massDiff-dist*j);
//		    if(minDiff>tmpDiff)
//			    minDiff = tmpDiff;
                    //int flag = 0;
                    //int z=0;
		double shiftM =p.getMass() + dist*j;
		double shiftMStart = shiftM-0.02;
		double shiftMEnd = shiftM+0.02;



		for(int k= 0; k<peakList.size();k++)
		{
			Peak eachP = peakList.get(k);

			if(eachP.getMass()>=shiftMStart && eachP.getMass()<=shiftMEnd) {
				matchCounter++;
				p.addIndex(k);
				totalIntensity+=eachP.getIntensity();
			}
		}
	    }

	p.setPeakMatch(matchCounter);
	p.setTotalIsoIntensity(totalIntensity);

	if(count<matchCounter) {
		maxPeakIndex = i;
		count = matchCounter;
	}

	if(maxInt<totalIntensity) {
		maxIntIndex = i;
		maxInt = totalIntensity;
	}


       //     matchCounts.add(matchCounter);
	//    if(minDiff<0.02) continue;  //roughly 20 ppm.  need to be corrected for ppm, instead of m/z
/*
            double tmp = maxIsoPeak - d;
            if(tmp==0) continue;
            if(tmp<diff)
                diff = tmp;
*/
        }



	//Peak matchPeak = peakList.get(maxPeakIndex);
	Peak matchPeak = peakList.get(maxIntIndex);


        float trueIntensitySum=0.0f;
	float falseIntensitySum=0.0f;

	Set matchSet = matchPeak.getIndexSet();
        for(int i=0;i<peakList.size();i++)
        {
            if(matchSet.contains(i)) {
                trueIntensitySum+=peakList.get(i).getIntensity();
	    }
            else {
                falseIntensitySum+=peakList.get(i).getIntensity();
	    }
        }


      //  System.out.print("\n-----" +sequence +"\t" + scanNum + "\t" + chargeState + "\t" + matchPeak.getIndexSet().size() + "\t" + (peakList.size() - matchPeak.getIndexSet().size()));
     //  System.out.println("\t" + trueIntensitySum +"\t" + falseIntensitySum);


        //diff is smallest gap comparing to between base peak
        double ratio = diff/maxIsoPeak;
        if(ratio<=0.8)
            return true;
        else
            return false;




	//return resultArr;
    }

    public static double[][] getSpectrumArr(IndexedFile iFile, int scanNum) throws IOException, CensusIndexOutOfBoundException
    {
	RandomAccessFile file = iFile.getFile();
        TIntLongHashMap index = iFile.getMsIndex();
	int[] keys = iFile.getKeys();
	int curIndex = Arrays.binarySearch(keys, scanNum);
	if(curIndex<0) //Cannot find index
	    curIndex=-(++curIndex); //Math.abs(++keyIndex);

	if(curIndex>=keys.length)
	    curIndex--;

	if(keys.length<=curIndex)
	    throw new CensusIndexOutOfBoundException();

        long startPos = index.get(keys[curIndex]);
        long endPos;

        file.seek(startPos);

        if( (curIndex+1)>=keys.length )
            endPos = file.length();
        else
            endPos = index.get(keys[curIndex+1]);

        int byteSize = (int)(endPos-startPos);

        byte[] bytes = new byte[byteSize];
        file.readFully(bytes);

        char ch;
        int pos=0;

        ch = (char)bytes[pos];

        //Remove Z, S, I, D lines
        while( (ch=(char)bytes[pos]) == 'Z' || ch == 'S' || ch == 'I' || ch == 'D')
        {
            while( ch != CARRIAGE_RETURN )
            {
                pos++;
                ch = (char)bytes[pos];
            }

            pos++;
        }

        StringBuilder mass = new StringBuilder(10);
        StringBuilder intensity = new StringBuilder(10);
	intensity.append('0');

        boolean isMass=true;
        boolean isInt=true;

        int arrSize=0;
        for(int j=pos;j<byteSize;j++)
        {
            if( CARRIAGE_RETURN == (char)bytes[j] )
                arrSize++;
        }

	double[] massArr = new double[arrSize];
	double[] intArr = new double[arrSize];
	double[][] resultArr = new double[2][arrSize];
	resultArr[0] = massArr;
	resultArr[1] = intArr;

        int massIndex=0;

	int spaceCount=0;

        for(int i=pos;i<byteSize;i++)
        {
            ch = (char)bytes[i];
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
        }


	return resultArr;
    }

        public static String getSpectrumString(IndexedFile iFile, int scanNum) throws IOException, CensusIndexOutOfBoundException
    {
	RandomAccessFile file = iFile.getFile();
        TIntLongHashMap index = iFile.getMsIndex();
	int[] keys = iFile.getKeys();
	int curIndex = Arrays.binarySearch(keys, scanNum);
	if(curIndex<0) //Cannot find index
	    curIndex=-(++curIndex); //Math.abs(++keyIndex);

	if(curIndex>=keys.length)
	    curIndex--;

	if(keys.length<=curIndex)
	    throw new CensusIndexOutOfBoundException();

        long startPos = index.get(keys[curIndex]);
        long endPos;

        file.seek(startPos);

        if( (curIndex+1)>=keys.length )
            endPos = file.length();
        else
            endPos = index.get(keys[curIndex+1]);

        int byteSize = (int)(endPos-startPos);

        byte[] bytes = new byte[byteSize];
        file.readFully(bytes);

        return new String(bytes);
            }


    public static long getBackGroundNoise (
           // double[] isoArr,
            long startPos,
            int diffPos,
            IndexedFile iFile
    ) throws IOException, CensusIndexOutOfBoundException, Exception
    {

        RandomAccessFile file = iFile.getFile();
        file.seek(startPos);

        int byteSize = diffPos;
        byte[] bytes = new byte[byteSize];
        file.readFully(bytes);

        String str = new String(bytes);
        String[] lines = str.split(System.getProperty("line.separator"));

        int count=0;

        for(String each:lines) {

            if( Character.isDigit(each.charAt(0)) )
                break;

            count++;
        }

        //double[] massArr = new double[lines.length-count];
        int[] intArr = new int[lines.length-count];
        // double[] rtArr = new double[lines.length-count];

        int massIndex=0;

        for(int i=count;i<lines.length;i++) {

          //  System.out.println(lines[i]);

            String[] arr = lines[i].split(" ");

      //  try {    //massArr[massIndex] = Double.parseDouble(arr[0]);
            intArr[massIndex] = (int) Double.parseDouble(arr[1]);
            massIndex++;
      //  } catch(Exception e) {
      //      System.out.println("aaa");
      //      throw new Exception("aa");
      //  }
        }

        System.out.println();
        System.out.print(">>>");
        for(double d: intArr)
        System.out.println(d+" ");
        System.out.println();

        Arrays.sort(intArr);

        int size = (int)(intArr.length*0.05);

        if(size == 0 ) size=1;
            //System.out.println("==========" + byteSize +  " " + diffPos + " " + intArr.length);


        double sum = 0;

        if(intArr.length<=0) return 1000;  //fixed background noise when there are no peaks

        for(int i=0;i<size;i++) {
            sum+=intArr[i];
        }

        return (long)sum/size;
    }


    public static long getBackGroundNoise (int scanNumber, SpectraDB spectraDB)
            throws SQLException {
        SpectraDB.Spectrum spectrum =  spectraDB.getSpectrumFromDB(scanNumber);
       // double [] tempIntArr = spectrum.getIntensityList().toNativeArray();
        int [] intArr = new int[spectrum.getIntensityList().size()];
        for(int i=0; i<spectrum.getIntensityList().size(); i++)
        {
            intArr[i] = (int)spectrum.getIntensityList().get(i);
        }


        Arrays.sort(intArr);

        int size = (int)(intArr.length*0.05);

        if(size == 0 ) size=1;
        //System.out.println("==========" + byteSize +  " " + diffPos + " " + intArr.length);

        //System.out.println("<<<>>>< "+intArr.length);
        double sum = 0;

        if(intArr.length<=0) return 1000;  //fixed background noise when there are no peaks

        for(int i=0;i<size;i++) {
            sum+=intArr[i];
        }

        return (long)sum/size;
    }


    public static double [][] readLabelfreeFullSpectrum(
            long startPos,
            int diffPos,
            IndexedFile iFile,
            int chargeState,
            String [] scanRetStrings,
            Configuration conf
    ) throws IOException {
        RandomAccessFile file = iFile.getFile();

        file.seek(startPos);

        int byteSize = diffPos;
        byte[] bytes = new byte[byteSize];
        file.readFully(bytes);


        String str = new String(bytes);
        String[] lines = str.split(System.getProperty("line.separator"));

        int count=0;
        double retTime=0;
        int scan=0;
        int oldii=0;
        int ii=0;

        for(String each:lines) {

            if( Character.isDigit(each.charAt(0)) )
                break;
            else if(each.startsWith("S\t")) {
                String tmpAr = each.split("\t")[1];
                scanRetStrings[0] = tmpAr;
                //scan = Integer.parseInt(tmpAr);
            }   else if(each.startsWith("I\tRetTime")) {
                String tmpAr = each.split("\t")[2];
                scanRetStrings[1] = tmpAr;
               // retTime = Double.parseDouble(tmpAr);
            }


            count++;
        }
        double [][] result = new double[2][];
        double[] massArr = new double[lines.length-count];
        double[] intArr = new double[lines.length-count];
        // double[] rtArr = new double[lines.length-count];
        result[0] = massArr;
        result[1] = intArr;

        int massIndex=0;

        for(int i=count;i<lines.length;i++) {
            String[] arr = lines[i].split(" ");
            if(conf.isLabelfreeCheckChargeState() && arr.length>2) {
                int cs = Integer.parseInt(arr[2]);
                if(chargeState==cs) {
                    massArr[massIndex] = Double.parseDouble(arr[0]);
                    intArr[massIndex] = Double.parseDouble(arr[1]);
                }
            } else {
                massArr[massIndex] = Double.parseDouble(arr[0]);
                intArr[massIndex] = Double.parseDouble(arr[1]);
            }

            massIndex++;
        }

        return result;
    }

    public static class SpectraArrays
    {
        public final double [] massArray;
        public final double [] intensityArray;
        public final int [] csArray;

        public SpectraArrays(double[] massArray, double[] intensityArray, int[] csArray) {
            this.massArray = massArray;
            this.intensityArray = intensityArray;
            this.csArray = csArray;
        }
    }

    public static SpectraArrays readLabelfreeFullSpectrum(
            long startPos,
            int diffPos,
            IndexedFile iFile,
            String [] scanRetStrings,
            Configuration conf
    ) throws IOException {
        RandomAccessFile file = iFile.getFile();

        file.seek(startPos);

        int byteSize = diffPos;
        byte[] bytes = new byte[byteSize];
        file.readFully(bytes);


        String str = new String(bytes);
        String[] lines = str.split(System.getProperty("line.separator"));

        int count=0;
        double retTime=0;
        int scan=0;
        int oldii=0;
        int ii=0;

        for(String each:lines) {

            if( Character.isDigit(each.charAt(0)) )
                break;
            else if(each.startsWith("S\t")) {
                String tmpAr = each.split("\t")[1];
                scanRetStrings[0] = tmpAr;
                //scan = Integer.parseInt(tmpAr);
            }   else if(each.startsWith("I\tRetTime")) {
                String tmpAr = each.split("\t")[2];
                scanRetStrings[1] = tmpAr;
                // retTime = Double.parseDouble(tmpAr);
            }


            count++;
        }
        double [][] result = new double[2][];
        double[] massArr = new double[lines.length-count];
        double[] intArr = new double[lines.length-count];
        int[] csArr = new int[lines.length-count];
        // double[] rtArr = new double[lines.length-count];
        result[0] = massArr;
        result[1] = intArr;

        int massIndex=0;

        for(int i=count;i<lines.length;i++) {
            String[] arr = lines[i].split(" ");
            if(conf.isLabelfreeCheckChargeState() && arr.length>2) {
                int cs = Integer.parseInt(arr[2]);
                massArr[massIndex] = Double.parseDouble(arr[0]);
                intArr[massIndex] = Double.parseDouble(arr[1]);
                csArr[massIndex]=Integer.parseInt(arr[2]);
            } else {
                massArr[massIndex] = Double.parseDouble(arr[0]);
                intArr[massIndex] = Double.parseDouble(arr[1]);
            }

            massIndex++;
        }

        return new SpectraArrays(massArr,intArr,csArr);
    }

    public static SpectrumModel   readLabelfreeFullSpectrum(

            PeakList peakList

    ) throws IOException, CensusIndexOutOfBoundException, Exception
    {

        TDoubleArrayList massList = new TDoubleArrayList();
        TDoubleArrayList intensityList = new TDoubleArrayList();
        TIntArrayList csList = new TIntArrayList();
        double retTime = peakList.getRetentionTime();
        int scanNum = peakList.getHiscan();
        for(ListIterator<edu.scripps.pms.util.spectrum.Peak> peakListIterator = peakList.getPeaks();
            peakListIterator.hasNext(); )
        {
            edu.scripps.pms.util.spectrum.Peak peak = peakListIterator.next();
            massList.add(peak.getM2z());
            intensityList.add(peak.getIntensity());
            csList.add(peak.getChargeState());
        }
        double [] massArr = massList.toNativeArray();
        double [] intArr = intensityList.toNativeArray();
        int [] csArr = csList.toNativeArray();

        edu.scripps.pms.census.labelFree.SpectrumModel spec = new edu.scripps.pms.census.labelFree.SpectrumModel();
        spec.setRetentionTime(retTime);
        spec.setScanNumber(scanNum);
        spec.setMass(massArr);
        spec.setIntensity(intArr);
        spec.setCsArray(csArr);

        return spec;
    }


    public static SpectrumModel   readLabelfreeFullSpectrum(
            long startPos,
            int diffPos,
            IndexedFile iFile

    ) throws IOException, CensusIndexOutOfBoundException, Exception
    {

        RandomAccessFile file = iFile.getFile();

        file.seek(startPos);

        int byteSize = diffPos;
        byte[] bytes = new byte[byteSize];
        file.readFully(bytes);
/*
        char ch;
        int pos=0;

        ch = (char)bytes[pos];

        //Remove Z, S, I, D lines
        while( (ch=(char)bytes[pos]) == 'Z' || ch == 'S' || ch == 'I' || ch == 'D')
        {
            while( ch != CARRIAGE_RETURN )
            {
                pos++;
                ch = (char)bytes[pos];
            }

            pos++;
        }*/


        String str = new String(bytes);
        String[] lines = str.split(System.getProperty("line.separator"));

        int count=0;
        double retTime=0;
        int scan=0;




        for(String each:lines) {

            if( Character.isDigit(each.charAt(0)) )
                break;
            else if(each.startsWith("S\t")) {
                String tmpAr = each.split("\t")[1];
                scan = Integer.parseInt(tmpAr);
            }   else if(each.startsWith("I\tRetTime")) {
                String tmpAr = each.split("\t")[2];
                retTime = Double.parseDouble(tmpAr);
            }


            count++;
        }

        double[] massArr = new double[lines.length-count];
        double[] intArr = new double[lines.length-count];
        int[] csArr = new int[lines.length-count];
        // double[] rtArr = new double[lines.length-count];

        int massIndex=0;

        for(int i=count;i<lines.length;i++) {
            String[] arr = lines[i].split(" ");


            // System.out.println(lines[i]);

            //check charge state



                massArr[massIndex] = Double.parseDouble(arr[0]);
                intArr[massIndex] = Double.parseDouble(arr[1]);
               if(arr.length==3) csArr[massIndex] = Integer.parseInt(arr[2]);
            //   System.out.println(massArr[massIndex] + "\t" + intArr[massIndex]);
            //    rtArr[massIndex] = Double.parseDouble(arr[2]);
            massIndex++;
        }

        // double[][] result = new double[2][2];
        // result[0] = massArr;
        //result[1] = intArr;

        edu.scripps.pms.census.labelFree.SpectrumModel spec = new edu.scripps.pms.census.labelFree.SpectrumModel();
        // spec.setMass(massArr);
        // spec.setIntensity(intArr);
        spec.setRetentionTime(retTime);
        spec.setScanNumber(scan);
        spec.setMass(massArr);
        spec.setIntensity(intArr);
        spec.setCsArray(csArr);
        file.close();

        return spec;
    }

    public static edu.scripps.pms.census.labelFree.SpectrumModel     readLabelfreeFullSpectrum(
            double[] isoArr,
            long startPos,
            int diffPos,
            double massTolerance,
            IndexedFile iFile,
            int chargeState,
            Configuration conf,
            double pepMass
            ) throws IOException, CensusIndexOutOfBoundException, Exception
    {

        RandomAccessFile file = iFile.getFile();

        file.seek(startPos);

        int byteSize = diffPos;
        byte[] bytes = new byte[byteSize];
        file.readFully(bytes);
/*
        char ch;
        int pos=0;

        ch = (char)bytes[pos];

        //Remove Z, S, I, D lines
        while( (ch=(char)bytes[pos]) == 'Z' || ch == 'S' || ch == 'I' || ch == 'D')
        {
            while( ch != CARRIAGE_RETURN )
            {
                pos++;
                ch = (char)bytes[pos];
            }

            pos++;
        }*/


        String str = new String(bytes);
        String[] lines = str.split(System.getProperty("line.separator"));

        int count=0;
        double retTime=0;
        int scan=0;




        for(String each:lines) {

            if( Character.isDigit(each.charAt(0)) )
                break;
            else if(each.startsWith("S\t")) {
                String tmpAr = each.split("\t")[1];
                scan = Integer.parseInt(tmpAr);
            }   else if(each.startsWith("I\tRetTime")) {
                String tmpAr = each.split("\t")[2];
                retTime = Double.parseDouble(tmpAr);
            }


            count++;
        }

	double[] massArr = new double[lines.length-count];
	double[] intArr = new double[lines.length-count];
       // double[] rtArr = new double[lines.length-count];

        int massIndex=0;

        for(int i=count;i<lines.length;i++) {
            String[] arr = lines[i].split(" ");


           // System.out.println(lines[i]);

            //check charge state
            if(conf.isLabelfreeCheckChargeState() && arr.length>2) {
                int cs = Integer.parseInt(arr[2]);
                if(chargeState==cs) {
                    massArr[massIndex] = Double.parseDouble(arr[0]);
                    intArr[massIndex] = Double.parseDouble(arr[1]);
                }
            } else {
                massArr[massIndex] = Double.parseDouble(arr[0]);
                intArr[massIndex] = Double.parseDouble(arr[1]);
            }
         //   System.out.println(massArr[massIndex] + "\t" + intArr[massIndex]);
        //    rtArr[massIndex] = Double.parseDouble(arr[2]);
            massIndex++;
        }

       // double[][] result = new double[2][2];
       // result[0] = massArr;
	//result[1] = intArr;

        edu.scripps.pms.census.labelFree.SpectrumModel spec = new edu.scripps.pms.census.labelFree.SpectrumModel();
       // spec.setMass(massArr);
       // spec.setIntensity(intArr);
        spec.setRetentionTime(retTime);
        spec.setScanNumber(scan);
        spec.setMass(massArr);
        spec.setIntensity(intArr);

	//result[0] = keys[curIndex];
	double[] tempArr = intensitySumWithIsotopeModeling(massArr, intArr, isoArr, massTolerance, pepMass);
        if(null != tempArr) {
            spec.setPrecursorPeakIntensity((long) tempArr[0]);
            //System.out.println(iFile.getFileName()+"\t"+scan+"\t"+(long) tempArr[0]);
            //System.exit(1);

        }
  /*      if(scan == 26391)
        {
            System.out.print(iFile.getFileName());
            System.out.println("Scan Number"+scan);
            for(int i=0; i<massArr.length;i++)
            {
                System.out.println(massArr[i]+"\t"+intArr[i]);
            }
        }*/
     //   System.out.println(tempArr[0]);

       // result[0] = tempArr[0];     //intensity
       // result[1] = isoArr.length;
       // result[2] = tempArr[1];      //Found iso num
       // result[3] = tempArr[2];	     //


	return spec;
    }


       //High resolution //lalebling data
    public static edu.scripps.pms.census.labelFree.SpectrumModel readLabelfreeFullSpectrumFast(
            double[] isoArr,
            long startPos,
            int diffPos,
            double massTolerance,
            IndexedFile iFile,
            BufferedRandomAccessFile file
            ) throws IOException, CensusIndexOutOfBoundException, Exception
    {

        //RandomAccessFile file = iFile.getFile();



        file.seek(startPos);

        int byteSize = diffPos;
        byte[] bytes = new byte[byteSize];
        file.readFully(bytes);
/*
        char ch;
        int pos=0;

        ch = (char)bytes[pos];

        //Remove Z, S, I, D lines
        while( (ch=(char)bytes[pos]) == 'Z' || ch == 'S' || ch == 'I' || ch == 'D')
        {
            while( ch != CARRIAGE_RETURN )
            {
                pos++;
                ch = (char)bytes[pos];
            }

            pos++;
        }*/


      //  System.out.println("");

        //String str = new String(bytes);
        //System.out.println(str);

        StringBuffer sb = new StringBuffer();
        int count=0;
        for(byte b : bytes) {
            sb.append((char)b);

            if(b == '\n') {
               // System.out.print("==" + sb.toString());
              //  System.out.println("===" + (char)bytes[count+1] + "===");

                if( Character.isDigit( (char)bytes[count+1] ) ) {
                    count++;
                    break;
                }
        //        sb.delete(0, sb.length());
            }

            count++;
        }


        String[] lines = sb.toString().split(System.getProperty("line.separator"));


        double retTime=0;
        int scan=0;

        for(String each:lines) {

            if( Character.isDigit(each.charAt(0)) )
                break;
            else if(each.startsWith("S\t")) {
                String tmpAr = each.split("\t")[1];
                scan = Integer.parseInt(tmpAr);
            }   else if(each.startsWith("I\tRetTime")) {
                String tmpAr = each.split("\t")[2];
                retTime = Double.parseDouble(tmpAr);
            }


        }


       // System.out.println( "===========" + (char)bytes[count] );

    //   System.out.println("===sb " + scan + " "  + retTime + " \n" + sb.toString());


        sb.delete(0, sb.length());


        StringBuffer massSb = new StringBuffer(16);
        StringBuffer intSb = new StringBuffer(16);
        StringBuffer retSb = new StringBuffer(16);

        int countSpace=0;
        boolean intSpace=false;
       // int scan=0;



        TDoubleArrayList massArr = new TDoubleArrayList();
        TIntArrayList intArr = new TIntArrayList();
        //TDoubleArrayList retArr = new TDoubleArrayList();
       // double[] massArr = new double[lines.length-count];
       // double[] intArr = new double[lines.length-count];
        //double[] rtArr = new double[lines.length-count];

        int massIndex=0;
        for(int i=count;i<bytes.length;i++) {

            if(bytes[i] == ' ') countSpace++;
            char ch = (char)bytes[i];
            System.out.print(ch);
            switch(countSpace) {
                case 0:
                    massSb.append(ch);
                    break;
                case 1:
                    if(ch == '.') intSpace = true;
                    if(!intSpace)
                        intSb.append(ch);
                    break;

               // case 2:
              //      retTimeSb.append((char)bytes[i]);
               //     break;
            }
             // sb.append((char)bytes[i]);

            if( bytes[i] == '\n') {
                countSpace=0;
                intSpace = false;

            //    System.out.println(massSb.toString() + "\t" + intSb.toString()); // + "\t" + retTimeSb.toString());

                //double dMass = (dotIndex>0)?mass/10000d:mass;
                double dMass = Double.parseDouble(massSb.toString());
               // System.out.println("corr==" + dMass + "\t" + mass + "==\t" + intSb.toString());
               // try {
                int intensity= Integer.parseInt(intSb.toString().trim());

                massArr.add(dMass);
                intArr.add(intensity);

          //      System.out.println(dMass + " " + intensity + " " + retTime);
                massSb.delete(0, massSb.length());
                intSb.delete(0, intSb.length());
               // retTimeSb.delete(0, retTimeSb.length());

              }
        }

      //  System.exit(0);
         /*
        String[] lines = str.split(System.getProperty("line.separator"));


       // double[][] result = new double[2][2];
       // result[0] = massArr;
	//result[1] = intArr;
        */
        edu.scripps.pms.census.labelFree.SpectrumModel spec = new edu.scripps.pms.census.labelFree.SpectrumModel();
       // spec.setMass(massArr);
       // spec.setIntensity(intArr);
        spec.setRetentionTime(retTime);
        spec.setScanNumber(scan);
        spec.setIntArrayList(intArr);
        spec.setMassArrayList(massArr);

    return spec;


        }
    public static String[] debugGetString(byte[] arr)
    {
        String str = new String(arr);
        String[] lines = str.split(System.getProperty("line.separator"));
        return lines;
    }


    public static edu.scripps.pms.census.labelFree.SpectrumModel labelFreeSpectrumReader(
            double[] isoArr,
            int scanNumber,
            double massTolerance,
           SpectraDB spectraDB,
            int chargeState,
            Configuration conf,
            double pepMass
    ) throws IOException, CensusIndexOutOfBoundException, Exception
    {

        SpectraDB.Spectrum spectrum = spectraDB.getSpectrumFromDB(scanNumber);
        double retTime = spectrum.retTime;
        double [] massArr = spectrum.getMzList().toNativeArray();
        double [] intArr = spectrum.getIntensityList().toNativeArray();

        edu.scripps.pms.census.labelFree.SpectrumModel spec = new edu.scripps.pms.census.labelFree.SpectrumModel();

        spec.setRetentionTime(retTime);
        spec.setScanNumber(scanNumber);

        double[] tempArr = intensitySumWithIsotopeModeling(massArr, intArr, isoArr, massTolerance, pepMass);
        if(null != tempArr) {
            spec.setPrecursorPeakIntensity((long) tempArr[0]);
        }
        return spec;
    }



    public static edu.scripps.pms.census.labelFree.SpectrumModel labelFreeSpectrumReader(
            double[] isoArr,
            long startPos,
            int diffPos,
            double massTolerance,
            IndexedFile iFile,
            int chargeState,
            Configuration conf,
            double pepMass
    ) throws IOException, CensusIndexOutOfBoundException, Exception
    {

        //RandomAccessFile file = iFile.getFile();
        RandomAccessFile file = iFile.getFile();


        file.seek(startPos);

        int byteSize = diffPos;
        byte[] bytes = new byte[byteSize];
        file.readFully(bytes);
        String [] debugString = debugGetString(bytes);




        // read header
        int startLinePos =0;
        int endLinePos =0;
        int scan =-1;
        double retTime=-1;
        StringBuilder sb = new StringBuilder(10);

        headerLoop: while(startLinePos<byteSize)
        {
            char firstChar = (char)bytes[startLinePos];

            int pos;
            char currentChar;
            sb.delete(0,10);
            switch (firstChar){
                case 'S':
                    //scan = readSline(bytes,startLinePos);
                    if(startLinePos!=0)
                    {
                        System.out.println(iFile.getFileName());
                        System.exit(1);
                    }
                     pos = startLinePos+2;
                     currentChar =(char)bytes[pos++] ;
                    while((currentChar!=WINDOW_CR && currentChar != CARRIAGE_RETURN && currentChar!=TAB ))
                    {
                        sb.append(currentChar);
                        currentChar = (char)bytes[pos++];
                    }
                    while((currentChar!=WINDOW_CR && currentChar != CARRIAGE_RETURN ))
                    {
                        currentChar = (char)bytes[pos++];
                    }
                    endLinePos = pos;
                    scan = Integer.parseInt(sb.toString());
                    break;
                case 'I':
                    //if(isRetTimeLine(bytes,startLinePos)) retTime = readILineRet(bytes,startLinePos);
                     pos = startLinePos;
                    boolean isRetString = true;
                    for(int i=0;i<ILINERETSTRING.length(); i++)
                    {
                        if(ILINERETSTRING.charAt(i)!=(char)bytes[pos++])
                            isRetString = false;
                    }
                    if(!isRetString)
                    {
                        endLinePos = findEndOfLine(bytes,startLinePos)+1;
                        break;
                    }
                    pos++;
                     currentChar =(char)bytes[pos++] ;
                    while((currentChar!=WINDOW_CR && currentChar != CARRIAGE_RETURN && currentChar!=TAB && currentChar!=' '))
                    {
                        sb.append(currentChar);
                        currentChar = (char)bytes[pos++];
                    }

                    while((currentChar!=WINDOW_CR && currentChar != CARRIAGE_RETURN ))
                    {
                        currentChar = (char)bytes[pos++];
                    }
                    endLinePos = pos;
                    retTime = Double.parseDouble(sb.toString());
                    break;
                default:
                    if(Character.isDigit(firstChar))
                        break headerLoop;
                    else
                        endLinePos = findEndOfLine(bytes,startLinePos);
                    break;
            }
            startLinePos = endLinePos;
        }

        int arrSize=0;
        for(int j=startLinePos;j<byteSize;j++)
        {
            if( CARRIAGE_RETURN == (char)bytes[j] )
                arrSize++;

        }
        char ch;
        double[] massArr = new double[arrSize];
        double[] intArr = new double[arrSize];
        //pos++;
        StringBuilder mass = new StringBuilder(10);
        StringBuilder intensity = new StringBuilder(10);
        intensity.append('0');

        boolean isMass=true;
        boolean isInt=true;

        int massIndex=0;
        //double tempMass;
        int spaceCount=0;
        int cs =0;
        for(int j=startLinePos;j<byteSize;j++)
        {
            ch = (char)bytes[j];
      //      System.out.print(ch);
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

                    //try {
                    if(conf.isLabelfreeCheckChargeState() && spaceCount>=2) {
                        if(chargeState==cs) {
                            massArr[massIndex] =Double.parseDouble(mass.toString());

try {
                            intArr[massIndex] = Long.parseLong(intensity.toString());
} catch (Exception e) {
//intensity can be negative.  assign zero value
intArr[massIndex] = 0;
}
                        }
                    } else {
                        massArr[massIndex] =Double.parseDouble(mass.toString());
                        intArr[massIndex] = Long.parseLong(intensity.toString());
                    }
                    spaceCount=0;

                    massIndex++;
                    mass.delete(0, mass.length());  //This is faster than creating new StringBuilder object
                    intensity.delete(0, intensity.length()).append('0');
                    cs =0;
                    //} catch(Exception e) { System.out.println("-->>" + mass.toString() + " " + intensity.toString()); System.exit(0); };

                    break;
                    /*
                case 'S':
                    System.out.println(iFile.getFileName());
                    System.out.println("Scan Number"+scan);
                    for(int i=0; i<massArr.length;i++)
                    {
                        System.out.println(massArr[i]+"\t"+intArr[i]);
                    }
                    System.exit(1);*/
                case DOT:
                    isInt=false;

                default:
                    if(spaceCount<2)
                    {
                    if (isMass)
                        mass.append(ch);
                    else if (isInt) //remove decimal value of intensity
                        intensity.append(ch);
                    }
                    else if (spaceCount ==2)
                    {
                        cs = ch - '0';
                    }
                    break;

            }
        }
/*
        System.out.print(iFile.getFileName());
        System.out.println("Scan Number"+scan);
        for(int i=0; i<massArr.length;i++)
        {
            System.out.println(massArr[i]+"\t"+intArr[i]);
        }
/*
        System.exit(1);*/

        edu.scripps.pms.census.labelFree.SpectrumModel spec = new edu.scripps.pms.census.labelFree.SpectrumModel();
        // spec.setMass(massArr);
        // spec.setIntensity(intArr);
        spec.setRetentionTime(retTime);
        spec.setScanNumber(scan);
/*
        if(scan == 26391)
        {
            System.out.println(iFile.getFileName());
            System.out.println("Scan Number"+scan);
            for(int i=0; i<massArr.length;i++)
            {
                System.out.println(massArr[i]+"\t"+intArr[i]);
            }
            System.exit(1);
        }

/*
        int scanCheck = Integer.parseInt(debugString[0].split("\t")[1]);
        double retCheck = Double.parseDouble(debugString[1].split("\t")[2]);
        //double massCheck = Double.parseDouble(debugString[4].split(" ")[0]);
        //int intensityCheck = Integer.parseInt(debugString[4].split(" ")[1]);

        if(scan!=scanCheck || (Double.compare(retCheck,retTime)!=0) )
        {
            System.out.print(iFile.getFileName());
            System.out.println("Scan Number"+scan);
            for(int i=0; i<massArr.length;i++)
            {
                System.out.println(massArr[i]+"\t"+intArr[i]);
            }
            System.exit(1);
        }
        else
        {
            for(int i=4;i<debugString.length; i++)
            {
                String [] arr = debugString[i].split(" ");
                double massCheck = Double.parseDouble(arr[0]);
                int intensityCheck = Integer.parseInt(arr[1]);
                if(Double.compare(massCheck,massArr[i-4])!=0 || intensityCheck!=intArr[i-4])
                {
                    System.out.print(iFile.getFileName());
                    System.out.println("Scan Number"+scan);
                    for(int j=0; j<massArr.length;j++)
                    {
                        System.out.println(massArr[j]+"\t"+intArr[j]);
                    }
                    System.exit(1);
                }
            }
        }
*/


        //result[0] = keys[curIndex];
        double[] tempArr = intensitySumWithIsotopeModeling(massArr, intArr, isoArr, massTolerance, pepMass);
        if(null != tempArr) {
            spec.setPrecursorPeakIntensity((long) tempArr[0]);
            //System.out.println(iFile.getFileName()+"\t"+scan+"\t"+(long) tempArr[0]);
            //System.exit(1);
        }
        return spec;
    }

    public static int findEndOfLine(byte [] arr, int start)
    {
        for(int i= start; i<arr.length; i++)
        {
            if(arr[i]==WINDOW_CR || arr[i]==CARRIAGE_RETURN)
                return i;
        }
        return arr.length-1;
    }


    public static int readSline(byte[] arr, int start)
    {
        int pos =0;
        char currentChar =(char)arr[pos++] ;
        StringBuilder sb = new StringBuilder();
        while(currentChar!=WINDOW_CR || currentChar != CARRIAGE_RETURN)
        {
            sb.append(currentChar);
            currentChar = (char)arr[pos++];
        }
        return Integer.parseInt(sb.toString());
    }
    public static boolean isRetTimeLine(byte[] arr, int start)
    {
        for(int i=0;i<ILINERETSTRING.length(); i++)
        {
            if(ILINERETSTRING.charAt(i)!=(char)arr[start++])
                return false;
        }
        return true;
    }
    public static double readILineRet(byte[] arr, int start)
    {
        int pos = start;
        char currentChar =(char)arr[pos++] ;
        StringBuilder sb = new StringBuilder();


        while(currentChar!=WINDOW_CR || currentChar != CARRIAGE_RETURN)
        {
            if(Character.isDigit(currentChar))
            {
                sb.append(currentChar);
            }
            currentChar = (char)arr[pos++];
        }
        return Double.parseDouble(sb.toString());
    }



    /*
    public static double[][] getSpectrumArr(int[] keys, int curIndex, TIntLongHashMap index, RandomAccessFile file) throws IOException, CensusIndexOutOfBoundException
    {

    }
*/

    public static double[] intensitySumWithIsotopeModeling(double[] massArr, double[] intArr, double[] isoArr, double massTolerance, double pepMass)
    {

        if(null == massArr || massArr.length<=0) return null;

        double[] resultArr = new double[3];

        double sumIntensity=0;
        int foundIsoNum=0;
        int totalPeakFound=0;
        double toleranceSum=0;

        double[] averagineArr = AssignMass.getIsotopeModel(pepMass);
        double[] peakArr = new double[averagineArr.length];

        for(int i=0;i<isoArr.length;i++)
        {
                //double tempTolerance = isoArr[i]/1000*massTolerance;
                //double tempTolerance = isoArr[i]/1000*massTolerance*0.001;
            double tempTolerance = isoArr[i]/1000*massTolerance;
            int start = Arrays.binarySearch(massArr, isoArr[i]-tempTolerance);
            if(start<0)
            start = -start -1;

            //int j=0;
            double small=100;
            //double massSmall=-1;
           // double highinten=0;

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

                    if(i<peakArr.length)
                        peakArr[i] += intArr[start];


                    totalPeakFound++;
                    toleranceSum += temp;

                    isFound = true;

                    double diff = (massArr[start] - isoArr[i]);
                    //if(diff<0)
                    //diff = -diff;
                    if(Math.abs(diff)<Math.abs(small))
                    {
                        small = massArr[start] - isoArr[i];
                //	massSmall = massArr[start];
                //	highinten = intArr[start];
                    }

                }

                if(massArr[start]>isoArr[i])
                    break;

                start++;
            }

                if(isFound)
                    foundIsoNum++;
        }

      //  LinearRegressionDouble

        SimpleRegression sr = new SimpleRegression();

     //   System.out.println("==============================================");
        for (int i = 0; i < averagineArr.length; i++) {
            sr.addData(averagineArr[i], peakArr[i]);

        //    System.out.println("corr====\t" + peakArr[i] + "\t\t" +  averagineArr[i]);
        }


   //     if( Double.isNaN(sr.getR()) || Double.isInfinite(sr.getR()) || sr.getR()<0.7)

            //System.out.println("uncomment this line...............................\t" + sr.getR() + "\t0");
     //       System.out.println("uncomment this line...............................\t" + "0");
    //    else
     //       System.out.println("uncomment this line...............................\t" + sumIntensity + "\t" + sr.getR());

        if( Double.isNaN(sr.getR()) || Double.isInfinite(sr.getR()) || sr.getR()<0.7)
            return resultArr;

            resultArr[0] = sumIntensity;
            resultArr[1] = foundIsoNum;
            resultArr[2] = (totalPeakFound>0)?(toleranceSum/totalPeakFound):(-1);

        return resultArr;
    }

    public static void main(String[] args) throws Exception {
               //High resolution //lalebling data
        /*
        public static edu.scripps.pms.census.labelFree.SpectrumModel readLabelfreeFullSpectrum(
            double[] isoArr,
            long startPos,
            int diffPos,
            double massTolerance,
            IndexedFile iFile
*/



        //Hashtable<String, IndexedFile> ht = createIndexedFiles("/Users/rpark/test_data", CensusConstants.MS2_FILE);
        Hashtable<String, IndexedFile> ht = createIndexedFiles("/home/rpark/test_data/sherry", CensusConstants.MS2_FILE);


        //IndexedFile iFile = ht.get("/Users/rpark/test_data/1993804.ms2");

        IndexedFile iFile = ht.get("/home/rpark/test_data/sherry/test.ms2");

         int[] keys = iFile.getKeys();


       // BufferedRandomAccessFile file = new BufferedRandomAccessFile("/Users/rpark/test_data/1993804.ms2", "r", 2000);
        BufferedRandomAccessFile file = new BufferedRandomAccessFile("/home/rpark/test_data/sherry/test.ms2", "r", 2000);

        //if(true) return;
//        Hashtable<String, IndexedFile> indexHt = new Hashtable<>();
        //Hashtable<String, IndexedFile> ht = ChroGenerator.createIndexedFiles("/data/2/rpark/ip2_data//xudong/AC_mousebrain_exp2_3_repeats_separate/sample1_3_2016_05_15_20_82232/search/upload2016_05_15_20_97593//", "ms2", false);
      //  File f = new File("/data/2/rpark/ip2_data//xudong/AC_mousebrain_exp2_3_repeats_separate/sample1_3_2016_05_15_20_82232/search/upload2016_05_15_20_97593//1993804.ms2.index");
//	IndexedFile iFile = new IndexedFile(f, "/data/2/rpark/ip2_data//xudong/AC_mousebrain_exp2_3_repeats_separate/sample1_3_2016_05_15_20_82232/search/upload2016_05_15_20_97593/1993804.ms2");
//"/Users/rpark/test_data"
	double[] isoArr = new double[4];
	isoArr[0] = 700.324;
	isoArr[1] = 701.324;
	isoArr[2] = 702.324;
	isoArr[3] = 703.324;

	CalcUtilGeneric cal = new CalcUtilGeneric();
	gnu.trove.TIntLongHashMap map = iFile.getScanPositionMap();

     //   for(int i=1;i<100000;i++) {
        for(int i=1;i<10;i++) {
            long start =map.get(keys[i]);
            long next = map.get(keys[i+1]);
            int diff = (int) (next- start);
//            System.out.print(".");
//System.out.println(start + " " + next + " " + diff);
/*
            int byteSize = diff;

            byte[] bytes = new byte[byteSize];
            file.readFully(bytes);
*/

           SpectrumModel spectrumModel = cal.readLabelfreeFullSpectrumFast(isoArr, start, diff, 0.05, iFile, file);

          //  cal.readLabelfreeFullSpectrumFast(isoArr, start, diff, 0.05, iFile);
        }

        System.out.println("");
        file.close();
       // iFile.get

    }

}
