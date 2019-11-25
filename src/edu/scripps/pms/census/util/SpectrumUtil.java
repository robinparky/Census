package edu.scripps.pms.census.util;

import edu.scripps.pms.census.ChroGenerator;
import edu.scripps.pms.census.exception.CensusGeneralException;
import edu.scripps.pms.census.exception.CensusIndexOutOfBoundException;
import edu.scripps.pms.census.hash.IndexedFile;
import edu.scripps.pms.util.sqlite.spectra.SpectraDB;
import gnu.trove.TDoubleArrayList;

import java.io.*;
import java.util.*;
import java.io.BufferedReader;

import static edu.scripps.pms.census.util.CalcUtil.intArr;


public class SpectrumUtil
{

	public static final double MZTHESHOLD = 0.10;
	public static void main(String args[]) throws Exception
	{

		//System.out.println( SpectrumUtil.getNumOfIsolationWindow(args[0]) );
		String file=args[0];
		String spath = file.substring(0, file.lastIndexOf(File.separator));
		String filename = file.substring(file.lastIndexOf(File.separator)+1);

		//List<XYPoint> massIntensity=getSpectrumMS1Arr(spath,filename,2093);
		List<XYPoint> massIntensity = getSpectrumMS1(spath,filename,90,470,480);
		//List<XYPoint> massIntensity=getSpectrum(spath,filename,1,300,400);
		for(XYPoint p: massIntensity)
		{
			System.out.println("Label: "+p.getLabel()+" x: "+p.getX()+" y: "+p.getY());
		}
		System.out.println("List size is: "+massIntensity.size());
	}

    public static int getNumOfIsolationWindow(String fileName) throws IOException
    {
        BufferedReader br = null;
	String lastLine;


	int i=1;
	try
	{
        	br = new BufferedReader(new FileReader(fileName));

        while( (lastLine=br.readLine()) != null && !lastLine.startsWith("S") );

	String start = lastLine.substring(lastLine.lastIndexOf("\t")+1);
//	System.out.println(start);


        	while( (lastLine=br.readLine()) != null )
		{
			if( lastLine.startsWith("S") )
			{
				i++;
			 //System.out.println(lastLine);
				if( start.equals(lastLine.substring(lastLine.lastIndexOf("\t")+1)) )
				break;
			}
		}			
	

	}
	catch(IOException e)
	{
		System.out.println("Error : " + e);
		throw new IOException(e.toString());
	
	}
	finally {
		
		if(null != br)
        		br.close();
	}


	return i;
    }


	public static List<edu.scripps.pms.census.model.XYPoint> getSpectrum(String spath, SpectraDB spectraDB, int scanNum, double startMass,
																		 double endMass)
			throws IOException, CensusGeneralException, CensusIndexOutOfBoundException {



		List<edu.scripps.pms.census.model.XYPoint> massInt = new ArrayList<>();

		try {
			SpectraDB.Spectrum spectrum  = spectraDB.getSpectrumFromDB(scanNum);
			intArr = spectrum.getIntensityList().toNativeArray();
			TDoubleArrayList massArr = spectrum.getMzList();
			edu.scripps.pms.census.model.XYPoint xy = null;

			for (int i = 0; i < intArr.length; i++) {
				if(massArr.get(i) >= startMass && massArr.get(i) <=endMass){
					xy =  new edu.scripps.pms.census.model.XYPoint();
					xy.setX(massArr.get(i));
					xy.setY(intArr[i]);
					massInt.add(xy);
				}
			}


		} catch (Exception e) {
			e.printStackTrace();
		}
		return massInt;
	}



	public static List<edu.scripps.pms.census.model.XYPoint> getSpectrum(String spath, String fileName, int scanNum, double startMass,
																		 double endMass)
			throws IOException, CensusGeneralException, CensusIndexOutOfBoundException {



		List<edu.scripps.pms.census.model.XYPoint> massInt = new ArrayList<>();

		try {

			String path = spath + File.separator + fileName;

			String spectrumPath = spath + File.separator + fileName + ".index";
			String spectrumPathMs2 = spath + File.separator + fileName;
			BufferedReader br = null;
			RandomAccessFile rfile = null;
			br = new BufferedReader(new FileReader(spectrumPath));

			String eachLine;
			while ((eachLine = br.readLine()) != null) {
				if (eachLine.startsWith(String.valueOf(scanNum)))
					break;
			}

			String[] arr = eachLine.split("\t");
			String nextLine = br.readLine();

			rfile = new RandomAccessFile(spectrumPathMs2, "r");

			long startPos = Long.parseLong(arr[1]);
			rfile.seek(startPos);
			int byteSize = -1;

			if (null == nextLine) {
				byteSize = (int) (rfile.length() - startPos);
			} else {
				String[] nextArr = nextLine.split("\t");
				long nextPos = Long.parseLong(nextArr[1]);
				byteSize = (int) (nextPos - startPos);
			}

			byte[] bytes = new byte[byteSize];
			rfile.readFully(bytes);

			int pos = 0;
			char ch = (char) bytes[pos];
			while ((ch = (char) bytes[pos]) == 'Z' || ch == 'S' || ch == 'I' || ch == 'D') {
				while (ch != '\n') {
					pos++;
					ch = (char) bytes[pos];
				}

				pos++;
			}
			int arrSize = 0;
			for (int j = pos; j < byteSize; j++) {
				if ('\n' == (char) bytes[j])
					arrSize++;

			}

			double[] massArr = new double[arrSize];
			double[] intArr = new double[arrSize];

			StringBuilder mass = new StringBuilder(10);
			StringBuilder intensity = new StringBuilder(10);
			intensity.append('0');

			boolean isMass = true;
			boolean isInt = true;
			int countInt = 0;
			int massIndex = 0;
			for (int j = pos; j < byteSize; j++) {
				ch = (char) bytes[j];

				switch (ch) {
					case '\r':
						break;

					case ' ':
						countInt++;
						isMass = false;
						isInt = countInt <= 1 ? true : false;
						break;

					case '\n':
						isMass = true;
						isInt = false;
						intArr[massIndex] = Long.parseLong(intensity.toString());
						massArr[massIndex++] = Double.parseDouble(mass.toString());
						mass.delete(0, mass.length());
						intensity.delete(0, intensity.length()).append('0');
						countInt = 0;
						break;

					case '.':
						isInt = false;

					default:
						if (isMass)
							mass.append(ch);
						else if (isInt)
							intensity.append(ch);

						break;
				}

			}

			edu.scripps.pms.census.model.XYPoint xy = null;

			for (int i = 0; i < intArr.length; i++) {
				if(massArr[i] >= startMass && massArr[i] <=endMass){
					xy =  new edu.scripps.pms.census.model.XYPoint();
					xy.setX(massArr[i]);
					xy.setY(intArr[i]);
					massInt.add(xy);
				}
			}

			rfile.close();
			br.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		return massInt;
	}

	public static List<XYPoint> getSpectrum(String spath,String fileName,int scanNum,int start_pos,int end_pos) throws IOException, CensusGeneralException, CensusIndexOutOfBoundException {
		String path = spath + File.separator + fileName;


		String spectrumPath = spath +  File.separator + fileName + ".index";
		String spectrumPathMs2 = spath + File.separator  + fileName;
		BufferedReader br = null;
		RandomAccessFile rfile = null;
		br = new BufferedReader(new FileReader(spectrumPath));

		String eachLine;
		while( (eachLine = br.readLine()) != null) {
			if(eachLine.startsWith(String.valueOf(scanNum)))
				break;
		}


		String[] arr = eachLine.split("\t");
		String nextLine = br.readLine();

		rfile = new RandomAccessFile(spectrumPathMs2, "r");

		long startPos = Long.parseLong(arr[1]);
		rfile.seek(startPos);
		int byteSize=-1;

		if(null == nextLine) {
			byteSize = (int)(rfile.length() - startPos);
		} else {
			String[] nextArr = nextLine.split("\t");
			long nextPos = Long.parseLong(nextArr[1]);
			byteSize = (int)(nextPos - startPos);
		}

		byte[] bytes = new byte[byteSize];
		rfile.readFully(bytes);

		int pos = 0;
		char ch = (char)bytes[pos];
		while( (ch=(char)bytes[pos]) == 'Z' || ch == 'S' || ch == 'I' || ch == 'D')
		{
			while( ch != '\n' )
			{
				pos++;
				ch = (char)bytes[pos];
			}

			pos++;
		}
		int arrSize=0;
		for(int j=pos;j<byteSize;j++)
		{
			if( '\n' == (char)bytes[j] )
				arrSize++;

		}

		double[] massArr = new double[arrSize];
		double[] intArr = new double[arrSize];


		StringBuilder mass = new StringBuilder(10);
		StringBuilder intensity = new StringBuilder(10);
		intensity.append('0');

		boolean isMass=true;
		boolean isInt=true;
		int countInt=0;
		int massIndex=0;
		boolean done=false;
		int len=0;
		for(int j=pos;j<byteSize;j++)
		{
			ch = (char)bytes[j];

			switch(ch)
			{
				case '\r':
					break;

				case ' ':
					countInt++;
					isMass=false;
					isInt=countInt<=1?true:false;
					break;

				case '\n':
					isMass=true;
					isInt=false;
					if(Double.parseDouble(mass.toString())>start_pos && Double.parseDouble(mass.toString())<end_pos){
						intArr[massIndex] = Long.parseLong(intensity.toString());
						massArr[massIndex++] = Double.parseDouble(mass.toString());
						len++;
					}
					if(Double.parseDouble(mass.toString())>end_pos){
						done=true;
					}
					mass.delete(0, mass.length());
					intensity.delete(0, intensity.length()).append('0');
					countInt=0;
					break;

				case '.':
					isInt=false;

				default:
					if(isMass)
						mass.append(ch);
					else if(isInt)
						intensity.append(ch);

					break;
			}
			if(done){
				break;
			}

		}

		List<XYPoint> massInt=new ArrayList<>();

		for(int i=0;i<len;i++){
			XYPoint xy=new XYPoint();
			xy.setX(massArr[i]);
			xy.setY(intArr[i]);
			massInt.add(xy);
		}

		rfile.close();
		br.close();

		return massInt;
	}


	public static List<XYPoint> getSpectrum(String spath,String fileName,int scanNum) throws IOException, CensusGeneralException, CensusIndexOutOfBoundException {
		String path = spath + File.separator + fileName;

		String spectrumPath = spath +  File.separator + fileName + ".index";
		String spectrumPathMs2 = spath + File.separator  + fileName;

		BufferedReader br = null;
		RandomAccessFile rfile = null;
		br = new BufferedReader(new FileReader(spectrumPath));

		String eachLine;
		while( (eachLine = br.readLine()) != null) {
			if(eachLine.startsWith(String.valueOf(scanNum)))
				break;
		}


		String[] arr = eachLine.split("\t");
		String nextLine = br.readLine();

		rfile = new RandomAccessFile(spectrumPathMs2, "r");

		long startPos = Long.parseLong(arr[1]);
		rfile.seek(startPos);
		int byteSize=-1;

		if(null == nextLine) {
			byteSize = (int)(rfile.length() - startPos);
		} else {
			String[] nextArr = nextLine.split("\t");
			long nextPos = Long.parseLong(nextArr[1]);
			byteSize = (int)(nextPos - startPos);
		}

		byte[] bytes = new byte[byteSize];
		rfile.readFully(bytes);

		int pos = 0;
		char ch = (char)bytes[pos];
		while( (ch=(char)bytes[pos]) == 'Z' || ch == 'S' || ch == 'I' || ch == 'D')
		{
			while( ch != '\n' )
			{
				pos++;
				ch = (char)bytes[pos];
			}

			pos++;
		}
		int arrSize=0;
		for(int j=pos;j<byteSize;j++)
		{
			if( '\n' == (char)bytes[j] )
				arrSize++;

		}

		double[] massArr = new double[arrSize];
		double[] intArr = new double[arrSize];


		StringBuilder mass = new StringBuilder(10);
		StringBuilder intensity = new StringBuilder(10);
		intensity.append('0');

		boolean isMass=true;
		boolean isInt=true;
		int countInt=0;
		int massIndex=0;
		for(int j=pos;j<byteSize;j++)
		{
			ch = (char)bytes[j];

			switch(ch)
			{
				case '\r':
					break;

				case ' ':
					countInt++;
					isMass=false;
					isInt=countInt<=1?true:false;
					break;

				case '\n':
					isMass=true;
					isInt=false;
					intArr[massIndex] = Long.parseLong(intensity.toString());
					massArr[massIndex++] = Double.parseDouble(mass.toString());
					mass.delete(0, mass.length());
					intensity.delete(0, intensity.length()).append('0');
					countInt=0;
					break;

				case '.':
					isInt=false;

				default:
					if(isMass)
						mass.append(ch);
					else if(isInt)
						intensity.append(ch);

					break;
			}

		}

		List<XYPoint> massInt=new ArrayList<>();

		for(int i=0;i<intArr.length;i++){
			XYPoint xy=new XYPoint();
			xy.setX(massArr[i]);
			xy.setY(intArr[i]);
			massInt.add(xy);
		}

		rfile.close();
		br.close();

		return massInt;
	}

	public static List<XYPoint> getSpectrumMS1(String spath,String fileName,int scanNum) throws IOException, CensusGeneralException, CensusIndexOutOfBoundException {
		String path = spath + File.separator + fileName;

		String fileNameWithoutExtension=fileName.substring(0,fileName.lastIndexOf("."));
		String spectrumPath = spath +  File.separator + fileName + ".index";

		BufferedReader br = null;
		RandomAccessFile rfile = null;
		br = new BufferedReader(new FileReader(spectrumPath));
		Hashtable<String, IndexedFile> ht=ChroGenerator.createIndexedFiles(spath, fileName);
		int  prevScanfromMS1=getMS1Scan(ht,path,scanNum);

		String eachLine;
		while( (eachLine = br.readLine()) != null) {
			if(eachLine.startsWith(String.valueOf(prevScanfromMS1)))
				break;
		}


		String[] arr = eachLine.split("\t");
		String nextLine = br.readLine();

		rfile = new RandomAccessFile(path, "r");

		long startPos = Long.parseLong(arr[1]);
		rfile.seek(startPos);
		int byteSize=-1;

		if(null == nextLine) {
			byteSize = (int)(rfile.length() - startPos);
		} else {
			String[] nextArr = nextLine.split("\t");
			long nextPos = Long.parseLong(nextArr[1]);
			byteSize = (int)(nextPos - startPos);
		}

		byte[] bytes = new byte[byteSize];
		rfile.readFully(bytes);

		int pos = 0;
		char ch = (char)bytes[pos];
		while( (ch=(char)bytes[pos]) == 'Z' || ch == 'S' || ch == 'I' || ch == 'D')
		{
			while( ch != '\n' )
			{
				pos++;
				ch = (char)bytes[pos];
			}

			pos++;
		}
		int arrSize=0;
		for(int j=pos;j<byteSize;j++)
		{
			if( '\n' == (char)bytes[j] )
				arrSize++;

		}

		double[] massArr = new double[arrSize];
		double[] intArr = new double[arrSize];


		StringBuilder mass = new StringBuilder(10);
		StringBuilder intensity = new StringBuilder(10);
		intensity.append('0');

		boolean isMass=true;
		boolean isInt=true;
		int countInt=0;
		int massIndex=0;
		for(int j=pos;j<byteSize;j++)
		{
			ch = (char)bytes[j];

			switch(ch)
			{
				case '\r':
					break;

				case ' ':
					countInt++;
					isMass=false;
					isInt=countInt<=1?true:false;
					break;

				case '\n':
					isMass=true;
					isInt=false;
					intArr[massIndex] = Long.parseLong(intensity.toString());
					massArr[massIndex++] = Double.parseDouble(mass.toString());
					mass.delete(0, mass.length());
					intensity.delete(0, intensity.length()).append('0');
					countInt=0;
					break;

				case '.':
					isInt=false;

				default:
					if(isMass)
						mass.append(ch);
					else if(isInt)
						intensity.append(ch);

					break;
			}

		}

		List<XYPoint> massInt=new ArrayList<>();

		for(int i=0;i<intArr.length;i++){
			XYPoint xy=new XYPoint();
			xy.setX(massArr[i]);
			xy.setY(intArr[i]);
			massInt.add(xy);
		}

		rfile.close();
		br.close();

		return massInt;
	}

	public static List<XYPoint> getSpectrumMS1(String spath,String fileName,int scanNum,int start_pos,int end_pos) throws IOException, CensusGeneralException, CensusIndexOutOfBoundException {
		String path = spath + File.separator + fileName;

		String fileNameWithoutExtension=fileName.substring(0,fileName.lastIndexOf("."));
		String spectrumPath = spath +  File.separator + fileName + ".index";
		String spectrumPathMs2 = spath + File.separator  + fileName;

		BufferedReader br = null;
		RandomAccessFile rfile = null;
		br = new BufferedReader(new FileReader(spectrumPath));

		Hashtable<String, IndexedFile> ht=ChroGenerator.createIndexedFiles(spath, fileName);
		int  prevScanfromMS1=getMS1Scan(ht,path,scanNum);


		String eachLine;
		while( (eachLine = br.readLine()) != null) {
			if(eachLine.startsWith(String.valueOf(prevScanfromMS1)))
				break;
		}


		String[] arr = eachLine.split("\t");
		String nextLine = br.readLine();

		rfile = new RandomAccessFile(spectrumPathMs2, "r");

		long startPos = Long.parseLong(arr[1]);
		rfile.seek(startPos);
		int byteSize=-1;

		if(null == nextLine) {
			byteSize = (int)(rfile.length() - startPos);
		} else {
			String[] nextArr = nextLine.split("\t");
			long nextPos = Long.parseLong(nextArr[1]);
			byteSize = (int)(nextPos - startPos);
		}

		byte[] bytes = new byte[byteSize];
		rfile.readFully(bytes);

		int pos = 0;
		char ch = (char)bytes[pos];
		while( (ch=(char)bytes[pos]) == 'Z' || ch == 'S' || ch == 'I' || ch == 'D')
		{
			while( ch != '\n' )
			{
				pos++;
				ch = (char)bytes[pos];
			}

			pos++;
		}
		int arrSize=0;
		for(int j=pos;j<byteSize;j++)
		{
			if( '\n' == (char)bytes[j] )
				arrSize++;

		}

		double[] massArr = new double[arrSize];
		double[] intArr = new double[arrSize];


		StringBuilder mass = new StringBuilder(10);
		StringBuilder intensity = new StringBuilder(10);
		intensity.append('0');
        boolean done=false;
		boolean isMass=true;
		boolean isInt=true;
		int countInt=0;
		int massIndex=0;
		int len=0;
		for(int j=pos;j<byteSize;j++)
		{
			ch = (char)bytes[j];

			switch(ch)
			{
				case '\r':
					break;

				case ' ':
					countInt++;
					isMass=false;
					isInt=countInt<=1?true:false;
					break;

				case '\n':
					isMass=true;
					isInt=false;
					if(Double.parseDouble(mass.toString())>start_pos && Double.parseDouble(mass.toString())<end_pos){
						intArr[massIndex] = Long.parseLong(intensity.toString());
						massArr[massIndex++] = Double.parseDouble(mass.toString());
						len++;
					}
					if(Double.parseDouble(mass.toString())>end_pos){
						done=true;
				    }
					mass.delete(0, mass.length());
					intensity.delete(0, intensity.length()).append('0');
					countInt=0;
					break;

				case '.':
					isInt=false;

				default:
					if(isMass)
						mass.append(ch);
					else if(isInt)
						intensity.append(ch);

					break;
			}

			if(done){
				break;
			}

		}

		List<XYPoint> massInt=new ArrayList<>();

		for(int i=0;i<len;i++){
			XYPoint xy=new XYPoint();
			xy.setX(massArr[i]);
			xy.setY(intArr[i]);
			massInt.add(xy);
		}

		rfile.close();
		br.close();

		return massInt;
	}
	public static void getSpectrumMS1Arr(String spath, String fileName, int scanNum, double start_pos, double end_pos,
										 int prevScan, TDoubleArrayList massArr, TDoubleArrayList intArr) throws IOException, CensusGeneralException, CensusIndexOutOfBoundException {
		String path = spath + File.separator + fileName;
		String fileNameWithoutExtension=fileName.substring(0,fileName.lastIndexOf("."));
		String spectrumPath = spath +  File.separator + fileName + ".index";
		String spectrumPathMs2 = spath + File.separator  + fileName;

		BufferedReader br = null;
		RandomAccessFile rfile = null;
		br = new BufferedReader(new FileReader(spectrumPath));

		int  prevScanfromMS1=prevScan;

	//	System.out.println(">>>> " +path);
	//	System.out.println(">>>> " +prevScanfromMS1);
		String eachLine;
		while( (eachLine = br.readLine()) != null) {
			if(eachLine.startsWith(String.valueOf(prevScanfromMS1)))
				break;
		}
		if(eachLine==null)
		{
			massArr = null;
			intArr = null;
			return;
		}


		String[] arr = eachLine.split("\t");
		String nextLine = br.readLine();

		rfile = new RandomAccessFile(spectrumPathMs2, "r");

		long startPos = Long.parseLong(arr[1]);
		rfile.seek(startPos);
		int byteSize=-1;

		if(null == nextLine) {
			byteSize = (int)(rfile.length() - startPos);
		} else {
			String[] nextArr = nextLine.split("\t");
			long nextPos = Long.parseLong(nextArr[1]);
			byteSize = (int)(nextPos - startPos);
		}

		byte[] bytes = new byte[byteSize];
		rfile.readFully(bytes);

		int pos = 0;
		char ch = (char)bytes[pos];
		while( (ch=(char)bytes[pos]) == 'Z' || ch == 'S' || ch == 'I' || ch == 'D')
		{
			while( ch != '\n' )
			{
				pos++;
				ch = (char)bytes[pos];
			}

			pos++;
		}
		int arrSize=0;
		for(int j=pos;j<byteSize;j++)
		{
			if( '\n' == (char)bytes[j] )
				arrSize++;

		}

		//double []massArr = new double[arrSize];
		//double []intArr = new double[arrSize];


		StringBuilder mass = new StringBuilder(10);
		StringBuilder intensity = new StringBuilder(10);
		intensity.append('0');
		boolean done=false;
		boolean isMass=true;
		boolean isInt=true;
		int countInt=0;
		int massIndex=0;
		int len=0;
		long maxIntensity = 0;
		for(int j=pos;j<byteSize;j++)
		{
			ch = (char)bytes[j];

			switch(ch)
			{
				case '\r':
					break;

				case ' ':
					countInt++;
					isMass=false;
					isInt=countInt<=1?true:false;
					break;

				case '\n':
					isMass=true;
					isInt=false;



					///mass in m/z
					double mz = Double.parseDouble(mass.toString());
					if(mz>start_pos && mz<end_pos){
						long ii = Long.parseLong(intensity.toString());
						massArr.add(Double.parseDouble(mass.toString()));
						intArr.add(ii);
					//	intArr[massIndex] = ii;
					//	massArr[massIndex++] = Double.parseDouble(mass.toString());
						if(ii>maxIntensity) maxIntensity =ii;
						len++;
					}
					if(Double.parseDouble(mass.toString())>end_pos){
						done=true;
					}
					mass.delete(0, mass.length());
					intensity.delete(0, intensity.length()).append('0');
					countInt=0;
					break;

				case '.':
					isInt=false;

				default:
					if(isMass)
						mass.append(ch);
					else if(isInt)
						intensity.append(ch);

					break;
			}

			if(done){
				break;
			}

		}

		//List<XYPoint> massInt=new ArrayList<>();
		//	long threshold = (long)(maxIntensity*MZTHESHOLD);


		rfile.close();
		br.close();
		//return new double[][]{massArr,intArr};
		//return massInt;
	}


	public static List<XYPoint> getSpectrumMS1(String spath,String fileName,int scanNum,double start_pos,double end_pos,int prevScan) throws IOException, CensusGeneralException, CensusIndexOutOfBoundException {
		String path = spath + File.separator + fileName;
		String fileNameWithoutExtension=fileName.substring(0,fileName.lastIndexOf("."));
		String spectrumPath = spath +  File.separator + fileName + ".index";
		String spectrumPathMs2 = spath + File.separator  + fileName;

		BufferedReader br = null;
		RandomAccessFile rfile = null;
		br = new BufferedReader(new FileReader(spectrumPath));

		int  prevScanfromMS1=prevScan;
/*
		System.out.println(">>>> " +path);
		System.out.println(">>>> " +prevScanFromMS1);*/
		String eachLine;
		while( (eachLine = br.readLine()) != null) {
			if(eachLine.startsWith(String.valueOf(prevScanfromMS1)))
				break;
		}


		String[] arr = eachLine.split("\t");
		String nextLine = br.readLine();

		rfile = new RandomAccessFile(spectrumPathMs2, "r");

		long startPos = Long.parseLong(arr[1]);
		rfile.seek(startPos);
		int byteSize=-1;

		if(null == nextLine) {
			byteSize = (int)(rfile.length() - startPos);
		} else {
			String[] nextArr = nextLine.split("\t");
			long nextPos = Long.parseLong(nextArr[1]);
			byteSize = (int)(nextPos - startPos);
		}

		byte[] bytes = new byte[byteSize];
		rfile.readFully(bytes);

		int pos = 0;
		char ch = (char)bytes[pos];
		while( (ch=(char)bytes[pos]) == 'Z' || ch == 'S' || ch == 'I' || ch == 'D')
		{
			while( ch != '\n' )
			{
				pos++;
				ch = (char)bytes[pos];
			}

			pos++;
		}
		int arrSize=0;
		for(int j=pos;j<byteSize;j++)
		{
			if( '\n' == (char)bytes[j] )
				arrSize++;

		}

		double[] massArr = new double[arrSize];
		double[] intArr = new double[arrSize];


		StringBuilder mass = new StringBuilder(10);
		StringBuilder intensity = new StringBuilder(10);
		intensity.append('0');
		boolean done=false;
		boolean isMass=true;
		boolean isInt=true;
		int countInt=0;
		int massIndex=0;
		int len=0;
		long maxIntensity = 0;
		for(int j=pos;j<byteSize;j++)
		{
			ch = (char)bytes[j];

			switch(ch)
			{
				case '\r':
					break;

				case ' ':
					countInt++;
					isMass=false;
					isInt=countInt<=1?true:false;
					break;

				case '\n':
					isMass=true;
					isInt=false;



					///mass in m/z
					double mz = Double.parseDouble(mass.toString());
					if(mz>start_pos && mz<end_pos){
						long ii = Long.parseLong(intensity.toString());
						intArr[massIndex] = ii;
						massArr[massIndex++] = Double.parseDouble(mass.toString());
						if(ii>maxIntensity) maxIntensity =ii;
						len++;
					}
					if(Double.parseDouble(mass.toString())>end_pos){
						done=true;
					}
					mass.delete(0, mass.length());
					intensity.delete(0, intensity.length()).append('0');
					countInt=0;
					break;

				case '.':
					isInt=false;

				default:
					if(isMass)
						mass.append(ch);
					else if(isInt)
						intensity.append(ch);

					break;
			}

			if(done){
				break;
			}

		}

		List<XYPoint> massInt=new ArrayList<>();
	//	long threshold = (long)(maxIntensity*MZTHESHOLD);
		for(int i=0;i<len;i++){
		/*	if(intArr[i]<threshold)
			{
				//System.out.println(">>>> "+massArr[i]+"\t"+intArr[i]+"\tthreshold\t"+threshold+"\t max "+maxIntensity+" discard!");
				continue;
			}*/
			XYPoint xy=new XYPoint();
			xy.setX(massArr[i]);
			xy.setY(intArr[i]);
			massInt.add(xy);
		}

		rfile.close();
		br.close();

		return massInt;
	}

	public static int getMS1FixedScan(Hashtable<String, IndexedFile> ht, String filename,int scan) throws IOException, CensusIndexOutOfBoundException {

		int keyIndex=0;


		IndexedFile file=ht.get(filename);



		return file.getIdScanToMs1ScanMap().get(scan);
	}



	public static int getMS1Scan(Hashtable<String, IndexedFile> ht, String filename,int scan) throws IOException, CensusIndexOutOfBoundException {

		int keyIndex=0;
		int[] keys = null;


		IndexedFile file=ht.get(filename);
		keys = file.getKeys();
		//int temp = file.getPrecScanMap().get(scan);
		keyIndex = Arrays.binarySearch(keys, scan);

		if (keyIndex < 0) //Cannot find index
		{
			keyIndex = -(++keyIndex); //Math.abs(++keyIndex);
		}
		if (keyIndex >= keys.length) {
			keyIndex--;
		}

		if(keyIndex>0) keyIndex--;

		return keys[keyIndex];
	}

	public static void getMs1ScanPeaksArr(String path, String filename, int scannumber, double startMz, double endMz,
										  int preVScanNumber, TDoubleArrayList massArr, TDoubleArrayList intArr) throws CensusGeneralException, CensusIndexOutOfBoundException, IOException {
		 getSpectrumMS1Arr(path,filename,scannumber,startMz,endMz,preVScanNumber,massArr,intArr);
	}

	public static List<XYPoint> getMs1ScanPeaks(String path, String filename, int scannumber, double startMz, double endMz, int preVScanNumber) throws CensusGeneralException, CensusIndexOutOfBoundException, IOException {

		return getSpectrumMS1(path,filename,scannumber,startMz,endMz,preVScanNumber);
		/*List<XYPoint> preresult = getSpectrumMS1Arr(path,filename,scannumber);
		List<XYPoint> result = new ArrayList<>();
		for(XYPoint point: preresult)
		{
			if(point.getX()>startMz && point.getX()<endMz)
			{
				result.add(point);
			}
		}
		return  result;*/
	}

	public static List<XYPoint> getSpectrumMS1(String spath, String fileName, int scanNum,
											   Hashtable<String, IndexedFile> ht)
			throws IOException, CensusGeneralException, CensusIndexOutOfBoundException {

		if (!spath.endsWith(File.separator)) {
			spath = spath + File.separator;
		}

		String path = spath + fileName;
		String fileNameWithoutExtension = fileName.substring(0, fileName.lastIndexOf("."));
		String spectrumPath = spath + File.separator + fileName + ".index";

		BufferedReader br = null;
		RandomAccessFile rfile = null;
		br = new BufferedReader(new FileReader(spectrumPath));
		// Hashtable<String, IndexedFile> ht=
		// ChroGenerator.createIndexedFiles(spath, fileName);
		int prevScanfromMS1 = getMS1Scan(ht, path, scanNum);

		String eachLine;
		while ((eachLine = br.readLine()) != null) {
			if (eachLine.startsWith(String.valueOf(prevScanfromMS1)))
				break;
		}

		String[] arr = eachLine.split("\t");
		String nextLine = br.readLine();

		rfile = new RandomAccessFile(path, "r");

		long startPos = Long.parseLong(arr[1]);
		rfile.seek(startPos);
		int byteSize = -1;

		if (null == nextLine) {
			byteSize = (int) (rfile.length() - startPos);
		} else {
			String[] nextArr = nextLine.split("\t");
			long nextPos = Long.parseLong(nextArr[1]);
			byteSize = (int) (nextPos - startPos);
		}

		byte[] bytes = new byte[byteSize];
		rfile.readFully(bytes);

		int pos = 0;
		char ch = (char) bytes[pos];
		while ((ch = (char) bytes[pos]) == 'Z' || ch == 'S' || ch == 'I' || ch == 'D') {
			while (ch != '\n') {
				pos++;
				ch = (char) bytes[pos];
			}

			pos++;
		}
		int arrSize = 0;
		for (int j = pos; j < byteSize; j++) {
			if ('\n' == (char) bytes[j])
				arrSize++;

		}

		double[] massArr = new double[arrSize];
		double[] intArr = new double[arrSize];

		StringBuilder mass = new StringBuilder(10);
		StringBuilder intensity = new StringBuilder(10);
		intensity.append('0');

		boolean isMass = true;
		boolean isInt = true;
		int countInt = 0;
		int massIndex = 0;
		for (int j = pos; j < byteSize; j++) {
			ch = (char) bytes[j];

			switch (ch) {
				case '\r':
					break;

				case ' ':
					countInt++;
					isMass = false;
					isInt = countInt <= 1 ? true : false;
					break;

				case '\n':
					isMass = true;
					isInt = false;
					intArr[massIndex] = Long.parseLong(intensity.toString());
					massArr[massIndex++] = Double.parseDouble(mass.toString());
					mass.delete(0, mass.length());
					intensity.delete(0, intensity.length()).append('0');
					countInt = 0;
					break;

				case '.':
					isInt = false;

				default:
					if (isMass)
						mass.append(ch);
					else if (isInt)
						intensity.append(ch);

					break;
			}

		}

		List<XYPoint> massInt = new ArrayList<>();

		for (int i = 0; i < intArr.length; i++) {
			XYPoint xy = new XYPoint();
			xy.setX(massArr[i]);
			xy.setY(intArr[i]);
			massInt.add(xy);
		}

		rfile.close();
		br.close();

		return massInt;
	}
	

}

