package edu.scripps.pms.census.chroalign;

import edu.scripps.pms.census.util.io.SpectrumReader;
import java.io.*;
import java.util.*;
import java.text.DecimalFormat;
import edu.scripps.pms.util.spectrum.*;
import edu.scripps.pms.census.model.*;
import edu.scripps.pms.census.chroalign.Elements;
import edu.scripps.pms.census.util.BasePeakFinder;
import edu.scripps.pms.census.ChroProgressDialog;
import org.jdom.*;
import org.jdom.output.*;
import edu.scripps.pms.census.conf.Configuration;

import javax.swing.*;

public class chroalign {
        private double[][][]    dataarray;
        private int             bandConstraint;
	private Elements[][] diaArray;
        private int[][][]       pathArray;
        private double[][]       retArray;
        private boolean         useBasePeak;
        private int fileNumber=0;        
        private double[][][] binSpectra = null;        
        private double[][][] info = null;        
        private ChroProgressDialog progress;
        private int alignNumber = 0; //increase this number as alignment proceeds
        private String[]        targetMS1Files = null;
        private String          referenceMS1File;
        private double percentDone;
        private Object[] filenames;
        private int masterfile;
        private int[][] mappingArr;
        int[][] scanArr;
	private String workFolder;
        private Hashtable[] htArr;
        private Hashtable[] indexArr;
        private ChroProgressDialog chroProgress;
        private int BOUNCEFACTOR = 1;
        
        public void setConstraint (int constraint) {
            this.bandConstraint = constraint;
        }
        
        public int getConstraint () {
            return bandConstraint;
        }
        
        public String getWorkFolder() {
            return workFolder;
        }

        public void setWorkFolder(String work) {
            this.workFolder = work;
        }
        
        public void setAlignType (boolean bp) {
            this.useBasePeak = bp;
        }
        
        public void setMasterFile (int master) {
            this.masterfile = master;
        }
        
        public int getMasterFile() {
            return masterfile;
        }
        
        public boolean getAlignType () {
            return this.useBasePeak;
        }
        
        public String[] getTargetMS1Files() {
            return targetMS1Files;
        }

        public void setTargetMS1Files(String[] target) {
            this.targetMS1Files = target;
        }
        
        public String getReferenceMS1File() {
            return referenceMS1File;
        }
        
        public void setReferenceMS1File(String reference) {
            this.referenceMS1File = reference;
        }

        public static void main(String args[]) throws Exception {
                chroalign           appObject = new chroalign();
                ArrayList           tfiles = new ArrayList();
                ArrayList           tnames = new ArrayList();
                String[]            files = null;
                String[]            names = null;
                int[][][]           subArray;
                boolean             mode;
                int                 cnt = 0;
               
                if(args.length < 1) {
                        appObject.usage();
            		System.exit(0);
        	}
                
                //HEADER
                System.out.println("Chroalign v1.0.1");
                
                //READ IN COMMAND LINE ARGS
                for(int i = 0; i < args.length; i++) {
                    if ( (args[i].equals("--help")) ) {
                        appObject.usage();
                        System.exit(0);
                    }
                    else if (args[i].equals("-bc")) {
                        //SET BAND CONSTRAINT
                        appObject.setConstraint(new Integer(args[i+1]).intValue());
                        i++;
                    }
                    else if (args[i].equals("-mode")) {
                        //SET ALIGNMENT MODE
                        if (args[i+1].equals("true")) {
                            mode = true;
                        } else {
                            mode = false;
                        }
                        appObject.setAlignType(mode);
                        i++;
                    }
                    else if (args[i].equals("-F")) {
                        //SET FILE PATH
                        tfiles.add(args[i+1]);
                        //SET FILE NAME
                        tnames.add(args[i+2]);
                        cnt++;
                        i+=2;
                    }
                    else if (args[i].equals("-ma")) {
                        //SET MASTER FILE
                        appObject.setMasterFile(new Integer(args[i+1]).intValue());
                        i++;
                    }
                    else {
                        System.out.println("I don't understand this option " + args[i]);
                        appObject.usage();
                        System.exit(0);
                    }
                }
                appObject.setWorkFolder(System.getProperty("user.dir"));
                files = new String[tfiles.size()];
                names = new String[tfiles.size()];
                for(int j=0; j<tfiles.size(); j++) {
                    files[j] = tfiles.get(j).toString();
                }
                for(int j=0; j<tfiles.size(); j++) {
                    names[j] = tnames.get(j).toString();
                }
                subArray = appObject.alignChro(null, files, names, appObject.getMasterFile(), appObject.getAlignType(), appObject.getConstraint(), appObject.getWorkFolder());
        }

        public void usage() {
                System.out.println("\nA program for aligning chromatograms");
                System.out.println();
                System.out.println("Options:");
                System.out.println("-bc\tSet band constraint (maximum # scans shifted between chromatograms)");
                System.out.println("-mode\tSet alignment mode (true = base peak, false = spectral based correlation)");
                System.out.println("-F\tFilePath FileName");
                System.out.println("-ma\tIndex of master file (0 to number of files - 1)");
                System.out.println("--help #\tDisplays options");
		System.out.println("USAGE: chroalign -F 022707_BSAdig_SetA-01_ftms.ms1 reference -F 022707_BSAdig_SetC-01_ftms.ms1 target -bc 500 -mode false -ma 0");
	}
        
        public int[][][] noAlignChro(            
                ChroProgressDialog progress,
                Object[] args, 
                Object[] datasetnamesObj, 
                String workFolder                
                ) throws IOException 
        {        
            this.workFolder = workFolder; 
            this.progress = progress;
            this.chroProgress = chroProgress;
            this.filenames = args;
            targetMS1Files = new String[args.length];
            String[]        tempnames = new String[args.length];
            Element         chroRoot;
            Element         pathRoot;
            Element         alignedRoot;
            ArrayList<Element>  results = new ArrayList<Element>();
            SpectrumReader  ms1;
            int             cnt;
            int             refNum;
            int[]           tarNum = new int[args.length];
            this.fileNumber = args.length;
            referenceMS1File = args[0].toString();
            String[] datasetnames = new String[datasetnamesObj.length];
                
            try {
                    ArrayList<SpectrumReader>  data = new ArrayList(); 
                    SpectrumReader sReader;
                    
                    int size=0;

                    for(int i=0;i<args.length;i++) {
                        long start = System.currentTimeMillis();
                        targetMS1Files[i] = args[i].toString();
                        datasetnames[i] = datasetnamesObj[i].toString();
                        System.out.println("Reading Target MS1 File # " + (i+1) + " ...");
                        ChroProgressDialog.addMessageWithLine(progress, "Reading Target MS1 File # " + (i+1) + " ...");
                        System.out.println(targetMS1Files[i]);
                        ChroProgressDialog.addMessageWithLine(progress, targetMS1Files[i]);
                        sReader = new SpectrumReader(targetMS1Files[i], "ms1");
                        
                        tarNum[i] = sReader.getNumSpectra();
                        
                    //    System.out.println("size===" + size + " " + tarNum[i] + " " + args[i].toString());
                        
                        if(tarNum[i]>size)
                            size = tarNum[i];

                        ms1 = new SpectrumReader(targetMS1Files[i], "ms1");

                        data.add(ms1);
                     //   System.out.println("===>"+ (System.currentTimeMillis()-start));
                    }
                        
                    dataarray = new double[this.fileNumber][2][size];//stores bp chromatograms for output
                    pathArray = new int[this.fileNumber][2][2*size];
                    retArray = new double[this.fileNumber][size];
                    this.populateNoAlignDataArray(data, size);
                    
                    double maxRet =0;
                    
                    for(int i=0;i<retArray.length;i++)
                    {
                        double tmp = retArray[i][retArray[i].length-1];
                        
                        System.out.println("ret time== " + maxRet);
                        
                        int tmp1 = retArray[i].length-1;
                        
                        
                        while(true)
                        {
                            double tmpRet = retArray[i][tmp1--];
                            
                            if(tmp1<0)
                                break;
                            
                            if(tmpRet>0)
                            {
                                
                                if(tmpRet>maxRet)
                                    maxRet = tmpRet;
                               
                                break;
                            }
                            
                        }
                    }
                    
                    int[][] binArr = new int[this.fileNumber][(int)((maxRet+1)*100)];                    
                    pathArray = new int[this.fileNumber][2][binArr[0].length];
                    
                    
                    for(int i=0;i<retArray.length;i++)
                    {
                        for(int j=0;j<retArray[i].length;j++)
                        {
                            if(retArray[i][j]>0)
			    {
                                binArr[i][(int)(retArray[i][j]*100)] = j; //retArray[i][j]; //scan index;                       
			    }
                        }
                    }   
                    
                    for(int i=0;i<binArr.length;i++)
                    {
			int prev=-1;

                        for(int j=0;j<binArr[i].length;j++)
			{
			    if(binArr[i][j] == 0 && j>0)
			    {
				binArr[i][j] = binArr[i][j-1]; 
			    }
			}
                    }          
                    int tmpInd=0;
                   
		   int k=0; //last Index
                    for(int i=0;i<binArr[0].length;i++)
                    {
                        for(int j=0;j<binArr.length;j++)
                        {
                            if(0<binArr[j][i])
                            {
                                for(k=0;k<binArr.length;k++)
                                {                                    
                                    if(binArr[k][i]<=0 && tmpInd>0 && pathArray[k][1][tmpInd-1]>0) {
                                        pathArray[k][0][tmpInd] = pathArray[k][1][tmpInd-1];
                                        pathArray[k][1][tmpInd] = pathArray[k][1][tmpInd-1];   
                                    }
                                    else {
                                        pathArray[k][0][tmpInd] = binArr[k][i];
                                        pathArray[k][1][tmpInd] = binArr[k][i];
                                    }
                                }
                                
                                tmpInd++;       
                                
                                break;
                            }
                        }
                    }                    

		    System.out.println(pathArray.length +  " " + k + " " + tmpInd);
		    System.out.println(pathArray[0][0].length +  " " + k + " " + tmpInd);
                  
		    //fill up all zero values with last filled value 

		    for(int i=0;i<pathArray.length;i++)
		    {
			for(int j=tmpInd;j<pathArray[0][0].length;j++)
			{
			    pathArray[i][0][j] = pathArray[i][0][tmpInd-1];
			    pathArray[i][1][j] = pathArray[i][1][tmpInd-1];
			}
		    }


                    Configuration conf = Configuration.getInstance();
                    conf.setRetArr(retArray);
                    
                    System.out.println("Outputting results ...");
                    ChroProgressDialog.addMessageWithLine(progress, "Outputting results ...");
                    
                    //UNALLIGNED CHROMATOGRAMS (chro_out.xml)
                    chroRoot = createXmlChroHeader();
                    chroRoot = createChrodata(chroRoot, dataarray, datasetnames, data.size(), false);
                    writeOutputFile(chroRoot, "chro_out.xml");

                    //ALLIGNED CHROMATOGRAMS (aligned_out.xml)
                    alignedRoot = createXmlChroHeader();
                    alignedRoot = createChrodata(alignedRoot, dataarray, datasetnames, data.size(), true);
                    writeOutputFile(alignedRoot, "aligned_out.xml");
                    
                    //paths (path_out.xml)
                    pathRoot = createXmlPlotHeader();
                    pathRoot = createPathdata(pathRoot, dataarray, datasetnames, data.size());
                    writeOutputFile(pathRoot, "path_out.xml");
                    
            } catch(IOException failure) {
                    System.out.println("Error while reading MS1 files");
                    ChroProgressDialog.addMessageWithLine(progress, "Error while reading MS1 files");
                    System.out.println(failure.toString());
                    ChroProgressDialog.addMessageWithLine(progress, failure.toString());
            }
            return this.pathArray;
        
        }
        
	//**************MAIN ROUTINE IN PROGRAM - GENERATES LOCAL COST ARRAY AND OVERWRITES**********//
        //**************WITH CUMULATIVE COST ARRAY. USES TEMPORARY ARRAY TO STORE 2 COLUMNS**********//
        //**************WORTH OF LOCAL COST ARRAY BEFORE OVERWRITING. ALSO USES A ROTATED************//
        //**************ARRAY (DIAARRAY) TO SAVE MEMORY**********************************************//
        public int[][][] alignChro(ChroProgressDialog progress, Object[] filepaths, 
                Object[] datasetnames,int masterfile,boolean useBasePeak,int bandConstraint, 
                String workFolder) throws IOException, Exception {
                      
		this.workFolder = workFolder; 
                this.progress = progress;                    
                this.filenames = filepaths;
                this.masterfile = masterfile;
                this.bandConstraint = bandConstraint;
                this.fileNumber = filepaths.length;
                this.useBasePeak = useBasePeak;
		targetMS1Files = new String[filepaths.length];
                String[]        tempnames = new String[filepaths.length];
                Element         chroRoot;
                Element         pathRoot;
                Element         alignedRoot;
                ArrayList<Element>  results = new ArrayList<Element>();
                SpectrumReader  ms1;
                int             cnt;
                int             refNum;
                int[]           tarNum = new int[filepaths.length];
                
                //REASSIGNS PATHS OF MS1 FILES INTO ARRAY WITH THE MASTER FILE = REFERENCEMS1FILE AND
                //THE REMAINING FILES IN TARGETMS1FILES
                //ALSO PUTS DATASET NAMES INTO ARRAY CALLED DATASETNAMES
                referenceMS1File = filepaths[masterfile].toString();
                tempnames[0] = datasetnames[masterfile].toString();
                cnt = 0;
                for(int i=0;i<filepaths.length;i++) {
                    if(i != masterfile) {
                        targetMS1Files[cnt] = filepaths[i].toString();
                        cnt++;
                        tempnames[cnt] = datasetnames[i].toString();
                    }
                }
                datasetnames = tempnames;
                
                //READS IN MS1 FILES AND CREATES SPECTRUM READER OBJECTS - PUTS SR'S INTO
                //ARRAY LIST CALLED DATA
		try {
			ArrayList<SpectrumReader>  data = new ArrayList(); 
                        //READ REFERENCE FILE FIRST
                        System.out.println("Reading Reference MS1 ...");
                        System.out.println(referenceMS1File);

			if(null != progress)
			{
			    ChroProgressDialog.addMessageWithLine(progress, "Reading Reference MS1 ...");
			    ChroProgressDialog.addMessageWithLine(progress, referenceMS1File);
			}

                        SpectrumReader sReader = new SpectrumReader(referenceMS1File, "ms1");
                        //FIND # SPECTRA IN REFERENCE FILE
                        refNum = sReader.getNumSpectra();
                        ms1 = new SpectrumReader(referenceMS1File, "ms1");   
                        data.add(ms1);

                        //ITERATE THROUGH READING ALL TARGET FILES
                        int size=0;
                        for(int i=0;i<filepaths.length-1;i++) {
                            long start = System.currentTimeMillis();
                            System.out.println("Reading Target MS1 File # " + (i+1) + " ...");
                            ChroProgressDialog.addMessageWithLine(progress, "Reading Target MS1 File # " + (i+1) + " ...");
                            System.out.println(targetMS1Files[i]);
                            ChroProgressDialog.addMessageWithLine(progress, targetMS1Files[i]);
                            sReader = new SpectrumReader(targetMS1Files[i], "ms1");
                            //FIND # SPECTRA IN TARGET FILES
                            tarNum[i] = sReader.getNumSpectra();
                            //DETERMINE MAXIMUM SIZE OF TARGET FILES
                            if(tarNum[i]>size)
                                size = tarNum[i];
                            ms1 = new SpectrumReader(targetMS1Files[i], "ms1");
                            data.add(ms1);
                        }
                        
                        //DETERMINE WHETHER TARGET OR REFERENCE IS BIGGER AND SAVE BIGGEST IN SIZE
                        size = refNum>size?refNum:size;
                        
                        //ADD MARGIN TO SIZE
                        size += 2;
                        
                        //GENERATES ARRAYS TO STORE
                        dataarray = new double[this.fileNumber][2][size];   //BP CHROMATOGRAMS
                        pathArray = new int[this.fileNumber][2][2*size];    //ALIGNMENT PATH DATA
                        mappingArr = new int[this.fileNumber][size];        //
                        
                        //ITERATES THROUGH SPECTRUM READERS AND POPULATES DATAARRAY AND INFO ARRAY
                        this.populateDataArray(data, size);
                        
                        //************************************************************//
                        //******PROCESS EACH PAIR OF CHROMATOGRAMS SEQUENTIALLY*******//
                        //***********************************************************//
			
                        //INITIALIZE DIAGONAL ARRAY WHICH IS A ROTATED ARRAY TO SAVE MEMORY
			int height = bandConstraint*2;
			diaArray = new Elements[size][height];                        
			for(int i=0;i<size;i++) 
			{
			    for(int j=0;j<height;j++)
				diaArray[i][j] = new Elements();
			}
                        
                        //UPDATE PROGRESS BAR
                        int progressValue = (int)(10.0/data.size()) + 5;
                        ChroProgressDialog.updateProgress(progress,  progressValue + 1 );
                        
                        //PROCESS EACH CHROMATOGRAM PAIR SEQUENTIALLY
                        for (int i=1;i<data.size();i++) {
                            
                            //UPDATE PROGRESS BAR
                            System.out.println("Working on matrix #" + i);
                            ChroProgressDialog.addMessageWithLine(progress, "Working on matrix #" + i);
                            
                            //PRINT OUT ALIGNMENT MODE
                            if(useBasePeak) {                                
                                ChroProgressDialog.addMessageWithLine(progress, "Alignment Mode = BasePeak");
                                System.out.println("Alignment Mode = BasePeak");
                            } else {
                                ChroProgressDialog.addMessageWithLine(progress, "Alignment Mode = Pearson Correlation @ Spectral Level");                                        
                                System.out.println("Alignment Mode = Pearson Correlation @ Spectral Level");
                            }
                            
                            //GENERATES LOCAL COST MATRIX AND THEN OVERWRITES WITH CUMULATIVE COST MATRIX
                            alignNumber = i;
                            genCostMatrix(data.get(0), data.get(i), i, size-2, size-2, useBasePeak);
                            
                            //FINDS PATH THROUGH CUMULATIVE COST MATRIX USING DYNAMIC PROGRAMMING
                            findPath(i, size-2, size-2, this.bandConstraint);
                            
                            //PRINT OUT FINISHED MESSAGE
                            System.out.println("Finished with this matrix ...");
                            ChroProgressDialog.addMessageWithLine(progress, "Finished with this matrix ...");
                         
                            //UPDATE PROGRESS BAR
                            progressValue += (int)10.0/data.size();
                            ChroProgressDialog.updateProgress(progress,  progressValue );
                        }
                        
                        //OUTPUTS FILES
			System.out.println("Outputting results ....");
                        ChroProgressDialog.addMessageWithLine(progress, "Outputting results ....");
                        
                        //UNALLIGNED CHROMATOGRAMS (chro_out.xml)
                        chroRoot = createXmlChroHeader();
                        chroRoot = createChrodata(chroRoot, dataarray, datasetnames, data.size(), false);
			writeOutputFile(chroRoot, "chro_out.xml");
			System.out.println("chro_out.xml was done");
                        ChroProgressDialog.addMessageWithLine(progress, "chro_out.xml was done");
                        
                        //ALLIGNED CHROMATOGRAMS (aligned_out.xml)
                        alignedRoot = createXmlChroHeader();
                        alignedRoot = createChrodata(alignedRoot, dataarray, datasetnames, data.size(), true);
			writeOutputFile(alignedRoot, "aligned_out.xml");
			System.out.println("aligned_out.xml was done");
                        ChroProgressDialog.addMessageWithLine(progress, "aligned_out.xml was done");
                        
                        //ALIGNMENT PATH DATA (path_out.xml) - THIS IS THE MOST IMPORTANT FILE THAT GETS READ
                        //LATER BY CENSUS TO DETERMINE RELATIONSHIP BETWEEN FILES 
                        pathRoot = createXmlPlotHeader();
                        pathRoot = createPathdata(pathRoot, dataarray, datasetnames, data.size());
			writeOutputFile(pathRoot, "path_out.xml");
			System.out.println("path_out.xml was done");
                        ChroProgressDialog.addMessageWithLine(progress, "path_out.xml was done");
                        
		} catch(IOException failure) {
                        System.out.println("Error while reading MS1 files");
                        ChroProgressDialog.addMessageWithLine(progress, "Error while reading MS1 files");
                        ChroProgressDialog.addMessageWithLine(progress, failure.toString());
                        System.out.println(failure.toString());
                        throw new IOException();
		} catch(Exception e) {
                    if(null != progress)
                        ChroProgressDialog.addMessageWithLine(progress, "Error : " + e);
                    
                    System.out.println("Error : " + e);
		    e.printStackTrace();
                    throw new Exception();
                }
                return this.pathArray;
        }
        
        public void populateNoAlignDataArray(ArrayList<SpectrumReader> readerList, int size) throws IOException
        {
            System.out.println("Generating Chromatogram Arrays ...");

	    if(null != progress)
		ChroProgressDialog.addMessageWithLine(progress, "Generating Chromatogram Arrays ...");
            
            int listSize = readerList.size();
            info = new double[listSize][size][4];
            scanArr = new int[listSize][size];
            int readerIndex=0;
            htArr = new Hashtable[readerList.size()];
            
            for(int i=0;i<htArr.length;i++)
            {
                htArr[i] = new Hashtable<Integer, Double>();
            }
            
            for(Iterator<SpectrumReader> itr=readerList.iterator(); itr.hasNext(); )
            {
                SpectrumReader reader = itr.next();
                
                int cnt=0;    
                for(Iterator<PeakList> peakItr=reader.getSpectra(); peakItr.hasNext(); )
                {
                    PeakList pList = peakItr.next();
                    scanArr[readerIndex][cnt] = pList.getLoscan();
                    
                    List<String> rtList = pList.getIlines();
                    
                    double rt=-1;
                    
                    for(Iterator<String> rtItr=rtList.iterator(); rtItr.hasNext(); )
                    {
                        String eachI = rtItr.next();
                        
                        if(eachI.startsWith("I\tRetTime"))
                        {
                            String[] arr = eachI.split("\t");
                            String retTime = arr[2];
                        
                            rt  = Double.parseDouble(retTime);
                            this.retArray[readerIndex][cnt] = rt;
                            
                            break;
                        }
                    }
                    htArr[readerIndex].put(cnt, rt);
                    
                    cnt++;                    
                }
             
                readerIndex++;
            }

            
        }
        
        public void populateDataArray(ArrayList<SpectrumReader> readerList, int size) throws IOException
        {
            
            System.out.println("Generating Chromatogram Arrays ...");

	    if(null != progress)
		ChroProgressDialog.addMessageWithLine(progress, "Generating Chromatogram Arrays ...");
            
            int listSize = readerList.size();
            binSpectra = new double[listSize][size][2];
            info = new double[listSize][size][4];
            scanArr = new int[listSize][size];

            int readerIndex=0;
            
            int progressValue=1;
                        
            for(Iterator<SpectrumReader> itr=readerList.iterator(); itr.hasNext(); )
            {
                SpectrumReader reader = itr.next();
                
                int cnt=0;    
                for(Iterator<PeakList> peakItr=reader.getSpectra(); peakItr.hasNext(); )
                {
                    PeakList pList = peakItr.next();
                    scanArr[readerIndex][cnt] = pList.getLoscan();
                    
                    List<String> rtList = pList.getIlines();
                    String[] content = rtList.get(0).split("\t");
                    double rt = new Double(content[2]).doubleValue();
                    Iterator peaks = pList.getPeaks();
                    
                    //get info from spectra
                    binSpectra[readerIndex][cnt] = binSpectrum(peaks, 1, pList.getMaxM2z(), 2000, pList.numPeaks());
                    info[readerIndex][cnt] = getSpectrumInfo(binSpectra[readerIndex][cnt]);
                    dataarray[readerIndex][0][cnt] = rt;
                    dataarray[readerIndex][1][cnt] = info[readerIndex][cnt][0];//base peak
                    cnt++;                    
                }
             
                readerIndex++;
                
                progressValue += (int)(5.0/readerList.size());
                ChroProgressDialog.updateProgress(progress,  progressValue + 1 );
                
                ChroProgressDialog.addMessage(progress, ".");
                        
            }
            
        }
        
        public void genCostMatrix (SpectrumReader refms1, SpectrumReader tarms1, int filenum, int refNum, int tarNum, boolean useBasePeak) {

            BasePeakFinder  bpf;
            double          temp;
            double  dmin = 0;
            int     index = 0;
            double  tt;
            double	dd = 0;
            double  dd1 = 0;
            double	dd2 = 0;
            double	dd3 = 0;
            String tstr;
            String  TString = "0.###";
            DecimalFormat    formatter = new DecimalFormat(TString);
            double[][]	tempLocalCostArray = new double[2][tarNum+3];
	    System.out.println("Calculating Local Distances ...");
            ChroProgressDialog.addMessageWithLine(progress, "Calculating Local Distances ...");
            
            //initialize temporary array and first two columns of matrix array with local score
            Elements[][] accCostArray = new Elements[3][tarNum+3];
            for(int i=0;i<3;i++)
                for(int j=0;j<tarNum+3;j++)
                    accCostArray[i][j] = new Elements();
                                
	    for (int i=2;i<3;i++) {
		for (int j=2;j<tarNum+2;j++) {
                    if (AlignNode.checkWithinBound(this.bandConstraint, i, j))
		    {
			if (useBasePeak) {
			    tempLocalCostArray[0][j] = euclidianDis(dataarray[0][1][i-2], dataarray[filenum][1][j-2]);
                            tempLocalCostArray[1][j] = euclidianDis(dataarray[0][1][i-1], dataarray[filenum][1][j-2]);
                        } else {
			    tempLocalCostArray[0][j] = spectrumDis(binSpectra[0][i-2], binSpectra[filenum][j-2], info[0][i-2], info[filenum][j-2], false); 
                            tempLocalCostArray[1][j] = spectrumDis(binSpectra[0][i-1], binSpectra[filenum][j-2], info[0][i-1], info[filenum][j-2], false);
                        }
		    }
                    accCostArray[1][j].score = tempLocalCostArray[0][j];
                    accCostArray[2][j].score = tempLocalCostArray[1][j];
		}
	    }
            
            double[][]    tarr = new double[3][tarNum+1];
	    for (int i=2;i<refNum;i++) {
		percentDone = (double)i/(double)refNum/(this.fileNumber-1) + (double)(this.alignNumber-1)/(this.fileNumber-1);
                //new
                if (i==2) {
                    //do nothing
                } else {
                    //shift tempLocalCostArray and calc local cost
                    for (int j=2;j<tarNum+3;j++) {
                            tempLocalCostArray[0][j] = tempLocalCostArray[1][j];
                            //calculate local cost and assign to tempLocalCostArray
                            if (useBasePeak) {
                                tempLocalCostArray[1][j] = euclidianDis(dataarray[0][1][i-1], dataarray[filenum][1][j-2]);
                                accCostArray[2][j].score = tempLocalCostArray[1][j];
                            } else {
                                tempLocalCostArray[1][j] = spectrumDis(binSpectra[0][i-1], binSpectra[filenum][j-2], info[0][i-1], info[filenum][j-2], false);
                                accCostArray[2][j].score = tempLocalCostArray[1][j];
                            }
                    }
                }	
		for (int j=2;j<tarNum+2;j++) {
		    //then for each row calculate scores
                    dd = tempLocalCostArray[1][j+1];
                    dd1 = accCostArray[1][j].score+2*dd;
		    dd2 = accCostArray[0][j].score + 2*tempLocalCostArray[0][j+1] + dd;
		    dd3 = accCostArray[1][j-1].score + 2*tempLocalCostArray[1][j] + dd;
		    if (useBasePeak)
		    {
			if (dd2 < dd1 && dd2 < dd3) {
			    index = 2;
			    dmin = dd2;
			} else if (dd3 < dd1 && dd3 < dd2) {
			    index = 3;
			    dmin = dd3;
			} else {
			    index = 1;
			    dmin = dd1;
			}
		    }                                
		    else
		    {
			if (dd2 > dd1 && dd2 > dd3) {
			    index = 2;
			    dmin = dd2;
			} else if (dd3 > dd1 && dd3 > dd2) {
			    index = 3;
			    dmin = dd3;
			} else {                              
			    index = 1;
			    dmin = dd1;
			}
		    }
                    
                    if(AlignNode.checkWithinBound(this.bandConstraint, i-1, j))
                    {
                        int[] newCor = AlignNode.rotateOrigin(this.bandConstraint, i-1, j);
                        this.diaArray[newCor[0]][newCor[1]].score = dmin;
                        
                    }
                    if(AlignNode.checkWithinBound(this.bandConstraint, i-2, j-1))
                    {
                        int[] newCor = AlignNode.rotateOrigin(this.bandConstraint, i-2, j-1);
                        this.diaArray[newCor[0]][newCor[1]].index = index;                                               
                    }
                    accCostArray[2][j+1].score = dmin;
		}
                for (int j=2;j<tarNum+2;j++) {
                    accCostArray[0][j].score = accCostArray[1][j].score;
                    accCostArray[1][j].score = accCostArray[2][j].score;
                }
	    }
          
        }

        
        public boolean checkConstraint(int refNum, int tarNum, int i, int j) {
            
            boolean b = false;
            double  tt;
            
            tt = (double)((double)tarNum/(double)j);
            tt = (double)(refNum / tt);
            tt = i-tt;
            if (Math.abs(tt) > this.bandConstraint) {
                b = true;
            }
            return b;
        }
        
        public void findPath(int numfile, int refNum, int tarNum, int bconstraint) {
		int	index = 0;
		int	cnt = 0;
                int     ind;
		double	temp;
                boolean exportMatrix = false;
		
		//start at top right hand corner of matrix and find path down to bottom left corner
                System.out.println("Finding Path Through Matrix ...");
                ChroProgressDialog.addMessageWithLine(progress, "Finding Path Through Matrix ...");
                int i = refNum;
                int j = tarNum;
                ind = numfile - 1;
		while (i > 2 || j > 2 ) {
                        if (AlignNode.checkWithinBound(bconstraint, i, j)) {
                            int[] newCor = AlignNode.rotateOrigin(bconstraint, i, j);
                            index = this.diaArray[newCor[0]][newCor[1]].index;
                            temp = 0;
                            if (index == 1) {
                                    i = i-1;
                                    j = j-1;
                                    pathArray[ind][0][cnt] = (i<2)?2:i;
                                    pathArray[ind][1][cnt] = (j<2)?2:j;

                            } else if (index == 2) {
                                    i = i-2;
                                    j = j-1;
                                    pathArray[ind][0][cnt] = (i<2)?2:i+1;
                                    pathArray[ind][1][cnt] = (j<2)?2:j;
                                    cnt++;
                                    pathArray[ind][0][cnt] = (i<2)?2:i;
                                    pathArray[ind][1][cnt] = (j<2)?2:j;
                            } else if (index == 3) {
                                    i = i-1;
                                    j = j-2;
                                    pathArray[ind][0][cnt] = (i<2)?2:i;
                                    pathArray[ind][1][cnt] = (j<2)?2:j+1;
                                    cnt++;
                                    pathArray[ind][0][cnt] = (i<2)?2:i;
                                    pathArray[ind][1][cnt] = (j<2)?2:j;
                            } else {    
                                    System.out.println("error finding path!");
                                    ChroProgressDialog.addMessageWithLine(progress, "error finding path!");
                                    System.exit(0);
                            }
                            cnt++;
                        } else {
                            //BOUNCE OFF BY MOVING 2 INDEXES IN APPROPRIATE DIRECTION
                            //MAY NEED TO THINK ABOUT A BETTER WAY TO DO THIS
                            if (i<j) {
                                i = i - 1;
                                j = j - 1 - this.BOUNCEFACTOR;
                                pathArray[ind][0][cnt] = (i<2)?2:i;
                                pathArray[ind][1][cnt] = (j<2)?2:j+1;
                                cnt++;
                                pathArray[ind][0][cnt] = (i<2)?2:i;
                                pathArray[ind][1][cnt] = (j<2)?2:j;
                            } else {
                                i = i - 1 - this.BOUNCEFACTOR;
                                j = j - 1;
                                pathArray[ind][0][cnt] = (i<2)?2:i+1;
                                pathArray[ind][1][cnt] = (j<2)?2:j;
                                cnt++;
                                pathArray[ind][0][cnt] = (i<2)?2:i;
                                pathArray[ind][1][cnt] = (j<2)?2:j;
                            }
                        }
		}	
	}

        private static Element createXmlChroHeader()
    {
            
        Element root = new Element("plot");
        
        Element title = new Element("title");
        title.setText("Chromatogram\n");
        root.addContent(title);

        Element xLabel = new Element("xLabel");
        xLabel.setText("\nRetention Time");
        root.addContent(xLabel);

        Element yLabel = new Element("yLabel");
        yLabel.setText("\nIntensity");
        root.addContent(yLabel);

        Element grid = new Element("noGrid");
        root.addContent(grid);

        Element size = new Element("size");
        size.setAttribute("width", "1000");
        size.setAttribute("height", "300");
        root.addContent(size);

        return root;
    }
        
        private static Element createXmlPlotHeader()
    {
            
        Element root = new Element("plot");
        
        Element title = new Element("title");
        title.setText("Paths\n");
        root.addContent(title);

        Element xLabel = new Element("xLabel");
        xLabel.setText("\nReference Chromatogram");
        root.addContent(xLabel);

        Element yLabel = new Element("yLabel");
        yLabel.setText("\nTarget Chromatogram");
        root.addContent(yLabel);

        Element grid = new Element("noGrid");
        root.addContent(grid);

        Element size = new Element("size");
        size.setAttribute("width", "750");
        size.setAttribute("height", "750");
        root.addContent(size);

        return root;
    }
        //first file in array is master
        public Element createChrodata(Element root, double[][][] files, Object[] datasetnames, int numfiles, boolean aligned) {
            int             j;
            int             i;
	    int             k;
	    int             tickLength;
            String          TString = "0.##";
            DecimalFormat   formatter = new DecimalFormat(TString);
            double[][]      arr;
            int[][]         alignedarr;
            int             index = 0;
            
            //master file is first
            arr = files[0];
            tickLength = (int)(arr[1].length / 10);
            Element xTicks = new Element("" +
                    "xTicks");
            root.addContent(xTicks);
            k=0;
            for(j=0;j<arr[1].length;j++) {
                if(arr[1][j] > 0) {
                    if (k == tickLength) {
                        Element tick = new Element("tick");
                        tick.setAttribute("label", formatter.format(arr[0][j]));
                        tick.setAttribute("position",formatter.format(arr[0][j]));
                        root.addContent(tick);
			k = 0;
                    }
                }
                k++;
            }   
            while(index < numfiles) {
                
                Element dataset = new Element("dataset");
                dataset.setAttribute("name", filenames[index].toString());                
                dataset.setAttribute("marks", "none");
                dataset.setAttribute("connected", "yes");
                dataset.setAttribute("stems", "no");
                k=0;
                if (aligned) {
                    alignedarr = pathArray[(int)index/2];
                    for(j=0;j<alignedarr[1].length;j++) {
                        if(alignedarr[0][j] > 0) {
                            Element point = new Element("p");
                            point.setAttribute("x", formatter.format(files[0][0][alignedarr[0][j]]));
                            if (index == 0) {
                                point.setAttribute("y", formatter.format(arr[1][alignedarr[0][j]]));
                            } else {
                                point.setAttribute("y", formatter.format(arr[1][alignedarr[1][j]]));
                            }
                            dataset.addContent(point);
                        }
                    }
                } else {
                    for(j=0;j<arr[1].length;j++) {
                        if(arr[1][j] > 0) {
                            Element point = new Element("p");
                            point.setAttribute("x", formatter.format(arr[0][j]));
                            point.setAttribute("y", formatter.format(arr[1][j]));
                            dataset.addContent(point);
                        }
                    }
                    k++;
                }
                root.addContent(dataset);            
                index++;
                if (index < numfiles) {
                    arr = files[index];
                }
            }
            return root;
        }
        
        public Element createPathdata(Element root, double[][][] files, Object[] datasetnames, int numfiles) {
            int             j;
            int             i;
	    int             k;
	    int             tickLength;
            String          TString = "0.##";
            DecimalFormat   formatter = new DecimalFormat(TString);
            int[][]      arr;
            int             index = 0;
            
            //master file is first
            arr = pathArray[0];
            tickLength = (int)(arr[1].length / 10);
            Element xTicks = new Element("xTicks");
            root.addContent(xTicks);
            k=0;
            for(j=0;j<arr[1].length;j++) {
                if(arr[0][j] > 0) {
                    if (k == tickLength) {
                        Element tick = new Element("tick");
                        tick.setAttribute("label", formatter.format(arr[0][j]));
                        tick.setAttribute("position",formatter.format(arr[0][j]));
                        root.addContent(tick);
			k = 0;
                    }
                }
                k++;
            }
            Element yTicks = new Element("yTicks");
            root.addContent(yTicks);
            k=0;
            for(j=0;j<arr[1].length;j++) {
                if(arr[0][j] > 0) {
                    if (k == tickLength) {
                        Element tick = new Element("tick");
                        tick.setAttribute("label", formatter.format(arr[1][j]));
                        tick.setAttribute("position",formatter.format(arr[1][j]));
                        root.addContent(tick);
			k = 0;
                    }
                }
                k++;
            }
            
            while(index < numfiles-1) {
                
                Element dataset = new Element("dataset");
                if(index>=this.masterfile)
                    dataset.setAttribute("name", filenames[index+1].toString());
                else
                    dataset.setAttribute("name", filenames[index].toString());

                dataset.setAttribute("sample", datasetnames[index+1].toString());
                dataset.setAttribute("marks", "none");
                dataset.setAttribute("connected", "yes");
                dataset.setAttribute("stems", "no");
                
                k=0;
                for(j=0;j<arr[1].length;j++) {
                    if(arr[0][j] > 0) {
                        Element point = new Element("p");
                        point.setAttribute("x", formatter.format(arr[0][j]));
                        point.setAttribute("y", formatter.format(arr[1][j]));
                        dataset.addContent(point);
                    }
                    k++;
                }
                
                root.addContent(dataset);
                
                index++;
                arr = pathArray[index];
            }
            return root;
        }
        
        public void writeOutputFile(Element root, String filename) {
            //output xml file
            try {
                Document doc = new Document(root);
		
                OutputStream os = null;
		if( this.workFolder.endsWith("/") || this.workFolder.endsWith("\\") )
		    os = new FileOutputStream(this.workFolder + filename);
		else
		    os = new FileOutputStream(this.workFolder + File.separator + filename);

                XMLOutputter outputter = new XMLOutputter();
                outputter.setFormat(Format.getPrettyFormat());
                outputter.output(doc, os);
                os.close();
            } catch (IOException e)
            {
                System.out.println("IO Error while generating " + filename + " file : " + e);
                e.printStackTrace();
            }
        }
        
	public void writeOutputLine(String Filename, StringBuffer Data, boolean OverWrite) {
                try {
                        File             OutputFile = new File(Filename);
                        FileWriter       OutputFileWriter = new FileWriter(OutputFile, OverWrite);
                        BufferedWriter   Outgoing = new BufferedWriter(OutputFileWriter);
                        Outgoing.write(Data.toString());
                        Outgoing.close();
                } catch (IOException failure) {
                        System.out.println("IO Error while writing " + Filename);
                        System.exit(0);
                }
        }	

	public Peak basePeak(Iterator<Peak> peaks) {
		Peak	tpeak;
		Peak	bpeak = new Peak(0,0);
		double	maxint = 0.0;

		while (peaks.hasNext()) {
			tpeak = peaks.next();
			if (tpeak.getIntensity() > maxint) {
				maxint = tpeak.getIntensity();
				bpeak  = tpeak;
			}
		}
		return bpeak;
	}
        
        public double[] getSpectrumInfo(double[] peaks) {
		double	maxint = 0;
                double[]    info = new double[4];
                double      sum = 0;
                double      avg = 0;
                int         cnt = 0;
                double      sumOfDiff = 0;
                double      stdev = 0;

		//calculate average intensity and sum of intensities and base peak in binned spectrum
                for (int i=0;i<peaks.length;i++) {
                    sum += peaks[i];
                    if (peaks[i] > maxint) {
                            maxint = (double)peaks[i];
                    }
                    cnt++;
		}
                
                avg = sum / peaks.length;
                
                //calculate standard deviation of peak intensities for binned spectrum
                for (int i=0;i<peaks.length;i++) {
                    sumOfDiff += (peaks[i] - avg) * (peaks[i] - avg);
		}
                stdev = sumOfDiff / (peaks.length - 1);
                stdev = (double)Math.sqrt((double)stdev);
                
                //output info as 1D array 
                info[0] = maxint; //base peak intensity
                info[1] = sum;  //sum of intensities
                info[2] = avg;  //average intensity
                info[3] = stdev;  //stdev of intensity
                
		return info;
        }

	//pearson correlation between two binned spectra
	public double spectrumDis(double[] binRefSpectra, double[] binTarSpectra, double[] refinfo, double[] tarinfo,  boolean DotProduct) {
		int	i;
                double  sumOfDiff=0;
                double  answer = 1;
                
                if(binRefSpectra.length == binTarSpectra.length) {
                    for(i=0;i<binRefSpectra.length;i++) {
                        sumOfDiff += (binRefSpectra[i] - refinfo[2]) * (binTarSpectra[i] - tarinfo[2]);
                    } 
                    answer = sumOfDiff / (binRefSpectra.length * refinfo[3] * tarinfo[3]);
                }
                return answer;
        }
        
        double[] binSpectrum(Iterator<Peak> peaks, double precision, double max, int maxSize, int numpeaks) {
		double[]    arr;
                Peak        p;
                int        pcnt = 0;
                if (this.useBasePeak)
                    arr = new double[numpeaks+1];
                else    
                    arr = new double[(int)(maxSize/precision) + 1];
                while(peaks.hasNext()) {
                    p = peaks.next();
                    if (this.useBasePeak) {
                        arr[pcnt] = p.getIntensity();
                    } else {
                        if((int)(p.getM2z()/precision) < (int)(maxSize/precision)) {
                            arr[(int)(p.getM2z()/precision)] += p.getIntensity();
                        } else {
                            System.out.println("Peak " + p.getM2z() + " is bigger than " + maxSize);
                        }
                    }
                    pcnt++;
                }
		return arr;
	}
        
	public double euclidianDis(double q, double c) {
		double d;
                d = (q - c);
                d = d*d;
		return d;
	}
}
