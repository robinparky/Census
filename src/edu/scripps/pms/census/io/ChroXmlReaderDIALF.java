
/*
 * ChroReader.java
 *
 * Created on May 17, 2005, 12:00 PM
 */

package edu.scripps.pms.census.io;

import java.io.*;
import java.util.*;

import edu.scripps.pms.census.model.*;
import edu.scripps.pms.census.*;
import edu.scripps.pms.census.conf.Configuration;
import edu.scripps.pms.census.util.IsoData;

import edu.scripps.pms.census.util.LinearRegression;
import org.jdom.*;
import org.jdom.output.*;
import org.jdom.input.*;
import rpark.statistics.GaussianFitting;
import rpark.statistics.model.GaussianPeakModel;

import java.net.URI;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * @author Robin Park
 * @version $Id: ChroXmlReader.java,v 1.43 2014/08/06 05:29:53 rpark Exp $
 */

public class ChroXmlReaderDIALF {

  private BufferedReader br;
  private String lastLine = "";
  private ArrayList list = new ArrayList();
  private ChroData data = null;
  private URI uri = null;
  private boolean isDataDependent;

  private SAXBuilder builder = new SAXBuilder();
  private String fileName;
  private File file = null;
  private Element rootEle;
  private Configuration conf;
  private int quantLevel = 1; //default quant is from MS
  private boolean labeled = true;
  private String quantType = null;
  private int expType = -1;

  private ArrayList<String> fileList = new ArrayList<String>();
  private ArrayList<String> sampleList = new ArrayList<String>(); //sample list
  private Hashtable<String, String> fileSampleHt = new Hashtable<String, String>();  //filename(key) sample(value)
  private Hashtable<String, Sample> sampleObjList = new Hashtable<String, Sample>(); //sample array
  private ArrayList<Sample> sampleExpList = new ArrayList<Sample>();

  public Sample getSampleObj(String sampleName) {
    return this.sampleObjList.get(sampleName);
  }

  public void addSampleObj(Sample sample) {
    this.sampleObjList.put(sample.getSampleName(), sample);
  }

  public static class Sample {
    private ArrayList intensityArr = new ArrayList();
    private ArrayList expList = new ArrayList();

    private String sampleName;

    private ArrayList<Double> normalizedSpecC = new ArrayList<Double>();

    public void addNormalizedSpecC(double d) {
      this.normalizedSpecC.add(d);
    }

    //average of normalized values
    public double getNormalizedSpecC() {
      double sum = 0;

      for (Iterator<Double> itr = normalizedSpecC.iterator(); itr.hasNext(); ) {
        double each = itr.next().doubleValue();
        sum += each;
      }

      return sum / normalizedSpecC.size();
    }

    public void addExperiment(String expName) {
      this.expList.add(expName);
    }

    public ArrayList getExpList() {
      return this.expList;
    }

    public Sample(String sampleName) {
      this.sampleName = sampleName;
    }


    public void addIntensity(long intensity) {
      this.intensityArr.add(intensity);
    }

    public long getSumIntensity() {
      long total = 0;

      for (Iterator<Long> itr = intensityArr.iterator(); itr.hasNext(); ) {
        total += itr.next();
      }

      return total;
    }

    public long getAverage() {
      if (this.intensityArr.size() <= 0)
        return 0;

      return this.getSumIntensity() / this.intensityArr.size();
    }

    public String getSampleName() {
      return sampleName;
    }

    public void setSampleName(String sampleName) {
      this.sampleName = sampleName;
    }


    public ArrayList getIntensityArr() {
      return intensityArr;
    }

    public void setIntensityArr(ArrayList intensityArr) {
      this.intensityArr = intensityArr;
    }
  }

  public ChroXmlReaderDIALF(String fileName) throws IOException, JDOMException, Exception {
    this(new File(fileName), true);
  }

  public ChroXmlReaderDIALF(String fileName, boolean doInit) throws IOException, JDOMException, Exception {
    this(new File(fileName), doInit);
  }

  public ChroXmlReaderDIALF(File file) throws IOException, JDOMException, Exception {
    this(file, true);
  }

  public ChroXmlReaderDIALF(File file, boolean doInit) throws IOException, JDOMException, Exception {
    this.file = file;

    if (doInit)
      init();
  }

  public ChroXmlReaderDIALF(java.net.URI uri) throws IOException, JDOMException, Exception {
    this.uri = uri;
    this.file = null;

    init();
  }

  private void init() throws IOException, JDOMException, Exception {
    Document doc = null;

    if (null != this.uri) {
      InputStream in = uri.toURL().openStream();
      doc = builder.build(in);
    } else {
      doc = builder.build(file);
    }

    rootEle = doc.getRootElement();
    String expTypeStr = rootEle.getChildText("experiment_type");


    if (null == conf)
      conf = Configuration.getInstance();


    if (null != expTypeStr)
      setExpType(Integer.parseInt(expTypeStr));

    String labelType = rootEle.getChildText("label_type");

    this.labeled = false;

    conf.setLabeling(this.labeled);
    String quantLevelString = rootEle.getChildText("quantLevel");
    quantType = rootEle.getChildText("quantType");

    if (null != quantLevelString)
      this.quantLevel = Integer.parseInt(quantLevelString);

    conf.setQuantLevel(this.quantLevel);
    Element each = null;
    Element eachFile = null;
    Sample s = null;

      /*
        for(Iterator<Element> itr=rootEle.getChildren("sample").iterator(); itr.hasNext(); )
        {
            each = itr.next();

            this.sampleList.add( each.getChildText("name") );
            s = new Sample(each.getChildText("name"));

            for(Iterator<Element> fitr=each.getChild("ms_files").getChildren("file").iterator(); fitr.hasNext(); )
            {
                eachFile = fitr.next();
                fileList.add(eachFile.getText());
                s.addExperiment(eachFile.getText());

                fileSampleHt.put(eachFile.getText(), each.getChildText("name"));

            }

            this.sampleExpList.add(s);
        }
        */

    this.getProteinList(); //initialize all parameters, especially for stand-alone report generation
    //Move line to parameters.  We assume parameters start after carriage return
  }

  public ArrayList<ChroProtein> getMrmCrvProteinList() throws IOException, JDOMException {
    Document doc = builder.build(file);
    rootEle = doc.getRootElement();

    ArrayList<ChroProtein> list = new ArrayList<ChroProtein>();

    ChroProtein protein;
    ChroPeptide peptide;
    ChroData data;

    List proList = rootEle.getChildren("protein");
    for (Iterator<Element> itrPro = proList.iterator(); itrPro.hasNext(); ) {
      Element eachProtein = itrPro.next();

      protein = new ChroProtein();
      protein.setLocus(eachProtein.getAttributeValue("locus"));
      protein.setSeqCount(eachProtein.getAttributeValue("seq_ct"));
      protein.setSpectrumCount(eachProtein.getAttributeValue("spec_ct"));
      protein.setSeqCoverage(eachProtein.getAttributeValue("seq_cov"));
      protein.setLength(eachProtein.getAttributeValue("length"));
      protein.setMolWt(eachProtein.getAttributeValue("molwt"));
      protein.setPI(eachProtein.getAttributeValue("pi"));
      protein.setValidation(eachProtein.getAttributeValue("val"));
      protein.setDescription(eachProtein.getAttributeValue("desc"));

      List pepList = eachProtein.getChildren("peptide");

      for (Iterator<Element> itrPep = pepList.iterator(); itrPep.hasNext(); ) {
        Element eachPeptide = itrPep.next();

        peptide = new ChroPeptide();
        peptide.setUnique(!"".equals(eachPeptide.getAttributeValue("unique")));
        peptide.setFileName(eachPeptide.getAttributeValue("file"));
        peptide.setScanNum(Integer.parseInt(eachPeptide.getAttributeValue("scan")));
        peptide.setSequence(eachPeptide.getAttributeValue("seq"));
        peptide.setXCorr(eachPeptide.getAttributeValue("xcorr"));
        peptide.setChargeState(eachPeptide.getAttributeValue("charge"));
        peptide.setDeltCN(eachPeptide.getAttributeValue("deltaCN"));
        peptide.setMhPlus(eachPeptide.getAttributeValue("MHPlus"));
        peptide.setDeltMass(eachPeptide.getAttributeValue("deltaMass"));
        peptide.setSpRank(eachPeptide.getAttributeValue("spRank"));
        peptide.setSpScore(eachPeptide.getAttributeValue("spScore"));

        peptide.setDtaStartRange(Integer.parseInt(eachPeptide.getAttributeValue("start_scan"))); //For dtaselect range
        peptide.setDtaEndRange(Integer.parseInt(eachPeptide.getAttributeValue("end_scan"))); //for dtaselect range

        addScores(peptide, eachPeptide);

        String[] dataArr = eachPeptide.getChildText("chro").split(";");
        String[] peakArr = dataArr[0].split(" ");

        data = new ChroData(); //i, Long.parseLong(intArr[1]), Long.parseLong(intArr[2]));

        data.setBsIntensity(dataArr[1].split(" "));
        data.setYsIntensity(dataArr[2].split(" "));
        //data.setYsIntensityReverse( dataArr[2].split(" ") );

        Element fragEle = eachPeptide.getChild("frag");
        String bsText = fragEle.getChildText("bs");
        String ysText = fragEle.getChildText("ys");
        String brText = fragEle.getChildText("br");
        String yrText = fragEle.getChildText("yr");

        peptide.setBsText(bsText);
        peptide.setYsText(ysText);
        peptide.setBrText(brText);
        peptide.setYrText(yrText);

        String[] bsArr = bsText.split(",");
        String[] ysArr = ysText.split(",");

        double[] bsStartMassArr = new double[bsArr.length];
        double[] bsEndMassArr = new double[bsArr.length];
        double[] ysStartMassArr = new double[ysArr.length];
        double[] ysEndMassArr = new double[ysArr.length];

        for (int i = 0; i < bsArr.length; i++) {
          String[] strArr = bsArr[i].split(" ");
          bsStartMassArr[i] = Double.parseDouble(strArr[1]);
          bsEndMassArr[i] = Double.parseDouble(strArr[2]);

          strArr = ysArr[i].split(" ");
          ysStartMassArr[i] = Double.parseDouble(strArr[1]);
          ysEndMassArr[i] = Double.parseDouble(strArr[2]);
        }

        data.setBsStartMass(bsStartMassArr);
        data.setBsEndMass(bsEndMassArr);
        data.setYsStartMass(ysStartMassArr);
        data.setYsEndMass(ysEndMassArr);

        peptide.addData(data);


        protein.addPeptide(peptide);
      }

      if (protein.getPeptideList().size() == 0)
        protein.setRedundant(true);

      list.add(protein);
    }

    //populate all redundant proteins, which has no peptides
    int listSize = list.size();
    for (int i = 0; i < listSize; i++) {
      if (i >= listSize - 1)
        break;

      ChroProtein p = list.get(i);

      if (p.getPeptideList().size() <= 0)
        p.setPeptideList(list.get(i + 1).getPeptideList());
    }

    return list;

  }

  public ArrayList<ChroProtein> getProteinList() throws IOException, Exception {
    return getProteinList(null, 0);
  }

  public ArrayList<ChroProtein> getProteinList(ProgressTask task, int initSize) throws IOException, Exception {

    ArrayList<ChroProtein> list = new ArrayList<ChroProtein>();

    ChroProtein protein;
    ChroPeptide peptide;
    ChroData data;

    List proList = rootEle.getChildren("protein");

    int proSize = proList.size();
    int leftSize = 100 - initSize;
    float seg = leftSize / (float) proSize;
    float progress = (float) initSize;

    for (Iterator<Element> itrPro = proList.iterator(); itrPro.hasNext(); ) {
      Element eachProtein = itrPro.next();

      if (null != task) {
        progress += seg;
        task.updateProgress((int) progress);
      }

      protein = new ChroProtein();
      protein.setLocus(eachProtein.getAttributeValue("locus"));
      protein.setSeqCount(eachProtein.getAttributeValue("seq_ct"));
      protein.setSpectrumCount(eachProtein.getAttributeValue("spec_ct"));


      String ltmpSpec = eachProtein.getAttributeValue("lspec_ct");
      //System.out.println("======" + ltmpSpec);

      if (null != ltmpSpec)
        protein.setLspectrumCount(ltmpSpec);

      String htmpSpec = eachProtein.getAttributeValue("hspec_ct");
      if (null != htmpSpec)
        protein.setHspectrumCount(htmpSpec);

      protein.setSeqCoverage(eachProtein.getAttributeValue("seq_cov"));
      protein.setLength(eachProtein.getAttributeValue("length"));
      protein.setMolWt(eachProtein.getAttributeValue("molwt"));
      protein.setPI(eachProtein.getAttributeValue("pi"));
      protein.setValidation(eachProtein.getAttributeValue("val"));
      protein.setDescription(eachProtein.getAttributeValue("desc"));


      Element redEle = eachProtein.getChild("redundant");
      if (null != redEle) {
        List redunList = redEle.getChildren("protein");

        for (Iterator<Element> itrRPro = redunList.iterator(); itrRPro.hasNext(); ) {

          Element eachRPro = itrRPro.next();

          ChroProtein rprotein = new ChroProtein();
          rprotein.setLocus(eachRPro.getAttributeValue("locus"));
          rprotein.setSeqCount(eachRPro.getAttributeValue("seq_ct"));
          rprotein.setSpectrumCount(eachRPro.getAttributeValue("spec_ct"));
          rprotein.setSeqCoverage(eachRPro.getAttributeValue("seq_cov"));
          rprotein.setLength(eachRPro.getAttributeValue("length"));
          rprotein.setMolWt(eachRPro.getAttributeValue("molwt"));
          rprotein.setPI(eachRPro.getAttributeValue("pi"));
          rprotein.setValidation(eachRPro.getAttributeValue("val"));
          rprotein.setDescription(eachRPro.getAttributeValue("desc"));

          protein.addRedunProtein(rprotein);
        }
      }

      List pepList = eachProtein.getChildren("peptide");

      long freeMemory = Runtime.getRuntime().freeMemory();

      //System.out.println(freeMemory);
      if (freeMemory < 10000) //virtual memory is available less than 1M
        throw new Exception("Memory is not enough. Increase virtual memory in the bat file in Windows or command in Linux. " + freeMemory + " memory left");

      for (Iterator<Element> itrPep = pepList.iterator(); itrPep.hasNext(); ) {
        Element eachPeptide = itrPep.next();

        peptide = new ChroPeptide();
        peptide.setUnique(!"".equals(eachPeptide.getAttributeValue("unique")));
        peptide.setFileName(eachPeptide.getAttributeValue("file"));
        peptide.setScanNum(Integer.parseInt((null != eachPeptide.getAttributeValue("scan")) ? eachPeptide.getAttributeValue("scan") : eachPeptide.getAttributeValue("start_scan")));
        peptide.setSequence(eachPeptide.getAttributeValue("seq"));
        peptide.setXCorr(eachPeptide.getAttributeValue("xcorr"));
        peptide.setChargeState(eachPeptide.getAttributeValue("charge"));
        peptide.setDeltCN(eachPeptide.getAttributeValue("deltaCN"));
        peptide.setSpecCount(eachPeptide.getAttributeValue("spC"));
        peptide.setDeltMass(eachPeptide.getAttributeValue("deltaMass"));
        peptide.setSpRank(eachPeptide.getAttributeValue("spRank"));
        peptide.setSpScore(eachPeptide.getAttributeValue("spScore"));

        addScores(peptide, eachPeptide);

        String lMass = eachPeptide.getAttributeValue("lightAvgMass");
        String hMass = eachPeptide.getAttributeValue("heavyAvgMass");
        if (null != lMass)
          peptide.setLightMass(Double.parseDouble(lMass));
        if (null != hMass)
          peptide.setHeavyMass(Double.parseDouble(hMass));


        //we will use exp type only in the future.  No more many if else.. checking quantlevel or labeled check.


        String[] dataArr;
        String[] massMonitorArr;
        String[] tempArr;
        String[] peakArr;

//		if(true)
//		    continue;
        peptide.setDtaStartRange(Integer.parseInt(eachPeptide.getAttributeValue("start_scan"))); //For dtaselect range
        peptide.setDtaEndRange(Integer.parseInt(eachPeptide.getAttributeValue("end_scan"))); //for dtaselect range

        dataArr = eachPeptide.getChildText("chro").split(";");
        peakArr = dataArr[0].split(" ");
        int startScan = Integer.parseInt(peakArr[1]);
        int endScan = Integer.parseInt(peakArr[2]);


        peptide.setStartRange(peakArr[1]);
        peptide.setEndRange(peakArr[2]);


        int ionCount = dataArr[1].split(" ").length-3;


        int totalScan = 0;
        for (int i = 1; i < dataArr.length; i++) {
          tempArr = dataArr[i].split(" ");
          int currentScan = (int)Double.parseDouble(tempArr[0]);

          if (currentScan >= startScan && currentScan <= endScan)
            totalScan++;
        }


        long[][] intensityArr = new long[ionCount][totalScan];
        double[] retArr = new double[totalScan];
        int tmpIndex=0;

        //dataArr[0]: scan number
        //dataArr[1]: ret time
        //dataArr[2]: total intensity

        for (int i = 1; i < dataArr.length; i++) {
          tempArr = dataArr[i].split(" ");

          int currentScan = (int)Double.parseDouble(tempArr[0]);

          if(currentScan>=startScan && currentScan<=endScan) {
         //   System.out.println(currentScan + " " + startScan + " " + endScan + " " + tempArr[1]);
            retArr[tmpIndex] = Double.parseDouble(tempArr[1]);

            for(int j=3;j<tempArr.length;j++) {
       //     System.out.println(j + " " + tmpIndex+  " " + intensityArr.length);

         //   System.out.println(intensityArr[tmpIndex].length);

        //      try {

                intensityArr[j-3][tmpIndex] = (long)Double.parseDouble(tempArr[j]);

           // System.out.print(" " + tmpIndex);
            //  } catch(Exception e) {
               // e.printStackTrace();
              //}


            }


            tmpIndex++;
/*
            for(int j=0;j<intensityArr[tmpIndex-1].length;j++) {
              for(int k=0;k<intensityArr[tmpIndex-1].length;k++) {
                System.out.print(intensityArr[tmpIndex-1][k] + " ");
              }

              System.out.println(" ");
            } */


          }

        }


        double[] regArr = new double[intensityArr.length];


        HashSet<Integer> zeroIntData = new HashSet();



        for(int i=0;i<intensityArr.length;i++) {


          //gaussian
          double[] tmpPeakArr = new double[intensityArr[i].length];
          double intSum = 0;
          int peakCount=0;
          for(int j=0;j<intensityArr[i].length;j++) {
            tmpPeakArr[j] = (double)intensityArr[i][j];

            if(tmpPeakArr[j]>0) peakCount++;

            intSum += tmpPeakArr[j];
         //   System.out.print("== " + tmpPeakArr[j]);
          }

          if(intSum<=0 || peakCount<3) {
            zeroIntData.add(i);
            peptide.addGaussianPeakModel(null);
            continue;
          }

          GaussianPeakModel gm = GaussianFitting.getGaussianPeakRangeIndex(retArr, tmpPeakArr, 0, intensityArr[i].length-1);
          peptide.addGaussianPeakModel(gm);
         // System.out.println(peptide.getSequence() + " " + gm);
          if(gm.getSigma()<=0 || gm.getPeakArea()<=0) continue;



         //   System.out.print("== " + gm.getPeakArea() + " " + gm.getSigma());
          /*
          boolean peakFound = false;
          for(double d:intensityArr[i]) {
            if (d > 0) {
              peakFound = true;
              break;
            }
          }

          if(!peakFound) {
            regArr[i] = 0;
            continue;
          }
          */
         // GaussianPeakModel gm = GaussianFitting.getGaussianPeakRangeIndex(intensityArr[i], yArr, 0, yArr.length-1);

       //   int start = gm.getPeakStartIndex();
       //   int end = gm.getPeakEndIndex();



          int validRegCount=0;
          for(int j=0;j<intensityArr.length;j++) {
            if(j==i)
              continue;

         //   System.out.println(i + " " + j + " " + intensityArr.length);
            LinearRegression reg = new LinearRegression(intensityArr[i], intensityArr[j], 0, intensityArr.length, 0, true);

            if(reg.getCorr()<=0)
              continue;

            regArr[i] += reg.getCorr();

/*
            if(!zeroIntData.contains(j)) {
              for (long l : intensityArr[i])
                System.out.print(l + " ");
              System.out.println("");

              for (long l : intensityArr[j])
                System.out.print(l + " ");
              System.out.println("");


              System.out.println(reg.getCorr());


            } */
            validRegCount++;
          }


          regArr[i] /= validRegCount;


        }
/*
        for(int i=0;i<regArr.length;i++) {
       //   if(regArr[i]>0)
       //     System.out.println(regArr[i]); // + " " + validRegCount + " " + regArr[i]/validRegCount);

          regArr[i] /= validRegCount;
        }
*/


        peptide.setDiaFragRegressionArr(regArr);

        /*
        for(int i=0;i<regArr.length;i++) {
          if(regArr[i]<=0) continue;

          System.out.println(regArr[i] + " " + i);
          for(double d:intensityArr[i]) {
            System.out.print(d + " " );
          }

          System.out.println();
        }
        */


        protein.addPeptide(peptide);

      }

      if (protein.getPeptideList().size() == 0)
        protein.setRedundant(true);

      list.add(protein);
    }

    //populate all redundant proteins, which has no peptides
    int listSize = list.size();
    for (int i = 0; i < listSize; i++) {
      if (i >= listSize - 1)
        break;

      ChroProtein p = list.get(i);

      if (p.getPeptideList().size() <= 0)
        p.setPeptideList(list.get(i + 1).getPeptideList());
    }

    return list;
  }

  public boolean isDataDependent() {
    return this.isDataDependent;
  }
    /*
    public Iterator<ChroProtein> getProteins() throws IOException {
        return new Iterator<ChroProtein>() {
            private ChroProtein protein;
            private ChroPeptide peptide;

            public boolean hasNext() {
                return lastLine != null; // && !lastLine.startsWith("\tProteins\t");
            }

            public ChroProtein next() {

                try {
                    protein = getProtein();
                }
                catch (IOException e) {
                    e.printStackTrace();
                }

                return protein;
            }

            public void remove() {
                throw new UnsupportedOperationException("Not supported");
            }

            private ChroProtein getProtein() throws IOException {

                //String[] strArr = lastLine.split("\t");
                //protein = new Protein(strArr);
                //System.out.println(lastLine);
                lastLine=br.readLine();

                protein = new ChroProtein(lastLine);
                peptide = null;

                while( (lastLine=br.readLine())!=null && !lastLine.startsWith("[PRO") )
                {
                    if(lastLine.startsWith("[PEP"))
                    {
                        if(null != peptide)
                            protein.addPeptide(peptide);

                        peptide = new ChroPeptide(lastLine=br.readLine());

                        lastLine=br.readLine();//Read away [CHROMATOGRAMS] line
                        String[] peakRange = (lastLine=br.readLine()).split("\t"); //read away Peak line;
                        if(peakRange.length==3)
                        {
                            peptide.setStartRange(peakRange[1]);
                            peptide.setEndRange(peakRange[2]);
                        }
                        //Read P line here
                        //lastLine=br.readLine();//Read the first data line

                    }
                    else// if(!lastLine.startsWith("[CHR"))
                    {
                        String[] str=lastLine.split("\t");
                        data = new ChroData(Integer.parseInt(str[0]), Integer.parseInt(str[1]), Integer.parseInt(str[2]) );
                        peptide.addData(data);
                    }
                }

                protein.addPeptide(peptide);

                return protein;
            }
        };
    }
       */

  public static void main(String args[]) throws Exception {
    ChroXmlReader cr = new ChroXmlReader(args[0]);

    ArrayList<ChroProtein> list = cr.getProteinList();

    System.out.println(list.size());
    System.out.println(list.get(0).getProteinLine());

//        List pepList = list.get(0).getPeptideList();

    ChroPeptide peptide;
    for (int i = 0; i < list.size(); i++) {
      ChroProtein protein = (ChroProtein) list.get(i);
      List pepList = protein.getPeptideList();
      for (int j = 0; j < pepList.size(); j++) {
        peptide = (ChroPeptide) pepList.get(j);
        //peptide.getF
        //System.out.println( peptide.getDataList().get(0) );
        System.out.println(peptide);


      }
    }

    ArrayList sampleList = cr.getSampleList();
    for (Iterator<String> itr = sampleList.iterator(); itr.hasNext(); ) {
      String each = itr.next();
    }
  }

  public int getQuantLevel() {
    return quantLevel;
  }

  public void setQuantLevel(int quantLevel) {
    this.quantLevel = quantLevel;
  }

  public boolean isLabeled() {
    return labeled;
  }

  public void setLabeled(boolean labeled) {
    this.labeled = labeled;
  }

  public ArrayList<String> getFileList() {
    return fileList;
  }

  public void setFileList(ArrayList<String> fileList) {
    this.fileList = fileList;
  }

  public String getSampleName(String fileName) {
    return this.fileSampleHt.get(fileName);
  }

  public String getSampleName(int fileIndex) {
    String curDir = getFileList().get(fileIndex);
    return getSampleName(curDir);
  }


  public String getFileName(int index) {
    String tmpFileName = this.fileList.get(index);

    return tmpFileName.substring(tmpFileName.lastIndexOf(File.separator) + 1);
  }

  public ArrayList getSampleList() {
    return sampleList;
  }

  public void setSampleList(ArrayList sampleList) {
    this.sampleList = sampleList;
  }

  public Hashtable<String, Sample> getSampleObjList() {
    return sampleObjList;
  }

  public void setSampleObjList(Hashtable<String, Sample> sampleObjList) {
    this.sampleObjList = sampleObjList;
  }

  public ArrayList<Sample> getSampleExpList() {
    return sampleExpList;
  }

  public void setSampleExpList(ArrayList<Sample> sampleExpList) {
    this.sampleExpList = sampleExpList;
  }

  public int getExpType() {
    return expType;
  }

  public void setExpType(int expType) {
    this.expType = expType;
  }

  public void addScores(ChroPeptide chroPeptide, Element peptide) {

    List scrList = peptide.getChildren("search_score");
    for (Iterator<Element> scrItr = scrList.iterator(); scrItr.hasNext(); ) {
      Element eachScr = scrItr.next();
      chroPeptide.addScore(eachScr.getAttributeValue("name"), eachScr.getAttributeValue("value"));
    }
  }

  public String getQuantType() {
    return quantType;
  }

  public void setQuantType(String quantType) {
    this.quantType = quantType;
  }

  public static HashMap<String, ChroPeptide> getPeptideMap(String fileName) {//key is the Sequence_cs
    ChroXmlReader cr = null;
    HashMap<String, ChroPeptide> peptideMap = new HashMap<>();
    try {
      cr = new ChroXmlReader(fileName);
      for (ChroProtein currentPRotein : cr.getProteinList()) {
        for (ChroPeptide currentPeptide : (List<ChroPeptide>) currentPRotein.getPeptideList()) {
          peptideMap.put(currentPeptide.getSequence() + "_" + currentPeptide.getChargeState(), currentPeptide);

        }
      }

    } catch (JDOMException ex) {
      Logger.getLogger(ChroXmlReader.class.getName()).log(Level.SEVERE, null, ex);
    } catch (Exception ex) {
      Logger.getLogger(ChroXmlReader.class.getName()).log(Level.SEVERE, null, ex);
    }
    return peptideMap;
  }


}


