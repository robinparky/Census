/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.scripps.pms.census.tools;

import edu.scripps.pms.census.ElementComposition;
import edu.scripps.pms.census.conf.Configuration;
import edu.scripps.pms.census.hash.IndexedFile;
import edu.scripps.pms.census.io.IsotopeReader;
import edu.scripps.pms.census.labelFree.LabelFreeCalcUtil;
import edu.scripps.pms.census.labelFree.SpectrumModel;
import edu.scripps.pms.census.model.IrisDataModel;
import edu.scripps.pms.census.util.CalcUtil;
import edu.scripps.pms.census.util.IsotopeDist;
import edu.scripps.pms.census.util.io.SpectrumReader;
import edu.scripps.pms.util.spectrum.Hline;
import edu.scripps.pms.util.spectrum.Peak;
import edu.scripps.pms.util.spectrum.PeakList;
import gnu.trove.TDoubleIntHashMap;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.jdom.Document;
import org.jdom.Element;
import org.jdom.input.SAXBuilder;

/**
 *
 * @author Harshil
 */
public class IsotopeTool {

    private final static double PROTON_MASS = 1.00728;

    public static void main(String[] args) throws Exception {
        String configFileName = "e:\\census_config.xml";
//        String configFileName = "/data/1/rpark/deploy/census_config_repository/census_config_silac_new.xml";
        //String sequence = "LPSSLDQNVPQYKR";
        String sequence = "LPSSLDQNVPQYK(114.042927)R";

        /*
         970.3814434
         971.3847973999999
         972.3881514
         973.3915053999999
         974.3948594
         975.3982133999999
         976.4015674
         977.4049213999999
         978.4082754
         979.4116293999999
         */
        double[] distArr = IsotopeTool.getIsotopePeaks(configFileName, sequence, 3);
        
        
        

        //586.9759606
                
    //            586.966
        //586.9758606
     //   586.9760606

//        double[] distArr = {586.9759606, 587.3104119333333, 587.6448632666667, 587.9793146000001};

        
/*
        for (double d : distArr) {
            if (d <= 0) {
                break;
            }
	
            System.out.println(d);
        }

*/
//        double[] temp = {350.2410,1987.2134};
//114.042927

	 String ms1FileName = "/home/rpark/rpark_on_data/project/data_analysis/khatereh/sample6505_projects2014_05_26_12_1818/Sample6505_1.ms1";



//        List<IsotopeModel> isotopeModel = getIntensity(20.0, 55.0, distArr, 0, 0.005, ms1FileName);

        //String fileName = "e:\\isoTopeResult.txt";
//        write(ms1FileName + "_quant_out.txt",isotopeModel);

        
        /*
         IsotopeReader isotopeReader = getIsotopePeaks(configFileName, sequence);
         for(String key  : isotopeReader.getIsotope().keySet())
         {
         System.out.println(key + "-->"+Arrays.toString(isotopeReader.getIsotope().get(key)));

         }*/
    }

    public static void write(String fileName, List<IsotopeModel> isotopeList)
    {
        BufferedWriter bw = null;
        try {
            bw = new BufferedWriter(new FileWriter(fileName));
            bw.write("ScanNumber\tr_time\tionInjection_time\tmassList\tintensityList\n");
            for(IsotopeModel isotopeModel : isotopeList)
            {
                bw.write(isotopeModel.getScanNumber()+"\t");
                bw.write(isotopeModel.getRetentionStartTime()+"\t");
                bw.write(isotopeModel.getIonInjectionTime()+"\t");
                for(double number : isotopeModel.getIsoArr())
                {
                    bw.write(number+",");
                }
                bw.write("\t");
                for(double number : isotopeModel.getIntensityArr())
                {
                    bw.write(number+",");
                }
                bw.write("\n");
//                bw.write(Arrays.toString(isotopeModel.getIsoArr())+"\t");
//                bw.write(Arrays.toString(isotopeModel.getIntensityArr())+"\n");
            }
        } catch (IOException ex) {
            Logger.getLogger(IsotopeTool.class.getName()).log(Level.SEVERE, null, ex);
        }
        finally{
            try {
                bw.close();
            } catch (IOException ex) {
                Logger.getLogger(IsotopeTool.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        
    }
    
    public static double[] getIsotopePeaks(String configFileName, String sequence, int cs) throws Exception {

        double[] isoArr = IsotopeTool.getIsotopePeaks(configFileName, sequence);

        for (int i = 0; i < isoArr.length; i++) {
            isoArr[i] = (isoArr[i] + cs * PROTON_MASS) / cs;
        }

        return isoArr;
    }

/*    public static List<IsotopeModel> getIntensity(double retentionStartTime, double retentionEndTime, double[] isoArr, double modMass, double tolerence, String ms1FileName) 
    {
        List<IsotopeModel> isotopeList = new ArrayList();
        IsotopeModel isotopeModel = null;

	if(modMass>0) {
		for(int i=0;i<isoArr.length;i++) {
			isoArr[i] += modMass;
		}
	}

        try 
        {
            String line = null;
            //SpectrumReader sr = new SpectrumReader("e:\\Sample6503_1.ms1", "ms1");
            SpectrumReader sr = new SpectrumReader(ms1FileName, "ms1");
            Iterator<PeakList> it = sr.getSpectra();
            while (it.hasNext()) 
            {
                double[] intensityArr = new double[isoArr.length];
                int counter = 0;

                PeakList list = it.next();
                if (retentionStartTime <= list.getRetentionTime() && retentionEndTime>= list.getRetentionTime() ) 
                {
//			System.out.println(list.getHiscan());
//                    System.out.println(list.getRetentionTime());
                    List<Peak> peakList = list.getPeakList();
                    List<Double> massList = new ArrayList<>();
                    List<Double> intensityList = new ArrayList<>();

                    
                    for (Peak currentPeak : peakList) 
                    {
                        massList.add(currentPeak.getM2z());
                        intensityList.add(currentPeak.getIntensity());
                    }
                 
                    for(double massValue : isoArr)
                    {
                        if(massList.contains(massValue))
                        {
                            intensityArr[counter] = intensityList.get(massList.indexOf(massValue));
                            counter++;
                        }
                        else
                        {
                            for(int i =0; i<massList.size();i++)
                            {
                                double actualMassValue = massList.get(i);
          //                      if(actualMassValue>586)
          //                          System.out.println("586....");
                                
                               if(actualMassValue <=massValue+tolerence && actualMassValue >= massValue-tolerence)
                               {
//				System.out.println(intensityArr.length + "\t" + counter);
//                                   if( counter < intensityArr.length && intensityArr[counter] == 0)
//                                       intensityArr[counter++] += intensityList.get(i);
//                                   else 
//                                       intensityArr[counter] += intensityList.get(i);
//                                   if((int)(intensityList.get(i)- 189312.9) == 0)
//                                           System.out.println("");
                                intensityArr[counter] += intensityList.get(i);
                               }
                            }
                            counter++;
                        }
                    }
                    isotopeModel = new IsotopeModel(list.getRetentionTime(), retentionEndTime, isoArr,list.getInjectionTime() , intensityArr,list.getHiscan());
                    isotopeList.add(isotopeModel);
                }
                else if(retentionEndTime < list.getRetentionTime())
                    break;
                
            }
            

        } catch (IOException ex) {
            Logger.getLogger(IsotopeTool.class.getName()).log(Level.SEVERE, null, ex);
        }
        return isotopeList;
    }
    */

    public static double[] getIsotopePeaks(String configFileName, String sequence) throws Exception {
        
//        SAXBuilder builder = new SAXBuilder();
//
//        Document doc = builder.build(new File(configFileName));
//        
//        Element rootEle = doc.getRootElement();
//        
        //           IsotopeReader isoReader = new IsotopeReader(rootEle);

        char[] ch = sequence.toCharArray();
        IsotopeReader isoReader = null;
        
        Configuration conf = Configuration.getInstance();
        if(!conf.isReadConfigFile())
            conf.readXMLParam(configFileName);
        isoReader = conf.getIsotopeReader();
            
        ElementComposition element = new ElementComposition(ch, 0, ch.length, isoReader.getIsotope());
        element.calculateSample();
       // int[] isoArr = element.getElementSampleArr();

        IsotopeDist dist = new IsotopeDist(element.getElementSampleArr(), element.getModShift(), true);

        //return dist.getMasslist();
        return dist.getHighMassList();
        //return isoReader;

    }
    
    /*Harshil : WOkred Till MS1Read in file
    public static List<IsotopeModel> getIntensityChro(double retentionStartTime, double retentionEndTime, double[] isoArr, double modMass, double tolerence, List<PeakList> peaks)
    {

        List<IsotopeModel> isotopeList = new ArrayList();
        IsotopeModel isotopeModel = null;

        System.out.println("New Entry-> Time ->"+ (retentionEndTime-retentionStartTime));
        Date startData = new Date();
        if(modMass>0) {
                for(int i=0;i<isoArr.length;i++) {
                        isoArr[i] += modMass;
                }
        }

        int  startIndex =getBinarySearchIndex(peaks, retentionStartTime);
        int endIndex  =getBinarySearchIndex(peaks, retentionEndTime);
        String line = null;
//        Iterator<PeakList> it = peaks.iterator();
//        while (it.hasNext())
        for(int z=startIndex;z<=endIndex;z++)
        {
            double[] intensityArr = new double[isoArr.length];
            int counter = 0;

            PeakList list = peaks.get(z);
            if (retentionStartTime <= list.getRetentionTime() && retentionEndTime>= list.getRetentionTime() )
            {
//                      System.out.println(list.getHiscan());
//                    System.out.println(list.getRetentionTime());

                List<Peak> peakList = list.getPeakList();
                List<Double> massList = new ArrayList<>();
                List<Double> intensityList = new ArrayList<>();


                for (Peak currentPeak : peakList)
                {
                    massList.add(currentPeak.getM2z());
                    intensityList.add(currentPeak.getIntensity());
                }

                for(double massValue : isoArr)
                {if(massList.contains(massValue))
                    {
                        intensityArr[counter] = intensityList.get(massList.indexOf(massValue));
                        counter++;
                    }
                    else
                    {
                        for(int i =0; i<massList.size();i++)
                        {
                            double actualMassValue = massList.get(i);
                            //                      if(actualMassValue>586)
                            //                          System.out.println("586....");

                            if(actualMassValue <=massValue+tolerence && actualMassValue >= massValue-tolerence)
                            {
//                              System.out.println(intensityArr.length + "\t" + counter);
//                                   if( counter < intensityArr.length && intensityArr[counter] == 0)
//                                       intensityArr[counter++] += intensityList.get(i);
//                                   else 
//                                       intensityArr[counter] += intensityList.get(i);
//                                   if((int)(intensityList.get(i)- 189312.9) == 0)
//                                           System.out.println("");
                                intensityArr[counter] += intensityList.get(i);
                            }
                        }
                        counter++;
                    }
                }
                isotopeModel = new IsotopeModel(list.getRetentionTime(), retentionEndTime, isoArr,list.getInjectionTime() , intensityArr,list.getHiscan());
                isotopeList.add(isotopeModel);
            }
            else if(retentionEndTime < list.getRetentionTime())
                break;

        }
        Date endData = new Date();
        double difference = endData.getTime() - startData.getTime();
        System.out.println("runTime : " + difference);
        return isotopeList;
    }
    
    
    
private static int getBinarySearchIndex(List<PeakList> aList,double time)
    {
        int keyIndex=-1;
        Comparator<PeakList> peakComaprator = new Comparator<PeakList>() {

            @Override
            public int compare(PeakList o1, PeakList o2) {
                return ((Double)o1.getRetentionTime()).compareTo((Double)o2.getRetentionTime());
            }
        };

        keyIndex = Collections.binarySearch(aList, new PeakList(time), peakComaprator);
        if(keyIndex<0) //Cannot find index
            keyIndex=-(++keyIndex); //Math.abs(++keyIndex);

        if(keyIndex>=aList.size())
            keyIndex--;

        return keyIndex;

    }

    */
    /**
     * 
     * @param retentionStartTime
     * @param retentionEndTime
     * @param isoArr
     * @param modMass
     * @param iFile
     * @return 
     */
    public static List<IsotopeModel> getIntensityChro(double retentionStartTime, double retentionEndTime, double[] isoArr, double modMass, IndexedFile iFile)
    {

        List<IsotopeModel> isotopeList = new ArrayList();
        IsotopeModel isotopeModel = null;
        //DEBUG purpose
//        System.out.println("New Entry: Time ->"+retentionEndTime +"-" + retentionStartTime + " = " + (retentionEndTime-retentionStartTime));
        
//        System.out.print("Search window: ->" + (retentionEndTime-retentionStartTime));
        if(modMass>0) {
                for(int i=0;i<isoArr.length;i++) {
                        isoArr[i] += modMass;
                }
        }

        int  startScanNumber =iFile.getRetentonToScanMap().get(retentionStartTime);
        int endScanNumber  =iFile.getRetentonToScanMap().get(retentionEndTime);
        
        if(startScanNumber>endScanNumber)
        {
            int temp = startScanNumber;
            startScanNumber = endScanNumber;
            endScanNumber = temp;
        }
        return getIntensityChroByScanNumber(startScanNumber, endScanNumber, isoArr, modMass, iFile);
    }
    
/**
 * 
 * @param startScanNumber
 * @param endScanNumber
 * @param isoArr
 * @param modMass
 * @param iFile
 * @return 
 */
public static List<IsotopeModel> getIntensityChroByScanNumber(int  startScanNumber, int  endScanNumber, double[] isoArr, double modMass, IndexedFile iFile)
    {
        double massTolerance = Configuration.getInstance().getMassTolerance();
        List<IsotopeModel> isotopeList = new ArrayList();
        IsotopeModel isotopeModel = null;
        int[] keys = iFile.getKeys();
        
        int startIndex =iFile.getIndexFromScanNumber(startScanNumber);
        int endIndex = iFile.getIndexFromScanNumber(endScanNumber);

        Date startDate = new Date();
        if(modMass>0) {
                for(int i=0;i<isoArr.length;i++) {
                        isoArr[i] += modMass;
                }
        }

  
//
        String line = null;
//        Iterator<PeakList> it = peaks.iterator();
//        while (it.hasNext())
        for(int z=startIndex;z<=endIndex;z++)
        {
            double[] intensityArr = new double[isoArr.length];
            int counter = 0;
            List<Double> massList = new ArrayList<>();
            List<Double> intensityList = new ArrayList<>();
            
            SpectrumModel spectrumModel = null;
            try {
                 spectrumModel = LabelFreeCalcUtil.readSpectrumPeak(iFile, z, iFile.getKeys());
            } catch (IOException ex) {
                Logger.getLogger(IsotopeTool.class.getName()).log(Level.SEVERE, null, ex);
            }
                for(int i =0; i<spectrumModel.getMass().length;i++)
                {
                    massList.add(spectrumModel.getMass()[i]);
                    intensityList.add(spectrumModel.getIntensity()[i]);
                }
                for(double massValue : isoArr)
                {
                    double tolerence = massValue/1000 *  massTolerance ;
//                    System.out.println( tolerence +" ---" + massValue);
                    if(massList.contains(massValue))
                    {
                        intensityArr[counter] = intensityList.get(massList.indexOf(massValue));
                        counter++;
                    }
                    else
                    {
                        for(int i =0; i<massList.size();i++)
                        {
                            double actualMassValue = massList.get(i);
                            //                      if(actualMassValue>586)
                            //                          System.out.println("586....");
                            
                            if(actualMassValue <=massValue+tolerence && actualMassValue >= massValue-tolerence)
                            {
//                              System.out.println(intensityArr.length + "\t" + counter);
//                                   if( counter < intensityArr.length && intensityArr[counter] == 0)
//                                       intensityArr[counter++] += intensityList.get(i);
//                                   else 
//                                       intensityArr[counter] += intensityList.get(i);
//                                   if((int)(intensityList.get(i)- 189312.9) == 0)
//                                           System.out.println("");
                                intensityArr[counter] += intensityList.get(i);
                            }
                        }
                        counter++;
                    }
                }
                isotopeModel = new IsotopeModel(isoArr,iFile.getRetentionTimeMap().get(spectrumModel.getScanNumber()),
                        spectrumModel.getIonInjectionTime(), intensityArr,spectrumModel.getScanNumber());
                isotopeList.add(isotopeModel);
           
            

        }
        Date endData = new Date();
        
        double difference = endData.getTime() - startDate.getTime();
        //DEBUG purpose
        
//        System.out.print(" ->runTime : " + difference +"ms isotopes: " +isotopeList.size());
//        if(isotopeList.size()==0)
//            System.out.println("error: isotopeModeLlist is empty....");
        return isotopeList;
    }
    
                    
    
    
}
