
/*
* Copyright (c) 2008 The Scripps Research Institute, Yates Lab.  All rights reserved.
*/

package edu.scripps.pms.census.util;

/**
 *
 * @author rohan
 * @author Sung Kyu, Robin, Park
 * @email rpark@scripps.edu
 * Created on May 23, 2008
 * $Revision: 1.18 $
 * $Date: 2014/08/27 18:00:35 $
 */

import java.util.*;
import java.io.*;

import org.dom4j.*;

import edu.scripps.pms.census.model.*;
import edu.scripps.pms.census.conf.Configuration;
import edu.scripps.pms.census.hash.IndexedFile;
import edu.scripps.pms.util.seq.Fasta;

import edu.scripps.pms.census.util.io.FastaReader;
import edu.scripps.pms.census.model.ChroPeptide;
import edu.scripps.pms.census.tools.Formatter;
import org.dom4j.Element;
import org.jdom.*;
import org.jdom.Document;
import org.jdom.input.SAXBuilder;

public class MergeLabelFreeChro {

	public static void readConfig(Configuration conf, String configFile) throws Exception {


		SAXBuilder sb = new SAXBuilder();
		Document doc = sb.build(configFile);
		final org.jdom.Element root = doc.getRootElement();
		conf.setRootConfEle(root);

		if(null == conf.getFilePath() || "".equals(conf.getFilePath()))
		{
			int fIndex = -1;

			if(File.separator.equals("/"))  //unix
			{
				fIndex = configFile.lastIndexOf("/");
			}
			else //windows
			{
				fIndex = configFile.lastIndexOf("\\");
			}

			if( fIndex <0 ) //current folder
				conf.setFilePath(".");
			else
				conf.setFilePath(configFile.substring(0, fIndex));
		}



		List sampleEleList = root.getChildren("sample");

		//Configuration.
		String refFileName = "";
		org.jdom.Element refEle = root.getChild("ref");
		if(null != refEle)
		{
			//String refSamName = refEle.getChildText("sample_name");
			refFileName = refEle.getChildText("file_name");

			conf.setRefFileName(refFileName);
		}


		int index=0;
		int refIndex=0;

		Vector<String> fileNameList = new Vector<String>();
		Vector<String> sampleNameList = new Vector<String>();
		Vector<String> pathFileNameList = new Vector<String>();

		//confSam.setRefFileName(refFileName);
		//Hashtable<String, Hashtable> masterHt = new Hashtable<String, Hashtable>();
		HashSet set = new HashSet();

		for(Iterator<org.jdom.Element> itr = sampleEleList.iterator(); itr.hasNext(); ) {
			org.jdom.Element sam = itr.next();
			//String samName = sam.getChildText("name");
			String samName = sam.getAttributeValue("group");
			//populate configuration class
			Configuration.Sample confSam = new Configuration.Sample();
			confSam.setName(samName);

			List<String> filesList = new ArrayList<String>();

			List<org.jdom.Element> sampleList = sam.getChildren("each_sample");
			for (Iterator<org.jdom.Element> eachSampleItr = sampleList.iterator(); eachSampleItr.hasNext(); ) {
				org.jdom.Element eachSample = eachSampleItr.next();
				List<org.jdom.Element> fileList = eachSample.getChild("ms_files").getChildren("file");


				SampleModel sampleModel = new SampleModel(samName);

				if (fileList.size() > 0) {
					String firstFileName = fileList.get(0).getText();
					String filePath = firstFileName.substring(0, firstFileName.lastIndexOf(File.separator));
					firstFileName = firstFileName.substring(firstFileName.lastIndexOf(File.separator) + 1);

					if (firstFileName.startsWith("*")) {
						String extension = firstFileName.substring(firstFileName.lastIndexOf(".") + 1);

						edu.scripps.pms.census.util.RelExFileFilter fFilter = new edu.scripps.pms.census.util.RelExFileFilter(extension);

						File specFile = new File(filePath);
						String[] splist = specFile.list(new edu.scripps.pms.census.util.RelExFileFilter(extension));

						for (String eachFileName : splist) {
							eachFileName = filePath + File.separator + eachFileName;
							set.add(filePath);
							filesList.add(filePath);
							sampleModel.addPath(filePath);

							confSam.addFile(eachFileName);
							pathFileNameList.add(eachFileName);
							fileNameList.add(eachFileName);
							sampleNameList.add(samName);
						}
					} else {
						for (Iterator<org.jdom.Element> itr1 = fileList.iterator(); itr1.hasNext(); ) {
							org.jdom.Element eachFile = itr1.next();
							String fileName = eachFile.getText();

							confSam.addFile(fileName);
							pathFileNameList.add(fileName);

							set.add(fileName.substring(0, fileName.lastIndexOf(File.separator)));
							filesList.add(fileName.substring(0, fileName.lastIndexOf(File.separator)));
							sampleModel.addPath(fileName.substring(0, fileName.lastIndexOf(File.separator)));

							if (fileName.endsWith("ms2")) {
								fileName = fileName.substring(0, fileName.length() - 3);
								fileName += "ms1";
							}

							fileNameList.add(fileName);
							sampleNameList.add(samName);

							if (null != refEle && refFileName.equals(fileName))
								refIndex = index;

							index++;
						}
					}

					conf.addExp(confSam);
					conf.addSample(sampleModel);
					//List<Element> fileList = sam.getChild("ms_files").getChildren("file");
				}
			}
			conf.getNonlabelFilenameGroupMap().put(samName, filesList);
		}
	}

	public static void main(String[] args) throws Exception {
		Configuration conf = Configuration.getInstance();
		//conf.setFilePath("/data/2/rpark/ip2_data/rpark/Xianyin_Lai/labelfree_quant/labelfree_12279");
		//String confFile = "/data/2/rpark/ip2_data/rpark/Xianyin_Lai/labelfree_quant/labelfree_12279/census_config_labelfree_12279.xml";
		conf.setFilePath("/home/rpark/test_data/carol");
		String confFile = "/home/rpark/test_data/carol/census_config_labelfree_16383.xml";

		MergeLabelFreeChro.readConfig(conf, confFile);

		System.out.println(conf.getSampleList());
		System.out.println(conf.getRootConfEle());


		MergeLabelFreeChro.mergeLabelFreeChro(conf.getSampleList(), "census_chro_temp.xml");

	}

    public MergeLabelFreeChro() {

    }
//, Hashtable<String,IndexedFile> ht
    public static void mergeLabelFreeChro(List<SampleModel> sampleList, String chroTmpFileName) throws DocumentException, Exception {

        Configuration conf = Configuration.getInstance();

	String outputFilename = conf.getOutputFilename();
	if(null == outputFilename)
	    outputFilename = "census_labelfree_out.txt";

	String outputFilenameTmp = outputFilename + "tmp";

        HashMap<String, IndexedFile> ht = new HashMap<>();
        for(SampleModel sModel : sampleList)
        {
            for(String path : sModel.getPathList())
            {
                File f = new File(path);
                String[] list = f.list(new RelExFileFilter("ms1"));
                for(String value : list)
                {
                    File  indexFile = new File(path+File.separator+value+".index");
					//if(!indexFile.exists()) {
						//create index file
					//}
                    ht.put(value, new IndexedFile(indexFile, path+File.separator+value));
                }

            }

        }
        List<ChroLabelfreeMergeProteinModel> allList = Collections.EMPTY_LIST; //new Vector<MergeProteinModel>();

	Hashtable<String, MultipleChroProtein> proteinHt = new Hashtable<String, MultipleChroProtein>();

	int expSize=0;
        int index=0;

	//BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream("census_label_free_result_with_id_temp.txt"));
	BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(outputFilenameTmp));
	PrintStream p = new PrintStream(out);

	List<org.jdom.Element> samGroupEleList = conf.getRootConfEle().getChildren("sample");
	for(Iterator<org.jdom.Element> samgItr=samGroupEleList.iterator(); samgItr.hasNext(); )
	{
	    org.jdom.Element groupEle = samgItr.next();
            p.print("H\tGROUP_SAMPLE\t");
            p.print(groupEle.getAttributeValue("group"));
            p.print("\t");

            List<org.jdom.Element> sampleEleList = groupEle.getChildren("each_sample");

            for(Iterator<org.jdom.Element> samItr=sampleEleList.iterator(); samItr.hasNext(); )
            {
                org.jdom.Element eachSample = samItr.next();
                p.print(eachSample.getAttributeValue("name"));
                p.print("\t");

                expSize++;
            }

            p.println();
	}

	p.print("PLINE\tACCESSION\tDESCRIPTION\t");
	for(int i=0;i<expSize;i++) {
	    p.print("SCOUNT_" + (i+1));
	    p.print("\t");
	}
//	for(int i=0;i<expSize;i++) {
//	    p.print("NORM_SCOUNT_" + (i+1));
//	    p.print("\t");
//	}
	p.print("PEP_COUNT\t");

/*
	for(int i=0;i<expSize;i++) {
	    p.print("AVG_INTENSITY_" + (i+1));
	    p.print("\t");
	}
*/

	for(int i=0;i<expSize;i++) {
	    p.print("NORM_INTENSITY_" + (i+1));
	    p.print("\t");
	}

	p.println("");

	p.print("SLINE\t");
	for(int i=0;i<expSize;i++) {
	    p.print("EXP_" + (i+1));
	    p.print("\t");
	    p.print("SEQUENCE\tFILENAME\tSCAN\tCSTATE\tINTENSITY\tCORRIONINJECTION_INTENSITY\tPROFILE_SCORE\tMHPLUS\tCALCMHPLUS\tTOTALINTENSITY\tXCORR\tDCN\tDMASS\tSPRANK\tSPSCORE\tREDUNDANCY\tSTARTRANGE\tENDRANGE\tRETENTIONTIME\tIONINJECTIONTIME\t");

	}
	/*
	for(int i=0;i<expSize;i++) {
	    p.print("NORM_AVG_INTENSITY" + (i+1));
	    p.print("\t");
	}*/

	for(int i=0;i<expSize;i++) {
	    p.print("NORM_INTENSITY" + (i+1));
	    p.print("\t");
	}

	p.println("");
	/****************************************
	    Hashtable<String, MultipleChroProtein>

	    MultipleChroProtein
	    |
	    |---- ChroProtein[]
	    |---- Hashtable<String, ChroPeptide[]> : key is (peptide_seq + charge_state)

	*****************************************/

        for(Iterator<SampleModel> samItr=sampleList.iterator(); samItr.hasNext(); )
	{
	    SampleModel sModel = samItr.next();

	    List<String> fList = sModel.getPathList();

	    if(fList.size()<=0)
	    {
		System.out.println("ms1 files are not found");
		System.exit(0);
	    }

//	    for(Iterator<String> fItr=fList.iterator(); fItr.hasNext(); ) {
	    String eachPath = fList.get(0);

	    //		String eachPath = fItr.next();

	    eachPath += File.separator;
	    eachPath += chroTmpFileName;

	    System.out.println("Reading "+ eachPath);
	    org.dom4j.Element rootEle = Dom4jUtil.getRootEle(new File(eachPath ));

	    //Hashtable<String, ChroPeptide> pepHt = new Hashtable<String, ChroPeptide>();

	    //treeWalk(ele);
	    List prolist = rootEle.selectNodes("protein");
	    List chroProteinList = new ArrayList();

	    for(Iterator<Element> proItr=prolist.iterator(); proItr.hasNext(); )
	    {
		Element proEle = proItr.next();
		String locus = proEle.attributeValue("locus");
		String desc = proEle.attributeValue("desc");
		if(locus.startsWith("Reverse"))
		    continue;

		MultipleChroProtein mprotein = proteinHt.get(locus);

		if(null == mprotein) {
		    mprotein = new MultipleChroProtein(expSize);
		    proteinHt.put(locus, mprotein);
		}

		ChroProtein cPro = new ChroProtein();
		cPro.setLocus(locus);
		cPro.setSeqCount(proEle.attributeValue("seq_ct"));
		cPro.setSpectrumCount(proEle.attributeValue("spec_ct"));
		//cPro.setDescription(proEle.attributeValue("desc"));
		cPro.setDescription(desc);

		mprotein.addChroProtein(cPro, index);
		mprotein.setDescription(desc);

		List pepList = proEle.selectNodes("peptide");

		for(Iterator<Element> pepItr=pepList.iterator(); pepItr.hasNext(); )
		{
		    Element pepEle = pepItr.next();
		    ChroPeptide peptide = new ChroPeptide();
		    peptide.setSequence(pepEle.attributeValue("seq"));
		    peptide.setFileName(pepEle.attributeValue("file"));
		    peptide.setScanNum(Integer.parseInt(pepEle.attributeValue("scan")));
		    peptide.setXCorr(pepEle.attributeValue("xcorr"));
		    peptide.setDeltCN(pepEle.attributeValue("deltaCN"));

                    peptide.setSpRank(pepEle.attributeValue("spRank"));
                    peptide.setSpScore(pepEle.attributeValue("spScore"));
                    //System.out.println("--------->>" + pepEle.attributeValue("spScore"));
                    //System.out.println("--------->>" + pepEle.attributeValue("redundancy"));
                    peptide.setDeltMass(pepEle.attributeValue("deltaMass"));
		    peptide.setChargeState(pepEle.attributeValue("charge"));

                    peptide.setGaussianPeakArea((long)Double.parseDouble(pepEle.attributeValue("peak_area")));


		    if(null == pepEle.attributeValue("peak_areaIonInjectionCorrection"))
			    peptide.setGaussianPeakAreaIonInjectionCorrection(0);
		    else
			    peptide.setGaussianPeakAreaIonInjectionCorrection((long) Double.parseDouble(pepEle.attributeValue("peak_areaIonInjectionCorrection")));
                 //   peptide.setDtaStartRange( Integer.parseInt(pepEle.attributeValue("start_scan")) ); //For dtaselect range


		 //   peptide.setDtaStartRange( Integer.parseInt(pepEle.attributeValue("start_scan")) ); //For dtaselect range
		 //   peptide.setDtaEndRange( Integer.parseInt(pepEle.attributeValue("end_scan")) ); //for dtaselect range
                    peptide.setMhPlus(pepEle.attributeValue("MHplus"));
                    peptide.setCalcMHplus(pepEle.attributeValue("calcMHplus"));
                    peptide.setRedundancy(pepEle.attributeValue("redundancy"));
  //                  System.out.println("====" + locus + " " + pepEle.attributeValue("redundancy") + " " + pepEle.attributeValue("redundancy"));
                    peptide.setTotalIntensity( Double.parseDouble(pepEle.attributeValue("totalIntensity")));
                    //robin

                    String 		str = pepEle.attributeValue("rt");
                    if(null != str)
                        peptide.setRetentionTime( Double.parseDouble(str));

                    str = pepEle.attributeValue("iit");
                    if(null != str)
                        peptide.setIonInjectionTime(Double.parseDouble(str));

                    //peptide.setSpecCount(pepEle.attributeValue("redundancy"));

//                    String reteintionTime[] = pepEle.elementText("retention-time").split(";");

                    Iterator<String> indexItr =  ht.keySet().iterator();
                    while(indexItr.hasNext())
                    {
                        String currentIndexFile = indexItr.next();
                        if(! currentIndexFile.contains(pepEle.attributeValue("file")))
                        {
                            continue;
                        }



                        int[] keys = ht.get(currentIndexFile).getKeys();
                        int keyIndex = Arrays.binarySearch(keys,peptide.getScanNum());

                        if(keyIndex<0) //Cannot find index
                             keyIndex=-(++keyIndex); //Math.abs(++keyIndex);

                        if(keyIndex>=keys.length)
                            keyIndex--;

                        //peptide.setRetentionTime( ht.get(currentIndexFile).getRetentionTimeMap().get(keys[keyIndex]));
                      //  peptide.setIonInjectionTime(ht.get(currentIndexFile).getIonInjectionMap().get(keys[keyIndex]));
//                        sdafasdfsdfsasdf
//                                Still to generate the retention time for the specififc peptide.....

                    }
//                    int[] keys =new int[reteintionTime.length];
//                    for(int i=0;i<reteintionTime.length;i++)
//                    {
//                        String data[] =  reteintionTime[i].split(" ");
//                        keys[i] = Integer.parseInt(data[0]);
//                    }
//                    int keyIndex = Arrays.binarySearch(keys,peptide.getScanNum());
//
//                    if(keyIndex<0) //Cannot find index
//                        keyIndex=-(++keyIndex); //Math.abs(++keyIndex);
//
//                    if(keyIndex>=keys.length)
//                        keyIndex--;


//                    String ionInjectionTime[] = pepEle.elementText("ionInjection-time").split(";");
//                        for(int i=0;i<reteintionTime.length;i++)
//                    {
//                        String data[] =  ionInjectionTime[i].split(" ");
//                        String retentiondata[] =  reteintionTime[i].split(" ");
//                        if(keyIndex == Integer.parseInt(data[0]))
//                        {
//                            peptide.setIonInjectionTime(Double.parseDouble(data[1]));
//                            peptide.setRetentionTime(Double.parseDouble(retentiondata[1]));
//                            break;
//                        }
//
//                    }


		    String chroValue = pepEle.elementText("chro");

		    if(null != chroValue) {
			String[] dataArr = chroValue.split(";");
			String[] peakArr = dataArr[0].split(" ");
			peptide.setStartRange( peakArr[1] );
			peptide.setEndRange( peakArr[2] );

			int pepScan = peptide.getScanNum();
			int prevDataScan=-1;

			for(int i=1;i<dataArr.length;i++)
			{
			    String[] eachArr = dataArr[i].split(" ");

			    int dataScan = Integer.parseInt(eachArr[0]);

			    ChroData data;
			    if(eachArr.length>5)
			    {
				data = new ChroData(
					Integer.parseInt(eachArr[0]),
					(long)Double.parseDouble(eachArr[1]),
					(long)Double.parseDouble(eachArr[2]),
					Double.parseDouble(eachArr[7]),
					Double.parseDouble(eachArr[8]),
					Integer.parseInt(eachArr[3]), //total iso sam
					Integer.parseInt(eachArr[5]), //found iso sam
					Integer.parseInt(eachArr[4]), //total iso ref
					Integer.parseInt(eachArr[6]) //found iso ref
					);
			    }
			    else
			    {
				data = new ChroData(
					Integer.parseInt(eachArr[0]),
					(long)Double.parseDouble(eachArr[1]),
					(long)Double.parseDouble(eachArr[2])
					);
			    }



			    if(!"single".equals(conf.getPeakCount()))
				peptide.addData(data);
			    else {
				if(pepScan < dataScan) {
				    peptide.addData(data);

				    break;
				}
			    }

			    prevDataScan = dataScan;

			}

		    }

//		    edu.scripps.pms.census.util.AllNoneUtil.getANScore(peptide);

		    //if(("K.IHEQTDKETIEQGVK.D2").equals(peptide.getSequence() + peptide.getChargeState()) && peptide.getFileName().contains("3M")) {
		    //if(("K.IHEQTDKETIEQGVK.D2").equals(peptide.getSequence() + peptide.getChargeState())) {
		   // 	System.out.println("+++++++>>" +  peptide.getFileName() + "\t" + peptide.getScanNum() + "\t" + peptide.getSamIntensity());
		     // }

		    //protein.addPeptide(peptide);
		    //
		    mprotein.putPeptideHt(peptide, index);



		}
		//chroProteinList.add(protein);
		//System.out.println(proEle);
	    }
		index++;
        }


	gnu.trove.TDoubleArrayList[] proteinIntensityArr = new gnu.trove.TDoubleArrayList[expSize];
	gnu.trove.TDoubleArrayList[] peptideIntensityArr = new gnu.trove.TDoubleArrayList[expSize];
	double[] intSumArr = new double[expSize];
	for(int i=0;i<expSize;i++) {
	    proteinIntensityArr[i] = new gnu.trove.TDoubleArrayList();
	    peptideIntensityArr[i] = new gnu.trove.TDoubleArrayList();
	}

	for(Iterator<String> itr=proteinHt.keySet().iterator(); itr.hasNext(); )
	{
	    String htKey = itr.next();

	    MultipleChroProtein pro = proteinHt.get(htKey);

	    ChroProtein[] proteinArr = pro.getProteinArr();

	    p.print("P\t" + htKey + "\t");
	    p.print(pro.getDescription() + "\t");
            int count=0;
	    int maxSpecCount=-1;
	    int curSpecCount=-1;
	    int maxSpecCountIndex=-1;
	    for(ChroProtein cPro : proteinArr)
	    {
		if(null != cPro) {
		    p.print(cPro.getSpectrumCount() + "\t");
		    curSpecCount = cPro.getSpectrumCountAsInt();
		    if(curSpecCount>maxSpecCount) {
			maxSpecCount = curSpecCount;
			maxSpecCountIndex = count;
		    }

                }
		else
		    p.print("NA\t");
		count++;
	    }



/*
	    for(int i=0;i<proteinArr.length;i++) {
		if(null == proteinArr[i])
		{
		    p.print("0\t");
		    continue;
		}

		p.print( Formatter.formatThreeDecimal( (double)proteinArr[i].getSpectrumCountAsInt()/proteinArr[maxSpecCountIndex].getSpectrumCountAsInt()) );

		p.print("\t");
	    }*/


	    StringBuffer pepLineSb = new StringBuffer();
            Hashtable<String, ChroPeptide[]> peptideHt = pro.getPeptideHt();
	    gnu.trove.TDoubleArrayList[] pepRatioListArr = new gnu.trove.TDoubleArrayList[proteinArr.length];
	    double[] pepIntensiltySuml = new double[proteinArr.length];

	    for(int i=0;i<proteinArr.length;i++)
		pepRatioListArr[i] = new gnu.trove.TDoubleArrayList();

            double[] ratioArr = new double[peptideHt.size()];

            for(Iterator<String> itrPep=peptideHt.keySet().iterator(); itrPep.hasNext(); ) {
		String pepKey = itrPep.next();
                ChroPeptide[] pepArr = peptideHt.get(pepKey);

		pepLineSb.append("S\t");

                for(int i=0;i<pepArr.length;i++) {
                    if(pepArr[i] != null) {


/*
		    if(("K.IHEQTDKETIEQGVK.D2").equals(pepArr[i].getSequence() + pepArr[i].getChargeState())) {
		    	System.out.println("---+++++++>>");

		    	System.out.println("---+++++++>>" +  pepArr[i].getFileName() + "\t" + pepArr[i].getScanNum() + "\t" + pepArr[i].getMultiIntensitySum());
		      }*/

                                                /*
			if("single".equals(conf.getPeakCount())) {
			    pepLineSb.append("[" + (i+1) + "]\t"
                                    + pepArr[i].getSequence() + "\t"
                                    + pepArr[i].getFileName() + "\t"
                                    + pepArr[i].getScanNum() + "\t"
                                    + pepArr[i].getChargeState() + "\t"
                                    + pepArr[i].getMultiAveIntensity() + "\t"
                        //            + pepArr[i].getAverageIntensity() + "\t"
                                    + pepArr[i].getAnCompositeScore() + "\t");
		//	    if(!peptideNull) peptideIntensityArr[i].add( pepArr[i].getMultiAveIntensity() );
			}
			else {			    */

//System.out.println("------------------" + pepArr[i].getMultiIntensity() + " " + pepArr[i].getGaussianPeakArea());

double mergedPeakArea=pepArr[i].getGaussianPeakArea(); //when there are multiple peptides with same seq and cs, get average intensity
int multiplePepSize = pepArr[i].getMultiIntensity().size() + 1;
for(double d:pepArr[i].getMultiIntensity()){

 mergedPeakArea += d;
}
double averageIntensity = mergedPeakArea/multiplePepSize;

/*
double gmergedPeakArea=pepArr[i].getGaussianPeakAreaIonInjectionCorrection(); //when there are multiple peptides with same seq and cs, get average intensity
int gmultiplePepSize = pepArr[i].getGaussianPeakAreaIonInjectionCorrection().size() + 1;
for(double d:pepArr[i].getGaussianPeakAreaIonInjectionCorrection()){

 gmergedPeakArea+= d;
}
double gaverageIntensity = gmergedPeakArea/gmultiplePepSize;
*/



//System.out.println("====" + multiplePepSize + " " + mergedPeakArea + " " + mergedPeakArea/multiplePepSize + " " + multiplePepSize);

			    pepLineSb.append("[").append((i+1)).append("]\t")
                                    .append(pepArr[i].getSequence()).append("\t")
                                    .append(pepArr[i].getFileName()).append("\t")
                                    .append(pepArr[i].getScanNum()).append("\t")
                                    .append(pepArr[i].getChargeState()).append("\t")
                                    //.append(pepArr[i].getMultiAveIntensity()).append("\t")
                                    //.append(pepArr[i].getGaussianPeakArea()).append("\t")
                                    .append(averageIntensity).append("\t")
                                    .append(pepArr[i].getGaussianPeakAreaIonInjectionCorrection()).append("\t")
                                    .append(pepArr[i].getAnCompositeScore()).append("\t")
                                    .append(pepArr[i].getMhPlus()).append("\t")
                                    .append(pepArr[i].getCalcMHplus()).append("\t")
                                    .append(pepArr[i].getTotalIntensity()).append("\t")
                                    .append(pepArr[i].getXCorr()).append("\t")
                                    .append(pepArr[i].getDeltCN()).append("\t")
                                    .append(pepArr[i].getDeltMass()).append("\t")
                                    .append(pepArr[i].getSpRank()).append("\t")
                                    .append(pepArr[i].getSpScore()).append("\t")

                                    .append(pepArr[i].getRedundancy()).append("\t")
                                    //.append(pepArr[i].getSpecCount()).append("\t")
                                    .append(pepArr[i].getStartRange()).append("\t")
                                    .append(pepArr[i].getEndRange()).append("\t")
                                    .append(pepArr[i].getRetentionTime()).append("\t")
                                    .append(pepArr[i].getIonInjectionTime()).append("\t");

                      //      System.out.println("====--===" + pepArr[i].getSpScore());
           //             }
		    }
                    else {
			pepLineSb.append("[").append((i+1)).append("]\t");

			for(int k=0;k<20;k++)
				pepLineSb.append("NA\t");

		    }

//NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t");
                }

		int intensityMaxIndex=-1;
		double intensityCurrentValue=0;
		double intensityMaxValue=-1;
		double totalIntensity=0;

                for(int i=0;i<pepArr.length;i++) {

		    if(null == pepArr[i])
			continue;

		    //intensityCurrentValue = pepArr[i].getMultiIntensitySum();
		    intensityCurrentValue = pepArr[i].getMultiAveIntensity();
		    if(intensityCurrentValue>intensityMaxValue) {
			intensityMaxIndex = i;
			intensityMaxValue = intensityCurrentValue;
		    }

		    totalIntensity += intensityCurrentValue;
		}

		double averageIntensity = totalIntensity/expSize;

//		System.out.println(intensityMaxIndex  +"\t" + pepArr[intensityMaxIndex].getSamIntensity() + "\t" + pepArr[intensityMaxIndex].getMultiAveIntensity());

//if(!"IPI00469221.5".equals(locus))
//continue;



                for(int i=0;i<pepArr.length;i++) {
	//	    if(!peptideNull) pepIntensitySum[i] += pepArr[i].getSamIntensity();

		    //if(!"single".equals(conf.getPeakCount())) {
		    /*
			if(pepArr[i] != null)
			{
			    //double pepRatio = pepArr[i].getMultiAveIntensity()/averageIntensity;
			    double pepRatio = pepArr[i].getMultiIntensitySum()/averageIntensity;
			    pepLineSb.append(pepRatio + "\t");
			    pepRatioListArr[i].add(pepRatio);
			}
			else {
			    pepLineSb.append("0\t");
			    pepRatioListArr[i].add(0);
			}
			*/
		    //}
		    //else
		    {
			if(null !=pepArr[i])
			{
			    //double pepRatio = pepArr[i].getMultiAveIntensity()/averageIntensity;
			    //double pepRatio = Math.log10(pepArr[i].getMultiAveIntensity())/Math.log10(averageIntensity);
			    double pepRatio = pepArr[i].getMultiAveIntensity()/averageIntensity;
	//		    pepIntensitySum[i] += pepArr[i].getMultiAveIntensity();

//			    System.out.println("----" +  pepArr[i].getMultiAveIntensity() + " " + averageIntensity + " " + pepRatio + "\t");
			    pepLineSb.append(pepRatio + "\t");
			    pepRatioListArr[i].add(pepRatio);
			}
			else {
			    pepLineSb.append("0\t");
			    pepRatioListArr[i].add(0);
			}
		    }
                }

		pepLineSb.append("\n");
	    }

	    //p.print(pepRatioListArr[0].toNativeArray().length);
	    p.print(peptideHt.keySet().size()+"\t");
	    //for NORM_INTENSITY later value will be added here

            for(ChroProtein cPro : proteinArr)
	    {
                p.print("NA\t");
            }
	    p.print("\n");
	    p.print(pepLineSb.toString());
	}



	if(null != p)
	    p.close();

	if(null != out)
	    out.close();



	int minIdPerReplicates = 0;
        org.jdom.Element paramEle = conf.getRootConfEle().getChild("params");
	org.jdom.Element minIdEle = paramEle.getChild("min_id_in_replicates");
	if(null != minIdEle)
	    minIdPerReplicates = Integer.parseInt(minIdEle.getText());

	scripts.labelfree.NormalizeLabelFreeIntensityBase norm = new scripts.labelfree.NormalizeLabelFreeIntensityBase();

	System.out.println("normalize...");
	double[] correctArr = norm.getNormalizeFactor(outputFilenameTmp); //"census_label_free_result_with_id_temp.txt");
	norm.normalize(outputFilename, outputFilenameTmp, correctArr, minIdPerReplicates);

	System.out.println("Done.  Please see census_label_free_result_with_id.txt");
	System.out.println("Done.  Please see " + outputFilename);
	/*
	MultipleChroProtein pro = proteinHt.get("contaminant_KERATIN21");

	ChroPeptide[] arr = pro.getPeptideHt().get("K.AQYEEIAQR.S2");

	for(ChroPeptide cp:arr)
	    System.out.println(cp);
	*/

    }
}
