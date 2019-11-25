/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package edu.scripps.pms.census.labelFree;

/**
 *
 * @author Harshil
 */
import java.util.*;
import java.io.*;
import org.jdom.input.SAXBuilder;
import org.jdom.*;
//import com.ipa.ip2.util.Formatter;
//
//import com.ipa.ip2.model.db.Quantitation;
//import com.ipa.ip2.model.Replicate;
//import com.ipa.ip2.model.Replicates;

import org.apache.commons.math.distribution.*;
import org.apache.commons.math.stat.inference.*;
public class LabelfreeUtil {
    


    private List<LabelfreeQuantitationModel>  quantModelList = new ArrayList<LabelfreeQuantitationModel>();
    private List<Replicates>  sampleList = new ArrayList<Replicates>();
       

    public LabelfreeUtil(String path, String id) throws Exception {
	readLabelfreeResult(path,id);
	

    }

public static void printUsage()
    {
            System.out.println("Usage: LabelfreeUtil /home/harshil/data/jolene/labelfree_quant 7841");
            System.out.println("[1] : path where the config file is stored ");
            System.out.println("[2] : id");
    }
    
    public static void main(String[] args) throws Exception {
//        Quantitation q = new Quantitation();
//        q.setPath("/ip2_data/hseol/Hutu_secretomes/labelfree_quant");
//        q.setId(712);
        if(args.length <2 )
        {
            printUsage();
            return;
        }
//        String path = "/home/harshil/data/jolene/labelfree_quant";
        String path = args[0];
        
        LabelfreeUtil lf = new LabelfreeUtil(path,args[1]);
	//System.out.println("-------------------------");// + lf.getQuantModelList());
        //System.out.println("-------------------------"+ lf.getQuantModelList().size());

//	lf.readLabelfreeResult(q);
    }


 public void readLabelfreeResult(String path,String id) throws IOException, Exception {
 		ChiSquaredDistribution dist = new ChiSquaredDistributionImpl(1);
    String confFilename = path + File.separator + "census_config_labelfree_" + id + ".xml";
//        String confFilename = path;
        
    Document doc = new SAXBuilder().build(new File(confFilename));
    Element root = doc.getRootElement();

//        int expSize = root.getChildren("sample").size();

    String filename = path + File.separator + "census_labelfree_out_" + id + "_stat.txt";

    BufferedReader br = null;

		try {
				br = new BufferedReader(new FileReader(filename));

	
				String eachLine;
	
				List<String> sNameList = new ArrayList<String>();
				List<Integer> sSizeList = new ArrayList<Integer>();
	
				int expSize=0;
	
				while( (eachLine = br.readLine()) != null ) {
					if(eachLine.startsWith("PLINE\tLOCUS"))
						break;
	
					if(eachLine.startsWith("H\tGROUP_SAMPLE") || eachLine.startsWith("H\tsample_group")) {
						String[] arr = eachLine.split("\t");
	
						if(arr.length<3)
							break;
		
						Replicates repList = new Replicates();
						repList.setName(arr[2]);
		
						for(int i=3;i<arr.length;i++) {
							repList.addReplicate(new Replicate(arr[i]));
						}
		
						this.sampleList.add(repList);
		
		
						sNameList.add(arr[2]);
						sSizeList.add(arr.length-3);
						expSize += arr.length-3;
	
					}//end if
				}//end while

				String[] tmpArr = eachLine.split("\t");


				int accessionIndex=-1;
				int descriptionIndex=-1;
				int normCorrectAnovaFvalueIndex = -1;
				int normCorrectAnovaPvalueIndex = -1;


				//        System.out.println(sNameList);
				//System.out.println(sSizeList);
		
				int sampleSize = sNameList.size();
				int[] specCountIndexArr = new int[sampleSize];
				int[] normSpecCountIndexArr = new int[sampleSize];
				int[] normIntensityIndexArr = new int[sampleSize];
				int[] normCorrectIntensityIndexArr = new int[sampleSize];
				int[] normCorrectIntensityAvgIndexArr = new int[sampleSize];
				int[] ratioPvalueIndexArr = new int[sampleSize-1];
				int[] logRatioChangeIndexArr = new int[sampleSize-1];
			
				double[] totalSCAvg=new double[sampleSize];
				double[] correctSC=new double[sampleSize];//corrected Scpect Count
				double[] eachSpecCountAvg=new double[sampleSize];
				   double[] eachNormSpecCountAvg=new double[sampleSize];
			
		
				int tmpScountIndex=0;
				int tmpNormScountIndex=0;
				int tmpNormIntensityIndex=0;
				int tmpNormCorrectIntensityIndex=0;
				int ratioPvalueIndex=0;
				int tmpLogRatioChangeIndex=0;
				int tmpNormCorrectIntensityAvgIndex=0;
                                int sortedCodeIndex = -1;
                                int sortedCodeSCFillerIndex = -1;

				for(int i=0;i<tmpArr.length;i++) {
					if(tmpArr[i].startsWith("LOCUS")) {
						accessionIndex = i;
					} else if(tmpArr[i].startsWith("NORM_CORRECT_ANOVA_FVALUE")) {
						normCorrectAnovaFvalueIndex = i;
					} else if(tmpArr[i].startsWith("NORM_CORRECT_ANOVA_PVALUE")) {
						normCorrectAnovaPvalueIndex = i;
					} else if(tmpArr[i].startsWith("DESCRIPTION")) {
						descriptionIndex = i;
					} else if(tmpArr[i].startsWith("SPEC_COUNT")) {
						specCountIndexArr[tmpScountIndex++] = i;
					} else if(tmpArr[i].startsWith("NORM_SPEC_COUNT")) {
						normSpecCountIndexArr[tmpNormScountIndex++] = i;
					} else if(tmpArr[i].startsWith("NORM_CORRECT_INTENSITY")) {
						normCorrectIntensityIndexArr[tmpNormCorrectIntensityIndex++] = i;
					} else if(tmpArr[i].startsWith("NORM_INTENSITY")) {
						normIntensityIndexArr[tmpNormIntensityIndex++] = i;
					} else if(tmpArr[i].startsWith("LOG_RATIO_CHANGE")) {
						logRatioChangeIndexArr[tmpLogRatioChangeIndex++] = i;
					} else if(tmpArr[i].startsWith("RATIO_PVALUE")) {
						ratioPvalueIndexArr[ratioPvalueIndex++] = i;
					} else if(tmpArr[i].startsWith("NORM_CORRECT_AVERAGE_RATIO")) {
						normCorrectIntensityAvgIndexArr[tmpNormCorrectIntensityAvgIndex++] = i;
					}
                                        else if(tmpArr[i].startsWith("SORTED_CODE")) {
						sortedCodeIndex = i;
					}
                                        else if(tmpArr[i].startsWith("SORTED_SCFILLERED_COD")) {
						sortedCodeSCFillerIndex = i;
					}     
				}//end for
				LabelfreeQuantitationModel protein ;
				while( (eachLine = br.readLine()) != null) {

						if(!eachLine.startsWith("P\t"))
							continue;
			
						String[] arr = eachLine.split("\t");
						protein= new LabelfreeQuantitationModel();
						protein.setAccession(arr[accessionIndex]);
						protein.setDescription(arr[descriptionIndex]);
						if("NA".equals(arr[normCorrectAnovaFvalueIndex]))
							protein.setAnovaFvalue(1000000.0);
						else {
							try { protein.setAnovaFvalue(Double.parseDouble(arr[normCorrectAnovaFvalueIndex])); } 
							catch(Exception pe) { protein.setAnovaFvalue(1000000.0); }
						}
			
						if("NA".equals(arr[normCorrectAnovaPvalueIndex]))
							protein.setAnovaPvalue(1000000.0);
						else {
							try { protein.setAnovaPvalue(Double.parseDouble(arr[normCorrectAnovaPvalueIndex])); }
							catch(Exception pe) { protein.setAnovaPvalue(1000000.0); }
						}
			
						protein.setAnovaPvalueLog10(-Math.log10(protein.getAnovaPvalue()));
						Double sum=0.0;
						int jj=0;
				  
						for(int in:specCountIndexArr) {				
							protein.addSpecCount(arr[in]);
							String[] tmp=arr[in].split(",");					
							sum=0.0;
							for(int i=0;i<tmp.length;i++){						
								sum+=Double.parseDouble(tmp[i]);	
							}				
							if(tmp.length!=0)
											eachSpecCountAvg[jj]=sum/tmp.length;
							else eachSpecCountAvg[jj]=0;
			
							protein.addSpecCountAvgList(eachSpecCountAvg[jj]);
							totalSCAvg[jj]+=	eachSpecCountAvg[jj];
							jj++;	    
					  }//end for
	//	for(int in:normSpecCountIndexArr) {
		//	protein.addNormSpecCount(arr[in]);
		//}
				  	double[] d_tmp1 =new double[arr[normSpecCountIndexArr[0]].split(",").length];
					double[] d_tmp2 =new double[arr[normSpecCountIndexArr[1]].split(",").length];
					jj=0;
		
					for(int in:normSpecCountIndexArr) {
						//protein.addNormSpecCount(arr[in]);
						sum=0.0;
						String[] tmp=arr[in].split(",");					
						String tmpOut="";
						for(int i=0;i<tmp.length;i++){						
							if(jj==0){
								d_tmp1[i]=Double.parseDouble(tmp[i]);		
							}
							else if(jj==1){
								d_tmp2[i]=Double.parseDouble(tmp[i]);		 
							}
							tmpOut+=""+Formatter.sciRound(Double.parseDouble(tmp[i]))+",";
							
							sum+=Double.parseDouble(tmp[i]);	
						}	
				
						protein.addNormSpecCountList(tmpOut);

						if(tmp.length!=0)
											eachNormSpecCountAvg[jj]=sum/tmp.length;
						else eachNormSpecCountAvg[jj]=0;
		
						protein.addNormSpecCountAvgList(eachNormSpecCountAvg[jj]);
							
					//   System.out.println(arr[1]+"===="+tmpOut+"\t==="+eachNormSpecCountAvg[jj]);
						jj++;	
					}//end for		
/*
		  //calculating NormSpecCount log 2 ratio change	
						for(int i=1;i<sampleSize;i++){			
					//System.out.println(protein.getNormSpecCountAvgList().get(0));
								 if(protein.getNormSpecCountAvgList().get(0)==0)		
									}
						}//end for	
	*/		
			
					  //calculate only first two norm scpectCount for t-test	
			
						if((d_tmp1.length==1)||(d_tmp2.length==1)){
							protein.setSpecCountTscore(100000.0);
						}
			
						else{
							List classes=new ArrayList();
							classes.add(d_tmp1);
							classes.add(d_tmp2);
	
							/*try{			protein.setSpecCountTscore(Formatter.round(TestUtils.pairedTTest(d_tmp1, d_tmp2),3));		  
							}catch(Exception e){
							
							protein.setSpecCountTscore(Formatter.round(TestUtils.oneWayAnovaPValue(classes),3));
							System.out.println("else"+protein.getSpecCountTscore());
							}
								}
							System.out.println("==33333" + eachLine);
							*/
							try{
							protein.setSpecCountTscore(Formatter.round(TestUtils.oneWayAnovaPValue(classes),3));
							}catch(Exception e){
								System.out.println("error"+e);
							}
			     		  	}		

						for (int in:normIntensityIndexArr){
							protein.addNormIntensityList(arr[in]);
			  			}		
						for(int in:normCorrectIntensityAvgIndexArr){
							if("NA".equals(arr[in]))
								protein.addCorrectedIntensityAverageList(0);
							else
								protein.addCorrectedIntensityAverageList(Double.parseDouble(arr[in]));
						}
						/*		for(int i=1;i<sampleSize;i++){
									double x=0.0;
									if(protein.getCorrectedIntensityAverageList().get(0)==0){
											System.out.println("=====NC");
									}	
									else if(protein.getCorrectedIntensityAverageList().get(i)==0)		
							System.out.println("=====NC");
									else {
										x=protein.getCorrectedIntensityAverageList().get(0)/protein.getCorrectedIntensityAverageList().get(i);
						System.out.println("====="+(Math.log(x)/Math.log(2)));
									}
								}
								*/
						for(int in:logRatioChangeIndexArr) {
							if("NA".equals(arr[in]))
								protein.addLogRatioChangeList(1000000);
								
							else {
								try { protein.addLogRatioChangeList(Double.parseDouble(arr[in])); } catch (Exception e) {protein.addLogRatioChangeList(1000000);}
							}
						}

						for(int in:ratioPvalueIndexArr) {
							if("NA".equals(arr[in]))
								protein.addRatioPvalueList(1000000);
							else {
								try { protein.addRatioPvalueList(Double.parseDouble(arr[in])); } catch (Exception e) {protein.addRatioPvalueList(1000000);}
							}
						}
						for(int in:normCorrectIntensityIndexArr) {
							protein.addNormIntensityCorrectList(arr[in]);
						}
                                                
                                                if(sortedCodeIndex != -1)
                                                    protein.setSortCode(arr[sortedCodeIndex]);
                                                if(sortedCodeSCFillerIndex != -1)
                                                    protein.setSortSCFillered(arr[sortedCodeSCFillerIndex]);
						quantModelList.add(protein);
				}//end while
		
				//to raise up Norm spect count values
				double ratio=0;
				for(double aTotalSC:totalSCAvg){
					ratio+=aTotalSC;					
				}
				ratio=ratio/sampleSize;
						
			  //calculate Gscore and Pvalue//
			  // calculate only first two spec count score//
		    double avg_total=(totalSCAvg[0]+totalSCAvg[1])/2;
			boolean anovaExceptionOccurred = false;
		  	for(Iterator<LabelfreeQuantitationModel> itr=this.quantModelList.iterator(); itr.hasNext(); ) {
		
		         LabelfreeQuantitationModel qmodel = itr.next();
		               
						 correctSC[0]=qmodel.getSpecCountAvgList().get(0);//C
						 correctSC[1]=qmodel.getSpecCountAvgList().get(1);//D
		         double avg=(correctSC[0]+ correctSC[1])/2;
		         int sign=0;
		         if((avg-correctSC[0])>0)sign=1;
		         else sign=-1;                 
		         double E=correctSC[0]+0.5*sign;//E
		                 
		         if((avg-correctSC[1])>0)sign=1;
		         else sign=-1;
		         double F=correctSC[1]+0.5*sign;//F
		                
		          double normSC1= E*avg_total/totalSCAvg[0];   //G
		          double normSC2= F*avg_total/totalSCAvg[1];  //H
		         double avg_norm=(normSC1+normSC2)/2;//I or J
		    	   double ratioSC1=normSC1/avg_norm;//K
		      	 double ratioSC2=normSC2/avg_norm;//L
		         double gscore=(normSC1*Math.log(ratioSC1)+normSC2*Math.log(ratioSC2))*2;
		         double pvalue = 1-dist.cumulativeProbability(gscore);       
		              		
		         qmodel.setGscore(Formatter.round(gscore , 3));
		         qmodel.setGscorePvalue(Formatter.round(pvalue , 3)); 
		  
		        //to raise value of norm spec counter
// System.out.println("--------------"+ratio+"\t"+sampleSize);
		   
			
					  for(int i=0;i<sampleSize;i++){
					  	    try{   
								//protein.addNormSpecCount(arr[in]);
								double sum=0.0;
								String[] tmp=qmodel.getNormSpecCountList().get(i).split(",");					
//	System.out.print("old="+qmodel.getNormSpecCountList().get(i)+"\t");
							 String tmpOut="";
								for(int j=0;j<tmp.length;j++){
							
									double dt=Double.parseDouble(tmp[j])*ratio;		
									tmpOut+=""+Formatter.sciRound(dt)+",";
									
									sum+=Double.parseDouble(tmp[j]);	
								}							
								qmodel.editNormSpecCountList(i,tmpOut);
										
							if(tmp.length!=0)
								qmodel.editNormSpecCountAvgList(i,sum*ratio/tmp.length);
							else 
								qmodel.editNormSpecCountAvgList(i,0);
			
	//System.out.println("");
	//System.out.println("Avg new= i="+i+"="+qmodel.getNormSpecCountAvgList().get(i));
		//System.out.println("new= i="+i+"="+qmodel.getNormSpecCountList().get(i));
		
				    				//calculating NormSpecCount log 2 ratio change		
				    				if(i==0) continue;
				    				else {
										if(qmodel.getNormSpecCountAvgList().get(0)==0)	{
											qmodel.addNormSpecCountlogRatioChangeList(-9999.999);
											qmodel.addNormSpecCountRatioChangeList(-9999.999);
										}
										else if(qmodel.getNormSpecCountAvgList().get(i)==0)		{
											qmodel.addNormSpecCountlogRatioChangeList(-9999.999);
											qmodel.addNormSpecCountRatioChangeList(-9999.999);
										}
										else {
											double x=qmodel.getNormSpecCountAvgList().get(0)/qmodel.getNormSpecCountAvgList().get(i);
											qmodel.addNormSpecCountlogRatioChangeList(Formatter.round((Math.log(x)/Math.log(2)),3));
											qmodel.addNormSpecCountRatioChangeList(Formatter.round(x,3));
										}
										}
									}catch (Exception e) {
									e.printStackTrace();
									}
							
						}//end for	
		        
//						if(qmodel.getAccession().startsWith("T11B7")){
							int num = qmodel.getNormSpecCountList().size();
							int len = 1;
							List classes = new ArrayList();
							for(int x=0; x<num; x++){
								String scList = qmodel.getNormSpecCountList().get(x);
								String[] ss = scList.split(",");
								len = (len<ss.length)?ss.length:len;
								double[] dd = new double[ss.length];
								for(int y = 0; y<ss.length; y++)
									dd[y] = Double.parseDouble(ss[y]);
								classes.add(dd);
							}
							try{
								qmodel.setSpecCountTscore(Formatter.round(TestUtils.oneWayAnovaPValue(classes),3));
							}catch(Exception e){
								if(!anovaExceptionOccurred){
									System.out.println("error for anova calculation: "+e);
									anovaExceptionOccurred=true;
								}
							}
		   	} //end iterator for              
		
		} catch (Exception e) {
					try {   if(null != br) br.close(); }
					catch(IOException ie) { }
		}
        //System.out.println("============" + l);
 }//end method readLabelfreeResult

    public void setQuantModelList(List<LabelfreeQuantitationModel> quantModelList) {
	this.quantModelList = quantModelList;
    }

    public List<LabelfreeQuantitationModel> getQuantModelList() {
	return quantModelList;
    }

    public void setSampleList(List<Replicates> sampleList) {
	this.sampleList = sampleList;
    }

    public List<Replicates> getSampleList() {
	return sampleList;
    }

}
