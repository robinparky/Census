 /*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.scripps.pms.census.tmtFilter;

import edu.scripps.pms.stats.GrubbsTest;
import java.io.*;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
//import sun.security.krb5.internal.crypto.Des;

import static edu.scripps.pms.util.isotope.RatioIntensityUtils.checkRatio;


 /**
 *
 * @author Harshil
 * @author 2018 Titus Jung
 * Its a generic version at peptide Level Analysis.................
 * There is a Main HashMap "peptideData" Contains the Key(Peptide ) to the Value(object of "Peptide" type which contains the data related to the peptide(Ex. the whole row).)
 *  
 * To run this file type 

java peptide_process.pepseq_v3 /data/6/harshil/peptide_process/census-out.txt

In CMD so java classfile pathofthecensus directory......
 * 
 */
public class TmtfilterPeptide extends AbstractTmt
{
//    Map<String, List> groupNameMap = new LinkedHashMap<>();
//    Map<String, List<String>> groupKeyMap = new LinkedHashMap<>();//GruoupName -> List<Key(header_Line_name)>
//    List<String> groupNameList = new ArrayList<>();
//    public static ArrayList<String> intensityNames = new ArrayList<>();
    String words[] = new String[10000];
//    String headerSline[] = new String[150];
    ArrayList<String> peptideNameList = new ArrayList<>();
    ArrayList<String> normintensity_name = new ArrayList<>();
    Map<String, Peptide> peptideData = new HashMap<>(); // The main map for contains the peptide sequence name and all the datails of each peptide sequence.
    int count_intensity = 0;
    boolean checkzero = false;
    int position_spc = 0, position_scn = 0, position_cstate = 0, position_file = 0;
        
    public boolean checkValidEntryForpeptide(Peptide pep_detail,int newScanNumber,String newFileName)
      {
          ArrayList scanNumber = pep_detail.getScannum();
          Set fileName =  pep_detail.getFilename();
          if(scanNumber.contains(newScanNumber) && fileName.contains(newFileName))
              return false;
          else
              return true;
      }
     
      public void parseCensusFile(String inputPath){
        BufferedReader br = null;

       
	System.out.println("Please wait...Parsing the file.......");
        try
        {
//          File f = new File(path + File.separator + "census-out-small.txt");
            File f = new File(inputPath + File.separator + "census-out.txt");
//            File f = new File(path + File.separator + "census-out-navin.txt");
            
            br = new BufferedReader(new FileReader(f));
            String line;
            line=br.readLine();
            Peptide pep_details = new Peptide();        
           while(line != null)
            {
                checkzero=false;
                pep_details=new Peptide();
                if(line.startsWith("H") || line.startsWith("h")) 
                {
                   // System.out.println("Its a header line.");
                    setHeaderLine(line);
                }
                if(line.contains("PLINE"))
                    header_line= line.split("\t");
                if(line.startsWith("P") || line.startsWith("p"))
                    words=line.split("[\t]");
                if(line.startsWith("S") || line.startsWith("s"))
                 {
//                     parsePeptideLine(line);
                     //System.out.print("Its a Peptide line.");
                     String[] words_peptide=line.split("[\t]");
                     //below logic is to check single intensity is 0 or not if its zero then skip it..
//                     for(int temp=0;temp<count_intensity;temp++)
//                     {
//                          if(Float.parseFloat(words_peptide[temp+4-1]) == 0)
//                            {
//                                checkzero=true;
//                                break;
//                            }
//                     }
//                     if(checkzero)
//                         {
//                            line=br.readLine();
//                            continue;
//                         }   
                     if(peptideData.containsKey(words_peptide[2]))
                     {
                         pep_details=peptideData.get(words_peptide[2]);
                         pep_details.setPeptideCounter();
                         if(!checkValidEntryForpeptide(pep_details,Integer.parseInt(words_peptide[position_scn]),words_peptide[position_file]))
                         {
                             line=br.readLine();
                             continue;
                         }
                         for(int temp=0;temp<count_intensity;temp++)
                        {
                            double value = Double.parseDouble(words_peptide[temp+4-1]);
                            if(value<=0)//if intensity is negative then make it zero..
                                value = 1;
                            double value2 = Double.parseDouble(words_peptide[temp+4]);
                            if(value2<=0)//if intensity is negative then make it zero..
                                value2 = 1;
                            pep_details.setallIntensity(headerSline[temp+4], value);
                            pep_details.setAllNormIntensity(headerSline[temp+4],value2);


                        }
                         setPeptideDetails(pep_details, words_peptide);
                     }
                     else
                     {
//                         for (int temp = 0; temp < count_intensity; temp++) {
//                             if (Float.parseFloat(words_peptide[temp + 4 - 1]) == 0) {
//                                 checkzero = true;
//                                 break;
//                             }
//                         }
//                         if (checkzero)
//                         {
//                            line=br.readLine();
//                            continue;
//                         }   
                         if(!peptideNameList.contains(words_peptide[2]))
                            peptideNameList.add(words_peptide[2]);
//                         else
//                           System.out.println("Error in addinf peptide at the peptide list");
                         for(int temp=0;temp<count_intensity;temp++)
                        {
                            double value = Double.parseDouble(words_peptide[temp+4-1]) ;
                            double value2 = Double.parseDouble(words_peptide[temp+4]);
                            if(value<1)//if intensity is negative then make it zero..
                                value = 1;
                            pep_details.setallIntensity(headerSline[temp+4], value);
                            pep_details.setAllNormIntensity(headerSline[temp+4],value2);
                         }
                         pep_details.setPeptideCounter();
                          setPeptideDetails(pep_details, words_peptide);
                          pep_details.setSequence(words_peptide[2]);
                     }
                         pep_details.addProteinLocus(words[1]);
                        String proteinDescription = words[1]+" "+words[words.length-1];
                         pep_details.addProteinDescription(proteinDescription);
                         String geneName = TmtfilterPeptideReader.extractGeneName(proteinDescription);
                        pep_details.addGeneName(geneName);
                         peptideData.put(words_peptide[2],pep_details);
//                         hash.add_data(words_peptide[2], words[1]);
                 }
                line=br.readLine();
            }
          //  findSum(peptideData);
//            setRatio(peptideData, normintensity_name);
//           System.out.println();
           // hash.display(peptideData,ar,path);
        }
        catch (FileNotFoundException ex)
		 {        
            Logger.getLogger(TmtfilterPeptide.class.getName()).log(Level.SEVERE, null, ex);
       	 }
		catch (IOException ex)
		 {
            Logger.getLogger(TmtfilterPeptide.class.getName()).log(Level.SEVERE, null, ex);
      		  }        finally
        {
           
          }
    }

    public void setPeptideDetails(Peptide pep_details,String[] words_peptide) {
        pep_details.setSpc(Integer.parseInt(words_peptide[position_spc]));
        pep_details.addScannum(Integer.parseInt(words_peptide[position_scn]));
        pep_details.setCstate(Integer.parseInt(words_peptide[position_cstate]));
     //   pep_details.setratio(normintensity_name);
        pep_details.setFilename(words_peptide[position_file]);
    }
    
    public void setHeaderLine(String line) 
    {
        headerSline = line.split("[\t]");
        if (headerSline[1].equals("SLINE"))
        {
            for (int temp = 4; temp < headerSline.length; temp++) 
            {
                String splitter[] = new String[10];
                splitter = headerSline[temp].split("[_()]");
                if (splitter[0].equals("norm") && !splitter[1].equals("ratio"))
                    normintensity_name.add(headerSline[temp]);
                if (headerSline[temp].contains("_") && !headerSline[temp].contains("norm")) 
                    intensityNames.add(headerSline[temp]);
                if (headerSline[temp].equals("SpC")) 
                {
                    count_intensity = temp - 4; // count total number of intensity.
                    position_spc = temp - 1;
                    position_scn = position_spc + 1;
                    position_cstate = position_scn + 1;
                    position_file = position_cstate + 1;
                    break;
                }
//                ar.add(headerSline[temp]);
            }
        }
    }
    
//    public void groupNameToKeys(Map groupNameMap)
//    {
//        if(!headerSline[1].equalsIgnoreCase("SLINE"))
//            return;
//        List keys = new ArrayList();
//        keys.addAll(groupNameMap.keySet());
//        for(int z =0;z<groupNameMap.size();z++)
//        {
//           List ar1= (ArrayList<String>) groupNameMap.get(keys.get(z));
//           List keyArray = new ArrayList();
//            for(int j=0;j<ar1.size();j++)
//            {
//                
//                for(int i =0;i<headerSline.length;i++)
//                {
//                     if(headerSline[i].toLowerCase().contains("norm") 
//                             && headerSline[i].toLowerCase().contains("m/z") 
//                             &&  headerSline[i].toLowerCase().contains((String)ar1.get(j)))
//                            keyArray.add(headerSline[i]);
//                }
//            }
//            groupKeyMap.put((String) keys.get(z), keyArray);
//        }
//        System.out.println("Group MAped successfully..............");
//            
//    }
    
    
    
    
    
    public void generateGroupInfo()
    {
        //first find Sum its done in the find sum method////
        Iterator itr = peptideData.keySet().iterator();
        List groupNameList = new ArrayList(groupKeyMap.keySet());
        while(itr.hasNext())
        {
            String currentPeptideName = (String) itr.next();
            Peptide currentPeptide = peptideData.get(currentPeptideName);
//            if(currentPeptide.getSequence().equals("K.TAVNALWGK.V"))
//                System.out.println("");
            Map<String, List> allIntensity = currentPeptide.getallIntensity();
            Map<String,List<Double>> normAllIntensity = currentPeptide.getNormAllIntensity();
            for(int z=0;z<allIntensity.get(normintensity_name.get(0)).size();z++)
            {
                List statList = new ArrayList();
                List normStatList = new ArrayList();
                for(int i=0;i<groupNameList.size()-1;i++)
                {
                    List currrentGroupKey = groupKeyMap.get(groupNameList.get(i));
                    for(int j=i+1; j<groupNameList.size();j++)
                    {
                        DescriptiveStatistics stat = new DescriptiveStatistics();
                        DescriptiveStatistics normStat = new DescriptiveStatistics();
                          List nextGroupKey = groupKeyMap.get(groupNameList.get(j));
                          for(int t1 = 0;t1<currrentGroupKey.size();t1++)
                          {
                              for(int t2 =0;t2<nextGroupKey.size();t2++)
                              {
                                  double value1 = (double)((List)allIntensity.get(currrentGroupKey.get(t1))).get(z);
                                  double value2=  (double)((List)allIntensity.get(nextGroupKey.get(t2))).get(z);
                                  double ratio = value1/value2;

                                  double nValue1 = (normAllIntensity.get(currrentGroupKey.get(t1)).get(z));
                                  double nValue2 = (normAllIntensity.get(nextGroupKey.get(t2)).get(z));
                                  double normRatio = nValue1/nValue2;

                                  ratio = checkRatio(ratio);

                                  if(!Double.isNaN(ratio))
                                  {
                                      double logRatio = Math.log(ratio) / Math.log(2);
                                      stat.addValue(logRatio);
                                  }


                                  normRatio = checkRatio(normRatio);

                                  if(!Double.isNaN(normRatio))
                                  {
                                      double normLogRatio = Math.log(normRatio)/ Math.log(2);
                                      normStat.addValue(normLogRatio);
                                  }

                              }
                          }

                            statList.add(stat);
                        normStatList.add(normStat);
                    }
                }
                currentPeptide.comparedGroup.add(statList);
                currentPeptide.normComparedGroup.add(normStatList);
            }
            // Add the grouping to the stat object from which we can have the detaisl ...........
        }
//        System.out.println("Generating the descriptive stats done..........");

    }
    public void filterValues(double pValue)
    {//This will generate filterValues and and remove all those which are below t
        // specified pValue passed in the argument...
        
        System.out.println("Filtering values......");
//        GrubbsTest test = new GrubbsTest();
        List<Peptide> peptideDataList= new ArrayList();
        Iterator itr= peptideNameList.iterator();
        int counter =0;
        while(itr.hasNext())
        {
            counter=0;
            String currentSequence = (String) itr.next();
//            System.out.println(currentSequence);
//            if(currentSequence.equals("K.C(272.205135)IIVEEGKEILVGDVGVTITDPFK.H"))
//                System.out.println("");
             counter=GrubbsTest.filterAndRemovePeptideNoLog(peptideData.get(currentSequence),pValue);
//             if(peptideData.get(currentSequence) == null)
//                 System.out.println("NUll peptide found.........");
//             if(counter!=0)
//                System.out.println("-------"+counter);
        }
//        System.out.println("Total peptide removed: " + counter);
    }
    public void generateRatio()
    {
        Iterator itr= peptideNameList.iterator();
        boolean status;
        List<String> removingPeptide = new ArrayList();
        List<List<String>> headerName = new ArrayList<>();
        for(String currentGroupKey : groupKeyMap.keySet())
        {
            headerName.add(groupKeyMap.get(currentGroupKey));
        }
        while(itr.hasNext())
        {
            String currentPeptideName = (String) itr.next();
            status = peptideData.get(currentPeptideName).generateRatio();
            peptideData.get(currentPeptideName).generateIntensityAverage(headerName);
//            gro
            if(status == false)
                removingPeptide.add(currentPeptideName);
        }
        itr= removingPeptide.iterator();
        while(itr.hasNext())
        {
            String currentPeptideName = (String) itr.next();
            peptideData.remove(currentPeptideName);
            peptideNameList.remove(currentPeptideName);
        }
    }
    public void printToFile(String inputPath,String directory,String fileName)
    {
         NumberFormat formatter = new DecimalFormat("##.###");
        StringBuffer output = new StringBuffer();
        int compareGroupSize = peptideData.get(peptideNameList.get(0)).getRatioAvg().size();
        output.append("H\t" + new Date() + "\n");
        output.append("H\tTMT Peptide Level\n");
        output.append("H\tPath\t"+inputPath +"\n");

        List groupNameList = new ArrayList(groupKeyMap.keySet());
        for(int i=0;i<groupNameList.size();i++)
        {   
            output.append("H\tSample\t"+ groupNameList.get(i));
            for(int j=0;j<groupKeyMap.get(groupNameList.get(i)).size();j++)
                 output.append("\t" +groupKeyMap.get(groupNameList.get(i)).get(j));
            output.append("\n");
        }
        int counter = 0;
        for(int i=0;i<groupNameList.size();i++)
        {
            for(int j=i+1;j<groupNameList.size();j++)
            {
                counter++;
                output.append("H\tCompare" + counter +"\t" + groupNameList.get(i) + "/" + groupNameList.get(j)+"\n");
            }
        }
        output.append("H\tpep_sequence");
        output.append("\tscan_num");
        
        int totalGroup = ((Peptide)((Collection<Peptide>)peptideData.values()).toArray()[0]).getpValueList().size();
        for(int i =0 ;i<totalGroup;i++)
        {
            output.append("\tpvalue_"+i);
        }
        
        for(int i =0 ;i<totalGroup;i++)
        {
            output.append("\tpvalueRatio_"+i);
        }
        
        for(int i =0 ;i<totalGroup;i++)
        {
            output.append("\tpvalueRatioData_"+i);
        }
        
        for(int i =0 ;i<totalGroup;i++)
        {
            output.append("\tBHCorrection_"+i);
        }
        
        for(int i =0 ;i<totalGroup;i++)
        {
            output.append("\tBHCorrectionRatio_"+i);
        }
        
        for(int i =0;i<groupNameList.size();i++)
        {
            output.append("\t" + groupNameList.get(i) + "_intensity_avg");
       }
        
        for(int i =0;i<compareGroupSize;i++)
        {
            output.append("\tcom_" + (i+1) + "_ratio_avg");
            output.append("\tcom_" + (i+1) + "_ratio_list");            
            output.append("\tcom_" + (i+1) + "_log_ratio_avg");
            output.append("\tcom_" + (i+1) + "_log_ratio_STDEV");
            output.append("\tcom_" + (i+1) + "_log_ratio_RSD");
            output.append("\tcom_" + (i+1) + "_log_ratio_STDEV_paired");
            output.append("\tcom_" + (i+1) + "_log_ratio_RSD_paired");
        }

        for(int i =0 ;i<totalGroup;i++)
        {
            output.append("\tNorm_pvalue_"+i);
        }

        for(int i =0 ;i<totalGroup;i++)
        {
            output.append("\tNorm_pvalueRatio_"+i);
        }

        for(int i =0 ;i<totalGroup;i++)
        {
            output.append("\tNorm_pvalueRatioData_"+i);
        }

        for(int i =0 ;i<totalGroup;i++)
        {
            output.append("\tNorm_BHCorrection_"+i);
        }

        for(int i =0 ;i<totalGroup;i++)
        {
            output.append("\tNorm_BHCorrectionRatio_"+i);
        }

        for(int i =0;i<groupNameList.size();i++)
        {
            output.append("\t" + groupNameList.get(i) + "_norm_intensity_avg");
        }

        for(int i =0;i<compareGroupSize;i++)
        {
            output.append("\tcom_" + (i+1) + "_norm_ratio_avg");
            output.append("\tcom_" + (i+1) + "_norm_ratio_list");
            output.append("\tcom_" + (i+1) + "_norm_log_ratio_avg");
            output.append("\tcom_" + (i+1) + "_norm_log_ratio_STDEV");
            output.append("\tcom_" + (i+1) + "_norm_log_ratio_RSD");
            output.append("\tcom_" + (i+1) + "_norm_log_ratio_STDEV_paired");
            output.append("\tcom_" + (i+1) + "_norm_log_ratio_RSD_paired");
        }

        output.append("\tSpC");
        output.append("\toutlier_count");
        output.append("\tfileName");
        output.append("\taccession");
        output.append("\tGene_Name");
        output.append("\tproteinDescription");
        output.append("\n");
        
        Iterator itr= peptideNameList.iterator();
        while(itr.hasNext())
        {
            Peptide currentPeptide= peptideData.get(itr.next());
            
            output.append("S\t"+ currentPeptide.getSequence());
            output.append("\t" + currentPeptide.getScannum() );
//            output.append("\t" + formatter.format(currentPeptide.getpValue()));
            
            for(double value : currentPeptide.getpValueList())
            {
                output.append("\t" +formatter.format(value));
            }
            
            for(double value : currentPeptide.getpValueRatioList())
            {
                output.append("\t" +formatter.format(value));
            }
            
            
                
            for(List<DescriptiveStatistics> currentGroup : (List<List<DescriptiveStatistics>> )currentPeptide.ratioListForPvalue)
            {
                output.append("\t");
                for(DescriptiveStatistics value : currentGroup)
                 {
                     double v = value.getMean();
                     if(Double.isNaN(v))
                     {
                         output.append("1.0,");
                     }
                     else
                     {
                         output.append(formatter.format(value.getMean())+",");
                     }
                 }
            }
            

            for(double value : currentPeptide.getbHCorrection())
            {
                output.append("\t" +formatter.format(value));
            }
            
            for(double value : currentPeptide.getbHCorrectionRatioList())
            {
                output.append("\t" +formatter.format(value));
            }
            
            
//            for(double value : currentPeptide.getIntensityAverage())
//            {
//                output.append(formatter.format(value)+",");
//            }
            
            
            for(List<Double> currentIntensityGroup : currentPeptide.getIntensityAverage())
            {
                output.append("\t");
                for(Double value : currentIntensityGroup)
                {
                    output.append(formatter.format(value)+",");
                }
            }
            
            
            
            for(int i=0;i<currentPeptide.getRatioAvg().size();i++)
            {
                double value = (double) currentPeptide.getRatioAvg().get(i);
                if(Double.isNaN(value))
                    output.append("\tNA");
                else
                    output.append("\t" + formatter.format(value));
                
                output.append("\t");
                for(int j=0;j<currentPeptide.comparedGroup.size();j++)
                {
                    value = Math.pow(2, currentPeptide.comparedGroup.get(j).get(i).getMean());
//                      value = currentPeptide.comparedGroup.get(j).get(i).getMean();
                    if(Double.isNaN(value))
                        output.append("NA,");
                    else
                        output.append(formatter.format(value)+",");
//                    output.append(currentPeptide.comparedGroup.get(j).get(i).getMean() + ",");
                }
                
                value = (double) currentPeptide.getLogRatioAvg().get(i);
                if(Double.isNaN(value))
                    output.append("\tNA");
                else
                    output.append("\t" + formatter.format(value));
                 
                 value = (double) currentPeptide.getLogRatioStdDev().get(i);
                if(Double.isNaN(value) || Double.isInfinite(value))
                    output.append("\tNA");
                else
                    output.append("\t" + formatter.format(value));
                
                value = (double) currentPeptide.getLogRatioRSD().get(i);
                if(Double.isNaN(value))
                    output.append("\tNA");
                else
                    output.append("\t" + formatter.format(value));
                
                value = (double) currentPeptide.getStdevPaired().get(i);
                if(Double.isNaN(value) || Double.isInfinite(value))
                    output.append("\tNA");
                else
                    output.append("\t" + formatter.format(value));
                
                value = (double) currentPeptide.getRsdPaired().get(i);
                if(Double.isNaN(value) || Double.isInfinite(value))
                    output.append("\tNA");
                else
                    output.append("\t" + formatter.format(value));
                
            }



            for(double value : currentPeptide.getNormPValueList())
            {
                output.append("\t" +formatter.format(value));
            }

            for(double value : currentPeptide.getpValueNormRatioList())
            {
                output.append("\t" +formatter.format(value));
            }



            for(List<DescriptiveStatistics> currentGroup : (List<List<DescriptiveStatistics>> )currentPeptide.normRatioListForPValue)
            {
                output.append("\t");
                for(DescriptiveStatistics value : currentGroup)
                {
                    double v = value.getMean();
                    if(Double.isNaN(v))
                    {
                        output.append("1.0,");
                    }
                    else
                    {
                        output.append(formatter.format(value.getMean())+",");
                    }
                }
            }


            for(double value : currentPeptide.getbHCorrectionNorm())
            {
                output.append("\t" +formatter.format(value));
            }

            for(double value : currentPeptide.getbHCorrectionNormRatioList())
            {
                output.append("\t" +formatter.format(value));
            }


//            for(double value : currentPeptide.getIntensityAverage())
//            {
//                output.append(formatter.format(value)+",");
//            }


            for(List<Double> currentIntensityGroup : currentPeptide.getNormIntensitySum())
            {
                output.append("\t");
                for(Double value : currentIntensityGroup)
                {
                    output.append(formatter.format(value)+",");
                }
            }



            for(int i=0;i<currentPeptide.getRatioNormAvg().size();i++)
            {
                double value = (double) currentPeptide.getRatioNormAvg().get(i);
                if(Double.isNaN(value))
                    output.append("\tNA");
                else
                    output.append("\t" + formatter.format(value));

                output.append("\t");
                for(int j=0;j<currentPeptide.normComparedGroup.size();j++)
                {
                    value = Math.pow(2, currentPeptide.normComparedGroup.get(j).get(i).getMean());
//                      value = currentPeptide.comparedGroup.get(j).get(i).getMean();
                    if(Double.isNaN(value))
                        output.append("1.0,");
                    else
                        output.append(formatter.format(value)+",");
//                    output.append(currentPeptide.comparedGroup.get(j).get(i).getMean() + ",");
                }

                value = (double) currentPeptide.getLogRatioNormAvg().get(i);
                if(Double.isNaN(value))
                    output.append("\tNA");
                else
                    output.append("\t" + formatter.format(value));

                value = (double) currentPeptide.getLogRatioNormStdDev().get(i);
                if(Double.isNaN(value) || Double.isInfinite(value))
                    output.append("\tNA");
                else
                    output.append("\t" + formatter.format(value));

                value = (double) currentPeptide.getLogRatioNormRSD().get(i);
                if(Double.isNaN(value))
                    output.append("\tNA");
                else
                    output.append("\t" + formatter.format(value));

                value = (double) currentPeptide.getNormStdevPaired().get(i);
                if(Double.isNaN(value) || Double.isInfinite(value))
                    output.append("\tNA");
                else
                    output.append("\t" + formatter.format(value));

                value = (double) currentPeptide.getNormRsdPaired().get(i);
                if(Double.isNaN(value) || Double.isInfinite(value))
                    output.append("\tNA");
                else
                    output.append("\t" + formatter.format(value));

            }
            output.append("\t" + currentPeptide.getSpc() );
            output.append("\t" + currentPeptide.getDiscaredCounter() );
            output.append("\t" + currentPeptide.getFilename() );                               
            output.append("\t" + currentPeptide.getProteinLocus().toString() );
            output.append("\t" + currentPeptide.getGeneNameSet().toString() );
            output.append("\t" + currentPeptide.getProteinDescription().toString() );

            
            output.append("\n");
        }
        
        BufferedWriter bufferedWriter = null;
        try{
            if(directory == "")
                bufferedWriter = new BufferedWriter(new FileWriter(fileName));
            else
                bufferedWriter = new BufferedWriter(new FileWriter(directory+ File.separator+fileName));
            bufferedWriter.write(output.toString());
            bufferedWriter.close();
            System.out.println("See the output File at " + directory + File.separator +fileName );
        }
        catch(Exception e)
        {
            Logger.getLogger(TmtfilterPeptide.class.getName()).log(Level.SEVERE, null, e);
        }
    }
    public void runTmtPeptideAnalysis(Map groupMap,String inputPath,String outputPath, String outPutFileName,boolean normalizedValue) throws Exception {
        isNormalized = normalizedValue;
        parseCensusFile(inputPath);
//        generateGroupList();
        groupNameToKeys(groupMap);
        generateGroupInfo();
        filterValues(0.05);
        generateRatio();
//        generatePValues();

        
        printToFile(inputPath,outputPath,outPutFileName);
    }
    public Map getPeptideMap(Map groupMap,String path) {
        parseCensusFile(path);
        groupNameToKeys(groupMap);
        generateGroupInfo();
        filterValues(0.05);
        generateRatio();
        
        return peptideData;
    }
public static void main(String[] args) throws Exception {
/*
      TmtfilterPeptide tmtfilter = new TmtfilterPeptide();
      Map<String, List> groupNameMap = new LinkedHashMap<>();
       List groupList= new ArrayList();
        groupList.add("126.1277");
        groupList.add("127.1248");
       groupNameMap.put("one",groupList);
       
       groupList = new ArrayList();
       groupList.add("128.1344");
       groupList.add("129.1315");
        groupNameMap.put("two", groupList);
       
       groupList = new ArrayList();
       groupList.add("130.1411");
       groupList.add("131.1382");
        groupNameMap.put("three", groupList);
       
      tmtfilter.runTmtPeptideAnalysis(groupNameMap,"C:\\Users\\Harshil\\Documents\\NetBeansProjects\\TmtFilter\\input_data", "census-out-peptideAnalysis.txt");
      tmtfilter.printAvg();//Just a tesing method...........
      //Give the directory of the Census-out.txt file.....
*/
      TmtfilterPeptide tmtfilter = new TmtfilterPeptide();
      Map<String, List> groupNameMap = new LinkedHashMap<>();
       List groupList= new ArrayList();
       groupList.add("126.1277");
        groupNameMap.put("control1", groupList);
        groupList = new ArrayList();
        
        groupList.add("127.1311");
        groupNameMap.put("control2", groupList);
        groupList = new ArrayList();
        
        groupList.add("128.1344");
        groupNameMap.put("control3", groupList);
        groupList = new ArrayList();
        
        groupList.add("129.1378");
        groupNameMap.put("control4", groupList);
        groupList = new ArrayList();
        
        groupList.add("130.1411");
        groupNameMap.put("control5", groupList);
        groupList = new ArrayList();
        
        groupList.add("131.1382");
        groupNameMap.put("control6", groupList);
  /*      
        groupList = new ArrayList();
        groupList.add("128.1344");
        groupList.add("131.1382");
        groupNameMap.put("test2", groupList);
*/
//      tmtfilter.runTmtPeptideAnalysis(groupNameMap,"C:\\Users\\Harshil\\Documents\\NetBeansProjects\\TmtFilter\\input_data","C:\\Users\\Harshil\\Documents\\NetBeansProjects\\TmtFilter\\input_data", "census-out-peptideAnalysis.txt",true);
        tmtfilter.runTmtPeptideAnalysis(groupNameMap,"C:\\Users\\Harshil\\Desktop","C:\\Users\\Harshil\\Desktop", "census-out-peptideAnalysisSimple.txt",true);        
//      tmtfilter.runTmtPeptideAnalysis(groupNameMap,"/home/rpark", "/home/rpark", "census-out-peptideAnalysis.txt",true);
//      tmtfilter.runTmtPeptideAnalysis(groupNameMap,args[0], "TMT_peptide_comare.txt");
//        tmtfilter.runTmtPeptideAnalysis(groupNameMap,"/home/rpark/test_data/meha_tmt", "TMT_peptide_comare.txt");
        
      //Give the directory of the Census-out.txt file.....
//      tmtfilter.printAvg();//Just a tesing method...........

    }
    
    public  void printAvg()//Just a tesing method...........
    {
        System.out.println("Printing the average values for testing .........");
        Iterator itr = peptideData.keySet().iterator();
        while(itr.hasNext())
        {
            Peptide currentPeptide = peptideData.get(itr.next());
            System.out.print(currentPeptide.getSequence());
            System.out.print("\t");
            
            for(int i=0;i<normintensity_name.size();i++)
            {
                DescriptiveStatistics stat = new DescriptiveStatistics();
                int limit =currentPeptide.getallIntensity(normintensity_name.get(0)).size();
                for(int j=0;j<limit;j++)
                {
                
                    stat.addValue( Double.parseDouble(currentPeptide.getallIntensity(normintensity_name.get(i)).get(j).toString()));
                }
                    System.out.print("\t"+ stat.getMean());
            }
        }
    }
}
