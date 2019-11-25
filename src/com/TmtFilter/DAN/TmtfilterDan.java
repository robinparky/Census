 /*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package com.TmtFilter.DAN;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

/**
 *
 * @author Harshil
 * There is a Main HashMap "peptide_data" Contains the Key(Peptide ) to the Value(object of "pep_details" type which contains the data related to the peptide(Ex. the whole row).)
 *  
 * To run this file type 

java peptide_process.pepseq_v3 /data/6/harshil/peptide_process/census-out.txt

In CMD so java classfile pathofthecensus directory......
 * 
 */
public class TmtfilterDan 
{
  //  public static List<Float> sumForAverage = new ArrayList<Float>();
      public  Map allIntensitySum = new HashMap();
      public  ArrayList<String> intensityNames=new ArrayList();

   
      
    public static void main(String[] args) 
    {
        TmtfilterDan dan = new TmtfilterDan();
        dan.parseCensusFile("C:\\Users\\Harshil\\Desktop\\DAn-tmt\\census-out.txt");
        
    }
    
    public void parseCensusFile(String fullPath)
    {
        BufferedReader br = null;
        String protine,peptide,directory="";
        String words[]=new String[10000];
        String words_peptide[]=new String[10000];
        String header_line[]=new String[150];
        multihash hash=new multihash();
        ArrayList<String> ar=new ArrayList();
        ArrayList<String> peptide_list=new ArrayList();
        ArrayList<String> normintensity_name=new ArrayList();
        
        Map<String,pep_details> peptide_data=new HashMap(); // The main map for contains the peptide sequence name and all the datails of each peptide sequence.
        DescriptiveStatistics stats = new DescriptiveStatistics();
        
  //      ratioOnnorm rn = new ratioOnnorm();
        

        int count_intensity = 0;
        boolean checkzero=false;
        int position_spc = 0,position_scn = 0,position_cstate = 0,position_file = 0;
	

	System.out.println("Please wait...");

        try
        {
            //   C:\\testing_data\\test072213_1tot1\\census-out.txt
                    fullPath += File.separator +"census-out.txt";

            File f= new File(fullPath);
          
	     if(fullPath.contains(File.separator))        
                directory=f.getParent();
            else
                directory = "";
	    directory = f.getParent();
             
//          File f = new File("C:\\Users\\Harshil\\Desktop\\DAn-tmt\\census-out.txt");
            br = new BufferedReader(new FileReader(f));
            String line;
            line=br.readLine();
            pep_details pep_details = new pep_details();        
           while(line != null)
            
          // for(int i=0;i<250;i++)
            {
                checkzero=false;
                pep_details=new pep_details();
                stats = new DescriptiveStatistics();
                if(line.startsWith("H") || line.startsWith("h")) 
                {
                   // System.out.println("Its a header line.");
                    header_line= line.split("[\t]");
                        if(header_line[1].equals("SLINE"))
                        {
                           for(int temp=4;temp<header_line.length;temp++)
                           {
                               String splitter[] =new String[10];
                               splitter=header_line[temp].split("[_()]");
                               if(splitter[0].equals("norm") && !splitter[1].equals("ratio"))
                               {
                                   normintensity_name.add(header_line[temp]);
                                  
                               }
                                 if(header_line[temp].contains("_") && !header_line[temp].contains("norm"))
                               {
                                   intensityNames.add(header_line[temp]);
                               }
                               
                               if(header_line[temp].equals("SpC"))
                               {    
                                   count_intensity=temp-4; // count total number of intensity.
                                   position_spc=temp-1;
                                   position_scn=position_spc+1;
                                   position_cstate=position_scn+1;
                                   position_file=position_cstate+1;  
                                   break;
                               }
                               ar.add(header_line[temp]);
                            }
                    // System.out.println("The sline header is " + line);                
                        }
                        
//                for(int j=0;j<count_intensity ;j++) {
//                        sumForAverage.add(j, (float)0);
//                    }
                }
                if(line.startsWith("P") || line.startsWith("p"))
                {
               //    System.out.print("Its a Protine line. \t The name is ");
                    words=line.split("[\t]");
            //        System.out.println(words[1]);
                }
                if(line.startsWith("S") || line.startsWith("s"))
                 {
                     //System.out.print("Its a Peptide line.");
                     words_peptide=line.split("[\t]");
                     for(int temp=0;temp<count_intensity;temp++)
                     {
                          if(Float.parseFloat(words_peptide[temp+4-1]) == 0)
                            {
                                checkzero=true;
                                break;
                            }
                     }
                     if(checkzero)
                         {
                            line=br.readLine();
                            continue;
                         }   
                     if(peptide_data.containsKey(words_peptide[2]))
                     {
                         pep_details=peptide_data.get(words_peptide[2]);
                         int peptide_count=pep_details.getPeptide_counter(words_peptide[2]);
                         pep_details.setPeptide_counter(words_peptide[2]);
                        
                         
                         for(int temp=0;temp<count_intensity;temp++)
                        {
                           
                            float z=pep_details.getIntensity(header_line[temp+4]) ;
//                            if(words_peptide[2].equals("R.ENEMDENLEQVSGIIGNLR.H"))
//                                 System.out.println(z+ "*"+  peptide_count + "+ " + words_peptide[temp+4-1] +"/" + (peptide_count+1) );
                            float new_average=(float) ((z*peptide_count + Float.parseFloat(words_peptide[temp+4-1]))/(peptide_count+1.0));
                           
                            pep_details.setIntensity(header_line[temp+4], new_average );
                            
                            pep_details.setallIntensity(header_line[temp+4], Float.parseFloat(words_peptide[temp+4-1]));
                            
                            // Adding the data to the standard deviation to find it.
                            stats=pep_details.getstat(header_line[temp+4]);
                            stats.addValue(Float.parseFloat(words_peptide[temp+4-1]));
                            pep_details.setstat(header_line[temp+4],stats);
                            
                         }
                         pep_details.setSpc(Integer.parseInt(words_peptide[position_spc]));
                         pep_details.addScannum(Integer.parseInt(words_peptide[position_scn]));
                         pep_details.setCstate(Integer.parseInt(words_peptide[position_cstate]));
                         pep_details.setratio(normintensity_name);
                         pep_details.setFilename(words_peptide[position_file]);
                        
                     }
                     else
                     {
                         for (int temp = 0; temp < count_intensity; temp++) {
                             if (Float.parseFloat(words_peptide[temp + 4 - 1]) == 0) {
                                 checkzero = true;
                                 break;
                             }
                         }
                         if (checkzero)
                         {
                            line=br.readLine();
                            continue;
                         }   
                         if(!peptide_list.contains(words_peptide[2]))
                            peptide_list.add(words_peptide[2]);
                         else
                           System.out.println("Error in addinf peptide at the peptide list");
                         for(int temp=0;temp<count_intensity;temp++)
                        {
                            stats = new DescriptiveStatistics();
                                
                            pep_details.setIntensity(header_line[temp+4],Float.parseFloat(words_peptide[temp+4-1] ));
                           
                            stats.addValue(Float.parseFloat(words_peptide[temp+4-1]));
                            pep_details.setstat(header_line[temp+4],stats);
                            pep_details.setallIntensity(header_line[temp+4], Float.parseFloat(words_peptide[temp+4-1]));

                         //   System.out.println("The Peptide detais of : " + words_peptide[2] + " " + pep_details.getstat(header_line[temp+4]).getStandardDeviation());
                         }
                         
                         pep_details.setPeptide_counter(words_peptide[2]);
                         pep_details.setSpc(Integer.parseInt(words_peptide[position_spc]));
                         pep_details.addScannum(Integer.parseInt(words_peptide[position_scn]));
                         pep_details.setCstate(Integer.parseInt(words_peptide[position_cstate]));
                         pep_details.setratio(normintensity_name);

                         pep_details.setFilename(words_peptide[position_file]);
                         
                     }
//                            pep_details.setratio_avg();
                         
                          peptide_data.put(words_peptide[2],pep_details);
                         hash.add_data(words_peptide[2], words[1]);
                 }
                line=br.readLine();
            }
          // System.out.println(peptide_data.get("R.RADQLADESLESTRR.M").getSpc());
          //  System.out.println(normintensity_name);
            
             findSum(peptide_data);   
            
            setRatio(peptide_data, normintensity_name);
            
               hash.display(peptide_data,ar,directory);
        }
        catch (FileNotFoundException ex)
		 {        
            Logger.getLogger(TmtfilterDan.class.getName()).log(Level.SEVERE, null, ex);
       	 }
		catch (IOException ex)
		 {
            Logger.getLogger(TmtfilterDan.class.getName()).log(Level.SEVERE, null, ex);
      		  }        finally
        {
           
          }
    }
    public  void findSum(Map peptide_data)
    {
        Iterator itr = peptide_data.keySet().iterator();
        while(itr.hasNext())
        {
            String currentPeptideName = (String) itr.next();
            pep_details currnetPeptideData= (pep_details) peptide_data.get(currentPeptideName);
            Map allIntensity = currnetPeptideData.getallIntensity();
            Iterator itr1 = allIntensity.keySet().iterator();
            while(itr1.hasNext())
            {
                String currentIntensityName = (String) itr1.next();
                List currentIntensityData = (List) allIntensity.get(currentIntensityName);
                float total = 0;
                for(Object z : currentIntensityData)
                {
                    total+= Float.parseFloat(z.toString());
                }
                if(allIntensitySum.containsKey(currentIntensityName))
                {
                    float temp = (float) allIntensitySum.get(currentIntensityName);
                    temp += total;
                    allIntensitySum.put(currentIntensityName, temp);
                }
                else
                {
                    allIntensitySum.put(currentIntensityName, total);
                }
            }
        }
        
    }   
    
    
    public static void setRatio(Map peptide_data, List norm_details)
    {
        Iterator keys = peptide_data.keySet().iterator();
        pep_details pep_details= new pep_details();
        while(keys.hasNext())
        {
            String currentPeptide = (String) keys.next();
            pep_details = (pep_details) peptide_data.get(currentPeptide);
            float sum[] = new float[3];
            List n1,n2;
            int i;
            for( i=0;i<pep_details.getallIntensity((String)norm_details.get(0)).size();i++)
            {
                
                sum[0]+= ((float)pep_details.getallIntensity((String)norm_details.get(0)).get(i)) / (float)pep_details.getallIntensity((String)norm_details.get(3)).get(i) * 1.0;
                sum[1]+= ((float)pep_details.getallIntensity((String)norm_details.get(1)).get(i)) / (float)pep_details.getallIntensity((String)norm_details.get(4)).get(i) *1.0;
                sum[2]+= ((float)pep_details.getallIntensity((String)norm_details.get(2)).get(i)) / (float)pep_details.getallIntensity((String)norm_details.get(5)).get(i)*1.0;
//                System.out.println(sum[0]+  " "+ sum[1] + " " + sum[2]); 
            } 
             sum[0] /= i*1.0;sum[1] /=i*1.0;sum[2] =(float) (sum[2]/i*1.0);
//             sum[0] = (sum[0]/((float)allIntensitySum.get( (String)norm_details.get(0)) )) * ((float)allIntensitySum.get( (String)norm_details.get(3)) );
//            sum[1] = (sum[1]/((float)allIntensitySum.get( (String)norm_details.get(1)) )) * ((float)allIntensitySum.get( (String)norm_details.get(4)) );
//            sum[2] =( sum[2]/((float)allIntensitySum.get( (String)norm_details.get(2)) )) * ((float)allIntensitySum.get( (String)norm_details.get(5)) );
//            if(currentPeptide.equals("K.GETPVNSTMSIGQAR.K"))
//                System.out.println(sum[0]+  " "+ sum[1] + " " + sum[2]);
            pep_details.setratio_avg(sum);
            
        }
    
    }}
