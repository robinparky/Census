/*
 * ReportGenerator.java
 *
 * Created on April 25, 2007, 3:02 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package edu.scripps.pms.census.util;

import javax.swing.*;
import java.io.*;
import java.util.*;

import edu.scripps.pms.census.conf.*;
import edu.scripps.pms.census.*;
import edu.scripps.pms.census.model.*;
import edu.scripps.pms.census.io.ChroXmlReader;

/**
 *
 * @author rpark
 */
public class ReportGenerator {
    
    /** Creates a new instance of ReportGenerator */
    public ReportGenerator() {
    }
 
    public static void exportITRAQMultipleReport(String chroFileName, final File file) throws Exception
    {
        exportITRAQMultipleReport(chroFileName, file, null);
    }
    
    public static void exportITRAQMultipleReport(final String chroFileName, final File file, JFrame parentFrame) throws Exception
    {
        try
        {
            final ChroXmlReader cr = new ChroXmlReader(new File(chroFileName)); //should be called before using configuration
            Configuration conf = Configuration.getInstance();
                        
            final List<ReportIon> massMonitorList = conf.getReportIonList();
            
                        
            Thread t = new Thread() {
                    private boolean isSuccessful=true;
                    private String errorMessage="";
                    private PrintStream p = new PrintStream( new BufferedOutputStream(new FileOutputStream(file)));    
                    public void run() 
                    {
                        
                        
			RelExMainFrame.printHeader(p);
                        p.println("H\tCensus msms analysis");                        
                        p.print("H\tcreated date\t"); p.println(new Date());                                             
                        p.print("H\tPLINE\tLOCUS\tSPEC_COUNT\t");

                        StringBuffer result = new StringBuffer();
                        
                        
			for(Iterator<ReportIon> itr=massMonitorList.iterator(); itr.hasNext(); )
			{
                            ReportIon ri = itr.next();
                            double each = ri.getMass();

			    p.print("m/z_");
                            p.print(each);
			    p.print("_total_int\t");                                            
			}                                               

                        p.println("DESCRIPTION");
                        
                        p.print("H\tSLINE\tUNIQUE\tSEQUENCE\t"); //RATIO\tREGRESSION_FACTOR\tDETERMINANT_FACTOR\tSAM_INT\tREF_INT\tSIGNAL_TO_NOISE_RATIO\tFILE_NAME");
                        
                        for(Iterator<ReportIon> itr=massMonitorList.iterator(); itr.hasNext(); )
			{
                            ReportIon ri = itr.next();
                            double each = ri.getMass();

			    p.print("area_m/z_");
                            p.print(each);
                            p.print("");
			    p.print("\t");
                        }

                        String massStr = massMonitorList.get(0).toString();
                        for(int i=1;i<massMonitorList.size();i++)
                        {
                            p.print("ratio(" + massStr + "/" + massMonitorList.get(i) + ")");                            
                        }
                        
                        p.println();


                        int totalCount=0;
                        int quantifiedCount=0;
                        int tempaaa=0;

                        try {
                            
                            ArrayList<ChroProtein> list = cr.getProteinList();
                                                        
                            for(Iterator<ChroProtein> proItr = list.iterator(); proItr.hasNext(); )
                            {                    
                                ChroProtein protein = proItr.next();            
                            
                                //**  for normalization **
                                //int totoalSpecAndProteinLength = 0; 
                                int proLength = Integer.parseInt(protein.getLength());
                                String[] specArr = protein.getSpectrumCount().split(",");

                                StringBuffer proteinSb = new StringBuffer();
                                
                                StringBuffer peptideLines = new StringBuffer();
                                List<ChroPeptide> peptideList = protein.getPeptideList();
                                //List<ChroPeptide> tempPepList = new Vector<ChroPeptide>();

                                for(Iterator<ChroPeptide> pepItr = peptideList.iterator(); pepItr.hasNext(); )
                                {   
                                    ChroPeptide peptide = pepItr.next();                                
                                    
                                    totalCount++;

                                    
                                    List dataList = peptide.getDataList();                                                                               

                                    long[] intensitySumArr = new long[massMonitorList.size()];

                                    int startRange = Integer.parseInt(peptide.getStartRange());
                                    int endRange = Integer.parseInt(peptide.getEndRange());

                                    //LinearRegression reg = null;
                                    for(Iterator<ChroiTRAQLabelData> dataItr=dataList.iterator(); dataItr.hasNext(); )
                                    {
                                        ChroiTRAQLabelData eachData = dataItr.next();

                                        int scanNum = eachData.getScanNum();
                                        long[] intenArr = eachData.getIntensityArr();

                                        for(int i=0;i<intenArr.length;i++)
                                        {                                    
                                            if(scanNum>=startRange && scanNum<=endRange)
                                                intensitySumArr[i] += intenArr[i];
                                        }

                                    }
                                    
                                    peptide.setTotalIntArr(intensitySumArr);
                                    
                                    
                                    
                                    
                                    //tempPepList.add(peptide);
                                    
                                }
               
                                proteinSb.append("P\t"); 
                                proteinSb.append(protein.getLocus()).append("\t");
                                proteinSb.append(protein.getSpectrumCount()).append("\t");
                 
                                double devSum=0;

                                StringBuffer pepSb = new StringBuffer();
 
				int tmpCount=0;

				for(Iterator<ChroPeptide> tempItr=peptideList.iterator(); tempItr.hasNext(); )
				{
				    ChroPeptide each = tempItr.next();

				    pepSb.append("S\t");
				    pepSb.append(each.isUnique()?"U":"");
				    pepSb.append("\t");
				    pepSb.append(each.getSequence());
				    pepSb.append("\t");

                                    long[] larr = each.getTotalIntArr();
                                    
                                    for(long l : larr)
                                        pepSb.append(l).append("\t");
                                    
                                    
                                    for(int i=1;i<larr.length;i++)
                                    {                                        
                                        pepSb.append( CensusHelper.format.format( (double)larr[0]/larr[i]) ).append("\t");
                                    }
                                    
    
				    pepSb.append(each.getFileName()).append("\n");

				    tmpCount++;
				}

                                
                                proteinSb.append(protein.getDescription());
                                proteinSb.append(peptideLines.toString());                
                                proteinSb.append("\n");

                                if( pepSb.length()<=0 )
                                    continue;

                                result.append(proteinSb.toString());
                                result.append(pepSb.toString());		
                                
                            }
                        }
                        catch(Exception e)  
                        {
                            isSuccessful=false;
                            errorMessage = e.getMessage();
                            e.printStackTrace();
                        }

                        SwingUtilities.invokeLater(new Runnable() {
                                public void run() {
                                    System.out.println("census-out.txt was generated successfully.");
                            }

                        }                    
                        );

                        
                        p.print(result.toString());

                        if(null != p)
                            p.close();  


                    }

                };
                
        
            t.start();            
            
        } catch(IOException e)
        {
            System.out.println("Failed to write file" + e);
            e.printStackTrace();
            
        }

    }
    
}
