/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.scripps.pms.census.io;

import edu.scripps.pms.census.conf.Configuration;
import edu.scripps.pms.census.model.ChroProtein;
import edu.scripps.pms.census.model.ChroPeptide;
import java.util.ArrayList;
import org.jdom.*;
import org.jdom.output.*;
import org.jdom.input.*;

import java.util.*;

/**
 *
 * @author rpark
 */
public class Ms2ReportGenerator {
    
    public static void main(String[] args) throws Exception {
        

        ChroXmlReader cr = new ChroXmlReader(args[0]);
        ArrayList<ChroProtein> proteinList = cr.getProteinList();
        
        double intThreshold = 0;
        if(args[0].length()>1 || args[1].equals("i"))
            intThreshold = Double.parseDouble(args[2]);
            
        StringBuffer sb = new StringBuffer();
        for(Iterator<ChroProtein> itr = proteinList.iterator(); itr.hasNext(); )
        {
            ChroProtein protein = itr.next();

            StringBuffer proSb = new StringBuffer();
            proSb.append(protein.getLocus()).append("\t");
            proSb.append(protein.getSeqCount()).append("\t");
            proSb.append(protein.getSeqCoverage()).append("\t");
            proSb.append(protein.getSpectrumCount()).append("\t");
            proSb.append(protein.getLength()).append("\t");
            proSb.append(protein.getMolWt()).append("\t");
            proSb.append(protein.getPI()).append("\t");
            proSb.append(protein.getValidation()).append("\t");
            proSb.append(protein.getDescription()).append("\t");
            
            List<ChroPeptide> pepList = protein.getPeptideList();
            for(Iterator<ChroPeptide> itr2 = pepList.iterator(); itr2.hasNext(); )
            {
                sb.append(proSb.toString());
                ChroPeptide pep = itr2.next();
                
                sb.append(pep.isUnique()).append("\t");
                sb.append(pep.getFileName()).append("\t");
                sb.append(pep.getScanNum()).append("\t");
                sb.append(pep.getSequence()).append("\t");
                sb.append(pep.getXCorr()).append("\t");
                sb.append(pep.getDeltCN()).append("\t");
                sb.append(pep.getCalcMHplus()).append("\t");
                sb.append(pep.getMhPlus()).append("\t");
                sb.append(pep.getTotalIntensity()).append("\t");
                sb.append(pep.getSpRank()).append("\t");
                sb.append(pep.getSpScore()).append("\t");
                sb.append(pep.getChargeState()).append("\t");
                sb.append(pep.getDeltMass()).append("\t");
                sb.append(pep.getSpecCount()).append("\t");
                
                String bsText = pep.getBsText();
                String brText = pep.getBrText();
                String ysText = pep.getYsText();
                String yrText = pep.getYrText();
                
                String[] bsArr = bsText.split(" ");
                String[] brArr = brText.split(" ");
                String[] ysArr = ysText.split(" ");
                String[] yrArr = yrText.split(" ");
                
                for(int i=0;i<bsArr.length;i+=2) {
		    if("".equals(bsArr[i])) continue;
                    double mass = Double.parseDouble(bsArr[i]);
                    if(mass<=0) continue;
                    
                    try {
                        double bs = Double.parseDouble(bsArr[i+1]);
                        double bf = Double.parseDouble(brArr[i+1]);
                        
                        double sum = bs + bf;
                        //if(sum<intThreshold)
                            sb.append("b").append(mass).append("\t").append(bsArr[i+1]).append("/").append(brArr[i+1]).append("\t");
                        //else {
                        //    double ratio = bs/bf;                        
                        //    sb.append("b").append(mass).append("\t").append(ratio).append("\t");                            
                        //}
                            
                    } catch(Exception e) {
                        sb.append("b").append(mass).append("\t").append(bsArr[i+1]).append("/").append(brArr[i+1]).append("\t");
                    }
                }
                   
                for(int i=0;i<ysArr.length;i+=2) {
		    if("".equals(ysArr[i])) continue;
                    double mass = Double.parseDouble(ysArr[i]);
                    if(mass<=0) continue;
                    
                    try {
                        double ys = Double.parseDouble(ysArr[i+1]);
                        double yf = Double.parseDouble(yrArr[i+1]);
                        
                        
                        double sum = ys + yf;
			sb.append("y").append(mass).append("\t").append(ysArr[i+1]).append("/").append(yrArr[i+1]).append("\t");

/*
                        if(sum<intThreshold)
                            sb.append("y").append(mass).append("\t").append(ysArr[i+1]).append("/").append(yrArr[i+1]).append("\t");
                        else {
                            double ratio = ys/yf;                        
                            sb.append("y").append(mass).append("\t").append(ratio).append("\t");                            
                        }
*/
                    } catch(Exception e) {
                        sb.append("y").append(mass).append("\t").append(ysArr[i+1]).append("/").append(yrArr[i+1]).append("\t");
                    }
                }

                for(Iterator<ChroProtein> redItr = protein.getRedunList().iterator(); redItr.hasNext(); ) {
                    ChroProtein redPro = redItr.next();
                    sb.append(redPro.getLocus()).append("\t").append(redPro.getDescription()).append("\t");
                }
                            
                /*
                System.out.println(bsText);
                System.out.println(brText);
                System.out.println(ysText);
                System.out.println(yrText);
                */
                sb.append("\n");
                
            }
                
                
            
            
            
            
           
        }
        
        System.out.println(sb.toString());
                        
    }
    
}
