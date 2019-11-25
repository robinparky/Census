
/*
* Copyright (c) 2008 The Scripps Research Institute, Yates Lab.  All rights reserved.  
*/

package edu.scripps.pms.census.model;

import java.util.*;
import edu.scripps.pms.census.util.*;

/**
 *
 * @author Sung Kyu, Robin, Park
 * @email rpark@scripps.edu
 * Created on Jun 10, 2008 
 * $Revision: 1.5 $
 * $Date: 2013/07/30 20:21:46 $
 */
public class MultipleChroProtein extends ChroProtein {

    private Hashtable<String, ChroPeptide[]> peptideHt = null;

    private ChroProtein[] proteinArr = null;
    private int expSize=0;
   
    public ChroProtein[] getProteinArr() {
	return proteinArr;
    }

    public MultipleChroProtein(int expSize) {
	this.expSize = expSize;
	proteinArr = new ChroProtein[expSize];
	peptideHt = new Hashtable<String, ChroPeptide[]>();
    }
   
    public void addChroProtein(ChroProtein cPro, int index) {
	proteinArr[index] = cPro;
    }
    
    public void putPeptideHt(ChroPeptide pep, int index) {
        String key = pep.getSequence() + pep.getChargeState();
        ChroPeptide[] arr = getPeptideHt().get(key);

	//if(pep.getAnCompositeScore()<0.5)
	//    return;
/*
	int peakStart = Integer.parseInt( pep.getStartRange() );
	int peakEnd = Integer.parseInt( pep.getEndRange() );
	List dataList = pep.getDataList();
	AllNoneUtil.getANScore(pep, dataList, peakStart, peakEnd);
	System.out.println("====" + pep.getSequence() + " " + pep.getAnCompositeScore());
*/


	//pep.addMultiIntensity(pep.getSamIntensity());

        if(null == arr) {
            arr = new ChroPeptide[expSize];            

            arr[index] = pep;            
            peptideHt.put(key, arr);


//if(pep.getSequence().equals("K.ILKEDLQPSPVCR.N"))
//		    System.out.println("1***===" + arr + " " + index + " " + key + " " + pep.getSamIntensity() + " " + pep.getMultiIntensity() + " " + pep.getMultiAveIntensity());
        } else {
	    ChroPeptide tmpPep = arr[index];

	    if(null != tmpPep) {
                tmpPep.addMultiIntensity(pep.getGaussianPeakArea());
		//tmpPep.addMultiIntensity(pep.getSamIntensity());

//if(pep.getSequence().equals("K.ILKEDLQPSPVCR.N"))
//		    System.out.println("===" + index + " " + key + " " + tmpPep.getMultiIntensity().size() + " " + tmpPep.getMultiIntensity());
		
		//arr[index] = tmpPep;
	    } else {
		arr[index] = pep;
	    }

//if(pep.getSequence().equals("K.ILKEDLQPSPVCR.N"))
	    //System.out.println("==" + aa.getSequence());
        }
/*       
	if(key.equals("K.AQYEEIAQR.S2") && pep.getFileName().contains("3M")) {
		ChroPeptide tmpPep1 = arr[index];
		
		    System.out.println("***");
		ArrayList<Double> mi = tmpPep1.getMultiIntensity();

		for (Double d:mi) {
		    System.out.println("***====" + d + "\t" + index);
		}

		System.out.println("22====" + index + "\t" + key + "\t" +  arr[index].getSamIntensity() + "\t" + arr[index].getMultiIntensitySum() + "\t");

	    }
*/
        getPeptideHt().put(pep.getSequence() + pep.getChargeState(), arr);
	/*
	if(key.equals("K.AQYEEIAQR.S2") && pep.getFileName().contains("3M")) {
	    System.out.println("22====" + index + "\t" + key + "\t" +  arr[index].getSamIntensity() + "\t" + arr[index].getMultiAveIntensity() + "\t");
	}*/
    }

    public String toString() {
	StringBuffer sb = new StringBuffer();

	for(ChroProtein cPro : proteinArr)
	{
	    if(null == cPro)
		sb.append("null");
	    else
		sb.append(cPro.getLocus());
	
	    sb.append("\t");
	}

	return sb.toString();
    }

    public Hashtable<String, ChroPeptide[]> getPeptideHt() {
        return peptideHt;
    }

    public void setPeptideHt(Hashtable<String, ChroPeptide[]> peptideHt) {
        this.peptideHt = peptideHt;
    }
}
