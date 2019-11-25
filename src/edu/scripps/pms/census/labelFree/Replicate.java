/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package edu.scripps.pms.census.labelFree;

/*
* Copyright (c) 2008 Integrated Proteomics Applications.  All rights reserved.  
*/
/**
 *
 * @author Sung Kyu, Robin, Park
 * @email robinparky@yahoo.com
 * Created on Mar 30, 2009 
 * $Revision: 1.5 $
 * $Date: 2013/05/17 05:02:53 $
 */
import java.util.Date;

public class Replicate implements Comparable {


    private Date expDate;
    private Date searchDate;
    private String sampleName;
    private String sampleDescription;
    private String searchName;

    private int searchId;
    private String path;
    private String labelType;
    private String labelFolder;
    private boolean labelReverse;
    
    public Replicate() {

    }
    
    public Replicate(String sampleName) {
	this.sampleName = sampleName;

    }

    public Date getExpDate() {
        return expDate;
    }

    public void setExpDate(Date expDate) {
        this.expDate = expDate;
    }

    public Date getSearchDate() {
        return searchDate;
    }

    public void setSearchDate(Date searchDate) {
        this.searchDate = searchDate;
    }

    public String getSampleName() {
        return sampleName;
    }

    public void setSampleName(String sampleName) {
        this.sampleName = sampleName;
    }

    public String getSampleDescription() {
        return sampleDescription;
    }

    public void setSampleDescription(String sampleDescription) {
        this.sampleDescription = sampleDescription;
    }

    public int getSearchId() {
        return searchId;
    }

    public void setSearchId(int searchId) {
        this.searchId = searchId;
    }

    public String getPath() {
        return path;
    }

    public void setPath(String path) {
        this.path = path;
    }

	public void setLabelType(String lt){
		labelType = lt;
	}
	public String getLabelType(){
		return labelType;
	}
	public void setLabelReverse(boolean lr){
		labelReverse = lr;
	}
	public boolean getLabelReverse(){
		if(this.labelType==null){
			return false;
		}
		else{
			if(labelType.contains("_HL") || labelType.contains("_ML") || labelType.contains("_HM") || labelType.contains("_REV")){
				return true;
			}
			else{
				return false;
			}
		}
	}
    public void setLabelFolder(String lf){
            this.labelFolder = lf;
    }
    public String getLabelFolder(){
            if(labelType==null || labelType.contains("FOR") || labelType.contains("REV")){
                    return "";
            }
            else{
                    String lt = labelType.substring(labelType.lastIndexOf("_")+1);
                    if(getLabelReverse()){
                            StringBuffer sbf = new StringBuffer();
                            sbf.append(lt.charAt(1));
                            sbf.append(lt.charAt(0));
                            return sbf.toString();
                    }
                    else{
                            return lt;
                    }
            }
    }

        
    public int compareTo(Object o) throws ClassCastException {

        Replicate r = (Replicate)o;
        
        //System.out.println("===" + sampleName + " " + this.getLabelReverse() + " " + this.labelType + " " +  r.getSampleName() + " " + r.getLabelReverse() + " " + r.getLabelType());
        
        
        //if(this.getLabelReverse() == false && r.getLabelReverse() == true)
        if(r.getLabelReverse() == true) {
            if(this.getLabelReverse()==true)
                return this.getSampleName().compareTo(r.getSampleName());
            else
                return -1;                    
        }
        else 
            return 1;
        //else            
          //  return 

    }

    public String getSearchName() {
        return searchName;
    }

    public void setSearchName(String searchName) {
        this.searchName = searchName;
    }


    
}
