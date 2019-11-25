/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package edu.scripps.pms.census.labelFree;
/**
 * @author  Robin Park
 * @version $Id: Replicates.java,v 1.3 2013/05/17 05:02:53 cvsuser Exp $
 */
import java.util.*;
import java.io.Serializable;

public class Replicates implements Serializable 
{

    private String name;
//    private ArrayList<String> pathList;
//    private ArrayList<Integer> searchIdList;
    
    private ArrayList<Replicate> replicateList;

    private int[] columnIndexArr;
    private int[] regIndexArr;

    private String tandemMass;
    
    public Replicates() {
//	pathList = new ArrayList<String>();
//	searchIdList = new ArrayList<Integer>();
        replicateList = new ArrayList<Replicate>();
    }

    public void addReplicate(Replicate replicate) {
            getReplicateList().add(replicate);
    }

    
    /*
    public void addPath(String path) {
	this.pathList.add(path);
    }

    public void addSearchId(Integer searchId) {
	this.searchIdList.add(searchId);
    }
    */

    /**
     * @return the name
     */
    public String getName() {
	return name;
    }

    /**
     * @param name the name to set
     */
    public void setName(String name) {
	this.name = name;
    }

    public ArrayList<Replicate> getReplicateList() {
        return replicateList;
    }

    public void setReplicateList(ArrayList<Replicate> replicateList) {
        this.replicateList = replicateList;
    }


    public int[] getColumnIndexArr() {
        return columnIndexArr;
    }

    public void setColumnIndexArr(int[] columnIndexArr) {
        this.columnIndexArr = columnIndexArr;
    }

    public String getTandemMass() {
	return tandemMass;
    }

    public void setTandemMass(String tandemMass) {
	this.tandemMass = tandemMass;
    }

    public int[] getRegIndexArr() {
        return regIndexArr;
    }

    public void setRegIndexArr(int[] regIndexArr) {
        this.regIndexArr = regIndexArr;
    }

    
    public void sort() {
        Collections.sort(this.replicateList);        
    }
}
