/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.scripps.pms.census.labelFree;

import java.util.AbstractList;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author Rohan
 */
public class LabelfreePeptide {

    private String sequence;
    private int chargestate;
    private ArrayList intensity1 = new ArrayList<>();
    private ArrayList intensity2 = new ArrayList<>();
    private List<String> filename1 = new ArrayList<>();
    private List<String> filename2 = new ArrayList<>();
    private ArrayList xCorr1 = new ArrayList<>();
    private ArrayList xCorr2 = new ArrayList<>();
    private ArrayList dcn1 = new ArrayList<>();
    private ArrayList dcn2 = new ArrayList<>();
    private ArrayList dmass1 = new ArrayList<>();
    private ArrayList dmass2 = new ArrayList<>();
    private ArrayList sprank1 = new ArrayList<>();
    private ArrayList sprank2 = new ArrayList<>();
    private ArrayList spscore1 = new ArrayList<>();
    private ArrayList spscore2 = new ArrayList<>();
    private ArrayList redundancy1 = new ArrayList<>();
    private ArrayList redundancy2 = new ArrayList<>();
    private ArrayList startrange1 = new ArrayList<>();
    private ArrayList startrange2 = new ArrayList<>();
    private ArrayList endrange1 = new ArrayList<>();
    private ArrayList endrange2 = new ArrayList<>();
    private ArrayList retentiontime1 = new ArrayList<>();
    private ArrayList retentiontime2 = new ArrayList<>();
    private ArrayList ioninjection1 = new ArrayList<>();
    private ArrayList ioninjection2 = new ArrayList<>();
    private String protein;
    private String description;

    public String getSequence() {
        return sequence;
    }

    public void setSequence(String sequence) {
        this.sequence = sequence;
    }

    public int getChargestate() {
        return chargestate;
    }

    public void setChargestate(int chargestate) {
        this.chargestate = chargestate;
    }

    public ArrayList getIntensity1() {
        return intensity1;
    }

    public void setIntensity1(ArrayList intensity1) {
        this.intensity1 = intensity1;
    }

    public ArrayList getIntensity2() {
        return intensity2;
    }

    public void setIntensity2(ArrayList intensity2) {
        this.intensity2 = intensity2;
    }

    public List<String> getFilename1() {
        return filename1;
    }

    public void setFilename1(List<String> filename1) {
        this.filename1 = filename1;
    }

    public List<String> getFilename2() {
        return filename2;
    }

    public void setFilename2(List<String> filename2) {
        this.filename2 = filename2;
    }

    public ArrayList getxCorr1() {
        return xCorr1;
    }

    public ArrayList getxCorr2() {
        return xCorr2;
    }

    public ArrayList getDcn1() {
        return dcn1;
    }

    public void setDcn1(ArrayList dcn1) {
        this.dcn1 = dcn1;
    }

    public ArrayList getDcn2() {
        return dcn2;
    }

    public void setDcn2(ArrayList dcn2) {
        this.dcn2 = dcn2;
    }

    public ArrayList getDmass1() {
        return dmass1;
    }

    public void setDmass1(ArrayList dmass1) {
        this.dmass1 = dmass1;
    }

    public ArrayList getDmass2() {
        return dmass2;
    }

    public void setDmass2(ArrayList dmass2) {
        this.dmass2 = dmass2;
    }

    public List<Integer> getSprank1() {
        return sprank1;
    }



    public List<Integer> getSprank2() {
        return sprank2;
    }

    public ArrayList getSpscore1() {
        return spscore1;
    }

    public void setSpscore1(ArrayList spscore1) {
        this.spscore1 = spscore1;
    }

    public ArrayList getSpscore2() {
        return spscore2;
    }

    public void setSpscore2(ArrayList spscore2) {
        this.spscore2 = spscore2;
    }

    public List<Integer> getRedundancy1() {
        return redundancy1;
    }


    public List<Integer> getRedundancy2() {
        return redundancy2;
    }

    public ArrayList getStartrange1() {
        return startrange1;
    }

    public void setStartrange1(ArrayList startrange1) {
        this.startrange1 = startrange1;
    }

    public ArrayList getStartrange2() {
        return startrange2;
    }

    public void setStartrange2(ArrayList startrange2) {
        this.startrange2 = startrange2;
    }

    public ArrayList getEndrange1() {
        return endrange1;
    }

    public void setEndrange1(ArrayList endrange1) {
        this.endrange1 = endrange1;
    }

    public ArrayList getEndrange2() {
        return endrange2;
    }

    public void setEndrange2(ArrayList endrange2) {
        this.endrange2 = endrange2;
    }

    public ArrayList getRetentiontime1() {
        return retentiontime1;
    }

    public void setRetentiontime1(ArrayList retentiontime1) {
        this.retentiontime1 = retentiontime1;
    }

    public ArrayList getRetentiontime2() {
        return retentiontime2;
    }

    public void setRetentiontime2(ArrayList retentiontime2) {
        this.retentiontime2 = retentiontime2;
    }

    public ArrayList getIoninjection1() {
        return ioninjection1;
    }

    public void setIoninjection1(ArrayList ioninjection1) {
        this.ioninjection1 = ioninjection1;
    }

    public ArrayList getIoninjection2() {
        return ioninjection2;
    }

    public void setIoninjection2(ArrayList ioninjection2) {
        this.ioninjection2 = ioninjection2;
    }

    public String getProtein() {
        return protein;
    }

    public void setProtein(String protein) {
        this.protein = protein;
    }

    public String getDescription() {
        return description;
    }

    public void setDescription(String description) {
        this.description = description;
    }
    public void addIntensity1(Object intensity){
        intensity1.add(intensity);
    }
    public void addIntensity2(Object intensity){
        intensity2.add(intensity);
    }
    public void addFilename1(String filename){
        filename1.add(filename);
    }
    public void addFilename2(String filename){
        filename2.add(filename);
    }
    public void addXCorr1(Object xcorr){
        xCorr1.add(xcorr);
    }
    public void addXCorr2(Object xcorr){
        xCorr2.add(xcorr);
    }
    public void addrt1(Object ret){
        retentiontime1.add(ret);
    }
    public void addrt2(Object ret){
        retentiontime2.add(ret);
    }
    public void addion1(Object ion){
        ioninjection1.add(ion);
    }
    public void addion2(Object ion){
        ioninjection2.add(ion);
    }
    public void adddcn1(Object dcn){
        dcn1.add(dcn);
    }
    public void adddcn2(Object dcn){
        dcn2.add(dcn);
    }

    public void adddmass1(Object dmass) {
        dmass1.add(dmass);
    }

    public void adddmass2(Object dmass) {
        dmass2.add(dmass);
    }
    public void addsprank1(Object sprank) {
        sprank1.add(sprank);
    }

    public void addsprank2(Object sprank) {
        sprank2.add(sprank);
    }
    public void addstart1(Object start){
        startrange1.add(start);
    }
    public void addstart2(Object start){
        startrange2.add(start);
    }
    public void addend1(Object end){
        endrange1.add(end);
    }
    public void addend2(Object end){
        endrange2.add(end);
    }
    public void addred1(Object red){
        redundancy1.add(red);
    }
    public void addred2(Object red){
        redundancy2 .add(red);
    }
    public void addsp1(Object sp){
        spscore1.add(sp);
    }
    public void addsp2(Object sp){
        spscore2 .add(sp);
    }
    
}
