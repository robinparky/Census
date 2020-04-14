package edu.scripps.pms.census.labelFree.util;

import java.sql.Connection;

public class LabelfreeTempDB {

    public final String path;
    private Connection conn = null;



    public LabelfreeTempDB(String path) {
        this.path = path;
    }


    public static class Peptide
    {
        String file;
        int scan;
        String seq;
        double xcorr;
        double calcMHPlus;
        double mhPlus;
        double totalIntensity;
        int spRange;
        double redundancy;
        double spScore;
        double deltaCN;
        double deltaMass;
        int charge;
        double spC;
        double rt;
        double lit;
        double lightStartMass;
        double lightAvgMass;
        double peakSigma;
        double peakX;
        double peakY;
        double startScan;
        double endScan;
        double startRt;
        double endRT;
        double peakArea;
        double peakSigmaIonInjectionCorrection;
        double peakXIonInjectionCorrection;
        double peakYIonInjectionCorrection;
        double peakAreaIonInjectionCorrection;
        int peakDetected;
        String chro;
        String Peaks;
        String expPath;
    }

    public static class Protein
    {
        String locus;
        int seqCount;
        int specCount;
        double seqCoverage;
        int length;
        double molWt;
        double pi;
        double val;
        String desc;
        String expPath;
    }



}
