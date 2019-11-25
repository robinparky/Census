package edu.scripps.pms.census.model;

import edu.scripps.pms.census.util.XYPoint;
import edu.scripps.pms.census.util.dtaselect.Peptide;

import java.util.*;

/**
 * Created by Titus Jung titusj@scripps.edu on 1/29/18.
 */
public class DisplayData {


    public static class DisplayProtein {

        private String accession;
        private double averageRatio;
        private double reverseAverageRatio;
        private double compositeRatio;
        private double stDev;
        private double reverseStDev;
        private double weightedAverageRatio;
        private int pepNum;
        private int spectrumCount;
        private String description;

        private List<DisplayPeptide> displayPeptideList = new ArrayList<DisplayPeptide>();

        public DisplayProtein(String accession, double averageRatio, double reverseAverageRatio, double compositeRatio,
                              double stDev, double reverseStDev, double weightedAverageRatio, int pepNum, int spectrumCount,
                              String description, List<DisplayPeptide> displayPeptideList) {
            this.accession = accession;
            this.averageRatio = averageRatio;
            this.reverseAverageRatio = reverseAverageRatio;
            this.compositeRatio = compositeRatio;
            this.stDev = stDev;
            this.reverseStDev = reverseStDev;
            this.weightedAverageRatio = weightedAverageRatio;
            this.pepNum = pepNum;
            this.spectrumCount = spectrumCount;
            this.description = description;
            this.displayPeptideList = displayPeptideList;
        }

        public String getAccession() {
            return accession;
        }

        public void setAccession(String accession) {
            this.accession = accession;
        }

        public double getAverageRatio() {
            return averageRatio;
        }

        public void setAverageRatio(double averageRatio) {
            this.averageRatio = averageRatio;
        }

        public double getReverseAverageRatio() {
            return reverseAverageRatio;
        }

        public void setReverseAverageRatio(double reverseAverageRatio) {
            this.reverseAverageRatio = reverseAverageRatio;
        }

        public double getCompositeRatio() {
            return compositeRatio;
        }

        public void setCompositeRatio(double compositeRatio) {
            this.compositeRatio = compositeRatio;
        }

        public double getStDev() {
            return stDev;
        }

        public void setStDev(double stDev) {
            this.stDev = stDev;
        }

        public double getReverseStDev() {
            return reverseStDev;
        }

        public void setReverseStDev(double reverseStDev) {
            this.reverseStDev = reverseStDev;
        }

        public double getWeightedAverageRatio() {
            return weightedAverageRatio;
        }

        public void setWeightedAverageRatio(double weightedAverageRatio) {
            this.weightedAverageRatio = weightedAverageRatio;
        }

        public int getPepNum() {
            return pepNum;
        }

        public void setPepNum(int pepNum) {
            this.pepNum = pepNum;
        }

        public int getSpectrumCount() {
            return spectrumCount;
        }

        public void setSpectrumCount(int spectrumCount) {
            this.spectrumCount = spectrumCount;
        }

        public String getDescription() {
            return description;
        }

        public void setDescription(String description) {
            this.description = description;
        }

        public List<DisplayPeptide> getDisplayPeptideList() {
            return displayPeptideList;
        }

        public void setDisplayPeptideList(List<DisplayPeptide> displayPeptideList) {
            this.displayPeptideList = displayPeptideList;
        }
    }

    public static class DisplayChroDataXY {

        private int x;
        private long y;

        public DisplayChroDataXY(int x, long y) {
            this.x = x;
            this.y = y;
        }

        public int getX() {
            return x;
        }

        public void setX(int x) {
            this.x = x;
        }

        public long getY() {
            return y;
        }

        public void setY(long y) {
            this.y = y;
        }
    }

    public static class DisplayChroData {

        private List<DisplayChroDataXY> data1 = new ArrayList<DisplayChroDataXY>();
        private List<DisplayChroDataXY> data2 = new ArrayList<DisplayChroDataXY>();
        private List<edu.scripps.pms.census.model.XYPoint> xyPoints = new ArrayList<>();
        private String thoMass;

        private int winSize = 5;

        private long maxIntensity;
        private int startRange;
        private int endRange;
        private double massAccuracy;

        public long getMaxIntensity() {
            return maxIntensity;
        }

        public void setMaxIntensity(long maxIntensity) {
            this.maxIntensity = maxIntensity;
        }

        public void smooth() {

            MovingAverage ma = new MovingAverage(winSize);
            for (Iterator<DisplayChroDataXY> itr = data1.iterator(); itr.hasNext(); ) {
                DisplayChroDataXY each = itr.next();
                ma.newNum(each.getY());
                each.setY((long) ma.getAvg());
            }

            ma = new MovingAverage(winSize);
            for (Iterator<DisplayChroDataXY> itr = data2.iterator(); itr.hasNext(); ) {
                DisplayChroDataXY each = itr.next();
                ma.newNum(each.getY());
                each.setY((long) ma.getAvg());
            }
            /*
			 * for(int i=0;i<data1.size()-5;i++) { DisplayChroDataXY dataA =
			 * data1.get(i); DisplayChroDataXY dataB = data1.get(i+5);
			 * dataA.setX(dataB.getX()); dataA.setY(dataB.getY()); }
			 *
			 * for(int i=0;i<data2.size()-5;i++) { DisplayChroDataXY dataA =
			 * data2.get(i); DisplayChroDataXY dataB = data2.get(i+5);
			 * dataA.setX(dataB.getX()); dataA.setY(dataB.getY()); }
			 */

        }

        public List<DisplayChroDataXY> getData1() {
            return data1;
        }

        public void setData1(List<DisplayChroDataXY> data1) {
            this.data1 = data1;
        }

        public List<DisplayChroDataXY> getData2() {
            return data2;
        }

        public void setData2(List<DisplayChroDataXY> data2) {
            this.data2 = data2;
        }

        public int getStartRange() {
            return startRange;
        }

        public void setStartRange(int startRange) {
            this.startRange = startRange;
        }

        public int getEndRange() {
            return endRange;
        }

        public void setEndRange(int endRange) {
            this.endRange = endRange;
        }

        public List<edu.scripps.pms.census.model.XYPoint> getXyPoints() {
            return xyPoints;
        }

        public void setXyPoints(List<edu.scripps.pms.census.model.XYPoint> xyPoints) {
            this.xyPoints = xyPoints;
        }

        public String getThoMass() {
            return thoMass;
        }

        public void setThoMass(String thoMass) {
            this.thoMass = thoMass;
        }

        public double getMassAccuracy() {
            return massAccuracy;
        }

        public void setMassAccuracy(double massAccuracy) {
            this.massAccuracy = massAccuracy;
        }
    }

    public static class DisplayPeptide {

        private boolean unique;
        private String sequence;
        //private String proteinSequence;
        private String ptmIndexInProtein;
        private String fileName;
        private String xCorr;
        private String deltaCN;
        private String ratio;
        private String ratioRev;
        private String regressionFactor;
        private String determinantFactor;
        private String pValue;
        private String samInt;
        private String refInt;
        private String areaRatio;
        private String profScore;
        private String scan;
        private String cs;
        private String enrichment;
        private String proteins;
        private String proteinDescription;
        private DisplayChroData displayChroData = new DisplayChroData();

        private String singleton;

        public DisplayPeptide(boolean unique, String sequence, String fileName, String xCorr, String deltaCN, String ratio,
                       String ratioRev, String regressionFactor, String determinantFactor, String pValue, String samInt,
                       String refInt, String areaRatio, String profScore, String scan, String cs, String enrichment,
                       String proteins, String proteinDescription, DisplayChroData displayChroData, String singleton) {
            this.unique = unique;
            this.sequence = sequence;
            this.fileName = fileName;
            this.xCorr = xCorr;
            this.deltaCN = deltaCN;
            this.ratio = ratio;
            this.ratioRev = ratioRev;
            this.regressionFactor = regressionFactor;
            this.determinantFactor = determinantFactor;
            this.pValue = pValue;
            this.samInt = samInt;
            this.refInt = refInt;
            this.areaRatio = areaRatio;
            this.profScore = profScore;
            this.scan = scan;
            this.cs = cs;
            this.enrichment = enrichment;
            this.proteins = proteins;
            this.proteinDescription = proteinDescription;
            this.displayChroData = displayChroData;
            this.singleton = singleton;
        }

        public DisplayPeptide(Peptide peptide) {
            this.unique = peptide.isUnique();
            this.sequence = peptide.getSequence();
            this.fileName = peptide.getFileName();
            this.xCorr = peptide.getXCorr();
            this.deltaCN = peptide.getDeltCN();
            this.ratio = "";
            this.ratioRev = "";
            this.regressionFactor = "";
            this.determinantFactor = "";
            this.pValue = "";
            this.samInt = samInt;
            this.refInt = refInt;
            this.areaRatio = areaRatio;
            this.profScore = profScore;
            this.scan = scan;
            this.cs = cs;
            this.enrichment = enrichment;
            this.proteins = proteins;
            this.proteinDescription = proteinDescription;
            this.displayChroData = displayChroData;
            this.singleton = singleton;

        }

        public boolean isUnique() {
            return unique;
        }

        public void setUnique(boolean unique) {
            this.unique = unique;
        }

        public String getSequence() {
            return sequence;
        }

        public void setSequence(String sequence) {
            this.sequence = sequence;
        }

        public String getFileName() {
            return fileName;
        }

        public void setFileName(String fileName) {
            this.fileName = fileName;
        }

        public String getxCorr() {
            return xCorr;
        }

        public void setxCorr(String xCorr) {
            this.xCorr = xCorr;
        }

        public String getDeltaCN() {
            return deltaCN;
        }

        public void setDeltaCN(String deltaCN) {
            this.deltaCN = deltaCN;
        }

        public String getRatio() {
            return ratio;
        }

        public void setRatio(String ratio) {
            this.ratio = ratio;
        }

        public String getRatioRev() {
            return ratioRev;
        }

        public void setRatioRev(String ratioRev) {
            this.ratioRev = ratioRev;
        }

        public String getRegressionFactor() {
            return regressionFactor;
        }

        public void setRegressionFactor(String regressionFactor) {
            this.regressionFactor = regressionFactor;
        }

        public String getDeterminantFactor() {
            return determinantFactor;
        }

        public void setDeterminantFactor(String determinantFactor) {
            this.determinantFactor = determinantFactor;
        }

        public String getpValue() {
            return pValue;
        }

        public void setpValue(String pValue) {
            this.pValue = pValue;
        }

        public String getSamInt() {
            return samInt;
        }

        public void setSamInt(String samInt) {
            this.samInt = samInt;
        }

        public String getRefInt() {
            return refInt;
        }

        public void setRefInt(String refInt) {
            this.refInt = refInt;
        }

        public String getAreaRatio() {
            return areaRatio;
        }

        public void setAreaRatio(String areaRatio) {
            this.areaRatio = areaRatio;
        }

        public String getProfScore() {
            return profScore;
        }

        public void setProfScore(String profScore) {
            this.profScore = profScore;
        }

        public String getScan() {
            return scan;
        }

        public void setScan(String scan) {
            this.scan = scan;
        }

        public String getCs() {
            return cs;
        }

        public void setCs(String cs) {
            this.cs = cs;
        }

        public String getEnrichment() {
            return enrichment;
        }

        public void setEnrichment(String enrichment) {
            this.enrichment = enrichment;
        }

        public String getProteins() {
            return proteins;
        }

        public void setProteins(String proteins) {
            this.proteins = proteins;
        }

        public String getProteinDescription() {
            return proteinDescription;
        }

        public void setProteinDescription(String proteinDescription) {
            this.proteinDescription = proteinDescription;
        }

        public void setDisplayChroData(DisplayChroData displayChroData) {
            this.displayChroData = displayChroData;
        }

        public DisplayChroData getDisplayChroData() {
            return displayChroData;
        }

        public String getSingleton() {
            return singleton;
        }

        public void setSingleton(String singleton) {
            this.singleton = singleton;
        }
/*

        public String getProteinSequence() {
            return proteinSequence;
        }

        public void setProteinSequence(String proteinSequence) {
            this.proteinSequence = proteinSequence;
        }
*/

        public String getPtmIndexInProtein() {
            return ptmIndexInProtein;
        }

        public void setPtmIndexInProtein(String ptmIndexInProtein) {
            this.ptmIndexInProtein = ptmIndexInProtein;
        }
    }

    public static class MovingAverage {
        private final Queue<Double> window = new LinkedList<Double>();
        private final int period;
        private double sum;

        public MovingAverage(int period) {
            assert period > 0 : "Period must be a positive integer";
            this.period = period;
        }

        public void newNum(double num) {
            sum += num;
            window.add(num);
            if (window.size() > period) {
                sum -= window.remove();
            }
        }

        public double getAvg() {
            if (window.isEmpty()) return 0; // technically the average is undefined
            return sum / window.size();
        }


    }



}
