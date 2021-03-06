/**
 * @file PeakList.java
 * This is the source file for edu.scripps.pms.util.spectrum.PeakList
 * @author Tao Xu
 * @date $Date: 2014/06/07 00:01:38 $
 */



package edu.scripps.pms.util.spectrum;

import edu.scripps.pms.mspid.MassSpecConstants;
import java.util.ArrayList;
import java.util.ListIterator;
import java.util.Iterator;
import java.util.List;
import java.util.Comparator;
import java.util.Collections;

public class PeakList {
    public static final String RETENTIONTIME = "I\tRetTime";
    public static final String INJECTIONTIME = "I\tIonInjectionTime";
    public static int DEFAULTNUMPEAKS = 1000;
    public static final int DEFAULTSPECTRUMSIZE = 100000;   
    public static boolean SORTBYM2Z = false;
    public static boolean SORTBYINTENSITY = true;
    private double totalIntensity = 0;
    private int loscan;
    private int hiscan;
    private double precursorMass; // this is actually precursor m/z, not M+H
    private List<String> hlines;
    private List<String> ilines = new ArrayList<String>();
    private ArrayList<Zline> zlines = new ArrayList<Zline>();
    private ArrayList<Peak> peaks = new ArrayList<Peak>(DEFAULTNUMPEAKS);
    private String listType;
    private int index = 0;
    private double retTime =-1;

    private ArrayList<Peak> peaksSortedByIntensity = null;
    private ArrayList<Peak> peaksSortedByM2z = null;

    public PeakList()
    {
    }

    public PeakList(double retentime)
    {
        this.retTime = retentime;
    }

    
    public int getNumZlines() {
        return zlines.size();
    }
    public Peak getPeak(double m2z, double massTolerance) {
        int index = getIndex(m2z, massTolerance);
        if(index > -1) {
            return getSortedPeaks(SORTBYM2Z).get(index);
        } else {
            return null;
        }
    }
    

    // return -1 if the peak not fund
    private int getIndex(double m2z, double massTolerance) {
        ArrayList<Peak> sortedPeaks = getSortedPeaks(SORTBYM2Z);
        int index = -1;
        int max = peaks.size() - 1;
        int min = 0;
        int mid = max/2;
        double maxM2z = m2z + massTolerance;
        double minM2z = m2z - massTolerance;

        while(max > (min+1)) {
            double m = sortedPeaks.get(mid).getM2z();
            if(m >= minM2z && m <= maxM2z) {
                return mid;
            } 
           
            if(m > minM2z) {
                max = mid;
                mid = (min + mid)/2;
            } else if(m < minM2z) {
                min = mid;
                mid = (mid + max)/2;
            }
        } 
        return index;
    }
//    public double getRetentionTime() {
//        double retentionTime = -1;
//        for(Iterator<String> it = ilines.iterator(); it.hasNext();) {
//            String i = it.next();
//            if(i.startsWith(RETENTIONTIME)) {
//                String [] contents = i.split("\t");
//                retentionTime = Double.parseDouble(contents[2]);
//                break;
//            }
//        
//        }
//        return retentionTime;
//    }
    public double getRetentionTime() {
        double retentionTime = this.retTime;
        for(Iterator<String> it = ilines.iterator(); it.hasNext();) {
            String i = it.next();
            if(i.startsWith(RETENTIONTIME)) {
                String [] contents = i.split("\t");
                retentionTime = Double.parseDouble(contents[2]);
                break;
            }

        }
        return retentionTime;
    }
    
    public double getInjectionTime() {
        double inejctionTime = -1;
        for(Iterator<String> it = ilines.iterator(); it.hasNext();) {
            String i = it.next();
            if(i.startsWith(INJECTIONTIME)) {
                String [] contents = i.split("\t");
                inejctionTime = Double.parseDouble(contents[2]);
                break;
            }
        
        }
        return inejctionTime;
    }
    
    public void setHlines(List<String> hlines) {
       this.hlines = hlines;
    }
    
    public Peak getPeakByM2z(double m2z) {
        double minDiff = 1000;
        Peak prcPeak = null;
        for(Iterator<Peak> it = peaks.iterator(); it.hasNext();) {
            Peak currPeak = it.next();
            double m2zCurrent = currPeak.getM2z(); 
            double diff = m2z > m2zCurrent? m2z-m2zCurrent : m2zCurrent-m2z; // avoid abs()
            if(diff == 0) {
                return currPeak;
            }
            if(diff < minDiff) {
                minDiff = diff;
                prcPeak = currPeak;
            }
        }
        return prcPeak;
    }
    public int getPrecursorScan() {
        for(Iterator<String> it = ilines.iterator(); it.hasNext();) {
            String s = it.next();
            if(s.startsWith("I\tPrecursorScan")) {
                String [] elements = s.split("\t");
                return Integer.parseInt(elements[2]);
            }
        }
        return -1; // I line for PrecursorScan was not found
    }
    public void setZlines(ArrayList<Zline> zlines) {
       this.zlines = zlines;
    }
    public void setLoscan(int loscan) {
        this.loscan = loscan;
    }
    public void setHiscan(int hiscan) {
        this.hiscan = hiscan;
    }
    public List<String> getIlines() {
        return ilines;
    }
    public List<String> getHlines() {
        return hlines;
    }
    public void setPrecursorMass(double precursorMass) {
        this.precursorMass = precursorMass;
    }
    public void addIline(String l) {
        ilines.add(l);
    }
    public void addZline(Zline z) {
        zlines.add(z);
    }

    public void addPeak(Peak p) {
        p.setIndex(index++);
        peaks.add(p);
        totalIntensity += p.getIntensity();
    }

    public int numPeaks() {
        return peaks.size();
    }
    public double getTotalIntensity() {
        return totalIntensity;
    }
    public Iterator<Zline> getZlines() {
        return zlines.listIterator();
    }
    public void sortPeaks(boolean sortByIntensity) {
        peaks = getSortedPeaks(sortByIntensity);
    }

    /**
     * Return a sorted list (incremental by m2z or intensity)
     * for the Peaks in this PeakList
     * @param sortByIntensity - indicate how the list should be sorted,
     *                          true for sort by intensity,
     *                          false for sort by M2z
     * @note user must not modify the List returned. 
     * 
     */
    public synchronized ArrayList<Peak> getSortedPeaks(boolean sortByIntensity) {
        ArrayList<Peak> sortedPeaks =  new ArrayList(DEFAULTNUMPEAKS);
        if (sortByIntensity) {
            if (peaksSortedByIntensity != null) {
                return peaksSortedByIntensity;
            } else {
                peaksSortedByIntensity = sortedPeaks;
            }
        } else {
            if (peaksSortedByM2z != null) {
                return peaksSortedByM2z;
            } else {
                peaksSortedByM2z = sortedPeaks;
            }
        }
        for (Peak p : peaks) {
            sortedPeaks.add(p);
        }

        Collections.sort(sortedPeaks, new PeakComparator(sortByIntensity));
        return sortedPeaks; 
    }
    public synchronized ArrayList<Peak> getSortedPeaks(int numPeaks, boolean sortByIntensity) {
        ArrayList topPeaks = new ArrayList<Peak>(numPeaks);
 
        List<Peak> sortedList = getSortedPeaks(sortByIntensity);
        int totalPeaks = sortedList.size()-1; 
        while(numPeaks > 0 && totalPeaks >= 0) {
            Peak p = sortedList.get(totalPeaks--);
            //if(p.getM2z() > getMaxChargeState()*getPrecursorMass()){
            topPeaks.add(p);
            numPeaks--;
             
        }
        return topPeaks;
    }
    public ListIterator<Peak> getPeaks() {
        return peaks.listIterator();
    }

    public String getListType()
    {
        return listType;
    }

    public void setListType(String listType)
    {
        this.listType = listType;
    }
    /**
     * Return the low scan number of this
     */
    public int getLoscan() {
        return loscan;
    }

    public int getHiscan() {
        return hiscan;
    }

    // return the m/z of the precursor ion, not the M+H value
    public double getPrecursorMass() {
        return precursorMass;
    }
    public double getMaxPrecursorMass() {
        double maxPrecMass = 0;
        for(Zline z : zlines) {
            double precMass = z.getM2z();
            maxPrecMass = maxPrecMass > precMass? maxPrecMass : precMass;
        }
        return maxPrecMass;
    } 
    public double getMaxM2z() {
        List<Peak> sortedPeaks = getSortedPeaks(SORTBYM2Z);
        return sortedPeaks.get(numPeaks()-1).getM2z();
    } 
    public double getMinM2z() {
        List<Peak> sortedPeaks = getSortedPeaks(SORTBYM2Z);
        return sortedPeaks.get(0).getM2z();
    } 
    public double getMaxIntensity() {
        List<Peak> sortedPeaks = getSortedPeaks(SORTBYINTENSITY);
        return sortedPeaks.get(sortedPeaks.size()-1).getIntensity();
    }
    public double getMinIntensity() {
        List<Peak> sortedPeaks = getSortedPeaks(SORTBYINTENSITY);
        return sortedPeaks.get(0).getIntensity();
    }
    public ArrayList<Peak> getLeastIntensePeaks(int percentage) {
        List<Peak> sortedPeaks = getSortedPeaks(SORTBYINTENSITY);
//System.out.println("NumPeaks: " + peaks.size() + "\tnumLeastIntensePeaks: " + numPeaks);
        int numPeaks = (int)(peaks.size()*percentage/100 + 0.5);
//System.out.println("NumPeaks: " + peaks.size() + "\tnumLeastIntensePeaks: " + numPeaks);
        numPeaks = numPeaks > 0? numPeaks : 1;
        ArrayList<Peak> leastIntensePeaks = new ArrayList<Peak>(numPeaks); 
        int lastIndex = peaks.size() - 1;
        for(int i = 0; i < numPeaks; i++) {
//System.out.println(sortedPeaks.get(i).getIntensity());
            leastIntensePeaks.add(sortedPeaks.get(i));
        }
        return leastIntensePeaks;
    }
    // assume the peaks are sorted by m2z
    public double getMaxIntensity(double minM2z, double maxM2z) {
        double maxIntensity = -10000;
        for(Iterator<Peak> it = peaks.iterator(); it.hasNext();) {
            Peak p = it.next();
            double intensity = p.getIntensity();
            double m2z = p.getM2z();
            if(m2z >= minM2z) {
                if(m2z <= maxM2z) {
                    maxIntensity = maxIntensity > intensity? maxIntensity : intensity;
                } else {
                    break;
                }
            }
        }
        return maxIntensity;    
    }
    public double getAvgIntensityOfLeastIntensePeaks(int percentage) {
        List<Peak> sortedPeaks = getSortedPeaks(SORTBYINTENSITY);
//System.out.println("NumPeaks: " + peaks.size() + "\tnumLeastIntensePeaks: " + numPeaks);
        int numPeaks = (int)(peaks.size()*percentage/100 + 0.5);

        numPeaks = numPeaks > 0? numPeaks : 1;
//System.out.println("NumPeaks: " + peaks.size() + "\tnumLeastIntensePeaks: " + numPeaks);
        //int lastIndex = peaks.size() - 1;
        double totalIntens = 0;
        for(int i = 0; i < numPeaks; i++) {
            totalIntens += sortedPeaks.get(i).getIntensity();
//System.out.println(i + "\t" + sortedPeaks.get(i).getIntensity() +"\ttotalIntens: " + totalIntens + "\tnumPeaks: " + numPeaks);
        }
//System.out.println("totalIntensity: " + totalIntens + "\tumPeaks: " + numPeaks + "\tavgNoise: " + totalIntens/numPeaks);
        return totalIntens/numPeaks;
    }
    public String getSpectrumWithHlines() {
        StringBuffer sb = new StringBuffer(DEFAULTSPECTRUMSIZE);
        for(String h : hlines) {
            sb.append(h);
            sb.append("\n");
        }
        getSpectrumWithoutHlines(sb); 
        return sb.toString();       
    }
    public String getSpectrumWithoutHlines() {
        StringBuffer sb = new StringBuffer(DEFAULTSPECTRUMSIZE);
        getSpectrumWithoutHlines(sb);
        return sb.toString();        
    }
    private void getSpectrumWithoutHlines(StringBuffer sb) {
        addSline(sb);
        addIlines(sb);
        addZlines(sb);
        addSpectrum(sb);
    }
    private void addIlines(StringBuffer sb) {
        for(String i : ilines) {
            sb.append(i);
            sb.append("\n");
        }
    }
    private void addSline(StringBuffer sb) {
        sb.append("S");
        sb.append("\t");
        sb.append(loscan);
        sb.append("\t");
        sb.append(hiscan);
        sb.append("\t");
        sb.append(precursorMass);
        sb.append("\n");
    }
    private void addZlines(StringBuffer sb) {
        for(Zline z : zlines) {
            sb.append("Z");
            sb.append("\t");
            sb.append(z.getChargeState());
            sb.append("\t");
            sb.append(z.getM2z());
            sb.append("\n");
            for(String d : z.getDlines()) {
                sb.append(d);
                sb.append("\n");
            }
        }
    }
    private void addSpectrum(StringBuffer sb) {
        for(Peak p : peaks) {
            sb.append(p.getM2z());
            //sb.append("\t");
            sb.append(" ");
            sb.append(p.getIntensity());
            sb.append("\n");
        }        
    }
    public int getMaxChargeState() {
        int maxChargeState = 0;
        for(Zline z : zlines) {
            if(z.getChargeState() > maxChargeState) {
                maxChargeState = z.getChargeState();
            }
        }
        return maxChargeState;
    } 
    private List<Peak> getMostIntensPeaks(int numPeaks, double lowM2zLimit, double highM2zLimit) {
        ArrayList<Peak> topPeaks = new ArrayList<Peak>(numPeaks);
 
        List<Peak> sortedList = getSortedPeaks(true);
        int totalPeaks = sortedList.size()-1; 
        while(numPeaks > 0 && totalPeaks >= 0) {
            Peak p = sortedList.get(totalPeaks--);
            //if(p.getM2z() > getMaxChargeState()*getPrecursorMass()){
            double m2z = p.getM2z();
            if(m2z > lowM2zLimit && m2z < highM2zLimit){
                topPeaks.add(p);
                numPeaks--;
            }
             
        }
        return topPeaks;
    }
    public PointList [] calcQCorrs(int numPeaks, int accuracyFactor, double massWindow) {
        int massH = (int)(MassSpecConstants.MASSH*accuracyFactor + 0.5);
        PointList [] points = new PointList[getNumZlines()];
        int zlineCounter = 0;
        
        for(Zline z : zlines) {
            int chargeState = z.getChargeState();
            double lowLimit = massWindow*chargeState;
            double noShiftMass = precursorMass*chargeState;
            double highLimit = noShiftMass - lowLimit;
            List<Peak> topList = getMostIntensPeaks(numPeaks, lowLimit, highLimit);
            int numBins = (int)(precursorMass*chargeState*accuracyFactor + 0.5);
            double [] intens = new double[numBins];    
            for(Peak p : topList) {
                int index = (int)(p.getM2z()*accuracyFactor+0.5);
                intens[index] += p.getIntensity();
            } 
            //double [] temp = intens;
            intens = reverse(intens); // get the revered list for correlation
            int numShift = (int)massWindow*chargeState*accuracyFactor;
            int offSet = numShift/2; 
//System.out.println("scan#: " + getLoscan() +"\tprecuMass: " + precursorMass +  "\toffSet: " + offSet);
            double [] qcorr = new double [numShift];
        
            for(Peak p : topList) {
                int m2z = (int)(p.getM2z()*accuracyFactor + 0.5);
                double intensity = p.getIntensity();  
                int startIndex = m2z-offSet;
                for(int i = 0; i < numShift; i++) {
                    qcorr[i] += intensity*intens[startIndex + i];
                }
            } 
            PointList pointList = new PointList();
            double firstMass = noShiftMass+massWindow/2*chargeState;
            double dAccFactor = accuracyFactor; // convert to accuracyFactor to double
System.out.println("FirstMass: " + firstMass + "\tnoShiftMass: " + noShiftMass+ "lowLimit: " + lowLimit + "\tMax: " + (numBins+offSet)/dAccFactor + "\tMin: " + (numBins-offSet)/dAccFactor);
            for(int i = 0; i < qcorr.length; i++) {
                pointList.addPoint(new Point(firstMass-i/dAccFactor, qcorr[i]));
            }
            points[zlineCounter++] = pointList;
            int indexMax = getMaxIndex(qcorr);
            //double mass = (numBins-offSet-indexMax-1.5*massH)/(accuracyFactor+0.0); 
            double mass = firstMass-indexMax/(accuracyFactor+0.0); 
//System.out.println("chargeState: " + z.getChargeState()+"\tnumBins: " + numBins +"\tm2z: " + z.getM2z() + "\tmass: " + mass + "\tindexMax: " + indexMax + "\tqcorr: " + qcorr[indexMax]);
//System.out.println(z.getChargeState() + "\t"+qcorr[indexMax] + "\t" + (z.getM2z()-mass));
            //String dline = "D\t" + mass + "\t" + qcorr[indexMax];
            //z.addDline(dline);
        } 
        return points;
    }
    private double [] reverse(double [] intens) {
        double [] reversed = new double[intens.length];
        int lastIndex = intens.length - 1;
        for(int i = 0; i <= lastIndex; i++) {
            reversed[lastIndex-i] = intens[i];
        }
        return reversed;
    }
    private int getMaxIndex(double [] qcorr) {
        
        int indexMax = 0;
        double max = 0;
        for(int i = 0; i < qcorr.length; i++) {
            if(qcorr[i] > max) {
                max = qcorr[i];
                indexMax = i;
            }            

        }
//System.out.println("max: " + max);
        return indexMax;
    } 

    public ArrayList<Peak> getPeakList() {
	return peaks;
    }
}
