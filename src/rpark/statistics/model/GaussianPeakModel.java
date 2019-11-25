/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rpark.statistics.model;

/**
 *
 * @author rpark
 * @author rohan
 */
public class GaussianPeakModel {
    private double lowBound;
    private double highBound;
    private double y;
    private double x;
    private double sigma;
    private double sigmaIonInjectCorrection;
    private double yIonInjectionCorrection;
    private double xIonInjectionCorrection;
    private double peakArea=-1;
    private double peakAreaIonInjectionCorrection=-1;
    private int peakStartIndex=-1;
    private int peakEndIndex=-1;
    private String chroData;
    private double[] peakArr;
    private double[] retArr;
    private int[] scanArr;
    private double[] gaussianXArr=null;
    private double[] gaussianYArr=null;
    private boolean hasPeak;
    private int chargeState;
    private boolean peakDetected=true;
    private double[] isoArr;
    //public final static double SQRT_PI=1.772454;  //sqrt(PI)
    public final static double SQRT_PI_60=106.3472;  //sqrt(PI) x 60 seconds

  private double maxIntensity=0;

    public GaussianPeakModel() {

    }

    public double getSigmaIonInjectCorrection() {
        return sigmaIonInjectCorrection;
    }

    public void setSigmaIonInjectCorrection(double sigmaIonInjectCorrection) {
        this.sigmaIonInjectCorrection = sigmaIonInjectCorrection;
    }

    public double getyIonInjectionCorrection() {
        return yIonInjectionCorrection;
    }

    public void setyIonInjectionCorrection(double yIonInjectionCorrection) {
        this.yIonInjectionCorrection = yIonInjectionCorrection;
    }

    public double getxIonInjectionCorrection() {
        return xIonInjectionCorrection;
    }

    public void setxIonInjectionCorrection(double xIonInjectionCorrection) {
        this.xIonInjectionCorrection = xIonInjectionCorrection;
    }


    public GaussianPeakModel(double low, double high) {
        lowBound = low;
        highBound = high;
    }
    public void setLowBound(double low) {
        lowBound = low;
    }
    public void setHighBound(double high) {
        highBound = high;
    }
    public double getLowBound() {
        return lowBound;
    }
    public double getHighBound() {
        return highBound;
    }
    public boolean isInRange(double value) {
        return value >= lowBound && value <= highBound;
    }

    public double getY() {
        return y;
    }

    public void setY(double y) {
        this.y = y;
    }

    public double getX() {
        return x;
    }

    public void setX(double x) {
        this.x = x;
    }

    public double getSigma() {
        return sigma;
    }

    public void setSigma(double sigma) {
        this.sigma = sigma;
    }

    public double getPeakArea() {

      //check if peak modeled max value is bigger than max real intensity, discard it

      if(this.y>maxIntensity*1.5) {
          this.peakArea = maxIntensity * this.sigma * SQRT_PI_60;
      }
        //peakArea = 0;
      else if(peakArea<0) {

          //full width at half maximum(FWHM).  If FWHM is bigger than peak height, we assume peak modeling fails.  return zero.
          //double fwhm = this.sigma * 2.35;

          //if(fwhm>y*2) peakArea = 0;
          //else
            this.peakArea = this.y * this.sigma * SQRT_PI_60;

          //System.out.println(">>>------" + y + " " + sigma + " " + peakArea);

        }

        return peakArea;
    }

    public double getPeakAreaTargeted() {

        if(peakArea<0)
        this.peakArea = this.y * this.sigma * SQRT_PI_60;


        return peakArea;
    }


    public void setPeakArea(double peakArea) {
        this.peakArea = peakArea;
    }

    public double getPeakAreaIonInjectionCorrection() {
        if(peakAreaIonInjectionCorrection<0)
            this.peakAreaIonInjectionCorrection = this.yIonInjectionCorrection*this.sigmaIonInjectCorrection*SQRT_PI_60;
        return peakAreaIonInjectionCorrection;
    }

    public void setPeakAreaIonInjectionCorrection(double peakAreaIonInjectionCorrection) {
        this.peakAreaIonInjectionCorrection = peakAreaIonInjectionCorrection;
    }


    public int getPeakStartIndex() {
        return peakStartIndex;
    }

    public void setPeakStartIndex(int peakStartIndex) {
        this.peakStartIndex = peakStartIndex;
    }

    public int getPeakEndIndex() {
        return peakEndIndex;
    }

    public void setPeakEndIndex(int peakEndIndex) {
        this.peakEndIndex = peakEndIndex;
    }

    public String getChroData() {
        return chroData;
    }

    public void setChroData(String chroData) {
        this.chroData = chroData;
    }

    public double[] getPeakArr() {
        return peakArr;
    }

    public void setPeakArr(double[] peakArr) {
        this.peakArr = peakArr;
    }

    public double[] getRetArr() {
        return retArr;
    }

    public void setRetArr(double[] retArr) {
        this.retArr = retArr;
    }

    public int[] getScanArr() {
        return scanArr;
    }

    public void setScanArr(int[] scanArr) {
        this.scanArr = scanArr;
    }

    public double[] getGaussianXArr() {
        return gaussianXArr;
    }

    public void setGaussianXArr(double[] gaussianXArr) {
        this.gaussianXArr = gaussianXArr;
    }

    public double[] getGaussianYArr() {
        return gaussianYArr;
    }

    public void setGaussianYArr(double[] gaussianYArr) {
        this.gaussianYArr = gaussianYArr;
    }

    public String getGaussianPeakString() {

        if(gaussianXArr == null) return "";

        StringBuffer sb = new StringBuffer();
        for(int i=0;i<gaussianXArr.length;i++)
            sb.append(gaussianXArr[i]).append(" ").append((int)gaussianYArr[i]).append(";");

        return sb.toString();
    }

    public boolean isHasPeak() {
        return hasPeak;
    }

    public void setHasPeak(boolean hasPeak) {
        this.hasPeak = hasPeak;
    }

    public int getChargeState() {
        return chargeState;
    }

    public void setChargeState(int chargeState) {
        this.chargeState = chargeState;
    }

    public boolean isPeakDetected() {
        return peakDetected;
    }

    public void setPeakDetected(boolean peakDetected) {
        this.peakDetected = peakDetected;
    }

    public double[] getIsoArr() {
        return isoArr;
    }

    public void setIsoArr(double[] isoArr) {
        this.isoArr = isoArr;
    }

  public double getMaxIntensity() {
    return maxIntensity;
  }

  public void setMaxIntensity(double maxIntensity) {
    this.maxIntensity = maxIntensity;
  }
}

