/*
 * LinearRegression2.java
 *
 * Created on July 25, 2005, 11:09 PM
 */

package scripts.emily;

/**
 *
 * @author rpark
 */
public class LinearRegression2 {
    
    private double[] xArr;    
    private double[] yArr;
    private int startIndex;
    private int endIndex;
    private int shift=0;
    
    private double corr=-1000;
    private int bestShift;

    private double slope;
    private double intercept;
    private double xAve;
    private double yAve;

    private double areaRatio;
    
    
    /** Creates a new instance of LinearRegression2 */
    public LinearRegression2(double[] yArr, double[] xArr, int startIndex, int endIndex, int shift) {
        this.xArr = xArr;
        this.yArr = yArr;
        this.startIndex = startIndex;
        this.endIndex = endIndex;
        this.shift = shift;
        
        calc();
    }
    
    private void calc()
    {
        startIndex = (startIndex>=0)?startIndex:0;
        endIndex = (endIndex<xArr.length)?endIndex:xArr.length-1;
        int size = xArr.length;        
        double curCorr;
        int dataEntry = endIndex - startIndex + 1;
        
        for(int move=-shift;move<=shift;move++)
        {
            double xSum=0;
            double ySum=0;
            double xxSum=0;
            double yySum=0;
            double xySum=0;
            for(int i=startIndex;i<=endIndex;i++)
            {
                int temp = i+move;
                if(temp<0)
                    temp=0;
                else if(temp>=size)
                    temp = size-1;
                
                xSum += xArr[i];
                ySum += yArr[temp];
                xxSum = xxSum + (double)xArr[i]*xArr[i];
                yySum += (double)yArr[temp]*yArr[temp];
                xySum += (double)xArr[i]*yArr[temp];
            }

            curCorr = (xySum-(double)xSum*ySum/dataEntry)/Math.sqrt((xxSum-(double)xSum*xSum/dataEntry)*(yySum-(double)ySum*ySum/dataEntry));

            if(this.corr<curCorr)
            {
                this.corr = curCorr;
                this.bestShift = move;
                
                slope = (xySum-(double)xSum*ySum/dataEntry)/(xxSum-xSum*xSum/dataEntry);
                yAve = ySum/dataEntry;
                xAve = xSum/dataEntry;

                intercept = yAve - slope*xAve;
                
                this.areaRatio = ySum/xSum;
            }
        }
    }

    public double getCorr() {
	if(corr<0 || corr>1)
	    return -1;

        return corr;
    }

    public void setCorr(double corr) {
        this.corr = corr;
    }

    public int getBestShift() {
        return bestShift;
    }

    public void setBestShift(int bestShift) {
        this.bestShift = bestShift;
    }

    public double getSlope() {
        return slope;
    }

    public void setSlope(double slope) {
        this.slope = slope;
    }

    public double getIntercept() {
        return intercept;
    }

    public void setIntercept(double intercept) {
        this.intercept = intercept;
    }

    public double getAreaRatio() {
        return areaRatio;
    }

    public void setAreaRatio(double areaRatio) {
        this.areaRatio = areaRatio;
    }
    
}
