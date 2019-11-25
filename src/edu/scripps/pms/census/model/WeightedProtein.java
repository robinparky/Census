/*
 * WeightedProtein.java
 *
 * Created on October 10, 2006, 5:39 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package edu.scripps.pms.census.model;

import java.util.*;
import org.apache.commons.lang3.ArrayUtils;
import gnu.trove.TDoubleArrayList;
/**
 *
 * @author rpark
 */
public class WeightedProtein {
    
    
    private Hashtable<String, ProteinModel> ht = new Hashtable<String, ProteinModel>();
    
    
    /** Creates a new instance of WeightedProtein */
    public WeightedProtein() {
        //add(stdev, intAvg, sampleName);
    }
    
    public void add(double stdev, long intAvg, String sampleName) {
        ProteinModel model = ht.get(sampleName);
        
        if(null == model)
        {
            ht.put(sampleName, new ProteinModel(stdev, intAvg));
        }
        else
        {
            model.add(stdev, intAvg);
        }
    }
    
    public double getWeightedAverage(String sampleName)
    {
        ProteinModel model = ht.get(sampleName);
        return model.getWeightedAverage();
        
    }
    
    
    public static class ProteinModel {
        
        private ArrayList<Double> stdevList = new ArrayList<Double>();
        private ArrayList<Long> intAvgList = new ArrayList<Long>();
        private ArrayList<Double> ratioList= new ArrayList<Double>();
        private TDoubleArrayList enrichRatioList= new TDoubleArrayList();
        

        private long totalIntensity=0;
        private double totalRatio=0;

        public TDoubleArrayList getEnrichRatioList() {
            return enrichRatioList;
        }

        public void setEnrichRatioList(TDoubleArrayList enrichRatioList) {
            this.enrichRatioList = enrichRatioList;
        }
    
        public void addEnrichRatio(double d) {
            this.enrichRatioList.add(d);
        }
        
        public double getAverageEnrichRatio(boolean useOutlier, double pValueThreshold) {
            
            if(pValueThreshold<0)
                pValueThreshold=0.1;
                
            double sum = 0;
            
            if(useOutlier && this.enrichRatioList.size()>=3) {
                
                double[] logArr = new double[enrichRatioList.size()];
                for(int i=0;i<enrichRatioList.size();i++) {
                    logArr[i] = Math.log(enrichRatioList.get(i))/Math.log(2);    
                    
                  //  System.out.println("========" + enrichRatioList.get(i) + " " + logArr[i]);
                }
                
                logArr = edu.scripps.pms.stats.GrubbsTest.filterAndRemove(logArr, pValueThreshold);
                
                
                //Math.exp(logArr[li]*Math.log(2));
                
                for(double d:logArr) {
                                        
                    sum += d;
                //    System.out.println("========" + d);
                }
                
                double avg = sum / logArr.length;
                avg = Math.exp(avg*Math.log(2));
                
              //  System.out.println("--------:" + avg);
                return avg;
                
                
                // ArrayUtils.toObject(enrichRatioList);
            } else {
                
                for(int i=0;i<enrichRatioList.size();i++) {
                    double d = enrichRatioList.get(i);
                    sum += d;
                }
                
                return sum/this.enrichRatioList.size();
            }
                      
         //   return 0;
        }
        
	public ProteinModel()
	{
	}

        public ProteinModel(double stdev, long intAvg)
        {
            add(stdev, intAvg);
        }
       
	
        public void add(double stdev, long intAvg)
        {
            totalIntensity += intAvg;
            
            stdevList.add(stdev);
            intAvgList.add(intAvg);
        }
       
        public void add(double stdev, double ratio)
        {
            totalRatio += ratio;
            
            stdevList.add(stdev);
            ratioList.add(ratio);
        }

	//weighed average for standard weight of the labeled
	public double getStandardWeightedAverage()
	{
            double[] weightArr = new double[stdevList.size()];
            
            double averageAll = totalRatio / weightArr.length;

            
            double sum=0;
            
            for(int i=0;i<weightArr.length;i++)
            {
                //this is for weight mean from encyclopedia
                double stdev = stdevList.get(i); 
                
                //relStdev = relStdev>0?relStdev:0.0;
                double temp;
                
                if(stdev>0)
                    temp = 1/(stdev*stdev);
                else
                    temp = 0;
                
                sum += temp;
                weightArr[i] = temp;                
            }
            
            double weightSum = 0;
            double weightAvgSum = 0;
            //System.out.println("==>>" + stdevList.size());
            
            for(int i=0;i<weightArr.length;i++)
            {                
                weightArr[i] /= sum;                                         
                
                weightSum += weightArr[i];                
	 //  System.out.println("-->>" +  weightArr[i] + " " + ratioList.get(i) + " " + sum);
                weightAvgSum += (weightArr[i]*ratioList.get(i));
            }
           
          //  System.out.println("temp..." + sum + " " + weightSum + " " + (weightAvgSum/weightSum));
            
            return weightAvgSum/weightSum;
        }
        
        
	//wikipedia version of weight.  Used for census paper for nonlabeling calculation
        public double getWeightedAverage()
        {
            double[] weightArr = new double[stdevList.size()];
            
            double averageAll = totalIntensity / weightArr.length;
            
            double sum=0;
            
            for(int i=0;i<weightArr.length;i++)
            {
                //relative standard deviation for nomalization
                //double relStdev = stdevList.get(i) / intAvgList.get(i);
                
                //this is for weight mean from encyclopedia
                double relStdev = Math.abs(averageAll-intAvgList.get(i));
                
                
                //relStdev = relStdev>0?relStdev:0.0;
                double temp;
                
                if(relStdev>0)
                    temp = 1/(relStdev*relStdev);
                else
                    temp = 0;
                
                //System.out.println(relStdev + " -----" + stdevList.get(i) + " " + intAvgList.get(i) + " " + temp);
                sum += temp;
                weightArr[i] = temp;                
            }
            
            double weightSum = 0;
            double weightAvgSum = 0;

            for(int i=0;i<weightArr.length;i++)
            {                
                weightArr[i] /= sum;           
                
                weightSum += weightArr[i];                
                weightAvgSum += (weightArr[i]*intAvgList.get(i));
            }

            return weightAvgSum/weightSum;
        }

        public static double getWeightedAverage(double [] list )
        {
            double totalIntensity = 0;
            for(double d: list)
            {
                totalIntensity+=d;
            }
            double[] weightArr = new double[list.length];

            double averageAll = totalIntensity / weightArr.length;

            double sum=0;

            for(int i=0;i<weightArr.length;i++)
            {
                //relative standard deviation for nomalization
                //double relStdev = stdevList.get(i) / intAvgList.get(i);

                //this is for weight mean from encyclopedia
                double relStdev = Math.abs(averageAll-list[i]);


                //relStdev = relStdev>0?relStdev:0.0;
                double temp;

                if(relStdev>0)
                    temp = 1/(relStdev*relStdev);
                else
                    temp = 0;

                //System.out.println(relStdev + " -----" + stdevList.get(i) + " " + intAvgList.get(i) + " " + temp);
                sum += temp;
                weightArr[i] = temp;
            }

            double weightSum = 0;
            double weightAvgSum = 0;

            for(int i=0;i<weightArr.length;i++)
            {
                weightArr[i] /= sum;

                weightSum += weightArr[i];
                weightAvgSum += (weightArr[i]*list[i]);
            }

            return weightAvgSum/weightSum;
        }

        public static double getWeightedAverage(List<Double> list)
        {
            double totalIntensity = 0;
            for(double d: list)
            {
                totalIntensity+=d;
            }
            double[] weightArr = new double[list.size()];

            double averageAll = totalIntensity / weightArr.length;

            double sum=0;

            for(int i=0;i<weightArr.length;i++)
            {
                //relative standard deviation for nomalization
                //double relStdev = stdevList.get(i) / intAvgList.get(i);

                //this is for weight mean from encyclopedia
                double relStdev = Math.abs(averageAll-list.get(i));


                //relStdev = relStdev>0?relStdev:0.0;
                double temp;

                if(relStdev>0)
                    temp = 1/(relStdev*relStdev);
                else
                    temp = 0;

                //System.out.println(relStdev + " -----" + stdevList.get(i) + " " + intAvgList.get(i) + " " + temp);
                sum += temp;
                weightArr[i] = temp;
            }

            double weightSum = 0;
            double weightAvgSum = 0;

            for(int i=0;i<weightArr.length;i++)
            {
                weightArr[i] /= sum;

                weightSum += weightArr[i];
                weightAvgSum += (weightArr[i]*list.get(i));
            }

            return weightAvgSum/weightSum;
        }
    }
}
