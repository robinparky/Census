/*          
 * 
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.scripps.pms.census.compare;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.math.BigDecimal;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.math.MathException;
import org.apache.commons.math.stat.inference.TestUtils;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

/**
 *
 * @author Harshil
 */
public class ProteinDetails
{
   private List locus = new ArrayList();
   private List specCount = new ArrayList();
   private List pepCount = new ArrayList();
   private List experimentGroup = new ArrayList();
   private List descriptions = new ArrayList();
   private String headerLine;
  private List experimentGroupNameIndexs = new ArrayList();
  private List experimentGroupAttribute = new ArrayList();
  private List pValues = new ArrayList();
  private List fValues = new ArrayList();
  private List stdDev = new ArrayList();
  private List average = new ArrayList();
  private List rsdValue = new ArrayList();
  private List variancePvalueBH = new ArrayList();
  
// Index
   private int locusIndex;
   private int specCountIndex;
   private int pepCountIndex;
   private int descriptionIndex;
   private   int startExperimentIndex;
   private     int endExperimentIndex;
 

    
   void ProteinDeails()
   {
       
   }
   public void add(ProteinDetails deletedEntries)
   {
       for(int i=0;i<deletedEntries.locus.size();i++)
       {
             locus.add(deletedEntries.locus.get(i));
               pepCount.add(deletedEntries.pepCount.get(i));
               specCount.add(deletedEntries.specCount.get(i));
               experimentGroup.add(deletedEntries.experimentGroup.get(i));
               descriptions.add(deletedEntries.descriptions.get(i));
//               pValues.add(deletedEntries.pValues.get(i));
//               fValues.add(deletedEntries.fValues.get(i));
               pValues.add(-1);
               fValues.add(-1);
               stdDev.add(deletedEntries.stdDev.get(i));
               average.add(deletedEntries.average.get(i));
               rsdValue.add(deletedEntries.rsdValue.get(i));
               variancePvalueBH.add(-1);
       }
   }
   public void removeNaNEntries(ProteinDetails deletedEntries)
   {
       for(int i=0;i<pValues.size();i++)
       {
           if(Double.isNaN((double) pValues.get(i)))
           {
               deletedEntries.locus.add(locus.get(i));
               deletedEntries.pepCount.add(pepCount.get(i));
               deletedEntries.specCount.add(specCount.get(i));
               deletedEntries.experimentGroup.add(experimentGroup.get(i));
               deletedEntries.descriptions.add(descriptions.get(i));
               deletedEntries.pValues.add(pValues.get(i));
               deletedEntries.fValues.add(fValues.get(i));
               deletedEntries.stdDev.add(stdDev.get(i));
               deletedEntries.average.add(average.get(i));
               deletedEntries.rsdValue.add(rsdValue.get(i));
               
               locus.remove(i);pepCount.remove(i);specCount.remove(i);
               experimentGroup.remove(i);
                descriptions.remove(i);
                pValues.remove(i);
                fValues.remove(i);
                stdDev.remove(i);
                average.remove(i);
                rsdValue.remove(i);
                               
           }
       }
   }
    void parseHeaderLine(String header) {
        int flag = 0;
        String words[] = header.split("\t");

        for (int i = 0; i < words.length; i++) 
        {
            if (words[i].equals("LOCUS")) {
                locusIndex = i - 1;
            }
            if (words[i].equals("SPEC_COUNT")) {
                specCountIndex = i - 1;
            }
            if (words[i].equals("PEP_NUM")) {
                pepCountIndex = i - 1;
            }
            if (words[i].equals("DESCRIPTION")) {
                descriptionIndex = i - 1;
            }
            if (flag == 0 && words[i].contains("m/z")) {
                startExperimentIndex = i - 1;
                flag = 1;
            }
            if (flag == 1 && !words[i].contains("m/z")) {
                endExperimentIndex = i - 1;
                flag = -1;
            }
        }
    }
  
    public List getExperimentGroup() {
        return experimentGroup;
    }
    public void removeCurrentEntry()
    {
        locus.remove(locus.size() -1 );
        specCount.remove(specCount.size() -1 );
        pepCount.remove(pepCount.size() -1 );
        descriptions.remove(descriptions.size() -1 );
        
    }
    public void setExperimentGroup(String[] words)
    {
        List proteinGroup= new ArrayList();
//        int zeroFlag=0;
        for(int i=0;i<getExperimentGroupNameIndexs().size();i++)
        {
            int[] indexes = (int[]) getExperimentGroupNameIndexs().get(i);
            double[] currentGroup = new double[indexes.length];
            
            for(int j=0;j<indexes.length;j++ )
            {
                currentGroup[j] =(double) Double.parseDouble(words[indexes[j]]);
//                if(currentGroup[j]!=0)
//                    zeroFlag =1;
            }
            
                proteinGroup.add(currentGroup);
        }
//        if(zeroFlag==0)
//                removeCurrentEntry();
//        else
            experimentGroup.add(proteinGroup);
       
    }

   
    public List getDescriptions() {
        return descriptions;
    }

    public String getHeaderLine() {
        return headerLine;
    }

    public void setHeaderLine(String headerLine) {
        this.headerLine = headerLine;
    }

    public List getLocus() {
        return locus;
    }
    public void setLocus(String[] words) {
        locus.add(words[locusIndex]);
    }

    public List getSpecCount() {
        return specCount;
    }

    public void setSpecCount(String[] words) {
        specCount.add(words[specCountIndex]);
    }

    public List getPepCount() {
        return pepCount;
    }

    public void setPepCount(String[] words) {
       pepCount.add(words[pepCountIndex]);
    }

    public void setDescriptions(String[] words) {
       descriptions.add(words[descriptionIndex]);
    }

    public void setExperimentGroupNameIndexs(List experimentGroupNames)
    {
        for(int i=0;i<experimentGroupNames.size();i++)
        {
            double[] currentGroup = (double[])(experimentGroupNames.get(i));
            int[] currentGroupIndex = new int[currentGroup.length];
            String words[]=headerLine.split("\t");
            for(int k=0;k< currentGroup.length;k++)
            {
                for(int j=startExperimentIndex;j<=endExperimentIndex;j++)
                {
                  if(words[j].contains(Double.toString(currentGroup[k])) &&  words[j].contains("norm"))
                  {
                      currentGroupIndex[k]=j-1;
                      break;
                  }
                }
            }
            experimentGroupNameIndexs.add(currentGroupIndex);
                               
        }
    }
    
    public void generateAnovaValues() {
        List classes = getExperimentGroup();
        try {
            for (int i = 0; i < classes.size(); i++) 
            {
                fValues.add(TestUtils.oneWayAnovaFValue((List) classes.get(i)));
                pValues.add(TestUtils.oneWayAnovaPValue((List) classes.get(i)));
             }
        } catch (IllegalArgumentException ex) {
            Logger.getLogger(ProteinDetails.class.getName()).log(Level.SEVERE, null, ex);
        } catch (MathException ex) {
            Logger.getLogger(ProteinDetails.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

    public List getExperimentGroupAttribute() {
        return experimentGroupAttribute;
    }

    public void setExperimentGroupAttribute(List experimentGroupAttribute) {
        this.experimentGroupAttribute = experimentGroupAttribute;
    }

    public List getStdDev() {
        return stdDev;
    }

    public void setStdDev() 
    {
        for(int i=0;i< experimentGroup.size();i++)
        {
            List CurrentGroup = (List) experimentGroup.get(i);
            List localStdev = new ArrayList();
            List localAvg = new ArrayList();
            List localRsd = new ArrayList();
            for(int j=0;j<CurrentGroup.size();j++)
            {
                DescriptiveStatistics stdevValue = new DescriptiveStatistics();
                for(int z = 0;z<((double[])CurrentGroup.get(j)).length ;z++)
                {
                    stdevValue.addValue((double) ((double[])CurrentGroup.get(j))[z]);
                }
                 localStdev.add(stdevValue.getStandardDeviation());
                 localAvg.add(stdevValue.getMean());
                 localRsd.add(Math.sqrt(stdevValue.getVariance()) * 100 / stdevValue.getMean());
            }
            stdDev.add(localStdev);
            average.add(localAvg);
            getRsdValue().add(localRsd);
        }
    }

    public List getExperimentGroupNameIndexs() {
        return experimentGroupNameIndexs;
    }

    public List getpValues() {
        return pValues;
    }

    public List getfValues() {
        return fValues;
    }

    public List getAverage() {
        return average;
    }

    public List getRsdValue() {
        return rsdValue;
    }
   
    
     public void setVariancePvalueBH()
    {
         int len = pValues.size();
        double[] orderedPValues = new double[len];
        double[] adjustedpValues = new double[len];
        int[] indexOfValues = new int[len];
        final int RESULT_SCALE = 10;

        for (int i = 0; i < len; i++) {
            orderedPValues[i] = (double) pValues.get(i);
            System.out.println(pValues.get(i) + "  " + i);
        }
        // sort the values
        java.util.Arrays.sort(orderedPValues);
        for (int i = 0; i < len; i++) {
//                              indexOfValues[i]=getIndexOf(orderedPValues, pValueDoubles[i]);
            indexOfValues[i] = getIndexOf(orderedPValues, (double)pValues.get(i));
//                              indexOfValues[i]=getIndexOf(pValueDoubles, orderedPValues[i]);
        }
                        // calculate the post hoc adjustment
                        BigDecimal min = new BigDecimal("" + 1);
                        BigDecimal mkprk;
                        for (int i = len; i > 0; i--) {
                         
                                mkprk = (new BigDecimal("" + len).multiply(new BigDecimal(orderedPValues[i - 1]))).divide(new BigDecimal("" + i), RESULT_SCALE, BigDecimal.ROUND_HALF_UP);
                                if (mkprk.compareTo(min) < 0) {
                                        min = mkprk;
                                }
                          
                                adjustedpValues[i - 1] = min.doubleValue();
//                                System.out.println( (i-1)+" " +adjustedpValues[i-1]);
                        }
                       
                        // adjust the sequence
                        len = pValues.size();
                        int j = 0;
                        for(int i=0; i<len; i++){
                                try {
                                        //      double tmp = Double.parseDouble(proteinList.get(i).getLogRatioVariance2());
//                                     if(Double.isNaN((double) pValues.get(i)))
//                                     {
//                                           getVariancePvalueBH().add(-1.0);
//                                           continue;
//                                     }
                                        double tmp = (double) pValues.get(i);
 
                                        //if(tmp>-100000.0){
                                        NumberFormat formatter = new DecimalFormat("##.#####");
                                        String apvString = formatter.format(adjustedpValues[indexOfValues[j]]);
                //                      proteinList.get(i).setPostHocP(apvString);
                                        getVariancePvalueBH().add(adjustedpValues[indexOfValues[j]]);
                                        j++;
                                                /*
                                        }
                                        else{
                        //                      proteinList.get(i).setPostHocP("-1");
                                                getVariancePvalueBH().add(-1.0);
                                        }*/
                                } catch(Exception e) {
                                        getVariancePvalueBH().add(-1.0);
                                }
                        }
        
    }
     public static int getIndexOf(double[] orderedPvalues, double value)
     {
         int index = -1;
         for(int i=0;i<orderedPvalues.length;i++)
         {
             if(orderedPvalues[i]== value)
             {
                return i;
             }
         }
         return index;
     }

    /**
     * @return the variancePvalueBH
     */
    public List getVariancePvalueBH() {
        return variancePvalueBH;
    }
    public static void main(String args[]) throws FileNotFoundException, IOException
    {
        ProteinDetails pd = new ProteinDetails();
        List a= new ArrayList();
//              a.add(0.12185);a.add(0.74764);a.add(0.32332);
        File f = new File ("C:\\Users\\Harshil\\Documents\\NetBeansProjects\\a1.csv");
        BufferedReader br = new BufferedReader(new FileReader(f));
        String line = br.readLine();
        while(line!= null)
        {
           try{
            if(!line.contains("NA"))
//                a.add((double) 0);
                
//            else 
                a.add(Double.parseDouble(line));
           }
           catch(Exception e)
                   {
//                       a.add((double) 0);
                   }
            line = br.readLine();
        }
        br.close();
        pd.setVariancePvalueBH(a);
//        System.out.println(pd.getVariancePvalueBH());
    }
    
     public void setVariancePvalueBH(List pValues)
    {
         int len = pValues.size();
        double[] orderedPValues = new double[len];
        double[] adjustedpValues = new double[len];
        int[] indexOfValues = new int[len];
        final int RESULT_SCALE = 10;

        for (int i = 0; i < len; i++) {
            orderedPValues[i] = (double) pValues.get(i);
        }
        // sort the values
        java.util.Arrays.sort(orderedPValues);
        for (int i = 0; i < len; i++) {
//                              indexOfValues[i]=getIndexOf(orderedPValues, pValueDoubles[i]);
            indexOfValues[i] = getIndexOf(orderedPValues, (double)pValues.get(i));
//                              indexOfValues[i]=getIndexOf(pValueDoubles, orderedPValues[i]);
        }
                        // calculate the post hoc adjustment
                        BigDecimal min = new BigDecimal("" + 1);
                        BigDecimal mkprk;
                        for (int i = len; i > 0; i--) {
                                mkprk = (new BigDecimal("" + len).multiply(new BigDecimal(orderedPValues[i - 1]))).divide(new BigDecimal("" + i), RESULT_SCALE, BigDecimal.ROUND_HALF_UP);
                                if (mkprk.compareTo(min) < 0) {
                                        min = mkprk;
                                }
                                adjustedpValues[i - 1] = min.doubleValue();
//                                System.out.println( (i-1)+" " +adjustedpValues[i-1]);
                        }
                        //System.exit(0);
                        // adjust the sequence
                        len = pValues.size();
                        int j = 0;
                        for(int i=0; i<len; i++){
                                try {
                                        //      double tmp = Double.parseDouble(proteinList.get(i).getLogRatioVariance2());
                                        double tmp = (double) pValues.get(i);
 
                                        if(tmp>-100000.0){
                                                NumberFormat formatter = new DecimalFormat("##.#####");
                                                String apvString = formatter.format(adjustedpValues[indexOfValues[j]]);
                        //                      proteinList.get(i).setPostHocP(apvString);
                                                getVariancePvalueBH().add(adjustedpValues[indexOfValues[j]]);
                                                j++;
                                        }
                                        else{
                        //                      proteinList.get(i).setPostHocP("-1");
                                                getVariancePvalueBH().add(-1.0);
                                        }
                                } catch(Exception e) {
                                        getVariancePvalueBH().add(-1.0);
                                }
                        }
        
    }
}
