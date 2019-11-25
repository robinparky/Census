/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rpark.statistics;

import java.math.BigDecimal;
import java.util.*;

/**
 *
 * @author Rohan Rampuria <rampuria@scripps.edu>
 */
public class BHCorrection {

  //get p values and get corresponding q values.  comprese
    public static void main (String [] args){
        List<Double> pValueOrigList = new ArrayList();

        pValueOrigList.add(0.7);
        pValueOrigList.add(0.8);
        pValueOrigList.add(0.85);
        pValueOrigList.add(0.850000001);
        pValueOrigList.add(0.8500000001);
        pValueOrigList.add(0.9);
        pValueOrigList.add(0.9);
        pValueOrigList.add(0.9);

        pValueOrigList.add(0.0);

        pValueOrigList.add(0.0);
        pValueOrigList.add(0.5);
        pValueOrigList.add(0.5);
        pValueOrigList.add(0.5);
        pValueOrigList.add(0.5);
        pValueOrigList.add(0.6);
        pValueOrigList.add(0.0);
        List<Double> qValueOrigList = BHCorrection.runBHCorrectionRemoveRedundant(pValueOrigList);

        for(int i=0; i<pValueOrigList.size(); i++)
        {
            System.out.println(pValueOrigList.get(i)+"\t"+qValueOrigList.get(i));
        }



        /*
        try {
            BufferedReader br = new BufferedReader(new FileReader(args[0]));
    //      BufferedReader br = new BufferedReader(new FileReader("/ip2_data4/rpark/libsearch/DDA/pvalue1.txt"));
            String line;
            List<Double> pvalueList = new ArrayList<>();
          String lastLine = null;
            while((line=br.readLine())!=null)
            {
              pValueOrigList.add(line);

              if(line.equals(lastLine)) continue;

              pvalueList.add(Double.parseDouble(line));
            }

            List<Double> qvalues =runBhCorrection(pvalueList);



            BufferedWriter bw = new BufferedWriter(new FileWriter(args[1]));

          //for(String each:pValueOrigList)
          for(int i=0;i<pValueOrigList.size();i++)
          {
            bw.write(pValueOrigList.get(i) +"\t"+qvalues.get(i));
            bw.newLine();
          }

            bw.close();
            br.close();
        }
        catch(IOException ioe)
        {

        }
/*
        pValues.add(0.2678678678E10);
        pValues.add(0.64556665656986767E9);
        pValues.add(0.4894768678676677667E11);

        System.out.println(runBhCorrection(pValues));*/
    }


    public static class ValueContainer
    {
        public final double pvalue;
        public final int index;
        private double qvalue=-1;
        private boolean toEvaluate  =false;
        private ValueContainer refContainer = null;

        public ValueContainer(double pvalue, int index) {
            this.pvalue = pvalue;
            this.index = index;
        }

        public double getPvalue() {
            return pvalue;
        }

        public int getIndex() {
            return index;
        }

        public double getQvalue() {
            if(refContainer!=null)
                return  refContainer.qvalue;
            return qvalue;
        }

        public void setQvalue(double qvalue) {
            this.qvalue = qvalue;
        }

        public boolean isToEvaluate() {
            return toEvaluate;
        }

        public void setToEvaluate(boolean toEvaluate) {
            this.toEvaluate = toEvaluate;
        }


        public void setRefContainer(ValueContainer refContainer) {
            this.refContainer = refContainer;
        }
        public String toString()
        {
            return pvalue+"\t"+index+"\t"+toEvaluate+"\t"+qvalue;
        }
    }



    public static class IndexComparator implements Comparator<Integer>
    {
        private double[] array;
        private Integer index[];
        public IndexComparator(double [] array)
        {
            this.array = array;
            index = new Integer[array.length];
            for(int i=0; i<index.length; i++)
            {
                index[i] = i;
            }
        }


        @Override
        public int compare(Integer integer, Integer t1) {
            return Double.compare(array[integer],array[t1]);
        }

        public Integer[] getIndex() {
            return index;
        }
    }

    public static List<Double> runBHCorrectionRemoveRedundant(List<Double> pValues)
    {
        List<ValueContainer> containerList = new ArrayList<>();
        for(int i=0; i<pValues.size(); i++)
        {
            ValueContainer valueContainer = new ValueContainer(pValues.get(i),i);
            containerList.add(valueContainer);
        }
        Collections.sort(containerList, new Comparator<ValueContainer>() {
            @Override
            public int compare(ValueContainer valueContainer, ValueContainer t1) {
                return Double.compare(valueContainer.getPvalue(),t1.getPvalue());
            }
        });
        ValueContainer refContainer = containerList.get(0);
        refContainer.setToEvaluate(true);
        for(int i=1; i<containerList.size(); i++)
        {
            ValueContainer currentContainer = containerList.get(i);
            if(Double.compare(currentContainer.getPvalue(),refContainer.getPvalue()) != 0 || currentContainer.getPvalue() == 0)
            {
                currentContainer.setToEvaluate(true);
                refContainer = currentContainer;
            }
            else
            {
                currentContainer.setToEvaluate(false);
                currentContainer.setRefContainer(refContainer);
            }
        }
        BigDecimal min = new BigDecimal("" + 1);
        BigDecimal mkprk;
        final int RESULT_SCALE = 10;
        int len = containerList.size();
        int count = containerList.size();
        for (int i = containerList.size(); i > 0; i--) {
            if (containerList.get(i - 1).isToEvaluate()) {
                mkprk = (new BigDecimal("" + len).multiply(new BigDecimal(containerList.get(i - 1).getPvalue())))
                        .divide(new BigDecimal("" + count), RESULT_SCALE, BigDecimal.ROUND_HALF_UP);
                if (mkprk.compareTo(min) < 0) {
                    min = mkprk;
                }
                count--;
                containerList.get(i - 1).setQvalue(min.doubleValue());
            }

        }
        Collections.sort(containerList, new Comparator<ValueContainer>() {
            @Override
            public int compare(ValueContainer valueContainer, ValueContainer t1) {
                return valueContainer.index - t1.index;
            }
        });
        List<Double> BHcorrection = new ArrayList<>();
        for(int i=0; i<containerList.size(); i++)
        {
            BHcorrection.add(containerList.get(i).getQvalue());
        }
        return BHcorrection;
    }



    public static List<Double> runBhCorrection(List<Double> pValues) {
        int len = pValues.size();
        double[] orderedPValues = new double[len];
        double[] adjustedpValues = new double[len];
        List<Double> BhCorrection = new ArrayList<>();
        int[] indexOfValues = new int[len];
        final int RESULT_SCALE = 10;

        for (int i = 0; i < len; i++) {
            orderedPValues[i] = (double) pValues.get(i);
        }
        // sort the values
        java.util.Arrays.sort(orderedPValues);
        for (int i = 0; i < len; i++) {
            indexOfValues[i] = getIndexOf(orderedPValues, (double) pValues.get(i));
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
        }
        // adjust the sequence
        len = pValues.size();
        int j = 0;
        for (int i = 0; i < len; i++) {
            try {
                double tmp = (double) pValues.get(i);

                if (tmp > -100000.0) {
                    //NumberFormat formatter = new DecimalFormat("##.#####");
                    //String apvString = formatter.format(adjustedpValues[indexOfValues[j]]);
                    //                      proteinList.get(i).setPostHocP(apvString);
                    BhCorrection.add(adjustedpValues[indexOfValues[j]]);
                    j++;
                } else {
                    //                      proteinList.get(i).setPostHocP("-1");
                    BhCorrection.add(-1.0);
                }
            } catch (Exception e) {
                BhCorrection.add(-1.0);
            }
        }
        //System.out.println(""+BhCorrection);
        return BhCorrection;
    }

    private static int getIndexOf(double[] orderedPvalues, double value) {
        int index = -1;
        for (int i = 0; i < orderedPvalues.length; i++) {
            if (orderedPvalues[i] == value) {
                return i;
            }
        }
        return index;
    }

}
