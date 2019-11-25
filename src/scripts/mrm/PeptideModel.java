
/*
* Copyright (c) 2008 Integrated Proteomics Applications.  All rights reserved.  
*/

package scripts.mrm;

/**
 *
 * @author Sung Kyu, Robin, Park
 * @email robinparky@yahoo.com
 * Created on Feb 1, 2010 
 * $Revision:$
 * $Date:$
 */
import java.util.*;

public class PeptideModel {

    private List<Precursor> plist = new ArrayList<Precursor>();
    private String name;

    public PeptideModel(String name) {
        this.name = name;
    }

    public void addPrecursor(Precursor precursor) {
        this.plist.add(precursor);
    }

    public List<Precursor> getPlist() {
        return plist;
    }

    public void setPlist(List<Precursor> plist) {
        this.plist = plist;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public void calculateRatio() {
        if(plist.size() < 2)
            return;

        Precursor p1 = plist.get(0);
        Precursor p2 = plist.get(1);

        int intArrSize=0;
        for(int i=0;i<p1.getDmassList().size();i++) {
            Daughter d1 = p1.getDaughter(i);
            Daughter d2 = p2.getDaughter(i);

            if(intArrSize<d1.getIntensityList().size())
                intArrSize = d1.getIntensityList().size();

            if(intArrSize<d2.getIntensityList().size())
                intArrSize = d2.getIntensityList().size();
        }

        int[] intSum = new int[intArrSize];


        for(int i=0;i<p1.getDmassList().size();i++) {

            Daughter d1 = p1.getDaughter(i);
            Daughter d2 = p2.getDaughter(i);

            List<Double> intList = d1.getIntensityList();

            int count=0;
            for(Iterator<Double> intItr=intList.iterator(); intItr.hasNext(); ) {



                intSum[count++] += (int)intItr.next().doubleValue();

                //if(count==1)
                  //  System.out.println("==>>" + intSum[0]);
            }

            intList = d2.getIntensityList();

            System.out.println(intList.size());

            count=0;
            for(Iterator<Double> intItr=intList.iterator(); intItr.hasNext(); ) {
                intSum[count++] += (int)intItr.next().doubleValue();

            //    if(count==1)
              //      System.out.println("=="+ intSum[0]);
            }

        }


        /*
        for(int i:intSum) {
            System.out.println(i);
        }
        */
    }

    /*
    public int findpeak() {
        int intArrSize=0;
        int daughterSize=0;

        for(Iterator<Precursor> itr=this.plist.iterator(); itr.hasNext(); ) {
            Precursor pre = itr.next();
            daughterSize += pre.getDmassList().size();

            for(int i=0;i<pre.getDmassList().size();i++) {
                Daughter d = pre.getDaughter(i);

                if(intArrSize<d.getIntensityList().size())
                    intArrSize = d.getIntensityList().size();
            }
        }

        int[][] intArr = new int[daughterSize][intArrSize];

        int daughterIndex=0;

        for(Iterator<Precursor> itr=this.plist.iterator(); itr.hasNext(); ) {
            Precursor pre = itr.next();
            daughterSize += pre.getDmassList().size();

            for(int i=0;i<pre.getDmassList().size();i++) {
                Daughter d = pre.getDaughter(i);

                List<Double> intList = d.getIntensityList();
                for(int j=0;j<intList.size();j++) {
                    intArr[daughterIndex][j] = (int)intList.get(j).doubleValue();

                }

                daughterIndex++;
            }
        }

	int sum=0;
	int peakIndex=0;

        for(int i=0;i<intArr[0].length;i++) {

	    int tmpSum=0;

            for(int j=0;j<intArr.length;j++) {
		tmpSum += intArr[j][i];
            }


	    if(sum<tmpSum) {
		sum = tmpSum;
		peakIndex = i;
	    }
        }

	//System.out.println("==" + peakIndex);
	return peakIndex;
    }
    */

    public void print() {

        int intArrSize=0;
        int daughterSize=0;

        for(Iterator<Precursor> itr=this.plist.iterator(); itr.hasNext(); ) {
            Precursor pre = itr.next();
            daughterSize += pre.getDmassList().size();

            for(int i=0;i<pre.getDmassList().size();i++) {
                Daughter d = pre.getDaughter(i);

                if(intArrSize<d.getIntensityList().size())
                    intArrSize = d.getIntensityList().size();
            }
        }

        int[][] intArr = new int[daughterSize][intArrSize];

        int daughterIndex=0;

        for(Iterator<Precursor> itr=this.plist.iterator(); itr.hasNext(); ) {
            Precursor pre = itr.next();
            daughterSize += pre.getDmassList().size();

            for(int i=0;i<pre.getDmassList().size();i++) {
                Daughter d = pre.getDaughter(i);


                List<Double> intList = d.getIntensityList();
                for(int j=0;j<intList.size();j++) {
                    intArr[daughterIndex][j] = (int)intList.get(j).doubleValue();

                }

                daughterIndex++;
            }
        }

        for(int i=0;i<intArr[0].length;i++) {

	    System.out.print(i + ":\t");
            for(int j=0;j<intArr.length;j++) {
                System.out.print(intArr[j][i] + "\t");
            }


            System.out.println("");
            
        }
        
    }

    public void print(int index) {

        Precursor p = plist.get(index);
        List<Double> l = p.getDmassList();

        for(Iterator<Double> mitr=l.iterator(); mitr.hasNext(); ) {
            double dmass = mitr.next();
            Daughter daughter = p.getDaughter(dmass);

            List<Double> intList = daughter.getIntensityList();

            for(Iterator<Double> iitr=intList.iterator(); iitr.hasNext(); ) {
                double intensity = iitr.next();
                System.out.println(p.getMass() + "\t" + dmass + "\t" + intensity);
            }

        }

    }

}
