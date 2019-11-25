
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
import edu.scripps.pms.util.spectrum.*;

public class Precursor {

    //private List<Daughter> dlist = new ArrayList<Daughter>();
    private Hashtable<Double, Daughter> ht = new Hashtable<Double, Daughter>();
    private List<Double> dmassList = new ArrayList<Double>();

    private double mass;
    

    public Precursor(double d) {
        this.mass = d;

    }

    public void addDaughter(double dmass) {

        ht.put(dmass, new Daughter(dmass));
        this.dmassList.add(dmass);
    }

    public double getMass() {
        return mass;
    }

    public void setMass(double mass) {
        this.mass = mass;
    }

    //public void addIntensity(double dmass, double intensity, int scan) {
    public void addIntensity(ArrayList<Peak> peakList, int scan, double retTime) {

	Set<Double> set = new HashSet<Double>();

        //System.out.println("-----" + peakList.size());
	for(Iterator<Peak> itr=peakList.iterator(); itr.hasNext(); ) {
	    Peak p = itr.next();

	    Daughter d = ht.get(p.getM2z());            
            d.addIntensity(p.getIntensity(), scan, retTime);
	    set.add(p.getM2z());
	}

	for(Iterator<Double> itr=ht.keySet().iterator(); itr.hasNext(); ) {
	    double dmass = itr.next();
	    if(!set.contains(dmass)) {
		Daughter d = ht.get(dmass);
		d.addIntensity(0, scan, -1);
	    }
	}
    }

    public Hashtable<Double, Daughter> getHt() {
        return ht;
    }

    public Daughter getDaughter(double dmass) {

        return this.ht.get(dmass);
    }

    public Daughter getDaughter(int index) {


        return this.ht.get(this.dmassList.get(index));
    }


    public void setHt(Hashtable<Double, Daughter> ht) {
        this.ht = ht;
    }

    public List<Double> getDmassList() {
        return dmassList;
    }

    public void setDmassList(List<Double> dmassList) {
        this.dmassList = dmassList;
    }

   
    
    public int findpeak() {
        int daughterSize=this.getDmassList().size();
	int intArrSize=0;

	for(int i=0;i<getDmassList().size();i++) {
	    Daughter d = getDaughter(i);

	    if(intArrSize<d.getIntensityList().size())
		intArrSize = d.getIntensityList().size();
	}

        int[][] intArr = new int[daughterSize][intArrSize];

        int daughterIndex=0;

	for(int i=0;i<getDmassList().size();i++) {
	    Daughter d = getDaughter(i);

	    List<Double> intList = d.getIntensityList();
	    for(int j=0;j<intList.size();j++) {
		intArr[daughterIndex][j] = (int)intList.get(j).doubleValue();

	    }

	    daughterIndex++;
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

	return peakIndex;
    }


    
}
