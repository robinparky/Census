/*
 * NonLabelMappingModel.java
 *
 * Created on August 24, 2006, 11:43 AM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package edu.scripps.pms.census.model;

import java.util.*;
import java.io.*;
import edu.scripps.pms.census.ChroGenerator;

import gnu.trove.TIntDoubleHashMap;

/**
 *
 * @author rpark
 */
public class NonLabelMappingModel extends Hashtable<String, Hashtable> {
    
    private Vector<String> pathFileNameList;
    private int[][][] pathArray;
    private int refIndex;
    
    //contains maximum scan index
    private Hashtable<String, Integer> fileNameHt = new Hashtable<String, Integer>();
    private Hashtable<String, Hashtable> ms2ms1Ht= new Hashtable<String, Hashtable>();
    private Hashtable<String, Hashtable> ms1ms2Ht= new Hashtable<String, Hashtable>();

    private int[] ms1KeyArr = null;
    private int[] ms2KeyArr = null;
            
    public Hashtable<String, Hashtable> getMs2ms1Ht()
    {
	return ms2ms1Ht;
    }

    /** Creates a new instance of NonLabelMappingModel */
    public NonLabelMappingModel(int[][][] pathArray, Vector<String> pathFileNameList, int refIndex) {
       
        this.pathFileNameList = pathFileNameList;
        this.pathArray = pathArray;
        this.refIndex = refIndex;        
        
        init();
    }

    public Hashtable<Integer, Hashtable> getParentHt(String fileName)
    {
        return this.get(fileName);
    }
    
    public int getMaxScanIndex(String fileName)
    {
        return this.fileNameHt.get(fileName);        
    }
        
    public Hashtable getChildHashtable(String fileName, int scanIndex)
    {           
        return this.getChildHashtable(fileName, Integer.valueOf(scanIndex));
    }    
    
    public Hashtable getChildHashtableByMS2(String ms2FileName, Integer scanIndex)
    {   
     
        Hashtable<Integer, Integer> ht = this.ms2ms1Ht.get(ms2FileName);

        /*
	 if(scanIndex==16)
         {
           System.out.println("-->>" + ht.keySet() + " " +ms2FileName + " " + scanIndex);
           System.out.println("-->>" + ht.get(scanIndex));
         }
        */
        
        Integer indexValue = ht.get(scanIndex);
        
        if(null == indexValue)
            return null;
        
        int index = indexValue;
        
       // System.out.println(scanIndex + " " + index);
        
     //   if(null == this.getChildHashtable(ms2FileName, Integer.valueOf(index)))
            //System.out.println(this.get(ms2FileName).keySet());
       
	// return ms1 mapping info
//	ms2FileName = ms2FileName.replace(".ms2", ".ms1");
        return this.getChildHashtable(ms2FileName, Integer.valueOf(index));
    }
        
    public Hashtable getChildHashtable(String fileName, Integer scanIndex)
    {           
        return (Hashtable)this.get(fileName).get(scanIndex);        
    }

    public void reinitializeMaxIndexByRet()
    {
        for(int i=0;i<pathArray.length;i++)
	{
	    String tempFileName = pathFileNameList.get(i);

            int tempIndex = 0;
            
            
            for(int ii=pathArray[i][0].length-1;ii>=0;ii--)
            {
                if(pathArray[i][1][ii] > 0)
                {
                    tempIndex = pathArray[i][1][ii];
                    break;
                }
            }

            fileNameHt.put(tempFileName, tempIndex);
        }
                
    }

    
    public void setMsmsMap(Hashtable<String, edu.scripps.pms.census.hash.IndexedFile> ht1, Hashtable<String, edu.scripps.pms.census.hash.IndexedFile> ht2) throws java.io.IOException
    {        
        for(Iterator<String> itr=ht1.keySet().iterator(); itr.hasNext(); )
        {
            String ms1File = itr.next();
            
	    String ms2File = ms1File.replace(".ms1", ".ms2");
            
            edu.scripps.pms.census.hash.IndexedFile eachIndex2 = ht2.get(ms2File);
	    if(null == eachIndex2)
	    {
		System.out.println("Error: cannot find " + ms2File);
		continue;
	    }

            this.fileNameHt.put(ms2File, ht2.get(ms2File).getKeys().length-1);
            
	    Hashtable ht = new Hashtable();
	    Hashtable<Integer, ArrayList> ms1PrecurHt = new Hashtable<Integer, ArrayList>();

            edu.scripps.pms.census.hash.IndexedFile eachIndex1 = ht1.get(ms1File);
            ms1KeyArr = eachIndex1.getKeys();
            ms2KeyArr = eachIndex2.getKeys();

	    TIntDoubleHashMap preMap = eachIndex2.getPrecursorMap();

	    int ms2Index=0;

	    //move to the first valid msms scan#
	    while(true)
	    {
		if(ms2KeyArr[ms2Index]>ms1KeyArr[0])
		    break;

		ms2Index++;
	    }
	  
            for(int i=0;i<ms1KeyArr.length-1;i++)
            {
//		System.out.println(ms1KeyArr[i] + " " + ms1KeyArr[i+1]);

//                int val1 = Integer.parseInt(ms1KeyArr[i].toString());
  //              int val2 = Integer.parseInt(ms1KeyArr[i+1].toString());
   
		while(true)
		{
		    if(ms2Index>=ms2KeyArr.length)
			break;

		    int value = ms2KeyArr[ms2Index];

		    if(ms1KeyArr[i]<value && ms1KeyArr[i+1]>value)
		    {
			//ht.put(value, ms1KeyArr[i]);
			ht.put(ms2Index, i); //save index instead of scan number

			ArrayList list = ms1PrecurHt.get(i);
			if(null == list)
			{
			    list = new ArrayList();
			    list.add(ms2Index);
			    ms1PrecurHt.put(i, list);
			}
			else
			    list.add(ms2Index);

//			ms1PrecurHt.put("" + i+preMap.get(ms2Index), ms2Index);
//			System.out.println(i+preMap.get(ms2Index), ms2Index);

			//System.out.println(ms2KeyArr[ms2Index] + " " + preMap.get(ms2KeyArr[ms2Index]));
		    }
		    else
			break;
			
		    ms2Index++;
		}		
            }

//	    System.out.println(ht);
	    ms2ms1Ht.put(ms2File, ht);
	    ms1ms2Ht.put(ms2File, ms1PrecurHt);
        }
        //System.exit(0);
//	    System.out.println(ms2ms1Ht);

    }
 
    public Integer getMs2Index(String ms2FileName, int ms1Index, double precursor, double massTolerance, edu.scripps.pms.census.hash.IndexedFile iFile)
    {
        
        //System.out.println("------------->>" + ms2FileName + " " + ms1Index + " " + precursor + " " + massTolerance  + " " + iFile);
        
	Hashtable<Integer, ArrayList> tmpHt = ms1ms2Ht.get(ms2FileName);


//	for(Iterator itr= tmpHt.keySet().iterator(); itr.hasNext(); )
//	    System.out.println( itr.next() );

        //if(null == tmpHt.get(ms1Index))
//	System.out.println(ms1Index + " " + tmpHt.keySet());
        
	int keys[] = iFile.getKeys();

        while(true)
        {
            ArrayList aList = tmpHt.get(ms1Index);
            
            for(Iterator<Integer> itr = aList.iterator(); itr.hasNext(); )
            {
                int eachIndex = itr.next();
          //      System.out.println("--->>" + keys[eachIndex] +  " " + eachIndex);

                double d = iFile.getPrecursorMap().get(keys[eachIndex]);

                if( precursor>(d-massTolerance) && precursor<(d+massTolerance) )
                {
                    return eachIndex;
                }
            }
            
            ms1Index--;
            
            if(ms1Index<0)
                break;
            
        }
        
        //System.out.println(ms1KeyArr[ms1Index] + " ");
        
        //for(int i=0;i<ms1KeyArr.length;i++)
          //  System.out.print(ms1KeyArr[i] + " ");
        
        
        /*
        System.out.println("------------->>" + ms1ms2Ht.get(ms2FileName));
        System.out.println("------------->>" + ms1Index + precursor);
        
        ArrayList tmpList = (ArrayList)ms1ms2Ht.get(ms2FileName).get(new Integer(ms1Index));
        
        for(Iterator<Integer> itr=tmpList.iterator(); itr.hasNext(); )
        {
            Integer each = itr.next();
            
            System.out.println(iFile.getPrecursorMap().get(keys[each]));
        }
        

	return (Integer)ms1ms2Ht.get(ms2FileName).get("" + ms1Index + precursor);
         */
        
        
        return -1;
    }
    
    private void init()
    {
        //Hashtable<String, Hashtable> masterHt = new Hashtable<String, Hashtable>();
        
        for(Iterator itr=pathFileNameList.iterator(); itr.hasNext(); )
        {
            this.put(itr.next().toString(), new Hashtable());
        }
        
	Hashtable<Integer, Hashtable> refHt = new Hashtable<Integer, Hashtable>();
        int refMinIndex = 2;  //minumum is always zero
        int refMaxIndex = 0; //pathArray[0][0][pathArray[0][0].length-1];        

	for(int i=0;i<pathArray.length-1;i++)
	{
	    String tempFileName = null;

	    if(i>=refIndex)
	    {
		tempFileName = pathFileNameList.get(i+1);

		int tempIndex = 0;
		for(int ii=0;ii<pathArray[i][0].length;ii++)
		{
		    if(0 != pathArray[i][1][ii])
		    {
			tempIndex = pathArray[i][1][ii];
			break;
		    }
		}

		fileNameHt.put(tempFileName, tempIndex);
                
                //System.out.println("111==" + tempFileName + " " + tempIndex + " " + i + " " + refIndex);
	    }
	    else
	    {
		tempFileName = pathFileNameList.get(i);    

		int tempIndex = 0;
		for(int ii=0;ii<pathArray[i][0].length;ii++)
		{
		    if(0 != pathArray[i][1][ii])
		    {
			tempIndex = pathArray[i][1][ii];
			break;
		    }
		}

		fileNameHt.put(tempFileName, tempIndex);
	    }

	    for(int j=0;j<pathArray[i][0].length;j++)
	    {
		Hashtable<String, Set> refChildHt = refHt.get(pathArray[i][0][j]);
                
                if(refMaxIndex<pathArray[i][0][j])                
                    refMaxIndex = pathArray[i][0][j];

		if( null == refChildHt )
		{
		    refChildHt = new Hashtable<String, Set>();
		    Set<Integer> tempSet = new HashSet<Integer>();
		    tempSet.add(pathArray[i][1][j]);
		    refChildHt.put( tempFileName,  tempSet);	

                    refHt.put(pathArray[i][0][j], refChildHt);
		}
		else
		{
		    Set<Integer> tempSet = refChildHt.get(tempFileName);

		    if(null == tempSet)
			tempSet = new HashSet<Integer>();                                      

		    tempSet.add(pathArray[i][1][j]);
		    refChildHt.put( tempFileName,  tempSet);	    

		}

		Set refChildSet = refChildHt.get(pathFileNameList.get(refIndex));

		if(null == refChildSet)
		{
		    refChildSet = new HashSet();
		    refChildSet.add(pathArray[i][0][j]);

		    refChildHt.put(pathFileNameList.get(refIndex), refChildSet);
		}

	    }

	}

	fileNameHt.put(pathFileNameList.get(refIndex), refMaxIndex);

	this.put(pathFileNameList.get(refIndex), refHt);
	
	for(Iterator<Integer> itr = refHt.keySet().iterator(); itr.hasNext(); )
	{
	    Integer refScanNumIndex = itr.next();

	    Hashtable<String, Set> refChildHt = (Hashtable)refHt.get( refScanNumIndex ); //retrieve hashtable by ref scan # index

	    for(Iterator<String> fItr=pathFileNameList.iterator(); fItr.hasNext(); )
	    {
		String tmpFileName = fItr.next();

		if( tmpFileName.equals(pathFileNameList.get(refIndex)) )
		    continue;

		Set tmpSet = refChildHt.get(tmpFileName);

		if(null == tmpSet)
		    continue;

		//each target hashtable
		Hashtable targetHt = this.get(tmpFileName);

		for(Iterator setItr=tmpSet.iterator(); setItr.hasNext(); )
		{
		    int ii = Integer.parseInt( setItr.next().toString() );

		    Hashtable<String, Set> targetChildHt = (Hashtable)targetHt.get(ii);

		    if(null == targetChildHt)
			targetChildHt = new Hashtable();

		    //for(Iterator<String> fItr2=pathFileNameList.iterator(); fItr2.hasNext(); )

		    for(int fIndex=0;fIndex<pathFileNameList.size();fIndex++)
		    {
			String tmpFile2 = pathFileNameList.get(fIndex); //fItr2.next();

			Set mappingSet = refChildHt.get(tmpFile2);
			Set targetSet = (Set)targetChildHt.get(tmpFile2);

			if(null == targetSet)
			    targetSet = new HashSet();

                        //if mappingSet is null, search closest scan index in refht
			if(null == mappingSet)
                        {
                            int tempIdx = refScanNumIndex.intValue()-1;
                            
                            while(true)
                            {
                                Hashtable<String, Set> tmpRefChildHt = (Hashtable)refHt.get( tempIdx ); 
                                
                                //well... if tempIdx becomes 2, tmpRefChildHt should not be null
                                if(null != tmpRefChildHt) { // || tempIdx<=2) {  

                                    mappingSet = tmpRefChildHt.get(tmpFile2);

				    if(null != mappingSet)
				    {
					targetSet.addAll(mappingSet);
					targetChildHt.put(tmpFile2, targetSet);   
				    }
                        
                                    break;
                                }                                    
                                
                                tempIdx--;                                
                            }
                            
			    continue;
                        }

			targetSet.addAll(mappingSet);
			targetChildHt.put(tmpFile2, targetSet);   

		    }

		    targetHt.put(ii, targetChildHt);
		}

		//if(true) continue;

	    }

	}
    }

    public Vector<String> getPathFileNameList() {
        return pathFileNameList;
    }

    public void setPathFileNameList(Vector<String> pathFileNameList) {
        this.pathFileNameList = pathFileNameList;
    }

    public int[][][] getPathArray() {
        return pathArray;
    }

    public void setPathArray(int[][][] pathArray) {
        this.pathArray = pathArray;
    }

    public int getRefIndex() {
        return refIndex;
    }

    public void setRefIndex(int refIndex) {
        this.refIndex = refIndex;
    }

    public Hashtable<String, Integer> getFileNameHt() {
        return fileNameHt;
    }
}
