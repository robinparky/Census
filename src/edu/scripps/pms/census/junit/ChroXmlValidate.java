package edu.scripps.pms.census.junit;

import edu.scripps.pms.census.io.*;
import edu.scripps.pms.census.model.*;
import java.util.*;

public class ChroXmlValidate 
{
    public static void main(String args[]) throws Exception
    {
        ChroXmlReader reader = new ChroXmlReader(args[0]);
        
        List list = reader.getProteinList();

        System.out.print("parsing...");
        for(Iterator<ChroProtein> itr = list.iterator(); itr.hasNext(); )
        {
            ChroProtein protein = itr.next();
            List peptideList = protein.getPeptideList();

            for(Iterator<ChroPeptide> pepItr = peptideList.iterator(); pepItr.hasNext(); )
            {
                ChroPeptide peptide = pepItr.next();

                List l = peptide.getDataList();

                long[] samArr = new long[l.size()];
                long[] refArr = new long[samArr.length];

                ChroData data1 = (ChroData)l.get(0);
                ChroData data2 = (ChroData)l.get(1);
                int scanDiff = data2.getScanNum()-data1.getScanNum();
                for(int i=0;i<l.size()-1;i++)
                {
                    data1 = (ChroData)l.get(i);
                    data2 = (ChroData)l.get(i+1);
                    
                    if(scanDiff != (data2.getScanNum()-data1.getScanNum()))
                    {
                        System.out.println(protein.getLocus());
                        System.out.println(peptide.getSequence());
                        System.out.println("data : " + data1.getScanNum() + " " + data2.getScanNum());

                        System.out.println("Data seems to be wrong");

                    }

                }
            }
        }

        System.out.println("\nvalidation completed");
    }

    
}
