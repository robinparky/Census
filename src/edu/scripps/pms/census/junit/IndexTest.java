package edu.scripps.pms.census.junit;

import edu.scripps.pms.census.io.*;
import edu.scripps.pms.census.model.*;
import java.util.*;

import javax.swing.*;
import java.io.*;
import java.awt.Color;
import java.util.*;
import java.text.DecimalFormat;

import edu.scripps.pms.census.util.SimpleFileNameFilter;
import edu.scripps.pms.census.io.ChroReader;
import edu.scripps.pms.census.io.parse.ChroXMLParser;

import edu.scripps.pms.census.util.RelExFileFilter;
import edu.scripps.pms.census.util.PostCalculation;
import edu.scripps.pms.census.util.CalcUtil;

import edu.scripps.pms.census.model.ChroProtein;
import edu.scripps.pms.census.model.ChroPeptide;
import edu.scripps.pms.census.model.ChroData;

import javax.swing.table.DefaultTableModel;

import ptolemy.plot.*;
import ptolemy.plot.plotml.PlotBoxMLParser;
import ptolemy.plot.plotml.PlotMLParser;
import edu.scripps.pms.census.plot.*;
import edu.scripps.pms.census.conf.*;

import edu.scripps.pms.census.util.LinearRegression;
import gnu.trove.TDoubleArrayList;
import gnu.trove.TIntDoubleHashMap;
import gnu.trove.TIntLongHashMap;

import edu.scripps.pms.census.hash.*;

public class IndexTest 
{
    private static String path;
    private Hashtable ht = new Hashtable();
    private static float startRange;
    private static float endRange;
    private static double massDiff;
    private static int numIsoWindows;
    private static double isoWin;
    private static double lastPrecur;
    private static double[] arr;
    private static TDoubleArrayList precursorList = new TDoubleArrayList();

    public static void main(String args[]) throws Exception
    {
        double isoWindow=15;

        Configuration conf = Configuration.getInstance();
        startRange = 700;
        conf.setStartMassRange(startRange);
        endRange = 1210;
        conf.setEndMassRange(endRange);
        conf.setIsolationWindow( isoWindow );
        double startPrecursor = startRange + isoWindow/2;
        for(double i=startPrecursor; i<=endRange; i+=isoWindow)
        {
            precursorList.add(i);
        }

        path = args[0];

        arr = precursorList.toNativeArray();
        conf.setPrecursorArr(arr);

        
        numIsoWindows = conf.getNumOfIsolationWindow();
        isoWin = conf.getIsolationWindow();
        lastPrecur = arr[arr.length-1];

        massDiff = isoWin - conf.getMassRange(); //this is negative value

        for( Iterator fileItr = getAllFiles(args[0]); fileItr.hasNext(); )
        {
            String eachFile = fileItr.next().toString();

            if(eachFile.endsWith("index"))
            {
                System.out.println("reading file " + eachFile + "...");

                IndexTest test = new IndexTest(eachFile);
                System.out.println("");

                //              System.out.println(reader.getHeader());
            }
        }

    }

    public IndexTest(String fileName) throws IOException, Exception
    {

        MSIndexBuilder builder = new MSIndexBuilder(path + File.separator + fileName);
        TIntLongHashMap posMap = builder.readIndexFile();
        TIntDoubleHashMap precursorMap = builder.getPrecursorMap();

        System.out.println(posMap.size());
        System.out.println(precursorMap.size());

        if(posMap.size() != precursorMap.size())
            System.out.println("Error : size is wrong");


    /*
        BufferedReader br = new BufferedReader(new FileReader(path + File.separator + fileName));

        String curLine;

        String lastLine = br.readLine();
        String[] str;
        str = lastLine.split("\t");

        int prevSNum = Integer.parseInt(str[0]);
        double prevPrecur = Double.parseDouble(str[2]);
        double curPrecur;
        posMap.put(prevSNum, Long.parseLong(str[1]));
        precursorMap.put(prevSNum, prevPrecur);

        while( (lastLine=br.readLine()) != null)
        {
            str = lastLine.split("\t");

            curPrecur = Double.parseDouble(str[2]);
            ht.put(str[0], str[2]);

            double curDiff = curPrecur-prevPrecur;

            int curSNum = Integer.parseInt(str[0]);

            int scanDiff = curSNum-prevSNum;
            long curPos = Long.parseLong(str[1]);
            double mDiff = curPrecur - prevPrecur;

            if( scanDiff==1 || (scanDiff==2 && curDiff==massDiff) )
                //if( massDiff==isoWin || (scanDiff==2 && curDiff==massDiff) )
            {
                posMap.put(curSNum, curPos);
                precursorMap.put(curSNum, curPrecur);
            }
            else if(scanDiff>numIsoWindows || mDiff<0)
            {
                for(double i=prevPrecur+isoWin;i<=lastPrecur;i+=isoWin)
                {
                    prevSNum++;



                    posMap.put(prevSNum, -1);
                    precursorMap.put(prevSNum, i); //prevPrecur+isoWin*(i+1));
                }

                prevSNum++; //skip scan num for the ms1

                for(double i=arr[0];i<curPrecur;i+=isoWin)
                {
                    prevSNum++;

                    posMap.put(prevSNum, -1);
                    precursorMap.put(prevSNum, i); //prevPrecur+isoWin*(i+1));
                }

                posMap.put(curSNum, curPos);
                precursorMap.put(curSNum, curPrecur);
            }
            else
            {
                mDiff = curPrecur - prevPrecur - isoWin;


                for(int i=0;i<mDiff;i+=isoWin)
                {
                    prevSNum++;

                    double temp = prevPrecur+isoWin+i;
                    //System.out.println("diff==>>" + mDiff + " " + i + " " + temp + " " + lastPrecur + " "  + isoWin + " " + isoWin*(i+1));
                    //+" " +  prevPrecur+isoWin*(i+1) + indexedFile.getName());
                    if(temp>lastPrecur)
                    {
                        prevSNum++;
                        continue;
                    }

                    posMap.put(prevSNum, -1);
                    precursorMap.put(prevSNum, temp); //prevPrecur+isoWin*(i+1));
                }

                posMap.put(curSNum, curPos);
                precursorMap.put(curSNum, curPrecur); //prevPrecur+isoWin*(i+1));

                prevPrecur = curPrecur;
                prevSNum = curSNum;
            }

        }

        int[] keys=posMap.keys();
        java.util.Arrays.sort(keys);

        for(int i=0;i<keys.length-1;i++)
        {
            double d = precursorMap.get(keys[i+1]) - precursorMap.get(keys[i]);
            if(d != isoWin && d != massDiff)
                System.out.println("error");

            Object precur = ht.get(String.valueOf(keys[i]));

        }
*/
    }

    public static Iterator getAllFiles(String path)
    {
        File f = new File(path);
        return Arrays.asList( f.list() ).iterator();
    }
}
