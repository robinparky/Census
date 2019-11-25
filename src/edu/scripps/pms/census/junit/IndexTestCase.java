package edu.scripps.pms.census.junit;

import junit.framework.*;

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
import java.util.*;
import java.io.*;

public class IndexTestCase extends TestCase {

    private String path="/home/rpark/rpark_on_data/relex_run/BDCHN15";
    private Hashtable ht = new Hashtable();
    private float startRange;
    private float endRange;
    private double massDiff;
    private int numIsoWindows;
    private double isoWin;
    private double lastPrecur;
    private double[] arr;
    private TDoubleArrayList precursorList = new TDoubleArrayList();

    public IndexTestCase(String name)
    {
        super(name);
    }

    protected void setUp() throws Exception
    {
        super.setUp();
        /**@todo verify the constructors*/
    }


    public void indexValidate() throws IOException, Exception
    {
        System.out.println("Validating Relex indexing...");

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

        arr = precursorList.toNativeArray();
        conf.setPrecursorArr(arr);

        numIsoWindows = conf.getNumOfIsolationWindow();
        isoWin = conf.getIsolationWindow();
        lastPrecur = arr[arr.length-1];

        massDiff = isoWin - conf.getMassRange(); //this is negative value

        for( Iterator fileItr = getAllFiles(path); fileItr.hasNext(); )
        {
            String eachFile = fileItr.next().toString();

            if(eachFile.endsWith("index"))
            {
                System.out.println("reading file " + eachFile + "...");

                indexTest(eachFile);
                System.out.println("");

                //              System.out.println(reader.getHeader());
            }
        }

    /*
        try
        {
            System.out.println(TestConfiguration.getUrl() + " , " + Integer.parseInt(TestConfiguration.getPort()) );
            Socket socket = new Socket( TestConfiguration.getUrl(), Integer.parseInt(TestConfiguration.getPort()) );
            Assert.assertTrue(socket.isConnected());
            socket.close();
        }
        catch (UnknownHostException e)
        {
            Assert.fail("Failed : " + e);
        }
        catch(IOException e)
        {
            Assert.fail("Failed " + e);
        }
        catch(Exception e)
        {
            Assert.fail("Failed " + e);
        }
    */
    }

    public void indexTest(String fileName) throws IOException
    {
        MSIndexBuilder builder = new MSIndexBuilder(path + File.separator + fileName);
        TIntLongHashMap posMap = builder.readIndexFile();
        TIntDoubleHashMap precursorMap = builder.getPrecursorMap();

        if(posMap.size() != precursorMap.size())
            Assert.fail("Error : index size is wrong");

        Hashtable indexHash = getIndexData(fileName);

        int[] keys=posMap.keys();
        Arrays.sort(keys);

        for(int i=0;i<keys.length-1;i++)
        {
            double d = precursorMap.get(keys[i+1]) - precursorMap.get(keys[i]);
            if(d != isoWin && d != massDiff)
                Assert.fail("Error " + precursorMap.get(keys[i+1]) + " " + precursorMap.get(keys[i]));

            Object precur = indexHash.get(String.valueOf(keys[i]));
            if(precur==null)
            {
                if(-1!= posMap.get(keys[i]))
                    System.out.println(keys[i] + " " + posMap.get(keys[i]) + " " + precursorMap.get(keys[i]));
                Assert.assertEquals(-1, posMap.get(keys[i]));
            }
            else
            {
                double tempPre = Double.parseDouble(precur.toString());
                Assert.assertTrue(tempPre>0);
            }

            int nextEle = i+numIsoWindows;

            if(nextEle<keys.length)
                Assert.assertEquals( precursorMap.get(keys[i]), precursorMap.get(keys[nextEle]) );
        }

    }

    public Hashtable getIndexData(String fileName) throws IOException
    {
        BufferedReader br = new BufferedReader(new FileReader(path + File.separator + fileName));

        String[] arr;
        Hashtable ht = new Hashtable();
        String lastLine;

        while(null != (lastLine=br.readLine()))
        {
            arr = lastLine.split("\t");
            ht.put(arr[0], arr[1]);
        }
        
        return ht;
    }

    public Iterator getAllFiles(String path)
    {
        File f = new File(path);
        return Arrays.asList( f.list() ).iterator();
    }

    protected void tearDown() throws Exception
    {
        super.tearDown();
    }
}
