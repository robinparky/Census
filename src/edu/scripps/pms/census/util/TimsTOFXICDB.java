package edu.scripps.pms.census.util;

import edu.scripps.pms.util.sqlite.spectra.SpectraDB;
import gnu.trove.TDoubleArrayList;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.lang3.tuple.Triple;
import org.sqlite.SQLiteConfig;
import rpark.statistics.model.GaussianPeakModel;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.sql.*;
import java.util.*;

import static rpark.statistics.GaussianFitting.getGaussianPeakRangeIndex;
import static rpark.statistics.model.GaussianPeakModel.getGaussianPeakArea;

public class TimsTOFXICDB implements Closeable {

    public final String dbPath;
    private Map<Integer,Integer> ms2TimstofPrecursorMap = null;
    private Connection conn = null;
    private PreparedStatement queryPrecursor  =null;
    private TimsTOFIndex index = null;
    public static final String DB_NAME = "precursor_xics.sqlite";

    public TimsTOFXICDB(String path) throws SQLException {
        dbPath = path;
        String url = "jdbc:sqlite:" +path;
        SQLiteConfig config = SpectraDB.GetDefaultConfig();
        config.setReadOnly(true);
        config.setTempStore(SQLiteConfig.TempStore.FILE);
        //config.setPageSize(100_000_000);
        config.setCacheSize(-50_000);
        config.setJournalMode(SQLiteConfig.JournalMode.OFF);
        conn = DriverManager.getConnection(url,config.toProperties());
    }

    public static class TimstofQueryResult
    {
        private List<Triple<String, Double,Double>> initialResult;
        public final double retTime;
        private List<Pair<Double,Double>> summedList = new ArrayList<>();

        public TimstofQueryResult(List<Triple<String, Double, Double>> initialResult, double retTime) {
            this.initialResult = initialResult;
            this.retTime = retTime;
        }
        private void sumPeaks()
        {
            Map<String, Pair<Double,Double>> map = new HashMap<>();
            for(Triple<String, Double, Double> r: initialResult)
            {
                Pair<Double,Double> p = map.get(r.getLeft());
                if(p!=null)
                {
                    map.put(r.getLeft(), Pair.of(r.getMiddle(), p.getRight() + r.getRight()));
                }
                else
                {
                    map.put(r.getLeft(), Pair.of(r.getMiddle(), r.getRight()));
                }
            }
            summedList = new ArrayList<>();
            for(Map.Entry<String, Pair<Double,Double>> entry: map.entrySet() )
            {
                summedList.add(entry.getValue());
            }
            summedList.sort(Comparator.comparingDouble(Pair::getLeft));
        }

        public List<Pair<Double, Double>> getSummedList() {
            return summedList;
        }
    }


    public TimstofQueryResult queryParentId(int id) throws SQLException {
        if(queryPrecursor == null)
        {
            queryPrecursor = conn.prepareStatement("select val.Time, val.Area, pre.Time as PrecTime, pre.Intensity" +
                    " as PreIntensity from XICValues val JOIN XICIsotopeTraces tr ON val.XicTrace = tr.Id" +
                    "  JOIN XICMapping mp on mp.XicId=tr.XIC" +
                    "  JOIN Precursors pre on pre.Id=mp.Precursor" +
                    "  Where pre.Id = ?");
        }
        queryPrecursor.setInt(1,id);
        ResultSet rs = queryPrecursor.executeQuery();
        List<Triple<String, Double,Double>> resultList = new ArrayList<>();

        double prectime =0;
        while(rs.next())
        {
            String s = rs.getString(1);

            double time = rs.getDouble(1);
            double area = rs.getDouble(2);
             prectime = rs.getDouble(3);
            resultList.add(Triple.of(s, time,area));
        }
        return new TimstofQueryResult(resultList, prectime);
    }




    public TimstofQueryResult queryAndSumParentId(int id) throws SQLException {
        TimstofQueryResult queryResult = queryParentId(id);
        queryResult.sumPeaks();
        return queryResult;
    }



    public TimstofQueryResult queryAndSumMS2(int ms2id) throws SQLException {
        int id = index.getParentId(ms2id);
        return queryAndSumParentId(id);
    }

    @Override
    public void close() throws IOException {
        try {
            if(conn != null && !conn.isClosed())
            {
                conn.close();
            }
        } catch (SQLException throwables) {
            throwables.printStackTrace();
        }
    }

    public static double getPeakAreaEstimate(double [] xarr, double[] yarr)
    {
        double areaSum = 0;

        double firstHeight =0;
        double firstWidth =0 ;

        double lastHeight =0;
        double lastWidth =0 ;
        for(int i=1; i<xarr.length; i++)
        {

            double width = xarr[i] - xarr[i-1];
            double height = (yarr[i] + yarr[i-1])/2;
            areaSum += width*height;
            if(i==1)
            {

                firstWidth = width;
            }
            else if(i==2)
            {
                firstWidth = firstWidth -( width - firstWidth);
            }
            else if(i==xarr.length-2)
            {
                lastWidth = width;
            }
            else if(i==xarr.length-1)
            {
                lastWidth = width + (lastWidth-width);
            }
        }

        double startWidth = firstWidth;
        double startHeight =yarr[0]/2;
        areaSum += startHeight*startWidth;

        double endWidth = lastWidth;
        double endHeight =yarr[yarr.length-1]/2;
        areaSum += endHeight*endWidth;
        return areaSum;
    }


    public static void main(String [] args) throws SQLException, IOException {
        String path = args[0];
        String ms2path = args[1];
        int id = Integer.parseInt(args[2]);
        TimsTOFIndex index = new TimsTOFIndex(ms2path);
        int precursorId = index.getParentId(id);
        TimsTOFXICDB timsTOFXICDB = new TimsTOFXICDB(path);
        List<Pair< Double,Double>> resutlt = timsTOFXICDB.queryAndSumParentId(precursorId).summedList;
        TDoubleArrayList xarrayList = new TDoubleArrayList();
        TDoubleArrayList yarrayList = new TDoubleArrayList();
        double max = Double.MIN_VALUE;
        double sum = 0;
        for(Pair< Double,Double> r: resutlt)
        {
            xarrayList.add(r.getLeft());
            yarrayList.add(r.getRight());
            if(r.getRight()> max)
                max = r.getRight();
            sum+=r.getRight();
            //System.out.println(r.getLeft()+"\t"+r.getRight());
        }
        double [] xarr = xarrayList.toNativeArray();
        double [] yarr = yarrayList.toNativeArray();
        double predictedArea = getPeakAreaEstimate(xarr, yarr);
        GaussianPeakModel model =  getGaussianPeakRangeIndex(xarr, yarr, -1, resutlt.size()-1);
        model.setMaxIntensity(max);
        double peakHeight = model.getY();
        double sigma = model.getSigma();
        double area = getGaussianPeakArea(peakHeight, sigma);
        System.out.println("<<>> rectangular peak area:\t"+predictedArea);
        System.out.println("<<>> gaussian peak area:\t"+area);

        timsTOFXICDB.close();

    }

    public TimsTOFIndex getIndex() {
        return index;
    }

    public void setIndex(TimsTOFIndex index) {
        this.index = index;
    }

    public static Map<String, TimsTOFXICDB> createDBIndexMap(String path) throws IOException, SQLException {

        if (!path.endsWith(File.separator)) {
            path += File.separator;
        }

        File f = new File(path);
        String[] list = f.list(new RelExFileFilter(".ms2"));
        Map<String,TimsTOFXICDB> timsTOFXICDBMap = new HashMap<>();

        for(String ms2Path :list)
        {
            String rawName = ms2Path.replace("_nopd.ms2",".d");
            TimsTOFIndex index = new TimsTOFIndex(path+ms2Path);
            TimsTOFXICDB db = new TimsTOFXICDB(path+rawName+File.separator+DB_NAME);
            db.setIndex(index);
            timsTOFXICDBMap.put(ms2Path, db);

        }
        return timsTOFXICDBMap;


    }



}
