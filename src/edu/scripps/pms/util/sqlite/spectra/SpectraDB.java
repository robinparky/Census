package edu.scripps.pms.util.sqlite.spectra;

import gnu.trove.TDoubleArrayList;
import org.sqlite.SQLiteConfig;

import java.io.*;
import java.sql.*;
import java.util.*;

public class SpectraDB implements Closeable {

    public final String path;
    private Connection conn = null;
  //  private Map<Integer, Spectrum> scanSpectraMap = new HashMap<>(1_000_000);
    private Map<Integer, Spectrum> scanSpectraMap = null;
    private boolean savePrevSpectra =true;
    private boolean useCache = false;

    public static SpectraDB connectToDBReadOnly(String path) throws SQLException {
        return connectToDBReadOnly(path, false, true, true);
    }

    public static SpectraDB connectToDBReadOnly(String path, boolean useCache, boolean savePrevSpectra, boolean preFetch)
            throws SQLException {
        SpectraDB db = new SpectraDB(path);
        db.useCache = useCache;
        db.savePrevSpectra = savePrevSpectra;
        db.scanSpectraMap = useCache ? new Cache<>(100_000_000, 0.75f, true)
                : new HashMap<>(1_000_000);
        db.connectReadOnly();;
        if(preFetch)
        {
            db.generateScanSpectraMap();
        }
        return db;
    }


    public Set<Integer> getScanSet() throws SQLException {
        if(scanSpectraMap == null)
        {
            generateScanSpectraMap();
        }
        return scanSpectraMap.keySet();
    }

    public void generateScanSpectraMap() throws SQLException {
        scanSpectraMap = new HashMap<>(1_000_000);
        String retrieveScans = "SELECT scan, prcMass, retTime, charge, prcMass_z FROM spectra;";
        PreparedStatement statement = conn.prepareStatement(retrieveScans);
        ResultSet result = statement.executeQuery();
        while((result.next()))
        {
            int scan = result.getInt(1);
            double prcMass =result.getDouble(2);
            double retTime = result.getDouble(3);
            int charge = result.getInt(4);
            double prcMassZ = result.getDouble(5);
            Spectrum spectra = new Spectrum(scan,prcMass, prcMassZ,retTime,charge);
            scanSpectraMap.put(scan,spectra);
        }
    }

    public String getSpectrumStrFromDB(int scan) throws SQLException {
        String retrieveSql = "SELECT spectrum FROM spectra WHERE scan = ?";
        PreparedStatement statement = conn.prepareStatement(retrieveSql);
        statement.setInt(1,scan);
        ResultSet resultSet = statement.executeQuery();
        String spectrum =null;
        while((resultSet.next()))
        {
            spectrum = resultSet.getString(1);
        }
        return  spectrum;
    }

    public Spectrum getSpectrumFromDB(int scan) throws SQLException {
        Spectrum spectrum = savePrevSpectra ?  scanSpectraMap.get(scan) : null ;
        if(spectrum == null)
        {
            spectrum = getEntryFromDB(scan);
            if(savePrevSpectra)
                scanSpectraMap.put(scan,spectrum);
        }
        if(spectrum == null )
            return null;
        if(spectrum.getSpectrum()==null)
        {
            spectrum.setSpectrum(getSpectrumStrFromDB(scan));
        }
        return  spectrum;
    }

    private Spectrum getEntryFromDB(int scan) throws SQLException {
        String retrieveScans = "SELECT scan, prcMass, retTime, charge, prcMass_z, spectrum FROM spectra WHERE scan= ?;";
        PreparedStatement statement = conn.prepareStatement(retrieveScans);
        statement.setInt(1,scan);
        ResultSet result = statement.executeQuery();
        while((result.next()))
        {
            double prcMass =result.getDouble(2);
            double retTime = result.getDouble(3);
            int charge = result.getInt(4);
            double prcMassZ = result.getDouble(5);
            String spectrum = result.getString(6);
            Spectrum spectra = new Spectrum(scan,prcMass, prcMassZ,retTime,charge);
            spectra.setSpectrum(spectrum);
           // scanSpectraMap.put(scan,spectra);
            return spectra;
        }
        return null;
    }


    private SpectraDB(String path) {
        this.path = path;
    }
    private void connectReadOnly() throws SQLException {
        String url = "jdbc:sqlite:" +path;
        SQLiteConfig config = SpectraDB.GetDefaultConfig();
        config.setReadOnly(true);
        config.setJournalMode(SQLiteConfig.JournalMode.OFF);

        conn = DriverManager.getConnection(url,config.toProperties());
    }




    @Override
    public void close() throws IOException {
        if(conn!=null)
        {
            try {
                if(!conn.isClosed())
                    conn.close();
            } catch (SQLException e) {
                e.printStackTrace();
            }

        }
    }


    public static SQLiteConfig GetDefaultConfig()
    {
        final int cacheSize = 100000 / 6;
        final int pageSize = 4096;
        SQLiteConfig config = new SQLiteConfig();
        //optimize for multiple connections that can share data structures
        config.setSharedCache(true);
        config.setCacheSize(cacheSize);
        config.setPageSize(pageSize);
        config.setJournalMode(SQLiteConfig.JournalMode.TRUNCATE);
        config.enableFullSync(false);
        config.enableRecursiveTriggers(false);
        config.setLockingMode(SQLiteConfig.LockingMode.NORMAL);
        config.setSynchronous(SQLiteConfig.SynchronousMode.OFF); //TODO may be dangerous on some systems to have off
        return config;
    }

    public static class Spectrum
    {
        private static int ID_COUNTER =0;
        public final int id;
        public final int scanNumber;
        public final double prcMass;
        public final double prcMassZ;
        public final double retTime;
        public final int charge;
        private String spectrum =null;
        private TDoubleArrayList mzList = null;
        private TDoubleArrayList intensityList = null;


        public Spectrum(int scanNumber, double prcMass, double prcMassZ, double retTime, int charge) {
            this.scanNumber = scanNumber;
            this.prcMassZ = prcMassZ;
            this.prcMass = prcMass;
            this.retTime = retTime;
            this.charge = charge;
            this.id = ID_COUNTER++;
        }

        private void fillArrays() {
            String [] lineArr = spectrum.split("\n");
            mzList = new TDoubleArrayList();
            intensityList = new TDoubleArrayList();
            for(String line: lineArr)
            {
                if(Character.isDigit(line.charAt(0)))
                {
                    String [] arr = line.split(" ");
                    String mzStr = arr[0];
                    String intStr = arr[1];
                    double mz = Double.parseDouble(mzStr);
                    double intensity = Double.parseDouble(intStr);
                    mzList.add(mz);
                    intensityList.add(intensity);
                }
            }
        }

        public String getSpectrum() throws SQLException {
            return spectrum;
        }

        private void setSpectrum(String spectrum) {
            this.spectrum = spectrum;
            fillArrays();
        }


        public TDoubleArrayList getMzList() throws SQLException {
            return mzList;
        }

        public TDoubleArrayList getIntensityList() throws SQLException {
            return intensityList;
        }


    }



    public static class Cache<K , V> extends LinkedHashMap<K , V>
    {
        public final int MAX_SIZE;
        public Cache(int initialCapacity, float loadFactor, boolean accessOrder) {
            super(initialCapacity, loadFactor, accessOrder);
            MAX_SIZE = initialCapacity;
        }

        @Override
        protected boolean removeEldestEntry(Map.Entry eldest) {
            return size() > MAX_SIZE;
        }
    }
}
