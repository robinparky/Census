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


    private  static int num_db_open =0;
    private static int max_spectra_size =2000;
    public  enum StorageMode
    {
        CACHE, REGULAR, NO_STORAGE
    }


    public static SpectraDB connectCreateDB(String filePath) throws Exception {
        File msFile = new File(filePath);
        File msDir = msFile.getParentFile().getAbsoluteFile();
        File spectraDir = new File( msDir.getParentFile().getAbsolutePath() + "/../../spectra/");
        String sqliteDBPath = filePath+".sqlite";
        File sqliteDB = new File(sqliteDBPath);
        if(sqliteDB.exists())
        {
            SpectraDB db = SpectraDB.connectToDBReadOnly(sqliteDBPath);
            return db;
        }
        else if(spectraDir.exists() && spectraDir.isDirectory()){
            String spectraSqliteDBPath = spectraDir.getAbsolutePath() +File.separatorChar + msFile.getName() + ".sqlite";
            File spectraSqliteDB = new File(spectraSqliteDBPath);
            if(spectraSqliteDB.exists())
            {
                SpectraDB spectraDB = SpectraDB.connectToDBReadOnly(spectraSqliteDB.getAbsolutePath());
                return spectraDB;
            }
            else
            {
                CreateDb.createNewDatabase(spectraDir.getAbsolutePath(),msFile.getName()+".sqlite",msFile.getName());
                SpectraDB spectraDB = SpectraDB.connectToDBReadOnly(spectraSqliteDB.getAbsolutePath());
                return  spectraDB;
            }
        }
        else if(!spectraDir.exists())
        {
            CreateDb.createNewDatabase(msDir.getAbsolutePath(),sqliteDB.getName(),msFile.getName());
            SpectraDB spectraDB = SpectraDB.connectToDBReadOnly(sqliteDBPath);
            return  spectraDB;
        }
        return null;
    }

    public static SpectraDB connectToDBReadOnly(String path) throws SQLException {
        return connectToDBReadOnly(path, true, true, false);
    }

    public static SpectraDB connectToDBReadOnly(String path, boolean useCache, boolean savePrevSpectra, boolean preFetch)
            throws SQLException {
        SpectraDB db = new SpectraDB(path);
        db.useCache = useCache;
        db.savePrevSpectra = savePrevSpectra;
        db.scanSpectraMap = useCache ? new Cache<>( 0.75f, true)
                : new HashMap<>(1_000_000);
        db.connectReadOnly();;

        if(preFetch)
        {
            if(useCache)
            {
                db.generateScanSpectraMapLarge();
            }
            else
            {
                db.generateScanSpectraMap();
            }
        }
        return db;
    }

    private static synchronized void editNumSpectaDBOpen(int diff)
    {
        num_db_open += diff;
    }

    public Set<Integer> getScanSet() throws SQLException {
        if(scanSpectraMap == null)
        {
            generateScanSpectraMap();
        }
        return scanSpectraMap.keySet();
    }

    public void generateScanSpectraMap() throws SQLException {
       // scanSpectraMap = new HashMap<>(1_000_000);
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
        result.close();
    }

    public void generateScanSpectraMapLarge() throws SQLException {
        // scanSpectraMap = new HashMap<>(1_000_000);
        String retrieveScans = "SELECT scan, spectrum FROM spectra;";
        PreparedStatement statement = conn.prepareStatement(retrieveScans);
        ResultSet result = statement.executeQuery();
        while((result.next()))
        {
            int scan = result.getInt(1);
            String spectrum = result.getString(2);
            Spectrum spectra = new Spectrum(scan,-1,-1,-1,-1);
            spectra.setSpectrum(spectrum);
            scanSpectraMap.put(scan,spectra);

        }
        result.close();
    }


    public synchronized String getSpectrumStrFromDB(int scan) throws SQLException {
        String retrieveSql = "SELECT spectrum FROM spectra WHERE scan = ?";
        PreparedStatement statement = conn.prepareStatement(retrieveSql);
        statement.setInt(1,scan);
        ResultSet resultSet = statement.executeQuery();
        String spectrum =null;
        while((resultSet.next()))
        {
            spectrum = resultSet.getString(1);
        }
        resultSet.close();
        return  spectrum;
    }

    public synchronized Spectrum getSpectrumFromDB(int scan) throws SQLException {
        Spectrum spectrum = savePrevSpectra ?  scanSpectraMap.get(scan) : null ;
        if(spectrum == null)
        {
            spectrum = getEntryFromDBLight(scan);
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

    private PreparedStatement getEntryFromDBListStatement = null;

    private synchronized Spectrum getEntryFromDBLight(int scan) throws SQLException {
        if(getEntryFromDBListStatement == null)
        {
            String retrieveScans = "SELECT spectrum FROM spectra WHERE scan= ?;";
            getEntryFromDBListStatement = conn.prepareStatement(retrieveScans);
        }
        getEntryFromDBListStatement.setInt(1,scan);
        ResultSet result = getEntryFromDBListStatement.executeQuery();
        while((result.next()))
        {

            String spectrum = result.getString(1);
            Spectrum spectra = new Spectrum(scan,-1, -1,-1,-1);
            spectra.setSpectrum(spectrum);
            // scanSpectraMap.put(scan,spectra);
            return spectra;
        }
        result.close();
        return null;
    }

    private PreparedStatement getEntriesFromDBListStatement = null;

    private List<Spectrum> getEntriesFromDBLight(int scanStart, int scanEnd) throws SQLException {
        if(getEntriesFromDBListStatement ==null)
        {
            String retrieveScans = "SELECT scan, spectrum FROM spectra WHERE scan  between ? and ?;";
             getEntriesFromDBListStatement = conn.prepareStatement(retrieveScans);
        }
        getEntriesFromDBListStatement.setInt(1,scanStart);
        getEntriesFromDBListStatement.setInt(2,scanEnd);
        List<Spectrum> specList = new ArrayList<>();
        ResultSet result = getEntriesFromDBListStatement.executeQuery();
        while((result.next()))
        {
            int scan = result.getInt(1);
            String spectrum = result.getString(2);
            Spectrum spectra = new Spectrum(scan,-1, -1,-1,-1);
            spectra.setSpectrum(spectrum);
            // scanSpectraMap.put(scan,spectra);
            specList.add(spectra);
        }
        result.close();
        return specList;
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
        result.close();
        return null;
    }


    private SpectraDB(String path) {
        this.path = path;
    }
    private void connectReadOnly() throws SQLException {
        String url = "jdbc:sqlite:" +path;
        SQLiteConfig config = SpectraDB.GetDefaultConfig();
        config.setReadOnly(true);
        config.setTempStore(SQLiteConfig.TempStore.FILE);
        //config.setPageSize(100_000_000);
        config.setCacheSize(-50_000);
        config.setJournalMode(SQLiteConfig.JournalMode.OFF);
        System.out.println("<<> " +url);
        conn = DriverManager.getConnection(url,config.toProperties());
        editNumSpectaDBOpen(1);
    }




    @Override
    public void close() throws IOException {

        if(conn!=null)
        {
            try {
                if(!conn.isClosed())
                {
                    editNumSpectaDBOpen(-1);
                    conn.close();
                }
            } catch (SQLException e) {
                e.printStackTrace();
            }

        }
    }


    public static SQLiteConfig GetDefaultConfig()
    {
        final int cacheSize = 100_000 ;
        final int pageSize = 4096;
        SQLiteConfig config = new SQLiteConfig();
        //optimize for multiple connections that can share data structures
        config.setSharedCache(true);
        config.setCacheSize(cacheSize);
        config.setPageSize(pageSize);
        config.setJournalMode(SQLiteConfig.JournalMode.OFF);
        config.enableFullSync(false);
        config.enableRecursiveTriggers(false);
        config.setLockingMode(SQLiteConfig.LockingMode.EXCLUSIVE);
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
            setMaxSpectraSize(mzList.size());
        }

        public String getSpectrum() throws SQLException {
            return spectrum;
        }

        private void setSpectrum(String spectrum) {
            this.spectrum = spectrum;
            fillArrays();
        }


        public TDoubleArrayList getMzList() {
            return mzList;
        }

        public TDoubleArrayList getIntensityList()  {
            return intensityList;
        }


    }



    public static class Cache<K , V> extends LinkedHashMap<K , V>
    {
        public static final int MAX_SIZE = 500_000;
        public Cache(float loadFactor, boolean accessOrder) {
            super(50_000, loadFactor, accessOrder);
            //MAX_SIZE = initialCapacity;
        }

        @Override
        protected boolean removeEldestEntry(Map.Entry eldest) {
            return size() > getMaxSize();
        }

        private static int getMaxSize()
        {
            int ratio = 2000/max_spectra_size;
            int tempMax = MAX_SIZE * ratio /getNumDBOpen();
            if(tempMax < 100 )
                return 100;
            else
                return MAX_SIZE * ratio /getNumDBOpen();
        }

    }




    private synchronized static int getNumDBOpen()
    {
        if(num_db_open>0)
        {
            return num_db_open;
        }
        else
        {
            return 1;
        }
    }

    private synchronized static void setMaxSpectraSize(int size)
    {
        if(size>max_spectra_size)
        {
            max_spectra_size = size;
        }
    }




}
