package edu.scripps.pms.util.sqlite.spectra;

import edu.scripps.pms.util.FileFilterUtil;
import org.sqlite.SQLiteConfig;

import java.io.*;
import java.sql.*;
import java.util.List;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadPoolExecutor;


public class CreateDb {

    public static void main(String[] args) throws Exception {

        /*
        String path = "/home/bsuwirjo/programming/scrippsResearch/testData/";
        String name = "a.ms2.sqlite";
        String result = getSpectrum(path, name, 13);
        System.out.println(result);
        */

        //String path = "/data/2/rpark/ip2_data/rpark/Demo_Project_C_elegant_label_free_test/D2a_2010_07_14_15_1521/spectra/temp";

        String path=args[0];


        if(!path.endsWith(File.separator))
            path += File.separator;

        List<String> fileList = FileFilterUtil.getFilesBySuffix(path, ".ms2");

        for(String each:fileList) {



            createNewDatabase(path, each + ".sqlite", each);
        }



        fileList = FileFilterUtil.getFilesBySuffix(path, ".ms1");

        for(String each:fileList) {



            createNewDatabase(path, each + ".sqlite", each);
        }


        fileList = FileFilterUtil.getFilesBySuffix(path, ".ms3");

        for(String each:fileList) {



            createNewDatabase(path, each + ".sqlite", each);
        }


      //  System.out.println("aaa");

    //      createNewDatabase(path, "test.db", "test");*/

    }

    public static void createNewDatabase(String path, String dbFilename, String spectralFile) throws IOException {
        createNewDatabase(path,dbFilename,spectralFile,false);
    }
        /**
         * Connect to a sample database
         *
         * @param dbFilename the database file name
         */
        public static void createNewDatabase(String path, String dbFilename, String spectralFile, boolean deleteMode) throws IOException {
            //String url = "jdbc:sqlite:/home/rpark/temp/" + fileName;
            if(!path.endsWith(File.separator))
                path+=File.separator;
            String filePath = path +dbFilename;
            String journal = path + dbFilename+".sqlite-journal";
            File f = new File(filePath);
            File journalF = new File(journal);
            if(!f.exists())
            {
                f.createNewFile();
            }
            else {
                if (deleteMode || f.length()==0)
                {
                    f.delete();
                    f.createNewFile();
                }
                else
                {
                    return;
                }
            }
            if(journalF.exists())
            {
                journalF.delete();
            }
            System.out.println("edu.scripps.pms.util.sqlite.spectra.CreateDb: creating " + dbFilename);

            SQLiteConfig config = SpectraDB.GetDefaultConfig();


            String url = "jdbc:sqlite:" + path + dbFilename;

            // SQL statement for creating a new table
            //String sql = "CREATE TABLE IF NOT EXISTS spectra (\n"
            String sql = "CREATE TABLE IF NOT EXISTS spectra (\n"
                    + "id integer PRIMARY KEY,\n"
                    + "scan integer NOT NULL,\n"
                    + "    filename text,\n"
                    + "    scanFilename text,\n"
                    + "spectrum text NON NULL,\n"
                    + "prcMass REAL ,\n"
                    + "prcMass_z REAL ,\n"
                    + "retTime REAL ,\n"
                    + "charge integer "
                    + ");";
            String sqlIndex = "CREATE INDEX scan_index ON spectra(scan);\n";
            String sqldrop = "DROP TABLE IF EXISTS spectra\n";
         //   ThreadPoolExecutor executor =  (ThreadPoolExecutor) Executors.newFixedThreadPool(1);

            PreparedStatement pstmt = null;
            try (Connection conn = DriverManager.getConnection(url,config.toProperties());
                 Statement stmt = conn.createStatement()) {
                // create a new table
                stmt.execute(sqldrop);
                stmt.execute(sql);
           //     conn.setAutoCommit(false);



                BufferedReader br = new BufferedReader(new FileReader(path + spectralFile));
                String eachLine="";
                int scanNum=0;
                int previousScan=0;
                StringBuffer spectrum = new StringBuffer();

                while( null != (eachLine = br.readLine()) && !eachLine.startsWith("S\t") );
                spectrum.append(eachLine).append("\n");

                String[] arr = eachLine.split("\t");
                scanNum = Integer.parseInt(arr[1]);
                double prcMass = -1;
                if(arr.length>3)
                    prcMass = Double.parseDouble(arr[3]);
                //insert sql

                String insertSql = "INSERT INTO spectra(scan, spectrum, prcMass,prcMass_z ,retTime, charge) VALUES(?,?,?,?,?,?);";
                pstmt = conn.prepareStatement(insertSql);
                double retTime = -1;
                int cs =-1;
                int batchCount =0;
                double prcMass_z = -1;
                while( null != (eachLine = br.readLine()) ) {

                    if(eachLine.startsWith("S\t")) {

                    //    System.out.println("========================");
                   //     System.out.println(scanNum);
                   //     System.out.println(spectrum.toString());


                        pstmt.setInt(1, scanNum);

                        pstmt.setString(2, spectrum.toString());
                        pstmt.setDouble(3, prcMass);
                        pstmt.setDouble(4, prcMass_z);
                        pstmt.setDouble(5, retTime);
                        pstmt.setInt(6, cs);



                        pstmt.addBatch();
                        batchCount++;
                        if(batchCount>=1_000)
                        {
                            batchCount =0;/*
                            PreparedStatement finalPstmt = pstmt;
                            executor.submit(() ->{
                                finalPstmt.executeBatch();
                                return null;
                            });*/
                            pstmt.executeBatch();

                            pstmt = conn.prepareStatement(insertSql);
                        }


                   //     System.out.println("========================>>");

                        spectrum.delete(0, spectrum.length());
                        arr = eachLine.split("\t");
                        scanNum = Integer.parseInt(arr[1]);
                        if(arr.length>3)
                            prcMass = Double.parseDouble(arr[3]);
                    }
                    if(eachLine.startsWith("I\tRetTime"))
                    {
                        arr = eachLine.split("\t");
                        String toConvert = arr[2];
                        retTime = Double.parseDouble(toConvert);
                    }
                    if(eachLine.startsWith("Z\t"))
                    {
                        arr = eachLine.split("\t");
                        String toConvert = arr[1];
                        cs = Integer.parseInt(toConvert);
                        toConvert = arr[2];
                        prcMass_z = Double.parseDouble(toConvert);
                    }


                    spectrum.append(eachLine).append("\n");





                    //System.out.println(eachLine);
                }
              //  executor.shutdown();
                pstmt.setInt(1, scanNum);
                pstmt.setString(2, spectrum.toString());
                pstmt.setDouble(3, prcMass);
                pstmt.setDouble(4, prcMass_z);
                pstmt.setDouble(5, retTime);
                pstmt.setInt(6, cs);
                pstmt.addBatch();
                pstmt.executeBatch();

                stmt.execute(sqlIndex);
                //conn.commit();

                stmt.execute("VACUUM;\n");
                //Statement vacuumStatement = conn.createStatement();
                //vacuumStatement.execute("VACUUM");
               // conn.commit();

            } catch (SQLException e) {
                System.out.println(e.getMessage());
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }

  //          pstmt.close();
        }

        //Connect to database
        private static Connection connect(String url) {
            // SQLite connection string

            Connection conn = null;
            try {
                conn = DriverManager.getConnection(url);
            } catch (SQLException e) {
                System.out.println(e.getMessage());
            }
            return conn;
        }


        public static String getSpectrum(String path, String dbFilename, int scanNum){

            //Create url and retrieve cmd
            String url = "jdbc:sqlite:" + path + dbFilename;
            String retrieveSql = "SELECT * FROM spectra WHERE scan = ?";

            //System.out.println(url);

            try (Connection conn = connect(url);
                 PreparedStatement stmt  = conn.prepareStatement(retrieveSql);){
                //Set value to inputted scan number
                stmt.setInt(1, scanNum);
                ResultSet rs  = stmt.executeQuery();

                String spectrum = rs.getString("spectrum");
                System.out.println(rs.getInt("scan"));
                return spectrum;

            } catch (SQLException e) {
                System.out.println(e.getMessage());
                return e.getMessage();
            }
        }
}
