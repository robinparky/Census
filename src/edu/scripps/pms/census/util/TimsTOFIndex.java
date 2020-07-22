package edu.scripps.pms.census.util;

import java.io.*;
import java.util.HashMap;
import java.util.Map;

public class TimsTOFIndex {

    public final String ms2Path;
    public final String indexpath;
    public final static String SUFFIX = ".timsstofindex";
    private Map<Integer,Integer> scan2ParentMap = new HashMap<>();


    public TimsTOFIndex(String path) throws IOException {
        if(path.endsWith(".ms2"))
        {
            this.ms2Path = path;
            this.indexpath = ms2Path + SUFFIX;
        }
        else
        {
            this.ms2Path = path.replace(SUFFIX, "");
            this.indexpath =path;
        }
        init();
    }

    private void init() throws IOException {
        File indexFile = new File(indexpath);
        if(indexFile.exists())
        {
            System.out.println(">>>> IndexFile exists; reading exisitng index file");
            scan2ParentMap = readIndex(indexpath);
        }
        else
        {
            System.out.println(">>>> IndexFile does not exist; reading ms2 file and writing index file");
            scan2ParentMap = readMs2(ms2Path);
            writeIndex(indexpath, scan2ParentMap);
        }
    }

    private static Map<Integer, Integer> readMs2(String ms2Path) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(ms2Path));
        String line;
        Map<Integer,Integer> scan2PrecusorMap = new HashMap<>();
        int scanNumber =-1;
        while((line = br.readLine())!=null)
        {
            if(line.startsWith("S\t"))
            {
                String [] arr = line.split("\t");
                scanNumber = Integer.parseInt(arr[1]);
            }
            else if(line.startsWith("I\tTIMSTOF_Precursor_ID"))
            {
                String [] arr = line.split("\t");
                int prec = Integer.parseInt(arr[2]);
                scan2PrecusorMap.put(scanNumber,prec);
            }
        }
        br.close();
        return scan2PrecusorMap;
    }

    private static Map<Integer,Integer> readIndex(String indexPath) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(indexPath));
        String line;
        Map<Integer,Integer> scan2PrecusorMap = new HashMap<>();
        while((line = br.readLine())!=null)
        {
            String [] arr = line.split("\t");
            int scanNumber = Integer.parseInt(arr[0]);
            int precNumber = Integer.parseInt(arr[1]);
            scan2PrecusorMap.put(scanNumber, precNumber);
        }
        br.close();;
        return scan2PrecusorMap;
    }

    private static void writeIndex(String indexpath, Map<Integer,Integer> scan2PrecursorMap) throws IOException {
        BufferedWriter bw  = new BufferedWriter(new FileWriter(indexpath));
        for(Map.Entry<Integer,Integer> entry: scan2PrecursorMap.entrySet())
        {
            String s1 = Integer.toString(entry.getKey());
            String s2 = Integer.toString(entry.getValue());
            bw.append(s1).append("\t").append(s2);
            bw.newLine();
        }
        bw.close();
    }

    public Map<Integer, Integer> getScan2ParentMap() {
        return scan2ParentMap;
    }

    public int getPrecursorID(int ms2ScanNumber)
    {
        return scan2ParentMap.get(ms2ScanNumber);
    }


    public static void main(String [] args) throws IOException {
        String path = args[0];
        TimsTOFIndex index = new TimsTOFIndex(path);
    }


}
