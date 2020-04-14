package edu.scripps.pms.census.hash;

import edu.scripps.pms.util.FileFilterUtil;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

/**
 * Created by rpark on 10/8/16.
 */
public class IndexUtil {

    public static void main(String[] args) throws Exception {

     //   buildMS2toMS1ScanMapFiles("/data/2/rpark/ip2_data//rpark/an_chi/Onbeads_Lung_DMSO_CS_fixed_2016_09_12_16_84798/spectra/");


        HashMap<String, HashMap<Integer, Integer>> ms2ToMs1Map = IndexUtil.buildMS2toMS1ScanMapFiles("/data/2/rpark/ip2_data//rpark/Indiana_University_Lisa_Jones_BlindPTM/FLOW_CELL_NL1_SubSearch_2016_08_07_14_83921/search/projects2016_08_07_14_101413/../../spectra/");


        if(true) return;
        IndexUtil.buildMS2toMS1ScanMap(
                "/home/rpark/a.ms1.index",
                "/home/rpark/a.ms2.index"
        );

            //"/data/2/rpark/ip2_data/rpark/an_chi/Onbeads_Lung_100nM_CS_fixed_2016_09_12_16_84796/spectra/Onbeads_Lung_100nM.ms1.index",
            //"/data/2/rpark/ip2_data/rpark/an_chi/Onbeads_Lung_100nM_CS_fixed_2016_09_12_16_84796/spectra/Onbeads_Lung_100nM.ms2.index");
    }

    public static HashMap<String, HashMap<Integer, Integer>> buildMS2toMS1ScanMapFiles(String path) throws IOException {
        HashMap<String, HashMap<Integer, Integer>> map = new HashMap<>();
        if(!path.endsWith("/")) path += "/";

        ArrayList<String> list = FileFilterUtil.getFilesBySuffix(path, "ms1");
        for(Iterator<String> itr=list.iterator(); itr.hasNext(); ) {
            String ms1File = itr.next();
            String ms2File = ms1File.substring(0, ms1File.length() - 1) + "2";

            if(new File(path + ms2File).exists()) {
                File ms1Index = new File(path + ms1File +".index");
                File ms2Index = new File(path + ms2File+".index");
                if(!ms1Index.exists())
                    MSIndexFileCreator.createIndexFile(path + ms1File);
                if(!ms2Index.exists())
                    MSIndexFileCreator.createIndexFile(path + ms2File);

                HashMap<Integer, Integer> eachMap = buildMS2toMS1ScanMap(path + ms1File +".index", path + ms2File+".index");
                map.put(path + ms1File, eachMap);
            }
        }

        return map;
    }

    public static HashMap<Integer, Integer> buildMS2toMS1ScanMapWithNewFormat(String ms2IndexFile) throws IOException {

        HashMap<Integer, Integer> map = new HashMap<>();

        BufferedReader br = new BufferedReader(new FileReader(ms2IndexFile));
        String eachLine;

        while(null != (eachLine=br.readLine())) {
            String[] arr = eachLine.split("\t");

            map.put(Integer.parseInt(arr[0]), Integer.parseInt(arr[5]));

        }

        br.close();

        return map;

    }

    public static HashMap<Integer, Integer> buildMS2toMS1ScanMap(String ms1IndexFile, String ms2IndexFile) throws IOException {

        HashMap<Integer, Integer> map = new HashMap<>();

        BufferedReader br = new BufferedReader(new FileReader(ms2IndexFile));
/*
        System.out.println("------------------->>" + ms2IndexFile);
        System.out.println(br);
        System.out.println(br.readLine());
        System.out.println(br.readLine().split("\t"));
*/

        String[] arr= br.readLine().split("\t");
        boolean isNewIndexFormat=false;
        if(arr.length>=6) {
               map = buildMS2toMS1ScanMapWithNewFormat(ms2IndexFile);
        }

        br.close();

        if(isNewIndexFormat)
            return map;


        List<Integer> ms1ScanList = getScanList(ms1IndexFile);
        List<Integer> ms2ScanList = getScanList(ms2IndexFile);

        int count=0;
        int curMs1Scan=ms1ScanList.get(count);
        for(Iterator<Integer> itr=ms2ScanList.iterator(); itr.hasNext(); ) {
            int ms2Scan = itr.next();

            while(curMs1Scan<ms2Scan) {
                count++;
                if(ms1ScanList.size()<=count) break;

               // if(count>=ms1ScanList.size())
                //    System.out.println("out");
                curMs1Scan = ms1ScanList.get(count);
            }

            if(ms1ScanList.size()<=count)
                map.put(ms2Scan, ms1ScanList.get(ms1ScanList.size()-1));
            else
                map.put(ms2Scan, ms1ScanList.get(count-1));
                //System.out.println(ms2Scan + " " + ms1ScanList.get(count-1) );

        }
        return map;

    }

    public static List<Integer> getScanList(String file) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(file));
        String lastLine;

        List<Integer> scanList = new ArrayList<>();
        while( null != (lastLine=br.readLine()) ) {
            String[] arr = lastLine.split("\t");
            scanList.add(Integer.parseInt(arr[0]));
           // System.out.println(arr[0]);
        }

        br.close();

        return scanList;

    }

}
