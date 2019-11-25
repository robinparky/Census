package edu.scripps.pms.census;

import edu.scripps.pms.census.conf.Configuration;
import edu.scripps.pms.census.exception.InvalidAAException;
import edu.scripps.pms.census.hash.IndexUtil;
import edu.scripps.pms.census.hash.IndexedFile;
import edu.scripps.pms.census.io.IsotopeReader;
import edu.scripps.pms.census.model.SampleModel;
import edu.scripps.pms.census.util.dtaselect.Peptide;
import edu.scripps.pms.census.util.dtaselect.Protein;
import edu.scripps.pms.census.util.io.DTASelectFilterReader;
import edu.scripps.pms.util.seq.Fasta;
import gnu.trove.TDoubleIntHashMap;
import rpark.statistics.BinarySearch;
import rpark.statistics.model.GaussianPeakModel;
import scripts.MSSplitFolderCreation;

import java.io.*;
import java.util.*;

import static edu.scripps.pms.census.labelFree.LabelfreeMissingPeptideBuilderSplit.isotopeCalc;
import static edu.scripps.pms.census.labelFree.LabelfreeMissingPeptideBuilderSplit.redoNormalization;

/**
 * Created by Titus Jung titusj@scripps.edu on 5/28/19.
 */
public class AccurateMassTag {

    private Map<String,IndexedFile> indexedMs1FileMap = new HashMap<>();
    private Map<String, String> splitSpectraMap = new HashMap<>();
    private Map<String, IndexedFile> splitMs1FileHt = new HashMap<>();
    private Map<String, HashMap<Integer, Integer>> ms2ToMs1Map = new HashMap<>();
    private List<Map<String, List<PeptideContainer>>> locusPeptideMap = new ArrayList<>();
    private Configuration conf = null;
    private         IsotopeReader isoReader = null;
    private String path;

    public AccurateMassTag(String configurationPath ) throws Exception {
        conf = Configuration.getInstance();
        path = configurationPath.substring(0,configurationPath.lastIndexOf(File.separator));
        if (!conf.isReadConfigFile()) {
            conf.setLabelfree(true);
            conf.readXMLParam(configurationPath);
        }
        init();
    }


    public static class PeptideContainer {
        public final String key;
        public final String sequence;
        public final String cs;
        public final int csInt;
        public final int scanNumber;
        public final String ms2Name;
        public final String dtaPath;
        private double retentionTime;
        private List<Protein> proteinList = new ArrayList<>();

        public PeptideContainer(String key, String sequence, String cs, int scanNumber, String ms2Name, String dtaPath) {
            this.key = key;
            this.sequence = sequence;
            this.cs = cs;
            csInt = Integer.parseInt(cs);
            this.scanNumber = scanNumber;
            this.ms2Name = ms2Name;
            this.dtaPath = dtaPath;
        }

        public double getRetentionTime() {
            return retentionTime;
        }

        public void setRetentionTime(double retentionTime) {
            this.retentionTime = retentionTime;
        }

        public List<Protein> getProteinList() {
            return proteinList;
        }

        public void setProteinList(List<Protein> proteinList) {
            this.proteinList = proteinList;
        }

        public void addProtein(Protein protein) {
            proteinList.add(protein);
        }

        public int hashCode() {
            return key.hashCode();
        }

    }



    public static void fillRetentionTimes(String path, Map<String,  Map<Integer,PeptideContainer>> ms2ContainerMap) throws IOException {
        for(Map.Entry<String, Map<Integer,PeptideContainer>> entry: ms2ContainerMap.entrySet())
        {
            String filename = entry.getKey();
            Map<Integer,PeptideContainer> containerMap = entry.getValue();
            String fullpath = path+ File.separator+filename+".ms2.index";
            BufferedReader br = new BufferedReader(new FileReader(fullpath));
            String line;
            PeptideContainer containerToFill = null;
            while((line = br.readLine())!=null)
            {
                String [] arr = line.split("\t");
                int scanNo = Integer.parseInt(arr[0]);
                containerToFill = containerMap.get(scanNo);
                if(containerToFill!=null)
                {
                    double retTime = Double.parseDouble(arr[3]);
                    containerToFill.setRetentionTime(retTime);
                }
            }
            br.close();
        }

    }


    public  Map<String,PeptideContainer> readDTASelect(String dtaPath) throws IOException {
        DTASelectFilterReader reader = new DTASelectFilterReader(dtaPath);
        Map<String, PeptideContainer> peptideContainerMap = new HashMap<>();
        Map<String, Map<Integer, PeptideContainer> > ms2ScanMap = new HashMap<>();
        for(Iterator<Protein> proteinIterator = reader.getProteins(); proteinIterator.hasNext(); )
        {
            Protein protein = proteinIterator.next();
            String locus = protein.getLocus();
            if(!locus.startsWith("Reverse_"))
            {

                for(Iterator<Peptide> peptideIterator = protein.getPeptides(); peptideIterator.hasNext(); )
                {
                    Peptide peptide = peptideIterator.next();
                    String sequence = peptide.getSequence();
                    String cs = peptide.getChargeState();
                    int scanNumber = peptide.getScanNumber();
                    String ms2Filename = peptide.getFileName();
                    String key = sequence+cs;


                    PeptideContainer container = peptideContainerMap.getOrDefault(key,
                            new PeptideContainer(key,sequence,cs,scanNumber,ms2Filename,dtaPath));
                    container.addProtein(protein);
                    peptideContainerMap.put(key,container);
                    Map<Integer, PeptideContainer> containerSet = ms2ScanMap.getOrDefault(ms2Filename,new HashMap<>());
                    containerSet.put(scanNumber, container);
                    ms2ScanMap.put(ms2Filename, containerSet);
                }
            }

        }
        reader.close();

        String dir = dtaPath.replaceAll("DTASelect-filter.txt","../../spectra/");
        fillRetentionTimes(dir, ms2ScanMap);
        return peptideContainerMap;
    }


    private void init() throws Exception {
        MSSplitFolderCreation msp = new MSSplitFolderCreation();
        isoReader = new IsotopeReader(conf.getRootConfEle());
        List<SampleModel> sampleModelList = conf.getSampleList();
        for (SampleModel eachSample : sampleModelList) {
            List<String> pathList = new ArrayList<>();
            for (String path : eachSample.getPathList()) {

                if(!path.endsWith("/"))
                    path += "/";

                String spectraPath = path + "../../spectra/";
                String splitSpectraPath = path + "../../spectra/split/";

                splitSpectraMap.putAll( msp.splitMS1Files(spectraPath, CensusConstants.LABELFREE_MS1_SPLIT_SCAN_NUM, true, conf.isLabelfreeCheckChargeState()) );
                splitMs1FileHt.putAll( ChroGenerator.createIndexedFiles(splitSpectraPath, CensusConstants.MS1_FILE) );
                ms2ToMs1Map.putAll( IndexUtil.buildMS2toMS1ScanMapFiles(spectraPath) );

                Hashtable<String, IndexedFile> ht = ChroGenerator.createIndexedFiles(spectraPath, "ms1", true,true);
                indexedMs1FileMap.putAll(ht);
                pathList.addAll(ht.keySet());
            }
            eachSample.setSpectraList(pathList);
        }
    }



    public  void test1Peptide() throws Exception {
        String configFile = "/home/yateslab/project_data/census/2805AccurateMassTag/census_config_lfree.xml";
        Configuration conf = Configuration.getInstance();

        if (!conf.isReadConfigFile()) {
            conf.setLabelfree(true);
            conf.readXMLParam(configFile);
        }
        //Hashtable<String, IndexedFile> indexHt = conf.getIndexHt();
        String [] arr  = {"/data/2/rpark/ip2_data/rpark/Jolene_Hela_loading/HeLa_100ng_BEH60_140_35_IC_DE5_5e3_2_2015_03_19_17_30921/search/projects2015_03_23_12_78271/../../spectra/20141026_HeLa_100ng_BEH60_140min_35ms_IC_DE5_5e3_2.ms1",
        "/data/2/rpark/ip2_data/rpark/Jolene_Hela_loading/HeLa_100ng_BEH60_140_35_IC_DE5_5e3_3_2015_03_19_17_30922/search/projects2015_03_23_12_78273/../../spectra/20141026_HeLa_100ng_BEH60_140min_35ms_IC_DE5_5e3_3.ms1"};
        double startRt = 96.566-conf.getRetentionTimeWindow();
        double endRt = 96.566+conf.getRetentionTimeWindow();
        String sequence = "K.MITGDSQETAVAIASR.L";
        int cs = 2;
        IsotopeReader isoReader = new IsotopeReader(conf.getRootConfEle());
        for(String eachKey: arr)
        {
            IndexedFile origIFile = indexedMs1FileMap.get(eachKey);
            TDoubleIntHashMap retentonToScanMap = origIFile.getRetentonToScanMap();

            int startScan = retentonToScanMap.get(startRt);
            int endScan = retentonToScanMap.get(endRt);

            double[] retKeys = origIFile.getRtArr();
            int startIndex = BinarySearch.binarySearch(retKeys, startRt);

            if (startScan <= 0) {

                double rtTime = retKeys[startIndex];

                startScan = retentonToScanMap.get(rtTime);
            }
            int endIndex = BinarySearch.binarySearch(retKeys, endRt);
            if (endScan <= 0) {

                double rtTime = retKeys[endIndex];
                endScan = retentonToScanMap.get(rtTime);
            }

            GaussianPeakModel peakModel = isotopeCalc(startScan, endScan, startIndex, endIndex, isoReader, sequence,
                    cs, origIFile, splitSpectraMap, splitMs1FileHt);
            int[] scanArr = peakModel.getScanArr();
            double[] retArr = peakModel.getRetArr();
            double[] peakArr = peakModel.getPeakArr();
        }

    }

    private List<String> getDTAsFromSampleList()
    {
        List<SampleModel> sampleModelList = conf.getSampleList();
        List<String> pathList = new ArrayList<>();
        for(SampleModel model : sampleModelList)
        {
            String path = model.getPathList().get(0)+File.separator+"DTASelect-filter.txt";
            pathList.add(path);
        }
        return pathList;
    }


    public GaussianPeakModel generatePeakModel(String spectraPath, double retTimeLow, double retTimeHigh, String sequence, int cs ) throws Exception {
        IndexedFile origIFile = indexedMs1FileMap.get(spectraPath);
        TDoubleIntHashMap retentonToScanMap = origIFile.getRetentonToScanMap();

        int startScan = retentonToScanMap.get(retTimeLow);
        int endScan = retentonToScanMap.get(retTimeHigh);

        double[] retKeys = origIFile.getRtArr();
        int startIndex = BinarySearch.binarySearch(retKeys, retTimeLow);

        if (startScan <= 0) {

            double rtTime = retKeys[startIndex];

            startScan = retentonToScanMap.get(rtTime);
        }
        int endIndex = BinarySearch.binarySearch(retKeys, retTimeHigh);
        if (endScan <= 0) {

            double rtTime = retKeys[endIndex];
            endScan = retentonToScanMap.get(rtTime);
        }

        GaussianPeakModel peakModel = isotopeCalc(startScan, endScan, startIndex, endIndex, isoReader, sequence,
                cs, origIFile, splitSpectraMap, splitMs1FileHt);
        return  peakModel;
    }

    public void run() throws Exception {

        List<SampleModel> sampleModelList = conf.getSampleList();
        IsotopeReader isoReader = new IsotopeReader(conf.getRootConfEle());

        List<Integer> countList = new ArrayList<>();
        List<Map<String,PeptideContainer>> peptideContainerMapList = new ArrayList<>();
        List<Set<String>> setList = new ArrayList<>();
        List<Map<String,List<PeptideContainer>>> newPeptideMapList = new ArrayList<>();
        Set<String> peptideUnionSet = new HashSet<>();

        List<String> dtaPaths = getDTAsFromSampleList();
        for(String path : dtaPaths)
        {
            Map<String,PeptideContainer> peptideContainerMap = readDTASelect(path);
            peptideContainerMapList.add(peptideContainerMap);
            countList.add(peptideContainerMap.size());
            peptideUnionSet.addAll(peptideContainerMap.keySet());
            setList.add(new HashSet<>(peptideContainerMap.keySet()));
            newPeptideMapList.add(new HashMap<>());
        }


        int count =0;
        double oldPercentDone = 0.0;
        for(String peptide : peptideUnionSet)
        {
            double retTimeLow = Double.MAX_VALUE;
            double retTimeHigh = Double.MIN_VALUE;
            String sequence = null;
            int cs =-1;
            List<Protein> proteinList = null;
            List<Integer> sampleIndexToSearch = new ArrayList<>();
            for(int i=0; i<peptideContainerMapList.size(); i++)
            {
                PeptideContainer container = peptideContainerMapList.get(i).get(peptide);
                if(container!=null)
                {
                    sequence = container.sequence;
                    cs = container.csInt;
                    double retTime = container.getRetentionTime();
                    double low = retTime-conf.getRetentionTimeWindow();
                    double high = retTime+conf.getRetentionTimeWindow();
                    retTimeLow = low < retTimeLow ? low : retTimeLow;
                    retTimeHigh = high > retTimeLow ? high: retTimeHigh;
                    proteinList = container.getProteinList();
                }
                else
                {
                    sampleIndexToSearch.add(i);
                }
            }
            for(int index: sampleIndexToSearch)
            {
                String name = sampleModelList.get(index).getSpectraList().get(0);
                GaussianPeakModel peakModel = generatePeakModel(name,retTimeLow,retTimeHigh,sequence,cs);
                if(peakModel.isHasPeak())
                {
                    countList.set(index,countList.get(index)+1);

                    setList.get(index).add(peptide);
                    PeptideContainer temp = new PeptideContainer(peptide, sequence, Integer.toString(cs), 0, name, name);
                    temp.setProteinList(proteinList);
                    for(Protein p: proteinList)
                    {
                        List<PeptideContainer> list = newPeptideMapList.get(index).getOrDefault(p.getLocus(),new ArrayList<>());
                        list.add(temp);
                        newPeptideMapList.get(index).put(p.getLocus(), list);
                    }
         //           System.out.println(">>>>\t"+index+"\t" +peptide+"\t"+name);
                }
            }
            double percentDone = (double)count/(double)peptideUnionSet.size()*100;
            if(percentDone - oldPercentDone > 5.0 )
            {
                oldPercentDone = percentDone;
                System.out.println(">>peptide number\t"+count+"\t"+percentDone+"% done");
            }
            count++;
        }
        Set<String> intersection = new HashSet<>(setList.get(0));
        for(int i=1; i<setList.size(); i++)
        {
            intersection.retainAll(setList.get(i));
        }

        int sizeIntersection = intersection.size();



        printSet(peptideUnionSet, path+File.separator+"peptideUnionSet.txt");

        System.out.println(">>> union size is: "+ peptideUnionSet.size());
        System.out.println(">>> intersection size is:" +sizeIntersection);
        for(int i=0 ; i<setList.size(); i++)
        {
            printSet(setList.get(i), path+File.separator+"peptideSetGroup_"+i+".txt");

        }

        for(int i=0; i<setList.size(); i++)
        {
            String dtaPath = dtaPaths.get(i);
            Map<String,List<PeptideContainer>> peptideSet = newPeptideMapList.get(i);
            writeNewDTASelect(dtaPath,peptideSet);
        }

    }

    private void printSet(Set<String> set, String path) throws IOException {
        BufferedWriter bw = new BufferedWriter(new FileWriter(path));
        for(String s : set)
        {
            bw.write(s);
            bw.newLine();
        }
        bw.close();
    }

    private void writeNewDTASelect(String path, Map<String,List<PeptideContainer>> containerMap) throws IOException {
        String outputPath = path +".filled";
        BufferedReader br = new BufferedReader(new FileReader(path));
        BufferedWriter bw = new BufferedWriter(new FileWriter(outputPath));
        boolean readMode = false;
        boolean addMode = false;
        String line;
        List<PeptideContainer> containerList = new ArrayList<>();
        while((line = br.readLine())!=null)
        {
            if(!readMode && line.startsWith("Unique\t"))
            {
                readMode = true;
            }
            else if(readMode)
            {
                if(line.startsWith("\tProteins\t"))
                {
                    readMode = false;
                    if(containerMap.size()>0)
                    {

                    }
                }
                else if(line.matches("^[A-z].*"))
                {
                    if(addMode)
                    {
                        for(PeptideContainer s: containerList)
                        {
                            bw.write("\t");
                            String newName = cleanMs2Name(s.ms2Name);
                            bw.write(s.ms2Name+".XXXX.XXXX."+s.cs);
                            bw.write("\t1.0\t1.0\tx\tx\tx\tx\tx\tx\tx\tx\tx\tx\t");
                            bw.write(s.sequence);
                            bw.newLine();
                        }
                        addMode = false;
                    }
                    String arr [] = line.split("\t");
                    String locus = Fasta.getAccession( arr[0]);
                    System.out.println(">>>>"+locus);
                    if(containerMap.containsKey(locus))
                    {
                        System.out.println(">>>>matches : "+locus);
                        containerList = new ArrayList<>();
                        containerList.addAll(containerMap.get(locus));
                        addMode = true;
                        containerMap.remove(locus);
                    }

                }
            }
            bw.write(line);
            bw.newLine();
        }
        br.close();
        bw.close();

    }

    public static String cleanMs2Name(String ms2name)
    {
        int index = ms2name.lastIndexOf("/");
        String temp = index>0 ? ms2name.substring(0,index) : ms2name;
        String temp2 = temp.replaceAll("\\.ms[12]","" );
        return temp2;
    }






    public static void main(String [] args) throws Exception {
        String configFile = args[0];
        AccurateMassTag massTag = new AccurateMassTag(configFile);
        massTag.run();
        //AccurateMassTag massTag = new AccurateMassTag(configFile);
        //massTag.test1Peptide();;

    }
}
