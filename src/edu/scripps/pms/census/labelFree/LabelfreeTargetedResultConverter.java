package edu.scripps.pms.census.labelFree;

import edu.scripps.pms.census.conf.Configuration;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;

import java.io.*;
import java.util.*;

/**
 * Created by Titus Jung titusj@scripps.edu on 8/17/18.
 */
public class LabelfreeTargetedResultConverter {

    private List<Integer> aucIndices = new ArrayList<>();
    private List<Integer> peakArrayIndices = new ArrayList<>();
    private List<Integer> isotopeArrayIndices = new ArrayList<>();
    private int seqLocation;
    private int csLocation;
    private int numColumns;


    private void readHeader(String header)
    {
        int num=0;
        String [] arr = header.split("\t");
        for(String s: arr)
        {
            if(s.startsWith("AUC"))
            {
                aucIndices.add(num);
            }
            else if(s.startsWith("isotope_Array"))
            {
                isotopeArrayIndices.add(num);
            }
            else if(s.startsWith("peaks"))
            {
                peakArrayIndices.add(num);
            }
            else if(s.startsWith("sequence"))
            {
                seqLocation = num;
            }
            else if(s.startsWith("chargeState"))
            {
                csLocation = num;
            }
            num++;
        }
        numColumns = arr.length;
    }

    public void convert(String input, String output) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(input));
        BufferedWriter bw = new BufferedWriter(new FileWriter(output));
        String line;
        boolean readHeader = false;
        while ((line = br.readLine()) != null)
        {
            if(readHeader)
            {
                String [] arr = line.split("\t");

                List<String> columns = new ArrayList<>();
                for(String s: arr)
                {
                    columns.add(s);
                }
                if(columns.size()==numColumns-1)columns.add("\t");
                List<String> aucList = new ArrayList<>();
                List<String> peakList = new ArrayList<>();
                List<String> isotopeList = new ArrayList<>();
                String seq = arr[seqLocation];
                String cs = arr[csLocation];
                for(int i: aucIndices)
                {
                    aucList.add(columns.get(i));
                }
                for(int i: peakArrayIndices)
                {
                    peakList.add(columns.get(i));
                }
                for(int i: isotopeArrayIndices)
                {
                    isotopeList.add(columns.get(i));
                }
                bw.append(">Seq= ").append(seq);
                bw.newLine();
                bw.append("Z= ").append(cs);
                bw.newLine();
                for(int i=0;i< aucList.size(); i++)
                {
                    bw.append("AUC").append(Integer.toString(i)).append("= ").append(aucList.get(i));
                    bw.newLine();
                }
                for(int i=0;i< isotopeList.size(); i++)
                {
                    bw.append("Isotope_array_").append(Integer.toString(i)).append("= ").append(isotopeList.get(i));
                    bw.newLine();
                }
                for(int i=0;i< peakList.size(); i++)
                {
                    bw.append("Peaks").append(Integer.toString(i)).append("= ").append(peakList.get(i));
                    bw.newLine();
                }
            }
            else if(line.startsWith("H\tsequence"))
            {
                readHeader(line);
                readHeader = true;
            }
            else
            {
                bw.write(line);
                bw.newLine();
            }
        }
        br.close();
        bw.close();
    }

    public void convertJSON(String inputPath, String outputPath, String configFile) throws Exception {

        JSONParser parser=new JSONParser();
        Object obj = parser.parse(new FileReader(inputPath));
        JSONObject jsonObject =  (JSONObject) obj;
        JSONArray peptideLst = (JSONArray) jsonObject.get("peptideList");
        BufferedWriter bw = new BufferedWriter(new FileWriter(outputPath));

        for(int j=0; j<peptideLst.size(); j++)
        {

                JSONArray lst=(JSONArray)(peptideLst.get(j));
            String sequence = "";
            Map<String,List<String>> aucMap = new HashMap<>();
            Map<String,List<String>> fileNameMap = new HashMap<>();
            Map<String,List<String>> rtMap = new HashMap<>();
            Map<String,List<String>> peakMap = new HashMap<>();
            Map<String,List<String>> isoMap = new HashMap<>();
            List<String> isotopeList;

            List<String> aucList ;
            List<String> fileNamelist;
            Set<String> sampleSet = new TreeSet<>();
            List<String> rtList;
            String cs = "";
            List<String> peakList;
            for(int k=0;k<lst.size();k++){
                JSONObject lstObj= (JSONObject) lst.get(k);
                JSONArray isotopeArray=(JSONArray) lstObj.get("isotopeArr");
                String auc = (String)lstObj.get("auc");

                sequence = (String)lstObj.get("seq");

                String fileName = (String)lstObj.get("file");
                String sampleName = (String)lstObj.get("sampleName");
                String groupName = (String)lstObj.get("sampleName");

                sampleSet.add(sampleName);
                cs = (String)lstObj.get("charge");

                String key =groupName+sampleName+fileName;


                String startRt = (String)lstObj.get("startRt");
                String endRt = (String)lstObj.get("endRt");
                String rtRange = startRt+":"+endRt;

                // if(j%2==0) {
                StringBuilder isoArrBuilder = new StringBuilder();
                Iterator<java.lang.Double> iterator = isotopeArray.iterator();
                while (iterator.hasNext()) {
                    isoArrBuilder.append(iterator.next() + ",");
                }
                String isoArrayString = isoArrBuilder.toString();
                String peaks = (String)lstObj.get("peaks");


                if(peakMap.get(key) == null)
                {
                    List<String> temp = new ArrayList<String>();
                    temp.add(peaks);
                    peakMap.put(key,temp);

                    temp = new ArrayList<String>();
                    temp.add(auc);
                    aucMap.put(key,temp);

                    temp = new ArrayList<String>();
                    temp.add(fileName);
                    fileNameMap.put(key,temp);

                    temp = new ArrayList<String>();
                    temp.add(rtRange);
                    rtMap.put(key,temp);

                    temp = new ArrayList<String>();
                    temp.add(rtRange);
                    rtMap.put(key,temp);

                    temp = new ArrayList<String>();
                    temp.add(isoArrayString);
                    isoMap.put(key,temp);
                }
                else
                {
                    peakMap.get(key).add(peaks);
                    aucMap.get(key).add(auc);
                    fileNameMap.get(key).add(fileName);
                    rtMap.get(key).add(rtRange);
                    isoMap.get(key).add(isoArrayString);
                }
            }
           // for(int i=0; )
            if(j==0)
            {
                Configuration conf = Configuration.getInstance();
                // conf.setLabelfreeCheckChargeState(true);

                if (!conf.isReadConfigFile()) {
                    conf.setLabelfree(true);
                    conf.readXMLParam(configFile);
                }

                List<org.jdom.Element> samGroupEleList1 = conf.getRootConfEle().getChildren("sample");
                for(int ii=0;ii<samGroupEleList1.size();ii++){
                    bw.append("H\tGROUP_NAME\t"+samGroupEleList1.get(ii).getAttributeValue("group")+"\t");
                    List<org.jdom.Element> sampleList = samGroupEleList1.get(ii).getChildren("each_sample");
                    for(org.jdom.Element sample : sampleList)
                    {
                        bw.append(sample.getAttributeValue("name"));
                        bw.append("\n");
                    }

                }
            }
            bw.append(">Seq= ").append(sequence);
            bw.newLine();
            bw.append("Z= ").append(cs);
            bw.newLine();
            StringBuilder aucSB = new StringBuilder();
            StringBuilder isotopeSB = new StringBuilder();
            int i=0;

            List<String> keyList = new ArrayList<>();
            keyList.addAll(aucMap.keySet());
            Collections.sort(keyList);
            for(String key: keyList)
            {
                //System.out.println(key);
                aucList = aucMap.get(key);
                peakList = peakMap.get(key);
                isotopeList = isoMap.get(key);
                aucSB.append("AUC").append(i).append("=");
                for(int k=0;k< aucList.size(); k++)
                {
                    aucSB.append(aucList.get(k)).append("[").append(peakList.get(k)).append("],");
                }
                aucSB.append("\n");

                isotopeSB.append("Isotope_array_").append(i).append("=");
                for(int k=0;k< aucList.size(); k++)
                {
                    isotopeSB.append(isotopeList.get(k));
                }
                isotopeSB.append("\n");

                i++;
            }


            bw.append(aucSB.toString());
            bw.append(isotopeSB.toString());



    /*
            for(int i=0;i< isotopeList.size(); i++)
            {
                bw.append("Isotope_array_").append(Integer.toString(i)).append("= ").append(isotopeList.get(i));
                bw.newLine();
            }
            for(int i=0;i< peakList.size(); i++)
            {
                bw.append("Peaks").append(Integer.toString(i)).append("= ").append(peakList.get(i));
                bw.newLine();
            }
*/

        }
        bw.close();




    }



    public static void main(String [] args) throws Exception {
        String input = args[0];
        String output = args[1];
        LabelfreeTargetedResultConverter converter = new LabelfreeTargetedResultConverter();

        String inputExtension = args[0].substring(args[0].indexOf("."),args[0].length()).toLowerCase();

        if(args.length>2)
        {
            String configPath = args[2];
            converter.convertJSON(input,output,configPath);
        }
        else
        {
            converter.convert(input,output);
        }

    }


    public LabelfreeTargetedResultConverter()
    {

    }

}
