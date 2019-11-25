/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package scripts;

import edu.scripps.pms.census.hash.MSIndexFileCreator;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;

/**
 *
 * @author rampuria
 */
public class MSSplitFolderCreation {

    public static void main(String[] args) throws Exception {
        String filename = "/data/2/rpark/ip2_data/yjwang/Eval_0203_human/A1_2016_11_15_09_86094/search/projects2016_11_17_05_105379/";
        //String filename = "/data/2/rpark/ip2_data/rpark/an_chi/Onbeads_Lung_DMSO_CS_fixed_2016_09_12_16_84798/spectra";
        MSSplitFolderCreation ms = new MSSplitFolderCreation();
        HashMap<String,String> map = ms.splitMS1Files(filename, 4);
        System.out.println(map);
    }

    public HashMap<String,String> splitMS1Files(String folder, int splitSize) throws Exception {

        return splitMS1Files(folder, splitSize, false, false);
    }

    public HashMap<String,String> splitMS1Files(String folder, int splitSize, boolean checkExistingIndex) throws Exception {

      return splitMS1Files(folder, splitSize, checkExistingIndex, false);
    }
    /*
    public HashMap<String,String> splitMS1Files(String folder, int splitSize, boolean checkExistingIndex) throws Exception {

        return splitMS1Files(folder, splitSize, checkExistingIndex);
    }
*/

    public HashMap<String,String> splitMS1Files(String folder, int splitSize, boolean checkExistingIndex, boolean chargeStateCheck) throws Exception {
        File[] fileList = null;
        if(folder.endsWith("/"))
            folder = folder.substring(0, folder.length()-1);

        HashMap<String,String> scantofileMap = new HashMap<>();
       if (new File(folder).isDirectory()) {

           if(folder.endsWith("spectra"))
                fileList = new File(folder).listFiles();


           File splitFolder = new File(folder + File.separator + "split");
           if(!splitFolder.exists()) splitFolder.mkdir();

           File versionFile = new File(folder + File.separator + "split/version2");
           File csFile = new File(folder + File.separator + "split/cs");  //charge state check file

           if(!checkExistingIndex || !versionFile.exists()) {
             if(chargeStateCheck) {
               if (!csFile.exists()) {
                 for (File file : fileList) {
                   if (file.getName().endsWith(".ms1")) {
                     run(file.getPath(), splitSize, scantofileMap);
                   }
                 }
                 System.out.println("generating split index files...");
                 MSIndexFileCreator.createIndexFileByPath(folder + File.separator + "split");

                 versionFile.createNewFile();
		 if(chargeStateCheck) csFile.createNewFile();
		   else csFile.delete();

               }
             } else {
               if (csFile.exists()) {
                 for (File file : fileList) {
                   if (file.getName().endsWith(".ms1")) {
                     run(file.getPath(), splitSize, scantofileMap);
                   }
                 }
                 System.out.println("generating split index files...");
                 MSIndexFileCreator.createIndexFileByPath(folder + File.separator + "split");

                 versionFile.createNewFile();

               } else if(splitFolder.listFiles().length<=0) {
                 for (File file : fileList) {
                   if (file.getName().endsWith(".ms1")) {
                     run(file.getPath(), splitSize, scantofileMap);
                   }
                 }

                   System.out.println("generating split index files...");
                   MSIndexFileCreator.createIndexFileByPath(folder + File.separator + "split");
                   versionFile.createNewFile();
		   if(chargeStateCheck) csFile.createNewFile();
		   else csFile.delete();

               }
             }


           }  else if( (csFile.exists() && !chargeStateCheck) || (!csFile.exists() && chargeStateCheck) ) {
                 for (File file : fileList) {
                   if (file.getName().endsWith(".ms1")) {
                     run(file.getPath(), splitSize, scantofileMap);
                   }
                 }

                   System.out.println("generating split index files...");
                   MSIndexFileCreator.createIndexFileByPath(folder + File.separator + "split");
                   versionFile.createNewFile();
		   if(chargeStateCheck) csFile.createNewFile();
		   else csFile.delete();

	   }

                fileList = new File(folder+ File.separator + "split").listFiles();
                BufferedReader br = null;
                for (File file : fileList) {
                    if (file.getName().endsWith(".index")) {
                        try{
                             br = new BufferedReader(new FileReader(file));
                            String eachLine = null;
                            while((eachLine=br.readLine()) != null){
                                String [] words = eachLine.split("\t");
                                String splitfileName = FilenameUtils.removeExtension(file.getName());
                                String originalFilename = splitfileName.substring(splitfileName.indexOf("-")+1);
                                scantofileMap.put(words[0]+"\t"+folder+File.separator+originalFilename, folder+File.separator+"split"+File.separator+splitfileName);
                            }

                            br.close();

                        }catch(Exception e){
                            e.printStackTrace();
                        }


                    }
                }
             //   return scantofileMap;


           // }
        } else {
            System.out.println("Need spectra directory");
        }

        System.out.println("indexing completed");
        return scantofileMap;
    }

  //  public HashMap<String,String> run(String filename, int splitnum, HashMap<String,String> scantofileMap) {
  public void run(String filename, int splitnum, HashMap<String,String> scantofileMap) {

/* try {

throw new Exception();
} catch(Exception e) {
e.printStackTrace();
System.exit(0);
} */
        File file = new File(filename);
      System.out.println("split folder creating on " + filename);
            try {
                int count = 0;
                BufferedReader br = new BufferedReader(new FileReader(file));
                BufferedWriter bw = null;
                String eachLine = br.readLine();
                List<Integer> scanList = new ArrayList<>();
                int firstscan =0;
                int lastscan =0;

                File testFile = new File(file.getParent() + File.separator +"split"+File.separator+"test.txt");
                if(!testFile.exists()) testFile.createNewFile();

                while (eachLine != null) {
                    if (eachLine.startsWith("H\t")) {
                        eachLine = br.readLine();
                    }
                    else if (eachLine.startsWith("S\t")) {
                        String [] words = eachLine.split("\t");
                        if(count == 0){
                            firstscan = Integer.parseInt(words[1]);
                            File splitFolder = new File(file.getParent() + File.separator +"split");
                            if(!splitFolder.exists()) splitFolder.mkdir();
                            else splitFolder.delete();

                            bw = new BufferedWriter(new FileWriter(file.getParent() + File.separator +"split"+File.separator+"test.txt"));
                        }
                        else{
                            lastscan= Integer.parseInt(words[1]);
                        }
                        scanList.add(Integer.parseInt(words[1]));
                        count++;
                        bw.write(eachLine+"\n");
                        while ((eachLine = br.readLine()) != null) {


                            if (eachLine.startsWith("S\t")) {
                                if(count == splitnum){
                                    for(int scan:scanList){

                                        scantofileMap.put(scan+"\t"+filename,file.getParent() + File.separator +"split"+File.separator+firstscan+"_"+lastscan+"-"+file.getName());
                                    }
                                    scanList.clear();
                                    count=0;
                                    new File(file.getParent() +File.separator+"split"+File.separator+"test.txt").renameTo(new File(file.getParent() + File.separator +"split"+File.separator+firstscan+"_"+lastscan+"-"+file.getName()));
                                    bw.close();
                                }
                                break;
                            }
                            bw.write(eachLine+"\n");

                        }
                        if(eachLine == null){
                            for(int scan:scanList){
                                        scantofileMap.put(scan+"\t"+filename,file.getParent() + File.separator +"split"+File.separator+firstscan+"_"+lastscan+"-"+file.getName());
                                    }
                                    scanList.clear();
                           new File(file.getParent() + File.separator+"split"+File.separator+"test.txt").renameTo(new File(file.getParent() + File.separator +"split"+File.separator+firstscan+"_"+lastscan+"-"+file.getName()));
                           bw.close();
                        }

                    } else {
                        eachLine = br.readLine();
                    }
                }
                br.close();

            } catch (Exception e) {
                e.printStackTrace();
            }

    }

}

