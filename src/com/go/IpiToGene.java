package com.go;


import edu.scripps.pms.util.seq.Fasta;
//import edu.scripps.pms.mspid.ProteinDatabase;
import edu.scripps.pms.census.util.io.FastaReader;
import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Iterator;
import java.util.List;
import java.util.LinkedList;
import java.util.HashSet;
import java.io.IOException;
import java.io.FileReader;
import java.io.FileInputStream;
import java.io.PrintStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;


/**
 * @author  Tao Xu
 * @version $Id
 */
public class IpiToGene {
    public static String USAGE = "Usage: java -Xmx500M IPIId2GeneSymbol IPI_accession_file fasta_file";
    private static HashSet<String> ipiacs = new HashSet<String>(1000000); // swissprot of trembl id
//    private static String fastadb = "/lustre/people/applications/yates/dbase/EBI-IPI_mouse_3.52_11-21-2008_reversed.fasta";
    private static String fastadb = "C:\\Users\\Harshil\\Desktop\\Go\\EBI-IPI_rat_3.30_06-28-2007_reversed.fasta";
    private static String fileName = "C:\\Users\\Harshil\\Desktop\\Go\\DTASelect-filter-ipi.txt";
    private static String fileType;
    private static HashMap<String, String> ipiac2genesymbol = new HashMap<String, String>();
    private static HashMap<String, String> ipiac2genedescription = new HashMap<String, String>();
    private static HashMap<String, String> ipiTogene = new HashMap<String, String>();
    
    public IpiToGene (String fileName)
    {
        this.fileName = fileName;
        this.fileType = checkFilename();
        parseFile();
    }
    public static void parseFile()
    {
        if(fileType.equals("IPI"))
            parseIpiToGene(fileName);
//        else if(fileType.equals("UniProt"))
            
    }
    public static String checkFilename()
    {
        BufferedReader br = null;
        try {
            String filePath;
            String[] currentLine = new String[50];
            String lastWord;
            br = new BufferedReader(new FileReader(new File(fileName)));
            for(int i=0;i<3;i++)
            {
                filePath = br.readLine();
                currentLine = filePath.split("/");
                lastWord = currentLine[currentLine.length-1];
                if(lastWord.contains("EBI-IPI"))
                {
//                    fastadb = filePath;
                    return "IPI";
                }
                else if (lastWord.contains("UniProt"))
                {
//                    fastadb = filePath;
                    return "UniProt";
                }
                                        
            }
        } 
        catch (FileNotFoundException ex) {
            Logger.getLogger(IpiToGene.class.getName()).log(Level.SEVERE, null, ex);
        } 
        catch (IOException ex) {
            Logger.getLogger(IpiToGene.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            try {
                br.close();
              
            } catch (IOException ex) {
                Logger.getLogger(IpiToGene.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
          return "error";
    }
    public static void parseIpiToGene(String ipiFile)
    {
        try{ 
            getIdentifiedIds(); 
            int numgenesymbol = 0;
            int numipiac = 0;
            readDatabase();
            
            String currentlIne;
            StringBuffer sb= new StringBuffer();
            BufferedReader br = new BufferedReader(new FileReader(new File(fileName)));
            while((currentlIne = br.readLine() )!= null)
            {
                if(currentlIne.startsWith("IPI"))
                {
                    String[] temp= currentlIne.split(":");
                    String[] extra= temp[1].split("|");
                }
            }
            
  /*        Iterator<Fasta> itr = FastaReader.getFastas(new FileInputStream(fastadb));
            
            while(itr.hasNext()){
           // for (Iterator<Fasta> itr = FastaReader.getFastas(new FileInputStream(fastadb)); itr.hasNext();) { 
                Fasta f = itr.next();
                String ipiac = f.getAccession();
                
                String ipiacnoversion = ipiac.split("\\.")[0];
                
                
                System.out.println("ipiac: " + ipiac + "\tipiacnoversion: " + ipiacnoversion);
                if(ipiacs.contains(ipiac) || ipiacs.contains(ipiacnoversion)) {
System.out.println("found ipiac: " + ipiac + "\tipiacnoversion: " + ipiacnoversion);
                    String defline = f.getDefline();
                    String [] arr = defline.split("Gene_Symbol="); 
System.out.println("arr.lenth: " + arr.length +"\tdefline: " + defline);
                    if(arr.length > 1) {
                        String [] arr2 = arr[1].split(" ");
                        String genesymbol = arr2[0];

                        //System.out.println("Gene_Symbol for " + ipiac + " is " + genesymbol);
                        ipiTogene.put(ipiac, genesymbol);
//                        ipiac2genesymbol.put(ipiacnoversion, genesymbol);
                        StringBuffer sb = new StringBuffer();
                        for(int i = 1; i < arr2.length; i++) {
                            String s = arr2[i];
                            if(s.startsWith("similar") || s.startsWith("Similar")) {
                                i++; i++; // ignore similar to
                            }
                            sb.append(arr2[i]);
                            sb.append(" ");
                        }
                        ipiac2genedescription.put(ipiac, sb.toString());
                        ipiac2genedescription.put(ipiacnoversion, sb.toString());
                    }
            
                }
            }

            output();   
 */           
        }
        catch (IOException ex) {
            Logger.getLogger(IpiToGene.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    public static void readDatabase()
    {
        try {
            Iterator<Fasta> itr = FastaReader.getFastas(new FileInputStream(fastadb));

            while (itr.hasNext()) {
                Fasta f = itr.next();
                String ipiac = f.getAccession();
                String defline = f.getDefline();
                String[] arr = defline.split("Gene_Symbol=");
                System.out.println("arr.lenth: " + arr.length + "\tdefline: " + defline);
                if (arr.length > 1) {
                    String[] arr2 = arr[1].split(" ");
                    String genesymbol = arr2[0];
                    ipiTogene.put(ipiac, genesymbol);
                }
                
            }
        } catch (FileNotFoundException ex) {
            Logger.getLogger(IpiToGene.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(IpiToGene.class.getName()).log(Level.SEVERE, null, ex);
        }
          
            
        
    
    }
    public static void main(String args[]) throws IOException 
    {
        
          IpiToGene obj = new IpiToGene("C:\\Users\\Harshil\\Desktop\\Go\\DTASelect-filter-ipi.txt");
          obj.checkFilename();
          System.out.println("success");
    }

    public static void output() throws IOException {
       
        String prefix = fileName.substring(0, fileName.lastIndexOf("."));
        PrintStream ps = new PrintStream(prefix + "_with_gene_symbol.txt");
        BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(fileName)));
        String line = null;
        while ((line = br.readLine()) != null) 
        {
            
            String ipiac = line.trim();
            
//System.out.println("line: " + ipiac + "\tipiac2genesymbol.szie: " + ipiac2genesymbol.size());
            String genesymbol = ipiac2genesymbol.get(ipiac);
            String description = ipiac2genedescription.get(ipiac); 
//System.out.println("genesymbol: " + genesymbol + "\tdesc: " + description);
            if(genesymbol != null) {
//System.out.println("Found genesymbol: " + genesymbol + "\tdesc: " + description);
                String [] arr = genesymbol.split(";");
                
                for(int i = 0; i < arr.length; i++) {
                    ps.println(ipiac + "\t" + arr[i] + "\t" + description);
                }
            } else {
                ps.println(ipiac + "\t" + "" + "\t" + description);
            } 
        }
        br.close();
        ps.close();
//System.out.println("Number of ids: " + ipiacs.size());
    }
    public static void getIdentifiedIds() throws IOException {


        BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(fileName)));
        //br.readLine(); // remove the first line
        String line = null;
        while ((line = br.readLine()) != null) 
        {
            
            ipiacs.add(line.trim());
            
        }
        br.close();
//System.out.println("Number of ids: " + ipiacs.size());
    }
}
