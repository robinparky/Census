/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.scripps.pms.go;

import edu.scripps.pms.census.util.dtaselect.Protein;
import edu.scripps.pms.census.util.io.DTASelectFilterReader;
import edu.scripps.pms.util.seq.Fasta;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Harshil
 */
public class MapUniprotToGo 
{
    private Map uniprotToGo= new HashMap<>();
    public static void main(String arg[])
    {
        String path= "C:\\Users\\Harshil\\Documents\\NetBeansProjects\\inputData";
        MapUniprotToGo mapUniproToGo = new MapUniprotToGo();
//        mapUniproToGo.startMapping(path);
        mapUniproToGo.readGoFile(path);   
         System.out.println()       ;
            
    }
    public void startMapping(String path)
    {
        try {
            DTASelectFilterReader reader = new DTASelectFilterReader(path + File.separator + "DTASelect-filter.txt");
            System.out.println("plese Wait.....");
            Iterator<Protein> pitr = reader.getProteins();
            String accession;
            while (pitr.hasNext()) {
                Protein p = pitr.next();
                accession = p.getAccession().split("\\|")[1];
//                 readGoFile(path);

            }
        } 
        catch (IOException ex) {
            Logger.getLogger(MapUniprotToGo.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    public void readGoFile(String path)
    {
        BufferedReader br = null;
        try 
        {
            File f = new File(path+ File.separator + "gene_association.goa_human" );
            br = new BufferedReader(new FileReader(f));
            String currentLine = br.readLine();
            String words[] = null;
            String uniprotName1=null;
            String uniprotName2=null;
            String goName=null;
            while(currentLine != null)
            {
                if(currentLine.startsWith("UniProt"))
                {
                   words = currentLine.split("\t");
                   uniprotName1=words[1];
                   if(words[7].contains(":"))
                       uniprotName2=words[7].split(":")[1];
                   else
                       uniprotName2="";
//                    System.out.println(uniprotName2);
                    goName = words[4].split(":")[1];
                    if(uniprotName1.equals(uniprotName2))
                    {
                        uniprotToGo.put(uniprotName1, goName);
                    }
                    else if(!uniprotName2.equals(""))
                    {
                        uniprotToGo.put(uniprotName1, goName);
                        uniprotToGo.put(uniprotName2, goName);
                    }
                    
                            
                }
                currentLine = br.readLine();
            }
            System.out.println("Success");
        } catch (FileNotFoundException ex) {
            Logger.getLogger(MapUniprotToGo.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(MapUniprotToGo.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            try {
                br.close();
            } catch (IOException ex) {
                Logger.getLogger(MapUniprotToGo.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        
    }

    /**
     * @return the uniprotToGo
     */
    public Map getUniprotToGo() {
        return uniprotToGo;
    }

    /**
     * @param uniprotToGo the uniprotToGo to set
     */
    public void setUniprotToGo(Map uniprotToGo) {
        this.uniprotToGo = uniprotToGo;
    }
}
