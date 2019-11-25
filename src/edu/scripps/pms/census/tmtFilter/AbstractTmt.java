/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.scripps.pms.census.tmtFilter;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 *
 * @author Harshil
 */
public class AbstractTmt {
    Map<String, List<String>> groupKeyMap = new LinkedHashMap();//GruoupName -> List<Key(header_Line_name)>
    public static ArrayList<String> intensityNames = new ArrayList();
    String header_line[] = new String[150];
    boolean isNormalized = true;
    String headerSline[] = new String[150];

    public void groupNameToKeys(Map groupNameMap) {

        if (!headerSline[1].equalsIgnoreCase("SLINE")) {
            return;
        }
        List keys = new ArrayList();
        keys.addAll(groupNameMap.keySet());
        for (int z = 0; z < groupNameMap.size(); z++) {
            List ar1 = (ArrayList<String>) groupNameMap.get(keys.get(z));
            Set keyArray = new LinkedHashSet();
            for (int j = 0; j < ar1.size(); j++) 
            {
                for (int i = 0; i < headerSline.length; i++) 
                {
                    String key = (String) ar1.get(j);
//                    if (key.contains(".")) 
//                        key = key.substring(0, key.indexOf("."));

                    if (isNormalized) 
                    {
                        if (headerSline[i].toLowerCase().contains("norm")
                                && headerSline[i].toLowerCase().contains("m/z")
                                && headerSline[i].toLowerCase().contains(key)) 
                        {
                            keyArray.add(headerSline[i]);
                        }
                    }
                    else {
/*
//                        if (!headerSline[i].toLowerCase().contains("norm")
//                                && headerSline[i].toLowerCase().contains("m/z")
//                                && headerSline[i].toLowerCase().contains((String) ar1.get(j))) {
//
//                        }
//                        else
//                        {
*/
                            if (!headerSline[i].toLowerCase().contains("norm")
                                    && headerSline[i].toLowerCase().contains("m/z")
                                    //&&  headerSline[i].toLowerCase().contains((String)ar1.get(j)+".")  )
                                    && headerSline[i].toLowerCase().contains(key)) {
                                keyArray.add(headerSline[i]);
                            }
 //                       }
                    }
                }
            }
            groupKeyMap.put((String) keys.get(z), new ArrayList<String>(keyArray));
        }
        System.out.println("Group MAped successfully..............");
    }

    public boolean checkValidEntryForpeptide(Peptide pep_detail, int newScanNumber, String newFileName) {
        ArrayList scanNumber = pep_detail.getScannum();
        Set fileName = pep_detail.getFilename();
        if (scanNumber.contains(newScanNumber) && fileName.contains(newFileName)) {
            return false;
        } else {
            return true;
        }
    }
    
    

}
