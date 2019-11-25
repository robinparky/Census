
package edu.scripps.pms.census;

import java.io.*;
import java.util.*;

public class CensusChroProgress {

    private String filePath;
    private Properties props;
    private OutputStream output;

    CensusChroProgress(String filePath, Map<String, Integer> groupCompletionProcess){
        this.filePath = filePath;
        init(groupCompletionProcess);
    }

    public void init(Map<String, Integer> groupCompletionProcess){
        try {
            output = new FileOutputStream(new File(filePath));
            props = new Properties();
            Set<String> groupNameSet = groupCompletionProcess.keySet();
            for (String groupName : groupNameSet) {
                Integer progress = groupCompletionProcess.get(groupName);
                props.setProperty(groupName, progress.toString());
            }
            synchronized (output) {
                props.store(output, null);
            }
        } catch (IOException io) {
            io.printStackTrace();
        } finally {
            if (output != null) {
                try {
                    output.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }

        }
    }

    public void updateProcess(String groupName, Integer progress){
        try {
            FileInputStream in = new FileInputStream(new File(filePath));
            props = new Properties();
            props.load(in);
            in.close();

            output = new FileOutputStream(new File(filePath));
            props.setProperty(groupName, progress.toString());
            synchronized (output) {
                props.store(output, null);
            }
        } catch (IOException io) {
            io.printStackTrace();
        } finally {
            if (output != null) {
                try {
                    output.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }

        }
    }

}