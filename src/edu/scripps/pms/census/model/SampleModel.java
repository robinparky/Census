/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package edu.scripps.pms.census.model;

import java.util.*;

/**
 *
 * @author rpark
 */

//one sample represents this object containing multiple files including path
public class SampleModel {

    private String sampleName;
    private List<String> pathList = new ArrayList<String>();
    private List<String> labelfreeFilenameList = new ArrayList<>();
    private List<String> spectraList = new ArrayList<>();

    public SampleModel(String sampleName) {
        this.sampleName = sampleName;
    }
    public String getSampleName() {
        return sampleName;
    }

    public void setSampleName(String sampleName) {
        this.sampleName = sampleName;
    }

    public void addPath(String path) {
        this.getPathList().add(path);
    }

    public List<String> getPathList() {
        return pathList;
    }

    public void setPathList(List<String> pathList) {
        this.pathList = pathList;
    }

  /*  public String getLabelfreeFilename() {
        return labelfreeFilename;
    }

    public void setLabelfreeFilename(String labelfreeFilename) {
        this.labelfreeFilename = labelfreeFilename;
    }
*/

    public List<String> getLabelfreeFilenameList() {
        return labelfreeFilenameList;
    }

    public void setLabelfreeFilenameList(List<String> labelfreeFilenameList) {
        this.labelfreeFilenameList = labelfreeFilenameList;
    }

    public void addLabelfreeFilename(String name) {
        this.labelfreeFilenameList.add(name);
    }

    public List<String> getSpectraList() {
        return spectraList;
    }

    public void setSpectraList(List<String> spectraList) {
        this.spectraList = spectraList;
    }
}
