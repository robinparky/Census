/*
 * SelectFileModel.java
 *
 * Created on July 24, 2006, 2:34 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package edu.scripps.pms.census.model;

/**
 *
 * @author rpark
 */
public class SelectFileModel
{
    private String sampleName;
    private String spectraFileName;

    public SelectFileModel(String sampleName, String spectraFileName)
    {
        this.sampleName = sampleName;
        this.spectraFileName = spectraFileName;
    }

    public String getSampleName() {
        return sampleName;
    }

    public void setSampleName(String sampleName) {
        this.sampleName = sampleName;
    }

    public String getSpectraFileName() {
        return spectraFileName;
    }

    public void setSpectraFileName(String spectraFileName) {
        this.spectraFileName = spectraFileName;
    }

    public String toString()
    {
        return this.spectraFileName;
    }

}
