package edu.scripps.pms.census.model;

import java.util.*;

/**
 * Created by rpark on 12/28/16.
 */
public class SampleGroup {
    private List<SampleModel> sampleModelList = null;
    private String name;
    public final int ID;

    public SampleGroup()
    {
        ID = IDCOUNTER++;
    }


    public List<SampleModel> getSampleModelList() {
        return sampleModelList;
    }

    public void setSampleModelList(List<SampleModel> sampleModelList) {
        this.sampleModelList = sampleModelList;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public void addSample(SampleModel sampleModel) {
        if(sampleModelList == null)
            sampleModelList = new ArrayList<>();

        sampleModelList.add(sampleModel);
    }

    private static int IDCOUNTER = 0;

}
