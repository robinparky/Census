package edu.scripps.pms.census.util;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by rampuria on 1/31/17.
 */
public class JSONData {

    private String key = "";

    private String color = "";

    private List<XYPoint> values = new ArrayList<>();

    public JSONData() {

    }

    public JSONData(String color) {
        this.color = color;
    }

    public String getKey() {
        return key;
    }

    public void setKey(String key) {
        this.key = key;
    }

    public String getColor() {
        return color;
    }

    public void setColor(String color) {
        this.color = color;
    }

    public List<XYPoint> getValues() {
        return values;
    }

    public void setValues(List<XYPoint> values) {
        this.values = values;
    }
}
