package edu.scripps.pms.census.util;

import java.util.Hashtable;

/**
 * Created by rampuria on 1/31/17.
 */
public class XYPoint {

    private double x;

    private double y;

    private double size;

    private boolean isMatched;

    private String label = "";

    /* Store key-value data pairs*/
    private Hashtable<String, String> data = new Hashtable<String, String>();

    public XYPoint(int x, int y, String label) {
        this.x = x;
        this.y = y;
        this.label = label;
    }
    public XYPoint(){}

    public XYPoint(double x, int y, String label) {
        this.x =x;
        this.y = y;
        this.label = label;
    }


    public XYPoint(double x, int y) {
        this.x =x;
        this.y = y;
    }

    public boolean isIsMatched() {
        return isMatched;
    }

    public void setIsMatched(boolean isMatched) {
        this.isMatched = isMatched;
    }

    public double getX() {
        return x;
    }

    public void setX(double x) {
        this.x = x;
    }

    public double getY() {
        return y;
    }

    public void setY(double y) {
        this.y = y;
    }

    public double getSize() {
        return size;
    }

    public void setSize(double size) {
        this.size = size;
    }

    public void setLabel(String label) {
        this.label = label;
    }

    public String getLabel() {
        return label;
    }

    public String getData(String key) {
        return data.get(key);
    }

    public void setData(String key, String value){
        data.put(key, value);
    }
}
