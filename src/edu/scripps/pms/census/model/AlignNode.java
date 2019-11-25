/*
 * AlignNode.java
 *
 * Created on May 8, 2006, 6:21 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package edu.scripps.pms.census.model;

import java.util.*;
/**
 *
 * @author $Id: AlignNode.java,v 1.1 2006/10/02 22:00:08 rpark Exp $ 
 */
public class AlignNode {
    
    private int x;
    private int y;
    private int i;
    private int j;
    private double localScore;   
    private double cumScore;
    private int pathIndex;
    private int halfBandWidth;
    
    /** Creates a new instance of AlignNode */
    public AlignNode(int x, int y, int halfBandWidth) {
        
        // i, j for matrix array coordinate
        // x, y for chromatograph cooridate
        
        this.setX(x);
        this.setY(y);

        if(y>=x)
        {
            i=y;
            j=y-x+halfBandWidth;
        }
        else
        {
            i=x;
            j=y-x+halfBandWidth;
        }
        
        this.halfBandWidth = halfBandWidth;
    }

    public int getX() {
        return x;
    }

    public void setX(int x) {
        this.x = x;
    }

    public int getY() {
        return y;
    }

    public void setY(int y) {
        this.y = y;
    }

    public double getLocalScore() {
        return localScore;
    }

    public void setLocalScore(double localScore) {
        this.localScore = localScore;
    }

    public double getCumScore() {
        return cumScore;
    }

    public void setCumScore(double cumScore) {
        this.cumScore = cumScore;
    }

    public int getPathIndex() {
        return pathIndex;
    }

    public void setPathIndex(int pathIndex) {
        this.pathIndex = pathIndex;
    }
    
    
    //convert from chromatograph coordinate to matrix array coordinate
    // n x n matrix -> n x hw array
    public static int[] rotateOrigin(int hw, int x, int y)
    {
	int i,j;

	if(y>=x)
	{
	    i=y;
	    j=y-x+hw;
	}
	else
	{
	    i=x;
	    j=y-x+hw;
	}
	
        
        int[]   temp = new int[2];
        

	temp[0] = i;

	if(hw*2<=j)
	    temp[1] = hw*2-1;
	else
	    temp[1] = j;

        return temp;
        
    }
    
    //convert from matrix array coordinate to chromatograph coordinate
    // n x hw array -> n x n matrix 
    public static int[] rotateNew(int hw, int i, int j)
    {
	
	int x;
	int y;
	int[]   temp = new int[2];
        
	if(hw>=j)
	{
	    x=i;
	    y=(i+j)-hw;
	}
	else
	{
	    //X=2*hw-1-j;
	    x=i-j+hw;
	    y=i;
	}

	temp[0] = x;
        temp[1] = y;
        
        return temp;
    }
    public static boolean checkWithinBound(int hw, int x, int y)
    {
       
	int j;
	//if(y>=x)
        {
	    j= (y-x) + hw;
	}
        /*
	else
	{
	    j = (y-x) + hw;
	}
         */
        
        //System.out.println(x + " " + y + " " + j + " " + hw + " "  + (j<0||j>=2*hw));
	if(j < 0 || j >= 2*hw)
	    return false;
	else
	    return true;
    }
    

    public int getI() {
        return i;
    }

    public void setI(int i) {
        this.i = i;
    }

    public int getJ() {
        return j;
    }

    public void setJ(int j) {
        this.j = j;
    }
}
