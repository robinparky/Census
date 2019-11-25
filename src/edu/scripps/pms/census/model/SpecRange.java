package edu.scripps.pms.census.model;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: Yates Lab</p>
 *
 * @author Robin Park
 * @version $Id: SpecRange.java,v 1.2 2008/09/09 22:29:42 rpark Exp $
 */
 
public class SpecRange
{
    private int min;
    private int max;

    public SpecRange(int min, int max)
    {
        this.min = min;
        this.max = max;
    }

    public int getMin()
    {
        return min;
    }

    public int getMax()
    {
        return max;
    }

    public void setMin(int min)
    {
        this.min = min;
    }

    public void setMax(int max)
    {
        this.max = max;
    }
}

