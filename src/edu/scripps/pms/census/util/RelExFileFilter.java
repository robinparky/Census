/*
 * RelaxFileFilter.java
 *
 * Created on March 21, 2005, 12:23 PM
 */

package edu.scripps.pms.census.util;

import java.io.FilenameFilter;
import java.io.File;

/**
 *
 * @author rpark
 */
public class RelExFileFilter implements FilenameFilter {
    
    /** Creates a new instance of RelaxFileFilter */
    private String format;
    
    public RelExFileFilter(String format) 
    {
        super();
        this.format = format;
    }
    
    public boolean accept(File dir, String fileName)
    {
        return fileName.endsWith(format);
    }
    
}
