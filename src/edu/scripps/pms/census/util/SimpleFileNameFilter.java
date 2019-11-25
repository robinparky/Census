/*
 * SimpleFileNameFilter.java
 *
 * Created on March 25, 2005, 11:10 AM
 */

package edu.scripps.pms.census.util;

/**
 *
 * @author rpark
 */


import java.io.File;
import javax.swing.filechooser.FileFilter;

public class SimpleFileNameFilter extends FileFilter
{

    private String extensions[];
    private String description;
    
    public SimpleFileNameFilter(String ext, String descr)
    {
        this(new String[] {
            ext
        }, descr);
    }

    public SimpleFileNameFilter(String exts[], String descr)
    {
        extensions = new String[exts.length];
        for(int i = exts.length - 1; i >= 0; i--)
            extensions[i] = exts[i].toLowerCase();

        description = descr != null ? descr : (new StringBuilder()).append(exts[0]).append(" files").toString();
    }

    public boolean accept(File f)
    {
        if(f.isDirectory())
            return true;
        String name = f.getName().toLowerCase();
        for(int i = extensions.length - 1; i >= 0; i--)
            if(name.endsWith(extensions[i]))
                return true;

        return false;
    }

    public String getDescription()
    {
        return description;
    }
}