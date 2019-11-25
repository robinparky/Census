/**
 * @file FileFormatUnknownException.java
 * This is the source file for edu.scripps.pms.util.spectrum.FileFormatUnknowException
 * @author Tao Xu
 * @date $Date: 2014/09/09 19:29:52 $
 */



package edu.scripps.pms.census.util.io;

import java.io.IOException;

public class FileFormatUnknownException extends IOException  {
    
    private static String messageStr = "Unknow file format"; 

    public FileFormatUnknownException () {

        super(messageStr);
    }
}
