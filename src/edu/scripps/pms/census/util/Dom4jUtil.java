/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package edu.scripps.pms.census.util;


import org.dom4j.Document;
import org.dom4j.DocumentException;
import org.dom4j.io.SAXReader;
import org.dom4j.Element;

import java.io.File;


/**
 *
 * @author rpark
 */
public class Dom4jUtil {           
    public static Document getDocument(File file) throws DocumentException {
        SAXReader reader = new SAXReader();
        Document document = reader.read(file);
        return document;
    }
    
    public static Element getRootEle(File file) throws DocumentException {
        
        Document document = getDocument(file);
        return document.getRootElement();
    }
}