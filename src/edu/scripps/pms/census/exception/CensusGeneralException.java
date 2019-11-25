/*
 * InvalidAAException.java
 *
 * Created on September 20, 2006, 1:59 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package edu.scripps.pms.census.exception;

/**
 *
 * @author rpark
 */
public class CensusGeneralException extends Exception {
    
    /** Creates a new instance of InvalidAAException */
    public CensusGeneralException() {
    }
    
    public CensusGeneralException(String message) {
	super(message);
    }
    
}
