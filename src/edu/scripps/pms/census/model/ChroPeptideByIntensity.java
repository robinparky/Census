/*
 * ChroPeptide.java
 *
 * Created on May 23, 2005, 11:38 AM
 */

package edu.scripps.pms.census.model;

import java.util.*;

/**
 *
 * @author rpark
 * @version $Id: ChroPeptideByIntensity.java,v 1.2 2014/08/27 18:00:35 rpark Exp $
 */
public class ChroPeptideByIntensity implements Comparator<ChroPeptide> {

    public int compare(ChroPeptide o1, ChroPeptide o2) {
	return (int)((o1.getTotalIntensity()) - (o2.getTotalIntensity()));        
    }
}
