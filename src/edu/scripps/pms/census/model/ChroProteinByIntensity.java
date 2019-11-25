/*
 * ChroProtein.java
 *
 * Created on May 23, 2005, 11:38 AM
 */

package edu.scripps.pms.census.model;

import java.util.*;

/**
 *
 * @author rpark
 * @version $Id: ChroProteinByIntensity.java,v 1.2 2014/08/27 18:00:35 rpark Exp $
 */
public class ChroProteinByIntensity implements Comparator<ChroProtein> {

    public int compare(ChroProtein o1, ChroProtein o2) {
	return o1.getIntensity().compareTo(o2.getIntensity());
    }
}
