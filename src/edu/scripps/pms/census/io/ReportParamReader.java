
/*
 * ChroReader.java
 *
 * Created on May 17, 2005, 12:00 PM
 */

package edu.scripps.pms.census.io;

import edu.scripps.pms.census.*;
import edu.scripps.pms.census.conf.Configuration;
import java.io.*;
import java.util.*;

/**
 *
 * @author  Robin Park
 * @version $Id: ReportParamReader.java,v 1.7 2014/02/21 22:10:52 rpark Exp $
 */

public class ReportParamReader {

    public static void main(String[] args) throws Exception {
	readParam(args[0]);
    }

    public static ReportParam readParam(String fileName) throws IOException
    {
	BufferedReader br = new BufferedReader(new FileReader(fileName));

	ParamTable ht = new ParamTable();

	String lastLine = null;

	while( (lastLine=br.readLine())!=null )
	{
	    ht.addParamLine(lastLine);
	}

        ReportParam rParam = new ReportParam();
	Object tmpValue = ht.get("discard_singleton_peptide");
	if(tmpValue!=null)
	{
		String tmp = tmpValue.toString();
		if("1".equals(tmp))
			rParam.setDiscardAN(true);
		else if("0".equals(tmp))
			rParam.setDiscardAN(false);
	}

	tmpValue = ht.get("nofilter");
	if(tmpValue!=null)
	{
		String tmp = tmpValue.toString();
		if("1".equals(tmp))
			rParam.setNoFilter(true);
		else if("0".equals(tmp))
			rParam.setNoFilter(false);
	}

	tmpValue = ht.get("use_outlier");
	if(tmpValue!=null)
	{
		String tmp = tmpValue.toString();
		if("1".equals(tmp))
			rParam.setPValueSelect(true);
		else if("0".equals(tmp))
			rParam.setPValueSelect(false);
	}


	Object tmpObject = ht.get("iterate_outlier");
	if(null == tmpObject) {
	    rParam.setIterateOutlier(true);
	} else {
		tmpValue = tmpObject.toString();
		if("1".equals(tmpValue))
			rParam.setIterateOutlier(true);
		else if("0".equals(tmpValue))
			rParam.setIterateOutlier(false);
	}


	tmpValue = ht.get("outlier_p_value");
	if(tmpValue!=null)
	    rParam.setPValue( Double.parseDouble(tmpValue.toString()) );

	tmpValue = ht.get("singleton_arearatio_lower_threshold");
	if(tmpValue!=null)    rParam.setAllNoneLowerBound( Double.parseDouble(tmpValue.toString()) );

	tmpValue = ht.get("singleton_arearatio_upper_threshold");
	if(tmpValue!=null)       rParam.setAllNoneUpperBound( Double.parseDouble(tmpValue.toString()) );

	tmpValue = ht.get("determinant_facter_threshold");
	if(tmpValue!=null)       rParam.setDetValue( Double.parseDouble(tmpValue.toString()) );

  	Object profileObj = ht.get("profile_score_threshold");
	if(null != profileObj) {
		String tmpsValue = profileObj.toString();
		rParam.setProfileScore(Double.parseDouble(tmpsValue));
	}

    Object minimumPepObj = ht.get("minimum_peptide_per_protein");
    if(null != minimumPepObj) {
      String tmpsValue = minimumPepObj.toString().trim();
      double tempD = Double.parseDouble(tmpsValue);
      int tempI = (int) tempD;
      rParam.setMinimumPeptidePerProtein(tempI);
    }
    Object minimumUniquePepObj = ht.get("minimum_unique_peptide_per_protein");
    if(null!=minimumUniquePepObj)
	{
		String tempValue = minimumUniquePepObj.toString().trim();
		double tempD = Double.parseDouble(tempValue);
		int tempI = (int) tempD;
		rParam.setMinimumNumberOfUniquePeptides(tempI);

	}

    String tmtPurity = (String)ht.get(CensusConstants.ISOBARIC_PURITY);
    if(tmtPurity!=null)
	{
		double purityDouble = Double.parseDouble(tmtPurity);
		rParam.setPurityThreshold(purityDouble);
	} else {
      rParam.setPurityThreshold(-1);

    }

	String signal = (String)ht.get(CensusConstants.ISOBARIC_SIGNAL_TO_NOISE);
	if(signal!=null)
	{
		rParam.setSignalToNoiseThreshold(Double.parseDouble(signal));
	} else {
    rParam.setSignalToNoiseThreshold(-1.0);
  }

	tmpValue = ht.get("unique_peptide_only");
	if(tmpValue!=null)
	{
		if("true".equals(tmpValue.toString()) || "1".equals(tmpValue.toString()) )
			rParam.setIsUniquePeptide(true);
		else if("false".equals(tmpValue))
			rParam.setIsUniquePeptide(false);
	}


	tmpValue = ht.get("remove_negative_determinant_factor");
	if(tmpValue!=null)
	{
		if("1".equals(tmpValue.toString()))
			rParam.setRemoveNegative(true);
		else if("0".equals(tmpValue))
			rParam.setRemoveNegative(false);
	}


	tmpValue = ht.get("discard_unlabeled_peptide");
	if(tmpValue!=null)
	{
		if("1".equals(tmpValue.toString()))
			rParam.setDiscardUnlabeledPeptide(true);
		else if("0".equals(tmpValue))
			rParam.setDiscardUnlabeledPeptide(false);
	}


	Object dRevProtein =  ht.get("discard_reverse_protein");

	if(null != dRevProtein) {
	    tmpValue = dRevProtein.toString();

	    if("1".equals(tmpValue))
		rParam.setDiscardReverseProtein(true);
	    else if("0".equals(tmpValue))
		rParam.setDiscardReverseProtein(false);
	}


	tmpValue = ht.get("correction_factor");
	 if(tmpValue!=null)   rParam.setCorrectFactorValue( Double.parseDouble(tmpValue.toString()) );

	tmpValue = ht.get("singleton_minimum_peptides");
		if(tmpValue!=null)    rParam.setAllNoneMinPeptide( Integer.parseInt(tmpValue.toString()) );

	tmpValue = ht.get("singleton_score_threshold");
		if(tmpValue!=null)    rParam.setAllNoneCompositeScore( Double.parseDouble(tmpValue.toString()) );

	tmpValue = ht.get("filter_fragment_ions");
	if(tmpValue!=null)
	{
		String tmpS = tmpValue.toString();
		if("1".equals(tmpS))
			rParam.setFilterFragmentIons(true);
		else if("0".equals(tmpS))
			rParam.setFilterFragmentIons(false);
	}


        Object smoothPeaks = ht.get("smooth_peaks");
        if(null != smoothPeaks) {
	tmpValue = smoothPeaks.toString();
            if("1".equals(tmpValue))
                rParam.setSmoothingPeaks(true);
            else if("0".equals(tmpValue))
                rParam.setSmoothingPeaks(false);
        }

	if(null == ht.get("profile_score")) {
		rParam.setUseProfileScore(false);
	} else {
		String tmpsValue = ht.get("profile_score").toString();
		rParam.setProfileScore( Double.parseDouble(tmpsValue) );
	}

	if(null != ht.get("max_spectrum_shift")) {
            String tmpsValue = ht.get("max_spectrum_shift").toString();
            rParam.setMaxSpectrumShift( Integer.parseInt(tmpsValue));
            Configuration conf = Configuration.getInstance();
            conf.setMaxSpectrumShift(rParam.getMaxSpectrumShift());

	}


        if(null != ht.get("intensity_threshold")) {
            String tmpsValue = ht.get("intensity_threshold").toString();

            rParam.setIntensityThreshold( Double.parseDouble(tmpsValue.substring(1)) );
	}

	return rParam;
    }

    static class ParamTable extends Hashtable {
	public void addParamLine(String paramLine) {
	    if(null == paramLine || paramLine.startsWith("H") || "".equals(paramLine.trim()) || paramLine.startsWith("#"))
		return;

	    String[] arr = paramLine.split("\t");
	    int index = arr[1].indexOf("#");

	    if(index<=0)
		put(arr[0].trim(), arr[1].trim());
	    else
		put(arr[0].trim(), arr[1].substring(0, index).trim());

	}


    }
}
