/*
 * CensusConstants.java
 *
 * Created on January 23, 2006, 2:33 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package edu.scripps.pms.census;

/**
 *
 * @author rpark
 */
public final class CensusConstants {

    /** Creates a new instance of CensusConstants */
    private CensusConstants() {
    }

    public static final String DATA_TYPE="data_type";
    public static final String ELEMENT_COMPOSITION_FILE="element_composition_file";
    public static final String ENRICHMENT="enrichment";
    public static final String RESOLUTION="resolution";
    public static final String START_RANGE="start_range";
    public static final String END_RANGE="end_range";
    public static final String ISOLATION_WINDOW="isolation_window";
    public static final String MASS_ACCURACY="mass_accuracy";
    public static final String MAX_WINDOW="max_window";
    public static final String WINDOW_MARGIN="window_margin";
    public static final String SEARCH_OUTPUT="DTASelect-filter.txt";
    public static final String ISOBARIC_PURITY="isobaric_purity_filter";
    public static final String ISOBARIC_SIGNAL_TO_NOISE="isobaric_signal_to_noise";
    public static final int LABELFREE_MS1_SPLIT_SCAN_NUM=100;


    public static final int BASIC_PEPTIDE_COLUMN_SIZE = 5;
    public static final String [] PEPTIDE_COLUMNS = {"Unique", "File Name", "Scan Num", "Sequence", "XCorr", "DeltCN", "Charge", };
//    public static final String [] PEPTIDE_COLUMNS = {"U", "File Name", "Scan Num", "Sequence", "ratio", "r", "r^2", "XCorr", "DeltCN", "Charge", };
    public static final String [] PROTEIN_COLUMNS = {"Locus", "Sequence Count", "Spectrum Count", "Sequence Coverage", "Length", "MolWt", "pI", "Description", };
    public static final String [] PEPTIDE_LABEL_FREE_COLUMNS = {"Unique", "File Name", "Scan Num", "Sequence", "XCorr", "DeltCN", "Charge", };

    public static final String [] PROTEIN_SIMPLE_COLUMNS = {"Locus", };
    public static final String [] PROTEIN_MRM_FRAG_COLUMNS = {"Protein/Peptide", "Frag Ion after Precursor", "Max m/z after Precursor", "max frag ion", "max m/z", "charge state", "Precursor", "Description", };

    public static final String [] NONLABEL_COLUMNS = {"Sample", "File name", "Path", "Intensity Area",};
    public static final String [] NONLABEL_SUMMARY_COLUMNS = {"Sample", "Intensity Sum", "Average", "StDev",};


    public final static double PROTON_MASS = 1.00728;


    //Experiment Types
    public static final int MSMS_SPECIFIC_SINGLE_MASS = 13; //iTRAQ
    public static final int MSMS_SPECIFIC_MULTIPLE_MASS = 14; //iTRAQ
    public static final int MSMS_DATA_INDEPENDENT = 15; //Data Independent
    public static final int MSMS_DATA_INDEPENDENT_LFREE = 17; //Data Independent
    public static final int MSMS_QUANT = 16; //For isotopologue
    public static final int MRM_EXPERIMENT = 20;

    public static final String MS1_FILE = "ms1";
    public static final String MS2_FILE = "ms2";
    public static final String MS3_FILE = "ms3";
    public static final String MZXML = "mzXML";

    public static final String[] getBasicColumn(int scoreNum) {
	String[] arr = new String[scoreNum + 5];
	arr[0] = "Unique";
	arr[1] = "File Name";
	arr[2] = "Scan Num";
	arr[3] = "Sequence";
	arr[4] = "Charge";

	return arr;
    }

}
