/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.scripps.pms.census.labelFree;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author Rohan
 */
public class LabelfreePeptideParser {
    private List<Integer> sequenceindex = new ArrayList<>();
    private List<Integer> chargestateindex = new ArrayList<>();
    private List<Integer> intensityindex = new ArrayList<>();
    private List<Integer> filenameindex = new ArrayList<>();
    private List<Integer> xCorrindex = new ArrayList<>();
    private List<Integer> dcnindex = new ArrayList<>();
    private List<Integer> dmassindex = new ArrayList<>();
    private List<Integer> sprankindex = new ArrayList<>();
    private List<Integer> spscoreindex = new ArrayList<>();
    private List<Integer> redundancyindex = new ArrayList<>();
    private List<Integer> startrangeindex = new ArrayList<>();
    private List<Integer> endrangeindex = new ArrayList<>();
    private List<Integer> retentiontimeindex = new ArrayList<>();
    private List<Integer> ioninjectionindex = new ArrayList<>();
    private List<Integer> proteinindex = new ArrayList<>();
    private List<Integer> descriptionindex = new ArrayList<>();
    private List<LabelfreePeptide> pepList = new ArrayList<>();

    public static void main(String[] args) {
        String inputpath = "C:/Users/Rohan/Documents/garbage_test/census_labelfree_out_11390_pep.txt";
        LabelfreePeptideParser l = new LabelfreePeptideParser();
        l.read(inputpath);
    }


    public List<LabelfreePeptide> read(String inputpath) {
        File file = new File(inputpath);
        List group1 = new ArrayList();
        List group2 = new ArrayList();
        int line = 0;
        int counter = 0;
        List<LabelfreePeptide> pepList = new ArrayList<>();
        try {
            FileReader fread = new FileReader(file);
            BufferedReader br = new BufferedReader(fread);
            String eachLine = null;
            while ((eachLine = br.readLine()) != null) {
                if (eachLine.startsWith("H\tGROUP_SAMPLE")) {
                    line++;
                    String words[] = eachLine.split("\t");
                    for (int i = 0; i < words.length; i++) {
                        if (i > 2 && line == 1) {
                            group1.add(counter);
                            counter++;
                        } else if (i > 2 && line == 2) {
                            group2.add(counter);
                            counter++;
                        }
                    }
                }
                if (eachLine.startsWith("SLINE\t")) {
                    readHeaderLine(eachLine);
                }
                if (eachLine.startsWith("S\t")) {
                    String[] words = eachLine.split("\t");
                    LabelfreePeptide peptide = new LabelfreePeptide();
                    for (int i = 0; i < words.length; i++) {
                        if (sequenceindex.contains(i)) {
                            if (!(words[i].contains("NA"))) {
                                peptide.setSequence(words[i]);
                            }
                        } else if (chargestateindex.contains(i)) {
                            if (!(words[i].contains("NA"))) {
                                peptide.setChargestate(Integer.parseInt(words[i]));
                            }
                        } else if (intensityindex.contains(i)) {

                            if (group1.contains(intensityindex.indexOf(i))) {
                                peptide.addIntensity1(words[i]);
                            } else if (group2.contains(intensityindex.indexOf(i))) {
                                peptide.addIntensity2(words[i]);
                            }

                        } else if (filenameindex.contains(i)) {
                            if (group1.contains(filenameindex.indexOf(i))) {
                                peptide.addFilename1(words[i]);
                            } else if (group2.contains(filenameindex.indexOf(i))) {
                                peptide.addFilename2(words[i]);
                            }
                        } else if (xCorrindex.contains(i)) {
                            if (group1.contains(xCorrindex.indexOf(i))) {
                                peptide.addXCorr1(words[i]);
                            } else if (group2.contains(xCorrindex.indexOf(i))) {
                                peptide.addXCorr2(words[i]);
                            }
                        } else if (descriptionindex.contains(i)) {
                            peptide.setDescription(words[i]);
                        } else if (proteinindex.contains(i)) {
                            peptide.setProtein(words[i]);
                        } else if (retentiontimeindex.contains(i)) {
                            if (group1.contains(retentiontimeindex.indexOf(i))) {
                                peptide.addrt1(words[i]);
                            } else if (group2.contains(retentiontimeindex.indexOf(i))) {
                                peptide.addrt2(words[i]);
                            }
                        } else if (ioninjectionindex.contains(i)) {
                            if (group1.contains(ioninjectionindex.indexOf(i))) {
                                peptide.addion1(words[i]);
                            } else if (group2.contains(ioninjectionindex.indexOf(i))) {
                                peptide.addion2(words[i]);
                            }
                        } else if (dcnindex.contains(i)) {
                            if (group1.contains(dcnindex.indexOf(i))) {
                                peptide.adddcn1(words[i]);
                            } else if (group2.contains(dcnindex.indexOf(i))) {
                                peptide.adddcn2(words[i]);
                            }
                        } else if (dmassindex.contains(i)) {
                            if (group1.contains(dmassindex.indexOf(i))) {
                                peptide.adddmass1(words[i]);
                            } else if (group2.contains(dmassindex.indexOf(i))) {
                                peptide.adddmass2(words[i]);
                            }
                        } else if (sprankindex.contains(i)) {
                            if (group1.contains(sprankindex.indexOf(i))) {
                                peptide.addsprank1(words[i]);
                            } else if (group2.contains(sprankindex.indexOf(i))) {
                                peptide.addsprank2(words[i]);
                            }
                        } else if (startrangeindex.contains(i)) {
                            if (group1.contains(startrangeindex.indexOf(i))) {
                                peptide.addstart1(words[i]);
                            } else if (group2.contains(startrangeindex.indexOf(i))) {
                                peptide.addstart2(words[i]);
                            }
                        } else if (endrangeindex.contains(i)) {
                            if (group1.contains(endrangeindex.indexOf(i))) {
                                peptide.addend1(words[i]);
                            } else if (group2.contains(endrangeindex.indexOf(i))) {
                                peptide.addend2(words[i]);
                            }
                        } else if (redundancyindex.contains(i)) {
                            if (group1.contains(redundancyindex.indexOf(i))) {
                                peptide.addred1(words[i]);
                            } else if (group2.contains(redundancyindex.indexOf(i))) {
                                peptide.addred2(words[i]);
                            }
                        } else if (spscoreindex.contains(i)) {
                            if (group1.contains(spscoreindex.indexOf(i))) {
                                peptide.addsp1(words[i]);
                            } else if (group2.contains(spscoreindex.indexOf(i))) {
                                peptide.addsp2(words[i]);
                            }
                        }
                    }

                    //        AnovaUtil.calculateAnovaPvalue(peptide.getIntensity1(),peptide.getIntensity2());
                    pepList.add(peptide);
                }
            }
            return pepList;

        } catch (Exception e) {
            e.printStackTrace();
        }
        return null;
    }

    private void readHeaderLine(String eachLine) {
        String[] words = eachLine.split("\t");
        for (int i = 1; i < words.length; i++) {
            if (words[i].startsWith("SEQUENCE")) {
                sequenceindex.add(i);
            } else if (words[i].startsWith("FILENAME")) {
                filenameindex.add(i);
            } else if (words[i].startsWith("CSTATE")) {
                chargestateindex.add(i);
            } else if (words[i].startsWith("INTENSITY")) {
                intensityindex.add(i);
            } else if (words[i].startsWith("XCORR")) {
                xCorrindex.add(i);
            } else if (words[i].startsWith("DCN")) {
                dcnindex.add(i);
            } else if (words[i].startsWith("DMASS")) {
                dmassindex.add(i);
            } else if (words[i].startsWith("SPRANK")) {
                sprankindex.add(i);
            } else if (words[i].startsWith("SPSCORE")) {
                spscoreindex.add(i);
            } else if (words[i].startsWith("REDUNDANCY")) {
                redundancyindex.add(i);
            } else if (words[i].startsWith("STARTRANGE")) {
                startrangeindex.add(i);
            } else if (words[i].startsWith("ENDRANGE")) {
                endrangeindex.add(i);
            } else if (words[i].startsWith("RETENTIONTIME")) {
                retentiontimeindex.add(i);
            } else if (words[i].startsWith("IONINJECTIONTIME")) {
                ioninjectionindex.add(i);
            } else if (words[i].contentEquals("PROTEIN")) {
                proteinindex.add(i);
            } else if (words[i].contains("PROTEIN DESCRIPTION")) {
                descriptionindex.add(i);
            }
        }
    }
}
