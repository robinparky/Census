package edu.scripps.pms.census.util.io;

import edu.scripps.pms.census.util.CSVHeaderIndex;
import org.jdom.Element;
import org.jfree.data.io.CSV;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class LuciphorReader {

    private CSVHeaderIndex specIDIndex = new CSVHeaderIndex("specId");
    private CSVHeaderIndex peptideSeqIndex = new CSVHeaderIndex("peptide");
    private CSVHeaderIndex predictedSeq1Index = new CSVHeaderIndex("predictedPep1");
    private CSVHeaderIndex predictedSeq2Index = new CSVHeaderIndex("predictedPep2");
    private CSVHeaderIndex numPPSIndex = new CSVHeaderIndex("numPPS");
    private CSVHeaderIndex numRPSIndex = new CSVHeaderIndex("numRPS");
    private CSVHeaderIndex psmScoreIndex = new CSVHeaderIndex("psmScore");
    private CSVHeaderIndex deltaScoreIndex = new CSVHeaderIndex("deltaScore");
    private CSVHeaderIndex pep1ScoreIndex = new CSVHeaderIndex("pep1score");
    private CSVHeaderIndex pep2ScoreIndex = new CSVHeaderIndex("pep2score");
    private CSVHeaderIndex globalFLRScoreIndex = new CSVHeaderIndex("globalFLR");
    private CSVHeaderIndex localFLRScoreIndex = new CSVHeaderIndex("localFLR");
    private CSVHeaderIndex [] headArr = {specIDIndex, peptideSeqIndex, predictedSeq1Index, predictedSeq2Index,
            numPPSIndex, numRPSIndex, psmScoreIndex, deltaScoreIndex, pep1ScoreIndex,pep2ScoreIndex,
            globalFLRScoreIndex, localFLRScoreIndex };


    public  class LuciphorPeptide
    {
        public final String specID;
        public final String peptide;
        public final String predictedPep1;
        public final String predictedPep2;
        public final String numPPS;
        public final String numRPS;
        public final String psmScore;
        public final String deltaScore;
        public final String pep1score;
        public final String pep2score;
        public final String globalFLR;
        public final String localFLR;

        public LuciphorPeptide(String [] arr)
        {
            specID = arr[specIDIndex.getIndex()];
            peptide = arr[peptideSeqIndex.getIndex()];
            predictedPep1 = arr[predictedSeq1Index.getIndex()];
            predictedPep2 = arr[predictedSeq2Index.getIndex()];
            numPPS = arr[numPPSIndex.getIndex()];
            numRPS = arr[numRPSIndex.getIndex()];
            psmScore = arr[psmScoreIndex.getIndex()];
            deltaScore = arr[deltaScoreIndex.getIndex()];
            pep1score = arr[pep1ScoreIndex.getIndex()];
            pep2score = arr[pep2ScoreIndex.getIndex()];
            globalFLR = arr[globalFLRScoreIndex.getIndex()];
            localFLR = arr[localFLRScoreIndex.getIndex()];
        }

        public LuciphorPeptide(String specID, String peptideSequence, org.dom4j.Element luciphorElement) {
            this.specID = specID;
            this.peptide = peptideSequence;
            predictedPep1= luciphorElement.attributeValue(predictedSeq1Index.headerName );
            predictedPep2= luciphorElement.attributeValue(predictedSeq2Index.headerName);
            numPPS = luciphorElement.attributeValue(numPPSIndex.headerName);
            numRPS = luciphorElement.attributeValue(numRPSIndex.headerName);
            psmScore= luciphorElement.attributeValue(psmScoreIndex.headerName);
            deltaScore= luciphorElement.attributeValue(deltaScoreIndex.headerName);
            pep1score= luciphorElement.attributeValue(pep1ScoreIndex.headerName);
            pep2score= luciphorElement.attributeValue(pep2ScoreIndex.headerName);
            globalFLR = luciphorElement.attributeValue(globalFLRScoreIndex.headerName);
            localFLR= luciphorElement.attributeValue(localFLRScoreIndex.headerName);
        }


        public Element generateElement()
        {
            Element luciphorElement = new Element("luciphor");
            luciphorElement.setAttribute(predictedSeq1Index.headerName, predictedPep1);
            luciphorElement.setAttribute(predictedSeq2Index.headerName, predictedPep2);
            luciphorElement.setAttribute(numPPSIndex.headerName, numPPS);
            luciphorElement.setAttribute(numRPSIndex.headerName, numRPS);
            luciphorElement.setAttribute(psmScoreIndex.headerName, psmScore);
            luciphorElement.setAttribute(deltaScoreIndex.headerName, deltaScore);
            luciphorElement.setAttribute(pep1ScoreIndex.headerName, pep1score);
            luciphorElement.setAttribute(pep2ScoreIndex.headerName, pep2score);
            luciphorElement.setAttribute(globalFLRScoreIndex.headerName, globalFLR);
            luciphorElement.setAttribute(localFLRScoreIndex.headerName, localFLR);
            return luciphorElement;
        }

    }

    private void setIndices(String [] arr)
    {
        resetIndices();
        for(int i=0; i<arr.length; i++)
        {
            String head = arr[i];
            for(CSVHeaderIndex index: headArr)
            {
                index.checkSetIndex(head, i);

            }
            /*
            specIDIndex.checkSetIndex(head, i);
            peptideSeqIndex.checkSetIndex(head, i);
            predictedSeq1Index.checkSetIndex(head, i);
            predictedSeq2Index.checkSetIndex(head, i);
            numPPSIndex.checkSetIndex(head, i);
            numRPSIndex.checkSetIndex(head, i);
            deltaScoreIndex.checkSetIndex(head, i);
            pep1ScoreIndex.checkSetIndex(head, i);
            globalFLRScoreIndex.checkSetIndex(head, i);
            pep2ScoreIndex.checkSetIndex(head, i);
            localFLRScoreIndex.checkSetIndex(head, i);*/
        }
    }

    private void resetIndices()
    {
        for(CSVHeaderIndex index: headArr)
        {
            index.reset();
        }
    }

    public Map<String,LuciphorPeptide> readLuciphorGetMap(String path) throws IOException {
        Map<String, LuciphorPeptide> result = new HashMap<>();
        List<LuciphorPeptide> peptideList = readLuciphor(path);
        for(LuciphorPeptide peptide: peptideList)
        {
            result.put(peptide.specID, peptide);
        }
        return result;
    }

    public  List<LuciphorPeptide> readLuciphor(String path) throws IOException {
        String line;
        BufferedReader br = new BufferedReader(new FileReader(path));
        List<LuciphorPeptide> luciphorPeptideList = new ArrayList<>();
        while((line = br.readLine()) != null)
        {
            String [] arr = line.split("\t");
            if(arr[0].equals("specId"))
            {
                setIndices(arr);
            }
            else
            {
                luciphorPeptideList.add(new LuciphorPeptide(arr));
            }
        }
        br.close();
        return luciphorPeptideList;
    }

    public static void main(String [] args) throws IOException {
        LuciphorReader reader = new LuciphorReader();
        List<LuciphorPeptide> peptides =  reader.readLuciphor(args[0]);

        for(LuciphorPeptide peptide: peptides)
        {
            System.out.println(peptide.specID);
        }


    }

    public LuciphorPeptide readLuciphorElement(String id, String sequence, org.dom4j.Element luciphorElement)
    {
        return new LuciphorPeptide(id, sequence, luciphorElement);
    }

}
