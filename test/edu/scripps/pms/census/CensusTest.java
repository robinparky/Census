package edu.scripps.pms.census;

import edu.scripps.pms.census.test.CensusTestConfig;
import edu.scripps.pms.census.util.io.FileUtil;
import org.apache.log4j.Logger;
import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;

public class CensusTest {

    private static final Logger LOGGER = Logger.getLogger(CensusTest.class);


    private String inputPath = CensusTestConfig.getInstance().getDataPath() + File.separator + "census_test";

    private String outputPath = CensusTestConfig.getInstance().getDataPath() + File.separator + "census_test_output";

    @BeforeClass
    public void beforeClass() throws IOException {
        LOGGER.info("inputPath ::" + inputPath);
        LOGGER.info("outputPath ::" + outputPath);
        LOGGER.info("working dir ::" + System.getProperty("user.dir"));


        LOGGER.info("Downloading  Data");

        //dataDownloader.downloadData("peptide_compare.zip");
        FileUtil.makeDir(outputPath);

        LOGGER.info("Done");

    }


    @AfterClass
    public void AfterTesting() {

        LOGGER.info("Remove ");


        LOGGER.info("Done ");
    }


    @Test
    public void TestTMTChro()
    {
        String path = inputPath+File.separator+"tmt_chro" + File.separator;
        String line = "-c "+path + "census_config.xml -i "+path +"DTASelect-filter.txt -f "+path;
        String [] args = line.split(" ");
        try {
            Census.main(args);
        } catch (Exception e) {
            e.printStackTrace();
            Assert.assertEquals(e.getMessage(), true);
            LOGGER.error("Error " + e.getMessage());
        }
    }




}


