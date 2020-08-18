package edu.scripps.pms.census.test;

import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;

public class CensusTestConfig {

    private static CensusTestConfig config = null;
    private String dataPath = null;
    private boolean downloadMode = false;
    private String url = null;
    private Properties properties = null;

    private CensusTestConfig()
    {}

    public static CensusTestConfig getInstance()
    {
        if(config ==null)
        {
            config = new CensusTestConfig();
            config.init();
        }
        return config;
    }

    public String getValue(String param) {
        return this.properties.getProperty(param);
    }


    private void init()
    {
        properties = new Properties();
        try (InputStream input = this.getClass().getResourceAsStream("/censusTest.properties")) {
            this.properties = new Properties();
            this.properties.load(input);
            this.url = this.properties.getProperty("test_data");
            this.dataPath = this.properties.getProperty("data_path");
        } catch (IOException ex) {
            ex.printStackTrace();
        }

    }


    public String getDataPath() {
        return dataPath;
    }

    public boolean isDownloadMode() {
        return downloadMode;
    }

    public String getUrl() {
        return url;
    }

    public Properties getProperties() {
        return properties;
    }
}
