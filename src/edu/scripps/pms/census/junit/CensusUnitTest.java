package edu.scripps.pms.census.junit;

import junit.framework.*;

public class CensusUnitTest 
    extends TestCase {

    public CensusUnitTest(String s)
    {
        super(s);
    }

    public static Test suite()
    {
        TestSuite suite = new TestSuite();

//        suite.addTest( new IndexTestCase("indexValidate") );
        suite.addTest( new ParamTestCase("paramTest") );


//suite.addTest( new TestQuery("testQuery") );

        return suite;
    }
}
