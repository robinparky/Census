package scripts;

import edu.scripps.pms.census.util.io.SpectrumReader;
import edu.scripps.pms.util.spectrum.Hline;
import edu.scripps.pms.util.spectrum.Peak;
import edu.scripps.pms.util.spectrum.PeakList;
import edu.scripps.pms.util.spectrum.Zline;

import java.io.File;
import java.util.Iterator;

/**
 * Created by rpark on 2/2/17.
 */
public class MSSpectrumReader {

    //private static final double REPORT_ION=127.13108;
    //private static final double REPORT_ION=128.12812;

    //private static final double REPORT_ION=128.13444;
    //private static final double REPORT_ION=126.12773;
    //private static final double REPORT_ION=127.12476;
    private static final double TOLERANCE=0.002;



    public static void main(String[] args) throws Exception {


        // 126.127726
        double[] inputArr = {126.127726, 127.124761, 127.131081, 128.128116, 128.134436, 129.131471, 129.13779,
                130.134825, 130.141145, 131.13818, 131.144499, };


    //    double REPORT_ION=126.130;
      //  double start = REPORT_ION-TOLERANCE;
        //double end = REPORT_ION+TOLERANCE;

        System.out.print("scan\t");
        for(int i=0;i<inputArr.length;i++)
            System.out.print(inputArr[i] + "\t");
        System.out.println("");

        int temp=0;
        int reportIonCount=0;

        //File path = new File("");

       // for(File f : path.listFiles()) {

           // if(!f.getName().endsWith("_ms3.ms2"))
           //     continue;




        File f = new File("/home/rpark/test_data/pfizer/tmt/PFMS_003913_22052019_JDL_900_1.ms3");
            String line = null;
            SpectrumReader sr = new SpectrumReader(f.getAbsolutePath(), "ms3");

            Iterator<PeakList> it = sr.getSpectra();
            int counter = 0;
            int numPeaks = 0;



            while (it.hasNext()) {
                PeakList list = it.next();

            //    System.out.println(list.getHiscan());

                Peak p;
                StringBuffer sb = new StringBuffer();

                double[] intArr = new double[inputArr.length];
                for(Iterator<Peak> itr=list.getPeaks(); itr.hasNext(); )
                {
                    p = itr.next();

                    //if(start<=p.getM2z() && end>=p.getM2z())
                     //   reportIonCount++;

                    for(int i=0;i<inputArr.length;i++) {

                        double start = inputArr[i] -TOLERANCE;
                        double end = inputArr[i]+TOLERANCE;

                        if(start<=p.getM2z() && end>=p.getM2z())
                            intArr[i] += p.getIntensity();
                    }

                    //System.out.println("Peak Lines: "+p.getM2z()+"\t"+p.getIntensity()+"\t"+p.getChargeState()+"\t"+p.getResolution());

                }

                System.out.print(list.getLoscan() + "\t");
                for(double d:intArr)
                    System.out.print(d + "\t");

                System.out.println("");
            }
      //      System.out.println("range start:\t"+start);
        //    System.out.println("range end:\t"+end);
          //  System.out.println("count: "+reportIonCount);

           // System.exit(0);




//        System.out.println("remove this line");

    }
}
