package scripts;

import edu.scripps.pms.census.util.io.SpectrumReader;
import edu.scripps.pms.util.spectrum.Hline;
import edu.scripps.pms.util.spectrum.Peak;
import edu.scripps.pms.util.spectrum.PeakList;
import edu.scripps.pms.util.spectrum.Zline;

import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.Iterator;

/**
 * Created by rpark on 3/3/17.
 */
public class WatersDDA {

    public static void main(String[] args) throws Exception {

        double intensityThreahold = 20.0;
        String path = "/data/2/rpark/ip2_data/rpark/merck_michael/BSA_DDA_2nd_sample_2017_03_02_08_227244/spectra/old/michael_20170215_mel_BSA_DDA_101.ms2_orig";
        String ext = "ms2";
        //SpectrumReader sr = new SpectrumReader(args[0], args[1]);

        PrintStream out = new PrintStream(new FileOutputStream(path +  ".newms2"));

        SpectrumReader sr = new SpectrumReader(path, ext);
        Hline h = new Hline(sr.getHlines());

        for(Iterator<String> itr = sr.getHlines(); itr.hasNext();  ) {
            out.println(itr.next());
        }


        Iterator<PeakList> it = sr.getSpectra();
        int counter = 0;
        int numPeaks = 0;
        //boolean sortByIntensity = true;

        while (it.hasNext()) {
            PeakList list = it.next();
//System.out.println("RetentionTime: " + list.getRetentionTime());
//          System.out.println("==" + list.getLoscan());
            //\t" + list.getHiscan() + "\t");

            Peak p;
            StringBuffer sb = new StringBuffer();
            sb.append("S\t");
            sb.append(list.getHiscan());
            sb.append("\t");
            sb.append(list.getHiscan());
            sb.append("\t");
            sb.append(list.getPrecursorMass());
            sb.append("\n");

            for(String each:list.getIlines()) {
                sb.append(each);
                sb.append("\n");
            }

            for(Iterator<Zline> itr = list.getZlines(); itr.hasNext();  ) {
                Zline z = itr.next();
                sb.append("Z\t" + z.getChargeState() + "\t" + z.getM2z()).append("\n");

            }


            StringBuffer peakSb = new StringBuffer();
            int peakCount=0;
            for(Iterator<Peak> itr=list.getPeaks(); itr.hasNext(); )
            {
                p = itr.next();

                if(p.getIntensity()>intensityThreahold) {
                    peakCount++;
                    peakSb.append(p.getM2z());
                    peakSb.append("\t");
                    peakSb.append((int)p.getIntensity());
                    peakSb.append("\n");
                }
                //System.out.println("Peak Lines: "+p.getM2z()+"\t"+p.getIntensity()+"\t"+p.getChargeState()+"\t"+p.getResolution());
            }

            if(peakCount>10) {
                out.print(sb.toString());
                out.print(peakSb.toString());
            }
        }
       // System.out.println("Total number of spectra processed: " + counter);
       // System.out.println("Average number of peaks per list: " + numPeaks/counter);
        sr.closeDataFile();
        out.close();

        System.out.println("Finished");
    }
}
