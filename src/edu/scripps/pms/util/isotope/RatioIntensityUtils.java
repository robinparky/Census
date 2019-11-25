package edu.scripps.pms.util.isotope;

/**
 * Created by Titus Jung titusj@scripps.edu on 4/25/18.
 */
public class RatioIntensityUtils {

    public static double checkRatio(double ratio)
    {
        if (ratio> 100 || Double.isInfinite(ratio)) {
            return 100;
        }
        else if(ratio <=0.01)
        {
            return 0.01;
        }
        return ratio;
    }

}
