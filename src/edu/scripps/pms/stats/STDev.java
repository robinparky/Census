package edu.scripps.pms.stats;

public class STDev 
{
    public static void main(String args[]) throws Exception
    {
	STDev dev = new STDev();

	//double[] arr = {4, 9, 11, 12, 17, 5, 8, 12, 14,};
	double[] arr = 
	{1.53,
	1.45,
	1.48,
	1.64,
	1.29,
	0.5,
	1.02,};

	double stdev = dev.getStdev(arr);
	double mean = dev.getMean(arr);

	for(int i=0;i<arr.length;i++)
	{
	    double diff = Math.abs(mean-arr[i]);
	    System.out.println( diff/stdev);
	}

	System.out.println("");

	int size=arr.length;
    
	for(int i=0;i<arr.length;i++)
	{
	    double diff = Math.abs(mean-arr[i]);
	    double z = diff/stdev;

	    double t = Math.sqrt( (size*(size-2)*z*z)/((size-1)*(size-1)-size*z*z) );
	    System.out.println(t);
	}
    }

    public static double getStdev(double[] arr)
    {
	double mean = getMean(arr);

	double devSum=0;
	for(int i=0;i<arr.length;i++)
	{
	    devSum += (mean-arr[i])*(mean-arr[i]);
	}

	return Math.sqrt(devSum/(arr.length-1));
    }

    //calculate stdev without considering negative values
    public static double getStdevWithoutNegative(double[] arr)
    {
	double sum =0;
	int negNum=0;
	for(int i=0;i<arr.length;i++)
	{
	    if(arr[i]<0)
	    {
		negNum++;
		continue;
	    }

	    sum+=arr[i];
	}

	//cannot caluclate stdev
	if(arr.length+1 <= negNum)
	    return -1;

	double mean = sum/(arr.length-negNum);

	double devSum=0;
	for(int i=0;i<arr.length;i++)
	{
	    if(arr[i]<0)
		continue;

	    devSum += (mean-arr[i])*(mean-arr[i]);
	}

	return Math.sqrt(devSum/(arr.length-1-negNum));
    }

    public static double getMean(double[] arr)
    {
	double sum =0;
	for(int i=0;i<arr.length;i++)
	{
	    sum+=arr[i];
	}

	double mean = sum/arr.length;

	return mean;
    }
    
    public static double getMeanExcludingNegative(double[] arr)
    {
	double sum =0;
	int size=0;
	for(int i=0;i<arr.length;i++)
	{
	    if(arr[i]<0)
		continue;

	    sum+=arr[i];
	    size++;
	}

	if(size==0)
	    return -1;

	//double mean = sum/arr.length;
	return sum/size;
    }
    
    //long
    public static double getStdev(long[] arr)
    {
	double mean = getMean(arr);

	double devSum=0;
	for(int i=0;i<arr.length;i++)
	{
	    devSum += (mean-arr[i])*(mean-arr[i]);
	}

	return Math.sqrt(devSum/(arr.length-1));
    }

    public static double getMean(long[] arr)
    {
	long sum =0;
	for(int i=0;i<arr.length;i++)
	{
	    sum+=arr[i];
	}

	double mean = (double)sum/arr.length;

	return mean;
    }
}
