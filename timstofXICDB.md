
## How to use TimstofXICDB class:

```
TimsTOFIndex index = new TimsTOFIndex(ms2path);
```
ms2path is path timstof ms2. It will either read index file or create new one if it does not exist
```
int precursorID = index.getPrecursorID(ms2_scanNumber);
```
put in ms2 scanNumber to get precursor ID. Precursor id will be used to access information from slqite


```
TimsTOFXICDB timsTOFXICDB = new TimsTOFXICDB(path);
```
Create TimsTOFXICDB, with path to sqlite file


```
TimstofQueryResult result = timsTOFXICDB.queryAndSumPrecursorID(precursorID)
List<Pair< Double,Double>> list = result.getSummedList();
double rettime = result.retTime
```
Use precursor Id from index file to query sqllite db. It will return TimstofQueryResult object. The object contains a list of summed peaks and retention time for that precursor 

```
        TDoubleArrayList xarrayList = new TDoubleArrayList();
        TDoubleArrayList yarrayList = new TDoubleArrayList();
        double max = Double.MIN_VALUE;
        double sum = 0;
        for(Pair< Double,Double> r: resutlt)
        {
            xarrayList.add(r.getLeft());
            yarrayList.add(r.getRight());
            if(r.getRight()> max)
                max = r.getRight();
            sum+=r.getRight();
            //System.out.println(r.getLeft()+"\t"+r.getRight());
        }
        double [] xarr = xarrayList.toNativeArray();
        double [] yarr = yarrayList.toNativeArray();
        double predictedArea = getPeakAreaEstimate(xarr, yarr);
        GaussianPeakModel model =  getGaussianPeakRangeIndex(xarr, yarr, -1, resutlt.size()-1);
        model.setMaxIntensity(max);
        double peakHeight = model.getY();
        double sigma = model.getSigma();
        double area = getGaussianPeakArea(peakHeight, sigma);

```
Use this formula to get Guassian peak are under curve