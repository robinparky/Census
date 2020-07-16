
## How to use TimstofXICDB class:

```
TimsTOFIndex index = new TimsTOFIndex(ms2path);
```
ms2path is path timstof ms2. It will either read index file or create new one if it does not exist
```
int parentID = index.getParentId(ms2_scanNumber);
```
put in ms2 scanNumber to get parent ID. Parent id will be used to access information from slqite


```
TimsTOFXICDB timsTOFXICDB = new TimsTOFXICDB(path);
```
Create TimsTOFXICDB, with path to sqlite file


```
TimstofQueryResult result = timsTOFXICDB.queryAndSumParentId(parentID)
List<Pair< Double,Double>> list = result.getSummedList();
double rettime = result.retTime
```
Use precursor Id from index file to query sqllite db. It will return TimstofQueryResult object. The object contains a list of summed peaks and retention time for that precursor 
