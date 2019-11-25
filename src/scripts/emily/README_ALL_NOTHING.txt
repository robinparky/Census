run_all_none.sh : a script to run all followings

1. filter all_nothing candidates from census-out-nofiltering.txt using CensusRatioFilter.java
Usage: CensusRatioFilter lower_bound_ratio upper_bound_ratio R_sqrt_threshold census_out.txt
ratio can be 0.5 and 2

2. create filtered_chro.xml from FilterLowQuanlityChro.java

3. run FindAllNothing.java on filtered_chro.xml file 
