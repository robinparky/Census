#java CensusRatioFilter 0.5 2.0 0.5 $* > low_quality_peptides.txt
#arg1 : no_filter_census-out.txt
#arg2 : census_chro.xml
#arg3 : threshold for all_none score

echo "generating low_quality_peptides.txt..."
#java CensusRatioFilter 0.5 2.0 0.5 $1 > low_quality_peptides.txt
java CensusRatioFilter 0.2 5.0 0.5 $1 > low_quality_peptides.txt
echo "done."

echo "generating census_chro_low_quality.xml..."
#Usage: FilterLowQuanityChro census_chro.xml low_quality_filtered_peptides.txt
java FilterLowQualityChro $2 low_quality_peptides.txt
echo "done."

echo "generating AN_candidates.txt..."
java FindAllNothing census_chro_low_quality.xml $3 > AN_candidates.txt
echo "done."

