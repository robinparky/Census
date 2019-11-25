grep "<chro" census_chro.xml | awk -F';' 'BEGIN{sum=0;} {sum = sum + NF -2;print sum" "(NF-2);} END {print sum}'
