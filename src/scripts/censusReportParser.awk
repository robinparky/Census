
cat "$1" | awk '{


    if($0 ~ /^P/)
    {
#	count=1;
	accession = $2;

#	split($6,sfile_name,".");
    }

    if($0 ~ /^S/)
    {
	print $0"\t"accession;
    }
}'

