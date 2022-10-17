#!/bin/bash
for ind in {1..28} ;
do
	echo $ind
	cp generic-job-adm-only.sh job-adm-only-$ind.sh
	sed -i -e 's/index_value/'"$ind"'/g' job-adm-only-$ind.sh
	qsub job-adm-only-$ind.sh
done
