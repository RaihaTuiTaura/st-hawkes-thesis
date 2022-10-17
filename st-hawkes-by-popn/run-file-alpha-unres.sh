#!/bin/bash
for ind in {1..28} ;
do
	echo $ind
	cp generic-job-alpha-unres.sh job-alpha-unres-$ind.sh
	sed -i -e 's/index_value/'"$ind"'/g' job-alpha-unres-$ind.sh
	qsub job-alpha-unres-$ind.sh
done
