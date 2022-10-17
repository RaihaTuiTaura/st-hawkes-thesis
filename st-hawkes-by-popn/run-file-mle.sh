#!/bin/bash
for ind in {1..4} ;
do
	echo $ind
	cp generic-job-mle.sh job-mle-$ind.sh
	sed -i -e 's/index_value/'"$ind"'/g' job-mle-$ind.sh
	qsub job-mle-$ind.sh
done
