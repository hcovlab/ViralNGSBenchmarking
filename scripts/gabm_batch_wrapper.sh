#!/bin/bash

#### go through every directory and do stuff
for d in */ ; do
    cd $d
    bash /home/hazaihiv/Analyses/myscripts/gabm_batch.sh
    cd ..
done

bash ~/Analyses/myscripts/gather_results.sh