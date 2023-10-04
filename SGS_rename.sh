for d in */ ; do
    cd ${d%/}
    file=$(ls)
    #file_r=$(ls | sed 's/\./_/g' | sed 's/_/\./10')
    file_r=${file%.*}_SGSref.fasta
    mv $file $file_r
    cd ..
done