BEGIN {
    OFS="\t";
}

{
    sj=$1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6;
    nSamples[sj]++;
    nU[sj]+=$7;
    nM[sj]+=$8;
    if (nO[sj]<$9) nO[sj]=$9;
};

END {
    for (sj in nSamples) print sj,nU[sj],nM[sj],nO[sj],nSamples[sj];
}
