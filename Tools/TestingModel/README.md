Delete Kenya, Uganda, and Russia studies. 



cat ChrAll_PC20_DPLT5 |  grep -v Kenya | grep -v "Uganda" | grep -v "Russia" > ChrAll_PC5_DPLT5.tsv


###Best Number of PCs

for i in {5,10,15,20}
do
  echo PC"$i"
  cat ChrAll_PC"$i"_DPLT5.tsv | awk -F '\t' '{sum+=$5/100*$2;} END{print "SVM: "sum/922*100}'
  cat ChrAll_PC"$i"_DPLT5.tsv | awk -F '\t' '{sum+=$4/100*$2;} END{print "NN: "sum/922*100}'
  cat ChrAll_PC"$i"_DPLT5.tsv | awk -F '\t' '{sum+=$3/100*$2;} END{print "RF: "sum/922*100}'
done






