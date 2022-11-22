cat FINAL_STRATA.txt | while read i;do
  sed "s/FINAL_STRATA/$i/g" InferAncestry.FINAL_STRATA.py > $i_InferAncestry.FINAL_STRATA.py
  snakemake -j 900 -k -s $i_InferAncestry.FINAL_STRATA.py --cluster "sbatch -t 02:00:00  -c 7 -p RM-shared"
  rm $i_InferAncestry.FINAL_STRATA.py
done
