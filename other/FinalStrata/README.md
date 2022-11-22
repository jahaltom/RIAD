```
cat FINAL_STRATA.txt | while read i;do
  sed "s/FINAL_STRATA/$i/g" InferAncestry.FINAL_STRATA.py > InferAncestry.FINAL_STRATA.$i.py
  snakemake -j 900 -k -s InferAncestry.FINAL_STRATA.$i.py --cluster "sbatch -t 02:00:00  -c 7 -p RM-shared"
  rm $i_InferAncestry.FINAL_STRATA.py
done
```
