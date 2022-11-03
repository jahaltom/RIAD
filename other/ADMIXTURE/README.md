
1KGP.pop much be in current working dir. Same as sample.py and admix.sh.

```
cat ids.txt | while read i;do
	cat admix.sh | sed "s/SAMPLE/$i/g" > output/$i/$i.sh
	cat sample.py | sed "s/SAMPLE/$i/g" > output/$i/$i.py
	sbatch output/$i/$i.sh
done
```

To have ADMIXTURE in RIA output format:
```
awk -F "\t" 'OFS="\t" {print $2,$2,$1,$7,$5,$4,$3,$6}' ADMIXTURE_Results.tsv > tail
more SuperpopulationChrAll.PC20SVMResults | head -1 > head
cat head tail  > SuperpopulationChrADMIXTURE
rm tail head
```
