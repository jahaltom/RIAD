
1KGP.pop much be in current working dir. Same as sample.py and admix.sh.

```
cat RAids.txt | while read i;do
	cat admix.sh | sed "s/SAMPLE/$i/g" > output/$i/$i.sh
	cat sample.py | sed "s/SAMPLE/$i/g" > output/$i/$i.py
	sbatch output/$i/$i.sh
done
```
