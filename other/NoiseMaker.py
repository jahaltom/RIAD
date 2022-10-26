#Grab all variants with a minor allele frequency =>0.05
#bcftools view -i 'AF>0 && AF<0.10' output/SRR6370304/sample.vcf.gz


#5% of list
percentualHead() {
  head -n "$(( "$(wc -l < "$2")" * "$1" / 100  ))" "$2"
}


#Gather single sample vcf for RNA-Seq sample
bcftools view -s SRR493967 --threads 7 --output-type z output/SRR6370304/Chr1.SRR6370304.vcf.gz > output/SRR6370304/sample.vcf.gz

#Gather cords and randomly suffel. Take 5% for targets.
bcftools query -f'[%CHROM:%POS:%REF:%ALT %SAMPLE %GT\n]' -i'GT="RR"' output/SRR6370304/sample.vcf.gz | awk '{print $1}'  > output/SRR6370304/HomRef
#bcftools query -f'%CHROM:%POS:%REF:%ALT %SAMPLE %GT\n]' -i'GT="het"' output/SRR6370304/sample.vcf.gz | awk '{print $1"_Het"}'  > output/SRR6370304/Het
#bcftools query -f'[%CHROM:%POS:%REF:%ALT %SAMPLE %GT\n]' -i'GT="AA"' output/SRR6370304/sample.vcf.gz | awk '{print $1"_HomAlt"}'  > output/SRR6370304/HomAlt
cat output/SRR6370304/HomRef | shuf  > output/SRR6370304/cords
#output/SRR6370304/Het output/SRR6370304/HomAlt  
rm output/SRR6370304/HomRef 
#output/SRR6370304/Het output/SRR6370304/HomAlt
percentualHead 5 output/SRR6370304/cords  > output/SRR6370304/targets
rm output/SRR6370304/cords


#create id in sample vcf 
bcftools annotate --set-id '%CHROM\:%POS\:%REF\:%ALT' --threads 7 output/SRR6370304/sample.vcf.gz > output/SRR6370304/All.vcf
rm output/SRR6370304/sample.vcf.gz


#Split targets into two files.
split -n 2 output/SRR6370304/targets

#Add in fake Hets
cat xaa | while read i;do
     sed -i "/$i/ s=0/0=0/1="  output/SRR6370304/All.vcf
done

#Add in fake HomAlts
cat xab | while read i;do
     sed -i "/$i/ s=0/0=1/1="  output/SRR6370304/All.vcf
done




 
############
bcftools query -f'[%AF %GT\n]' output/SRR6370304/sample.vcf.gz > AF_vs_GT
