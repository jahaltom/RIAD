source activate Ancestry

#5% of list
percentualHead() {
  head -n "$(( "$(wc -l < "$2")" * "$1" / 100  ))" "$2"
}


#Extract SAMPLE from 1KGP-SAMPLE vcf
bcftools view -s SAMPLE --threads THREADS --output-type z OUTPUT/SAMPLE/ChrAll_.SAMPLE.vcf.gz > OUTPUT/SAMPLE/SAMPLE.vcf.gz

#Gather coordinates for all HomRef calls and randomly suffel. Take 5% for targets.
bcftools query -f'[%CHROM:%POS:%REF:%ALT\n]' -i'GT="RR"' OUTPUT/SAMPLE/SAMPLE.vcf.gz   > OUTPUT/SAMPLE/HomRef
#bcftools query -f'%CHROM:%POS:%REF:%ALT\n]' -i'GT="het"' OUTPUT/SAMPLE/SAMPLE.vcf.gz  > OUTPUT/SAMPLE/Het
#bcftools query -f'[%CHROM:%POS:%REF:%ALT\n]' -i'GT="AA"' OUTPUT/SAMPLE/SAMPLE.vcf.gz  > OUTPUT/SAMPLE/HomAlt
cat OUTPUT/SAMPLE/HomRef | shuf  > OUTPUT/SAMPLE/cords
#OUTPUT/SAMPLE/Het OUTPUT/SAMPLE/HomAlt
rm OUTPUT/SAMPLE/HomRef
#OUTPUT/SAMPLE/Het OUTPUT/SAMPLE/HomAlt
percentualHead 5 OUTPUT/SAMPLE/cords  > OUTPUT/SAMPLE/targets
rm OUTPUT/SAMPLE/cords


#Targets file CHROM POS
cat  OUTPUT/SAMPLE/targets | sed 's/:/ /g' | awk '{print $1"\t"$2}'  > OUTPUT/SAMPLE/targets_file

#create id in SAMPLE vcf
bcftools annotate --set-id '%CHROM\:%POS\:%REF\:%ALT' --threads THREADS OUTPUT/SAMPLE/SAMPLE.vcf.gz > OUTPUT/SAMPLE/SAMPLE.vcf
rm OUTPUT/SAMPLE/SAMPLE.vcf.gz

#Gather targets and non target vcf.
bcftools view -T OUTPUT/SAMPLE/targets_file OUTPUT/SAMPLE/SAMPLE.vcf  > OUTPUT/SAMPLE/targets.vcf
bcftools view -T ^OUTPUT/SAMPLE/targets_file OUTPUT/SAMPLE/SAMPLE.vcf  > OUTPUT/SAMPLE/non-targets.vcf



#Split targets into two files.
split -n l/2  OUTPUT/SAMPLE/targets OUTPUT/SAMPLE/


#Add in fake Hets
cat OUTPUT/SAMPLE/aa | while read i;do
     sed -i "/$i/ s=0/0=0/1="  OUTPUT/SAMPLE/targets.vcf
done

#Add in fake HomAlts
cat OUTPUT/SAMPLE/ab | while read i;do
     sed -i "/$i/ s=0/0=1/1="  OUTPUT/SAMPLE/targets.vcf
done

grep "#" OUTPUT/SAMPLE/targets.vcf > OUTPUT/SAMPLE/header
cat OUTPUT/SAMPLE/targets.vcf  OUTPUT/SAMPLE/non-targets.vcf | grep -v "#" | sort -k1,1 -k2,2n > OUTPUT/SAMPLE/variants
cat OUTPUT/SAMPLE/header OUTPUT/SAMPLE/variants   | bgzip > OUTPUT/SAMPLE/Temp_ChrAll_Noise.SAMPLE.vcf.gz
#Remove INFO. Otherwise will mess up merge.
bcftools annotate -x ID,INFO --threads THREADS  --output-type z OUTPUT/SAMPLE/Temp_ChrAll_Noise.SAMPLE.vcf.gz  > OUTPUT/SAMPLE/ChrAll_Noise.SAMPLE.vcf.gz
bcftools index -t OUTPUT/SAMPLE/ChrAll_Noise.SAMPLE.vcf.gz
rm OUTPUT/SAMPLE/header OUTPUT/SAMPLE/variants OUTPUT/SAMPLE/targets.vcf  OUTPUT/SAMPLE/non-targets.vcf OUTPUT/SAMPLE/Temp_ChrAll_Noise.SAMPLE.vcf.gz


#Remove SAMPLE from 1KGP-SAMPLE vcf
bcftools view -s ^SAMPLE --threads THREADS  --output-type z OUTPUT/SAMPLE/ChrAll_.SAMPLE.vcf.gz  > OUTPUT/SAMPLE/ChrAll_.1KGP.SAMPLE.vcf.gz
bcftools index -t OUTPUT/SAMPLE/ChrAll_.1KGP.SAMPLE.vcf.gz

bcftools merge OUTPUT/SAMPLE/ChrAll_.1KGP.SAMPLE.vcf.gz OUTPUT/SAMPLE/ChrAll_Noise.SAMPLE.vcf.gz --output-type z --threads THREADS > OUTPUT/SAMPLE/ChrAll.SAMPLE.Noise.vcf.gz
bcftools index -t OUTPUT/SAMPLE/ChrAll.SAMPLE.Noise.vcf.gz
rm OUTPUT/SAMPLE/ChrAll_Noise.SAMPLE.vcf.gz* OUTPUT/SAMPLE/ChrAll_.1KGP.SAMPLE.vcf.gz*

