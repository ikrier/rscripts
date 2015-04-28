echo $1
echo $2
table_annovar.pl  $1  /data/software/annovar/humandb/ -vcfinput -buildver hg19 -outfile $2 -remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2014oct_all,1000g2014oct_afr,1000g2014oct_eas,1000g2014oct_eur,snp138,ljb26_all,cosmic70 -operation g,r,r,f,f,f,f,f,f,f,f -nastring .
