sample=$1
perl vcf2tab.pl $sample.mp.vcf /data/pipeline/NGS_Target_Sequencing/Organ_Transplant/config/pt2.SNP.GRCh38.bed $sample.mp.tab

#checking for uncomplete lines
less $sample.mp.tab|awk 'NF>=13' >t
mv t $sample.mp.tab
