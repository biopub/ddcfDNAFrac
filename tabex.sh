cfdnaID=$1
recipientID=$2

parameter="--bam mapped/${cfdnaID}.target.sort.bam --snptab ${cfdnaID}.mp.tab --recipientID genotype/${recipientID}.mp.tab"

perl /nfs1/public2/User/lxf/src/extract_organ_source_read.pl ${parameter}

Rscript /nfs1/public2/User/lxf/src/ddcfDNAFrac.MPCR.R ${cfdnaID}.mp.tab.ex ${cfdnaID} ${recipientID}