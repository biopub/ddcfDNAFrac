# ddcfDNAFrac

# install
git clone https://github.com/biopub/ddcfDNAFrac.git

# usage
perl ddcfDNAFrac.pl --gdnaID <STR> --gdnaR1 <file> --gdnaR2 [file] --cfdnaID <STR> --cfdnaR1 <file> --cfdnaR2 [file]
Params:
--gdnaID	STR	-	gDNA ID <required>
--gdnaR1	FILE	-	read1 of paired-end data for gdna <required>
--gdnaR2	FILE	-	read2 of paired-end data for gdna [optional]
--cfdnaID	STR	-	cfDNA ID <required>
--cfdnaR1	FILE	-	read1 of paired-end data for cfdna <required>
--cfdnaR2	FILE	-	read2 of paired-end data for cfdna [optional]
--config	FILE	- FILE of configuration.(defaut: cofig.txt)
--run		- Start running directly
--help		- Help document.
