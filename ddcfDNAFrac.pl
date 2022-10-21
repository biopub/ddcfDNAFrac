#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin);
use SOT;

&main();

sub main
{
	my $gdnaID;
	my $gdnaR1;
	my $gdnaR2;
	my $cfdnaID;
	my $cfdnaR1;
	my $cfdnaR2;
	my $config;
	my $run;
	my $help;
	GetOptions(
		"gdnaID=s" => \$gdnaID,
		"gdnaR1=s" => \$gdnaR1,
		"gdnaR2=s" => \$gdnaR2,
		"cfdnaID=s" => \$cfdnaID,
		"cfdnaR1=s" => \$cfdnaR1,
		"cfdnaR2=s" => \$cfdnaR2,
		"config=s" => \$config,
		"run" => \$run,
		"help" => \$help
	) or &usage();	


	$config ||= "$Bin/config.txt";
	&usage if $help;
	
	if ( $gdnaID & $gdnaR1 & $cfdnaID & $cfdnaR1 )
	{
		if( !$gdnaR2 & !$cfdnaR2 )
		{
			&start_SE($gdnaID,$gdnaR1,$cfdnaID,$cfdnaR1,$config);
		}else{
			&start_PE($gdnaID,$gdnaR1,$gdnaR2,$cfdnaID,$cfdnaR1,$cfdnaR2,$config);
		}
	}else{
		&usage;
	}

	#&run($gdnaID,$cfdnaID) if $run;
}

sub usage
{
	print<<EOF;
usage: perl $0 --gdnaID <STR> --gdnaR1 <file> --gdnaR2 [file] --cfdnaID <STR> --cfdnaR1 <file> --cfdnaR2 [file]
	   perl $0 --help
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

Eg:
   #quantification
   perl $0 --gdnaID gdnaID --gdnaR1 gdnaID.R1.fq.gz --gdnaR2 gdnaID.R2.fq.gz --cfdnaID cfdnaID --cfdnaR1 cfdnaID.R1.fq.gz --cfdnaR2 cfdnaID.R2.fq.gz --run
   #help
   perl $0 --help
EOF
exit(0)
}
