package NGS_PIP;

use warnings;
use Try::Tiny;
use NGS_TOOLS;
#use strict;
require Exporter;
@ISA = qw(Exporter);

@EXPORT = qw(init download_from_oss shell_for_SOT shell_for_TUMOR manual_start);
@EXPORT_OK = qw(%TASKS %CONFIG %STATUS %STATUS_FH);

my %TASKS = ();
my %TASKS2TARGET = ();
my %OSS_FILEs = ();
my %STATUS = ();
my %CONFIG = ();
my %STATUS_FH = ();
my %okstatu=(
	"00"=>"Both R1 and R2 haven't been downloaded",
	"12"=>"Both R1 and R2 have been downloaded",
	"10"=>"Only R1 have been downloaded",
	"02"=>"Only R2 have been downloaded"
);
my $background = 2;#how many QZ sample download for background analysis

my ($ossBucketP,$ossaccount) = ('','');
my $manufacturer = '';

my $regulex_QY='S\d+W?G?S?Q?\d?(QY)Q?[BWRXN]?(\d+)[HQG]\d+';
my $regulex_IF='S\d+W?G?S?Q?\d?Q?Y?(IF)X?(CS|T)(\d+)';
my $regulex_ZL='(ZL)+CS[XC]?(\d+)';
my $regulex_NIPT='S\d+(WS)(\d+)';
my $regulex_QZ='S\d+(QZ)(\d+)S';
my $regulex_YS='S\d+(YS)(\d+)';
my $regulex_QH='S\d+Q?\d?(QH)[N]?(\d+)[HQG]\S+';
sub init
{
	#use Getopt::Long;
	my $table=shift;
	my $osspath=shift;
	my $conf=shift;
	my $clean=shift;
	my $run=shift;
	my $pro=shift;
	my $bucket=shift;
	my $Ncores=shift;
	my $runlog=shift;
	my $docheck=shift;
	my $skipIFD=shift;

	&load_config($conf);	

	($ossBucketP,$ossaccount) = &ossRR($bucket);

	#my @months = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );
	#my @days = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
	#my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
	#$mon+=1;
	#my $output_prefix="SOT_$year\_$mon\_$mday\_$days[$wday]";
	my $output_prefix = &set_datelike_prefix();
	unless( -s "$CONFIG{SOT_PROJECT_DIR}/SHELL/$output_prefix")
	{
		mkdir("$CONFIG{SOT_PROJECT_DIR}/SHELL/$output_prefix",0755) or die $!;
	}
	unless( -s "$CONFIG{SOT_PROJECT_DIR}/report/$output_prefix")
	{
		mkdir("$CONFIG{SOT_PROJECT_DIR}/report/$output_prefix", 0755) or die $!;
	}

	my @savedir = &get_savedir($osspath,$output_prefix);

	if($table eq "")
	{
		&load_sample_from_oss(\@savedir,$output_prefix);
	}else{
		&load_sample_from_table($table,$output_prefix);
	}
	my $tasks_n = keys %TASKS;
	my $process_n = 0;
	foreach my $sample(keys %TASKS)
	{
		$process_n++;
		print STDERR "process:\t$process_n/$tasks_n\t";
		unless($clean)
		{
			&load_status($sample);
		}

		my $is_download_ok = &download_from_oss($sample,$runlog,$docheck);
		next unless($is_download_ok eq "12");
		unless(exists $STATUS{$sample}{'INITIAL'} && $STATUS{$sample}{'INITIAL'} eq 'DONE')
		{
			if($sample=~/QY/ || $sample=~/IF/)
			{
				&shell_for_SOT($sample,$TASKS2TARGET{$sample},$output_prefix,$Ncores,$skipIFD) if(exists $TASKS{$sample});
				unless(exists $STATUS{$sample}{'RUN'} && $STATUS{$sample}{'RUN'} eq 'DONE')
				{
					&run("$TASKS2TARGET{$sample}.script.sh",$pro) if($run && exists $TASKS{$sample} && (! -e "$TASKS2TARGET{$sample}.script.sh.statu"));
					$STATUS{$sample}{'RUN'} = 'DONE';
				}
			}elsif($sample=~/ZL/){
				&shell_for_TUMOR($sample) if($run && exists $TASKS{$sample} && (! -e "$sample.script.sh.statu"));
			}elsif($sample=~/QZ/){
				&shell_for_QZ($sample,$output_prefix,$Ncores,$skipIFD) if(exists $TASKS{$sample});
				unless(exists $STATUS{$sample}{'RUN'} && $STATUS{$sample}{'RUN'} eq 'DONE')
				{
					&run("$TASKS2TARGET{$sample}.script.sh",$pro) if($run && exists $TASKS{$sample} && (! -e "$TASKS2TARGET{$sample}.script.sh.statu"));
					$STATUS{$sample}{'RUN'} = 'DONE';
				}
			}
		}
		#====sleep================
		for(my $t=0;$t<1;$t++)
		{
			print STDERR "......... ";
			sleep(1);
		}
		print STDERR "\n";
	}

	#close statu handle
	foreach my $sp(keys %STATUS_FH)
	{
		my $statu_handle = $STATUS_FH{$sp};
		if(exists $STATUS{$sp})
		{
			my %tmp = %{$STATUS{$sp}};
			foreach my $k(keys %tmp)
			{
				print $statu_handle "$k\t$tmp{$k}\n";
			}
		}
		$STATUS_FH{$sp}->close if ($STATUS_FH{$sp}->opened);
	}
}

sub download_from_oss
{
	my ($sample,$runlog,$docheck)=@_;

	my $basedir= "";

	my @R1R2 = @{$OSS_FILEs{$sample}};
	#print "==========\n",join("\n",@R1R2),"\n";
	my @a = split("/",$R1R2[0]);
	
	my $dir=join("/",@a[3..$#a-1]);
	my $isok = "";

	my $downloadpath = join("/",@a[0..$#a-1]);
	if($sample=~/$regulex_QY/ || $sample=~/$regulex_IF/ || $sample=~/$regulex_QH/)
	{
		print STDERR "downloading $sample\n" if($runlog);
		$basedir = $CONFIG{'SOT_RAWDATA_DIR'};
		#no statu about download 
		unless(exists $STATUS{$sample}{'DOWNLOAD'} && $STATUS{$sample}{'DOWNLOAD'} eq 'DONE')
		{
			#but downloaded and md5check is ok
			$isok = &md5check($sample,$basedir,$dir,$downloadpath);
			chdir($CONFIG{'SOT_RAWDATA_DIR'});
			if($isok eq "12")
			{
				print STDERR "$sample md5check is ok\n" if($runlog);
				$STATUS{$sample}{'DOWNLOAD'} = 'DONE';
			}else{# no statu and md5check is not ok then
				print STDERR "isok=$isok means \"$okstatu{$isok}\"\n";
				if($isok eq "02")
				{
					system("ossutil64 cp -j 10 $ossaccount -rf $R1R2[0] --checkpoint-dir /data/raw_data/check_point . ");# or warn $!;
					system("chmod 444 $basedir/$dir/$sample*R1*gz");
				}elsif($isok eq "10"){
					system("ossutil64 cp -j 10 $ossaccount -rf $R1R2[1] --checkpoint-dir /data/raw_data/check_point . ");# or warn $!;
					system("chmod 444 $basedir/$dir/$sample*R2*gz");

				}elsif($isok eq "00"){
					#print "ossutil64 cp -j 10 $ossaccount -rf $R1R2[0] --checkpoint-dir /data/raw_data/check_point .\n";
					system("ossutil64 cp -j 10 $ossaccount -rf $R1R2[0] --checkpoint-dir /data/raw_data/check_point . ");# or warn $!;
					system("chmod 444 $basedir/$dir/$sample*R1*gz");
					#print "ossutil64 cp -j 10 $ossaccount -rf $R1R2[1] --checkpoint-dir /data/raw_data/check_point .\n";
					system("ossutil64 cp -j 10 $ossaccount -rf $R1R2[1] --checkpoint-dir /data/raw_data/check_point . ");# or warn $!;
					system("chmod 444 $basedir/$dir/$sample*R2*gz");
				}
				if($docheck)
				{
					$isok = &md5check($sample,$basedir,$dir,$downloadpath);
				}else{
					$isok = "12";
				}
				$STATUS{$sample}{'DOWNLOAD'} = 'DONE' if($isok eq "12") ;
			}
		}
		#}elsif($sample =~ /(ZL)[XC]?(\d+)/){
	}elsif($sample =~ /$regulex_ZL/){

		print STDERR "downloading $sample\n" if($runlog);
		$basedir = $CONFIG{'TUMOR_RAWDATA_DIR'};
		#no statu about download #
		unless(exists $STATUS{$sample}{'DOWNLOAD'} && $STATUS{$sample}{'DOWNLOAD'} eq 'DONE')	
		{
			#but downloaded and md5check is ok#		
			$isok = &md5check($sample,$basedir,$dir,$downloadpath);
			chdir($CONFIG{'TUMOR_RAWDATA_DIR'});
			if($isok eq "12")
			{
				$STATUS{$sample}{'DOWNLOAD'} = 'DONE';
			}else{# no statu and md5check is not ok then
				if($isok eq "02")
				{
					system("ossutil64 cp -j 10 $ossaccount -rf $R1R2[0] .");
					system("chmod 444 $basedir/$dir/$sample*R1*gz");
				}elsif($isok eq "10")
				{
					system("ossutil64 cp -j 10 $ossaccount -rf $R1R2[1] .");
					system("chmod 444 $basedir/$dir/$sample*R2*gz");
				}else{
					system("ossutil64 cp -j 10 $ossaccount -rf $R1R2[0] .");
					system("chmod 444 $basedir/$dir/$sample*R1*gz");
					system("ossutil64 cp -j 10 $ossaccount -rf $R1R2[1] .");
					system("chmod 444 $basedir/$dir/$sample*R2*gz");
				}
				$isok = &md5check($sample,$basedir,$dir,$downloadpath);
				$STATUS{$sample}{'DOWNLOAD'} = 'DONE' if($isok eq "12");

			}
		}
	}elsif($sample =~ /$regulex_NIPT/){
		print STDERR "downloading $sample\n" if($runlog);
		$basedir = $CONFIG{'NIPT_RAWDATA_DIR'};
		#no statu about download #
		unless(exists $STATUS{$sample}{'DOWNLOAD'} && $STATUS{$sample}{'DOWNLOAD'} eq 'DONE')	
		{
			#but downloaded and md5check is ok#		
			$isok = &md5check($sample,$basedir,$dir,$downloadpath);
			print "isok=$isok\n";#$R1R2[0]\n$R1R2[1]\n";
			chdir($CONFIG{'NIPT_RAWDATA_DIR'});
			if($isok eq "12")
			{
				$STATUS{$sample}{'DOWNLOAD'} = 'DONE';
			}else{# no statu and md5check is not ok then
				if($isok eq "02")
				{
					system("ossutil64 cp -j 10 $ossaccount -rf $R1R2[0] .");
					system("chmod 444 $basedir/$dir/$sample*R1*gz");
				}elsif($isok eq "10")
				{
					system("ossutil64 cp -j 10 $ossaccount -rf $R1R2[1] .");
					system("chmod 444 $basedir/$dir/$sample*R2*gz");
				}else{
					system("ossutil64 cp -j 10 $ossaccount -rf $R1R2[0] .");
					system("chmod 444 $basedir/$dir/$sample*R1*gz");
					system("ossutil64 cp -j 10 $ossaccount -rf $R1R2[1] .");
					system("chmod 444 $basedir/$dir/$sample*R2*gz");
				}
				if($docheck)
				{
					$isok = &md5check($sample,$basedir,$dir,$downloadpath);
				}else{
					$isok = "12";
				}
				$STATUS{$sample}{'DOWNLOAD'} = 'DONE' if($isok eq "12");
			}
		}
	}elsif($sample =~ /$regulex_QZ/){

		if($background>0)
		{
			print STDERR "downloading $sample\n" if($runlog);
			$background-=1;
			$basedir = $CONFIG{'QZ_RAWDATA_DIR'};
			#no statu about download #
			unless(exists $STATUS{$sample}{'DOWNLOAD'} && $STATUS{$sample}{'DOWNLOAD'} eq 'DONE')	
			{
				#but downloaded and md5check is ok#		
				$isok = &md5check($sample,$basedir,$dir,$downloadpath);
				chdir($CONFIG{'QZ_RAWDATA_DIR'});
				if($isok eq "12")
				{
					$STATUS{$sample}{'DOWNLOAD'} = 'DONE';
				}else{# no statu and md5check is not ok then
					if($isok eq "02")
					{
						system("ossutil64 cp -j 10 $ossaccount -rf $R1R2[0] .");
						system("chmod 444 $basedir/$dir/$sample*R1*gz");
					}elsif($isok eq "10")
					{
						system("ossutil64 cp -j 10 $ossaccount -rf $R1R2[1] .");
						system("chmod 444 $basedir/$dir/$sample*R2*gz");
					}else{
						system("ossutil64 cp -j 10 $ossaccount -rf $R1R2[0] .");
						system("chmod 444 $basedir/$dir/$sample*R1*gz");
						system("ossutil64 cp -j 10 $ossaccount -rf $R1R2[1] .");
						system("chmod 444 $basedir/$dir/$sample*R2*gz");
					}
					$isok = &md5check($sample,$basedir,$dir,$downloadpath);
					$STATUS{$sample}{'DOWNLOAD'} = 'DONE' if($isok eq "12");
				}
			}
		}else{
			delete $TASKS{$sample};
		}
	}#QZ as background
	return $isok;
}

sub start
{
	my ($gdnaID,$gdnaR1,$gdnaR2,$cfdnaID,$cfdnaR1,$cfdnaR2,$config)=@_;

	my $output_prefix = &set_datelike_prefix();
	&load_config($config);


	unless($r1 && $r2 && $sample)
	{
		print STDERR "program <read1> <read2> <sample> <pro>[SOT]\n";
		die "manual start analysis failed. Exiting...";
	}

	open SHELL, ">$sample.script.sh" or die $!;
	if($type eq "WGS")
	{
		my $r1fastp="$sample\_fastp_R1.fastq.gz";
		my $r2fastp="$sample\_fastp_R2.fastq.gz";

		print SHELL "#1 3G filter\n";
		print SHELL "$CONFIG{FASTP} -i $r1 -I $r2 -o $r1fastp -O $r2fastp -j $sample.json -h $sample.html  -l 31 -w 16 -Q\n\n";

		print SHELL "#1#2 40G kraken\n";
		print SHELL "export KRAKEN_DB_PATH=/data/soft/kraken/\n";
		print SHELL "if [ ! -f $sample.kraken ]\nthen\n";
		print SHELL "\t$CONFIG{KRAKEN} -db $CONFIG{KRAKEN_DB} --threads 16 --fastq-input --gzip-compressed --output $sample.kraken --preload $r1fastp $r2fastp\nfi\n\n";

		print SHELL "#3#4 5G BWAalign\n";
		print SHELL "$CONFIG{'PATHOGEN_REF_READ'} $sample.kraken.cdc $sample $sample\_fastp_R1.fastq.gz $sample\_fastp_R2.fastq.gz $CONFIG{'KRAKEN_PATH'}/$CONFIG{'KRAKEN_DB'} \n";
		#print SHELL "$CONFIG{'BWA'} index $sample.aligned_cdc.fa\n";
		#print SHELL "$CONFIG{'BWA'} mem $CONFIG{'CDC_TAXDB'} $sample\_CDC_R1.fastq $sample\_CDC_R2.fastq -t $CONFIG{'Ncores'}  -R \"\@RG:\\tID:$sample\\tSM:$sample\" -v 1 \| $CONFIG{'SAMTOOLS'} sort -T tmp -o $sample.cdc.sort.bam -\n";
		print SHELL "$CONFIG{'BWA'} mem $CONFIG{'CDC_TAXDB'} $sample\_CDC_R1.fastq $sample\_CDC_R2.fastq -t 16  -R \"\@RG:\\tID:$sample\\tSM:$sample\" -v 1 \| $CONFIG{'SAMTOOLS'} sort -T tmp -o $sample.cdc.sort.bam -\n";
		print SHELL "$CONFIG{'PICARD'} MarkDuplicates I=$sample.cdc.sort.bam O=$sample.cdc.picard.sort.bam R=$CONFIG{'CDC_TAXDB'} VALIDATION_STRINGENCY=SILENT QUIET=true REMOVE_DUPLICATES=true M=$sample.picard.txt\n";

		#print "Script made. Run with:\n$sample.script.sh\n\n";

	}elsif($type eq "HYB")
	{
		#set parameter
		my ($readLength,$barcodeLength,$repFilt)=cal_barcode_read_len($r1,$CONFIG{'SPACER'});
		print SHELL "#1 3G filter\n";
		print SHELL "barcodeLength=$barcodeLength\n";
		print SHELL "if [ ! -f $sample.withtag.R1.fq.gz ]\nthen\n";
		print SHELL "\tif [ \$barcodeLength -gt 2 ]\nthen\n";
		print SHELL "\t\tperl $CONFIG{tag_to_header} $r1 $r2 $sample\n\telse\n";
		print SHELL "\t\tln -s $r1 $sample.withtag.R1.fq.gz\n\t\tln -s $r2 $sample.withtag.R2.fq.gz\n\tfi\nfi\n\n";
		unless($skipIFD)
		{
			print SHELL "#1#2 40G kraken\n";
			print SHELL "export KRAKEN_DB_PATH=/data/soft/kraken/\n";
			print SHELL "if [ ! -f $sample.kraken ]\nthen\n";
			print SHELL "$CONFIG{'KRAKEN'} -db $CONFIG{'KRAKEN_DB'} --threads $CONFIG{'Ncores'} --fastq-input --output $sample.kraken --paired --preload --gzip-compressed $sample.withtag.R1.fq.gz $sample.withtag.R2.fq.gz\n";
			print SHELL "awk \'NR==FNR{a[\$1]=\$1;}NR!=FNR{if(a[\$3]){print;}}\' $CONFIG{'CDC_TAXID'} $sample.kraken >$sample.kraken.cdc\n";
			print SHELL "$CONFIG{'KRAKEN-FILTER'} --db $CONFIG{'KRAKEN_DB'} --threshold 0.3 $sample.kraken.cdc  > $sample.kraken.cdc.2\n";
			print SHELL "$CONFIG{'KRAKEN-REPORT'} --db $CONFIG{'KRAKEN_DB'} $sample.kraken.cdc.2 | sort > $sample.kraken.cdc.2.report.txt\nfi\n\n";


			print SHELL "#2#3 20G CDCalign\n";

			unless($sample =~ /QYX/)
			{
				print SHELL "if [ ! -f $sample\_CDC_R1.fastq ]\nthen\n";
				print SHELL "\t$CONFIG{'PATHOGEN_REF_READ'} $sample.kraken.cdc $sample $sample.withtag.R1.fq.gz $sample.withtag.R2.fq.gz $CONFIG{'KRAKEN_PATH'}/$CONFIG{'KRAKEN_DEFAULT_DB'}\nfi\n";
				print SHELL "if [ ! -f $sample.cdc.sort.bam ]\nthen\n";
				#print SHELL "\t$CONFIG{'BWA'} mem $CONFIG{'CDC_TAXDB'} $sample\_CDC_R1.fastq $sample\_CDC_R2.fastq -t $CONFIG{'Ncores'} -R \"\@RG\\tID:$sample\\tSM:$sample\" -v 1 | $CONFIG{'SAMTOOLS'} sort -T $sample.cdc.tmp -o $sample.cdc.sort.bam -\nfi\n";
				print SHELL "\t$CONFIG{'BWA'} mem $CONFIG{'CDC_TAXDB'} $sample\_CDC_R1.fastq $sample\_CDC_R2.fastq -t 16 -R \"\@RG\\tID:$sample\\tSM:$sample\" -v 1 | $CONFIG{'SAMTOOLS'} sort -T $sample.cdc.tmp -o $sample.cdc.sort.bam -\nfi\n";
				print SHELL "if [ ! -f $sample.cdc.sscs.sort.bam ]\nthen\n";
				if($barcodeLength >= 2)
				{
					print SHELL "\texport PYTHONPATH=/home/lxf/.PY2.7LIB/lib64/python2.7/site-packages/\n";
					print SHELL "\t/usr/bin/python $CONFIG{'DSpath'}/ConsensusMaker.v2.py --infile $sample.cdc.sort.bam --outfile $sample.cdc.sscs.bam --minmem $CONFIG{'minMem'} --maxmem $CONFIG{'maxMem'} --read_length $readLength --cut_off $CONFIG{'cutOff'} --Ncut_off $CONFIG{'nCutOff'} --read_type $CONFIG{'readTypes'} --filt $CONFIG{'filtersSet'} --isize $CONFIG{'iSize'} --rep_filt $repFilt\n";
					print SHELL "$CONFIG{'SAMTOOLS'} view -h $sample.cdc.sscs.bam|awk '{if (\$6!~/^[[:digit:]]+N[[:digit:]]+D|^[[:digit:]]+P|[[:digit:]]+P\$|^[[:digit:]]+I/) print \$0}' | $CONFIG{'SAMTOOLS'} sort -T tmp -o $sample.cdc.sscs.sort.bam -\@ 8 -\nfi\n\n";
				}else{
					print SHELL "\t$CONFIG{'PICARD'} MarkDuplicates I=$sample.cdc.sort.bam O=$sample.cdc.sort.picard.bam M=$sample.cdc.picard.txt R=$CONFIG{'hg38'} VALIDATION_STRINGENCY=SILENT QUIET=true REMOVE_DUPLICATES=true\n";
					print SHELL "ln -s $sample.cdc.sort.picard.bam $sample.cdc.sscs.sort.bam\nfi\n\n";
				}
			}
		}else{
			#print "echo non-neccessary CDC Aligning on $sample\n";
		}

		#======================XX=========================
		#======================XY=========================
		#genotyping
		print SHELL "#1#4 5G genotyping\n";
		print SHELL "if [ ! -f $sample.sort.bam ]\nthen\n";
		#print SHELL "\t$CONFIG{'BWA'} mem $CONFIG{'hg38'}  $sample.withtag.R1.fq.gz $sample.withtag.R2.fq.gz -t $CONFIG{'Ncores'} -R \"\@RG\\tID:$sample\\tSM:$sample\" -v 1 | $CONFIG{'SAMTOOLS'} sort -T $sample.tmp -o $sample.sort.bam -\nfi\n";
		print SHELL "\t$CONFIG{'BWA'} mem $CONFIG{'hg38'}  $sample.withtag.R1.fq.gz $sample.withtag.R2.fq.gz -t 16 -R \"\@RG\\tID:$sample\\tSM:$sample\" -v 1 | $CONFIG{'SAMTOOLS'} sort -T $sample.tmp -o $sample.sort.bam -\nfi\n";
		print SHELL "if [ ! -f $sample.sort.bam.bai ]\nthen\n\t$CONFIG{'SAMBAMBA'}  index $sample.sort.bam\nfi\n";
		print SHELL "if [ ! -f $sample.UT.sort.bam ]\nthen\n\t$CONFIG{'SAMBAMBA'} view -t 12 -L $CONFIG{'BED_SOT_MTB'} -h -f bam -F \"mapping_quality >= 30 and not (unmapped or secondary_alignment) and not ([XA] != null or [SA] != null) and ((not [MD]=~/\\^/ and [NM]<=2) or (cigar=~/[ID]/ and [NM]<8))\" $sample.sort.bam -o $sample.UT.sort.bam\nfi\n";
		print SHELL "if [ ! -f $sample.UT.sscs.sort.bam ]\nthen\n";
		if($barcodeLength >= 2)
		{
			print SHELL "\texport PYTHONPATH=/home/lxf/.PY2.7LIB/lib64/python2.7/site-packages/\n";
			print SHELL "\t/usr/bin/python $CONFIG{'DSpath'}/ConsensusMaker.v2.py --infile $sample.UT.sort.bam --outfile $sample.UT.sscs.bam --minmem $CONFIG{'minMem'} --maxmem $CONFIG{'maxMem'} --read_length $readLength --cut_off $CONFIG{'cutOff'} --Ncut_off $CONFIG{'nCutOff'} --read_type $CONFIG{'readTypes'} --filt $CONFIG{'filtersSet'} --isize $CONFIG{'iSize'} --rep_filt $repFilt\n";
			print SHELL "\t$CONFIG{'SAMTOOLS'} view -h $sample.UT.sscs.bam|awk '{if (\$6!~/^[[:digit:]]+N[[:digit:]]+D|^[[:digit:]]+P|[[:digit:]]+P\$|^[[:digit:]]+I/) print \$0}' | $CONFIG{'SAMTOOLS'} sort -T $sample.tmp -o $sample.UT.sscs.sort.bam -\@ 8 -\nfi\n";
		}else{
			print SHELL "\t$CONFIG{'PICARD'} MarkDuplicates I=$sample.UT.sort.bam O=$sample.UT.sort.picard.bam M=$sample.picard.txt R=$CONFIG{'hg38'} VALIDATION_STRINGENCY=SILENT QUIET=true REMOVE_DUPLICATES=true\nln -s $sample.UT.sort.picard.bam $sample.UT.sscs.sort.bam\nfi\n";
		}
		#call variant
		#print SHELL "if [ ! -f $sample.sscs.sort.bam.bai ]\nthen\n";
		#print SHELL "\t$CONFIG{'SAMTOOLS'} index $sample.sscs.sort.bam\nfi\n";
		print SHELL "if [ ! -f $sample.mp.vcf ]\nthen\n";
		print SHELL "\t$CONFIG{'SAMTOOLS'} mpileup -A -uvf $CONFIG{'hg38'} -l $CONFIG{'BED_SOT_MTB'} -t DP,AD $sample.UT.sscs.sort.bam | $CONFIG{'BCFTOOLS'} call --multiallelic-caller --keep-alts --targets-file $CONFIG{'BED_SOT_MTB'} > $sample.mp.vcf\nfi\n";
		print SHELL "if [ ! -f $sample.mp.tab ]\nthen\n";
		print SHELL "\tperl $CONFIG{'VCF2TAB'} $sample.mp.vcf $CONFIG{'BED_SOT'} $sample.mp.tab\n";
		print SHELL "\tperl $CONFIG{'VCF2TAB'} $sample.mp.vcf $CONFIG{'BED_MTB'} $sample.mp.mtb.tab\n";
		print SHELL "\tcp $sample.mp.mtb.tab $CONFIG{SOT_PROJECT_DIR}/report/$output_prefix/$sample.Pharmacogenomics.xls\n";
		print SHELL "\tcat $sample.mp.mtb.tab >> $CONFIG{SOT_PROJECT_DIR}/report/$output_prefix/$output_prefix.Pharmacogenomics.xls\nfi\n\n";
		#print SHELL "\tless $sample.mp.vcf | $CONFIG{'BIOVCF'} --skip-header --eval '[header.samples,r.chrom,r.pos,r.ref,r.alt.join(\",\")]' --seval '[s.dp, s.ad[0], s.ad[1]]' > $sample.mp.tab\nfi\n\n";
		#print SHELL "\tless $sample.mp.vcf | $CONFIG{'BIOVCF'} --skip-header --eval '[header.samples,r.chrom,r.pos,r.ref,r.alt.join(\",\")]' --seval '[s.dp, s.ad[0], s.ad[1]]' | awk -v s=$sample '{print s\"\\t\"\$0}' |awk 'BEGIN{OFS=\"\\t\";}{if(\$5==\".\"&&\$7<2&&\$8<2){\$7=\$6;print;}else{print;}}' > $sample.mp.tab\nfi\n\n";
		
		print SHELL "#4#5 500M cfDNA_info\n";
		if($sample !~ /QYX/)
		{
			# ddcfDNA quantification
			print SHELL "export LD_LIBRARY_PATH=/usr/local/lib64:/usr/local/lib:/data/soft/lib:/usr/local/lib64:/usr/local/lib:/data/soft/lib:/opt/gridengine/lib/lx-amd64:/opt/openmpi/lib\n";
			print SHELL "sh $CONFIG{'SOT_SCRIPT_DIR'}/tabex.sh $sample\n";

			#print SHELL "/data/soft/R-3.4.3/bin/Rscript $CONFIG{'SOT_SCRIPT_DIR'}/ddcfDNA_Percent.R $sample.mp.tab $sample\n";
			#print SHELL "/data/soft/R-3.4.3/bin/Rscript /data/Develop/ddcfDNA_quantification/ddcfDNA_Percent.ML.Sampling.v2.R $sample.mp.tab.ex $sample\n";
			print SHELL "cat $sample.ddcfDNA_percent.xls >> $CONFIG{'SOT_PROJECT_DIR'}/report/$output_prefix/$output_prefix.ddcfDNA_percent.xls\n";
			#print SHELL "if [ ! -f $CONFIG{'SOT_PROJECT_DIR'}/cfDNA_length/$sample.cfDNA_length.png ]\nthen\n";
			#print SHELL "\tsh $CONFIG{'SOT_SCRIPT_DIR'}/extract_organ_source_read.sh \$(pwd)\n\tmv *png $CONFIG{'SOT_PROJECT_DIR'}/cfDNA_length\nfi\n";
			#print SHELL "if [ ! -f $sample.sort.bam.IS ]\nthen\n";
			#print SHELL "\tsh /data/pipeline/NGS_Target_Sequencing/Organ_Transplant/IS.sh $sample.sort.bam\nfi\n";
			#print SHELL "if [ ! -f $sample.UT.sscs.sort.bam.IS ]\nthen\n";
			#print SHELL "\tsh /data/pipeline/NGS_Target_Sequencing/Organ_Transplant/IS.sh $sample.UT.sscs.sort.bam\nfi\n";
			print SHELL "if [[ -f $sample.mp.tab.xd-cfDNA.insertsize ]]\nthen\n";
			print SHELL "\t/data/soft/bin/Rscript /data/pipeline/NGS_Target_Sequencing/Organ_Transplant/IS.sh.R $sample\nfi\n";
			print SHELL "cp *pdf $CONFIG{'SOT_PROJECT_DIR'}/report/$output_prefix/\n\n";
		}

		print SHELL "#4#6 500M QC\n";
		print SHELL "export LD_LIBRARY_PATH=/usr/local/lib64:/usr/local/lib:/data/soft/lib:/usr/local/lib64:/usr/local/lib:/data/soft/lib:/opt/gridengine/lib/lx-amd64:/opt/openmpi/lib\n";
		print SHELL "if [ ! -f $sample.QC ]\nthen\n";
		print SHELL "\tsh /data/pipeline/OT.sh $sample.sort.bam $sample.UT.sort.bam $sample.UT.sscs.sort.bam $CONFIG{'BED_SOT'} > $sample.QC\nfi\n";
		print SHELL "if [ ! -f QC/$sample.QC.xls ]\nthen\n";
		print SHELL "\tsh $CONFIG{'SOT_SCRIPT_DIR'}/QC.sh $sample.withtag.R1.fq.gz $sample.withtag.R2.fq.gz $sample\nfi\n\n";
	}
	close SHELL;

	print STDERR "$sample.script.sh created\n\n";
#	system("perl /data/pipeline/NGS_Target_Sequencing/ngs_pip.pl --runscript $sample.script.sh") if($run);
}

sub shell_for_SOT
{
	my ($sample,$targetsp,$output_prefix,$Ncores,$skipIFD)=@_;

	my $type = "";

	unless($sample)
	{
		print STDERR "program <sample> <experiment type>[HYB]\n";
		print STDERR "Initial failed. Exiting...\n";
		exit;
	}
	my @R1R2 = @{$OSS_FILEs{$sample}};
	my @a = split("/",$R1R2[0]);
	my @b = split("/",$R1R2[1]);

	my $dir = join("/",@a[3..$#a-1]);
	my $R1 = $a[-1];
	my $R2 = $b[-1];

	chdir($CONFIG{'SOT_PROJECT_DIR'});
	unless( -s "$CONFIG{SOT_PROJECT_DIR}/$TASKS{$sample}")
	{
		mkdir("$CONFIG{SOT_PROJECT_DIR}/$TASKS{$sample}",0755) or die "mkdir $CONFIG{SOT_PROJECT_DIR}/$TASKS{$sample} Failed";
	}
	unless( -s "$CONFIG{SOT_PROJECT_DIR}/$TASKS{$sample}/$targetsp")
	{
		mkdir("$CONFIG{SOT_PROJECT_DIR}/$TASKS{$sample}/$targetsp", 0755) or die "mkdir $CONFIG{SOT_PROJECT_DIR}/$TASKS{$sample}/$targetsp Failed";
	}
	unless( -s "$CONFIG{SOT_PROJECT_DIR}/SHELL/$output_prefix")
	{
		mkdir("$CONFIG{SOT_PROJECT_DIR}/SHELL/$output_prefix", 0755) or die "mkdir $CONFIG{SOT_PROJECT_DIR}/SHELL/$output_prefix Failed";
	}

	chdir("$TASKS{$sample}/$targetsp");

	unless( -s "$targetsp\_combined_R1.fastq.gz" )
	{
		system("ln -s $CONFIG{'SOT_RAWDATA_DIR'}/$dir/$R1 $targetsp\_combined_R1.fastq.gz");
	}
	unless( -s "$targetsp\_combined_R2.fastq.gz" )
	{
		system("ln -s $CONFIG{'SOT_RAWDATA_DIR'}/$dir/$R2 $targetsp\_combined_R2.fastq.gz");
	}

	if($targetsp=~/WGS/ || $targetsp=~/IF/)
	{
		$type = "WGS";

	}elsif($targetsp=~/QY\d+[HQG]?/)
	{
		$type = "HYB";
	}elsif($targetsp=~/QY[XN]\d+/)
	{
		$type = "HYB";
	}elsif($targetsp=~/S\dQY[RNWB]\d+H\d+/)
	{
		$type = "WGS";#only analysis infection
	}elsif($targetsp=~/QY\d+G/)
	{
		$type = "HYB";
	}

	print "ProType\t$sample==>$targetsp\t$type\n";

	#my $r1 = $R1;
	#my $r2 = $R2;
	my $r1 = "$targetsp\_combined_R1.fastq.gz";
	my $r2 = "$targetsp\_combined_R2.fastq.gz";

	open SHELL, ">$CONFIG{SOT_PROJECT_DIR}/$TASKS{$sample}/$targetsp/$targetsp.script.sh";
	if($type eq "WGS")
	{
		my $r1fastp="$targetsp\_fastp_R1.fastq.gz";
		my $r2fastp="$targetsp\_fastp_R2.fastq.gz";

		print SHELL "#1 3G filter\n";
		print SHELL "$CONFIG{FASTP} -i $r1 -I $r2 -o $r1fastp -O $r2fastp -j $targetsp.json -h $targetsp.html  -l 31 -w 16 -Q\n\n";

		print SHELL "#1#2 40G kraken\n";
		print SHELL "export KRAKEN_DB_PATH=/data/soft/kraken/\n";
		print SHELL "if [ ! -f $targetsp.kraken ]\nthen\n";
		print SHELL "\t$CONFIG{KRAKEN} -db $CONFIG{KRAKEN_DB} --threads 16 --fastq-input --gzip-compressed --output $targetsp.kraken --preload $r1fastp $r2fastp\nfi\n\n";

		print SHELL "#2#3 500M CDCReport\n";
		print SHELL "awk 'NR==FNR{a[\$1]=\$1;}NR!=FNR{if(a[\$3]){print;}}' $CONFIG{'CDC_TAXID'} $targetsp.kraken >$targetsp.kraken.cdc\n";
		print SHELL "$CONFIG{'KRAKEN-FILTER'} --db $CONFIG{'KRAKEN_DB'} --threshold 0.3 $targetsp.kraken.cdc > $targetsp.kraken.cdc.2\n";
		print SHELL "#$CONFIG{'KRAKEN-REPORT'} --db $CONFIG{'KRAKEN_DB'} $targetsp.kraken.cdc |sort > $targetsp.kraken.cdc.report.txt\n";
		print SHELL "#$CONFIG{'KRAKEN-MPA-REPORT'} --db $CONFIG{'KRAKEN_DB'}  $targetsp.kraken.cdc | sort > $targetsp.kraken.cdc.mpa.txt\n";
		print SHELL "$CONFIG{'KRAKEN-REPORT'} --db $CONFIG{'KRAKEN_DB'} $targetsp.kraken.cdc.2 | sort > $targetsp.kraken.cdc.2.report.txt\n";
		print SHELL "$CONFIG{'KRAKEN-MPA-REPORT'} --db $CONFIG{'KRAKEN_DB'}  $targetsp.kraken.cdc.2 | sort > $targetsp.kraken.cdc.2.mpa.txt\n\n";

		print SHELL "#3#4 5G BWAalign\n";
		print SHELL "$CONFIG{'PATHOGEN_REF_READ'} $targetsp.kraken.cdc $targetsp $targetsp\_fastp_R1.fastq.gz $targetsp\_fastp_R2.fastq.gz $CONFIG{'KRAKEN_PATH'}/$CONFIG{'KRAKEN_DB'} \n";
		#print SHELL "$CONFIG{'BWA'} index $targetsp.aligned_cdc.fa\n";
		print SHELL "$CONFIG{'BWA'} mem $CONFIG{'CDC_TAXDB'} $targetsp\_CDC_R1.fastq $targetsp\_CDC_R2.fastq -t $Ncores  -R \"\@RG:\\tID:$targetsp\\tSM:$targetsp\" -v 1 \| $CONFIG{'SAMTOOLS'} sort -T tmp -o $targetsp.cdc.sort.bam -\n";
		print SHELL "$CONFIG{'PICARD'} MarkDuplicates I=$targetsp.cdc.sort.bam O=$targetsp.cdc.picard.sort.bam R=$CONFIG{'CDC_TAXDB'} VALIDATION_STRINGENCY=SILENT QUIET=true REMOVE_DUPLICATES=true M=$targetsp.picard.txt\n";

		#print "Script made.  Run with:\n$targetsp.script.sh\n\n";

	}elsif($type eq "HYB")
	{
		#set parameter
		my ($readLength,$barcodeLength,$repFilt)=cal_barcode_read_len($r1,$CONFIG{'SPACER'});

		print SHELL "#1 3G filter\n";
		print SHELL "barcodeLength=$barcodeLength\n";
		print SHELL "if [ ! -f $targetsp.withtag.R1.fq.gz ]\nthen\n";
		print SHELL "\tif [ \$barcodeLength -gt 2 ]\nthen\n";
		print SHELL "\t\tperl $CONFIG{tag_to_header} $r1 $r2 $targetsp\n\telse\n";
		print SHELL "\t\tln -s $r1 $targetsp.withtag.R1.fq.gz\n\t\tln -s $r2 $targetsp.withtag.R2.fq.gz\n\tfi\nfi\n\n";
		unless($skipIFD)
		{
			print SHELL "#1#2 40G kraken\n";
			print SHELL "export KRAKEN_DB_PATH=/data/soft/kraken/\n";
			print SHELL "if [ ! -f $targetsp.kraken ]\nthen\n";
			print SHELL "$CONFIG{'KRAKEN'} -db $CONFIG{'KRAKEN_DB'} --threads $CONFIG{'Ncores'} --fastq-input --output $targetsp.kraken --paired --preload --gzip-compressed $targetsp.withtag.R1.fq.gz $targetsp.withtag.R2.fq.gz\n";
			print SHELL "awk \'NR==FNR{a[\$1]=\$1;}NR!=FNR{if(a[\$3]){print;}}\' $CONFIG{'CDC_TAXID'} $targetsp.kraken >$targetsp.kraken.cdc\n";
			print SHELL "$CONFIG{'KRAKEN-FILTER'} --db $CONFIG{'KRAKEN_DB'} --threshold 0.3 $targetsp.kraken.cdc  > $targetsp.kraken.cdc.2\n";
			print SHELL "$CONFIG{'KRAKEN-REPORT'} --db $CONFIG{'KRAKEN_DB'} $targetsp.kraken.cdc.2 | sort > $targetsp.kraken.cdc.2.report.txt\nfi\n\n";


			print SHELL "#2#3 20G CDCalign\n";

			unless($targetsp =~ /QYX/)
			{
				print SHELL "if [ ! -f $targetsp\_CDC_R1.fastq ]\nthen\n";
				print SHELL "\t$CONFIG{'PATHOGEN_REF_READ'} $targetsp.kraken.cdc $targetsp $targetsp.withtag.R1.fq.gz $targetsp.withtag.R2.fq.gz $CONFIG{'KRAKEN_PATH'}/$CONFIG{'KRAKEN_DEFAULT_DB'}\nfi\n";
				print SHELL "if [ ! -f $targetsp.cdc.sort.bam ]\nthen\n";
				print SHELL "\t$CONFIG{'BWA'} mem $CONFIG{'CDC_TAXDB'} $targetsp\_CDC_R1.fastq $targetsp\_CDC_R2.fastq -t $Ncores -R \"\@RG\\tID:$targetsp\\tSM:$targetsp\" -v 1 | $CONFIG{'SAMTOOLS'} sort -T $targetsp.cdc.tmp -o $targetsp.cdc.sort.bam -\nfi\n";
				print SHELL "if [ ! -f $targetsp.cdc.sscs.sort.bam ]\nthen\n";
				if($barcodeLength >= 2)
				{
					print SHELL "\texport PYTHONPATH=/home/lxf/.PY2.7LIB/lib64/python2.7/site-packages/\n";
					print SHELL "\t/usr/bin/python $CONFIG{'DSpath'}/ConsensusMaker.v2.py --infile $targetsp.cdc.sort.bam --outfile $targetsp.cdc.sscs.bam --minmem $CONFIG{'minMem'} --maxmem $CONFIG{'maxMem'} --read_length $readLength --cut_off $CONFIG{'cutOff'} --Ncut_off $CONFIG{'nCutOff'} --read_type $CONFIG{'readTypes'} --filt $CONFIG{'filtersSet'} --isize $CONFIG{'iSize'} --rep_filt $repFilt\n";
					print SHELL "$CONFIG{'SAMTOOLS'} view -h $targetsp.cdc.sscs.bam|awk '{if (\$6!~/^[[:digit:]]+N[[:digit:]]+D|^[[:digit:]]+P|[[:digit:]]+P\$|^[[:digit:]]+I/) print \$0}' | $CONFIG{'SAMTOOLS'} sort -T tmp -o $targetsp.cdc.sscs.sort.bam -\@ 8 -\nfi\n\n";
				}else{
					print SHELL "\t$CONFIG{'PICARD'} MarkDuplicates I=$targetsp.cdc.sort.bam O=$targetsp.cdc.sort.picard.bam M=$targetsp.cdc.picard.txt R=$CONFIG{'hg38'} VALIDATION_STRINGENCY=SILENT QUIET=true REMOVE_DUPLICATES=true\n";
					print SHELL "\tif [ ! -f $targetsp.cdc.sscs.sort.bam ]\n\tthen\n";
					print SHELL "\t\tln -s $targetsp.cdc.sort.picard.bam $targetsp.cdc.sscs.sort.bam\n\tfi\nfi\n";
				}
			}
		}else{
			#print "echo non-neccessary CDC Aligning on $targetsp\n";
		}

		#======================XX=========================
		#======================XY=========================
		#genotyping
		print SHELL "#1#4 5G genotyping\n";
		print SHELL "if [ ! -f $targetsp.sort.bam ]\nthen\n";
		print SHELL "\t$CONFIG{'BWA'} mem $CONFIG{'hg38'}  $targetsp.withtag.R1.fq.gz $targetsp.withtag.R2.fq.gz -t $Ncores -R \"\@RG\\tID:$targetsp\\tSM:$targetsp\" -v 1 | $CONFIG{'SAMTOOLS'} sort -T $targetsp.tmp -o $targetsp.sort.bam -\nfi\n";
		print SHELL "if [ ! -f $targetsp.sort.bam.bai ]\nthen\n\t$CONFIG{'SAMBAMBA'}  index $targetsp.sort.bam\nfi\n";
		print SHELL "if [ ! -f $targetsp.UT.sort.bam ]\nthen\n\t$CONFIG{'SAMBAMBA'} view -t 12 -L $CONFIG{'BED_SOT_MTB'} -h -f bam -F \"mapping_quality >= 30 and not (unmapped or secondary_alignment) and not ([XA] != null or [SA] != null) and ((not [MD]=~/\\^/ and [NM]<=2) or (cigar=~/[ID]/ and [NM]<8))\" $targetsp.sort.bam -o $targetsp.UT.sort.bam\nfi\n";
		print SHELL "if [ ! -f $targetsp.UT.sscs.sort.bam ]\nthen\n";
		if($barcodeLength >= 2)
		{
			print SHELL "\texport PYTHONPATH=/home/lxf/.PY2.7LIB/lib64/python2.7/site-packages/\n";
			print SHELL "\t/usr/bin/python $CONFIG{'DSpath'}/ConsensusMaker.v2.py --infile $targetsp.UT.sort.bam --outfile $targetsp.UT.sscs.bam --minmem $CONFIG{'minMem'} --maxmem $CONFIG{'maxMem'} --read_length $readLength --cut_off $CONFIG{'cutOff'} --Ncut_off $CONFIG{'nCutOff'} --read_type $CONFIG{'readTypes'} --filt $CONFIG{'filtersSet'} --isize $CONFIG{'iSize'} --rep_filt $repFilt\n";
			print SHELL "\t$CONFIG{'SAMTOOLS'} view -h $targetsp.UT.sscs.bam|awk '{if (\$6!~/^[[:digit:]]+N[[:digit:]]+D|^[[:digit:]]+P|[[:digit:]]+P\$|^[[:digit:]]+I/) print \$0}' | $CONFIG{'SAMTOOLS'} sort -T $targetsp.tmp -o $targetsp.UT.sscs.sort.bam -\@ 8 -\nfi\n";
		}else{
			print SHELL "\t$CONFIG{'PICARD'} MarkDuplicates I=$targetsp.UT.sort.bam O=$targetsp.UT.sort.picard.bam M=$targetsp.picard.txt R=$CONFIG{'hg38'} VALIDATION_STRINGENCY=SILENT QUIET=true REMOVE_DUPLICATES=true\n";
			print SHELL "\tif [ ! -f $targetsp.UT.sscs.sort.bam ]\n\tthen\n";
			print SHELL "\t\tln -s $targetsp.UT.sort.picard.bam $targetsp.UT.sscs.sort.bam\n\tfi\nfi\n";
		}
		#call variant
		#print SHELL "if [ ! -f $targetsp.sscs.sort.bam.bai ]\nthen\n";
		#print SHELL "\t$CONFIG{'SAMTOOLS'} index $targetsp.sscs.sort.bam\nfi\n";
		print SHELL "if [ ! -f $targetsp.mp.vcf ]\nthen\n";
		print SHELL "\t$CONFIG{'SAMTOOLS'} mpileup -A -uvf $CONFIG{'hg38'} -l $CONFIG{'BED_SOT_MTB'} -t DP,AD $targetsp.UT.sscs.sort.bam | $CONFIG{'BCFTOOLS'} call --multiallelic-caller --keep-alts --targets-file $CONFIG{'BED_SOT_MTB'} > $targetsp.mp.vcf\nfi\n";
		print SHELL "if [ ! -f $targetsp.mp.tab ]\nthen\n";
		print SHELL "\tperl $CONFIG{'VCF2TAB'} $targetsp.mp.vcf $CONFIG{'BED_SOT'} $targetsp.mp.tab\n";
		print SHELL "\tperl $CONFIG{'VCF2TAB'} $targetsp.mp.vcf $CONFIG{'BED_MTB'} $targetsp.mp.mtb.tab\n";
		print SHELL "\tcp $targetsp.mp.mtb.tab $CONFIG{SOT_PROJECT_DIR}/report/$output_prefix/$targetsp.Pharmacogenomics.xls\n";
		print SHELL "\tcat $targetsp.mp.mtb.tab >> $CONFIG{SOT_PROJECT_DIR}/report/$output_prefix/$output_prefix.Pharmacogenomics.xls\nfi\n\n";
		#print SHELL "\tless $targetsp.mp.vcf | $CONFIG{'BIOVCF'} --skip-header --eval '[header.targetsps,r.chrom,r.pos,r.ref,r.alt.join(\",\")]' --seval '[s.dp, s.ad[0], s.ad[1]]' > $targetsp.mp.tab\nfi\n\n";
		#print SHELL "\tless $targetsp.mp.vcf | $CONFIG{'BIOVCF'} --skip-header --eval '[header.targetsps,r.chrom,r.pos,r.ref,r.alt.join(\",\")]' --seval '[s.dp, s.ad[0], s.ad[1]]' | awk -v s=$targetsp '{print s\"\\t\"\$0}' |awk 'BEGIN{OFS=\"\\t\";}{if(\$5==\".\"&&\$7<2&&\$8<2){\$7=\$6;print;}else{print;}}' > $targetsp.mp.tab\nfi\n\n";

		print SHELL "#4#5 500M cfDNA_info\n";
		if($targetsp !~ /QYX/)
		{
			# ddcfDNA quantification
			print SHELL "export LD_LIBRARY_PATH=/usr/local/lib64:/usr/local/lib:/data/soft/lib:/usr/local/lib64:/usr/local/lib:/data/soft/lib:/opt/gridengine/lib/lx-amd64:/opt/openmpi/lib\n";
			print SHELL "sh $CONFIG{'SOT_SCRIPT_DIR'}/tabex.sh $targetsp\n";

			#print SHELL "/data/soft/R-3.4.3/bin/Rscript $CONFIG{'SOT_SCRIPT_DIR'}/ddcfDNA_Percent.R $targetsp.mp.tab $targetsp\n";
			#print SHELL "/data/soft/R-3.4.3/bin/Rscript /data/Develop/ddcfDNA_quantification/ddcfDNA_Percent.ML.Sampling.v2.R $targetsp.mp.tab.ex $targetsp\n";
			print SHELL "cat $targetsp.ddcfDNA_percent.xls >> $CONFIG{'SOT_PROJECT_DIR'}/report/$output_prefix/$output_prefix.ddcfDNA_percent.xls\n";
			#print SHELL "if [ ! -f $CONFIG{'SOT_PROJECT_DIR'}/cfDNA_length/$targetsp.cfDNA_length.png ]\nthen\n";
			#print SHELL "\tsh $CONFIG{'SOT_SCRIPT_DIR'}/extract_organ_source_read.sh \$(pwd)\n\tmv *png $CONFIG{'SOT_PROJECT_DIR'}/cfDNA_length\nfi\n";
			#print SHELL "if [ ! -f $targetsp.sort.bam.IS ]\nthen\n";
			#print SHELL "\tsh /data/pipeline/NGS_Target_Sequencing/Organ_Transplant/IS.sh $targetsp.sort.bam\nfi\n";
			#print SHELL "if [ ! -f $targetsp.UT.sscs.sort.bam.IS ]\nthen\n";
			#print SHELL "\tsh /data/pipeline/NGS_Target_Sequencing/Organ_Transplant/IS.sh $targetsp.UT.sscs.sort.bam\nfi\n";
			print SHELL "if [[ -f $targetsp.mp.tab.xd-cfDNA.insertsize ]]\nthen\n";
			print SHELL "\t/data/soft/bin/Rscript /data/pipeline/NGS_Target_Sequencing/Organ_Transplant/IS.sh.R $targetsp\nfi\n";
			print SHELL "cp *pdf $CONFIG{'SOT_PROJECT_DIR'}/report/$output_prefix/\n\n";

		}
		print SHELL "#4#6 500M QC\n";
		print SHELL "export LD_LIBRARY_PATH=/usr/local/lib64:/usr/local/lib:/data/soft/lib:/usr/local/lib64:/usr/local/lib:/data/soft/lib:/opt/gridengine/lib/lx-amd64:/opt/openmpi/lib\n";
		print SHELL "if [ ! -f $targetsp.QC ]\nthen\n";
		print SHELL "\tsh /data/pipeline/OT.sh $targetsp.sort.bam $targetsp.UT.sort.bam $targetsp.UT.sscs.sort.bam $CONFIG{'BED_SOT'} > $targetsp.QC\nfi\n";
		print SHELL "if [ ! -f QC/$targetsp.QC.xls ]\nthen\n";
		print SHELL "\tsh $CONFIG{'SOT_SCRIPT_DIR'}/QC.sh $targetsp.withtag.R1.fq.gz $targetsp.withtag.R2.fq.gz $targetsp\nfi\n\n";
	}
	close SHELL;

	print STDERR "$targetsp.script.sh created\n";
	system("cp $targetsp.script.sh $CONFIG{SOT_PROJECT_DIR}/SHELL/$output_prefix");
	#report list
	open REPORTLIST,">>$CONFIG{SOT_PROJECT_DIR}/report/$output_prefix/$output_prefix.Pathogen.list" or die $!;
	print REPORTLIST "$CONFIG{SOT_PROJECT_DIR}/$TASKS{$sample}/$targetsp/$targetsp.kraken.cdc.2.report.txt\n";
	close REPORTLIST;
	$STATUS{$targetsp}{'INITIAL'} = "DONE";
}

sub shell_for_QZ
{
	my ($sample,$output_prefix,$Ncores,$skipIFD)=@_;

	my $type = "";

	unless($sample)
	{
		print STDERR "program <sample> <experiment type>[HYB]\n";
		print STDERR "Initial failed. Exiting...\n";
		exit;
	}
	my @R1R2 = @{$OSS_FILEs{$sample}};
	my @a = split("/",$R1R2[0]);
	my $dir = join("/",@a[3..$#a-1]);
	my $R1 = $a[-1];
	my $R2 = $R1; $R2 =~ s/R1/R2/;

	chdir($CONFIG{'QZ_PROJECT_DIR'});
	unless( -s $TASKS{$sample})
	{
		mkdir($TASKS{$sample},0755) or die "mkdir $TASKS{$sample} Failed in $CONFIG{'QZ_PROJECT_DIR'}";
	}
	unless( -s "$TASKS{$sample}/$sample")
	{
		mkdir("$TASKS{$sample}/$sample",0755)  or die "mkdir $TASKS{$sample}/$sample Failed in $CONFIG{'QZ_PROJECT_DIR'}";
	}
	chdir("$TASKS{$sample}/$sample");
	unless( -s $R1 )
	{
		system("ln -s $CONFIG{'QZ_RAWDATA_DIR'}/$dir/$R1 .");
	}
	unless( -s $R2 )
	{
		system("ln -s $CONFIG{'QZ_RAWDATA_DIR'}/$dir/$R2 .");
	}
	my $r1 = $R1;
	my $r2 = $R2;

	#set parameter
	my ($readLength,$barcodeLength,$repFilt)=cal_barcode_read_len($r1,$CONFIG{'SPACER'});
	open SHELL, ">$CONFIG{QZ_PROJECT_DIR}/$TASKS{$sample}/$sample/$sample.script.sh";		
	print SHELL "#1 3G filter\n";
	print SHELL "barcodeLength=$barcodeLength\n";
	print SHELL "if [ ! -f $sample.withtag.R1.fq.gz ]\nthen\n";
	print SHELL "\tif [ \$barcodeLength -gt 2 ]\nthen\n";
	print SHELL "\t\tperl $CONFIG{tag_to_header} $r1 $r2 $sample\n\telse\n";
	print SHELL "\t\tln -s $r1 $sample.withtag.R1.fq.gz\n\t\tln -s $r2 $sample.withtag.R2.fq.gz\n\tfi\nfi\n\n";

	print SHELL "#1#2 40G kraken\n";
	print SHELL "export KRAKEN_DB_PATH=/data/soft/kraken/\n";
	print SHELL "if [ ! -f $sample.kraken ]\nthen\n";
	print SHELL "$CONFIG{'KRAKEN'} -db $CONFIG{'KRAKEN_DB'} --threads $CONFIG{'Ncores'} --fastq-input --output $sample.kraken --paired --preload --gzip-compressed $sample.withtag.R1.fq.gz $sample.withtag.R2.fq.gz\n";
	print SHELL "awk \'NR==FNR{a[\$1]=\$1;}NR!=FNR{if(a[\$3]){print;}}\' $CONFIG{'CDC_TAXID'} $sample.kraken >$sample.kraken.cdc\n";
	print SHELL "$CONFIG{'KRAKEN-FILTER'} --db $CONFIG{'KRAKEN_DB'} --threshold 0.3 $sample.kraken.cdc  > $sample.kraken.cdc.2\n";
	print SHELL "$CONFIG{'KRAKEN-REPORT'} --db $CONFIG{'KRAKEN_DB'} $sample.kraken.cdc.2 | sort > $sample.kraken.cdc.2.report.txt\nfi\n\n";


	print SHELL "#2#3 20G CDCalign\n";
	print SHELL "if [ ! -f $sample\_CDC_R1.fastq ]\nthen\n";
	print SHELL "\t$CONFIG{'PATHOGEN_REF_READ'} $sample.kraken.cdc $sample $sample.withtag.R1.fq.gz $sample.withtag.R2.fq.gz $CONFIG{'KRAKEN_PATH'}/$CONFIG{'KRAKEN_DEFAULT_DB'}\nfi\n";
	print SHELL "if [ ! -f $sample.cdc.sort.bam ]\nthen\n";
	print SHELL "\t$CONFIG{'BWA'} mem $CONFIG{'CDC_TAXDB'} $sample\_CDC_R1.fastq $sample\_CDC_R2.fastq -t $Ncores -R \"\@RG\\tID:$sample\\tSM:$sample\" -v 1 | $CONFIG{'SAMTOOLS'} sort -T $sample.cdc.tmp -o $sample.cdc.sort.bam -\nfi\n";

	print SHELL "if [ ! -f $sample.cdc.sscs.sort.bam ]\nthen\n";
	if($barcodeLength >= 2)
	{
		print SHELL "\texport PYTHONPATH=/home/lxf/.PY2.7LIB/lib64/python2.7/site-packages/\n";
		print SHELL "\t/usr/bin/python $CONFIG{'DSpath'}/ConsensusMaker.v2.py --infile $sample.cdc.sort.bam --outfile $sample.cdc.sscs.bam --minmem $CONFIG{'minMem'} --maxmem $CONFIG{'maxMem'} --read_length $readLength --cut_off $CONFIG{'cutOff'} --Ncut_off $CONFIG{'nCutOff'} --read_type $CONFIG{'readTypes'} --filt $CONFIG{'filtersSet'} --isize $CONFIG{'iSize'} --rep_filt $repFilt\n";
		print SHELL "$CONFIG{'SAMTOOLS'} view -h $sample.cdc.sscs.bam|awk '{if (\$6!~/^[[:digit:]]+N[[:digit:]]+D|^[[:digit:]]+P|[[:digit:]]+P\$|^[[:digit:]]+I/) print \$0}' | $CONFIG{'SAMTOOLS'} sort -T tmp -o $sample.cdc.sscs.sort.bam -\@ 8 -\nfi\n\n";
	}else{
		print SHELL "\t$CONFIG{'PICARD'} MarkDuplicates I=$sample.cdc.sort.bam O=$sample.cdc.sort.picard.bam M=$sample.cdc.picard.txt R=$CONFIG{'hg38'} VALIDATION_STRINGENCY=SILENT QUIET=true REMOVE_DUPLICATES=true\n";
		print SHELL "ln -s $sample.cdc.sort.picard.bam $sample.cdc.sscs.sort.bam\nfi\n\n";
	}

	#======================XX=========================
	#======================XY=========================
	#genotyping
	print SHELL "#1#4 5G genotyping\n";
	print SHELL "if [ ! -f $sample.sort.bam ]\nthen\n";
	print SHELL "\t$CONFIG{'BWA'} mem $CONFIG{'hg38'}  $sample.withtag.R1.fq.gz $sample.withtag.R2.fq.gz -t $Ncores -R \"\@RG\\tID:$sample\\tSM:$sample\" -v 1 | $CONFIG{'SAMTOOLS'} sort -T $sample.tmp -o $sample.sort.bam -\nfi\n";
	print SHELL "if [ ! -f $sample.sort.bam.bai ]\nthen\n\t$CONFIG{'SAMBAMBA'}  index $sample.sort.bam\nfi\n";
	print SHELL "if [ ! -f $sample.UT.sort.bam ]\nthen\n\t$CONFIG{'SAMBAMBA'} view -t 12 -L $CONFIG{'BED_SOT_MTB'} -h -f bam -F \"mapping_quality >= 30 and not (unmapped or secondary_alignment) and not ([XA] != null or [SA] != null) and ((not [MD]=~/\\^/ and [NM]<=2) or (cigar=~/[ID]/ and [NM]<8))\" $sample.sort.bam -o $sample.UT.sort.bam\nfi\n";
	print SHELL "if [ ! -f $sample.UT.sscs.bam ]\nthen\n";
	if($barcodeLength >= 2)
	{
		print SHELL "\texport PYTHONPATH=/home/lxf/.PY2.7LIB/lib64/python2.7/site-packages/\n";
		print SHELL "\t/usr/bin/python $CONFIG{'DSpath'}/ConsensusMaker.v2.py --infile $sample.UT.sort.bam --outfile $sample.UT.sscs.bam --minmem $CONFIG{'minMem'} --maxmem $CONFIG{'maxMem'} --read_length $readLength --cut_off $CONFIG{'cutOff'} --Ncut_off $CONFIG{'nCutOff'} --read_type $CONFIG{'readTypes'} --filt $CONFIG{'filtersSet'} --isize $CONFIG{'iSize'} --rep_filt $repFilt\n";
		print SHELL "\t$CONFIG{'SAMTOOLS'} view -h $sample.UT.sscs.bam|awk '{if (\$6!~/^[[:digit:]]+N[[:digit:]]+D|^[[:digit:]]+P|[[:digit:]]+P\$|^[[:digit:]]+I/) print \$0}' | $CONFIG{'SAMTOOLS'} sort -T tmp -o $sample.UT.sscs.sort.bam -\@ 8 -\nfi\n";
	}else{
		print SHELL "\t$CONFIG{'PICARD'} MarkDuplicates I=$sample.UT.sort.bam O=$sample.UT.sort.picard.bam M=$sampl.UT.picard.txt R=$CONFIG{'hg38'} VALIDATION_STRINGENCY=SILENT QUIET=true REMOVE_DUPLICATES=true\n";
		print SHELL "\tif [ ! -f $sample.UT.sscs.sort.bam ]\n\tthen\n";
		print SHELL "\t\tln -s $sample.UT.sort.picard.bam $sample.UT.sscs.sort.bam\n\tfi\nfi\n";
	}

	#call variant
	#print SHELL "if [ ! -f $sample.sscs.sort.bam.bai ]\nthen\n";
	#print SHELL "\t$CONFIG{'SAMTOOLS'} index $sample.sscs.sort.bam\nfi\n";
	print SHELL "if [ ! -f $sample.mp.vcf ]\nthen\n";
	print SHELL "\t$CONFIG{'SAMTOOLS'} mpileup -A -uvf $CONFIG{'hg38'} -l $CONFIG{'BED_SOT'} -t DP,AD $sample.UT.sscs.sort.bam | $CONFIG{'BCFTOOLS'} call --multiallelic-caller --keep-alts --targets-file $CONFIG{'BED_SOT'} > $sample.mp.vcf\nfi\n";
	print SHELL "if [ ! -f $sample.mp.tab ]\nthen\n";
	print SHELL "\tperl $CONFIG{'VCF2TAB'} $sample.mp.vcf $CONFIG{'BED_SOT'} $sample.mp.tab\n";
	print SHELL "\tperl $CONFIG{'VCF2TAB'} $sample.mp.vcf $CONFIG{'BED_MTB'} $sample.mp.mtb.tab\n";
	print SHELL "\tcp $sample.mp.mtb.tab $CONFIG{SOT_PROJECT_DIR}/report/$output_prefix/$sample.Pharmacogenomics.xls\n";
	print SHELL "\tcat $sample.mp.mtb.tab >> $CONFIG{SOT_PROJECT_DIR}/report/$output_prefix/$output_prefix.Pharmacogenomics.xls\nfi\n\n";


	#print SHELL "\tless $sample.mp.vcf | $CONFIG{'BIOVCF'} --skip-header --eval '[header.samples,r.chrom,r.pos,r.ref,r.alt.join(\",\")]' --seval '[s.dp, s.ad[0], s.ad[1]]' > $sample.mp.tab\nfi\n\n";
	#print SHELL "\tless $sample.mp.vcf | $CONFIG{'BIOVCF'} --skip-header --eval '[header.samples,r.chrom,r.pos,r.ref,r.alt.join(\",\")]' --seval '[s.dp, s.ad[0], s.ad[1]]' | awk -v s=$sample '{print s\"\\t\"\$0}' |awk 'BEGIN{OFS=\"\\t\";}{if(\$5==\".\"&&\$7<2&&\$8<2){\$7=\$6;print;}else{print;}}' > $sample.mp.tab\nfi\n\n";

	print SHELL "#4#5 500M cfDNA_info\n";
	if($sample !~ /QYX/)
	{
		# ddcfDNA quantification
		print SHELL "export LD_LIBRARY_PATH=/usr/local/lib64:/usr/local/lib:/data/soft/lib:/usr/local/lib64:/usr/local/lib:/data/soft/lib:/opt/gridengine/lib/lx-amd64:/opt/openmpi/lib\n";
		print SHELL "sh $CONFIG{'SOT_SCRIPT_DIR'}/tabex.sh $sample\n";

		#print SHELL "/data/soft/R-3.4.3/bin/Rscript $CONFIG{'SOT_SCRIPT_DIR'}/ddcfDNA_Percent.R $sample.mp.tab $sample\n";
		#print SHELL "/data/soft/R-3.4.3/bin/Rscript /data/Develop/ddcfDNA_quantification/ddcfDNA_Percent.ML.Sampling.v2.R $sample.mp.tab.ex $sample\n";
		print SHELL "cat $sample.ddcfDNA_percent.xls >> $CONFIG{SOT_PROJECT_DIR}/report/$output_prefix/$output_prefix.ddcfDNA_percent.xls\n";
		#print SHELL "if [ ! -f $CONFIG{'SOT_PROJECT_DIR'}/cfDNA_length/$sample.cfDNA_length.png ]\nthen\n";
		#print SHELL "\tsh $CONFIG{'SOT_SCRIPT_DIR'}/extract_organ_source_read.sh \$(pwd)\n\tmv *png $CONFIG{'SOT_PROJECT_DIR'}/cfDNA_length\nfi\n";
		#print SHELL "if [ ! -f $sample.sort.bam.IS ]\nthen\n";
		#print SHELL "\tsh /data/pipeline/NGS_Target_Sequencing/Organ_Transplant/IS.sh $sample.sort.bam\nfi\n";
		#print SHELL "if [ ! -f $sample.UT.sscs.sort.bam.IS ]\nthen\n";
		#print SHELL "\tsh /data/pipeline/NGS_Target_Sequencing/Organ_Transplant/IS.sh $sample.UT.sscs.sort.bam\nfi\n";
		print SHELL "if [[ -f $sample.mp.tab.xd-cfDNA.insertsize ]]\nthen\n";
		print SHELL "\t/data/soft/bin/Rscript /data/pipeline/NGS_Target_Sequencing/Organ_Transplant/IS.sh.R $sample\nfi\n";
		print SHELL "cp *pdf $CONFIG{'SOT_PROJECT_DIR'}/report/$output_prefix/\n\n";
	}
	print SHELL "#4#6 500M QC\n";
	print SHELL "export LD_LIBRARY_PATH=/usr/local/lib64:/usr/local/lib:/data/soft/lib:/usr/local/lib64:/usr/local/lib:/data/soft/lib:/opt/gridengine/lib/lx-amd64:/opt/openmpi/lib\n";
	print SHELL "if [ ! -f $sample.QC ]\nthen\n";
	print SHELL "\tsh /data/pipeline/OT.sh $sample.sort.bam $sample.UT.sort.bam $sample.UT.sscs.sort.bam $CONFIG{'BED_SOT'} > $sample.QC\nfi\n";
	print SHELL "if [ ! -f QC/$sample.QC.xls ]\nthen\n";
	print SHELL "\tsh $CONFIG{'SOT_SCRIPT_DIR'}/QC.sh $sample.withtag.R1.fq.gz $sample.withtag.R2.fq.gz $sample\nfi\n\n";
	close SHELL;

	print STDERR "$sample.script.sh\tcreated\n";
	system("cp $sample.script.sh $CONFIG{SOT_PROJECT_DIR}/SHELL/$output_prefix");
	# print list 
	unless($skipIFD)
	{
		open REPORTLIST,">>$CONFIG{SOT_PROJECT_DIR}/report/$output_prefix/$output_prefix.Pathogen.list" or die $!;
		print REPORTLIST "$CONFIG{SOT_PROJECT_DIR}/$TASKS{$sample}/$sample/$sample.kraken.cdc.2.report.txt\n";
		close REPORTLIST;
	}
	$STATUS{$sample}{'INITIAL'} = "DONE";
}

sub shell_for_TUMOR
{
	my $sample=shift;
	unless($sample)
	{
		print STDERR "program <sample> <experiment type>[TUMOR:NORMAL]\n";
		print STDERR "Initial failed. Exiting...\n";
		exit;
	}
	my @R1R2 = @{$OSS_FILEs{$sample}};
	my @a = split("/",$R1R2[0]);
	my $dir = join("/",@a[3..$#a-1]);
	my $R1 = $a[-1];
	my $R2 = $R1; $R2 =~ s/R1/R2/;

	chdir($CONFIG{'TUMOR_PROJECT_DIR'});
	unless( -s $TASKS{$sample})
	{
		mkdir($TASKS{$sample},0755) or die "mkdir $TASKS{$sample} Failed in $CONFIG{'TUMOR_PROJECT_DIR'}";
	}
	unless( -s "$TASKS{$sample}/$sample")
	{
		mkdir("$TASKS{$sample}/$sample",0755)  or die "mkdir $TASKS{$sample}/$sample Failed in $CONFIG{'TUMOR_PROJECT_DIR'}";
	}
	unless( -s "$TASKS{$sample}/Fusion")
	{
		mkdir("$TASKS{$sample}/Fusion",0755)  or die "mkdir $TASKS{$sample}/Fusion Failed in $CONFIG{'TUMOR_PROJECT_DIR'}";
	}

	#unless( -s "$TASKS{$sample}/Mutect")
	#{
	#	mkdir("$TASKS{$sample}/Mutect",0755)  or die "mkdir $TASKS{$sample}/Mutect Failed in $CONFIG{'TUMOR_PROJECT_DIR'}";
	#}
	chdir("$TASKS{$sample}/$sample");

	unless( -s $R1 )
	{
		system("ln -s $CONFIG{'TUMOR_RAWDATA_DIR'}/$dir/$R1 .");
	}
	unless( -s $R2 )
	{
		system("ln -s $CONFIG{'TUMOR_RAWDATA_DIR'}/$dir/$R2 .");
	}

	my $type="";
	if($sample=~/ZLX/)
	{
		$type = "NORMAL";
	}elsif($sample=~/ZL/)
	{
		$type = "TUMOR";
	}

	my $r1 = $R1;
	my $r2 = $R2;

	#open SHELL,">run.$sample.sh" or die $!;
	if($type eq "TUMOR")
	{
		#
		system("echo sh /data/pipeline/NGS_Target_Sequencing/ctDNA/tumor.sh $r1 $r2 $sample|qsub -cwd -l vf=5G -P SOT -V -N $sample.ctdna.sh");
		#chdir("$CONFIG{'TUMOR_PROJECT_DIR'}/$TASKS{$sample}/$sample");
		system("echo sh /data/pipeline/NGS_Target_Sequencing/ctDNA/mutect_queue.sh $sample.sscs.sort.bam ../Normal/*.bam $sample > run.mutect.sh");
		chdir("$CONFIG{'TUMOR_PROJECT_DIR'}/$TASKS{$sample}/Fusion");
		system("echo sh /data/pipeline/NGS_Target_Sequencing/ctDNA/fusion.sh $sample.sort.bam $sample > run.fusion.sh");
	}elsif($type eq "NORMAL")
	{
		#
		system("echo sh /data/pipeline/NGS_Target_Sequencing/ctDNA/normal.sh $r1 $r2 $sample|qsub -cwd -l vf=5G -P SOT -V -N $sample.normal.sh");
	}

	#call variants

	#filter variants

	#annotate variants

	#report
	#close SHELL;
	$STATUS{$sample}{'INITIAL'} = "DONE";
}
sub load_config
{
	my $conf=shift;
	#my $hash=shift;
	open CFG, $conf or die $!;
	while(<CFG>)
	{
		chomp;
		next if(/#/);
		next if(/^$/);
		my @a = split /=/,$_;
		#$$hash{$a[0]} = $a[1];
		$CONFIG{$a[0]} = $a[1];
	}
	close CFG;
}
sub md5check
{
	my $sample=shift;
	my $basedir=shift;
	my $dir = shift;## unfinished 
	my $downloadpath=shift;

	my @d=split("/",$downloadpath);
	my $ok = '';
	my ($R1,$R2,$R1md5,$R2md5,$md5file) = ('','','','','');
	chdir($basedir);

	if($downloadpath=~/RRG/ || $downloadpath=~/s836/ || $downloadpath=~/s1293/)
	{
		$manufacturer = "wixinextcode";
		print STDERR "Manufacturer:\tNextcode\n";
		unless( -s "$basedir/$dir/$sample\_combined_R1.fastq.gz.md5" )
		{
			if(system("ossutil64 cp -j 1 $ossaccount -rf $downloadpath/$sample\_combined_R1.fastq.gz.md5 ."))
			{
				#die "Failed downloading $sample\_combined_R1.fastq.gz.md5 into $basedir/$dir\n";
				warn "Failed downloading $sample\_combined_R1.fastq.gz.md5, creating md5 file...";
				chdir("$basedir/$dir");
				system("md5sum $sample\_combined_R1.fastq.gz >$sample\_combined_R1.fastq.gz.md5");
				warn "done\n";
				chdir($basedir);
			}
		}
		unless( -s "$basedir/$dir/$sample\_combined_R2.fastq.gz.md5" )
		{
			if(system("ossutil64 cp -j 1 $ossaccount -rf $downloadpath/$sample\_combined_R2.fastq.gz.md5 ."))
			{
				#die "Failed downloading $sample\_combined_R2.fastq.gz.md5 into $basedir/$dir\n";
				warn "Failed downloading $sample\_combined_R2.fastq.gz.md5, creating md5 file...";
				chdir("$basedir/$dir");
				system("md5sum $sample\_combined_R2.fastq.gz >$sample\_combined_R2.fastq.gz.md5");
				warn "done\n";
				chdir($basedir);
			}
		}
		$R1md5 = "$basedir/$dir/$sample\_combined_R1.fastq.gz.md5";
		$R2md5 = "$basedir/$dir/$sample\_combined_R2.fastq.gz.md5";
		#print "R1md5==$R1md5\nR2md5==$R2md5\n";
	}elsif($downloadpath=~/SR\d+/)
	{
		$manufacturer = "shihe";
		print STDERR "Manufacturer:\tShihe\n";
		$md5file="md5sum.txt";
	}else{
		$manufacturer = "genwiz";
		print STDERR "Manufacturer:\tGenwiz\n";
		$md5file="md5.md5";
	}

	if( $md5file ne '' )
	{
		$R1md5 = "$basedir/$dir/$sample.R1.md5";
		$R2md5 = "$basedir/$dir/$sample.R2.md5";

		unless( -s "$basedir/$dir/$md5file" )
		{
			print STDERR "Downloading md5:\t$basedir/$dir/$md5file\n";
			if(system("ossutil64 cp -j 1 $ossaccount -rf $downloadpath/$md5file ."))
			{
				warn "Failed downloading $downloadpath/$md5file into $basedir/$dir\n";
			}
		}else{	
			system("grep $sample $basedir/$dir/$md5file |grep R1 > $R1md5");
			system("grep $sample $basedir/$dir/$md5file |grep R2 > $R2md5");
		}
	}
	#print "R1md5=$R1md5\n";
	$R1 = `cat $R1md5 |awk '{print \$2;}'`;chomp($R1);
	$R2 = `cat $R2md5 |awk '{print \$2;}'`;chomp($R2);

	chdir($dir);
	if ( -s $R1)
	{
		if(system("md5sum -c $R1md5"))
		{
			$ok = "0";
		}else{
			$ok = "1";
		}
	}else{
		$ok = "0";
	}
	if ( -s $R2)
	{
		if(system("md5sum -c $R2md5"))
		{
			$ok .= "0";
		}else{
			$ok .= "2";
		}
	}else{
		$ok .= "0";
	}
	return ($ok eq "") ? "00" : $ok;
}

sub load_sample_from_table
{
	my $table = shift;
	my $output_prefix=shift;
	#my $savedir = shift;
	print STDERR "load sample names from $table\n";
	my $workbook = $table;
	if($table=~/.xls$/)
	{
		use Spreadsheet::ParseExcel;
		#use Spreadsheet::ParseExcel::FmtUnicode;
		my $parser = Spreadsheet::ParseExcel->new();
		#my $formatter = Spreadsheet::ParseExcel::FmtUnicode->new(Unicode_Map=>"CP936");
		$workbook = $parser->parse($table);#,$formatter);
		if( ! defined $workbook )
		{
			print STDERR $parser->error(),"\n";
		}

		for my $worksheet ( $workbook->worksheets() ) 
		{

			my ( $row_min, $row_max ) = $worksheet->row_range();
			my ( $col_min, $col_max ) = $worksheet->col_range();

			for my $row ( $row_min .. $row_max ) 
			{
				my $cell = $worksheet->get_cell( $row, 2 );
				next unless $cell;
				#my $sp = $cell->unformatted();
				my $sp = $cell->value();
				next unless($sp =~ /S/);
				if($sp=~/$regulex_QY/)
				{
					next if($sp=~/DC/);
					$TASKS{$sp} = "$1$2";
					open my $SFH,">>$CONFIG{'STATUDIR'}/$sp.statu" or die $!; 
					$STATUS_FH{$sp} = $SFH;
				}elsif($sp=~/$regulex_ZL/)
				{
					$TASKS{$sp} = "$1$2";
					open my $SFH,">>$CONFIG{'STATUDIR'}/$sp.statu" or die $!; 
					$STATUS_FH{$sp} = $SFH;
				}elsif($sp=~/$regulex_NIPT/)
				{
					$TASKS{$sp} = "$1$2";
					open my $SFH,">>$CONFIG{'STATUDIR'}/$sp.statu" or die $!; 
					$STATUS_FH{$sp} = $SFH;
				}elsif($sp=~/$regulex_QZ/)
				{
					$TASKS{$sp} = "$1$2";
					open my $SFH,">>$CONFIG{'STATUDIR'}/$sp.statu" or die $!; 
					$STATUS_FH{$sp} = $SFH;
				}
			}
		}
	}else{
		print $table,"\n";
		open TB,$table or die $!;
		while(<TB>)
		{
			chomp;
			#print STDERR "TB:$_";
			my @a = split /\t/,$_;
			next if($a[2]=~/DC/);
			if($a[2]=~/$regulex_QY/)
			{
				next if($a[2]=~/DC/);
				$TASKS{$a[2]} = "$1$2";
				open my $SFH,">>$CONFIG{'STATUDIR'}/$a[2].statu" or die $!; 
				$STATUS_FH{$a[2]} = $SFH;				
			}elsif($a[2] =~ /$regulex_ZL/)
			{
				$TASKS{$a[2]} = "$1$2";
				open my $SFH,">>$CONFIG{'STATUDIR'}/$a[2].statu" or die $!; 
				$STATUS_FH{$a[2]} = $SFH;
			}elsif($a[2] =~ /$regulex_NIPT/)
			{
				$TASKS{$a[2]} = "$1$2";
				open my $SFH,">>$CONFIG{'STATUDIR'}/$a[2].statu" or die $!; 
				$STATUS_FH{$a[2]} = $SFH;
			}elsif($a[2] =~ /$regulex_QZ/)
			{
				$TASKS{$a[2]} = "$1$2";
				open my $SFH,">>$CONFIG{'STATUDIR'}/$a[2].statu" or die $!; 
				$STATUS_FH{$a[2]} = $SFH;
			}
		}
		close TB;
	}
	my @lines=`ossutil64 ls $ossBucketP $ossaccount`;
	chomp(@lines);
	for(my $i=0;$i<@lines;$i++)
	{
		my @a = split /\s+/,$lines[$i];
		my @b = split /\//,$a[-1];
		next if($a[-1]=~/\/$/ || $a[-1]!~/oss/);
		#next if($savedir ne join("/",@b[3..$#b-1]));
		my $sp = (split("_",$b[-1]))[0];
		if(exists $TASKS{$sp})
		{
			push(@{$OSS_FILEs{$sp}},$a[-1]);
			#print STDERR "$TASKS{$sp}\t$sp\t$b[-1]\n";
		}
	}
}
sub load_sample_from_oss
{
	my $savedir = shift; 
	my $output_prefix=shift;
	my @lines=();
	for my $path(@$savedir)
	{
		#print STDERR "$ossBucketP\/$path $ossaccount\n";
		my $osscmd = "ossutil64 ls $ossBucketP/$path $ossaccount";
		$osscmd =~ s/\/\//\//g;
		$osscmd =~ s/oss:\//oss:\/\//;
		my @tmp=`$osscmd`;
		push(@lines,@tmp);
	}
	chomp(@lines);
	for(my $i=0;$i<@lines;$i++)
	{
		#print STDERR "$lines[$i]\n";
		my @a = split /\s+/,$lines[$i];
		next if($a[-1]=~/\/$/ || $a[-1]!~/oss/);
		# $a[-1] oss file path
		my @b = split /\//,$a[-1];
		my $sp = "";
		my $targetsp = "";
		#  $b[-1] file name
		if($b[-1] =~ /([^_]+)_/)
		{
			$sp = $1;
			print STDERR "$sp\t$b[-1]\t$a[-1]\n";
		}
		if($b[-1] =~ /-(\w+-?\d+)_/)
		{
			$targetsp=$1;

		}else{
			$targetsp = $sp;
		}
		$TASKS2TARGET{$sp} = $targetsp;
		#if($manufacturer eq "wuxinextcode")
		#{
		#		my @c = split /-/,$sp;
		#		#$sp = "$c[2]-$c[3]";
		#		$sp = $c[-1];
		#}
		if($sp=~/$regulex_QY/)
		{
			next if($sp=~/DC/);
			print STDERR "$1$2\t$sp\t$b[-1]\n";
			$TASKS{$sp} = "$1$2";
			push (@{$OSS_FILEs{$sp}},$a[-1]) if ($a[-1] =~ /gz$/);
			open my $SFH,">>$CONFIG{'STATUDIR'}/$sp.statu" or die $!; 
			$STATUS_FH{$sp} = $SFH;
		}elsif($sp=~/$regulex_ZL/)
		{
			$TASKS{$sp} = "$1$2";
			push (@{$OSS_FILEs{$sp}},$a[-1]) if ($a[-1] =~ /gz$/);
			open my $SFH,">>$CONFIG{'STATUDIR'}/$sp.statu" or die $!; 
			$STATUS_FH{$sp} = $SFH;
		}elsif($sp=~/$regulex_NIPT/)
		{
			$TASKS{$sp} = "$1$2";
			push (@{$OSS_FILEs{$sp}},$a[-1]) if ($a[-1] =~ /gz$/);
			open my $SFH,">>$CONFIG{'STATUDIR'}/$sp.statu" or die $!;
			$STATUS_FH{$sp} = $SFH;
		}elsif($sp=~/$regulex_QZ/)
		{
			$TASKS{$sp} = "$1$2";
			push (@{$OSS_FILEs{$sp}},$a[-1]) if ($a[-1] =~ /gz$/);
			open my $SFH,">>$CONFIG{'STATUDIR'}/$sp.statu" or die $!; 
			$STATUS_FH{$sp} = $SFH;
		}elsif($sp=~/$regulex_IF/)
		{
			$TASKS{$sp} = "$1$2$3";
			#print "$sp\t$1$2$3\t$a[-1]\n";
			push(@{$OSS_FILEs{$sp}},$a[-1]) if ($a[-1] =~ /gz$/);
			open my $SFH,">>$CONFIG{'STATUDIR'}/$sp.statu" or die $!;
			$STATUS_FH{$sp} = $SFH;
		}elsif($sp=~/$regulex_YS/){
			$TASKS{$sp} = "$1$2";
			push (@{$OSS_FILEs{$sp}},$a[-1]) if ($a[-1] =~ /gz$/);
			open my $SFH,">>$CONFIG{'STATUDIR'}/$sp.statu" or die $!;
			$STATUS_FH{$sp} = $SFH;
		}elsif($sp=~/$regulex_QH/)
		{
			$TASKS{$sp} = "$1$2";
			push (@{$OSS_FILEs{$sp}},$a[-1]) if ($a[-1] =~ /gz$/);
			open my $SFH,">>$CONFIG{'STATUDIR'}/$sp.statu" or die $!;
			$STATUS_FH{$sp} = $SFH;
		}elsif($a[-1]=~/Total_QC|xlsx/i)
		{
			chdir($CONFIG{'SOT_RAWDATA_DIR'});
			print STDERR "Downloading QC:$a[-1]\n";
			system("ossutil64 cp -j 1 $ossaccount -rf $a[-1] .");
			chdir("$CONFIG{SOT_PROJECT_DIR}/report/$output_prefix");
			my $qc2=$a[-1];
			#$qc2=~s/oss:\/\/rrgene\///;
			$qc2=~s/$ossBucketP//;
			print STDERR "SeqQC:\t$a[-1]\nlnk $CONFIG{SOT_RAWDATA_DIR}/$qc2 to $CONFIG{SOT_PROJECT_DIR}/report/$output_prefix\n";
			try{
				system("ln -s $CONFIG{SOT_RAWDATA_DIR}/$qc2 .") or warn "$qc2 exists\n";;
			}catch{
				warn "$qc2 exists\n";
			};
			chdir($CONFIG{'SOT_RAWDATA_DIR'});
		}
	}
}

sub get_savedir
{
	my $osspath = shift;
	my $output_prefix = shift;
	my @savedir = ();
	my %date_data=();
	my $latest="";
	if($osspath ne "")
	{
		print STDERR "Getting savedir through user input path\n";
		my @b = split("/",$osspath);
		push(@savedir,join("/",@b[3..$#b-1]));
	}else{
		print STDERR "Getting savedir through latest date mark when no path input by user\n";
		#use Date::Manip;
		#&Date_Init("DateFormat=YYYY-MM-DD");


		my $predate="";
		my $postdate="";

		my @lines=`ossutil64 ls $ossBucketP $ossaccount`;
		chomp(@lines);

		for(my $i=1;$i<@lines;$i++)
		{	
			my @a = split /\s+/,$lines[$i];
			#print STDERR "$lines[$i]\n";
			next if($a[-1]=~/\/$/ || $a[-1]!~/oss/);
			my @b = split /\//,$a[-1];
			my $tmppath = "";
			if($a[-1]=~/s836/ || $a[-1]=~/s1293/)
			{
				$tmppath = join("/",@b[4..$#b-1]);
			}else{
				$tmppath=join("/",@b[3..$#b-1]);
			}
			$date_data{$a[0]}{$tmppath} = '';
			#print "$predate\t";
			if($predate eq "")
			{
				$predate = ParseDate($a[0]);
				$latest = $a[0];
				#print STDERR "test: $predate\t$latest\n";
			}else{
				$postdate = ParseDate($a[0]);
				my $flag = Date_Cmp($predate,$postdate);
				if($flag < 0)#predate is earlier date
				{
					$latest = $a[0];
					$predate = $postdate;
				}

			}
		}
		@savedir=keys %{$date_data{$latest}};
		print STDERR "Setting date : $latest\n";
	}

	for my $path(@savedir)
	{
		die "No Newest Data Detected\n" if($path eq "");
		#print STDERR "Setting savedir: $path\n";
	}
	return @savedir;
}

sub load_status
{
	my $sample=shift;
	print STDERR "load status of $sample\n";
	if ( -s "$CONFIG{'STATUDIR'}/$sample.statu" )
	{
		open IN,"$CONFIG{'STATUDIR'}/$sample.statu" or die $!;
		while(<IN>)
		{
			chomp;
			my @a = split /\s+/,$_;
			$STATUS{$sample}{$a[0]} = $a[1];
		}
		close IN;
	}
}
1;
