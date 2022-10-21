package SOT;

use warnings;
use Try::Tiny;
use Cwd 'abs_path';
use File::Basename;

require Exporter;

@ISA = qw(Exporter);

@EXPORT = qw(run start_SE start_PE);

sub start_SE
{
	my ($gdnaID,$gdnaR1,$cfdnaID,$cfdnaR1,$config)=@_;

	$gdnaR1 = abs_path($gdnaR1);
	$cfdnaR1 = abs_path($cfdnaR1);

	#loading dependencies
	my %CONFIG = &load_config($config);

	#create outdir
	my $outdir = $cfdnaID;
	unless (-e $outdir)
	{
		mkdir($outdir) or die "Failed: create $outdir\n";
	}

	unless (-e "$outdir/data")
	{
		mkdir("$outdir/data") or die "Failed: create $outdir/data\n";
	}

	#chdir($outdir) or die "Failed: enter $outdir\n";

	unless(-e "$outdir/data/$gdnaID.raw.R1.fastq.gz")
	{
		print "$outdir/data/$gdnaID.raw.R1.fastq.gz";
		symlink ($gdnaR1, "$outdir/data/$gdnaID.raw.R1.fastq.gz") or die "Failed link $gdnaR1 to $outdir/data/$gdnaID.raw.R1.fastq.gz\n";
	}
	
	unless(-e "$outdir/data/$cfdnaID.raw.R1.fastq.gz")
	{
		symlink ($cfdnaR1, "$outdir/data/$cfdnaID.raw.R1.fastq.gz") or die "Failed link $cfdnaR1 to $outdir/data/$cfdnaID.raw.R1.fastq.gz\n";
	}
	

	#generate shell script
	open SHELL, ">$outdir/Snakefile" or die $!;
	print SHELL "samples = [\"$gdnaID\", \"$cfdnaID\"]\n\n";

	my $is_PE = 0;

	&all();

	&filter(\%CONFIG,$is_PE);

	&aln2genome(\%CONFIG,$is_PE);

	&dedup(\%CONFIG);

	&genotyping(\%CONFIG);

	&ddcfDNA_frac(\%CONFIG);

	#&QC($cfdnaID);

	close SHELL;

}

sub start_PE
{
	my ($gdnaID,$gdnaR1,$gdnaR2,$cfdnaID,$cfdnaR1,$cfdnaR2,$config)=@_;

	$gdnaR1 = abs_path($gdnaR1);
	$gdnaR2 = abs_path($gdnaR2);
	$cfdnaR1 = abs_path($cfdnaR1);
	$cfdnaR2 = abs_path($cfdnaR2);

	#loading dependencies
	my %CONFIG = &load_config($config);

	#create outdir
	my $outdir = $cfdnaID;
	unless (-e $outdir)
	{
		mkdir($outdir) or die "Failed: create $outdir\n";
	}

	unless (-e "$outdir/data")
	{
		mkdir("$outdir/data") or die "Failed: create $outdir/data\n";
	}

	#chdir($outdir) or die "Failed: enter $outdir\n";

	unless(-e "$outdir/data/$gdnaID.raw.R1.fastq.gz")
	{
		print "$outdir/data/$gdnaID.raw.R1.fastq.gz";
		symlink ($gdnaR1, "$outdir/data/$gdnaID.raw.R1.fastq.gz") or die "Failed link $gdnaR1 to $outdir/data/$gdnaID.raw.R1.fastq.gz\n";
	}
	unless(-e "$outdir/data/$gdnaID.raw.R2.fastq.gz")
	{
		symlink ($gdnaR2, "$outdir/data/$gdnaID.raw.R2.fastq.gz") or die "Failed link $gdnaR2 to $outdir/data/$gdnaID.raw.R2.fastq.gz\n";
	}
	
	unless(-e "$outdir/data/$cfdnaID.raw.R1.fastq.gz")
	{
		symlink ($cfdnaR1, "$outdir/data/$cfdnaID.raw.R1.fastq.gz") or die "Failed link $cfdnaR1 to $outdir/data/$cfdnaID.raw.R1.fastq.gz\n";
	}
	unless(-e "$outdir/data/$cfdnaID.raw.R2.fastq.gz")
	{
		symlink ($cfdnaR2, "$outdir/data/$cfdnaID.raw.R2.fastq.gz") or die "Failed link $cfdnaR2 to$outdir/data/$cfdnaID.raw.R2.fastq.gz\n";
	}

	#generate shell script
	open SHELL, ">$outdir/Snakefile" or die $!;
	print SHELL "samples = [\"$gdnaID\", \"$cfdnaID\"]\n\n";

	my $is_PE = 1;

	&all();

	&filter(\%CONFIG,$is_PE);

	&aln2genome(\%CONFIG,$is_PE);

	&targetSelect(\%CONFIG);

	&dedup(\%CONFIG);

	&genotyping(\%CONFIG);

	&ddcfDNA_frac(\%CONFIG);

	#&QC($cfdnaID);

	close SHELL;

}

sub filter
{
	my ($CONFIG, $is_PE) = @_;
	if ($is_PE == 1)
	{
		my $rule_PE = <<"EOF";
rule filter:
	input:
		r1 = "data/{sample}.raw.R1.fastq.gz",
		r2 = "data/{sample}.raw.R2.fastq.gz"
	output:
		o1 = "data/{sample}.clean.R1.fastq.gz",
		o2 = "data/{sample}.clean.R2.fastq.gz",
		json = "data/{sample}.json",
		html = "data/{sample}.html"
	shell:
		"fastp -i {input.r1} -I {input.r2} -o {output.o1} -O {output.o2} -j {output.json} -h {output.html}  -l 31 -w 16 -Q"
EOF
		print SHELL $rule_PE,"\n";
	}else{
		my $rule_SE = << "EOF";
rule filter:
	input:
		"data/{sample}.raw.R1.fastq.gz"
	output:
		out = "data/{sample}.clean.R1.fastq.gz",
		json = "data/{sample}.json",
		html = "data/{sample}.html"
	shell:
		"fastp -i {input} -o {output.out} -j {output.json} -h {output.html}  -l 31 -w 16 -Q"
EOF
		print SHELL $rule_SE,"\n";
	}

}

sub aln2genome
{
	#genotyping
	my ($CONFIG,$is_PE) = @_;
	if ($is_PE == 1)
	{
		my $rule_PE = <<"EOF";
rule aln2genome:
	input:
		r1 = "data/{sample}.clean.R1.fastq.gz",
		r2 = "data/{sample}.clean.R2.fastq.gz"
	output:
		bam = "mapped/{sample}.sort.bam",
		bai = "mapped/{sample}.sort.bam.bai",
		targetbam = "mapped/{sample}.target.sort.bam"
	params:
		rg = r"\@RG\\tID:{sample}\\tSM:{sample}",
		fc = r"mapping_quality >= 30 and not (unmapped or secondary_alignment) and not ([XA] != null or [SA] != null) and ((not [MD]=~/\\^/ and [NM]<=2) or (cigar=~/[ID]/ and [NM]<8))"
	threads: 16
	shell:"""
		bwa mem $$CONFIG{'hg38'}  {input.r1} {input.r2} -t {threads} -R '{params.rg}' -v 1 | samtools sort -T mapped/{wildcards.sample}.tmp -o {output.bam} -;
		sambamba index {output.bam}
		sambamba view -t 10 -L $$CONFIG{'PANEL_BED'} -h -f bam -F '{params.fc}' -o {output.targetbam} {output.bam};
		"""
EOF
		print SHELL "$rule_PE\n\n";
	}else{
		my $rule_SE = <<"EOF";
rule aln2genome:
	input:
		"data/{sample}.clean.R1.fastq.gz"
	output:
		bam = "mapped/{sample}.sort.bam",
		bai = "mapped/{sample}.sort.bam.bai",
	params:
		rg = r"\@RG\\tID:{sample}\\tSM:{sample}",
	threads: 16
	shell:"""
		bwa mem $$CONFIG{'hg38'}  {input} -t {threads} -R '{params.rg}' -v 1 | samtools sort -T mapped/{wildcards.sample}.tmp -o {output.bam} -;
		sambamba index {output.bam}
		"""
EOF
		print SHELL "$rule_SE\n\n";
	}
	
}

sub targetSelect
{
	my $CONFIG = shift;
	my $rule = <<"EOF";
rule targetSelect:
	input:
		bam = "mapped/{sample}.sort.bam",
		bai = "mapped/{sample}.sort.bam.bai",
	output:
		targetbam = "mapped/{sample}.target.sort.bam"
	params:
		fc = r"mapping_quality >= 30 and not (unmapped or secondary_alignment) and not ([XA] != null or [SA] != null) and ((not [MD]=~/\\^/ and [NM]<=2) or (cigar=~/[ID]/ and [NM]<8))"
	threads: 10
	shell:"""
		sambamba view -t {threads} -L $$CONFIG{'PANEL_BED'} -h -f bam -F '{params.fc}' -o {output.targetbam} {output.bam};
		"""
EOF
		print SHELL "$rule\n\n";
}

sub dedup
{
	#remove duplicates
	my $CONFIG = shift;
	my $barcodeLength = 1;
	if($barcodeLength >= 2)
	{
		my $rule =<<"EOF";
rule dedup:
	input:
		"mapped/{sample}.target.sort.bam"
	output:
		dedup = "dedup/{sample}.target.dedup.bam",
		dedupsort = "dedup/{sample}.target.dedup.sort.bam"
	shell:"""
		$$CONFIG{'DSpath'}/ConsensusMaker.v2.py --infile {input} --outfile {output.dedup} --minmem $$CONFIG{'minMem'} --maxmem $$CONFIG{'maxMem'} --cut_off $$CONFIG{'cutOff'} --Ncut_off $$CONFIG{'nCutOff'} --read_type $$CONFIG{'readTypes'} --filt $$CONFIG{'filtersSet'} --isize $$CONFIG{'iSize'}  --read_length $readLength --rep_filt $repFilt
		samtools view -h {output.dedup}|awk '{if (\$6!~/^[[:digit:]]+N[[:digit:]]+D|^[[:digit:]]+P|[[:digit:]]+P\$|^[[:digit:]]+I/) print \$0}' | samtools sort -T dedup/{wildcards.sample}.tmp -o {output.dedupsort} -\@ 8 -
		"""
EOF
	print SHELL "$rule\n\n";
	}else{
		my $rule = <<"EOF";
rule dedup:
	input:
		"mapped/{sample}.target.sort.bam"
	output:
		Mtxt = "dedup/{sample}.picard.txt",
		dedupsort = "dedup/{sample}.target.dedup.sort.bam"
	shell:
		"picard MarkDuplicates I={input} O={output.dedupsort} M={output.Mtxt} R=$$CONFIG{'hg38'} VALIDATION_STRINGENCY=SILENT QUIET=true REMOVE_DUPLICATES=true"
EOF
	print SHELL "$rule\n\n";
	}
}

sub genotyping
{
	#call variant
	my $CONFIG = shift;
	my $rule =<<"EOF";
rule genotyping:
	input:
		"dedup/{sample}.target.dedup.sort.bam"
	output:
		vcf = "genotype/{sample}.mp.vcf",
		tab = "genotype/{sample}.mp.tab"
	shell:"""
		samtools mpileup -A -uvf $$CONFIG{'hg38'} -l $$CONFIG{'PANEL_BED'} -t DP,AD {input} | bcftools call --multiallelic-caller --keep-alts --targets-file $$CONFIG{'PANEL_BED'} > {output.vcf}
		perl $$CONFIG{'VCF2TAB'} {output.vcf} $$CONFIG{'PANEL_BED'} {output.tab}
		"""
EOF
	print SHELL "$rule\n\n";
}

sub ddcfDNA_frac
{
	# ddcfDNA quantification
	my $CONFIG = shift;
	my $rule = <<"EOF";
rule ddcfDNAFrac:
	input:
		gdna = "genotype/{sample[0]}.mp.tab",
		cfdna = "genotype/{sample[1]}.mp.tab"
	output:
		"ddcfDNA_percent.xls"
	shell:
		"sh $$CONFIG{'SOT_SCRIPT_DIR'}/tabex.sh {sample[0]} {sample[1]}"
EOF
	print SHELL "$rule\n\n";
}

sub all
{
	my $rule = <<"EOF";
rule all:
	input:
		"ddcfDNA_percent.xls"
EOF
	print SHELL "$rule\n";
}
sub QC
{
	print SHELL "#4#6 500M QC\n";
	print SHELL "export LD_LIBRARY_PATH=/usr/local/lib64:/usr/local/lib:/data/soft/lib:/usr/local/lib64:/usr/local/lib:/data/soft/lib:/opt/gridengine/lib/lx-amd64:/opt/openmpi/lib\n";
	print SHELL "if [ ! -f $cfdnaID.QC ]\nthen\n";
	print SHELL "\tsh /data/pipeline/OT.sh $cfdnaID.sort.bam $cfdnaID.UT.sort.bam $cfdnaID.UT.sscs.sort.bam $$CONFIG{'BED_SOT'} > $cfdnaID.QC\nfi\n";
	print SHELL "if [ ! -f QC/$sample.QC.xls ]\nthen\n";
	print SHELL "\tsh $$CONFIG{'SOT_SCRIPT_DIR'}/QC.sh $cfdnaID.withtag.R1.fq.gz $cfdnaID.withtag.R2.fq.gz $cfdnaID\nfi\n\n";
}

sub load_config
{
	my $conf=shift;
	#my $hash=shift;
	open CFG, $conf or die $!;
	my %CONFIG = ();
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
	return %CONFIG;
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
		##&Date_Init("DateFormat=YYYY-MM-DD");


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

1;
