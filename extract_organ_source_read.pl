use strict;
use warnings;
use Getopt::Long;
use RRG::SOT_app;

my ($bam,$snptab,$vcf,$bed,$donor,$recipient)=("","","","","","");
my $help;
my $features;
my $refresh;

GetOptions(
	"--bam=s" => \$bam,
	"--snptab=s" => \$snptab,
	"--vcf=s"=>\$vcf,
	"--bed=s"=>\$bed,
	"--donor=s" => \$donor,
	"--recipient=s" => \$recipient,
	"--features" => \$features,
	"--refresh" => \$refresh,
	"--help" => \$help
);
print STDERR $snptab,"\n";
if(($snptab eq "" && $vcf eq "") || $help)
{
	&usage();
}
if($features && $bam eq "")
{
	&usage();
}

$bed ||= "/data/pipeline/NGS_Target_Sequencing/Organ_Transplant/config/pt3.SNP.GRCh38.full.bed";

my %Organ_Snp_incfDNA = (); # store donor's SNP in cfDNA 
my %Donor_Snp = (); # store donor's SNP of white blood cell
my %Recipient_Snp = (); #store recipient's SNP of white blood cell
my %rs = ();


if($features)
{
	my $effectiveSNP = "$snptab.ex.informativeSNP.effective";
	if( -f $effectiveSNP)
	{
		MXddcfDNA_info_from_bam($bam,$snptab,$effectiveSNP);
	}else{
		die "no $effectiveSNP found\n";
	}
}else{
	
	#load rs of SNPs of population
	load_rs_from_bed($bed,\%rs);
	#load GT of donor's or recipient's from tab
	GT_from_TAB($donor,\%Donor_Snp,0) if($donor ne "");
	GT_from_TAB($recipient,\%Recipient_Snp,1) if($recipient ne "");

	#Donor SNP Genotyping
	if($refresh || !-f "$snptab.ex")
	{
		GT_from_TAB_MXcfDNA($snptab,\%rs,\%Recipient_Snp,\%Donor_Snp,\%Organ_Snp_incfDNA);
	
		#GT_from_VCF_MXcfDNA($vcf);
		
		open OUT,">$snptab.ex.informativeSNP" or die $!;
		foreach my $k(keys %Organ_Snp_incfDNA)
		{
			my $kk = $k;
			$kk=~s/:/\t/;
			print OUT "$kk\t$Organ_Snp_incfDNA{$k}\n";
		}
		close OUT;
		
		open OUT2,">$snptab.ex.informativeSNP.effective" or die $!;
		open EX,"$snptab.ex" or die $!;
		while(<EX>)
		{
			my @a = split;
			if(exists $Organ_Snp_incfDNA{"$a[1]:$a[2]"} || $a[0] eq "Sample")
			{
				print OUT2 $_;
			}
		}
		close EX;
		close OUT2;
	}
}



sub usage
{
	print "This software is apt to extract read from a alignment by input genotype as baits  .\n";
	print "Usage:\n";
	print "\t--bam\t<file>\talignment with bam format[required]\n";
	print "\t--snptab\t<file>\tsnp genotype of mixed sample[required]\n";
	print "\t--donor\t<file>\tsnp genotype of donor sample, for including\n";
	print "\t--recipient\t<file>\tsnp genotype of recipient sample,for excluding\n";
	print "\t--refresh\t\tupdate existed results of genotyping\n";
	print "\t--features\t\toutput features of dd-cfDNA and rd-cfDNA\n";
	print "\t--help\t\tprint help information\n";
	exit;
}
