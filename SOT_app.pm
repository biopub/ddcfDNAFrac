package SOT_app;

use File::Basename;

require Exporter;

@ISA = qw(Exporter);

@EXPORT = qw(load_rs_from_bed VCF_to_TAB GT_from_TAB GT_from_TAB_MXcfDNA MXddcfDNA_info_from_bam extract_info_from_bam dump_cfDNA_length dump_cfDNA_AB);
#@EXPORT_OK = qw();

sub VCF_to_TAB
{
	my $vcf = shift;
	my $Hrs = shift;
	my $output = shift;
	my ($ref , $alt , $geno , $genobase , $DP , $refdp , $altdp , $rs , $rf) = ("" , "" , "" , "" , 0 , 0 , 0 , 0 , 0);
	open VCF , $vcf or die $!;
	open OUT,">$output" or die $!;
	print OUT "Sample\tChr\tPos\tRef\tAlt\tAll.dp\tRef.dp\tAlt.dp\trf\tGT\tGT.Base\tRsID\tRef.AF\n";
	my $sample = "";
	while(<VCF>)
	{
		chomp;
		next if(/^##/);
		my @a = split /\t/ , $_;
		if(/^#CHROM/)
		{
			$sample = $a[-1];
			next;
		}
		next if(length($a[3])>1);

		next unless(exists $$Hrs{"$a[0]:$a[1]"});
		my @format = split /:/ , $a[-1];
		my $DP = $format[-2];
		my @dp = split /,/ , $format[-1];
		my ($alt) = $a[4] =~ /(\S)/;
		if($alt eq ".")
		{
			$alt = $a[3];
			$altdp = 0;
		}
		($geno) = $format[0];
		#next if($geno eq "1/2"); #ignore tri-polymorphism snp
		if($geno eq "0/0")
		{
			$ref = $a[3];
			#$alt = $a[4];
			$refdp = $dp[0];
			$altdp = ($dp[1]) ? join(",",@dp[1..$#dp]) : 0;
			$genobase = "$ref|$ref"
		}elsif($geno eq "0/1")
		{
			$ref = $a[3];
			#$alt = $a[4];
			$refdp = $dp[0];
			$altdp = join(",",@dp[1..$#dp]);
			$genobase = "$ref|$alt"
		}elsif($geno eq "1/1")
		{
			$ref = $a[3];
			#$alt = $a[4];
			$refdp = ($dp[1]) ? $dp[0] : 0;
			$altdp = ($dp[1]) ? join(",",@dp[1..$#dp]) : $dp[0];
			$genobase = "$alt|$alt";
		}elsif($geno eq "./.")
		{
			$ref = $a[3];
			#$alt = $a[4];
			$refdp = 0;
			$altdp = 0;
			$genobase = ".|.";
		}elsif($geno eq "1/2") ##ignore tri-polymorphism snp
		{
			$ref = $a[3];
			$refdp = $dp[0];
			$altdp = $dp[1]+$dp[2];
			$genobase = $a[4];$genobase =~ s/,/\|/;
		}

		if($DP == 0)
		{
			$rf = 0;
		}else{
			$rf = $refdp/$DP;
		}
		$rs_and_maf = (exists $$Hrs{"$a[0]:$a[1]"})  ?  $$Hrs{"$a[0]:$a[1]"}{$ref}  :  "RS\tRef.AF";

		print OUT "$sample\t$a[0]\t$a[1]\t$ref\t$a[4]\t$DP\t$refdp\t$altdp\t$rf\t$geno\t$genobase\t$rs_and_maf\n";
	}
	close VCF;
	close OUT;
}

sub load_rs_from_bed
{
	my $bed = shift;
	my $hash = shift;
	open BED , $bed or die $!;
	while(<BED>)
	{
		chomp;
		my @a  = split /\t/ , $_;
		$$hash{"$a[0]:$a[2]"}{$a[-3]} = $a[3]."\t".(1-$a[-1]);
		$$hash{"$a[0]:$a[2]"}{$a[-2]} = $a[3]."\t".$a[-1];
	}
	close BED;
}

sub GT_from_TAB_MXcfDNA
{
	my $snp=shift; #tab format
	my $rs=shift;#hash store RS information
	my $recipient_Snp=shift; #reference  of hash
	my $Donor_Snp=shift; #reference  of hash
	my $Organ_Snp_incfDNA=shift; #reference  of hash
	#my $RS=shift;

	my $has_recipient_data = keys %$recipient_Snp;
	my $has_donor_data = keys %$Donor_Snp;
	#print STDERR "recipient:$has_recipient_data\t-\tDonor:$has_donor_data\n";
	my $output = $snp;
	$output =~ s/mp.tab/mp.tab.ex/;	
	open TT , $snp or die $!;
	open OO , ">$output" or die $!;
	print OO "Sample\tChr\tPos\tRef\tAlt\tAll.dp\tRef.dp\tAlt.dp\tRecipient\tRecipient.dp\tDonor\tDonor.dp\tD.rf\tRsID\tRR.AF\tD.hete.AF\tD.homo.AF\n";
	while(<TT>)
	{
		chomp;
		my @a = split /\t/ , $_;
		next if($a[5] < 20);
		next if($a[4]!~/,/&&length($a[4])>1);
		next if($a[9]=~/2/);
		my (@b) = split(/,/,$a[4]);
		my (@dp_of_b) = split(/,/,$a[7]);


		my $donor_sites = ""; #get a specific donor sites from GT of donor
		my $recipient_sites = "";#get a specific recipient site from GT of recipient

		my $organ_source_site = "";#orgran site in cfDNA mixture sample
		my $organ_source_site_dp = "";
		my $recipient_source_site = "";#recipient site in cfDNA mixture sample
		my $recipient_source_site_dp = "";
		my $differentiate_donor_site = "";#record donor site which can differentiate to recipient site
 		my $both_dp_added = '';
		#filter error-like site 
		#?????
		if(@b>1 || $b[0] ne ".")
		{
			my @bb=();
			my @dp_of_bb=();

			for(my $i=0;$i<@b;$i++)
			{
				#if base exists in population database, it may be true genotype
				if(exists $$rs{"$a[1]:$a[2]"}{$b[$i]})
				{
					push(@bb,$b[$i]);
					push(@dp_of_bb,$==$dp_of_b[$i]);
					#print $_,"\n";
				}elsif($dp_of_b[$i]>2)#if depth of one base >2, it may be true genotype
				{
					#push(@bb,$b[$i]);
					#push(@dp_of_bb,$dp_of_b[$i]);
				}
			}
			@b=();@dp_of_b=();
			@b = @bb;
			@dp_of_b = @dp_of_bb;
			if(@b==0)
			{
				push(@b,".");
				push(@dp_of_b,0);
				#print "==========$_\n";
			}
		}

		my $altsum = &getsum(\@dp_of_b);
		if(($altsum+$a[6])==0)
		{
			#usually appeared when genotype is 1/2,and depth equal to zero
			#print STDERR $_,"\n";
			#print STDERR "b=",@b,"\n";
			#print STDERR "dpofb=",@dp_of_b,"\n";
			next;
		}
		my $rf = $a[6] / ($a[6] + $altsum);
		
		if($has_recipient_data < 1)
		{

			next if($rf>=0.2 && $rf<=0.8);
			if($a[6]==0)#ref dp ==0 ,   shift the first base of alt
			{

				$recipient_source_site = shift @b;
				$recipient_source_site_dp = shift @dp_of_b;
				#===============================================
				#need more check whether snp really exists (dbsnp 1000g)
				#===============================================
				#$organ_source_site = join("," ,  @b);
				#$organ_source_site_dp = &getsum(\@dp_of_b);
				$organ_source_site = $a[3];#slient signal
				$organ_source_site_dp = $a[6];#slient signal
			}elsif($a[7] == 0)
			{
				$recipient_source_site = $a[3];
				$recipient_source_site_dp = $a[6];
				$organ_source_site = $a[4];
				$organ_source_site_dp = $a[7];
			}else{
				#===============================================
				#need more check whether snp really exists (dbsnp 1000g)
				#===============================================
				if($rf<=0.2)
				{
					$recipient_source_site = $b[0];
					$recipient_source_site_dp = $dp_of_b[0];
					$organ_source_site = $a[3];
					$differentiate_donor_site = $a[3];
					$organ_source_site_dp = $a[6];
				}else{
					$recipient_source_site = $a[3];
					$recipient_source_site_dp = $a[6];
					$organ_source_site = $b[0];
					$differentiate_donor_site = $b[0];
					$organ_source_site_dp = $dp_of_b[0];					
				}
			}
		}

		if($has_recipient_data > 0)
		{
			# for excluding of recipient's snp
			next unless(exists $$recipient_Snp{"$a[1]:$a[2]"});# || ($rf>0&&$rf<0.2 || $rf>0.8&&$rf<1))
			$recipient_sites = $$recipient_Snp{"$a[1]:$a[2]"};
			next if($recipient_sites=~/,/); #heterozigous

			if($recipient_sites  eq  $a[3])
			{
				$recipient_source_site = $a[3];
				$recipient_source_site_dp = $a[6];
			}else{
				$organ_source_site = $a[3];# if($a[6]);
				$differentiate_donor_site = $a[3];
				$organ_source_site_dp = $a[6];# if($a[6]);
			}
			for(my $i=0;$i<@b;$i++)
			{
				if($recipient_sites  eq  $b[$i])
				{
					$recipient_source_site .= ",$b[$i]";
					##orignal depth information
					#$recipient_source_site_dp .= ",$dp_of_b[$i]";
					##accumulation
					$recipient_source_site_dp += $dp_of_b[$i];
				}else{
					#need to be checked by 1000g vcf
					$organ_source_site .= ",$b[$i]";
					$differentiate_donor_site .= ",$b[$i]";
					#$organ_source_site_dp .= ",$dp_of_b[$i]";
					$organ_source_site_dp += $dp_of_b[$i];
				}
			}
		}

		if($has_donor_data > 0 )
		{
			#print "=>$_\n";
			
			#for including snp of donor's snp according to chromosome location
			next unless(exists $$Donor_Snp{"$a[1]:$a[2]"});
			$donor_sites = $$Donor_Snp{"$a[1]:$a[2]"};
			
			$organ_source_site = $donor_sites;
			$organ_source_site_dp = 0;
			$recipient_source_site = $recipient_sites;
			$recipient_source_site_dp = 0;
            
           
			
			if($donor_sites=~/,/ && $donor_sites=~/$recipient_sites/)
			{
				if($donor_sites eq $recipient_sites)
				{
					#it's not likely to differentiate them
					next;
				}else{
					#donor sites belong to heterozygous
					if($a[3] eq $recipient_sites)
					{
						$recipient_source_site_dp = $a[6];
					}elsif($donor_sites=~/$a[3]/){
						$organ_source_site_dp = $a[6];
						$differentiate_donor_site = $a[3];
					}
					for(my $i=0;$i<@b;$i++)
					{
						if($b[$i] eq $recipient_sites)
						{
							$recipient_source_site_dp = $dp_of_b[$i];
						}elsif($donor_sites=~/$b[$i]/)
						{
                            if($organ_source_site_dp > 0)
                            {
                                $both_dp_added = "TRUE";
                            }
							$differentiate_donor_site .= ",$b[$i]";
							$organ_source_site_dp += $dp_of_b[$i];
						}
					}   
				}  		
			}else{
                if($donor_sites eq $recipient_sites)
                {
                    next;
                }
				#sites from donor and recipient are totally different
				if($a[3] eq $recipient_sites)
				{
					$recipient_source_site_dp = $a[6];
				}elsif($a[3] eq $donor_sites){
					$organ_source_site_dp = $a[6];
					$differentiate_donor_site = $a[3];
				}elsif($donor_sites=~/$a[3]/)
                {
                    $organ_source_site_dp = $a[6];
					$differentiate_donor_site = $a[3];
                }
				for(my $i=0;$i<@b;$i++)
				{
					if($b[$i] eq $recipient_source_site)
					{
						$recipient_source_site_dp = $dp_of_b[$i];
					}elsif($donor_sites eq $b[$i])
					{
						$organ_source_site_dp = $dp_of_b[$i];
						$differentiate_donor_site = $b[$i];
					}elsif($donor_sites=~/$b[$i]/)#donor sites are heterozygous sites
					{
                        if($organ_source_site_dp > 0)
                        {
                            $both_dp_added = "TRUE";
                        }
						$differentiate_donor_site .= ",$b[$i]"; 
						$organ_source_site_dp += $dp_of_b[$i];
					}
				}
			}

		}
		$organ_source_site=~s/^,//;
		$organ_source_site_dp=~s/^,//;
		$recipient_source_site=~s/^,//;
		$recipient_source_site_dp=~s/^,//;
		$differentiate_donor_site=~s/^,//;
		#print "$donor_sites:$organ_source_site\t$organ_source_site_dp\t$recipient_sites:$recipient_source_site\t$recipient_source_site_dp\n\n";
		
		next if($organ_source_site eq "");#no organ source site
		next if($recipient_source_site eq "");#no recipient source site
		next if($recipient_source_site =~ /,/); # heterozigous genotype of recipient
		
		my $donor_rf = $organ_source_site_dp/($organ_source_site_dp+$recipient_source_site_dp);
		
		if($differentiate_donor_site=~/[ATCG]/ && $organ_source_site_dp>0)
		{
			#$$Organ_Snp_incfDNA{"$a[1]:$a[2]"} = $organ_source_site if($organ_source_site);
			$$Organ_Snp_incfDNA{"$a[1]:$a[2]"} = $differentiate_donor_site if($organ_source_site);
			#print STDERR "$a[1]:$a[2]\t$organ_source_site\t",$$Organ_Snp_incfDNA{"$a[1]:$a[2]"},"\t$_\n";
		}
		my $RRAF = ($a[3] eq $recipient_source_site) ? $a[-1] : (1-$a[-1]);
		my $d_hete_AF = 0;
		my $d_homo_AF = 0;
		if($donor_sites ne "") # has donor
		{
			if($donor_sites =~ /,/)#hete
			{ 
				$d_hete_AF = 1;
				$d_homo_AF = 0;
				if($both_dp_added eq "TRUE")
				{
					$organ_source_site_dp = int($organ_source_site_dp/2); #get mean depth
				}
			}else{#homo
				$d_hete_AF = 0;
				$d_homo_AF = 1;
			}
		}else{
			$d_hete_AF = 2 * $a[-1] * (1-$a[-1]);
			$d_homo_AF = (1-$a[-1])**2;
		}
		
		#print STDERR "$a[0]\t$a[1]\t$a[2]\t$recipient_sites\t$donor_sites\t$RRAF\t$d_hete_AF\t$d_homo_AF\n";
		
		print OO "$a[0]\t$a[1]\t$a[2]\t$a[3]\t$a[4]\t$a[5]\t$a[6]\t$a[7]";
		print OO "\t$recipient_source_site\t$recipient_source_site_dp\t$organ_source_site\t$organ_source_site_dp\t$donor_rf\t$a[-2]\t$RRAF\t$d_hete_AF\t$d_homo_AF\n";
	}
	close TT;
	close OO;
}

sub GT_from_TAB
{
	#like donor derived
	my $snp=shift;
	my $hash = shift;
	my $is_recipient=shift;
	open  GT , $snp or die $!;
	#print "$snp\n";
	while(<GT>)
	{
		chomp;
		my @a = split /\t/ , $_;
		next if(length($a[4])>1 && $a[4]!~/,/); #ignore indel
		next if($a[5] < 10 || ($a[6]+$a[7]) < 10); # ignore snp which total depth less than 3
		my $rf = $a[6] / ($a[6] + $a[7]);
		my $site = &getGTbase($a[10]);

		next if((/0\/1/ || /1\/2/) && $is_recipient);#move to the next site if GT is heterozygous
	
		if($rf > 0.2 && $rf<0.8)
		{
			#heterozigous

			unless($is_recipient)
			{
				$$hash{"$a[1]:$a[2]"} = $site;
			}
			
		}else{
			#homozigous
			#$site = ($rf<=0.2) ? $a[4] : $a[3];
			#$site=~s/\.//;
			$$hash{"$a[1]:$a[2]"} = $site;
		}
	}
	close GT;
}

sub getGTbase
{
	# return base[s] of genotype on SNP loci
	my $base=shift;
	my @a=split(/\|/, $base);
	return ($a[0] eq $a[1])?$a[0]:"$a[0],$a[1]";
}
sub MXddcfDNA_info_from_bam
{
	#extract reads from donor in recipient cfDNA data
	my $bam = shift; #bam format file
	my $snptab = shift; #filename of snp table
	my $effectiveSNP=shift; #filename of effective SNP table

	my %Organ_Snp=(); #reference of a hash
	my %donor_cfdna_length=(); #reference of a hash
	my %donor_short_long=(); #reference of a hash
	my %donor_AB=(); #reference of a hash
	my %recipient_cfdna_length=(); #reference of a hash
	my %recipient_short_long=(); #reference of a hash
	my %recipient_AB=(); #reference of a hash
	my %depth=(); #reference of a hash
	my $length_cutoff=140;

	&load_effective_donor_sites($effectiveSNP,\%Organ_Snp);
	
	open BAM , "/data/soft/bin/samtools view --thread 5 $bam|" or die $!;
	open RD ,  ">$snptab.ex.read" or die $!;

	while(<BAM>)
	{
		chomp;
		my $line = $_;
		my @a = split /\t/ , $_;
		next if($a[5] eq '*');
		next if($a[6] ne "=");

		my %tmp = ();
		my (@cigar) = $a[5] =~ /(\d+\D)/g;
		#my (@cigar) = $line =~ /MC:Z:(\d+\w)/g;
		my $AB = "";
		if($cigar[0]=~/(\d+)S$/)
		{
			shift @cigar;
			$AB = substr($a[9] , $1 , 2);
		}else{
			$AB = substr($a[9] , 0 , 2);
		}

		my $ref_point = $a[3];
		my $base_point = -1;
		for(my $i=0;$i<@cigar;$i++)
		{
			
			if($cigar[$i]=~/(\d+)M/)
			{
				for(my $j=0;$j<$1;$j++)
				{
					$base_point += 1;
					#$ref_point += $j;
					next unless(exists $Organ_Snp{"$a[2]:$ref_point"});
					$tmp{$ref_point} = substr($a[9] , $base_point , 1);
					#print "$a[2]\t$ref_point\t",$Organ_Snp{"$a[2]:$ref_point"},"=$tmp{$ref_point}\n";
				}
				$ref_point += $1;
			}elsif($cigar[$i]=~/(\d+)D/)
			{
				$ref_point += $1;
				#$tmp{$ref_point} = substr($ref , $ref_point , $1);
			}elsif($cigar[$i]=~/(\d+)I/)
			{
				#$tmp{$ref_point} = substr($a[9] , $base_point , $1);
				$base_point += $1;
			}elsif($cigar[$i]=~/(\d+)S/)
			{
				$ref_point += $1;
				$base_point += $1;
			}
		}

		#test whether it is a polymorphsm site which mapping with table
		foreach my $pos( keys %tmp)
		{
		#	if(exists $Organ_Snp{"$a[2]:$pos"})
		#	{
				#print STDERR $Organ_Snp{"$a[2]:$pos"},"\t",$tmp{$pos},"\n";
				if($Organ_Snp{"$a[2]:$pos"} =~ /$tmp{$pos}/)
				{
					#yes organ
					print RD $a[2] , "\t" , $pos , "\t" , $Organ_Snp{"$a[2]:$pos"} , "\t" , $tmp{$pos} , "\t" , $line , "\n";
					$donor_cfdna_length{abs($a[8])}++;

					if(abs($a[8])<$length_cutoff)
					{
						$donor_short_long{"$a[2]:$pos"}{'SHORT'} ++;
					}else{
						$donor_short_long{"$a[2]:$pos"}{'LONG'} ++;
					}
					$donor_AB{"$a[2]:$pos"}{$AB}++ if($AB !~ /N/);
					$depth{"$a[2]:$pos"}{'donor'}++;
				}else{
					$recipient_cfdna_length{abs($a[8])}++;
					if(abs($a[8])<$length_cutoff)
					{
						$recipient_short_long{"$a[2]:$pos"}{'SHORT'} ++;
					}else{
						$recipient_short_long{"$a[2]:$pos"}{'LONG'} ++;
					}
					$recipient_AB{"$a[2]:$pos"}{$AB}++ if($AB !~ /N/);
					$depth{"$a[2]:$pos"}{'recipient'}++;
				}
		#	}
		}
	}
	close BAM;
	close RD;

	&dump_cfDNA_length($snptab,\%donor_cfdna_length,\%recipient_cfdna_length);
	&dump_cfDNA_AB($snptab,$donor_AB,\%recipient_AB,\%donor_short_long,\%recipient_short_long,\%depth);
}

sub dump_cfDNA_length
{
	my ($snptab , $hd , $hr) = @_;
	open HOUT , ">$snptab.xd-cfDNA.insertsize" or die $!;
    print HOUT "Type\tLength\tCount\n";
	foreach my $k(sort {$a<=>$b} keys %{$hd})
	{
		next if($k>1000||$k==0);
		print HOUT "Donor\t$k\t$$hd{$k}\n";
	}
	foreach my $k(sort {$a<=>$b} keys %{$hr})
	{
		next if($k>1000||$k==0);
		print HOUT "Recipient\t$k\t$$hr{$k}\n";
	}
	close HOUT;
}

sub dump_cfDNA_AB
{
	my ($snptab , $dAB , $rAB , $donor_short_long , $recipient_short_long , $depth) = @_;
	my ($short , $long , $shortR , $longR)=(0 , 0 , 0 , 0);
	my %ABS = ();
	my ($title , $subtitle , $subtitle2) = ("" , "" , "");
	my @base = ("A" , "T" , "C" , "G");
	for my $a(@base)
	{
		for my $b(@base)
		{
			$ABS{"$a$b"} = 1 unless(exists $ABS{"$a$b"});
		}
	}
	open ABOUT , ">$snptab.ex.AB" or die $!;
	$title = "Chr\tPos\tR_DP\tD_DP\tR_Scfdna\tR_Lcfdna\tD_Scfdna\tD_Lcfdna";
	foreach my $AB(sort {$a cmp $b} keys %ABS)
	{	
		$subtitle .= "\tR_$AB";
		$subtitle2 .= "\tD_$AB";
	}
	$title .= $subtitle.$subtitle2."\n";

	print ABOUT $title;
	foreach my $site(keys %$recipient_short_long)
	{
		$short = (exists $$recipient_short_long{$site}{'SHORT'}) ? $$recipient_short_long{$site}{'SHORT'} : 0;
		$long = (exists $$recipient_short_long{$site}{'LONG'}) ? $$recipient_short_long{$site}{'LONG'} : 0;
		$shortR = (($short + $long) > 0) ? $short*100/($short+$long):0;
		$longR = (($short + $long) > 0) ? $long*100/($short+$long):0;
		my $dp_r = (exists $$depth{$site}{'recipient'}) ? $$depth{$site}{'recipient'} : 0;
		my $dp_d = (exists $$depth{$site}{'donor'}) ? $$depth{$site}{'donor'} : 0; 
		my $site_out = $site; $site_out=~s/;/\t/;
		print ABOUT "$site_out";
		print ABOUT "\t$dp_r";
		print ABOUT "\t$dp_d";
		print ABOUT "\t$shortR\t$longR";

		$short = (exists $$donor_short_long{$site}{'SHORT'}) ? $$donor_short_long{$site}{'SHORT'} : 0;
		$long = (exists $$donor_short_long{$site}{'LONG'}) ? $$donor_short_long{$site}{'LONG'} : 0;
		$shortR = (($short + $long) > 0) ? $short*100/($short+$long):0;
		$longR = (($short + $long) > 0) ? $long*100/($short+$long):0;
		print ABOUT "\t$shortR\t$longR";

		my $sum_donor = 0;
		my $sum_recipient=0;
		foreach my $AB(keys %ABS)
		{
			$sum_recipient += $$rAB{$site}{$AB} if(exists $$rAB{$site}{$AB});
			$sum_donor += $$dAB{$site}{$AB} if(exists $$dAB{$site}{$AB}); 
		}

		foreach my $AB(keys %ABS)
		{
			if(exists $$rAB{$site}{$AB})
			{
				print ABOUT "\t" , $$rAB{$site}{$AB}/$sum_recipient*100;
			}else{
				print ABOUT "\t0";
			}

			if(exists $$dAB{$site}{$AB})
			{
				print ABOUT "\t" , $$dAB{$site}{$AB}/$sum_donor*100;
			}else{
				print ABOUT "\t0";
			}
		}
		print ABOUT "\n";
		delete $$donor_short_long{$site} if(exists $$donor_short_long{$site});
		delete $$recipient_short_long{$site} if(exists $$recipient_short_long{$site});
	}

	foreach my $site(keys %$donor_short_long)
	{
		$short = (exists $$recipient_short_long{$site}{'SHORT'}) ? $$recipient_short_long{$site}{'SHORT'} : 0;
		$long = (exists $$recipient_short_long{$site}{'LONG'}) ? $$recipient_short_long{$site}{'LONG'} : 0;
		$shortR = (($short + $long) > 0) ? $short*100/($short+$long):0;
		$longR = (($short + $long) > 0) ? $long*100/($short+$long):0;
		my $dp_r = (exists $$depth{$site}{'recipient'}) ? $$depth{$site}{'recipient'} : 0;
		my $dp_d = (exists $$depth{$site}{'donor'}) ? $$depth{$site}{'donor'} : 0; 
		print ABOUT "$site";
		print ABOUT "\t$dp_r";
		print ABOUT "\t$dp_d";
		print ABOUT "\t$shortR\t$longR";

		$short = (exists $$donor_short_long{$site}{'SHORT'}) ? $$donor_short_long{$site}{'SHORT'} : 0;
		$long = (exists $$donor_short_long{$site}{'LONG'}) ? $$donor_short_long{$site}{'LONG'} : 0;
		$shortR = (($short + $long) > 0) ? $short*100/($short+$long):0;
		$longR = (($short + $long) > 0) ? $long*100/($short+$long):0;
		print ABOUT "\t$shortR\t$longR";

		my $sum_donor = 0;
		my $sum_recipient=0;
		foreach my $AB(keys %ABS)
		{
			$sum_recipient += $$rAB{$site}{$AB} if(exists $$rAB{$site}{$AB});
			$sum_donor += $$dAB{$site}{$AB} if(exists $$dAB{$site}{$AB}); 
		}

		foreach my $AB(keys %ABS)
		{
			if(exists $$rAB{$site}{$AB})
			{
				print ABOUT "\t" , $$rAB{$site}{$AB}/$sum_recipient*100;
			}else{
				print ABOUT "\t0";
			}

			if(exists $$dAB{$site}{$AB})
			{
				print ABOUT "\t" , $$dAB{$site}{$AB}/$sum_donor*100;
			}else{
				print ABOUT "\t0";
			}
		}
		print ABOUT "\n";
		delete $$donor_short_long{$site} if(exists $$donor_short_long{$site});
		delete $$recipient_short_long{$site} if(exists $$recipient_short_long{$site});
	}

	close ABOUT;

}

sub getsum
{
	my $arr = shift;
	my $sum = 0;
	for(my $i=0;$i<@$arr;$i++)
	{
		$sum += $$arr[$i];
	}
	return $sum;
}

sub load_effective_donor_sites
{
	my($esnp,$hash)=@_;
	open ESNP,$esnp or die $!;
	while(<ESNP>)
	{
		next if(/^Sample/);
        my @a = split;
        $$hash{"$a[1]:$a[2]"} = $a[10];
	}
	close ESNP;
}

1;
