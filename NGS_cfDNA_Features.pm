package RRG::NGS_TOOLS;

use File::Basename;

require Exporter;

@ISA = qw(Exporter);

@EXPORT = qw();
#@EXPORT_OK = qw();

sub genotype_sites
{
	my ($snptab,$donor,$recipient)=@_;

	my %Organ_Snp = ();
	my %Donor_Snp = ();
	my %Recipient_Snp = ();
	my %donor_cfdna_length = ();
	my %recipient_cfdna_length = ();
	my %donor_short_long = ();
	my %recipient_short_long = ();
	my %donor_AB = ();
	my %recipient_AB = ();
	my %depth = ();

	&include($donor) if($donor ne "");
	&exclude($recipient) if($recipient ne "");
	&load_snp($snptab);
}

sub dump
{
	&extract_info_from_bam($bam,$snptab);
	&dump_cfDNA_length(\%donor_cfdna_length,\%recipient_cfdna_length);
	&dump_cfDNA_AB(\%donor_AB,\%recipient_AB,\%donor_short_long,\%recipient_short_long,\%depth);
}

sub load_snp
{
	my $snp=shift;
	open TT,$snp or die $!;
	while(<TT>)
	{
		chomp;
		my @a = split /\t/,$_;
		next if($a[5]<10 || $a[4] eq ".");
		my (@b) = $a[4] =~ /\w+/g;
		#my $flag = 0;
		#for(my $i=0;$i<@b;$i++)
		#{
		#		$flag = 1 if(length($b[$i])>1);
		#}
		#next if($flag);

		my $organ_source_site = "";
		if($recipient eq "")
		{
			my $rf = $a[6] / ($a[6] + $a[7]);
			next if($rf>0.2 && $rf <0.8);
			if($a[6]==0)
			{
				my @b = split /,/,$a[4];
				shift @b;
				$organ_source_site = join(",",@b);
			}else{
				$organ_source_site = ($rf<=0.2) ? $a[3] : $a[4];
			}
			print "$a[1]\t$a[2]\t$a[3]\t$a[4]\t$a[5]\t$a[6]\t$a[7]\n";
		}

		if($recipient ne "")
		{
			# for excluding of recipient's snp
			
			next unless(exists $Recipient_Snp{"$a[1]:$a[2]"});
			print "$a[1]\t$a[2]\t$a[3]\t$a[4]\t$a[5]\t$a[6]\t$a[7]\trecipient=",$Recipient_Snp{"$a[1]:$a[2]"},"\t";
			if($Recipient_Snp{"$a[1]:$a[2]"}!~/$a[3]/)
			{
				$organ_source_site = $a[3];
				print "\tOrgan=$a[3]";
			}
			for(my $i=0;$i<@b;$i++)
			{
				if($Recipient_Snp{"$a[1]:$a[2]"}!~/$b[$i]/)
				{
					$organ_source_site .= ",$b[$i]";
					print "\tOrgan=$b[$i]";
				}
			}
			print "\n";
		}

		if($donor ne "")
		{
			#for including snp of donor's snp
			next unless(exists $Donor_Snp{"$a[1]:$a[2]"});
			print "$a[1]\t$a[2]\t$a[3]\t$a[4]\t$a[5]\t$a[6]\t$a[7]\tdonor=",$Donor_Snp{"$a[1]:$a[2]"},"\n";
			if($a[3] =~ /$Donor_Snp{"$a[1]:$a[2]"}/)
			{
				$organ_source_site = "$a[3]";
			}	
			for(my $i=0;$i<@b;$i++)
			{
				if($b[$i]=~/$Donor_Snp{"$a[1]:$a[2]"}/)
				{
					$organ_source_site .= ",$b[$i]";
				}
			}
		}
		$organ_source_site=~s/^,//;
		$Organ_Snp{"$a[1]:$a[2]"} = $organ_source_site if($organ_source_site);
		#print STDERR $_,"\n";
	}
	close TT;
}

sub include
{
	#like donor derived
	my $include_snp=shift;
	open ICDSNP,$include_snp or die $!;
	while(<ICDSNP>)
	{
		chomp;
		my @a = split /\t/,$_;
		next if($a[5] == 0);
		my $rf = $a[6] / ($a[6] + $a[7]);
		my $donor_source_site = "";
		if($rf > 0.2 && $rf<0.8)
		{
			$donor_source_site = "$a[3],$a[4]";
		}else{
			$donor_source_site = ($rf<=0.2) ? $a[4] : $a[3];
		}
		$Donor_Snp{"$a[1]:$a[2]"} = $donor_source_site;				
	}
	close ICDSNP;
}

sub exclude
{
	#recipient's genotype
	my $exclude_snp=shift;
	open EXDSNP,$exclude_snp or die $!;
	while(<EXDSNP>)
	{
		chomp;
		my @a = split /\t/,$_;
		next if($a[5] < 1);
		my $rf = $a[6] / ($a[6] + $a[7]);
		my $recipient_source_site = "";
		if($rf > 0.2 && $rf<0.8)
		{
			$recipient_source_site = "$a[3],$a[4]";
		}else{
			$recipient_source_site = ($rf<=0.2) ? $a[4] : $a[3];
		}

		$Recipient_Snp{"$a[1]:$a[2]"} = $recipient_source_site;
	}
	close EXDSNP;
}
sub extract_info_from_bam
{
	my $bam = shift;
	my $snptab = shift;
	my $length_cutoff=140;
	open BAM,"/data/soft/bin/samtools view --thread 2 $bam|" or die $!;
	open READ,">$snptab.read" or die $!;

	while(<BAM>)
	{
		chomp;
		my $line = $_;
		&searching($line,$length_cutoff);
	}

	close BAM;
	close READ;
}

sub searching
{
	my $line = shift;
	my $length_cutoff=shift;

	my @a = split /\t/,$_;
	return 0 if($a[5] eq '*');
	return 0 if($a[6] ne "=");
	my %tmp = ();
	my (@cigar) = $a[5] =~ /(\d+\D)/g;
	#my (@cigar) = $line =~ /MC:Z:(\d+\w)/g;
	my $AB = "";
	if($cigar[0] eq "")
	{
		print STDERR "$line\n";
	}
	if($cigar[0]=~/(\d+)S$/)
	{
		shift @cigar;
		$AB = substr($a[9],$1,2);
	}else{
		$AB = substr($a[9],0,2);
	}
	#print STDERR join("-",@cigar),"\n";
	my $ref_point = $a[3];
	my $base_point = -1;
	for(my $i=0;$i<@cigar;$i++)
	{
		if($cigar[$i]=~/(\d+)M/)
		{
			for(my $j=0;$j<$1;$j++)
			{
				$base_point += 1;
				#print STDERR join(":",@cigar),"\ti=$i\tj=$j\tbase_point=$base_point\t$a[9]\n";
				$tmp{$ref_point+$j} = substr($a[9],$base_point,1);
			}
			$ref_point += $1;
		}elsif($cigar[$i]=~/(\d+)D/)
		{
			$ref_point += $1;
			#$tmp{$ref_point} = substr($ref,$ref_point,$1);
		}elsif($cigar[$i]=~/(\d+)I/)
		{
			#$tmp{$ref_point} = substr($a[9],$base_point,$1);
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
		if(exists $Organ_Snp{"$a[2]:$pos"})
		{
			if($Organ_Snp{"$a[2]:$pos"} =~ /$tmp{$pos}/)
			{
				#yes organ
				print READ $a[2],"\t",$pos,"\t",$Organ_Snp{"$a[2]:$pos"},"\t",$tmp{$pos},"\t",$line,"\n";
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
				#print STDERR "========$a[5]==============\n";
			}
		}
	}
}

sub dump_cfDNA_length
{
	my ($hd,$hr) = @_;
	open HDOUT,">$snptab.donor_cfDNA" or die $!;
	open HROUT,">$snptab.recipient_cfDNA" or die $!;
	foreach my $k(sort {$a<=>$b} keys %{$hd})
	{
		next if($k>1000||$k==0);
		print HDOUT "$k\t$$hd{$k}\n";
	}
	foreach my $k(sort {$a<=>$b} keys %{$hr})
	{
		next if($k>1000||$k==0);
		print HROUT "$k\t$$hr{$k}\n";
	}
	close HDOUT;
	close HROUT;
}
sub dump_cfDNA_AB
{
	my ($dAB,$rAB,$donor_short_long,$recipient_short_long,$depth) = @_;
	my ($short,$long,$shortR,$longR)=(0,0,0,0);
	my %ABS = ();
	my ($title,$subtitle,$subtitle2) = ("","","");
	my @base = ("A","T","C","G");
	for my $a(@base)
	{
		for my $b(@base)
		{
			$ABS{"$a$b"} = 1 unless(exists $ABS{"$a$b"});
		}
	}
	open ABOUT,">$snptab.AB" or die $!;
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
				print ABOUT "\t",$$rAB{$site}{$AB}/$sum_recipient*100;
			}else{
				print ABOUT "\t0";
			}
			
			if(exists $$dAB{$site}{$AB})
			{
				print ABOUT "\t",$$dAB{$site}{$AB}/$sum_donor*100;
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
				print ABOUT "\t",$$rAB{$site}{$AB}/$sum_recipient*100;
			}else{
				print ABOUT "\t0";
			}
			
			if(exists $$dAB{$site}{$AB})
			{
				print ABOUT "\t",$$dAB{$site}{$AB}/$sum_donor*100;
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

sub usage
{
	print "This software is apt to extract read from a alignment by input genotype as baits  .\n";
	print "Usage:\n";
	print "\t--bam\t<file>\talignment with bam format[required]\n";
	print "\t--snptab\t<file>\tsnp genotype of mixed sample[required]\n";
	print "\t--donor\t<file>\tsnp genotype of donor sample, for including\n";
	print "\t--recipient\t<file>\tsnp genotype of recipient sample,for excluding\n";
	print "\t--help\t\tprint help information\n";
	exit;
}
