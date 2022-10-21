package NGS_TOOLS;

use File::Basename;

require Exporter;

@ISA = qw(Exporter);

@EXPORT = qw(run select_micro_read rev_com cal_barcode_read_len ossRR set_datelike_prefix);
#@EXPORT_OK = qw();

my $MAX_TASKS_ONLINE_PER_USER = 200;

sub set_datelike_prefix
{
	my @months = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );
	my @days = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
	$mon+=1;
	my $output_prefix="SOT_$year\_$mon\_$mday\_$days[$wday]";
	return $output_prefix;
}
sub ossRR
{
	my $bucket=shift;
	my ($ossBucketP,$ossAccount) = ("","");
	if($bucket eq "rrgene")
	{
		$ossBucketP = "oss://rrgene";
		$ossAccount = "--endpoint=oss-cn-hangzhou.aliyuncs.com --access-key-id=LTAIcFJhm0HfkDHm --access-key-secret=QTYxS3xB2p711xfjE20hMRYbrCg5Mt";
	}elsif($bucket eq "rrgene2")
	{
		$ossBucketP = "oss://rrgene2";
		$ossAccount = "--endpoint=oss-cn-hangzhou.aliyuncs.com --access-key-id=LTAI926MxhIXzDHn --access-key-secret=5KSFJlsRqzWwWKJ6O0KMOAhv24KISE";
	}elsif($bucket eq "rrgene3")
	{
		$ossBucketP = "oss://rrgene3";
		$ossAccount = "--endpoint=oss-cn-hangzhou.aliyuncs.com --access-key-id=LTAIeboEMrBA9lYT --access-key-secret=fxvvaXqC7SUpxvl8YYcrVcXKYnUpRD";
	}elsif($bucket eq "allodx")
	{
		$ossBucketP = "oss://allodx";
		$ossAccount = "--endpoint=oss-cn-hangzhou.aliyuncs.com --access-key-id=LTAIZ2bznMdYr75H --access-key-secret=Iq6V9y0gk4Xie2OnHyZbYoznW0qeb6";
	}elsif($bucket eq "nextcode")
	{
		#$ossBucketP = "oss://delivery-data/s836/";
		#$ossAccount = "--endpoint=oss-cn-shanghai.aliyuncs.com --access-key-id=LTAI8FmldxfwkZxS --access-key-secret=wDamc7f29HPJYcndNBZh0QVgwMP0Zd";
		$ossBucketP = "oss://delivery-data/s1293/";
		$ossAccount = "--endpoint=oss-cn-shanghai.aliyuncs.com --access-key-id=LTAIqyIWO0IE15k1 --access-key-secret=c0VAxj7AtMvoSgGASJt3R8znm49Fay";
	}
	return ($ossBucketP,$ossAccount);
}
sub cal_barcode_read_len
{
	my $read=shift;
	my $spacer=shift;
	my $blen=0;
	my $rlen = 0;
	my %h=();
	my $lines = 0;
	open IN,"gzip -dc $read|" or die $!;
	while(<IN>)
	{
		my $seq = <IN>;chomp;
		<IN>;#sign
		<IN>;#quality
		$rlen=length($seq);
		$blen=index($seq,$spacer);
		$h{$blen}++;
		$lines++;
		if($lines == 400)
		{
			chomp($seq);
			$rlen = length($seq);
			last;
		}
	}
	close IN;
	my $max = 0;
	foreach my $k(keys %h)
	{
		if($h{$k}>$max)
		{
			$blen = $k;
			$max = $h{$k};
		}
	}
	$rlen = $rlen - $blen - length($spacer);
	$repFilt=$blen * 2 + 1;
	#print "rlen:$rlen\nblen:$blen\nrepFilt:$repFilt\n";
	return ($rlen,$blen,$repFilt);
}

sub run
{
	#input <sample.sh>  <dir>
	my $script = shift;
	my $project = shift;
	my $dir = dirname $script;
	print STDERR "run $script\n";
	if(! -e "$script.statu")
	{
		system("touch $script.statu");
	}
	my $task_running_count = 0;
	$task_running_count = `qstat -u $ENV{'USER'}|wc -l`;
	while($task_running_count > $MAX_TASKS_ONLINE_PER_USER)
	{
		sleep(300); #sleep 30 min to count
		$task_running_count = `qstat -u $ENV{'USER'}|wc -l`;
	}

	#entry dir
	chdir($dir);

	#load statu
	my %task_statu = ();
	if( -f "$script.statu" )
	{
		open SST,"$script.statu" or die $!;
		while(<SST>)
		{
			chomp;
			my @a = split;
			$task_statu{$a[0]} = $a[1];
		}
		close SST;
	}

	open SS,$script or die $!;
	my @steps = ();
	my @task = ();
	my %job_in_queue=();
	my $is_last_line_empty = 1;
	while(<SS>)
	{
		if(/#(\d+)/)
		{
			$is_last_line_empty = 0;
			chomp;
			@task = split;
			(@steps) = $_ =~ /#(\d+)/g;
			open SUBTASK,">$script.$task[-1].$steps[-1]" or die $!;
			print SUBTASK $_,"\n";

		}elsif(/^$/ || eof)
		{
			next if($is_last_line_empty);
			$is_last_line_empty = 1;
			print SUBTASK "$_\n";
			print SUBTASK "if [ ! \$\? -eq 0 ]\nthen\n\texit\nelse\n";
			print SUBTASK "\techo $task[-1] done >> $script.statu\nfi\n";
			close SUBTASK;
			
			my $hold_jobid = "";
			
			for(my $i=0;$i<$#steps;$i++)
			{
				next if($task_statu{$task[-1]} eq "done");
				$hold_jobid .= "$job_in_queue{$steps[$i]},";
			}
			$hold_jobid =~ s/,$//;

			#make qsub command
			my $cmd = "qsub -cwd -V -l vf=$task[-2] -P $project";
			$cmd .= " -q big.q" if($task[-1] =~ /kraken/i);
			$cmd .= ($hold_jobid ne "") ? " -hold_jid $hold_jobid" : "";
			$cmd .= " $script.$task[-1].$steps[-1]";

			#excution
			if(!exists $task_statu{$task[-1]} || $task_statu{$task[-1]} ne "done")
			{
				#clear unfinished task run log
				my @job_run_log_file = glob("$script.$task[-1].$steps[-1].*");
				for(my $i=0;$i<@job_run_log_file;$i++)
				{
					if( -f $job_run_log_file[$i] )
					{
						unlink($job_run_log_file[$i]);
					}
				}

				#excution
				#my $info = `$cmd`;
				my ($jobid) = $info =~ /Your\s+job\s+(\d+)/;
				print STDERR "Step$steps[-1]:$task[-1]\t\t\tstep:$steps[-1]\t<<--\t";
				print STDERR (@steps==1)?"start\t\t":"steps:@steps[0..$#steps-1]\t\t";
				print STDERR (exists $task_statu{$task[-1]})?"$task_statu{$task[-1]}":"qsubed";
				print STDERR "\t$jobid\n";
				$job_in_queue{$steps[-1]} = $jobid;
			}
			@steps = ();
			@task = ();

		}else{
			print SUBTASK $_;
			$is_last_line_empty = 0;
		}
	}
	close SS;
	print STDERR "\n\n";
}


sub rev_com 
{
	my ($seq) = @_;
	$seq = reverse($seq);
	$seq =~ tr/ACGT/TGCA/;
	return($seq);
}

