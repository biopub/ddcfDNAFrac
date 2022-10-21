use strict;
use warnings;

use SOT_app;
use NGS_PIP;


my $vcf=shift;
my $bed=shift;
my $output=shift;

my %rs = ();

#load_config($conf);


load_rs_from_bed($bed,\%rs);


VCF_to_TAB($vcf,\%rs,$output);

