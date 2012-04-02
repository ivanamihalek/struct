#! /usr/bin/perl -w

# infile is a list
# of the form
#  \@p \@center_of_mass
#  \@p \@center_of_mass
#  .....

$cc_bond_length = 1.54;
$max = 10;


(@ARGV) ||
    die "Usage: carbon_line.pl   <p and cm file>  [<half length>]\n";

$infile = $ARGV[0];
open (IF, "<$infile") ||
    die "Cno $infile: $!.\n";
( defined $ARGV[1] ) && ($max = $ARGV[1]);


$ascii = ord ("Z");

while ( <IF> ) {
    next if ( ! /\S/ );
    @aux =split;
    @p = @aux [0 .. 2];
    @cm = @aux [3 .. 5];
    $chain = chr($ascii);

    #printf " %d \n\n", 2*$max+1;
    $ctr = 0;
    for ( $step = -$max;  $step <= $max; $step++ ) {
	$x = $cm[0]+$step*$cc_bond_length*$p[0];
	$y = $cm[1]+$step*$cc_bond_length*$p[1];
	$z = $cm[2]+$step*$cc_bond_length*$p[2];
	$ctr++;
	$crap = sprintf  "ATOM  %5d  C   UNK $chain   1", 10000+$ctr;
	printf "%-30s%8.3f%8.3f%8.3f \n",
	$crap,  $x, $y, $z;
    }
    print "TER\n";
    $ascii--;
}


close IF;
