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
$fake_atom_type = "C";
while ( <IF> ) {
    next if ( ! /\S/ );
    next if (  /name/ );
    next if (  /residues/ );
    next if (  /helices/ );
    next if (  /strands/ );
    next if (  /^\#/ );
    @aux = split;
    @p   = @aux [5 .. 7];
    @cm  = @aux [8 .. 10];
    $chain = chr($ascii);
    if ( $aux[1] == 2 ) {
	$fake_atom_type = "C";
    } else {
	$fake_atom_type = "O";
    }

    #printf " %d \n\n", 2*$max+1;
    $ctr = 0;
    for ( $step = -$max;  $step <= $max; $step++ ) {
	$x = $cm[0]+$step*$cc_bond_length*$p[0];
	$y = $cm[1]+$step*$cc_bond_length*$p[1];
	$z = $cm[2]+$step*$cc_bond_length*$p[2];
	$ctr++;
	$crap = sprintf  "ATOM  %5d  $fake_atom_type   UNK $chain   1", 10000+$ctr;
	printf "%-30s%8.3f%8.3f%8.3f \n",
	$crap,  $x, $y, $z;
    }
    print "TER\n";
    $ascii--;
}


close IF;
