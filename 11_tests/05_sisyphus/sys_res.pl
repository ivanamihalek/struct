#!/usr/bin/perl -w 

$tfm_table = "sysphus_tfms.csv";
$pdbdir    = "/home/ivanam/databases/pdbfiles";
$pdbdown   = "/home/ivanam/perlscr/downloading/pdbdownload.pl";
$pdtfm     = "/home/ivanam/perlscr/pdb_manip/pdb_affine_tfm.pl";
$extr      = "/home/ivanam/perlscr/pdb_manip/pdb_extract_chain.pl";

foreach ($tfm_table, $pdbdir, $pdbdown, "fold", "homologous", "fragment") {
    (-e $_) || die "$_ not found.\n";
}

($alignment_id, $pdb_code, $pdb_chain, $mat11, $mat12, 
 $mat13, $mat21, $mat22, $mat23, $mat31, $mat32, $mat33, 
 $shift1, $shift2, $shift3) = ();


open ( LOG, ">log")  ||
    die "Cno log: $!\n";

open ( IF, "<$tfm_table") ||
    die "Cno $tfm_table: $!\n";

$home = `pwd`; chomp $home;


while (<IF>) {


    ($alignment_id, $alig_type, $pdb_code, $pdb_chain, $mat11, $mat12, 
     $mat13, $mat21, $mat22, $mat23, $mat31, $mat32, $mat33, 
     $shift1, $shift2, $shift3) = split ();

    next if ($alignment_id =~ /alignment_id/);

    $pdb_code =~ s/\"//g;
    $pdb_chain =~ s/\"//g;
    $alig_type =~ s/\"//g;


    if ( /1 0 0 0 1 0 0 0 1 0 0 0/ )  {
	$is_query    = 1;
	$current_qry = "$pdb_code$pdb_chain";
    } else {
	$is_query    = 0;
    }

    print "$pdb_code\n";

    if ( ! -e  "$pdbdir/$pdb_code.pdb") {
	$cmd = "$pdbdown $pdb_code\n";
	system $cmd;
    }

    if (! -e "$pdbdir/$pdb_code.pdb") {
	print LOG "$pdb_code.pdb not found\n";
	next;
    }
    


    chdir "$home/$alig_type";

    (-e "pdbchains") || `mkdir pdbchains`;
    (-e "sys_tfms")   || `mkdir sys_tfms`;


    $chainfile     = "pdbchains/$pdb_code$pdb_chain.pdb";
    next if (-e $chainfile);

    if ( $chain ) {

	$cmd = "$extr $pdbdir/$pdb_code.pdb $pdb_chain >  $chainfile";
	if (system $cmd) {
	    print LOG "Error running $cmd.\n";
	}

    } else {
	
	`cp $pdbdir/$pdb_code.pdb $chainfile`;
    }
   
  

    # extract chain
    next if ($is_query);

    $rot_chainfile = "pdbchains/$pdb_code$pdb_chain.to_$current_qry.pdb";

    # write tfm
    $filename = "$pdb_code.to_$current_qry.tfm";
    open (OF, ">sys_tfms/$filename") ||
	die "Cno $filename: $!.\n";
    printf OF "%8.4f   %8.4f   %8.4f   %8.4f \n",  
    $mat11, $mat12, $mat13, $shift1;
    printf OF "%8.4f   %8.4f   %8.4f   %8.4f \n",  
    $mat21, $mat22, $mat23, $shift2;
    printf OF "%8.4f   %8.4f   %8.4f   %8.4f \n",  
    $mat31, $mat32, $mat33, $shift3;
    close OF;    

    # pdb transform
    $cmd = "$pdtfm $chainfile sys_tfms/$filename > $rot_chainfile\n";
    if (system $cmd) {
	print LOG "Error running $cmd.\n";
    }

    
   print `pwd`;
    print $chainfile , " $alig_type \n";
    exit;

}

close IF;


