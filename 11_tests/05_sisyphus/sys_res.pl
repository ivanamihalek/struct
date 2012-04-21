#!/usr/bin/perl -w 

$tfm_table = "sysphus_tfms.csv";
$pdbdir    = "/Users/ivana/databases/pdbfiles";
$pdbdown   = "/Users/ivana/perlscr/downloading/pdbdownload.pl";
$pdtfm     = "/Users/ivana/perlscr/pdb_manip/pdb_affine_tfm.pl";
$extr      = "/Users/ivana/perlscr/pdb_manip/pdb_extract_chain.pl";
$struct    = "/Users/ivana/kode/03_struct/struct"; 
foreach ($tfm_table, $pdbdir, $pdbdown, "fold", "homologous", "fragment", "params",
	 $pdtfm, $extr, $struct) {
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

$qryfile = "";

$ctr = 0;

while (<IF>) {

    $ctr++;
    #($ctr==10) && exit;

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
    (-e "pdbchains/$current_qry") || `mkdir pdbchains/$current_qry`;
    (-e "sys_tfms")  || `mkdir sys_tfms`;


    $chainfile     = "pdbchains/$current_qry/$pdb_code$pdb_chain.pdb";
    
    if (!-e $chainfile) {
	# extract chain

	    $cmd = "$extr $pdbdir/$pdb_code.pdb $pdb_chain >  $chainfile";
	    if (system $cmd) {
		print LOG "Error running $cmd.\n";
	    }

    }
   
  

    if ($is_query) {
	$qryfile = $chainfile;
	next;
    }


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
    #$cmd = "$pdtfm $chainfile sys_tfms/$filename > $rot_chainfile\n";
    #if (system $cmd) {
#	print LOG "Error running $cmd.\n";
    #}

    $cmd = "$struct -in1 $chainfile -in2 $qryfile -p ../params";
    if  (system $cmd ) {
	print LOG "Error running $cmd.\n";
	next;
    }

    $rot_chainfile_orig    = "$pdb_code$pdb_chain.rot_onto_$current_qry.pdb";
    $rot_chainfile_renamed = "pdbchains/$current_qry/$pdb_code$pdb_chain.to_$current_qry.struct.pdb";

    `mv $rot_chainfile_orig $rot_chainfile_renamed`;
    `rm *.struct_out*`;

 
}

close IF;


