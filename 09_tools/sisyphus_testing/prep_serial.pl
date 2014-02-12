#!/usr/bin/perl -w 

$top_path  = "/home/ivanam";
$tfm_table = "sysphus_tfms.csv";
$pdbdir    = "$top_path/databases/pdbfiles";
$pdbdown   = "$top_path/perlscr/downloading/pdbdownload.pl";
$pdtfm     = "$top_path/perlscr/pdb_manip/pdb_affine_tfm.pl";
$extr      = "$top_path/perlscr/pdb_manip/pdb_extract_chain.pl";
$pdb2seq   = "$top_path/perlscr/pdb_manip/pdb2seq.pl";
$afa2msf   = "$top_path/perlscr/translation/afa2msf.pl";
$nossa     = "$top_path/kode/03_struct/09_tools/no_seq_sim_struct_superp_analyzer/nossa";



$struct      = "$top_path/kode/03_struct/struct"; 

foreach ($top_path, $tfm_table, $pdbdir, $pdbdown, "params",
	 $pdtfm, $extr, $struct,$pdb2seq, $afa2msf, $nossa) {
    (-e $_) || die "$_ not found.\n";
}

foreach ("fold", "homologous", "fragment") {
    (-e $_) || `mkdir $_`;
}





($alignment_id, $pdb_code, $pdb_chain, $mat11, $mat12, 
 $mat13, $mat21, $mat22, $mat23, $mat31, $mat32, $mat33, 
 $shift1, $shift2, $shift3) = ();


open ( LOG, ">struct.log")  ||
    die "Cno struct.log: $!\n";

open ( IF, "<$tfm_table") ||
    die "Cno $tfm_table: $!\n";

$home = `pwd`; chomp $home;

$qryfile = "";


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
	#$ctr++;
	#($ctr==11) && exit;
    }

    print "\n\n\n##############################################\n";
    print "$pdb_code\n";


    chdir $home;

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
    
    if (! -e $chainfile || -z $chainfile) {
	# extract chain

	$cmd = "$extr $pdbdir/$pdb_code.pdb $pdb_chain >  $chainfile";
	if (system $cmd) {
	    print LOG "Error running $cmd.\n";
	}
	if (! -e $chainfile || -z $chainfile) {
	    die " could not create $chainfile.\n";
	}
   }
   


    if ($is_query) {
	$qryfile = $chainfile;
	next;
    }


    # write tfm - note the negative T
    $filename = "$pdb_code.to_$current_qry.tfm";
    open (OF, ">sys_tfms/$filename") ||
	die "Cno $filename: $!.\n";
    printf OF "%8.4f   %8.4f   %8.4f   %8.4f \n",  
    $mat11, $mat12, $mat13, -$shift1;
    printf OF "%8.4f   %8.4f   %8.4f   %8.4f \n",  
    $mat21, $mat22, $mat23, -$shift2;
    printf OF "%8.4f   %8.4f   %8.4f   %8.4f \n",  
    $mat31, $mat32, $mat33, -$shift3;
    close OF;    

    ###########################################################
    ###########################################################
    ###########################################################
    # tfm according to sysiphus
    $sys_chainfile_renamed = "pdbchains/$current_qry/$current_qry.to_$pdb_code$pdb_chain.sys.pdb";
    if ( ! -e $sys_chainfile_renamed ) {
	$sys_chainfile_orig    = "$current_qry.to_$pdb_code$pdb_chain.sys.pdb";
	# pdb transform using Sysiphus matrix - note translation first
	$cmd = "$pdtfm $qryfile sys_tfms/$filename -trnsl_first > $sys_chainfile_orig\n";
	if (system $cmd) {
	    print LOG "Error running $cmd.\n";
	}
	`mv $sys_chainfile_orig $sys_chainfile_renamed`;
    }


}

close IF;


