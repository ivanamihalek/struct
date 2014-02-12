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



$mustang     = "$top_path/Downloads/MUSTANG_v3.2.1/bin/mustang-3.2.1";
$extr_chain  = "$top_path/perlscr/pdb_manip/pdb_extract_chain.pl";



foreach ($top_path, $tfm_table, $pdbdir, $pdbdown, "params",
	 $pdtfm, $extr, $pdb2seq, $afa2msf, $nossa, $mustang, 
	 $extr_chain) {
    (-e $_) || die "$_ not found.\n";
}

foreach ("fold", "homologous", "fragment") {
    (-e $_) || `mkdir $_`;
}



($alignment_id, $pdb_code, $pdb_chain, $mat11, $mat12, 
 $mat13, $mat21, $mat22, $mat23, $mat31, $mat32, $mat33, 
 $shift1, $shift2, $shift3) = ();


open ( LOG, ">mustang.log")  ||
    die "Cno mustang.log: $!\n";

open ( IF, "<$tfm_table") ||
    die "Cno $tfm_table: $!\n";

$home = `pwd`; chomp $home;

$qryfile = "";


while (<IF>) {

 
    ($alignment_id, $alig_type, $pdb_code, $pdb_chain, $mat11, $mat12, 
     $mat13, $mat21, $mat22, $mat23, $mat31, $mat32, $mat33, 
     $shift1, $shift2, $shift3) = split ();

    next if ($alignment_id =~ /alignment_id/);

    $pdb_code  =~ s/\"//g;
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

    (  -e  "$pdbdir/$pdb_code.pdb") || next;

    

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

    ###########################################################
    # tfm according to sysiphus
    $sys_chainfile_renamed = "pdbchains/$current_qry/$current_qry.to_$pdb_code$pdb_chain.sys.pdb";
    (-e $sys_chainfile_renamed) || next;

    ###########################################################
    ###########################################################
    ###########################################################
    # tfm according to MUSTANG
    $mustang_chainfile_renamed =  "pdbchains/$current_qry/$current_qry.to_$pdb_code$pdb_chain.mustang.pdb";
    $sys_mustang_score = "none";
    #if ( ! -e $mustang_chainfile_renamed ) {
	# run mustang
	$cmd = "time $mustang  -i  $chainfile  $qryfile > /dev/null \n";
	if (system $cmd) {
	    print LOG "Error running $cmd.\n";
	}
	# extract the matrix
	$mustang_chainfile_orig    = "$current_qry.to_$pdb_code$pdb_chain.mustang.pdb";
	$cmd = "$extr_chain results.pdb B > $mustang_chainfile_orig ";
	if (system $cmd) {
	    print LOG "Error running $cmd.\n";
	} else {
	    `mv $mustang_chainfile_orig $mustang_chainfile_renamed`;
	}
    `rm results.pdb results.html`;

    #}

}

close IF;


