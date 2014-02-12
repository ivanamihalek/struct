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

$runce     = "$top_path/Downloads/jfatcat/runCE.sh";
$runfatcat = "$top_path/Downloads/jfatcat/runFatcat.sh";

$runclick     = "$top_path/Downloads/Click/click";
$click_prms  = "$top_path/Downloads/Click/Parameters.inp";
$click_prms_local  = "Parameters.inp";

$struct      = "$top_path/kode/03_struct/struct"; 
$extr_matrix = "$top_path/perlscr/extractions/matrix_from_jce.pl";

foreach ($top_path, $tfm_table, $pdbdir, $pdbdown, "params",
	 $pdtfm, $extr, $struct,$pdb2seq, $afa2msf, $nossa,
         $runce, $runfatcat, $extr_matrix, $runclick, $click_prms) {
    (-e $_) || die "$_ not found.\n";
}

foreach ("fold", "homologous", "fragment") {
    (-e $_) || `mkdir $_`;
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

    ###########################################################
    ###########################################################
    ###########################################################
    # tfm according to CE
    $ce_chainfile_renamed =  "pdbchains/$current_qry/$current_qry.to_$pdb_code$pdb_chain.ce.pdb";
    $sys_ce_score = "none";
    if ( ! -e $ce_chainfile_renamed ) {
	# run ce
	$cmd = "$runce   -file1  $chainfile  -file2 $qryfile -printCE > tmp.out  \n";
	if (system $cmd) {
	    print LOG "Error running $cmd.\n";
	}
	# extract the matrix
	$cmd = "$extr_matrix tmp.out tmp";
	if (system $cmd) {
	    print LOG "Error running $cmd.\n";
	}


	# transform the qry file
	$ce_chainfile_orig    = "$current_qry.to_$pdb_code$pdb_chain.ce.pdb";
	$cmd = "$pdtfm $qryfile tmp.0.tfm > $ce_chainfile_orig \n";
	if (system $cmd) {
	    print LOG "Error running $cmd.\n";
	} else {
	    `mv $ce_chainfile_orig $ce_chainfile_renamed`;
	}
	
	`rm tmp.out tmp.0.tfm`;

    }

    if (  -e $ce_chainfile_renamed ) {
	@aux = split " ", `$nossa $sys_chainfile_renamed $ce_chainfile_renamed`; 
	chomp $aux[3];
	$sys_ce_score = $aux[3];
    } else {
	$sys_ce_score = "none";
    }
 
 


    ###########################################################
    ###########################################################
    ###########################################################
    # tfm according to Fatcat
    $fatcat_chainfile_renamed =  "pdbchains/$current_qry/$current_qry.to_$pdb_code$pdb_chain.fatcat.0.pdb";
    $max_score{"fatcat"}   = -1;
    if ( ! -e $fatcat_chainfile_renamed ) {
	# run fatcat
	$cmd = "$runfatcat   -file1  $chainfile  -file2 $qryfile -printCE -flexible true > tmp.out  \n";
	if (system $cmd) {
	    print LOG "Error running $cmd.\n";
	}
	# extract the matrix
	$cmd = "$extr_matrix tmp.out tmp";
	if (system $cmd) {
	    print LOG "Error running $cmd.\n";
	}


	# transform the qry file
	$fatcat_chainfile_orig    = "$current_qry.to_$pdb_code$pdb_chain.fatcat.pdb";
	$cmd = "$pdtfm $qryfile tmp.0.tfm > $fatcat_chainfile_orig \n";
	if (system $cmd) {
	    print LOG "Error running $cmd.\n";
	} else {
	    `mv $fatcat_chainfile_orig $fatcat_chainfile_renamed`;
	
	    @tfms_produced = split "\n", `ls tmp.*.tfm`;
	    $no_tfms = scalar @tfms_produced;
	
	    $oldfile = $fatcat_chainfile_renamed;
	    for $match_no ( 1 .. $no_tfms-1 ) {
		$newfile = "pdbchains/$current_qry/$current_qry.to_$pdb_code$pdb_chain.fatcat.$match_no.pdb";
		$cmd = "$pdtfm $oldfile tmp.$match_no.tfm > $newfile \n";
		if (system $cmd) {
		    print LOG "Error running $cmd.\n";
		    last;
		}

		$oldfile = $newfile;
	    }
	}
	`rm tmp.out tmp.*.tfm`;
    }

    $max_score{"fatcat"}    = -1;
    $max_match_no{"fatcat"} = -1;
    @fatcat_files = split "\n", `ls pdbchains/$current_qry/$current_qry.to_$pdb_code$pdb_chain.fatcat.*.pdb`;
    $match_no = 0;
    foreach $fatcat_match (@fatcat_files ) {
	@aux = split " ", `$nossa $sys_chainfile_renamed $fatcat_match`; 
	chomp $aux[3];

	if ( $max_score{"fatcat"} < $aux[3]) {
	    $max_score{"fatcat"}    = $aux[3];
	    $max_match_no{"fatcat"} = $match_no;
	}
	$match_no++;
    }
 		

    ###########################################################
    ###########################################################
    ###########################################################
    # tfm according to click
    chdir "pdbchains/$current_qry";
    $click_chainfile_renamed =  "$current_qry.to_$pdb_code$pdb_chain.click.0.pdb";
    $max_score{"click"}   = -1;
    if ( ! -e $click_chainfile_renamed ) {

	( -e $click_prms_local) || `cp  $click_prms  $click_prms_local`;
	# run click
	$cmd = "$runclick  $current_qry.pdb $pdb_code$pdb_chain.pdb  -a 3 \n";
	if (system $cmd) {
	    print LOG "Error running $cmd.\n";
	}

	# rename the files 
	for $match_no ( 1 .. 3 ) {
	    $click_chainfile_orig    = "$current_qry\-$pdb_code$pdb_chain.$match_no.pdb";
	    (-e $click_chainfile_orig) || next;
	    $i = $match_no-1;
	    $click_chainfile_renamed =  "$current_qry.to_$pdb_code$pdb_chain.click.$i.pdb";
	    `mv $click_chainfile_orig $click_chainfile_renamed`;
	}
	`rm *\-*.pdb *.clique $click_prms_local`;

	`rm $click_prms_local`;
    }
    $max_score   {"click"} = -1;
    $max_match_no{"click"} = -1;
    @click_files = split "\n", `ls $current_qry.to_$pdb_code$pdb_chain.click.*.pdb`;
    $match_no = 0;
    foreach $click_match (@click_files ) {
	@aux = split " ", `$nossa $current_qry.to_$pdb_code$pdb_chain.sys.pdb   $click_match`; 
	chomp $aux[3];

	if ( $max_score{"click"} < $aux[3]) {
	    $max_score{"click"}    = $aux[3];
	    $max_match_no{"click"} = $match_no;
	}
	$match_no++;
    }
    chdir "../../";



    ###########################################################
    ###########################################################
    ###########################################################
    # apply struct to the same problem

    $cmd = "time $struct  -in1 $qryfile -in2 $chainfile -p ../params";
    if  (system $cmd ) {
	print LOG "Error running $cmd.\n";
	print  "Error running $cmd.\n"; 
	next;
    }
    `rm *.struct_out*`;

    $max_score{"struct"} = -1;
    $max_match_no{"struct"} = -1;
    for $match_no ( 0 .. 2 ) {
	$struct_chainfile_orig    = "$current_qry.to_$pdb_code$pdb_chain.$match_no.pdb";
	$struct_chainfile_renamed = 
	    "pdbchains/$current_qry/$current_qry.to_$pdb_code$pdb_chain.struct.$match_no.pdb";

	if (-e $struct_chainfile_orig) {
	    `mv $struct_chainfile_orig $struct_chainfile_renamed`;

	} else {
	    printf "$struct_chainfile_orig not found after\n$cmd\n";
	    next;
	}

 
	@aux = split " ", `$nossa $sys_chainfile_renamed $struct_chainfile_renamed`; 
	chomp $aux[3];

	if ( $max_score{"struct"} < $aux[3]) {
	    $max_score{"struct"}    = $aux[3];
	    $max_match_no{"struct"} = $match_no;
	}
    }

    print "\n $alig_type  $current_qry  $pdb_code$pdb_chain  ".
	" sys-ce $sys_ce_score ".
	"   sys-fatcat (match_no: ".$max_match_no{"fatcat"}." ) ". $max_score{"fatcat"}.
	"   sys-click (match_no: ".$max_match_no{"click"}." ) ". $max_score{"click"}.
	"   sys-struct (match_no: ".$max_match_no{"struct"}." ) ". $max_score{"struct"}."\n";
    

}

close IF;


