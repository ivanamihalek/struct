#!/usr/bin/perl -w 

$top_path  = "/home/ivanam";
$tfm_table = "sysphus_tfms.csv";

$nossa     = "$top_path/kode/03_struct/09_tools/no_seq_sim_struct_superp_analyzer/nossa";


foreach ($top_path, $tfm_table,  $nossa) {
    (-e $_) || die "$_ not found.\n";
}


($alignment_id, $alig_type, $pdb_code, $pdb_chain, $mat11, $mat12, 
 $mat13, $mat21, $mat22, $mat23, $mat31, $mat32, $mat33, 
 $shift1, $shift2, $shift3) = ();



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
    }

    chdir $home;
    chdir "$home/$alig_type";

    (-e "pdbchains") ||  die "pdbchains dir not found\n";
    (-e "sys_tfms")  ||  die "sys_tfms not found\n";
    (-e "pdbchains/$current_qry") || die "pdbchains/$current_qry not found";
    

    my $chainfile     = "pdbchains/$current_qry/$pdb_code$pdb_chain.pdb";
    if (! -e $chainfile) {
	#print "$chainfile not found in ".`pwd`;
	next;
    }
    if (-z $chainfile) {
	#print  "$chainfile empty in ".`pwd`;
	next;
    }
   
    if ($is_query) {
	$qryfile = $chainfile;
	#print "\n\n\n##############################################\n";
	#print "$pdb_code\n";
	#(-e "pdbchains/$current_qry") || die "pdbchains/$current_qry not found\n"; 
	next;
    }

    ###########################################################
    ###########################################################
    ###########################################################
    # tfm according to sysiphus
    $sys_chainfile_renamed = "pdbchains/$current_qry/$current_qry.to_$pdb_code$pdb_chain.sys.pdb";
    (  -e $sys_chainfile_renamed ) || next;


    ###########################################################
    ###########################################################
    ###########################################################
    # tfm according to CE
    $ce_chainfile_renamed =  "pdbchains/$current_qry/$current_qry.to_$pdb_code$pdb_chain.ce.pdb";
    $sys_ce_score = "none";

    if (  -e $ce_chainfile_renamed ) {
	@aux = split " ", `$nossa $sys_chainfile_renamed $ce_chainfile_renamed`; 
	chomp $aux[3];
	$sys_ce_score = $aux[3];
    }
 
 
    ###########################################################
    ###########################################################
    ###########################################################
    # tfm according to mustang
    $mustang_chainfile_renamed =  "pdbchains/$current_qry/$current_qry.to_$pdb_code$pdb_chain.mustang.pdb";
    $sys_mustang_score = "none";

    if (  -e $mustang_chainfile_renamed ) {
	@aux = split " ", `$nossa $sys_chainfile_renamed $mustang_chainfile_renamed`; 
	chomp $aux[3];
	$sys_mustang_score = $aux[3];
    }
 
 


    ###########################################################
    ###########################################################
    ###########################################################
    # tfm according to Fatcat
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
    # tfm according to C implementation of ce
    $max_score{"ceC"}    = -1;
    $max_match_no{"ceC"} = -1;
    @ceC_files = split "\n", `ls pdbchains/$current_qry/$current_qry.to_$pdb_code$pdb_chain.ceC.*.pdb`;
    $match_no = 0;
    foreach $ceC_match (@ceC_files ) {
	@aux = split " ", `$nossa $sys_chainfile_renamed $ceC_match`; 
	chomp $aux[3];

	if ( $max_score{"ceC"} < $aux[3]) {
	    $max_score{"ceC"}    = $aux[3];
	    $max_match_no{"ceC"} = $match_no;
	}
	$match_no++;
    }
 		
 		
  		

    ###########################################################
    ###########################################################
    ###########################################################
    # struct

    $max_score{"struct"} = -1;
    $max_match_no{"struct"} = -1;
    for $match_no ( 0 .. 2 ) {
	$struct_chainfile_renamed = 
	    "pdbchains/$current_qry/$current_qry.to_$pdb_code$pdb_chain.struct.$match_no.pdb";

	if (! -e $struct_chainfile_renamed) {
	    next;
	}

 
	@aux = split " ", `$nossa $sys_chainfile_renamed $struct_chainfile_renamed`; 
	chomp $aux[3];

	if ( $max_score{"struct"} < $aux[3]) {
	    $max_score{"struct"}    = $aux[3];
	    $max_match_no{"struct"} = $match_no;
	}
    }

    print " $alig_type  $current_qry  $pdb_code$pdb_chain  ".
	" sys-ce $sys_ce_score ".
	" sys-mustang $sys_mustang_score ".
	"   sys-fatcat (match_no: ".$max_match_no{"fatcat"}." ) ". $max_score{"fatcat"}.
	"   sys-ceC    (match_no: ".$max_match_no{"ceC"}." ) ". $max_score{"ceC"}.
	"   sys-struct (match_no: ".$max_match_no{"struct"}." ) ". $max_score{"struct"}."\n";
    

}

close IF;


