#!/usr/bin/perl -w 

$top_path  = "/home/miles";
$tfm_table = "sysphus_tfms.csv";
#$tfm_table = "test.csv";


$nossa     = "$top_path/Projects/struct/09_tools/no_seq_sim_struct_superp_analyzer/nossa";


foreach ($top_path, $tfm_table, "params",
	  $nossa) {
    (-e $_) || die "$_ not found.\n";
}

foreach ("fold", "homologous", "fragment") {
    (-e $_) || `mkdir $_`;
}


($alignment_id, $pdb_code, $pdb_chain, $mat11, $mat12, 
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
	#$ctr++;
	#($ctr==11) && exit;
    }

    chdir $home;
    chdir "$home/$alig_type";

    (-e "pdbchains") || `mkdir pdbchains`;
    (-e "pdbchains/$current_qry") || `mkdir pdbchains/$current_qry`;
    (-e "sys_tfms")  || `mkdir sys_tfms`;


    $chainfile = "pdbchains/$current_qry/$pdb_code$pdb_chain.pdb";
    
    next if (! -e $chainfile || -z $chainfile);


    if ($is_query) {
	$qryfile = $chainfile;
	next;
    }

    $chainfile_chain_orig = "";
    $qry_chain_orig = "";

    @qry_name = split(/\//, $qryfile);
    $qry_pdb = substr ($qry_name[2], 0,4);
    $qry_chain = substr ($qry_name[2], 4, 1);
    if ($qry_chain eq ".") {
        $qry_chain = "A";
        $qry_chain_orig = "";
    } else {
        $qry_chain_orig = $qry_chain;
    }
    
    @chainfile_name = split(/\//, $chainfile);
    $chainfile_pdb = substr ($chainfile_name[2], 0, 4);
    $chainfile_chain = substr ($chainfile_name[2], 4, 1);
    if ($chainfile_chain eq ".") {
        $chainfile_chain = "A";
        $chainfile_chain_orig = "";
    } else {
        $chainfile_chain_orig = $chainfile_chain;
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
    # apply kaksi db to the same problem

#    $max_score{"kaksi"} = -1;
#    $max_match_no{"kaksi"} = -1;
#    for $match_no ( 0 .. 2 ) {
#	$struct_chainfile_renamed = 
#	    "pdbchains/$current_qry/$qry_pdb$qry_chain.to_$chainfile_pdb$chainfile_chain.kaksi.$match_no.pdb";

#	if (! -e $struct_chainfile_renamed) {
#		next;
#	}

# 
#	@aux = split " ", `$nossa $sys_chainfile_renamed $struct_chainfile_renamed`; 
#	chomp $aux[3];

#	if ( $max_score{"kaksi"} < $aux[3]) {
#	    $max_score{"kaksi"}    = $aux[3];
#	    $max_match_no{"kaksi"} = $match_no;
#	}
#    }
    
    ###########################################################
    ###########################################################
    ###########################################################
    # apply stride db to the same problem

#    $max_score{"stride"} = -1;
#    $max_match_no{"stride"} = -1;
#    for $match_no ( 0 .. 2 ) {
#	$struct_chainfile_renamed = 
#	    "pdbchains/$current_qry/$qry_pdb$qry_chain.to_$chainfile_pdb$chainfile_chain.stride.$match_no.pdb";

#	if (! -e $struct_chainfile_renamed) {
#		next;
#	}

# 
#	@aux = split " ", `$nossa $sys_chainfile_renamed $struct_chainfile_renamed`; 
#	chomp $aux[3];

#	if ( $max_score{"stride"} < $aux[3]) {
#	    $max_score{"stride"}    = $aux[3];
#	    $max_match_no{"stride"} = $match_no;
#	}
#    }
    
    ###########################################################
    ###########################################################
    ###########################################################
    # apply bii db to the same problem

    $max_score{"bii"} = -1;
    $max_match_no{"bii"} = -1;
    for $match_no ( 0 .. 2 ) {
	$struct_chainfile_renamed = 
	    "pdbchains/$current_qry/$qry_pdb$qry_chain.to_$chainfile_pdb$chainfile_chain.bii.$match_no.pdb";

	if (! -e $struct_chainfile_renamed) {
	    next;
	}

 
	@aux = split " ", `$nossa $sys_chainfile_renamed $struct_chainfile_renamed`; 
	chomp $aux[3];

	if ( $max_score{"bii"} < $aux[3]) {
	    $max_score{"bii"}    = $aux[3];
	    $max_match_no{"bii"} = $match_no;
	}
    }
    
    ###########################################################
    ###########################################################
    ###########################################################
    # apply comb db to the same problem

#    $max_score{"comb"} = -1;
#    $max_match_no{"comb"} = -1;
#    for $match_no ( 0 .. 2 ) {
#	$struct_chainfile_renamed = 
#	    "pdbchains/$current_qry/$qry_pdb$qry_chain.to_$chainfile_pdb$chainfile_chain.comb.$match_no.pdb";

#	if (! -e $struct_chainfile_renamed) {
#	    next;
#	}

# 
#	@aux = split " ", `$nossa $sys_chainfile_renamed $struct_chainfile_renamed`; 
#	chomp $aux[3];

#	if ( $max_score{"comb"} < $aux[3]) {
#	    $max_score{"comb"}    = $aux[3];
#	    $max_match_no{"comb"} = $match_no;
#	}
#    }
#    
    ###########################################################
    ###########################################################
    ###########################################################
    # apply comb db to the same problem
#    $alg = "bii_new";
#    $max_score{$alg} = -1;
#    $max_match_no{$alg} = -1;
#    for $match_no ( 0 .. 2 ) {
#	$struct_chainfile_renamed = 
#	    "pdbchains/$current_qry/$qry_pdb$qry_chain.to_$chainfile_pdb$chainfile_chain.bii_new.$match_no.pdb";

#	if (! -e $struct_chainfile_renamed) {
#	    next;
#	}

# 
#	@aux = split " ", `$nossa $sys_chainfile_renamed $struct_chainfile_renamed`; 
#	chomp $aux[3];

#	if ( $max_score{$alg} < $aux[3]) {
#	    $max_score{$alg}    = $aux[3];
#	    $max_match_no{$alg} = $match_no;
#	}
#    }
#    
    
    ###########################################################
    ###########################################################
    ###########################################################
    # apply comb db to the same problem
    $alg = "bii_python";
    $max_score{$alg} = -1;
    $max_match_no{$alg} = -1;
    for $match_no ( 0 .. 2 ) {
	$struct_chainfile_renamed = 
	    "pdbchains/$current_qry/$qry_pdb$qry_chain.to_$chainfile_pdb$chainfile_chain.bii_python.$match_no.pdb";

	if (! -e $struct_chainfile_renamed) {
	    next;
	}

 
	@aux = split " ", `$nossa $sys_chainfile_renamed $struct_chainfile_renamed`; 
	chomp $aux[3];

	if ( $max_score{$alg} < $aux[3]) {
	    $max_score{$alg}    = $aux[3];
	    $max_match_no{$alg} = $match_no;
	}
    }
    
    ###########################################################
    ###########################################################
    ###########################################################
    # apply struct to the same problem

#    $max_score{"struct"} = -1;
#    $max_match_no{"struct"} = -1;
#    for $match_no ( 0 .. 5 ) {
#	$struct_chainfile_renamed = 
#	    "pdbchains/$current_qry/$current_qry.to_$pdb_code$pdb_chain.struct.$match_no.pdb";

#	if (! -e $struct_chainfile_renamed) {
#	    next;
#	}

# 
#	@aux = split " ", `$nossa $sys_chainfile_renamed $struct_chainfile_renamed`; 
#	chomp $aux[3];

#	if ( $max_score{"struct"} < $aux[3]) {
#	    $max_score{"struct"}    = $aux[3];
#	    $max_match_no{"struct"} = $match_no;
#	}
#    }

#    print " $alig_type  $current_qry  $pdb_code$pdb_chain  ".
#    "   sys-kaksi (match_no: ".$max_match_no{"kaksi"}." ) ". $max_score{"kaksi"}.
#	"   sys-stride (match_no: ".$max_match_no{"stride"}." ) ". $max_score{"stride"}.
#	"   bii (match_no: ".$max_match_no{"bii"}." ) ". $max_score{"bii"}.
#	"   sys-comb (match_no: ".$max_match_no{"comb"}." ) ". $max_score{"comb"}.
#	"   fatcat (match_no: ".$max_match_no{"fatcat"}." ) ". $max_score{"fatcat"}.
#	"   struct (match_no: ".$max_match_no{"struct"}." ) ". $max_score{"struct"}.
#	"   sys-bii_new (match_no: ".$max_match_no{"bii_new"}." ) ". $max_score{"bii_new"}."\n"
    
    print " $alig_type  $current_qry  $pdb_code$pdb_chain  ".
#    "   kaksi " . $max_score{"kaksi"} .
#	"   stride " . $max_score{"stride"} .
	"   bii " . $max_score{"bii"}.
#	"   comb " . $max_score{"comb"}.
	"   fatcat " . $max_score{"fatcat"}.
#	"   struct " . $max_score{"struct"}.
#	"   bii_new " . $max_score{"bii_new"}.
	"   bii_python " . $max_score{"bii_python"}."\n"	
}

close IF;


