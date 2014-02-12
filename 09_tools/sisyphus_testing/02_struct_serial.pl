#!/usr/bin/perl -w 

#use strict;

my $top_path  = "/home/ivanam";
my $tfm_table = "sysphus_tfms.csv";
my $struct    = "$top_path/kode/03_struct/struct"; 

foreach ($top_path, $tfm_table,  $struct) {
    (-e $_) || die "$_ not found.\n";
}


my ($alignment_id, $alig_type, $pdb_code, $pdb_chain, $mat11, $mat12, 
 $mat13, $mat21, $mat22, $mat23, $mat31, $mat32, $mat33, 
 $shift1, $shift2, $shift3) = ();



open (IF, "<$tfm_table") ||
    die "Cno $tfm_table: $!\n";

(-e "sisyphus_alnmt_types") || die "sisyphus_alnmt_types\n";
chdir "sisyphus_alnmt_types";
my $home = `pwd`; chomp $home;
foreach ("fold", "homologous", "fragment") {
    (-e $_) || die "$_ not found in $home.\n";
}

my $qryfile = "";
my $is_query    = 0;
my $current_qry = "";
while (<IF>) {

 
    ($alignment_id, $alig_type, $pdb_code, $pdb_chain, $mat11, $mat12, 
     $mat13, $mat21, $mat22, $mat23, $mat31, $mat32, $mat33, 
     $shift1, $shift2, $shift3) = split ();

    next if ($alignment_id =~ /alignment_id/);

    $pdb_code =~ s/\"//g;
    $pdb_chain =~ s/\"//g;
    $alig_type =~ s/\"//g;

    $is_query = 0;
    if ( /1 0 0 0 1 0 0 0 1 0 0 0/ )  {
	$is_query    = 1;
	$current_qry = "$pdb_code$pdb_chain";
    } 



    chdir $home;
    chdir "$home/$alig_type";
    (-e "pdbchains") ||  die "pdbchains dir not found\n";
    (-e "sys_tfms")  || die "sys_tfms not found\n";


    my $chainfile     = "pdbchains/$current_qry/$pdb_code$pdb_chain.pdb";
    (-e $chainfile) || die "$chainfile not found in".`pwd`;
    (-z $chainfile)  && die "$chainfile empty in".`pwd`;
   
    if ($is_query) {
	$qryfile = $chainfile;
	print "\n\n\n##############################################\n";
	print "$pdb_code\n";
	(-e "pdbchains/$current_qry") || die "pdbchains/$current_qry not found\n"; 
	next;
    }

    print "$pdb_code\n";

    ###########################################################
    # tfm according to sysiphus
    my $sys_chainfile_renamed = "pdbchains/$current_qry/$current_qry.to_$pdb_code$pdb_chain.sys.pdb";
    (-e $sys_chainfile_renamed) || next;

    ###########################################################
    ###########################################################
    ###########################################################
    # apply struct to the same problem

    my $cmd = "time $struct  -from $qryfile -to $chainfile "; #-p ../params";
    print $cmd, "\n";
    if  (system $cmd ) {
	print  "Error running $cmd.\n"; 
	next;
    }
    my $name_root = "$current_qry\_to_$pdb_code$pdb_chain";


    for my  $match_no ( 0 .. 5 ) {
	my $struct_chainfile_orig    = "$name_root.$match_no.pdb";
	my $struct_chainfile_renamed = 
	    "pdbchains/$current_qry/$current_qry.to_$pdb_code$pdb_chain.struct.$match_no.pdb";

	if (-e $struct_chainfile_orig) {
	    `mv $struct_chainfile_orig $struct_chainfile_renamed`;
	    
	} else {
	    printf "$struct_chainfile_orig not found after\n$cmd\n";
	    #exit;
	}

    }
    # cleanup after ourselves
    `rm -f $name_root.struct_out`;
    `rm -f $name_root.*.aln`;
    `rm -f $name_root.*.pdb`;
    `rm -f *.db`;
}

close IF;


