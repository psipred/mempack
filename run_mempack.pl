#!/usr/bin/perl
##
## ********************************************************
## *   MEMPACK - Predicting transmembrane helix packing   *
## *       arrangements using residue contacts and        *
## *              a force-directed algorithm.             *
## * Copyright (C) 2009 Timothy Nugent and David T. Jones *
## ********************************************************
## 
## This program is copyright and may not be distributed without
## permission of the author unless specifically permitted under
## the terms of the license agreement.
## 
## THIS SOFTWARE MAY ONLY BE USED FOR NON-COMMERCIAL PURPOSES. PLEASE 
## CONTACT THE AUTHOR IF YOU REQUIRE A LICENSE FOR COMMERCIAL USE.

use strict;
use warnings;
use lib 'lib';
use Round qw(:all);
use Getopt::Long;

## NCBI / Database paths - these need to be set!
my $ncbidir = '';
my $dbname = '';

## Executable and input/output paths
my $mem_dir = '';
my $input_path = $mem_dir.'input/';
my $output_path = $mem_dir.'output/';
my $model_path = $mem_dir.'models/';
my $datadir = $mem_dir.'data/';
my $svm_classify = $mem_dir.'bin/svm_classify';
my $mtx = 0;
my $remove_files = 0;
my $globmem = 0;
my $graphics = 1;
my $format = 1;
my $cores = 1;
my $erase_previous = 0;
my $draw_rr_contacts = 1;

my (@mtx,$blast_out,$svm_all,%range,$header);
my ($system);
my %topology = ();

my ($sequence,@topology,%profile,@empty,$topology_string,$length,%aa,%aa3,@all_contact);
my $window = 7;
my $dp = 0.0001;
my $verbose = 0;
my $distance = 1;
my $def = 1;
my $kk_plot = $mem_dir.'bin/kk_plot';
my $font = 'lib/arial.ttf';

&main();
exit;

sub main{

	print "\n********************************************************\n";
	print "*   MEMPACK - Predicting transmembrane helix packing   *\n";
	print "*       arrangements using residue contacts and        *\n";
	print "*              a force-directed algorithm.             *\n";
	print "* Copyright (C) 2009 Timothy Nugent and David T. Jones *\n";
	print "********************************************************\n\n";

	&get_arguments();
	&get_normalisation_values();
	
	# If we've been passed .mtx files rather than fasta files
	if ($mtx){
		foreach (@ARGV){
			push @mtx,$_;
			$_ =~ s/\.mtx//g;
			$header = $_;
		}

		if ($erase_previous){
			my $erase_string = $input_path.$header."*";;
			if($header =~ /\//){
				my @split = split(/\//,$header);
				$erase_string = $input_path.$split[-1]."*";
			}
			$system = `rm $erase_string &> /dev/null`;
			$erase_string =~ s/^$input_path/$output_path/; 
			$system = `rm $erase_string &> /dev/null`;
		}	

	# Otherwise run PSI-BLAST to create .mtx files
	}else{
		my @fastas;	
		foreach (@ARGV){
			$header = $_;
			$header =~ s/\.fasta//g;
			$header =~ s/\.fa//g;			
			
			if ($erase_previous){
				my $erase_string = $input_path.$header."*";;
				if($header =~ /\//){
					my @split = split(/\//,$header);
					$erase_string = $input_path.$split[-1]."*";
				}
				$system = `rm $erase_string`;
				$erase_string =~ s/^$input_path/$output_path/; 
				$system = `rm $erase_string`;
			}	
			
			push @fastas,$_;
			&run_psiblast($_);
		}
		@ARGV = @fastas;
	}

	if($header =~ /\//){
		my @split = split(/\//,$header);
		$header = $split[-1];
	}

	$system = `rm $input_path/$header* &> /dev/null` if $erase_previous;
	$system = `rm $output_path/$header* &> /dev/null` if $erase_previous

	&load_mtx($header);

## Predict lipid exposure

	my $input_file = $input_path.$header."_LIPID_EXP.dat";

	if (-e $input_file){
		if (-s $input_file){
			`rm $input_file`;
		}else{
			exit;
		}
	}
	
	open (INPUT,">$input_file");

	my (%lipid_scores,@positions,@residues);
	
	## Loop through all the helices
	my $helix1_count = 1;
	for (my $h = 0; $h < scalar @topology; $h+=2){	
					
		#print "TMH $helix1_count :\t$topology[$h] - $topology[$h+1]\n";
		my $helix1_length = $topology[$h+1] - $topology[$h] + 1;

		## Loop through all positions in this helix
		for my $p1 ($topology[$h]..$topology[$h+1]){
						
			my ($array,$p1_seq) = &create_lipid_input($p1);								
				
			print INPUT "0 ";
			my $feature = 1;
			foreach my $v (@{$array}){
				print INPUT $feature.":".$v." ";
				$feature++;
			}
			
			push @positions,$p1;
			push @residues,substr($p1_seq,3,1);
			my $comment = $p1."_".$p1_seq;
			print INPUT "# $comment\n";
	
		}
		$helix1_count++;
	}

	close INPUT;

	my $model = $model_path."LIPID_EXPOSURE_ALL.model";
	my $prediction = $output_path.$header."_LIPID_EXPOSURE.predictions";
		
	if (-e $model){

		if (-e $prediction){
			print "$prediction exists!\n\n";
		}else{
			print "$svm_classify -v 0 $input_file $model $prediction\n";
			$system = `$svm_classify -v 0 $input_file $model $prediction`;	
		}
	}else{
		die "$model doesn't exist.\n";
	}


	my $lipid_out = $output_path.$header."_LIPID_EXPOSURE.results";
	
	
	if (-e $lipid_out){

		print "$lipid_out exists!\n\n";

		unless(open(LIPID,$prediction)){
			die "Couldn't open $prediction\n";	
		}else{

			my $count = 0;
			while(<LIPID>){
				my $value = $_;
				$value =~ s/\s+//g;
				$lipid_scores{$positions[$count]} = $value;
				$count++;
			}			
		}

	}else{

		unless(open(LIPID_OUT,">$lipid_out")){
			die "Couldn't write to $lipid_out\n";	
		}
		print LIPID_OUT "Pos\tRes\tScore\n";

		unless(open(LIPID,$prediction)){
			die "Couldn't open $prediction\n";	
		}else{

			if (scalar @positions != scalar @residues && scalar @residues != scalar keys %lipid_scores){
				die "Predictions do not equal number of residues!\n\n";
			}

			my $count = 0;
			while(<LIPID>){
				my $value = $_;
				$value =~ s/\s+//g;
				$lipid_scores{$positions[$count]} = $value;
				print  LIPID_OUT "$positions[$count]\t$residues[$count]\t$value\n";				
				$count++;
			}			
		}
		close LIPID_OUT;
		print "Written $lipid_out\n\n";
	
	}

## Generate SVM input files for residue-residue contact prediction

	$input_file = $input_path.$header."_CONTACT.dat";

	$helix1_count = 1;	
	for (my $h = 0; $h < scalar @topology; $h+=2){	
		## Loop through all positions in this helix
		my $helix1_pos = 1;
		for my $p1 ($topology[$h]..$topology[$h+1]){			
			## Loop through all the helices after current one
			my $helix2_count = $helix1_count + 1;				
			for (my $nh = $h+2; $nh < scalar @topology; $nh+=2){	
				my $helix2_pos = 1;
				for my $p2 ($topology[$nh]..$topology[$nh+1]){				
					my $tmp_res = $p1."-".$p2;
					push @all_contact,$tmp_res;
					$helix2_pos++;
				}				
				$helix2_count++;				
			}			
			$helix1_pos++;		
		}		
		$helix1_count++;
	}

	if (!-e $input_file){

	open (INPUT,">$input_file");

	## Loop through all the helices
	$helix1_count = 1;
	for (my $h = 0; $h < scalar @topology; $h+=2){	
					
		#print "TMH $helix1_count :\t$topology[$h] - $topology[$h+1]\n";
		my $helix1_length = $topology[$h+1] - $topology[$h] + 1;

		## Loop through all positions in this helix
		my $helix1_pos = 1;
		for my $p1 ($topology[$h]..$topology[$h+1]){
			
			## Loop through all the helices after current one
			my $helix2_count = $helix1_count + 1;			
			
			for (my $nh = $h+2; $nh < scalar @topology; $nh+=2){	
			
				my $helix2_length = $topology[$nh+1] - $topology[$nh] + 1;
			
				## Loop through all positions in next helix
				my $helix2_pos = 1;
				for my $p2 ($topology[$nh]..$topology[$nh+1]){
								
					print "$p1 (TMH ".$helix1_count.") --> $p2 (TMH ".$helix2_count.")\n" if $verbose;					
					my $distance = $p2 - $p1;
					
					# Get windows around these two residues
					my ($array,$p1_p2_seq) = &create_contact_input($p1,$p2);
					
					## Global features
					print "Helix $helix1_count length = $helix1_length\n" if $verbose;
					print "Residue 1 pos in helix = $helix1_pos\n" if $verbose;			
					my ($relative_pos1,$relative_pos2);
					
					if ($helix1_count % 2){
						$relative_pos1 = $helix1_pos/$helix1_length;
					}else{
						$relative_pos1 = 1 - ($helix1_pos/$helix1_length);
					}
					print "Residue 1 relative pos in helix = $relative_pos1\n" if $verbose;				
				
					print "Helix $helix2_count length = $helix2_length\n" if $verbose;
					print "Residue 2 pos in helix = $helix2_pos\n" if $verbose;
					
					if ($helix2_count % 2){
						$relative_pos2 = $helix2_pos/$helix2_length;											
					}else{
						$relative_pos2 = 1 - ($helix2_pos/$helix2_length);	
					}				
					print "Residue 2 relative pos in helix = $relative_pos2\n" if $verbose;					
					print "Distance = $distance\n\n" if $verbose;
					
					my $s1 = substr($sequence,$p1-1,1);
					my $s2 = substr($sequence,$p2-1,1);
					my $res1 = $aa3{$s1};
					my $res2 = $aa3{$s2};			
					my $label = 0;	
					
					print INPUT "$label ";

					my $feature = 1;
					foreach my $v (@{$array}){
						print INPUT $feature.":".$v." ";
						$feature++;
					}

					## Relative position in helix 1
					print  INPUT  $feature.":",nearest_ceil($dp,$relative_pos1)," ";
					$feature++;
					## Relative position in helix 2
					print  INPUT  $feature.":",nearest_ceil($dp,$relative_pos2)," ";
					$feature++;

					## Distance between residues
					if ($distance <= 25){
						print INPUT $feature.":1 ";
					}else{
						print INPUT $feature.":0 ";
					}
					$feature++;
					if ($distance > 25 && $distance <= 50){
						print INPUT $feature.":1 ";
					}else{
						print INPUT $feature.":0 ";
					}
					$feature++;					
					if ($distance > 50 && $distance <= 75){
						print INPUT $feature.":1 ";
					}else{
						print INPUT $feature.":0 ";
					}
					$feature++;
					if ($distance > 75 && $distance <= 100){
						print INPUT $feature.":1 ";
					}else{
						print INPUT $feature.":0 ";
					}
					$feature++;	
					if ($distance > 100 && $distance <= 125){
						print INPUT $feature.":1 ";
					}else{
						print INPUT $feature.":0 ";
					}
					$feature++;
					if ($distance > 125 && $distance <= 150){
						print INPUT $feature.":1 ";
					}else{
						print INPUT $feature.":0 ";
					}
					$feature++;					
					if ($distance > 150 && $distance <= 175){
						print INPUT $feature.":1 ";
					}else{
						print INPUT $feature.":0 ";
					}
					$feature++;
					if ($distance > 175 && $distance <= 200){
						print INPUT $feature.":1 ";
					}else{
						print INPUT $feature.":0 ";
					}
					$feature++;						
					if ($distance > 200){
						print INPUT $feature.":1 ";
					}else{
						print INPUT $feature.":0 ";
					}
					$feature++;							
						
					## Add lipid exposure scores					
					print INPUT $feature.":".$lipid_scores{$p1}." ";
					$feature++;
					print INPUT $feature.":".$lipid_scores{$p2}." ";
					$feature++;					
									
					my $comment = $p1."_".$p2."_".$p1_p2_seq."_".$res1."_".$res2;

					print INPUT "# $comment\n";

					$helix2_pos++;
				}				
				$helix2_count++;				
			}			
			$helix1_pos++;		
		}		
		$helix1_count++;
	}

	close INPUT;
	
	}else{
		print "$input_file exits!\n\n";
	}

	$model = $model_path."CONTACT_ALL_DEF1.model";
	$prediction = $output_path.$header."_CONTACT_DEF1.predictions";
	my $output = $output_path.$header."_CONTACT_DEF1.results";
	my $graph_out = $output_path.$header."_graph.out";
		
	if($def == 3){
		$model = $model_path."CONTACT_ALL_DEF3.model";
		$prediction = $output_path.$header."_CONTACT_DEF3.predictions";	
		$output = $output_path.$header."_CONTACT_DEF3.results";
	}elsif($def == 2){
		$model = $model_path."CONTACT_ALL_DEF2.model";
		$prediction = $output_path.$header."_CONTACT_DEF2.predictions";	
		$output = $output_path.$header."_CONTACT_DEF2.results";
	}		
		
	if (-e $model){
		if (-e $prediction){
			print "$prediction exists!\n\n";
		}else{
			print "$svm_classify -v 0 $input_file $model $prediction\n";
			$system = `$svm_classify -v 0 $input_file $model $prediction`;	
		}
	}else{
		die "$model doesn't exist.\n";
	}

## Construct results file

	my %residue_to_helix;

	my $helix = 1;
	for (my $h = 0;$h < scalar @topology; $h += 2){	
		my $r_start = $topology[$h];
		my $r_stop = $topology[$h+1]+1;
		for my $r ($r_start..$r_stop){		
			$residue_to_helix{$r} = $helix;		
		}
		$helix++;	
	}

	my %contact_scores;
	my @predicted_contact = ();
	my %pred_hh;

	if (-e $output){
	
		print "$output exists!\n\n";
		unless(open(OUTPUT,$output)){
			die "Couldn't open $output\n";	
		}else{
			while(<OUTPUT>){
				next if $_ =~ /Topology/;				
				my @split = split(/\s+/,$_);
				push @predicted_contact,$split[0];
				$pred_hh{$split[1]} = 1;
			}
			close OUTPUT;			
		}	
	}else{

		open(OUT,">$output");
		my $top_string = "";
		print OUT "# Topology:\t";
		foreach my $t (@topology){
			$top_string = $top_string.$t.",";
		
		}
		chop $top_string;
		print OUT "$top_string\n";

		unless(open(CONTACT,$prediction)){
			die "Couldn't open $prediction\n";	
		}else{
			my $count = 0;
			while(<CONTACT>){
				my $value = $_;
				$value =~ s/\s+//g;
				if($value > 0){

					my @split = split(/-/,$all_contact[$count]);		
					if (exists $residue_to_helix{$split[0]} && exists $residue_to_helix{$split[1]}){			
						print OUT "$all_contact[$count]\t",$residue_to_helix{$split[0]},"-",$residue_to_helix{$split[1]},"\t$value\n";
						push @predicted_contact,$all_contact[$count];
						my $hh = $residue_to_helix{$split[0]}."-".$residue_to_helix{$split[1]};
						$pred_hh{$hh} = 1;
					}			
				}							
				$count++;
			}
			close CONTACT;			
		}	
		close OUT;	
		print "Written $output\n\n";
	}

	exit unless $graphics;

## Draw layout

	unless (scalar keys %pred_hh){
		print "\nNo predicted contacts!\n\n";
		exit;
	}

	print "Generating layout...\n";
	my $system = `rm  $graph_out` if -e $graph_out;
	
	print "$kk_plot $output > $graph_out\n\n";
	$system = `$kk_plot $output > $graph_out`;
		
	if (!-e $graph_out){
		die "Couldn't plot layout!\n\n";
	}

	my @rotations = ();
	my $image_no = 1;
	my %kk_graph;

	if (-e $graph_out && $graphics){

		if(load_module('DrawHelicalWheel')){
			die "Couldn't load module DrawHelicalWheel (requires GD) - graphical output disabled.\n\n";
		}
		if(load_module('Image::Magick')){
			die "Couldn't load module Image::Magick - graphical output disabled.\n\n";
		}

		open(KK,$graph_out);
		my @kk = <KK>;
		close KK;	
		my $tag = 0;
		
		foreach my $k (@kk){
			
			if ($k =~ /===/){
				
				my @selected_helices = keys %kk_graph;
				my $im;

				if($draw_rr_contacts){

					$im = DrawHelicalWheel->new(-title=>$header,
			       	    		-sequence=>$sequence,
                                    		-helices=>\@topology,
                                    		-ttf_font=>$font,
				    		-hh_contacts=>\%pred_hh,
				    		-graph_centres=>\%kk_graph,
				    		-residue_contacts=>\@predicted_contact,
				    		-rotations=>\@rotations,
				    		-helix_diameter=>235,
				    		-only_these_helices=>,\@selected_helices				    
				    );

				}else{

					$im = DrawHelicalWheel->new(-title=>$header,
			       	    		-sequence=>$sequence,
                                    		-helices=>\@topology,
                                    		-ttf_font=>$font,
				    		-hh_contacts=>\%pred_hh,
				    		-graph_centres=>\%kk_graph,
				    		-rotations=>\@rotations,
				    		-helix_diameter=>235,
				    		-only_these_helices=>,\@selected_helices				    
				    	);				
				
				}

				my $svg = $output_path.$header.'_Kamada-Kawai.svg';
				open(OUTPUT, ">$svg");
				binmode OUTPUT;
				print OUTPUT $im->svg; 
				close OUTPUT;

				my $jpg = $output_path.$header.'_Kamada-Kawai_'.$image_no.'.jpg';
				print "Generating JPG image $jpg\n\n";
				my $p = Image::Magick->new(magick=>'jpg',quality=>70);
				$p->Read($svg);
				$p->Write($jpg);
				$system = `rm $svg`;
				$image_no++;
				@rotations = ();
				%kk_graph = ();	
			}
			
			
			next unless $k =~ /\(/;
			my @split = split(/\s+/,$k);
			my $x = 0;
			my $y = 0;
			if ($split[1] =~ /\((-?\d+\.?\d*e?-?\d*),(-?\d+\.?\d*e?-?\d*)\)/){
				$x = $1;
				$y = $2;
			}			
			$kk_graph{$split[0]} = [$x,$y];
			push @rotations,$split[2];
		}
	}else{
		if ($graphics){
			die "Couldn't create $graph_out!\n";
		}
	}

	# Clean up
	$system = `rm $input_path/$header* &> /dev/null` if $remove_files;
}

# Create input files for SVM classify
sub create_contact_input {

	my ($pos1,$pos2) = @_;
	my $p1_p2_seq = "";
	my $array = ();
	
	## Get window for 1st residude
		
	print "$pos1\t- central position in SW sequence\n" if $verbose;
	$pos1--; ## Since position 1 in @profile is 0

	## Push all the arrays onto @array
	
	for my $w ($pos1 - (($window- 1)/2) .. $pos1 + (($window- 1)/2)){
		print "$w\t- window positions required from profile\n" if $verbose;
		
		if (defined @{$profile{$w}} ){
			push @{$array},@{$profile{$w}};
		}else{
			push @{$array},@empty;
		}
	}
		
	## Generate sequence window for the comment		
	my $window_seq1 = "";		
	for my $j ($pos1 - (($window- 1)/2) .. $pos1 + (($window- 1)/2)){
			
		print "j = $j\t- substring position required from sequence\n" if $verbose;
		if ($j < 0 || $j >= $length){
			$window_seq1 .= 'X';
		}else{
			$window_seq1 .= substr($sequence,$j,1);
		}
	}
	print "$window_seq1\n\n" if $verbose;

	## Get window for 2md residude
		
	print "$pos2\t- central position in SW sequence\n" if $verbose;
	$pos2--; ## Since position 1 in @profile is 0

	## Push all the arrays onto @array
	
	for my $w ($pos2 - (($window- 1)/2) .. $pos2 + (($window- 1)/2)){
		print "$w\t- window positions required from profile\n" if $verbose;
		push @{$array},@{$profile{$w}}
	}
		
	## Generate sequence window for the comment		
	my $window_seq2 = "";		
	for my $j ($pos2 - (($window- 1)/2) .. $pos2 + (($window- 1)/2)){
			
		print "j = $j\t- substring position required from sequence\n" if $verbose;
		if ($j < 0 || $j >= $length){
			$window_seq2 .= 'X';
		}else{
			$window_seq2 .= substr($sequence,$j,1);
		}
	}
	print "$window_seq2\n\n" if $verbose;
	
	$p1_p2_seq = $window_seq1."-".$window_seq2;

	return($array,$p1_p2_seq);

}

# Create input files for SVM classify
sub create_lipid_input {

	my $pos1 = shift;
	my $array = ();

	$pos1--; ## Since position 1 in @profile is 0
	
	for my $w ($pos1 - (($window- 1)/2) .. $pos1 + (($window- 1)/2)){
		#print "$w\t- window positions required from profile\n";
		
		if (defined @{$profile{$w}} ){
			push @{$array},@{$profile{$w}};
		}else{
			push @{$array},@empty;
		}
	}
		
	## Generate sequence window for the comment		
	my $window_seq1 = "";		
	for my $j ($pos1 - (($window- 1)/2) .. $pos1 + (($window- 1)/2)){
			
		#print "j = $j\t- substring position required from sequence\n" if $verbose;
		if ($j < 0 || $j >= $length){
			$window_seq1 .= 'X';
		}else{
			$window_seq1 .= substr($sequence,$j,1);
		}
	}
	#print "$window_seq1\n\n" if $verbose;

	return($array,$window_seq1);

}

# Load PSI-BLAST mtx file
sub load_mtx{

	my $header = shift;

	my $mtx = $output_path.$header.".mtx";

	if (-e $mtx[0]){
		$mtx = $mtx[0];
	}

	if (-e $mtx){
		open (MTX,$mtx);
		@mtx = <MTX>;
		close MTX;	
	}else{
		die "Couldn't find $mtx\n";	
	}

	$length = $mtx[0];
	$sequence = $mtx[1];
	$length =~ s/\s+//g;
	$sequence =~ s/\s+//g;
	#print "$sequence\n";

	for (1..($length + $window- 1)){

		if (($_ <= ($window- 1)/2)||($_ > $length + (($window- 1)/2))){

			## Missing
			@{$profile{$_}} = @empty;

		}else{
			## Start on MTX line 15

			if (defined $mtx[$_+13-(($window- 1)/2)]){

				@{$profile{$_}} = split(/\s+/,$mtx[$_+13-(($window- 1)/2)]);
			
				## Now do the Z score normalisation
			
				## Not normalised
				#print "Not normalised: @{$profile{$_}}\n";
			
				my $z_count = 0;
				foreach my $z (@{$profile{$_}}){
			
					## Normalise to Z score
					$z= ($z - $range{$z_count}{'mean'})/$range{$z_count}{'sd'};
					## Scale
					$z = ($z + (-1 * $range{$z_count}{'lower'}))/$range{$z_count}{'range'};
					## 5 Decimal places
					$z = nearest_ceil($dp,$z);
					$z_count++;
				}

    				## Get 20 residues rather than 28
				my @aa_20 = (${$profile{$_}}[1],${$profile{$_}}[3],${$profile{$_}}[4],${$profile{$_}}[5],${$profile{$_}}[6],
				${$profile{$_}}[7],${$profile{$_}}[8],${$profile{$_}}[9],${$profile{$_}}[10],${$profile{$_}}[11],
				${$profile{$_}}[12],${$profile{$_}}[13],${$profile{$_}}[14],${$profile{$_}}[15],${$profile{$_}}[16],
				${$profile{$_}}[17],${$profile{$_}}[18],${$profile{$_}}[19],${$profile{$_}}[21],${$profile{$_}}[22]);
			
				@{$profile{$_}} = @aa_20;
		

			}
		}
	}

}

# Get arguments
sub get_arguments{

	# Process command line arguments
	if (!$ARGV[0]){
		&usage;
	}else{

		my $result = GetOptions ("n=s" => \$ncbidir,
                                 	"mtx=i" => \$mtx,
				 	"d=s" => \$dbname,
					"j=s" => \$output_path,
					"t=s" => \$topology_string,
					"w=s" => \$mem_dir,
				 	"e=i" => \$remove_files,
					"a=i" => \$def,
				 	"g=i" => \$graphics,
					"f=i" => \$erase_previous,
					"c=i" => \$cores,
					"r=i" => \$draw_rr_contacts,							
			         	"h"  => sub {&usage;});

		## Get rid of trailing slashes
		$ncbidir =~ s/\/+$//;
		
		unless($mtx){
		
			## Check the NCBI directory
			unless (-d $ncbidir){

				## Look for the NCBI directory
				my $system = `which blastpgp`;
				if ($system =~ /(.*)\/blastpgp$/){
					$ncbidir = $1;
					$ncbidir =~ s/\s+//g;	
				}

				unless (-d $ncbidir){
					print "NCBI directory $ncbidir doesn't exist. Please pass it using\n"; 
					print "the -n paramater or modifiy the value at the top of the script.\n\n";
					exit;
				}
			}

			## Make sure we can find blastpgp & makemat
			my $psiblast = $ncbidir."/blastpgp";
			my $makemat = $ncbidir."/makemat";
			unless (-e $psiblast){
				print "Can't find the program blastpgp in the NCBI directory $ncbidir\n"; 
				print "Please pass the correct NCBI location using the -n paramater or modifiy\n";
				print "the value at the top of the script..\n\n";
				exit;			
			}
			unless (-e $makemat){
				print "Can't find the program makemat in the NCBI directory $ncbidir\n"; 
				print "Please pass the correct NCBI location using the -n paramater or modifiy\n";
				print "the value at the top of the script..\n\n";
				exit;				
			}

			if ($dbname eq ''){
				print "The database name for PSI-BLAST searches has not been set. Please pass it using\n"; 
				print "the -d paramater or modifiy the value at the top of the script.\n\n";
				exit;		
			}
		}
		
		if (defined $topology_string){		
			@topology = split(/,/,$topology_string);			
			#shift @topology;						
			if(scalar @topology % 2){
				print "Uneven number of helix boundaries detected.\n\n"; 
				exit;	
			}		
		}else{
			print "Topology for this sequence was not provided. Please pass it using\n"; 
			print "the -t paramater.\n\n";
			exit;			
		}
				
		if (defined $mem_dir){
			$svm_classify = $mem_dir.'bin/svm_classify';
			$input_path = $mem_dir.'input/';
			my $new_lib = $mem_dir."lib";
			$font = $mem_dir.$font;
			unshift @INC, $new_lib;
		}
		
		if($output_path =~ /^output\/$/)
		{
			$output_path = $mem_dir.'output/';
		}
		$model_path = $mem_dir.'models/';
		$datadir = $mem_dir.'data/';

		unless (-e $svm_classify){
			print "Can't find $svm_classify - have you run make?\n\n"; 
			exit;			
		}
				
	}
}

# Usage
sub usage {

	print "Version 1.0\n\n";
	print "Usage: run_mempack.pl [options] <fasta file>\n\n";
	print "Options:\n\n";
	print "-a <1|2|3>     Distance between atoms in interacting residue pair. Default 1.\n";
	print "               1 = Less than 8 angstroms between C-beta atoms (C-alpha for glycine).\n";
	print "               2 = Less than the sum of their van der Waals radii plus a threshold of 0.6 angstroms.\n";
	print "               3 = Less than 5.5 angstroms between sidechain or backbone heavy atoms.\n";
	print "-mtx <0|1>     Process PSI-BLAST .mtx files instead of fasta files. Default 0.\n";
	print "-n <directory> NCBI binary directory (location of blastpgp and makemat)\n";
	print "-d <path>      Database for running PSI-BLAST.\n";
	print "-t <topology>  Transmembrane topology helix boundaries in the form: 20,35,50,70,91,108\n";
	print "-j <path>      Output path for all files. Default: output/\n";
	print "-w <path>      Directory that contains mempack. Default ''\n";
	print "-e <0|1>       Erase intermediate files. Default 0.\n";
	print "-f <0|1>       Erase files from previous runs. Default 0.\n";
	print "-g <0|1>       Draw schematic. Default 1.\n";
	print "-r <0|1>       Draw residue-residue contacts. Default 1.\n";
	print "-c <int>       Number of CPU cores to use for PSI-BLAST. Default 1.\n";
	print "-h <0|1>       Show help. Default 0.\n\n";
	exit;
}

# Run PSI-BLAST
sub run_psiblast {

	my $fasta = shift;
	my $mtx;
	if ($fasta =~ /\//){
		my @tmp = split(/\//,$fasta);
		$mtx = $tmp[-1].".mtx"; 
	}else{
		$mtx = $fasta.".mtx"; 
	}	
	
	if ($fasta =~ /mtx/){
		print "This looks like an .mtx file! It should be a fasta file, otherwise\n";
		print "pass the -mtx 1 flag.\n\n";
		exit;
	}	
	
	$mtx =~ s/\.fa//;
	$mtx =~ s/\.fasta//;
	my $out_mtx = $output_path.$mtx;

	my $blast_out = "mempack_tmp.out"; 

	die "Fasta file $fasta doesn't exist!\n" unless -e $fasta;

	unless (-e $out_mtx || -e $mtx){
	
		print "Running PSI-BLAST: $fasta\n";

		my $system = `cp -f $fasta mempack_tmp.fasta`;
		print "$ncbidir/blastpgp -a $cores -j 2 -h 1e-3 -e 1e-3 -b 0 -d $dbname -i mempack_tmp.fasta -C mempack_tmp.chk >& $blast_out\n\n";
		$system = `$ncbidir/blastpgp -a $cores -j 2 -h 1e-3 -e 1e-3 -b 0 -d $dbname -i mempack_tmp.fasta -C mempack_tmp.chk >& $blast_out`;
		
		unless (-e "mempack_tmp.chk"){
		
			print "There was an error running PSI-BLAST. Did you set the database path correctly?\n\n";
			open(ERROR,$blast_out);
			my @error = <ERROR>;
			close ERROR;
			foreach my $line (@error){
				print $line;
			}
			print "\n";
			exit;				
		}		
		
		$system = `echo mempack_tmp.chk > mempack_tmp.pn`;
		$system = `echo mempack_tmp.fasta > mempack_tmp.sn`;
		$system = `echo "$ncbidir/makemat -P mempack_tmp"`;		
		$system = `$ncbidir/makemat -P mempack_tmp`;
		$system = `cp mempack_tmp.mtx $mtx`;
		$system = `rm -f mempack_tmp*`;
		$system = `rm -f error.log` if -e "error.log";	

	}
	
	if (-e $mtx){
		$system = `mv $mtx $output_path`;
		$mtx = $output_path.$mtx;	
		push @mtx,$mtx;;
	}elsif (-e $out_mtx){	
		$mtx = $out_mtx;	
		push @mtx,$mtx;;		
	}else{
		die "Problem creating $mtx\n";
	}	
}

# Load graphics modules
sub load_module {
	eval "require $_[0]";
    	return(1) if $@;
    	$_[0]->import(@_[1 .. $#_]);
}

# Normalisation values
sub get_normalisation_values {

        $range{0}{'mean'} = -32768;
        $range{0}{'sd'} = 1;
        $range{0}{'lower'} = 0;
        $range{0}{'upper'} = 0;
        $range{0}{'range'} = 1;
        $range{1}{'mean'} = -65.7855796773707;
        $range{1}{'sd'} = 165.265945369071;
        $range{1}{'lower'} = -2.67577460882869;
        $range{1}{'upper'} = 4.54289344365968;
        $range{1}{'range'} = 7.21866805248836;
        $range{2}{'mean'} = -32768;
        $range{2}{'sd'} = 1;
        $range{2}{'lower'} = 0;
        $range{2}{'upper'} = 33266;
        $range{2}{'range'} = 33266;
        $range{3}{'mean'} = -263.751620684456;
        $range{3}{'sd'} = 123.206910338565;
        $range{3}{'lower'} = -3.37844994385232;
        $range{3}{'upper'} = 11.223003781888;
        $range{3}{'range'} = 14.6014537257403;
        $range{4}{'mean'} = -242.89740690487;
        $range{4}{'sd'} = 188.205050583637;
        $range{4}{'lower'} = -2.4234343960522;
        $range{4}{'upper'} = 5.73787687182693;
        $range{4}{'range'} = 8.16131126787913;
        $range{5}{'mean'} = -194.765716870194;
        $range{5}{'sd'} = 180.643157809653;
        $range{5}{'lower'} = -2.38168048182125;
        $range{5}{'upper'} = 5.31858349089852;
        $range{5}{'range'} = 7.70026397271977;
        $range{6}{'mean'} = -135.345507311925;
        $range{6}{'sd'} = 251.687254013587;
        $range{6}{'lower'} = -2.19182530656994;
        $range{6}{'upper'} = 4.06196774373301;
        $range{6}{'range'} = 6.25379305030295;
        $range{7}{'mean'} = -192.372795115332;
        $range{7}{'sd'} = 230.165335412046;
        $range{7}{'lower'} = -2.14466354805761;
        $range{7}{'upper'} = 4.12908743801041;
        $range{7}{'range'} = 6.27375098606802;
        $range{8}{'mean'} = -221.253316749585;
        $range{8}{'sd'} = 155.427737434726;
        $range{8}{'lower'} = -2.65555357147081;
        $range{8}{'upper'} = 8.40424841993316;
        $range{8}{'range'} = 11.059801991404;
        $range{9}{'mean'} = -92.5749283883612;
        $range{9}{'sd'} = 246.71482166236;
        $range{9}{'lower'} = -2.51474584068856;
        $range{9}{'upper'} = 3.33411232793338;
        $range{9}{'range'} = 5.84885816862194;
        $range{10}{'mean'} = -179.119704507764;
        $range{10}{'sd'} = 167.0612436708;
        $range{10}{'lower'} = -2.51931738471669;
        $range{10}{'upper'} = 5.86084290403797;
        $range{10}{'range'} = 8.38016028875466;
        $range{11}{'mean'} = -89.1816674204734;
        $range{11}{'sd'} = 236.538032097989;
        $range{11}{'lower'} = -2.59078139419735;
        $range{11}{'upper'} = 3.07426954122637;
        $range{11}{'range'} = 5.66505093542372;
        $range{12}{'mean'} = -88.3949193426805;
        $range{12}{'sd'} = 185.952348451384;
        $range{12}{'lower'} = -2.88033512412105;
        $range{12}{'upper'} = 5.65410938930711;
        $range{12}{'range'} = 8.53444451342817;
        $range{13}{'mean'} = -188.69587667722;
        $range{13}{'sd'} = 183.448068516649;
        $range{13}{'lower'} = -2.52008171610062;
        $range{13}{'upper'} = 5.73838626478471;
        $range{13}{'range'} = 8.25846798088533;
        $range{14}{'mean'} = -249.184494195688;
        $range{14}{'sd'} = 178.122836782516;
        $range{14}{'lower'} = -2.32320298328491;
        $range{14}{'upper'} = 6.35058656502781;
        $range{14}{'range'} = 8.67378954831272;
        $range{15}{'mean'} = -162.548130559325;
        $range{15}{'sd'} = 164.437104518497;
        $range{15}{'lower'} = -2.5690787409429;
        $range{15}{'upper'} = 6.32185864378553;
        $range{15}{'range'} = 8.89093738472843;
        $range{16}{'mean'} = -207.047225991256;
        $range{16}{'sd'} = 173.171717333037;
        $range{16}{'lower'} = -2.25182714599291;
        $range{16}{'upper'} = 5.99432310297486;
        $range{16}{'range'} = 8.24615024896777;
        $range{17}{'mean'} = -90.7070707070707;
        $range{17}{'sd'} = 163.603610103856;
        $range{17}{'lower'} = -2.83791371717405;
        $range{17}{'upper'} = 4.99198685278779;
        $range{17}{'range'} = 7.82990056996184;
        $range{18}{'mean'} = -86.9080732700136;
        $range{18}{'sd'} = 145.736036521422;
        $range{18}{'lower'} = -2.80021288125261;
        $range{18}{'upper'} = 5.76321473616139;
        $range{18}{'range'} = 8.56342761741401;
        $range{19}{'mean'} = -80.3608095884215;
        $range{19}{'sd'} = 211.306064426292;
        $range{19}{'lower'} = -2.69106895703854;
        $range{19}{'upper'} = 3.6599082552985;
        $range{19}{'range'} = 6.35097721233704;
        $range{20}{'mean'} = -278.93061209106;
        $range{20}{'sd'} = 192.05331159166;
        $range{20}{'lower'} = -2.45801222585864;
        $range{20}{'upper'} = 8.10676264412129;
        $range{20}{'range'} = 10.5647748699799;
        $range{21}{'mean'} = -100;
        $range{21}{'sd'} = 1;
        $range{21}{'lower'} = 0;
        $range{21}{'upper'} = 0;
        $range{21}{'range'} = 1;
        $range{22}{'mean'} = -167.07315694256;
        $range{22}{'sd'} = 185.875809713378;
        $range{22}{'lower'} = -2.53355637715103;
        $range{22}{'upper'} = 6.16042054481576;
        $range{22}{'range'} = 8.69397692196679;
        $range{23}{'mean'} = -32768;
        $range{23}{'sd'} = 1;
        $range{23}{'lower'} = 0;
        $range{23}{'upper'} = 33219;
        $range{23}{'range'} = 33219;
        $range{24}{'mean'} = -32768;
        $range{24}{'sd'} = 1;
        $range{24}{'lower'} = 0;
        $range{24}{'upper'} = 32687;
        $range{24}{'range'} = 32687;
        $range{25}{'mean'} = -401.752826775215;
        $range{25}{'sd'} = 3.56447071169638;
        $range{25}{'lower'} = -11.8522991607579;
        $range{25}{'upper'} = 3.29721513397469;
        $range{25}{'range'} = 15.1495142947326;
        $range{26}{'mean'} = -32768;
        $range{26}{'sd'} = 1;
        $range{26}{'lower'} = 0;
        $range{26}{'upper'} = 32687;
        $range{26}{'range'} = 32687;
        $range{27}{'mean'} = -32768;
        $range{27}{'sd'} = 1;
        $range{27}{'lower'} = 0;
        $range{27}{'upper'} = 33113;
        $range{27}{'range'} = 33113;
	
	%aa = (
	'GLY'=>'G',
	'PRO'=>'P',
	'ALA'=>'A',
	'VAL'=>'V',
	'LEU'=>'L',
	'ILE'=>'I',
	'MET'=>'M',
	'CYS'=>'C',
	'PHE'=>'F',
	'TYR'=>'Y',
	'TRP'=>'W',
	'HIS'=>'H',
	'LYS'=>'K',
	'ARG'=>'R',
	'GLN'=>'Q',
	'ASN'=>'N',
	'GLU'=>'E',
	'ASP'=>'D',
	'SER'=>'S',
	'THR'=>'T'
	);

	%aa3 = (
	'G' => 'GLY',
	'P' => 'PRO',
	'A' => 'ALA',
	'V' => 'VAL',
	'L' => 'LEU',
	'I' => 'ILE',
	'M' => 'MET',
	'C' => 'CYS',
	'F' => 'PHE',
	'Y' => 'TYR',
	'W' => 'TRP',
	'H' => 'HIS',
	'K' => 'LYS',
	'R' => 'ARG',
	'Q' => 'GLN',
	'N' => 'ASN',
	'E' => 'GLU',
	'D' => 'ASP',
	'S' => 'SER',
	'T' => 'THR'
	);
	
}
