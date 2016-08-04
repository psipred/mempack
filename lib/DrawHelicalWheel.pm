package DrawHelicalWheel;
use strict;
use warnings;
use Math::Trig;
use Data::Dumper;

our @ISA = "GD"; 
our @EXPORT = qw(new png svg);
our $VERSION = '1.0';

sub new {

	my $class = shift;
 	my %options = @_;
	
	my $self = {

		## general paramaters
		'title' => exists $options{-title} ? $options{-title} : '',
		'sequence' => exists $options{-sequence} ? uc $options{-sequence} : 0,
		'write_sequence' => exists $options{-write_sequence} ? uc $options{-write_sequence} : 0,
		'helices' => \@{$options{-helices}},
		'ttf_font' => exists $options{-ttf_font} ? $options{-ttf_font} : 0,
		'ttf_font_size' => exists $options{-ttf_font_size} ? $options{-ttf_font_size} : 8,
		'helix_diameter' => exists $options{-helix_diameter} ? $options{-helix_diameter} : 400,
		'helix_spacing' => exists $options{-helix_spacing} ? $options{-helix_spacing} : 140,
		'residue_diameter' => exists $options{-residue_diameter} ? $options{-residue_diameter} : 44,
		'vertical_padding' => exists $options{-vertical_padding} ? $options{-vertical_padding} : 80,
		'horizontal_padding' => exists $options{-horizontal_padding} ? $options{-horizontal_padding} : 80,
		'angle' => exists $options{-angle} ? $options{-angle} : 100,
		'draw_centre' => exists $options{-draw_centre} ? $options{-draw_centre} : 1,
		'draw_spiral' => exists $options{-draw_spiral} ? $options{-draw_spiral} : 1,
		'label_helix' => exists $options{-label_helix} ? $options{-label_helix} : 1,
		'label_residues' => exists $options{-label_residues} ? $options{-label_residues} : 1,
		'label_helix_offset_row1' => exists $options{-label_helix_offset_row1} ? $options{-label_helix_offset_row1} : 36,
		'label_helix_offset_row2' => exists $options{-label_helix_offset_row2} ? $options{-label_helix_offset_row2} : 2,
		'residue_contacts' => \@{$options{-residue_contacts}},
		'rotations' => \@{$options{-rotations}},
		'only_these_helices' => \@{$options{-only_these_helices}},
		'aa_type' => \%{$options{-aa_type}},
		'hh_contacts' => \%{$options{-hh_contacts}},
		'graph_centres' => \%{$options{-graph_centres}},
		'custom_colour' => exists $options{-custom_colour} ? $options{-custom_colour} : [255,255,255]
	};

  	bless ($self,$class);
	
	%{$self->{'centres'}} = ();
	%{$self->{'residue_positions'}} = ();
	$self->{'_png'} = 0;
	$self->{'_radius'} = $self->{'helix_diameter'} / 2;

	$self->{'aa_type'}{'A'} = 'non-polar' unless exists $self->{'aa_type'}{'A'};
	$self->{'aa_type'}{'R'} = 'charged-positive' unless exists $self->{'aa_type'}{'R'};
	$self->{'aa_type'}{'N'} = 'polar' unless exists $self->{'aa_type'}{'N'};
	$self->{'aa_type'}{'D'} = 'charged-negative' unless exists $self->{'aa_type'}{'D'};
	$self->{'aa_type'}{'C'} = 'polar' unless exists $self->{'aa_type'}{'C'};
	$self->{'aa_type'}{'Q'} = 'polar' unless exists $self->{'aa_type'}{'Q'};
	$self->{'aa_type'}{'E'} = 'charged-negative' unless exists $self->{'aa_type'}{'E'};
	$self->{'aa_type'}{'G'} = 'non-polar' unless exists $self->{'aa_type'}{'G'};
	$self->{'aa_type'}{'H'} = 'polar' unless exists $self->{'aa_type'}{'H'};
	$self->{'aa_type'}{'I'} = 'non-polar' unless exists $self->{'aa_type'}{'I'};
	$self->{'aa_type'}{'L'} = 'non-polar' unless exists $self->{'aa_type'}{'L'};
	$self->{'aa_type'}{'K'} = 'charged-positive' unless exists $self->{'aa_type'}{'K'};
	$self->{'aa_type'}{'M'} = 'non-polar' unless exists $self->{'aa_type'}{'M'};
	$self->{'aa_type'}{'F'} = 'non-polar-aromatic' unless exists $self->{'aa_type'}{'F'};
	$self->{'aa_type'}{'P'} = 'non-polar' unless exists $self->{'aa_type'}{'P'};
	$self->{'aa_type'}{'S'} = 'polar' unless exists $self->{'aa_type'}{'S'};
	$self->{'aa_type'}{'T'} = 'polar' unless exists $self->{'aa_type'}{'T'};
	$self->{'aa_type'}{'W'} = 'non-polar-aromatic' unless exists $self->{'aa_type'}{'W'};
	$self->{'aa_type'}{'Y'} = 'polar-aromatic' unless exists $self->{'aa_type'}{'Y'};
	$self->{'aa_type'}{'V'} = 'non-polar' unless exists $self->{'aa_type'}{'V'};

	return $self;

}

sub do_checks {

	my $self = shift;

	## check to make sure we have pairs of helix boundaries and that data is numeric otherwise quit
	if (scalar @{$self->{'helices'}} % 2){
		die "\nUneven number of helix boundaries.\n\n";
	}

	foreach (@{$self->{'helices'}}){
		if ($_ =~ /\D/){
			die "\nTopology data is not numeric. $_\n\n";
		}
	}

	## check to make sure the TTF font exists, otherwise use gdSmallFont
	if ($self->{'ttf_font'}){
		unless (-e $self->{'ttf_font'}){
			print "\nCan't find font ".$self->{'ttf_font'}.".\n";
			$self->{'ttf_font'} = 0;
		}
	}

	unless ($self->{'sequence'}){
		die "\nNo sequence data found.\n\n";
	}

	unless (scalar keys %{$self->{'graph_centres'}}){
		die "\nNo helix centres found.\n\n";
	}

	$self->{'helix_count'} = scalar @{$self->{'helices'}} / 2;

	unless ($self->{'helix_count'}){
		push @{$self->{'helices'}},1;
		push @{$self->{'helices'}},length $self->{'sequence'};
		$self->{'helix_count'} = scalar @{$self->{'helices'}} / 2;
	}


	@{$self->{'helices'}} = sort {$a <=> $b} @{$self->{'helices'}};

	#$self->{'width'} = ($self->{'horizontal_padding'} * 2) + ($self->{'helix_diameter'} * $self->{'helix_count'}) + ($self->{'helix_spacing'} * ($self->{'helix_count'} - 1));
	#$self->{'height'} = $self->{'helix_diameter'} + ($self->{'vertical_padding'} * 2) + ((scalar keys %{$self->{'hh_contacts'}}) * 10);

}

sub allocate_colours {

	my $self = shift;

	$self->{'black'} = $self->{'im'}->colorAllocate(0,0,0);
	$self->{'white'} = $self->{'im'}->colorAllocate(255,255,255);
	$self->{'dark_blue'} = $self->{'im'}->colorAllocate(40,89,177);
	$self->{'medium_blue'} = $self->{'im'}->colorAllocate(54,106,207);
	$self->{'light_blue'} = $self->{'im'}->colorAllocate(70,130,255);
	$self->{'lightest_blue'} = $self->{'im'}->colorAllocate(70,155,255);	
	$self->{'dark_red'} = $self->{'im'}->colorAllocate(210,0,0);
	$self->{'medium_red'} = $self->{'im'}->colorAllocate(255,0,0);
	$self->{'light_red'} = $self->{'im'}->colorAllocate(255,33,33);
	$self->{'lightest_red'} = $self->{'im'}->colorAllocate(255,96,0);
	$self->{'dark_green'} = $self->{'im'}->colorAllocate(28,150,100);
	$self->{'medium_green'} = $self->{'im'}->colorAllocate(53,107,126);
	$self->{'light_green'} = $self->{'im'}->colorAllocate(74,205,150);
	$self->{'lightest_green'} = $self->{'im'}->colorAllocate(83,224,166);	
	$self->{'dark_purple'} = $self->{'im'}->colorAllocate(110,20,165);
	$self->{'medium_purple'} = $self->{'im'}->colorAllocate(125,45,180);
	$self->{'light_purple'} = $self->{'im'}->colorAllocate(150,70,210);
	$self->{'lightest_purple'} = $self->{'im'}->colorAllocate(170,90,230);
	$self->{'yellow'} = $self->{'im'}->colorAllocate(250,250,25);
	$self->{'custom'} = $self->{'im'}->colorAllocate(@{$self->{'custom_colour'}});	

}

sub png {

	use GD;

 	my $self = shift;
	
	## This stops us drawing ttf fonts for svg
	$self->{'_png'} = 1;

	$self->do_checks;
	$self->transform_centres;
	
	## create a new image
	$self->{'im'} = new GD::Image($self->{'width'},$self->{'height'});
	
	&allocate_colours($self);
	
	$self->{'im'}->filledRectangle(0,0,$self->{'width'},$self->{'height'},$self->{'white'});
	
	$self->{'small_font'} = GD::Font->Small;
	
	if ($self->{'ttf_font'} && $self->{'_png'}){
		$self->{'im'}->stringFT($self->{'black'},$self->{'ttf_font'},$self->{'ttf_font_size'},0,4,12,$self->{'title'},{linespacing=>0.6,charmap  => 'Unicode',}) if $self->{'title'};
	}else{
		$self->{'im'}->string($self->{'small_font'},4,3,$self->{'title'},$self->{'black'}) if $self->{'title'};
	}

	$self->get_residue_positions;
	$self->draw_hh_contacts;
	$self->draw_rr_contacts;
	$self->draw_helices;
	
	$self->{'_png'} = 0;

	return $self->{'im'}->GD::Image::png;

}

sub svg {

 	my $self = shift;

	unless ($self->{'_png'}){
		use GD::SVG;
	}

	$self->do_checks;
	$self->transform_centres;
	
	## create a new image
	$self->{'im'} = new GD::SVG::Image($self->{'width'},$self->{'height'});

	&allocate_colours($self);	

	$self->{'im'}->filledRectangle(0,0,$self->{'width'},$self->{'height'},$self->{'white'});

	$self->{'small_font'} = gdSmallFont;
	
	# Not sure how to use TTF fonts with SVG yet so we won't be going here...
	if ($self->{'ttf_font'} && $self->{'_png'}){
		$self->{'im'}->stringFT($self->{'black'},$self->{'ttf_font'},$self->{'ttf_font_size'},0,4,12,$self->{'title'},{linespacing=>0.6,charmap  => 'Unicode',}) if $self->{'title'};
	}else{
		$self->{'im'}->string($self->{'small_font'},4,3,$self->{'title'},$self->{'black'}) if $self->{'title'};
	}

	$self->get_residue_positions;
	$self->draw_hh_contacts;
	$self->draw_rr_contacts;
	$self->draw_helices;

	return $self->{'im'}->svg();

}

sub draw_rr_contacts {

	my $self = shift;
	my $mod = 0;

	foreach (@{$self->{'residue_contacts'}}){		
		
		my ($r1,$r2) = split(/-/,$_);
		$r1 =~ s/\s+//g;
		$r2 =~ s/\s+//g;
		
		if (exists ${$self->{'residue_positions'}}{$r1} && exists ${$self->{'residue_positions'}}{$r2}){
		
			$self->{'im'}->line(
			${$self->{'residue_positions'}}{$r1}[0],
			${$self->{'residue_positions'}}{$r1}[1],
			${$self->{'residue_positions'}}{$r2}[0],
			${$self->{'residue_positions'}}{$r2}[1],$self->{'black'});
		
		}else{
			#print "Undefined residue centres for $r1 and $r2\n";
		}
			
	}

	#my $radius_mod = $self->{'residue_diameter'}*2;
	#foreach (sort {$a <=> $b} keys %{$self->{'graph_centres'}}){		
	#	$self->{'im'}->filledArc(($self->{'helix_diameter'}/2)+${$self->{'graph_centres'}}{$_}[0],($self->{'helix_diameter'}/2)+${$self->{'graph_centres'}}{$_}[1],$radius_mod+$self->{'helix_diameter'},$radius_mod+$self->{'helix_diameter'},0,360,$self->{'white'});
	#}

}

sub draw_hh_contacts {

	my $self = shift;
	my $mod = 0;

	foreach (sort {$a cmp $b} keys %{$self->{'hh_contacts'}}){		
		
		my ($h1,$h2) = split(/-/,$_);
		$h1 =~ s/\s+//g;
		$h2 =~ s/\s+//g;

		if (defined ${$self->{'graph_centres'}}{$h1}[0] && defined ${$self->{'graph_centres'}}{$h2}[0] && defined ${$self->{'graph_centres'}}{$h1}[1] && defined ${$self->{'graph_centres'}}{$h2}[1]){
			$self->{'im'}->line(($self->{'helix_diameter'}/2)+${$self->{'graph_centres'}}{$h1}[0],($self->{'helix_diameter'}/2)+${$self->{'graph_centres'}}{$h1}[1],($self->{'helix_diameter'}/2)+${$self->{'graph_centres'}}{$h2}[0],($self->{'helix_diameter'}/2)+${$self->{'graph_centres'}}{$h2}[1],$self->{'medium_red'});
		}else{
			#print "Undefined hh centres for $h1 and $h2\n";
		}			
	}

	my $radius_mod = $self->{'residue_diameter'}*2;
	foreach (sort {$a <=> $b} keys %{$self->{'graph_centres'}}){	
	
		if (defined ${$self->{'graph_centres'}}{$_}[0] && defined ${$self->{'graph_centres'}}{$_}[1]){
			$self->{'im'}->filledArc(($self->{'helix_diameter'}/2)+${$self->{'graph_centres'}}{$_}[0],($self->{'helix_diameter'}/2)+${$self->{'graph_centres'}}{$_}[1],$radius_mod+$self->{'helix_diameter'},$radius_mod+$self->{'helix_diameter'},0,360,$self->{'white'});
		}
	}

}

sub get_residue_positions {

	my $self = shift;

	my $rotation_counter = 0;
	for (1..$self->{'helix_count'}){

		my $tag = 0;
		foreach my $h (@{$self->{'only_these_helices'}}){	
			$tag++ if $_ == $h;
		}
		next unless $tag;		
		
		#print "Drawing helix $_\n";

		my $x = ${$self->{'graph_centres'}}{$_}[0] + $self->{'_radius'};
		my $y = ${$self->{'graph_centres'}}{$_}[1] + $self->{'_radius'};

		my $start = ($_ * 2) - 2;
		my $start_pos = ${$self->{'helices'}}[$start];
		my $stop  = ($_ * 2) - 1;
		my $stop_pos = ${$self->{'helices'}}[$stop];
		my $length = ($stop_pos - $start_pos) + 1;

		my $helix = substr($self->{'sequence'},$start_pos - 1,$length);
		my $helix_length = length $helix;
		my ($mod,$last_x,$last_y,$shift) = 0;
		my $cx = $x;
		my $cy = $y;
		
		${$self->{'centres'}}{$_} = [$cx,$cy];
		
		my $t = -90;
		if (scalar @{$self->{'rotations'}}){
			$t += ${$self->{'rotations'}}[$rotation_counter];
		}

		for my $aa (1..$helix_length){

			$t = ($t - 360) if $t > 360;

			my $aa_cx = $cx + (($self->{'_radius'} + $mod) * cos deg2rad($t));
			my $aa_cy = $cy + (($self->{'_radius'} + $mod) * sin deg2rad($t));

	
			my $residue = substr($helix,$aa-1,1);
			
			my $real_pos = $start_pos + $aa - 1;	
			
			${$self->{'residue_positions'}}{$real_pos} = [$aa_cx,$aa_cy];
			
			$t += $self->{'angle'};
			$mod += $self->{'residue_diameter'} + 5 if $aa == 18;
			$mod += $self->{'residue_diameter'} + 5 if $aa == 36;
			$mod += $self->{'residue_diameter'} + 5 if $aa == 54;
			
		}
		$rotation_counter++;
	}

	return $self;
}

sub draw_helices {

	my $self = shift;

	my $rotation_counter = 0;
	for (1..$self->{'helix_count'}){

		my $tag = 0;
		foreach my $h (@{$self->{'only_these_helices'}}){	
			$tag++ if $_ == $h;
		}
		next unless $tag;		
		
		#print "Drawing helix $_\n";

		my $x = ${$self->{'graph_centres'}}{$_}[0] + $self->{'_radius'};
		my $y = ${$self->{'graph_centres'}}{$_}[1] + $self->{'_radius'};

		my $start = ($_ * 2) - 2;
		my $start_pos = ${$self->{'helices'}}[$start];
		my $stop  = ($_ * 2) - 1;
		my $stop_pos = ${$self->{'helices'}}[$stop];
		my $length = ($stop_pos - $start_pos) + 1;

		my $helix = substr($self->{'sequence'},$start_pos - 1,$length);
		my $helix_length = length $helix;
		my ($mod,$last_x,$last_y,$shift) = 0;
		my $cx = $x;
		my $cy = $y;
		
		$self->{'im'}->arc($cx,$cy,5,5,0,360,$self->{'black'}) if $self->{'draw_centre'};
		
		${$self->{'centres'}}{$_} = [$cx,$cy];
		
		if ($self->{'label_helix'}){

			my $label = "Helix $_ : $start_pos - $stop_pos";
			my $rect_height = 7;
			my $rect_width = 40;

			$self->{'im'}->filledRectangle($cx-$rect_width-1,$cy-$rect_height-15,$cx+$rect_width+4,$cy+$rect_height-15,$self->{'white'});

			if ($self->{'ttf_font'} && $self->{'_png'}){
				$self->{'im'}->stringFT($self->{'black'},$self->{'ttf_font'},$self->{'ttf_font_size'},0,$cx - $self->{'label_helix_offset_row1'},$cy - 10,$label,{linespacing=>0.6,charmap  => 'Unicode',});
			}elsif($self->{'_png'}){
				$self->{'im'}->string($self->{'small_font'},$cx - $self->{'label_helix_offset_row1'} - 14,$cy - 20,$label,$self->{'black'});
			}else{
				$self->{'im'}->string($self->{'small_font'},$cx - $self->{'label_helix_offset_row1'},$cy - 20,$label,$self->{'black'});
			}
		
			# Write sequence
			if ($self->{'write_sequence'}){
				if ($length > 30){
		
					my $text_mod = ((30 * 6.5) / 2) - $self->{'label_helix_offset_row2'};

					if ($self->{'ttf_font'} && $self->{'_png'}){
						$self->{'im'}->stringFT($self->{'black'},$self->{'ttf_font'},$self->{'ttf_font_size'},0,$cx - $text_mod,$cy + 20,substr($helix,0,30),{linespacing=>0.6,charmap  => 'Unicode',});
					}else{
						$self->{'im'}->string($self->{'small_font'},$cx - $text_mod,$cy + 10,substr($helix,0,30),$self->{'black'});
					}

					$text_mod = (($length - 30) * 6.5) / 2;

					if ($self->{'ttf_font'} && $self->{'_png'}){
						$self->{'im'}->stringFT($self->{'black'},$self->{'ttf_font'},$self->{'ttf_font_size'},0,$cx - $text_mod,$cy + 30,substr($helix,30,$length - 30),{linespacing=>0.6,charmap  => 'Unicode',});
					}else{
						$self->{'im'}->string($self->{'small_font'},$cx - $text_mod,$cy + 20,substr($helix,30,$length - 30),$self->{'black'});
					}

				}else{
					my $text_mod = (($length * 6.5) / 2) - $self->{'label_helix_offset_row2'};

					if ($self->{'ttf_font'} && $self->{'_png'}){
						$self->{'im'}->stringFT($self->{'black'},$self->{'ttf_font'},$self->{'ttf_font_size'},0,$cx - $text_mod,$cy + 20,$helix,{linespacing=>0.6,charmap  => 'Unicode',});
					}else{
						$self->{'im'}->string($self->{'small_font'},$cx - $text_mod,$cy + 10,$helix,$self->{'black'});
					}		
				}	
			}	
		}

		my $t = -90;
		if (scalar @{$self->{'rotations'}}){
			#print "Adding ${$self->{'rotations'}}[$_-1] degrees to helix rotation.\n";
			$t += ${$self->{'rotations'}}[$rotation_counter];
		}
		
		if ($self->{'draw_spiral'}){
		
			for my $aa (1..$helix_length){
		
				$t = ($t - 360) if $t > 360;

				my $aa_cx = $cx + ($self->{'_radius'} * cos deg2rad($t));
				my $aa_cy = $cy + ($self->{'_radius'} * sin deg2rad($t));
			
				if ($last_x && $last_y){	

					$self->{'im'}->line($last_x,$last_y,$aa_cx,$aa_cy,$self->{'black'});

					$last_x = $aa_cx;
					$last_y = $aa_cy;
				
				}else{
					$last_x = $aa_cx;
					$last_y = $aa_cy;
				
				}

				$t += $self->{'angle'};
			}		
		}

		$t = -90;
		if (scalar @{$self->{'rotations'}}){
			$t += ${$self->{'rotations'}}[$rotation_counter];
		}

		for my $aa (1..$helix_length){

			$t = ($t - 360) if $t > 360;

			my $aa_cx = $cx + (($self->{'_radius'} + $mod) * cos deg2rad($t));
			my $aa_cy = $cy + (($self->{'_radius'} + $mod) * sin deg2rad($t));

			## Draw residues
			$self->{'im'}->filledArc($aa_cx,$aa_cy,$self->{'residue_diameter'},$self->{'residue_diameter'},0,360,$self->{'black'});
	
			my $residue = substr($helix,$aa-1,1);
			
			my $real_pos = $start_pos + $aa - 1;
	
			my $tmp_colour1 = 0;
			my $tmp_colour2 = 0;
			my $tmp_colour3 = 0;
			my $tmp_colour4 = 0;
			
			if (!exists $self->{'aa_type'}{$residue}){
				&draw_other($self,$aa_cx,$aa_cy);
			}elsif ($self->{'aa_type'}{$residue} eq 'polar'){
				&draw_polar($self,$aa_cx,$aa_cy);
			}elsif ($self->{'aa_type'}{$residue} eq 'non-polar'){
				&draw_non_polar($self,$aa_cx,$aa_cy);
			}elsif ($self->{'aa_type'}{$residue} eq 'charged-positive'){
				&draw_charged_positive($self,$aa_cx,$aa_cy);
			}elsif ($self->{'aa_type'}{$residue} eq 'charged-negative'){
				&draw_charged_negative($self,$aa_cx,$aa_cy);
			}elsif ($self->{'aa_type'}{$residue} eq 'polar-aromatic'){
				&draw_polar_aromatic($self,$aa_cx,$aa_cy);
			}elsif ($self->{'aa_type'}{$residue} eq 'non-polar-aromatic'){
				&draw_non_polar_aromatic($self,$aa_cx,$aa_cy);
			}elsif ($self->{'aa_type'}{$residue} eq 'custom-colour'){
				&draw_custom_colour($self,$aa_cx,$aa_cy);
			}else{
				&draw_other($self,$aa_cx,$aa_cy);
			}
				
			$self->{'im'}->filledArc($aa_cx,$aa_cy,$self->{'residue_diameter'}-26,$self->{'residue_diameter'}-26,0,360,$self->{'white'});
			
			if ($residue eq 'I'){
				$shift = 2;
			}elsif($residue eq 'V' || $residue eq 'C' || $residue eq 'G' || $residue eq 'M'){
				$shift = -1;
			}elsif($residue eq 'W'){
				$shift = -2;
			}elsif($residue eq 'L'){
				$shift = 1;
			}else{
				$shift = 0;
			}
	
			if ($self->{'label_residues'}){
	
				if ($self->{'ttf_font'} && $self->{'_png'}){
					$self->{'im'}->stringFT($self->{'black'},$self->{'ttf_font'},$self->{'ttf_font_size'},0,$aa_cx - 3 + $shift,$aa_cy + 4,$residue,{linespacing=>0.6,charmap  => 'Unicode',});
				}else{
					$self->{'im'}->string($self->{'small_font'},$aa_cx - 3 + $shift,$aa_cy - 6,$residue,$self->{'black'});
				}

				if ($real_pos < 10){

					if ($self->{'ttf_font'} && $self->{'_png'}){
						$self->{'im'}->stringFT($self->{'black'},$self->{'ttf_font'},$self->{'ttf_font_size'},0,$aa_cx - 2,$aa_cy + 19,$real_pos,{linespacing=>0.6,charmap  => 'Unicode',});
					}else{
						$self->{'im'}->string($self->{'small_font'},$aa_cx - 2,$aa_cy + 8,$real_pos,$self->{'black'});
					}
				
				}elsif ($real_pos < 100){
				
					if ($self->{'ttf_font'} && $self->{'_png'}){
						$self->{'im'}->stringFT($self->{'black'},$self->{'ttf_font'},$self->{'ttf_font_size'},0,$aa_cx - 5,$aa_cy + 19,$real_pos,{linespacing=>0.6,charmap  => 'Unicode',});
					}else{
						$self->{'im'}->string($self->{'small_font'},$aa_cx - 5,$aa_cy + 8,$real_pos,$self->{'black'});
					}
				}else{
					if ($self->{'ttf_font'} && $self->{'_png'}){
						$self->{'im'}->stringFT($self->{'black'},$self->{'ttf_font'},$self->{'ttf_font_size'},0,$aa_cx - 8,$aa_cy + 19,$real_pos,{linespacing=>0.6,charmap  => 'Unicode',});
					}else{
						$self->{'im'}->string($self->{'small_font'},$aa_cx - 8,$aa_cy + 8,$real_pos,$self->{'black'});
					}
				}		
			}
						
			#${$self->{'residue_positions'}}{$real_pos} = [$aa_cx,$aa_cy];
			
			$t += $self->{'angle'};
			$mod += $self->{'residue_diameter'} + 5 if $aa == 18;
			$mod += $self->{'residue_diameter'} + 5 if $aa == 36;
			$mod += $self->{'residue_diameter'} + 5 if $aa == 54;
			
		}
		$rotation_counter++;
	}

	return $self;
}

sub transform_centres{

	my $self = shift;
	my $x_low = 100000000;
	my $y_low = 100000000;
	my $x_high = -100000000;
	my $y_high = -100000000;
	my $x_mod = 0;
	my $y_mod = 0;
		
	foreach (keys %{$self->{'graph_centres'}}){	
		$x_low = ${$self->{'graph_centres'}}{$_}[0] if ${$self->{'graph_centres'}}{$_}[0] < $x_low;
		$y_low = ${$self->{'graph_centres'}}{$_}[1] if ${$self->{'graph_centres'}}{$_}[1] < $y_low;	
		$x_high = ${$self->{'graph_centres'}}{$_}[0] if ${$self->{'graph_centres'}}{$_}[0] > $x_high;
		$y_high = ${$self->{'graph_centres'}}{$_}[1] if ${$self->{'graph_centres'}}{$_}[1] > $y_high;
	}

	$x_mod = -$x_low;
	$y_mod = -$y_low;	
	
	my $max_dimension = 1500;
	my $max_x = 0;
	my $max_y = 0;
	foreach (sort {$a <=> $b} keys %{$self->{'graph_centres'}}){	
		${$self->{'graph_centres'}}{$_} = [${$self->{'graph_centres'}}{$_}[0]+$x_mod+$self->{'horizontal_padding'}+($self->{'helix_diameter'}/2),${$self->{'graph_centres'}}{$_}[1]+$y_mod+$self->{'vertical_padding'}+($self->{'helix_diameter'}/2)];
		if (${$self->{'graph_centres'}}{$_}[0] > $max_dimension){
			$max_x = ${$self->{'graph_centres'}}{$_}[0];
		}

		if (${$self->{'graph_centres'}}{$_}[1] > $max_dimension){
			$max_y = ${$self->{'graph_centres'}}{$_}[1];
		}
		#print "$_ => ${$self->{'graph_centres'}}{$_}[0],${$self->{'graph_centres'}}{$_}[1]\n";
	}	

	$self->{'width'} = ($x_high - $x_low) + ($self->{'horizontal_padding'} * 1) + ($self->{'helix_diameter'} * 2);
	$self->{'height'} = ($y_high - $y_low) + ($self->{'vertical_padding'} * 1) + ($self->{'helix_diameter'} * 2);
}

sub draw_polar{

	my ($self,$aa_cx,$aa_cy) = @_;
	$self->{'im'}->filledArc($aa_cx,$aa_cy,$self->{'residue_diameter'}-3,$self->{'residue_diameter'}-3,0,360,$self->{'dark_red'});
	$self->{'im'}->filledArc($aa_cx,$aa_cy,$self->{'residue_diameter'}-5,$self->{'residue_diameter'}-5,0,360,$self->{'medium_red'});
	$self->{'im'}->filledArc($aa_cx,$aa_cy,$self->{'residue_diameter'}-8,$self->{'residue_diameter'}-8,0,360,$self->{'light_red'});
	$self->{'im'}->filledArc($aa_cx,$aa_cy,$self->{'residue_diameter'}-10,$self->{'residue_diameter'}-10,0,360,$self->{'lightest_red'});

}

sub draw_non_polar{

	my ($self,$aa_cx,$aa_cy) = @_;
	$self->{'im'}->filledArc($aa_cx,$aa_cy,$self->{'residue_diameter'}-3,$self->{'residue_diameter'}-3,0,360,$self->{'dark_blue'});
	$self->{'im'}->filledArc($aa_cx,$aa_cy,$self->{'residue_diameter'}-5,$self->{'residue_diameter'}-5,0,360,$self->{'medium_blue'});
	$self->{'im'}->filledArc($aa_cx,$aa_cy,$self->{'residue_diameter'}-8,$self->{'residue_diameter'}-8,0,360,$self->{'light_blue'});
	$self->{'im'}->filledArc($aa_cx,$aa_cy,$self->{'residue_diameter'}-10,$self->{'residue_diameter'}-10,0,360,$self->{'lightest_blue'});

}

sub draw_charged_positive{

	my ($self,$aa_cx,$aa_cy) = @_;
	$self->{'im'}->filledArc($aa_cx,$aa_cy,$self->{'residue_diameter'}-3,$self->{'residue_diameter'}-3,0,360,$self->{'dark_purple'});
	$self->{'im'}->filledArc($aa_cx,$aa_cy,$self->{'residue_diameter'}-5,$self->{'residue_diameter'}-5,0,360,$self->{'medium_purple'});
	$self->{'im'}->filledArc($aa_cx,$aa_cy,$self->{'residue_diameter'}-8,$self->{'residue_diameter'}-8,0,360,$self->{'light_purple'});
	$self->{'im'}->filledArc($aa_cx,$aa_cy,$self->{'residue_diameter'}-10,$self->{'residue_diameter'}-10,0,360,$self->{'lightest_purple'});

}

sub draw_charged_negative{

	my ($self,$aa_cx,$aa_cy) = @_;
	$self->{'im'}->filledArc($aa_cx,$aa_cy,$self->{'residue_diameter'}-3,$self->{'residue_diameter'}-3,0,360,$self->{'dark_green'});
	$self->{'im'}->filledArc($aa_cx,$aa_cy,$self->{'residue_diameter'}-5,$self->{'residue_diameter'}-5,0,360,$self->{'medium_green'});
	$self->{'im'}->filledArc($aa_cx,$aa_cy,$self->{'residue_diameter'}-8,$self->{'residue_diameter'}-8,0,360,$self->{'light_green'});
	$self->{'im'}->filledArc($aa_cx,$aa_cy,$self->{'residue_diameter'}-10,$self->{'residue_diameter'}-10,0,360,$self->{'lightest_green'});

}

sub draw_polar_aromatic{

	my ($self,$aa_cx,$aa_cy) = @_;
	$self->{'im'}->filledArc($aa_cx,$aa_cy,$self->{'residue_diameter'}-3,$self->{'residue_diameter'}-3,0,360,$self->{'dark_red'});
	$self->{'im'}->filledArc($aa_cx,$aa_cy,$self->{'residue_diameter'}-5,$self->{'residue_diameter'}-5,0,360,$self->{'yellow'});
	$self->{'im'}->filledArc($aa_cx,$aa_cy,$self->{'residue_diameter'}-8,$self->{'residue_diameter'}-8,0,360,$self->{'light_red'});
	$self->{'im'}->filledArc($aa_cx,$aa_cy,$self->{'residue_diameter'}-10,$self->{'residue_diameter'}-10,0,360,$self->{'lightest_red'});

}

sub draw_non_polar_aromatic{

	my ($self,$aa_cx,$aa_cy) = @_;
	$self->{'im'}->filledArc($aa_cx,$aa_cy,$self->{'residue_diameter'}-3,$self->{'residue_diameter'}-3,0,360,$self->{'dark_blue'});
	$self->{'im'}->filledArc($aa_cx,$aa_cy,$self->{'residue_diameter'}-5,$self->{'residue_diameter'}-5,0,360,$self->{'yellow'});
	$self->{'im'}->filledArc($aa_cx,$aa_cy,$self->{'residue_diameter'}-8,$self->{'residue_diameter'}-8,0,360,$self->{'light_blue'});
	$self->{'im'}->filledArc($aa_cx,$aa_cy,$self->{'residue_diameter'}-10,$self->{'residue_diameter'}-10,0,360,$self->{'lightest_blue'});

}

sub draw_custom_colour{

	my ($self,$aa_cx,$aa_cy) = @_;
	$self->{'im'}->filledArc($aa_cx,$aa_cy,$self->{'residue_diameter'}-3,$self->{'residue_diameter'}-3,0,360,$self->{'custom'});

}

sub draw_other{

	my ($self,$aa_cx,$aa_cy) = @_;
	$self->{'im'}->filledArc($aa_cx,$aa_cy,$self->{'residue_diameter'}-3,$self->{'residue_diameter'}-3,0,360,$self->{'white'});

}

1;

=head1 NAME

LOCAL::lib::DrawHelicalWheel.pm - draw a helical wheel.

=head1 SYNOPSIS


use DrawHelicalWheel;

## Simple use - -sequence is the only option that is required

my $im = DrawHelicalWheel->new(-title=>'An amphipathic alpha-helix',
                               -sequence=>'EIVRSIVKTMQNVATQAP');


## Write to a .png file

open(OUTPUT, ">amphipathic_helix.png");

binmode OUTPUT;

print OUTPUT $im->png;

close OUTPUT;



## Write to a .svg file
open(OUTPUT, ">amphipathic_helix.svg");

binmode OUTPUT;

print OUTPUT $im->svg;

close OUTPUT;


## More advanced use - draw multiple helices

my $title = 'Aquaporin-4 (AQP-4) water channel';

my $sequence = 'MSDGAAARRWGKCGPPCSRESIMVAFKGVWTQAFWKAVTAEFLAMLIFVLLSVGSTINWGGSENPLPVDMVLISLCFGLSIATMVQCFGHISGGHINPAVTVAMVCTRKISIAKSVFYITAQCLGAIIGAGILYLVTPPSVVGGLGVTTVHGNLTAGHGLLVELIITFQLVFTIFASCDSKRTDVTGSVALAIGFSVAIGHLFAINYTGASMNPARSFGPAVIMGNWENHWIYWVGPIIGAVLAGALYEYVFCPDVELKRRLKEAFSKAAQQTKGSYMEVEDNRSQVETEDLILKPGVVHVIDIDRGDEKKGKDSSGEVLSSV';

my @helices = (34,56,70,88,136,178,189,203,231,252);

my $font = '/usr/share/fonts/msttcorefonts/arial.ttf';

my $im = DrawHelicalWheel->new(-title=>$title,
                               -sequence=>$sequence,
                               -helices=>\@helices,
                               -ttf_font=>$font);


## Write to a .png file

open(OUTPUT, ">aquaporin.png");

binmode OUTPUT;

print OUTPUT $im->png;

close OUTPUT;


=head1 DESCRIPTION


A module to draw a helical wheel. It uses GD and/or GD:SVG to write to a .png or .svg file.

The options are a set of tag/value pairs as follows:

  Option                      Value                                         Default
  ------                      -----                                         -------

  -title                      Title to add to the image                     none
  
  -sequence                   The protein sequence of the helix, or the     none
                              whole protein (if you give it the whole
                              protein you should also use the
                              -helices tag

  -write_sequence             Write the helix sequence in the helix         0
		      		      		      
  -helices                    An array containing pairs of helix            none
                              boundaries

  -ttf_font                   Path to TTF font, e.g.                        none 
                              /usr/share/fonts/msttcorefonts/arial.ttf
                              Only implemented for PNG files at the moment		      
		      
  -ttf_font_size              Default size for TTF font. Use 7-9 with       8
                              Arial for best results  
		      
  -helix_diameter             Diameter of each helix                        400
  
  -helix_spacing              Spacing between each helix                    140                           

  -residue_diameter           Diameter of each residue                      44

  -vertical_padding           Vertical padding                              130

  -horizontal_padding         Horizontal Padding                            80

  -angle                      Angle between each residue.                   100
                              e.g. 100 degrees = 3.6 residues/revolution  

  -draw_spiral                Draw the spiral in the centre of the helix    1
  
  -label_helix                Label the centre of the helix                 1
  
  -label_residues             Label the residues                            1
  
  -label_helix_offset_row1    Text offset for helix label row 1             40
  
  -label_helix_offset_row2    Text offset for helix label row 2             4

  -aa_type                    Hash ref containing residue and type, e.g.    none
                              %aa_type = ('A' = 'non-polar', 
			                  'R' = 'charged-positive',
                                          'N' = 'polar',
                                          'E' = 'charged-negative',
                                          'W' = 'non-polar-aromatic'
                                          'Y' = 'polar-aromatic',
                                          'K' = 'custom-colour');

                              Use this hash to change the type and thus
                              colour of a residue. Any other type will
                              be white. custom-colour type will be set
                              to the -custom_colour RGB code.

  -rotations                  Array ref containing degrees to rotate each   none
                              helix by.

  -custom_colour              Array ref containing RGB code of              [255,255,255]
                              custom-colour residues		      		      			      			      					  

	
=head1 AUTHOR

Tim Nugent E<lt>timnugent@gmail.comE<gt>

=cut
