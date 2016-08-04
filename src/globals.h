// ********************************************************
// *   MEMPACK - Predicting transmembrane helix packing   *
// *       arrangements using residue contacts and        *
// *              a force-directed algorithm.             *
// * Copyright (C) 2009 Timothy Nugent and David T. Jones *
// ********************************************************
// 
// This program is copyright and may not be distributed without
// permission of the author unless specifically permitted under
// the terms of the license agreement.
// 
// THIS SOFTWARE MAY ONLY BE USED FOR NON-COMMERCIAL PURPOSES. PLEASE CONTACT
// THE AUTHOR IF YOU REQUIRE A LICENSE FOR COMMERCIAL USE.
// Based on code by Timothy Nugent: November 5th 2009
// 
// Description: This program uses a force-directed algorithm
// to predict the helical packing arrangement of transmembrane 
// proteins based on multiple sequence profiles.
// 


#ifndef GLOBALS_H
#define GLOBALS_H

#include <map>
#include <vector>
#include <string>

using namespace std;

typedef map<const pair<int,int>, const pair<int,int> > rr_hh_map;
typedef map<const int, pair<double,double> > r_xy_map;
extern unsigned int total;
extern map<const int, pair<double,double> > helix_positions;
extern vector< map<const int, pair<double,double> > > all_helix_positions;
extern map<const string,int> helix_swaps_seen;
extern vector<int> boundaries;
extern vector<vector<int> > all_rotations;
extern rr_hh_map residues_to_helices;
extern vector<r_xy_map> all_residue_positions;

#endif
