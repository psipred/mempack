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

double distance_rr(double, double, double, double);
double optimise_rotation(int);
