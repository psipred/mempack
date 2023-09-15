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

#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/random_layout.hpp>
#include <boost/graph/circle_layout.hpp>
#include "kamada_kawai_spring_layout.h"
#include <boost/random/linear_congruential.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/simple_point.hpp>
#include <boost/graph/connected_components.hpp> 
#include "globals.h"
#include "paramopt.h"

using namespace boost;
using namespace std;

bool verbose = false;

typedef adjacency_list<vecS, vecS, undirectedS, no_property, property<edge_weight_t, double> > Graph;
typedef graph_traits<Graph>::vertex_descriptor Vertex;
typedef std::map<std::string, Vertex> NameToVertex;
typedef property_map<Graph, edge_weight_t>::type WeightMap;
typedef std::vector<simple_point<double> > PositionVec;
typedef iterator_property_map<PositionVec::iterator, property_map<Graph, vertex_index_t>::type> PositionMap;

static void tokenize(const string str, vector<int>& tokens){

    	const string delimiters = ",";
    
	// Skip delimiters at beginning
	string::size_type lastPos = str.find_first_not_of(delimiters, 0);
	
	// Find first non-delimiter
	string::size_type pos = str.find_first_of(delimiters, lastPos);

 	while (string::npos != pos || string::npos != lastPos){
	
        	// Found a token, add it to the vector
        	//tokens.push_back(str.substr(lastPos, pos - lastPos));
		int i;
		istringstream myStream(str.substr(lastPos, pos - lastPos));
		myStream >> i;
		tokens.push_back(i);

        	// Skip delimiters.  Note the not_of
        	lastPos = str.find_first_not_of(delimiters, pos);
		
        	// Find next non-delimiter
        	pos = str.find_first_of(delimiters, lastPos);
    	}
	
}

void get_residue_positions(double rotate, int helix, int x){

	double helix_x = all_helix_positions[x][helix].first;
	double helix_y = all_helix_positions[x][helix].second;
	int seq_start = boundaries[helix*2];
	int seq_stop = boundaries[(helix*2)+1];
	double radius = 5.0;
	double angle = -90.0 + rotate;

	for (int i = seq_start; i <= seq_stop; i++){	
		double aa_x = helix_x + radius * cos (angle * (M_PI/180));
		double aa_y = helix_y + radius * sin (angle * (M_PI/180));	
		//cout << i << ":\t" << aa_x << "," << aa_y << endl; 
		pair <double,double> pos (aa_x,aa_y);
		all_residue_positions[x][i] = pos;
		angle += 100.0; 		
		//cout << i << ":\t" << all_residue_positions[x][i].first << "," << residue_positions[i].second << endl; 		
	}
}

double distance_rr(double x1, double y1, double x2, double y2){
	double d1 = x1 - x2;
	double d2 = y1 - y2;
	return(sqrt(d1*d1+d2*d2));
}

double optimise_rotation(int x){

	double total_distance = 0;

	for (unsigned int i = 0; i < total; i++){
		//cout << "Rotating helix " << i+1 << " by " << all_rotations[x][i] << " degrees..." << endl;	
		get_residue_positions(all_rotations[x][i], i, x);
	}

    	for(rr_hh_map::const_iterator it = residues_to_helices.begin(); it != residues_to_helices.end(); ++it){
        	/*
		cout << "Residues: " << it->first.first << "-" << it->first.second << endl;
		cout << "r1 x: " << residue_positions[it->first.first].first << endl;
		cout << "r1 y: " << residue_positions[it->first.first].second << endl;
		cout << "r2 x: " << residue_positions[it->first.second].first << endl;
		cout << "r2 y: " << residue_positions[it->first.second].second << endl;
		*/

		double this_distance = distance_rr(all_residue_positions[x][it->first.first].first,all_residue_positions[x][it->first.first].second,all_residue_positions[x][it->first.second].first,all_residue_positions[x][it->first.second].second);

		total_distance += this_distance;

    	}	

	return(total_distance);

}

bool lineSegmentIntersection(double Ax, double Ay,double Bx, double By,double Cx, double Cy,double Dx, double Dy) {

	double  distAB, theCos, theSin, newX, ABpos ;

	// Fail if either line segment is zero-length.
	if (Ax==Bx && Ay==By || Cx==Dx && Cy==Dy){
		return false;
	}

  	// Fail if the segments share an end-point.
	if (Ax==Cx && Ay==Cy || Bx==Cx && By==Cy ||  Ax==Dx && Ay==Dy || Bx==Dx && By==Dy){
		return false;
	}

	// Translate the system so that point A is on the origin.
	Bx-=Ax; By-=Ay;
	Cx-=Ax; Cy-=Ay;
	Dx-=Ax; Dy-=Ay;

	//  Discover the length of segment A-B.
	distAB=sqrt(Bx*Bx+By*By);

	//  Rotate the system so that point B is on the positive X axis.
	theCos=Bx/distAB;
	theSin=By/distAB;
	newX=Cx*theCos+Cy*theSin;
	Cy  =Cy*theCos-Cx*theSin; Cx=newX;
	newX=Dx*theCos+Dy*theSin;
	Dy  =Dy*theCos-Dx*theSin; Dx=newX;

	// Fail if segment C-D doesn't cross line A-B.
	if (Cy<0. && Dy<0. || Cy>=0. && Dy>=0.){
		return false;
	}

	// Discover the position of the intersection point along line A-B.
	ABpos=Dx+(Cx-Dx)*Dy/(Dy-Cy);

	//  Fail if segment C-D crosses line A-B outside of segment A-B.
	if (ABpos<0. || ABpos>distAB){
		return false;
	}

	//  Success.
	return true;
}

int count_crossovers(unsigned int total, int x){

	vector<vector<double> > lines;
	int cross_overs = 0;

	//cout << "First side..." << endl;
	for (unsigned int i= 0; i < total; i+=2){	
		if ((i+1) < total){
			if(verbose) cout << "Loop between helix " << i+1 << " and " << i+2 << endl;	
			if(verbose) cout << all_helix_positions[x][i].first << "," << all_helix_positions[x][i].second << endl;
			if(verbose) cout << all_helix_positions[x][i+1].first << "," << all_helix_positions[x][i+1].second << endl;
			vector<double> line;
			line.push_back(all_helix_positions[x][i].first);
			line.push_back(all_helix_positions[x][i].second);
			line.push_back(all_helix_positions[x][i+1].first);
			line.push_back(all_helix_positions[x][i+1].second);
			lines.push_back(line);
		}
	}
	//cout << endl;
	for (unsigned int i= 0; i < lines.size(); i++){
		for (unsigned int j= i+1; j < lines.size(); j++){
			if(verbose) cout << "Line " << i+1 << "\t" << lines[i][0] << "," << lines[i][1]  << " - "<< lines[i][2] << "," << lines[i][3] << endl;
			if(verbose) cout << "Line " << j+1 << "\t" << lines[j][0] << "," << lines[j][1]  << " - "<< lines[j][2] << "," << lines[j][3] << endl;			
			if (lineSegmentIntersection(lines[i][0],lines[i][1],lines[i][2],lines[i][3],lines[j][0],lines[j][1],lines[j][2],lines[j][3])){
				if(verbose) cout << "Lines intersect." << endl;
				cross_overs++;			
			}else{
				if(verbose) cout << "Lines do not intersect!" << endl;
			}			
		}

	}
	lines.clear();
	//cout << endl;

	//cout << "Second side..." << endl;
	for (unsigned int i= 1; i < total; i+=2){
		if ((i+1) < total){
			if(verbose) cout << "Loop between helix " << i+1 << " and " << i+2 << endl;	
			if(verbose) cout << all_helix_positions[x][i].first << "," << all_helix_positions[x][i].second << endl;
			if(verbose) cout << all_helix_positions[x][i+1].first << "," << all_helix_positions[x][i+1].second << endl;
			vector<double> line;
			line.push_back(all_helix_positions[x][i].first);
			line.push_back(all_helix_positions[x][i].second);
			line.push_back(all_helix_positions[x][i+1].first);
			line.push_back(all_helix_positions[x][i+1].second);
			lines.push_back(line);
		}	
	}
	//cout << endl;
	for (unsigned int i= 0; i < lines.size(); i++){
		for (unsigned int j= i+1; j < lines.size(); j++){
			//cout << "Line " << i+1 << "\t" << lines[i][0] << "," << lines[i][1]  << " - "<< lines[i][2] << "," << lines[i][3] << endl;
			//cout << "Line " << j+1 << "\t" << lines[j][0] << "," << lines[j][1]  << " - "<< lines[j][2] << "," << lines[j][3] << endl;			
			if (lineSegmentIntersection(lines[i][0],lines[i][1],lines[i][2],lines[i][3],lines[j][0],lines[j][1],lines[j][2],lines[j][3])){
				if(verbose) cout << "Lines intersect." << endl;
				cross_overs++;			
			}else{
				if(verbose) cout << "Lines do not intersect!" << endl;
			}			
		}

	}
	lines.clear();
	//cout << endl;
	
	return(cross_overs);


}

void swap_helices(int total, vector<int>& h1, vector<int>& h2, int x){	

	// Loop through first helix
	for (int h = 0; h < total; h++){
		//cout << h1[i] << "--->" << h2[i] << endl;

		vector<int> h1_contact,h2_contact;

		for (unsigned int i = 0; i < h1.size(); i++){
			if(h1[i] == h) h1_contact.push_back(h2[i]);
			if(h2[i] == h) h1_contact.push_back(h1[i]);
		}	
			
		// Loop through second helix
		for (int j = h+1; j < total; j++){		

			for (unsigned int i = 0; i < h1.size(); i++){
				if(h1[i] == j) h2_contact.push_back(h2[i]);
				if(h2[i] == j) h2_contact.push_back(h1[i]);
			}			
		
			if (h1_contact.size() == h2_contact.size() && h1_contact.size()){

				// Make sure we skip these two helices when we recurse
				string seen("");
				stringstream n,m;
				n << h+1;
				seen.append(n.str());
				seen.append("-");
				m << j+1;
				seen.append(m.str());
				
				int required_score = h1_contact.size();
				int correct = 0;
				int seen_j = 0;
				int seen_h = 0;
				
				if(!helix_swaps_seen[seen]){
				
					for (unsigned int k = 0; k < h1_contact.size(); k++){
						for (unsigned int l = 0; l < h2_contact.size(); l++){
							if (h1_contact[k] == h2_contact[l]){
								correct++;
							}else if(h1_contact[k] == j){
								if(!seen_j){
									seen_j++;
									correct++;
								}							
							}else if(h2_contact[l] == h){
								if(!seen_h){
									seen_h++;
									correct++;	
								}					
							}						
						}
					}				
				}else{
					//cout << "Already seen these two helices " << seen << "..." << endl;
				}
				
				if (seen_h && seen_j) required_score++;
				if (seen_h && seen_j && h1_contact.size() == 1){
					correct = -1;
				}	
				if (correct == required_score){
			
					// Make sure we skip these two helices when we recurse
					helix_swaps_seen[seen] = 1;
					map<const string,int> helix_swaps_seen_copy = helix_swaps_seen;
		
					// These helices have intechangable positions
					// Work out same-side loop crossover scores					

					if(verbose) cout << "These helices have the same contacts:" << endl;

					if(verbose) cout << "Helix " << h+1 << ":\t";
					for (unsigned int k = 0; k < h1_contact.size(); k++){
						if(verbose) cout << h1_contact[k]+1 << " ";
					}
					if(verbose) cout << endl;
					if(verbose) cout << "Helix " << j+1 << ":\t";
					for (unsigned int k = 0; k < h2_contact.size(); k++){
						if(verbose) cout << h2_contact[k]+1 << " ";
					}
					if(verbose) cout << endl;						
					
					int cross_overs = count_crossovers(total,x);
					if(verbose) cout << "Number of loop crossovers in current helix positions:\t" << cross_overs << endl;
					
					pair<double,double> swap_h1,swap_h2;
					swap_h1 = make_pair(all_helix_positions[x][h].first,all_helix_positions[x][h].second);
					swap_h2 = make_pair(all_helix_positions[x][j].first,all_helix_positions[x][j].second);
					
					if(verbose) cout << "Swapping positions of helices " << h+1 << " and " << j+1 << endl;
					
					all_helix_positions[x][h] = swap_h2;
					all_helix_positions[x][j] = swap_h1;
					
					int cross_overs_swapped = count_crossovers(total,x);
					
					if(verbose) cout << "Number of loop crossovers in new helix positions:\t" << cross_overs_swapped << endl;
					
					if(cross_overs < cross_overs_swapped){
						if(verbose) cout << "Swapped conformation has more crossovers, reverting to original helix positions." << endl;
						if(verbose) cout << "Swapping positions of helices " << j+1 << " and " << h+1 << endl;
						// Swap back
						all_helix_positions[x][h] = swap_h1;
						all_helix_positions[x][j] = swap_h2;
					
					}else if(cross_overs == cross_overs_swapped){
						if(verbose) cout << "Crossovers are equal; helices " << h+1 << " and " << j+1 << " are interchangable." << endl;
						all_helix_positions.push_back(all_helix_positions[x]);
						// Swap back
						all_helix_positions[x][h] = swap_h1;
						all_helix_positions[x][j] = swap_h2;
						
						// Recurse
						if(verbose) cout << "Recursing..." << endl;
						swap_helices(total,h1,h2,x+1);
						if(verbose) cout << "Finished recursing..." << endl;						
						helix_swaps_seen = helix_swaps_seen_copy;
					}else{
						if(verbose) cout << "Swapped conformation has fewer crossovers, using new helix positions." << endl;
					}								
				}		
			}			
			h2_contact.clear();		
		}
		h1_contact.clear();		
	}
}

int main(int argc, char* argv[]){

  	double width = 2000;
  	double height = 2000;	

       // Exit unless filename given as argument
        if (argc < 2){
                cout << "File containing graph data is required." << endl;
		exit(1);
        }
	
       	ifstream is(argv[1]);

        // Check to make sure the stream is ok
        if(!is.good()){
                cout << "Cannot open file "<< argv[1] << endl;
                exit(1);
	}else{
		cout << "Reading input file " << argv[1] << "..." << endl;
	}

        string line,topology_string;
        char ch[150],topology[150];
	int a,b,c,d, edges = 0;
	vector<int> h1, h2;
	
        // Get line from stream
        while(getline(is,line)){

                // copy the string to a character array
                strcpy(ch, line.c_str());
			
                if(sscanf(ch, "%d%*c%d %d%*c%d %*f", &a,&b,&c,&d) == 4){
			//cout << "Residues:\t" << a << " ---> " << b << " Helices:\t" << c << " ---> " << d << endl;			
			
			const pair <int,int> rr (a,b);
			const pair <int,int> hh (c,d);

			residues_to_helices.insert(pair<const pair<int,int>, const pair<int,int> >(rr,hh));			
			
			int tag = 0;
			for (unsigned int i = 0; i < h1.size(); i++){
				if (h1[i] == c-1 && h2[i] == d-1){
					tag++;
				}
			}			
			if(!tag){
				h1.push_back(c-1);
				h2.push_back(d-1);
				edges++;
			}
		
                }else if (sscanf(ch, "# Topology: %s", topology) == 1){
			topology_string = topology;
			tokenize(topology_string, boundaries);
                       	//cout << "Topology:\t" << topology << endl << endl;
		} 
        }
	is.close();

	if(!edges){
		cout << "No contacts predicted!" << endl;
		exit(1);
	}
	if(boundaries.size() % 2){
		cout << "Uneven number of helix boundaries!" << endl;
		exit(1);
	}else{
		total = boundaries.size()/2;
	}
	if(!total){
		cout << "No helices found!" << endl;
		exit(1);
	}

	// Construct graph
	Graph g;

	WeightMap weightmap = get(edge_weight, g);
      	for (int j = 0; j < edges; ++j){
	 	graph_traits<Graph>::edge_descriptor e;
		bool inserted;  
         	tie(e, inserted) = add_edge(h1[j],h2[j],g);  
		weightmap[e] =  1.2;
	}	
 
 	// connected_components
    	vector<unsigned int> component(num_vertices(g));
    	unsigned int num = connected_components(g, &component[0]);
    
    	if (num > 1) cout << endl << "Graph contains " << num << " components." << endl << endl;
    	
	//for (i = 0; i != component.size(); ++i){
      	//	cout << "Vertex " << i <<" is in component " << component[i] << endl;
    	//}
	//cout << endl;

 	for (unsigned int i = 0; i < num; i++){
	
		if (num > 1) cout << "Processing graph component " << i+1 << "..." << endl << endl; 
	
		Graph G;
		WeightMap component_weightmap = get(edge_weight, G);
		map<int,int> original_component_vertices;	
		total = 0;

    		for (unsigned int j = 0; j != component.size(); ++j){				
			if (i == component[j]){
			
				if (num > 1) cout << "Helix " << j+1 <<" is in component " << component[j]+1 << endl;
				original_component_vertices[total] = j;
				
      				for (int k = 0; k < edges; ++k){
					if ((j == (unsigned int)h1[k] && component[h2[k]] == component[j])||(j == (unsigned int)h2[k] && component[h1[k]] == component[j])){
	 					graph_traits<Graph>::edge_descriptor  e;
						bool inserted;  
         					tie(e, inserted) = add_edge(h1[k],h2[k],G);  
						component_weightmap[e] =  1.2;
						//cout << "Edge: " << h1[k]+1 << " => " << h2[k]+1 << endl;
					}
				}
				total++;			
			}			
		}	
 
  		PositionVec position_vec(num_vertices(G));  	
  		PositionMap position(position_vec.begin(), get(vertex_index, G));
		minstd_rand gen;
		double radius = 300;
		if (total >= 5) radius = 500;
		if (total >= 10) radius = 800;

		// For some reason the KK algorithm seems to get stuck when
		// we have 5 helices arranged in a circle, so randomise the
		// initial arrangement rather than arrange in a circle
		if (total == 5){			
			//cout << "Random layout:" << endl;
			//random_graph_layout(G, position, 0, (int)width, 0, (int)height, gen);
			typedef boost::rectangle_topology<minstd_rand> RectTopology;
			typedef RectTopology::point_type point;
			RectTopology rect_top(gen, 0, (int)width, 0, int(height));
			random_graph_layout(G, position, rect_top);
  		}else{	
			//cout << "Circle layout:" << endl;
			circle_graph_layout(G, position, radius);
		}

		cout << endl << "Calculating layout..." << endl << endl;
		double edge_width = 590;

		if (total >= 4){	
				
			kamada_kawai_spring_layout(G, position, component_weightmap, edge_length(edge_width)); 		
			
			if (num == 1 && (unsigned int)edges == total-1){
				cout << "Linear helix arrangment detected. Using circular layout..." << endl;				
				double x_centre = position[0].x;
				double y_centre = position[0].y - radius;
				double t = -90;	
				map<double,int> helix_x_positions;
				map<double,int>::iterator it_h1;				
				// KK will put helices in a line, sort by x coord and use to make a circle
				// without any edge crossings
				for (unsigned int c = 0; c < total; c++){
					helix_x_positions[position[c].x] = c;
				}
				for (it_h1 = helix_x_positions.begin(); it_h1 != helix_x_positions.end(); it_h1++){
					if (t > 360) t -= 360;				
					position[(*it_h1).second].x = x_centre + radius*cos(t*M_PI/180);
					position[(*it_h1).second].y = y_centre + radius*sin(t*M_PI/180);				
					t += 360.0/total;			
				
				}
			}else{
				cout << "Using Kamada-Kawai spring layout..." << endl;
			}			
		}
	
		for (unsigned int i = 0; i < total; i++){
			pair <double,double> helix_xy (position[i].x, position[i].y);		
			helix_positions[i] = helix_xy;	
		}	

		all_helix_positions.push_back(helix_positions);	
		if(total > 4){
			cout << endl << "Generating loop crossover scores for alternate arrangements..." << endl;
			swap_helices(total,h1,h2,0);
			cout << endl;
		}	
		
		map<double,int> scores;	
		double score = 0;
	
		for (unsigned int a = 0; a < all_helix_positions.size(); a++){		
			vector<int> rotations;
			for (unsigned int i = 0; i < total; i++){
				rotations.push_back(0);	
			}
			all_rotations.push_back(rotations);
		
			r_xy_map residue_positions;
			all_residue_positions.push_back(residue_positions);
				
			for (unsigned int i = 0; i < total; i++){
				get_residue_positions(all_rotations[a][i], i, a);
			}		
			cout << "Optimising helix rotation..." << endl << endl;
			score = optimise_parameters(total,a)/h1.size();
			scores[score] = a;
			cout << endl;	
		}

		int count = 1;
		map<double,int>::iterator it_h1;
		for (it_h1 = scores.begin(); it_h1 != scores.end(); it_h1++){
			//cout << "----" << (*it_h1).first << "--->" << (*it_h1).second << endl;
			//cout << "Arrangement " << count << " (" << (*it_h1).second << ")" << endl;
			if (scores.size() > 1) cout << "Arrangement " << count << endl;
			cout << "Helix\tPosition\t\tRotation" << endl;
			for (unsigned int i = 0; i < total; i++){
				cout << original_component_vertices[i]+1 << "\t(" << all_helix_positions[(*it_h1).second][i].first << "," << all_helix_positions[(*it_h1).second][i].second << ")\t" << all_rotations[(*it_h1).second][i] << endl;
			}
			cout << "Score:\t" << (*it_h1).first << endl;
       			cout << "========================================" << endl;
			cout << endl;	
			count++;	
		}
 

		helix_positions.clear();
		all_helix_positions.clear();
		helix_swaps_seen.clear();
		all_rotations.clear();
		all_residue_positions.clear();	
		original_component_vertices.clear();
	}
	
  	return 1;
}

