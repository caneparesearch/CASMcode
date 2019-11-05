#include<iostream>
#include<fstream>
#include<istream>
#include<cstdlib>
#include<cstring>
#include<cmath>
#include<vector>
#include<stdlib.h>
#include<stdio.h>
#include<limits.h>
#include<fnmatch.h>
#include<sys/types.h>
#include<unistd.h>
#include<dirent.h>
#include<iomanip.h>
#include<time.h>
#include "clusters11.0.h"
#include "Array.h"

using namespace std;

//NOTE: This is licensed software. See LICENSE.txt for license terms.
//written by Anton Van der Ven, John Thomas, Qingchuan Xu, and Jishnu Bhattacharya
//please see CASMdocumentation.pdf for a tutorial for using the code. 
//Version as of May 26, 2010




int main(){

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //READ THE PRIM FILE
	
  structure prim;
  	
  ifstream in;
  if(scandirectory(".","prim")) in.open("prim");
  else if(scandirectory(".","PRIM")) in.open("PRIM");
  else cout << "No PRIM file in the current directory \n";	
	
  if(!in){
    cout << "cannot open file.\n";
    return 0;
  }
  prim.read_struc_prim(in);
  in.close();
  cout << " finished read_prim \n";
  prim.read_species();   // jishnu
  prim.collect_components();   // jishnu
  prim.get_trans_mat();
  prim.update_struc();
  prim.write_point_group();
  prim.write_factor_group();
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  cout << " Checkpoint 1: Crystallography routines complete. \n";

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //CLUSTER GENERATION PART

  int min_num_compon=2;
  int max_num_points;
  vector<double> max_radius;
  multiplet clustiplet;
  multiplet basiplet;

  read_cspecs(max_radius);
  max_num_points=max_radius.size()-1;
  generate_ext_clust(prim,min_num_compon,max_num_points,max_radius,clustiplet);
  generate_ext_basis(prim,clustiplet,basiplet);

  write_clust(clustiplet,"CLUST");

  write_fclust(clustiplet,"FCLUST");

  write_clust(basiplet,"BCLUST");

  write_fclust(basiplet,"FBCLUST");

  basiplet.get_hierarchy();
  //-----------------------------	
  int eci_in_choice =1;
  if(scandirectory(".","eci.in")){	
    cout << "generate (0) or keep (1) eci.in\n";
    cin >> eci_in_choice;
  } 
  if(eci_in_choice == 0) {
    string command="rm eci.in ";
    int s=system(command.c_str());
    if(s == -1){cout << "was unable to perform system command\n";}
  }
  if(!(eci_in_choice == 0) && !(eci_in_choice == 1)) {
    cout << "Invalid choice of eci.in generation option : exiting \n";
    exit(1);
  }
		

  if(!scandirectory(".","eci.in")){
    string eci_in_file="eci.in";
    ofstream out;
    out.open(eci_in_file.c_str());
    if(!out){
      cout << "cannot open file.\n";
      return 0;
    }
    basiplet.print_hierarchy(out);
    out.close();
  }
  //-----------------------------	
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  cout << " Checkpoint 2:  Cluster routines complete. \n";

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //SUPERCELL GENERATION PART
	
  vector<structure> suplat;
  int gr=-1;
  if(scandirectory(".","SCEL")){
    if((gr > 1) || (gr < 0)){
      cout << " generate (0) or read (1) supercells\n";  
      cin >> gr;    	
    }
  }

  if(gr==0 || !scandirectory(".","SCEL")){
    //generate the supercells
    int sup_vol=1;
    cout << "enter the maximum supercell volume\n";  
    cin >> sup_vol;
		
    if(sup_vol <= 0) {  // jishnu
      cout << "Invalid choice of supercell volume: exiting \n";
      exit(1);
    }
		
    int dimension=1;
    while(dimension >3 || dimension < 2){
      cout << "generate 3-dimensional or 2-dimensional supercells\n";  
      cout << "(enter the number 3 or 2)\n";
      cin >> dimension;
      if((dimension > 3)||(dimension < 2)){  // jishnu
	cout << "Invalid choice of supercell dimension: exiting \n";
	exit(1);
      }
			
    }
    if(dimension == 3){
      prim.generate_3d_supercells(suplat,sup_vol);
      for(int sc=0; sc<suplat.size(); sc++){
	suplat[sc].generate_3d_reduced_cell();
	suplat[sc].get_latparam();
      }
    }
    if(dimension == 2){
      int excluded_axis=0;
      while(excluded_axis > 3 || excluded_axis < 1){ 
	cout << "enter the excluded axis \n";
	cout << "(enter a number between 1 and 3)\n";
	cin >> excluded_axis;
      }
      excluded_axis=excluded_axis-1;
      prim.generate_2d_supercells(suplat,sup_vol,excluded_axis);
      for(int sc=0; sc<suplat.size(); sc++){
	suplat[sc].generate_2d_reduced_cell(excluded_axis);
	suplat[sc].get_latparam();
      }
    }
    write_scel(suplat);
  }
  else if(gr==1){
    //read in the supercells
    read_scel(suplat,prim);
  }
  else {  // jishnu
    cout << "Invalid choice of supercell generation option: exiting \n";
    exit(1);
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cout << " Checkpoint 3:  Supercell routines complete. \n";


  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //CONFIGURATION GENERATION PART

  configurations configs;
  configs.prim=prim;
  for(int np=0; np<basiplet.orb.size(); np++){
    configs.basiplet.orb.push_back(basiplet.orb[np]);
  }
  int conf_choice = 0;
  if(scandirectory(".","configuration")){	
    cout << "generate (0) or read (1) configuration and configuration.corr\n";
    cin >> conf_choice;
  }
  if(conf_choice == 0){	  
    configs.generate_configurations_fast(suplat);
    configs.print_con();
    configs.print_corr();
    // configs.print_con_old();  // old format for fitting code
    // configs.print_corr_old(); // old format for fitting code
  }  
  else if (conf_choice == 1) {    
    configs.reconstruct_from_read_files();
    //configs.print_con_old();
    //configs.print_corr_old();
  }
  else {
    cout << "Invalid choice of configuration generation option : exiting \n";
    exit(1);
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  cout << " Checkpoint 4:  Configuration routines complete. \n";

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //CREATE VASP FILES
  int make_dir_choice =1;
  if(scandirectory(".","make_dirs")){	
    cout << "generate (0) or read (1) make_dirs\n";
    cin >> make_dir_choice;
  } 
  if(make_dir_choice == 0) {
    string command="rm make_dirs ";
    int s=system(command.c_str());
    if(s == -1){cout << "was unable to perform system command\n";}
  }
  if(!(make_dir_choice == 0) && !(make_dir_choice == 1)) {
    cout << "Invalid choice of make_dirs generation option : exiting \n";
    exit(1);
  }
		
		
  configs.print_make_dirs();
  configs.read_make_dirs();
  configs.print_make_dirs();
  configs.generate_vasp_input_directories();
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cout << " Checkpoint 5:  VASP input routines complete. \n"; 

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //COLLECT VASP ENERGIES
  //CALCULATE FORMATION ENERGIES AND PRINT INFO INTO ener.in and corr.in
		
  int energy_choice = 0;
  if(scandirectory(".","energy")){
    cout << "generate (0) or read (1) energy and corr.in \n";
    cin >> energy_choice;
  }
  if(energy_choice == 0){
    configs.collect_energies_fast();
    // collect refernces                                                                                                                                                                                                                                                           
    int reference_choice =1;
    if(scandirectory(".","reference")){
      cout << "generate (0) or keep (1) reference\n";
      cin >> reference_choice;
    }
    if(reference_choice == 0) {
      string command="rm reference ";
      int s=system(command.c_str());
      if(s == -1){cout << "was unable to perform system command\n";}
    }
    if(!(reference_choice == 0) && !(reference_choice == 1)) {
      cout << "Invalid choice of reference generation option : exiting \n";
      exit(1);
    }

	cout << "The size of the configs vector is " << configs.superstruc.size()<< "\n";
	
//	for (int i=0; i<configs.superstruc.size();i++) {
//	  for (int j=0; j<configs.superstruc[i].conf.size();j++) {
	
//	  cout << "The name of the arrangement indexed by " << i << " " << j <<  " is " << configs.superstruc[i].conf[j].name << " \n";
//	  cout << "The energy of the arrangement is " << configs.superstruc[i].conf[j].energy << " \n";
//	  }
//	}
	
	
	cout << "before collect reference \n";
    configs.collect_reference();
    //---------------------------
	
	cout << "before calculate_formation_energy \n";	
    configs.calculate_formation_energy();
	
//	cout << "before assemble_hull \n";	
//    configs.assemble_hull();
	
	cout << "before print_eci \n";	
    configs.print_eci_inputfiles();
  }
  else if (energy_choice == 1) {
    configs.read_energy_and_corr();
    configs.assemble_hull();
  }
  else {
    cout << "Invalid choice of energy and corr.in files generation option : exiting \n";
    exit(1);
  }

  configs.print_make_dirs();
  cout << " Checkpoint 6:  VASP output interpretation routines completed. \n";  
  if(configs.chull.face.size()){
    configs.get_delE_from_hull();
    configs.chull.write_hull();
  }
  if(!scandirectory(".","eci.out")){
    cout << "No cluster expansion is present.  If eci.out file has been obtained, copy to current directory and rerun code.\n";
    cout << "Otherwise, CASM code has completed successfully.  Exiting...\n";
    return 0;
  }
  else if(configs.chull.face.size()) configs.CEfenergy_analysis();    
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cout << " Checkpoint 7:  First principles and cluster expansion energy analysis completed.  \n \n";
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
  cout << " Proceeding with MONTE CARLO tests... \n";

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //MONTE CARLO TESTS

  //ask if Monte Carlo should be made

  //read the eci and assign them to the correct clusters in basiplet

  string eci_file="eci.out";
  ifstream eci;
  eci.open(eci_file.c_str());
  if(!in){
    cout << "cannot open file.\n";
    return 0;
  }

  basiplet.read_eci(eci);

  structure groundstate=prim;
  //  string struc_file="POSCAR";


  //  in.open(struc_file.c_str());
  //  if(!in){
  //    cout << "cannot open " << struc_file << "\n";
  //    return 0;
  //  }
  //  groundstate.read_struc_poscar(in);
  //  in.close();

  string class_file="clusters11.0.h";

  //input size of the Monte Carlo cell

  Monte_Carlo monte(prim, groundstate, basiplet, 12, 12, 12);
  monte.write_monte_h(class_file);   
  cout << "Clusters code has executed successfully.  Exiting...\n";
  return 0;
}

