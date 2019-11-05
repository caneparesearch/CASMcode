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
#include "monte.h"
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
  string prim_file="PRIM";

  ifstream in;
  in.open(prim_file.c_str());
  if(!in){
    cout << "cannot open file.\n";
    return 0;
  }
  prim.read_struc_poscar(in);
  in.close();
  cout << " finished read_poscar \n";

  prim.get_trans_mat();
  prim.update_struc();
  prim.write_point_group();
  prim.write_factor_group();
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



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
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //MONTE CARLO TESTS


  //read the eci and assign them to the correct clusters in basiplet

  string eci_file="eci.out";
  ifstream eci;
  eci.open(eci_file.c_str());
  if(!in){
    cout << "cannot open file" << eci_file << ".\n";
    return 0;
  }  

  basiplet.read_eci(eci);

  structure groundstate;
  string struc_file="POSCAR";


  in.open(struc_file.c_str());
  if(!in){
    cout << "cannot open " << struc_file << "\n";
    return 0;
  }
  groundstate.read_struc_poscar(in);
  in.close();  

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //read in the conditions file

  concentration tconc;
  tconc.collect_components(prim);

  int n_pass,n_equil_pass, nx, ny, nz, output_poscar_step, corr_flag, temp_chem;
  double Tinit,Tmin,Tmax,Tinc;
  double Temp;
  chempot muinit,mu_min,mu_max;
  vector<chempot> muinc;  
  chempot mu;
  mu.initialize(tconc);
  cout << "JUST INITIALIZED mu \n";
  cout << "THE COMPONENTS OF MU ARE \n";
  mu.print_compon(cout);
  cout << "\n";
  muinit.initialize(tconc);
  //  muinc.initialize(tconc); //edited
  mu_min.initialize(tconc);
  mu_max.initialize(tconc);

  string cond_file="mc_input";
  if(!read_mc_input(cond_file, n_pass, n_equil_pass, nx, ny, nz, Tinit, Tmin, Tmax, Tinc, muinit, mu_min, mu_max, muinc, output_poscar_step, corr_flag, temp_chem)){
    exit(1);
  }
 
  cout << "Printing mu_min \n";
  for(int ii=0; ii<mu.compon.size(); ii++){
    for(int jj=0; jj<mu.compon[ii].size(); jj++){
      cout << mu_min.compon[ii][jj].name << ' ' << mu_min.m[ii][jj] << "\t";
    }
  }
  cout << "\n";
  cout << "Printing mu_max \n";
  for(int ii=0; ii<mu.compon.size(); ii++){
    for(int jj=0; jj<mu.compon[ii].size(); jj++){
      cout << mu_max.compon[ii][jj].name << ' ' << mu_max.m[ii][jj] << "\t";
    }
  }
  cout << "\n";

  cout << "Printing muinit \n";
  for(int ii=0; ii<mu.compon.size(); ii++){
    for(int jj=0; jj<mu.compon[ii].size(); jj++){
      cout << muinit.compon[ii][jj].name << ' ' << muinit.m[ii][jj] << "\t";
    }
  }
  cout << "\n";

  cout << "Printing muinc \n";
  for(int kk=0; kk<muinc.size(); kk++){
    cout << "muinc[" << kk << "]\n";
    for(int ii=0; ii<mu.compon.size(); ii++){
      for(int jj=0; jj<mu.compon[ii].size(); jj++){
	cout << muinc[kk].compon[ii][jj].name << ' ' << muinc[kk].m[ii][jj] << "\t";
      }
    }
    cout << "\n";
  }


  Monte_Carlo monte(prim, groundstate, basiplet, nx, ny, nz);
  monte.corr_flag=corr_flag;

  for(int i=0; i<monte.montiplet.size(); i++){
    string track;
    string mfile="MCLUST";
    int_to_string(i, track, 10);    
    mfile.append(track);
    write_fclust(monte.montiplet[i],mfile);
  }

  double energy;
  monte.calc_energy(energy);

  cout << " The energy of the cell is = " << energy/monte.nuc << "\n";

  if(temp_chem==0){

    mu=muinit;
    bool continue_mu_flag=true;  //edited
    int inc_count=0;
    while(continue_mu_flag){     //edited
      cout << "Working on chemical potential: ";
      mu.print(cout);
      Temp=Tinit;

      string temp_file="temp";
      string track;
      for(int s_m=0; s_m<muinc.size(); s_m++){
	for(int i_m=0; i_m<muinc[s_m].m.size(); i_m++){
	  for(int j_m=0; j_m<muinc[s_m].m[i_m].size()-1; j_m++){
	    if(abs(muinc[s_m].m[i_m][j_m])>1e-6){
	      string ttrack;
	      track.append(".");
	      int_to_string(int(round(mu.m[i_m][j_m]*1000)), ttrack, 10);
	      track.append(ttrack);
	    }
	  }
	}
      }
      temp_file.append(track);
      ofstream out(temp_file.c_str());

      if(!out){
        cout << "cannot open " << temp_file << "\n";
        return 0;
      }

      string corr_file="corr.";
      corr_file.append(track);
      ofstream corrstream;
      if(corr_flag){
	corrstream.open(corr_file.c_str());
	if(!corrstream){
	  cout << "cannot open " << corr_file << "\n";
	  return 0;
	}
      }
      out << "# chemical potential = ";
      mu.print(out);
      out << "\n";
      out << "# Temp  AVconc  AVenergy  heatcapacity  ";
      monte.AVSusc.print_elements(out);
      out << "  flipfreq \n";
      monte.initialize(groundstate);

      while(Temp >= Tmin && Temp <= Tmax){

        double beta=1.0/(kb*Temp);
        monte.grand_canonical(beta,mu,n_pass,n_equil_pass);

	//Print structure if dictated by mc_input
	if(output_poscar_step>0 && !(inc_count%output_poscar_step)){
	  string outstruc_file="struc.";
	  string track2;
	  int_to_string(int(Temp),track2,10);
	  outstruc_file.append(track2);
	  for(int s_m=0; s_m<muinc.size(); s_m++){
	    for(int i_m=0; i_m<muinc[s_m].m.size(); i_m++){
	      for(int j_m=0; j_m<muinc[s_m].m[i_m].size()-1; j_m++){
		if(abs(muinc[s_m].m[i_m][j_m])>tol){
		  track2.clear();
		  outstruc_file.append(".");
		  int_to_string(int(mu.m[i_m][j_m]*1000), track2,10);
		  outstruc_file.append(track2);
		}
	      }
	    }
	  }
	  outstruc_file.append(".POSCAR");
	  ofstream struc_stream(outstruc_file.c_str());
	  monte.write_monte_poscar(struc_stream);
	  struc_stream.close();
	}
	inc_count++;
	//\Edited code

        //Edited Code -- Print  correlations.  
	if(corr_flag){
	  for(int ii=0; ii<monte.AVcorr.size(); ii++)
	    corrstream << monte.AVcorr[ii] << "  ";
	  corrstream << "\n";
	}
	//\Edited Code 

        out << Temp << "  ";
        monte.AVconc.print_concentration_without_names(out);
        out << monte.AVenergy/monte.nuc << "  ";
        out << monte.heatcap << "  ";
	monte.AVSusc.print(out);
        out << "  " << monte.flipfreq << "\n";

        Temp=Temp+Tinc;

      }
      out.close();

      //Edited Code
      for(int ii=0; ii<muinc.size(); ii++){
	bool break_flag=false;
	for(int jj=0; jj<muinc[ii].m.size(); jj++){
	  for(int kk=0; kk<muinc[ii].m[jj].size(); kk++){
	    
	    if(abs(muinc[ii].m[jj][kk])>0.000001 && (mu.m[jj][kk]+muinc[ii].m[jj][kk]>=mu_min.m[jj][kk]) && (mu.m[jj][kk]+muinc[ii].m[jj][kk]<=mu_max.m[jj][kk])){
	      mu.increment(muinc[ii]);
	      break_flag=true;
	    }
	    else{
	      mu.m[jj][kk]=muinit.m[jj][kk];
	    }
	    if(break_flag)break;
	  }
	  if(break_flag)break;
	}
	if(break_flag)break;
	else if(ii==(muinc.size()-1)) continue_mu_flag=false;
      }

      //Edited Code

    }
  }
  else{
    Temp=Tinit;

    while(Temp >= Tmin && Temp <= Tmax){
      double beta=1.0/(kb*Temp);

      string chem_file="chem.";
      string track;
      int_to_string(Temp, track, 10);
      chem_file.append(track);
      ofstream out(chem_file.c_str());
      if(!out){
        cout << "cannot open " << chem_file << "\n";
        return 0;
      }

      string corr_file="corr.";
      corr_file.append(track);
      ofstream corrstream;
      if(corr_flag){
	corrstream.open(corr_file.c_str());
	if(!corrstream){
	  cout << "cannot open " << corr_file << "\n";
	  return 0;
	}
      }


      out << "# chemical potential = " << Temp << "\n";
      out << "# mu of ";
      mu.print_compon(out); 
      out << " AVconc  AVenergy  heatcapacity  ";
      monte.AVSusc.print_elements(out);
      out << "  flipfreq \n";

      cout << " about to assign mu = muinit \n";
      mu=muinit;

      cout << "just assigned mu = muinit \n";
      bool continue_mu_flag=true;   //edited
      int inc_count=0;
      monte.initialize(groundstate);
      while(continue_mu_flag){  //edited
        cout << "Working on mu = ";
        mu.print(cout);
        cout << "\n";


        monte.grand_canonical(beta,mu,n_pass,n_equil_pass);

	//Edited code -- Print structure
	if(output_poscar_step>0 && !(inc_count%output_poscar_step)){
	  string outstruc_file="struc.";
	  outstruc_file.append(track);
	  for(int s_m=0; s_m<muinc.size(); s_m++){
	    for(int i_m=0; i_m<muinc[s_m].m.size(); i_m++){
	      for(int j_m=0; j_m<muinc[s_m].m[i_m].size()-1; j_m++){
		if(abs(muinc[s_m].m[i_m][j_m])>tol){
		  string track2;
		  outstruc_file.append(".");
		  int_to_string(int(mu.m[i_m][j_m]*1000), track2,10);
		  outstruc_file.append(track2);
		}
	      }
	    }
	  }
	  outstruc_file.append(".POSCAR");
	  ofstream struc_stream(outstruc_file.c_str());
	  monte.write_monte_poscar(struc_stream);
	  struc_stream.close();
	}
	inc_count++;
	//\Edited code


        //Edited Code -- Print  correlations.  
	if(corr_flag){
	  for(int ii=0; ii<monte.AVcorr.size(); ii++)
	    corrstream << monte.AVcorr[ii] << "  ";
	  corrstream << "\n";
	}
	//\Edited Code 


        mu.print(out);
        out << "  ";
        monte.AVconc.print_concentration_without_names(out);
        out << monte.AVenergy/monte.nuc << "  ";
        out << monte.heatcap << "  ";
	monte.AVSusc.print(out);
	out << "  " << monte.flipfreq << "\n";

	//Edited Code
	bool break_flag=false;
	for(int ii=0; ii<muinc.size(); ii++){
	  for(int jj=0; jj<muinc[ii].m.size(); jj++){
	    for(int kk=0; kk<muinc[ii].m[jj].size(); kk++){
	      if(abs(muinc[ii].m[jj][kk])>0.000001 && (mu.m[jj][kk]+muinc[ii].m[jj][kk]>=mu_min.m[jj][kk]) && (mu.m[jj][kk]+muinc[ii].m[jj][kk]<=mu_max.m[jj][kk])){
		mu.increment(muinc[ii]);
		break_flag=true;
	      }
	      else if(muinc[ii].m[jj][kk]<-0.000001||muinc[ii].m[jj][kk]>0.000001){
		mu.m[jj][kk]=muinit.m[jj][kk];
	      }
	      if(break_flag){
		break;
	      }
	    }
	    if(break_flag){
	      break;
	    }
	  }
	  if(break_flag){
	    break;
	  }
	  else if(ii==(muinc.size()-1)){
	    continue_mu_flag=false;
	  }
	}
	
	//Edited Code
      }
      Temp=Temp+Tinc;
      out.close();
    }

  }




  return 0;
}

