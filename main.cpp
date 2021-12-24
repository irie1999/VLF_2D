/*
 * 2次元FDTD法 with electron density perturbation
 *
 * Usage:
 * main Lp[km] z[km] sigma[km]
 *
 * Lp[km]: The center of the perturbation measured from the current source
 * z_dec[km]: decreasing height of the electron density profile
 * sigma[km]: the standard deviation of the perturbation
 *
 */
#include <bits/stdc++.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>

#include <eigen3/Eigen/Dense>

#include <memory_allocate.h>
#include "fdtd2d.h"

std::string data_dir = "data/";


int main(int argc, char **argv){
  double h_prime = 85.0e3; /* [m] effective reflection height */
  double beta = 0.63; /* sharpness factor */

  double Lp = 0.0; /* [m] center_perturbation */
  double z_dec = 0.0; /* [m] height_decrease */
  double sig_per = 1.0; /* [m] sigma_perturbation*/
  //double *Ne = new double [80];
  double **Ei_tm = allocate_memory2d(3,  Nth, 0.0); /*前の時刻との電界強度の変化率*/  
  
  //input(Ne, t, s);
 for(int t = 0; t < 3; t++){ /*t0=6:10 t1=6:20 t2=6:30*/
    if(t == 0){  
    h_prime = 77.69128e3;
    beta = 0.493661;
    }

    if(t == 1){
      h_prime = 76.9801e3;
      beta = 0.511845;
    }

    if(t == 2){
      h_prime = 76.0717e3;
      beta = 0.53594;
    }


    double ***Dr  = allocate_memory3d(2, Nr,   Nth+1, 0.0);
    double ***Dth = allocate_memory3d(2, Nr+1, Nth, 0.0);
    double ***Dph = allocate_memory3d(2, Nr+1, Nth+1, 0.0);
    double ***Er  = allocate_memory3d(2, Nr,   Nth+1, 0.0);;
    double ***Eth = allocate_memory3d(2, Nr+1, Nth, 0.0);
    double ***Eph = allocate_memory3d(2, Nr+1, Nth+1, 0.0);
    double **Hr   = allocate_memory2d(Nr+1, Nth, 0.0);
    double **Hth  = allocate_memory2d(Nr,   Nth+1, 0.0);
    double **Hph  = allocate_memory2d(Nr,   Nth, 0.0);
     
    
  
    /* Coefficient matrix to update E taking account into the ionosphere */
    Eigen::Matrix3d **C = new Eigen::Matrix3d* [Nr_iono];
    Eigen::Matrix3d **F = new Eigen::Matrix3d* [Nr_iono];
    Eigen::Matrix3d *C1 = new Eigen::Matrix3d [Nr_iono*(Nth+1)];
    Eigen::Matrix3d *F1 = new Eigen::Matrix3d [Nr_iono*(Nth+1)];

    for(int i = 0; i < Nr_iono; i++){
      C[i] = C1 + i*Nth;
      F[i] = F1 + i*Nth;
      for(int j = 0; j <= Nth; j++){
        C[i][j] << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
        F[i][j] << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
      }
    }
  
    /*master*/
    initialize_conductivity(C, F, h_prime, beta, Lp, z_dec, sig_per);
    /*IRI*/
    //initialize_conductivity(C, F, h_prime, beta, Lp, z_dec, sig_per, Ne);
    /*betaとh*/
    //initialize_conductivity(C, F, h_prime, beta, Lp, z_dec, sig_per, t);
    
    /* PML only for +Theta-directed layer */
    double **Dr1 =    allocate_memory2d(Nr,   PML_L+1, 0.0);
    double **Dr2 =    allocate_memory2d(Nr,   PML_L+1, 0.0);
    double **Dph_r =  allocate_memory2d(Nr+1, PML_L+1, 0.0);
    double **Dph_th = allocate_memory2d(Nr+1, PML_L+1, 0.0);
    double **Hr1 =    allocate_memory2d(Nr+1, PML_L,   0.0);
    double **Hr2 =    allocate_memory2d(Nr+1, PML_L,   0.0);
    double **Hph_r  = allocate_memory2d(Nr,   PML_L,   0.0);
    double **Hph_th = allocate_memory2d(Nr,   PML_L,   0.0);
    
    double **Bph = allocate_memory2d(2,   PML_L,   0.0);
    double *Bph_r = new double [PML_L];
    double *Bph_th = new double [PML_L];
    for(int i = 0; i < PML_L; i++){
      Bph_r[i] = 0.0;
      Bph_th[i] = 0.0;
    }

    double *C01 = new double [PML_L+1];
    double *C02 = new double [PML_L+1];
    double *C11 = new double [PML_L];
    double *C12 = new double [PML_L];
    initialize_pml(C01, C02, C11, C12);

    /* Vertical E-field at Earth's surface (f kHz) */
    std::complex <double> zj { 0., 1. };
    std::complex <double> *Er0 = new std::complex <double> [Nth+1 - PML_L];
    for(int i = 0; i <= Nth - PML_L; i++){
      Er0[i] = std::complex <double> {0., 0.};
    }

    /* Surface impedance */
    double *Rs = new double [Nth + 1];
    double *Ls = new double [Nth + 1];
    initialize_surface_impedance(Rs, Ls);

    ///時間ループ///
    for(int n = 1; n <= Nt; n++){
      if ( n%100 == 0 ){
        std::cout << " " << n << " / " << Nt << " : " << Lp/1000 << "\n";
      }
      int NEW = (n+1) % 2;
      int OLD = n % 2;

      update_Dr(Dr, Hph, NEW, OLD);
      update_Dth(Dth, Hph, NEW, OLD);
      update_Dph(Dph, Hr, Hth, NEW, OLD);

      update_Dr_PML(Dr[NEW], Dr1, Dr2, Hph, C01, C02);
      update_Dph_PML(Dph[NEW], Dph_r, Dph_th, Hth, Hr, C01, C02);

      double t = (n - 0.5) * Dt;
      Dr[NEW][0][0] -= Dt * Jr(t); //θ=0, i=j=0  //Jr

      update_Er(Er, Eth, Eph, Dr, Dth, Dph, NEW, OLD, C, F);
      update_Eth(Er, Eth, Eph, Dr, Dth, Dph, NEW, OLD, C, F);
      update_Eph(Er, Eth, Eph, Dr, Dth, Dph, NEW, OLD, C, F);

      update_Hr(Hr, Eph[NEW]);
      update_Hth(Hth, Eph[NEW], Rs, Ls);
      update_Hph(Hph, Er[NEW], Eth[NEW], Rs, Ls);

      update_Hr_PML(Hr, Hr1, Hr2, Eph[NEW], C11, C12);
      update_Hph_PML(Hph, Hph_r, Hph_th, Er[NEW], Eth[NEW], C11, C12,
          Rs, Ls, Bph, Bph_r, Bph_th, NEW);

     //    if ( n%10 == 0 ) output(Er, NEW, n);

      t = n * Dt;
      for(int j = 0; j <= Nth - PML_L; j++){
        Er0[j] += Er[NEW][0][j] * std::exp( -1.0 * zj * OMG * t ) * Dt;
      }
    }

    for(int j = 0; j <= 2 * Nr_GA; j += 2){  /*Ei_tm[][0] = 100kmの電界　Ei_tm[][800] = 900kmの電界*/
      Ei_tm[t][j/2] = 20 * std::log10(std::abs( Er0[j + 2 * 100] )); /*Er0[セル数]*/
    }
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "Ei_tm[0][0]= " << Ei_tm[0][0] << std::endl
              << "Ei_tm[1][0]= " << Ei_tm[1][0] << std::endl;

    deallocate_memory3d(Dr);
    deallocate_memory3d(Dth);
    deallocate_memory3d(Dph);
    deallocate_memory3d(Er);
    deallocate_memory3d(Eth);
    deallocate_memory3d(Eph);
    deallocate_memory2d(Hr);
    deallocate_memory2d(Hth);
    deallocate_memory2d(Hph);
    delete [] C1;
    delete [] C;
    delete [] F1;
    delete [] F;
    delete [] Er0;
 }
  double time[3];
  time[0] = 6.16667, time[1] = 6.33333, time[2] = 6.5;
  double **Si_tm = allocate_memory2d(3,  Nth, 0.0);
  for(int t = 1; t < 3; t++){
    for(int j = 0; j <= Nr_GA; j++){
      Si_tm[t][j] = (Ei_tm[t][j] - Ei_tm[t-1][j]) / (time[t] - time[t-1]);
    }
  }

  /*100kmから900kmまで観測点*/
  std::ofstream ofs("../data/" "Si_tm_1.dat");
  for(int j = 0; j <= Nr_GA; j++){
    ofs << j << " " << Si_tm[1][j] << std::endl;
  }
  ofs.close();

  std::ofstream ofs1("../data/" "Si_tm_2.dat");
  for(int j = 0; j <= Nr_GA; j++){
    ofs1 << j << " " << Si_tm[2][j] << std::endl;
  }
  ofs1.close();
  

    
  return 0;
}
