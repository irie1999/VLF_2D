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

#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>

#include <eigen3/Eigen/Dense>

#include <memory_allocate.h>
#include "fdtd2d.h"

std::string data_dir = "data/";

int main(int argc, char **argv){

  double z_prime = 85.0e3; /* [m] effective reflection height */
  double beta = 0.63; /* sharpness factor */

  double Lp = 0.0; /* [m] center_perturbation */
  double z_dec = 0.0; /* [m] height_decrease */
  double sig_per = 1.0; /* [m] sigma_perturbation*/

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
  
  initialize_conductivity(C, F, z_prime, beta, Lp, z_dec, sig_per);
  
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
      std::cout << n << " / " << Nt << " : " << Lp/1000 << "\n";
    }
    int NEW = (n+1) % 2;
    int OLD = n % 2;
    
   std::cout<<n<<" "<<  Dr[NEW][0][0]   <<std::endl;

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

  std::ofstream ofs(data_dir + "Er_surface" + suffix(Lp, z_dec, sig_per) + ".dat");
  for(int j = 0; j <= Nth - PML_L; j++){
    ofs << j * R0 * dth * 1e-3 << " " << std::abs( Er0[j] ) << " "
        << std::arg( Er0[j] ) << "\n";
  }
  ofs.close();

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
    
  return 0;
}
