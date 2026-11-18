#ifndef SBN_CODE_EVENTGENERATOR_MEVPRTL_TOOLS_HIGGS_KMPIONFORMFACTOR_H
#define SBN_CODE_EVENTGENERATOR_MEVPRTL_TOOLS_HIGGS_KMPIONFORMFACTOR_H

#include <array>
#include <cmath>
#include <complex>

namespace evgen {
namespace ldm {
namespace KM {

using Complex = std::complex<double>;

constexpr double Fpi = 0.093;

constexpr double l3  =  2.9;
constexpr double l4  =  4.3;
constexpr double l7  =  1.0;
constexpr double l11 =  8.7;
constexpr double l12 = -0.4;
constexpr double l13 =  9.4;

constexpr double delta = -0.3;

inline double LoopFactor(){
  return 1.0 / std::pow(4.0 * M_PI * Fpi, 2);
}

inline Complex K(double s, double m){

  const double threshold = 4.0 * m * m;
  if (s < threshold) {
    const double rho = std::sqrt(threshold / s - 1.0);
    return {-2.0 + 2.0 * rho * std::atan(1.0 / rho), 0.0};
    }
  if (s == threshold) return {-2.0, 0.0};
  const double sigma = std::sqrt(1.0 - threshold / s);

  return {-2.0 + sigma * std::log((1.0 + sigma) / (1.0 - sigma)), -M_PI * sigma};
}

inline Complex GammaPiPM(double s, double mpi0, double mpipm){
  const double m02 = mpi0 * mpi0;
  const double logMass = std::log(mpipm / mpi0);
  return 1.0 + LoopFactor()*(
         s*(l4-1.0-logMass) + 0.5*(m02-s)*K(s, mpi0) - 0.5*s*K(s, mpipm));
}

inline Complex GammaPi0(double s, double mpi0, double mpipm){
  const double m02 = mpi0 * mpi0;
  const double logMass = std::log(mpipm / mpi0);
  return 1.0 + LoopFactor()*(
         s*(l4-1.0-2.0*logMass) + (m02-s)*K(s, mpipm) - 0.5*m02*K(s, mpi0));
}

inline double GammaTildePiPM0(double mpi0, double mpipm){
  const double m02 = mpi0*mpi0;
  const double mpm2 = mpipm*mpipm;
  return mpm2 + m02*m02*LoopFactor()*0.5*(1.0 - l3) - (mpm2 - m02);
}

inline double GammaTildePi00(double mpi0, double mpipm){
  const double m02 = mpi0 * mpi0;
  const double logMass = std::log(mpipm / mpi0);

  const double correction = m02*LoopFactor()*
                            (0.5*(1.0 - l3) - 2.0*l7*delta*delta + 2.0*logMass);

  return m02*(1.0 + correction);
}

inline Complex ThetaPiPM(double s, double mpi0, double mpipm){
  const double mc2 = mpipm * mpipm;
  const double m02 = mpi0 * mpi0;
  const double logMass = std::log(mpipm / mpi0);

  return 2.0*mc2 + s + 
         LoopFactor()* ((2.0 / 3.0)*(m02 - 2.0*s)*s 
                     + s*(3.0*s - 2.0*m02)*logMass
                     - 0.5*(s - m02)*(s + 2.0*m02)*K(s, mpi0)
                     - 0.5*s*(s + 2.0*mc2)*K(s, mpipm)
                     + s*(s - 2.0*m02)*l11
                     + 4.0*s*(s - m02)*l12
                     + 1.5*s*m02*l13);
}

inline Complex ThetaPi0(double s, double mpi0, double mpipm)
{
  const double m02 = mpi0 * mpi0;
  const double logMass = std::log(mpipm / mpi0);

  return 2.0*m02 + s + 
         LoopFactor()* ((2.0 / 3.0)*(m02 - 2.0*s)*s
                     - 2.0*s*(s - m02)*logMass
                     - 0.5*m02*(s + 2.0*m02)*K(s, mpi0)
                     - (s - m02)*(s + 2.0*m02)*K(s, mpipm)
                     + s*(s - 2.0*m02)*l11
                     + 4.0*s*(s - m02)*l12
                     + 1.5*s*m02*l13);
}

inline std::array<double, 2> PionFormFactors(double higgsMass, double mpi0, double mpipm){
  const double s = higgsMass * higgsMass;
  const Complex F_pi0 = (2.0/9.0)*ThetaPi0(s, mpi0, mpipm) +
                        (7.0/9.0)*GammaTildePi00(mpi0, mpipm)*GammaPi0(s, mpi0, mpipm);
          
  const Complex F_pipm = (2.0/9.0)*ThetaPiPM(s, mpi0, mpipm) +
                         (7.0/9.0)*GammaTildePiPM0(mpi0, mpipm)*GammaPiPM(s, mpi0, mpipm);

  return {std::abs(F_pi0),std::abs(F_pipm)};
}

}//namespace KM
}//namespace ldm
}//namespace evgen

#endif
