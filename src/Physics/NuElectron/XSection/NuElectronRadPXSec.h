//____________________________________________________________________________
/*!

\class    genie::NuElectronPXSec

\brief    nu/nubar + e- scattering differential cross section \n
          The cross section algorithm handles:
             - nue/nuebar + e- -> nue/nuebar + e- [CC + NC + interference]
             - numu/nutau + e- -> numu/nutau + e- [NC]
             - numubar/nutaubar + e- -> numubar/nutaubar + e- [NC]

          Is a concrete implementation of the XSecAlgorithmI interface. \n

\ref      Theory of elastic neutrino-electron scattering
          Oleksandr Tomalak and Richard J. Hill
          https://arxiv.org/pdf/1907.03379.pdf

\author   Brinden Carlson bcarlson1@ufl.edu

\created  June 1, 2023

\cpright  Copyright (c) 2003-2023, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _NU_ELECTRON_PARTIAL_XSEC_H_
#define _NU_ELECTRON_PARTIAL_XSEC_H_

#include "Physics/NuElectron/XSection/PXSecOnElectron.h"
#include <Math/IFunction.h>

namespace genie {

class NuElectronRadPXSec : public PXSecOnElectron {

public:
  NuElectronRadPXSec();
  NuElectronRadPXSec(string config);
  virtual ~NuElectronRadPXSec();

  //-- XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t k) const;

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:
  void LoadConfig (void);

  //Second order polylog
  double Li2(double z) const;

  //Vertex form factors
  double F1(double beta) const;
  double F2(double beta) const;


  //Leading order kinematic factors
  double IL0() const;
  double IR0( double omega, double omegap) const;
  double ILR0( double omega, double omegap) const;
  double IXX( double omega, double omegap) const;

  //EFT closed fermion loop factors - dep. on renormlization scale mu MS bar and strong coupling alphas
  double PiFermiLoop(double mf, double Q2) const;
  double Pi0FermiLoop(double mf) const;
  double Pi0rFermiLoop(double mf) const;

  //Soft photon and real photon radiation correction
  double deltasoft(double beta) const;
  double deltaemIR1(double beta,  double omega, double omegap) const;
  double deltaemIR2(double beta,  double omega, double omegap) const;

  //Kinematic corrections
  double ILfelectronexact(double beta,  double omega) const;
  double IRfelectronexact(double beta,  double omega) const;
  double ILRfelectronexact(double beta,  double omega) const;

  //Differential cross sections - do we need four functions or just update the couplings?
  double dsigmadE(double X, double omega, double cLe, double cLe_mu, double cRe) const;

  double fSin28w; // sin^2(theta-weinberg)
  double fSin48w;
  double me = genie::constants::kElectronMass;

  double alphas; 
  double alpha;
  double mu; //GeV
  double Nc;
  double lambda;
  double epsilon;
  double lambdac;
  double mcmu;
  double Pi_gg;
  //Couplings and charges
  double Qe;
  double Qd;
  double Qu;
  double cLl;
  double cLl_mu;
  double cRl;
  double cLu;
  double cRu;
  double cLd;
  double cRd;
};
}       // genie namespace

//____________________________________________________________________________
/*!
\class    genie::BardinIMDRadCorIntegrand

\brief    Auxiliary scalar function for the internal integration in Bardin's
          IMD d2xsec/dxdy cross section algorithm

\author   Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
          University of Liverpool & STFC Rutherford Appleton Laboratory

\created  February 20, 2006
*/
//____________________________________________________________________________

namespace genie {
 namespace utils {
  namespace gsl   {
   namespace wrap   {

    class NuElectronRadIntegrand : public ROOT::Math::IBaseFunctionOneDim
    {
     public:
       NuElectronRadIntegrand(double z);
      ~NuElectronRadIntegrand();
       // ROOT::Math::IBaseFunctionOneDim interface
       unsigned int                      NDim   (void)       const;
       double                            DoEval (double xin) const;
       ROOT::Math::IBaseFunctionOneDim * Clone  (void)       const;
     private:
       double fZ;
     };

   } // wrap namespace
  } // gsl namespace
 } // utils namespace
} // genie namespace


#endif  // _NU_ELECTRON_PARTIAL_XSEC_H_
