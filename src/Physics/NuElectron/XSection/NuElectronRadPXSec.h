//____________________________________________________________________________
/*!

\class    genie::NuElectronPXSec

\brief    nu/nubar + e- scattering differential cross section \n
          The cross section algorithm handles:
             - nue/nuebar + e- -> nue/nuebar + e- [CC + NC + interference]
             - numu/nutau + e- -> numu/nutau + e- [NC]
             - numubar/nutaubar + e- -> numubar/nutaubar + e- [NC]
             - numu/nutau + e- -> l- + nu_e [CC]

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

namespace genie {

class NuElectronPXSec : public PXSecOnElectron {

public:
  NuElectronPXSec();
  NuElectronPXSec(string config);
  virtual ~NuElectronPXSec();

  //-- XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t k) const;

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:
  void LoadConfig (void);

  //Second order polylog
  double NuElectronRadPXSec::Li2(double z) const

  //

  //const ElectronVelocity * fElectronVelocity;

  double fSin28w; // sin^2(theta-weinberg)
  double fSin48w;
};



}       // genie namespace
#endif  // _NU_ELECTRON_PARTIAL_XSEC_H_
