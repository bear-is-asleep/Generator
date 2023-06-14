//____________________________________________________________________________
/*
 Copyright (c) 2003-2023, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Full one loop corrections to nu+e scattering

 Brinden Carlson University of Florida
 bcarlson1@fufl.edu
*/
//____________________________________________________________________________

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/RefFrame.h"
#include "Physics/NuElectron/XSection/NuElectronRadPXSec.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Physics/NuElectron/XSection/PXSecOnElectron.h"
#include "Framework/Numerical/GSLUtils.h"


#include "TFile.h"
#include "TGraph.h"
#include "TMath.h"
#include <Math/Integrator.h>
#include <fstream>
#include <iterator>

using namespace genie;
using namespace genie::constants;
using namespace genie::controls;
using namespace TMath;

//____________________________________________________________________________
NuElectronRadPXSec::NuElectronRadPXSec() :
PXSecOnElectron::PXSecOnElectron("genie::NuElectronRadPXSec","Default")
{

}
//____________________________________________________________________________
NuElectronRadPXSec::NuElectronRadPXSec(string config) :
PXSecOnElectron::PXSecOnElectron("genie::NuElectronRadPXSec", config)
{

}
//____________________________________________________________________________
NuElectronRadPXSec::~NuElectronRadPXSec()
{

}
double NuElectronRadPXSec::Li2(double z) const
{
  double eps = 1e-2;
  double tmin = eps;
  double tmax = 1. - eps;

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("BardinIMD", pDEBUG)
    << "Summing NuElectronRadIntegrand in [" << tmin<< ", " << tmax<< "]";
#endif

  ROOT::Math::IBaseFunctionOneDim * integrand = new
          utils::gsl::wrap::NuElectronRadIntegrand(z);
  ROOT::Math::IntegrationOneDim::Type ig_type =
          utils::gsl::Integration1DimTypeFromString("adaptive");

  double abstol   = 1;    // We mostly care about relative tolerance
  double reltol   = 1E-4;
  int    nmaxeval = 100000;
  ROOT::Math::Integrator ig(*integrand,ig_type,abstol,reltol,nmaxeval);
  double li2 = ig.Integral(tmin, tmax);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("BardinIMD", pDEBUG) << "Li2(z = " << z << ")" << li2;
#endif

  delete integrand;

  return li2;
}

//____________________________________________________________________________
// Auxiliary scalar function for internal integration
//____________________________________________________________________________
utils::gsl::wrap::NuElectronRadIntegrand::NuElectronRadIntegrand(double z):
ROOT::Math::IBaseFunctionOneDim()
{
  fZ = z;
}
//____________________________________________________________________________
utils::gsl::wrap::NuElectronRadIntegrand::~NuElectronRadIntegrand()
{

}
//____________________________________________________________________________
unsigned int utils::gsl::wrap::NuElectronRadIntegrand::NDim(void) const
{
  return 1;
}
//____________________________________________________________________________
double utils::gsl::wrap::NuElectronRadIntegrand::DoEval(double xin) const
{
  if(xin<=0) return 0.;
  if(xin*fZ >= 1.) return 0.;
  double f = TMath::Log(1.-fZ*xin)/xin;
  return f;
}
//____________________________________________________________________________
ROOT::Math::IBaseFunctionOneDim *
  utils::gsl::wrap::NuElectronRadIntegrand::Clone(void) const
{
  return new utils::gsl::wrap::NuElectronRadIntegrand(fZ);
}
//____________________________________________________________________________

double NuElectronRadPXSec::F1(double beta) const{
  return (-1 + (((3 + Sqrt(1 - Power(beta,2)))*Log((1 + beta)/(1 - beta)))/8. - (Log((1 + beta)/(1 - beta))*Log((2*(1 + Sqrt(1 - Power(beta,2))))/Sqrt(1 - Power(beta,2))))/8. + 
      (beta - Log((1 + beta)/(1 - beta))/2.)*Log(me/lambda) + (Li2((1 + beta - Sqrt(1 - Power(beta,2)))/(2.*beta)) - Li2((-1 + beta + Sqrt(1 - Power(beta,2)))/(2.*beta)))/2.
      )/beta);
}
double NuElectronRadPXSec::F2(double beta) const{
  return ((Sqrt(1 - Power(beta,2))*Log((1 + beta)/(1 - beta)))/(4.*beta));
}
double NuElectronRadPXSec::IL0() const{
  return 1; //
}
double NuElectronRadPXSec::IR0( double omega, double omegap) const{
  return Power(omegap/omega,2);
}
double NuElectronRadPXSec::ILR0( double omega, double omegap) const{
  return -me/omega*(1-omegap/omega);
}
double NuElectronRadPXSec::IXX( double omega, double omegap) const{
  return ((-3 - Power(omegap,2)/Power(omega,2) + 2*(1 - omegap/omega) + (2*me*(1 - omegap/omega))/omega + Power(1 - omegap/omega,2))/4.);
}

//Coupling corrections
// double NuElectronRadPXSec::ax( double omega, double omegap) const{
//   return ((omega*(-7*me + 8*omega - 4*omegap)*Power(omegap,3))/(5.*Power(me,3)) + (omega*
//       (-15*Power(me,4) - (15*Power(me,5))/(me - 2*omegap) + (30*Power(me,4)*(me + 2*omega))/(me + 2*omega - 2*omegap) - 2*(me - 2*omega)*(11*Power(me,2) + 4*Power(omega,2))*omegap + 
//         4*(109*Power(me,2) + 78*me*omega + 2*Power(omega,2))*Power(omegap,2)))/(120.*Power(me,3)));
// }
// double NuElectronRadPXSec::bL( double omega, double omegap) const{
//   return (-2*Power(omega,2) - Power(omega,3)/(3.*me) + (2*Power(omega,5))/(15.*Power(me,3)));
// }

// double NuElectronRadPXSec::cL( double omega, double omegap) const{
//   return (Power(me,2)/60. - (9*me*omega)/8. - (5*Power(omega,2))/3. + Power(omega,3)/(3.*me) - (2*Power(omega,5))/(15.*Power(me,3)));
// }

// double NuElectronRadPXSec::dL( double omega, double omegap) const{
//   return ((2*Power(me,5) - 135*Power(me,4)*omega - 16*Power(omega - omegap,3)*(Power(omega,2) + 3*omega*omegap + 6*Power(omegap,2)))/(120.*Power(me,3)) + 
//    (3*Power(omega - omegap,2)*Power(omegap,2) + Power(me,2)*(-5*Power(omega,2) + 6*omega*omegap + 6*Power(omegap,2)) + 
//       me*(omega - omegap)*(Power(omega,2) + 7*omega*omegap + 13*Power(omegap,2)))/(3.*Power(me,2)));
// }

// double NuElectronRadPXSec::axLR( double omega, double omegap) const{
//   return ((4*omega*(-13*Power(me,2) + Power(-3*me + omega - omegap,2))*omegap)/(3.*me*(-me + 2*(me + omega - omegap))));
// }

// double NuElectronRadPXSec::bLR( double omega, double omegap) const{
//   return (2*(-Power(omega,2) + Power(omega,3)/(3.*me) - me*omegap + omega*omegap));
// }

// double NuElectronRadPXSec::cLR( double omega, double omegap) const{
//   return (2*Power(me,2))/3. - (2*Power(omega,3))/(3.*me) + 2*omega*(omega - omegap) + me*(3*omega - omegap);
// }

// double NuElectronRadPXSec::dLR( double omega, double omegap) const{
//   return (2*Power(me,2))/3. - (2*Power(omega,3))/(3.*me) + 2*omega*(omega - omegap) + me*(3*omega - omegap) + (2*omegap*(-3*(me + omega - omegap) + Power(omegap,2)/me))/3.;
// }

double NuElectronRadPXSec::ILfelectronexact(double beta, double omega) const{
  return ((Power(kPi,2)*((((-8*Power(beta,2))/(15.*(1 - Power(beta,2))) + ((-49 + 25*Power(beta,2))*(1 - 1/Sqrt(1 - Power(beta,2))))/(60.*Power(1 - Power(beta,2),1.5)))*Power(me,4) + 
          ((-105 + 38*beta + 51*Power(beta,2) - 20*Power(beta,3))/(60.*Power(1 - Power(beta,2),1.5)) - 
             (-105 - 82*beta + 54*Power(beta,2) + 55*Power(beta,3))/(60.*(1 + beta)*(1 - Power(beta,2))))*Power(me,3)*omega + 
          ((-23 + 8*beta + 7*Power(beta,2))/(30.*(1 + beta)*Sqrt(1 - Power(beta,2))) - (-23 + 6*beta + 15*Power(beta,2))/(30.*(1 - Power(beta,2))))*Power(me,2)*Power(omega,2) + 
          (-0.03333333333333333*(3 + 2*beta)/(1 + beta) + (3 - beta)/(30.*Sqrt(1 - Power(beta,2))))*me*Power(omega,3) + ((1 - Sqrt(1 - Power(beta,2))/(1 + beta))*Power(omega,4))/15.)/
        Power(me,2) + (((1 - Sqrt(1 - Power(beta,2)))*Power(me,2))/(2.*beta) + (beta*me*omega)/(2.*Sqrt(1 - Power(beta,2))) + 
          (-0.5*(1 + beta)/beta + Sqrt(1 - Power(beta,2))/(2.*beta*(1 + beta)))*Power(omega,2))*Log((1 + beta)/(1 - beta)) + 
       (((-22 + 14*Power(beta,2) - 7*Power(beta,4))/(15.*Power(1 - Power(beta,2),2)) + (22 - 25*Power(beta,2) + 15*Power(beta,4))/(15.*Power(1 - Power(beta,2),2.5)))*Power(me,2) + 
          (2/Power(1 - Power(beta,2),1.5) - (2 - Power(beta,2) + Power(beta,4))/Power(1 - Power(beta,2),2))*me*omega + 
          ((1 + 3*Power(beta,2))/(3.*Power(1 - Power(beta,2),1.5)) - (7 - 11*Power(beta,2) + 4*Power(beta,4))/(3.*Power(1 - Power(beta,2),2)))*Power(omega,2) + Power(omega,3)/(3.*me) - 
          (2*Power(omega,5))/(15.*Power(me,3)))*Log((2*(me - me/Sqrt(1 - Power(beta,2)) + omega))/(me*(1 - (1 + beta)/Sqrt(1 - Power(beta,2)) + (2*omega)/me))) + 
       (((44 + 43*beta + 14*Power(beta,2))/(60.*Power(1 + beta,2)) - (Sqrt(1 - Power(beta,2))*(22 + 36*beta + 17*Power(beta,2)))/(30.*Power(1 + beta,3)))*Power(me,2) + 
          ((-1 + beta + Power(beta,2))/((1 + beta)*Sqrt(1 - Power(beta,2))) + (1 + Power(1 + beta,-2))/2.)*me*omega + 
          ((4 + beta)/(6.*(1 + beta)) + ((-2 - beta)*Sqrt(1 - Power(beta,2)))/(3.*Power(1 + beta,2)))*Power(omega,2))*
        Log((1 - (1 + beta)/Sqrt(1 - Power(beta,2)))/((1 + beta)/(1 - beta) - ((1 + beta)*(1 + (2*omega)/me))/Sqrt(1 - Power(beta,2)))) + 
       (omega*(-me + omega)*Log((2*omega)/(me*(-1 + (Sqrt(1 - Power(beta,2))*(1 + (2*omega)/me))/(1 + beta)))))/2. + 
       ((Power(me,2)/2. + 2*me*omega + Power(omega,2))*(-0.16666666666666666*Power(kPi,2) + 
            Li2((1 + beta)/Sqrt(1 - Power(beta,2))) - Li2(1 + (2*omega)/me) + Li2((Sqrt(1 - Power(beta,2))*(1 + (2*omega)/me))/(1 + beta))))/2.))/Power(omega,3));
}

double NuElectronRadPXSec::IRfelectronexact(double beta, double omega) const{
  return ((Power(kPi,2)*((((-269 + 270*Power(beta,2))/(60.*Power(1 - Power(beta,2),1.5)) + (538 - 839*Power(beta,2) + 309*Power(beta,4))/(120.*Power(1 - Power(beta,2),2)))*Power(me,5) + 
          ((-439 + beta*(44 + (433 - 45*beta)*beta))/(30.*Power(1 - Power(beta,2),1.5)) + (439 + beta*(-29 + beta*(-730 + 9*beta*(3 + 33*beta))))/(30.*Power(1 - Power(beta,2),2)))*
           Power(me,4)*omega + ((851 + 778*beta - 825*Power(beta,2) - 760*Power(beta,3))/(60.*(1 + beta)*(1 - Power(beta,2))) - 
             (851 - 353*beta - 805*Power(beta,2) + 339*Power(beta,3))/(60.*Power(1 - Power(beta,2),1.5)))*Power(me,3)*Power(omega,2) + 
          ((-133 - 2*beta + 113*Power(beta,2))/(30.*(1 + beta)*Sqrt(1 - Power(beta,2))) - (-133 - 34*beta + 143*Power(beta,2))/(30.*(1 - Power(beta,2))))*Power(me,2)*Power(omega,3) + 
          ((8 - beta)/(15.*(1 + beta)) - 8/(15.*Sqrt(1 - Power(beta,2))))*me*Power(omega,4) + (2*(1 - Sqrt(1 - Power(beta,2))/(1 + beta))*Power(omega,5))/15.)/
        (Power(me,2)*(me + 2*omega)) + (((-1 + beta + 2*Power(beta,2) - Power(beta,3))/(2.*beta*(1 + beta)*Sqrt(1 - Power(beta,2))) + 
             (2 - 4*beta + Power(beta,2) + 4*Power(beta,3))/(4.*beta*(1 - Power(beta,2))))*Power(me,2) + 
          (2 - 1/(1 + beta) - (2 - beta)/(2.*Sqrt(1 - Power(beta,2))))*me*(me - me/Sqrt(1 - Power(beta,2)) + omega) + 
          (-0.5*(1 + beta)/beta + Sqrt(1 - Power(beta,2))/(2.*beta*(1 + beta)))*Power(me - me/Sqrt(1 - Power(beta,2)) + omega,2))*Log((1 + beta)/(1 - beta)) - 
       (Power(me - me/Sqrt(1 - Power(beta,2)) + omega,2)*(-30*Power(me,3) + 35*Power(me,2)*(me - me/Sqrt(1 - Power(beta,2)) + omega) - 
            10*me*Power(me - me/Sqrt(1 - Power(beta,2)) + omega,2) + 2*Power(me - me/Sqrt(1 - Power(beta,2)) + omega,3))*
          Log((2*(me - me/Sqrt(1 - Power(beta,2)) + omega))/(me*(1 - (1 + beta)/Sqrt(1 - Power(beta,2)) + (2*omega)/me))))/(15.*Power(me,3)) + 
       (((-2 + 41*beta + 14*Power(beta,2) - 23*Power(beta,3))/(30.*Power(1 + beta,2)*Sqrt(1 - Power(beta,2))) + 
             (2 + 43*Power(beta,2) - 28*beta*(1 - Power(beta,2)))/(30.*(1 + beta)*(1 - Power(beta,2))))*Power(me,2) + 
          ((1 - 5*beta + Power(beta,2))/(3.*(1 + beta)*Sqrt(1 - Power(beta,2))) + (-2 + 8*beta + 7*Power(beta,2))/(6.*Power(1 + beta,2)))*me*(me - me/Sqrt(1 - Power(beta,2)) + omega) + 
          ((4 + beta)/(6.*(1 + beta)) + ((-2 - beta)*Sqrt(1 - Power(beta,2)))/(3.*Power(1 + beta,2)))*Power(me - me/Sqrt(1 - Power(beta,2)) + omega,2))*
        Log((1 - (1 + beta)/Sqrt(1 - Power(beta,2)))/((1 + beta)/(1 - beta) - ((1 + beta)*(1 + (2*omega)/me))/Sqrt(1 - Power(beta,2)))) + 
       ((((-2 + Power(beta,2) + 2*Sqrt(1 - Power(beta,2)))*Power(me,4))/(1 - Power(beta,2)) + 
            ((-10 + 6*Power(beta,2) + 9*Sqrt(1 - Power(beta,2)))*Power(me,3)*omega)/(1 - Power(beta,2)) + 
            ((-16 + 12*Power(beta,2) + 11*Sqrt(1 - Power(beta,2)))*Power(me,2)*Power(omega,2))/(1 - Power(beta,2)) - 2*(5 - 1/Sqrt(1 - Power(beta,2)))*me*Power(omega,3) - Power(omega,4)
            )*Log((2*omega)/(me*(-1 + (Sqrt(1 - Power(beta,2))*(1 + (2*omega)/me))/(1 + beta)))))/Power(me + 2*omega,2) + 
       ((((Power(beta,2) + Sqrt(1 - Power(beta,2)))*Power(me,2))/(1 - Power(beta,2)) + Power(me - me/Sqrt(1 - Power(beta,2)) + omega,2))*
          (-0.16666666666666666*Power(kPi,2) + Li2((1 + beta)/Sqrt(1 - Power(beta,2))) - Li2(1 + (2*omega)/me) + 
              Li2((Sqrt(1 - Power(beta,2))*(1 + (2*omega)/me))/(1 + beta))))/2.))/Power(omega,3));
}
double NuElectronRadPXSec::ILRfelectronexact(double beta, double omega) const{
  return ((Power(kPi,2)*(((me - me/Sqrt(1 - Power(beta,2)) + omega - (Sqrt(1 - Power(beta,2))*omega)/(1 + beta))*(9*me + 2*(me - me/Sqrt(1 - Power(beta,2)) + omega)))/6. + 
       ((beta*Power(me,2))/Sqrt(1 - Power(beta,2)) + ((me - me/Sqrt(1 - Power(beta,2)))*((1 + beta/2.)*Power(me,2) - 2*Sqrt(1 - Power(beta,2))*me*omega + (1 - beta)*Power(omega,2)))/
           (beta*me))*Log((1 + beta)/(1 - beta)) + ((3*Power(me,2)*omega + 3*Power(me,2)*(me - me/Sqrt(1 - Power(beta,2)) + omega) - 
            3*me*Power(me - me/Sqrt(1 - Power(beta,2)) + omega,2) - 2*Power(me - me/Sqrt(1 - Power(beta,2)) + omega,3))*
          Log((2*(me - me/Sqrt(1 - Power(beta,2)) + omega))/(me*(1 - (1 + beta)/Sqrt(1 - Power(beta,2)) + (2*omega)/me))))/(3.*me) + 
       (((14 + 5*beta + (2*(-7 - 4*beta + 2*Power(beta,2)))/Sqrt(1 - Power(beta,2)))*Power(me,2))/(6.*(1 + beta)) + (1 + (1 - 2*Sqrt(1 - Power(beta,2)))/(1 + beta))*me*omega)*
        Log((1 - (1 + beta)/Sqrt(1 - Power(beta,2)))/((1 + beta)/(1 - beta) - ((1 + beta)*(1 + (2*omega)/me))/Sqrt(1 - Power(beta,2)))) + 
       (Power(me,2)*(1 - (Sqrt(1 - Power(beta,2))*(-(me*omega) + Power(me + 2*omega,2)))/(me*(me + 2*omega)))*
          Log((2*omega)/(me*(-1 + (Sqrt(1 - Power(beta,2))*(1 + (2*omega)/me))/(1 + beta)))))/Sqrt(1 - Power(beta,2)) + 
       (me*(-me + 2*(me - me/Sqrt(1 - Power(beta,2)) + omega))*(-0.16666666666666666*Power(kPi,2) + 
            Li2((1 + beta)/Sqrt(1 - Power(beta,2))) - Li2(1 + (2*omega)/me) + Li2((Sqrt(1 - Power(beta,2))*(1 + (2*omega)/me))/(1 + beta))))/2.))/Power(omega,3));
}


double NuElectronRadPXSec::PiFermiLoop(double mf, double Q2) const{
  return (5/9 - (4*Power(mf,2))/(3.*Q2) + Log(Power(mu,2)/Power(mf,2))/3. + 
   ((1 - (2*Power(mf,2))/Q2)*Sqrt(1 + (4*Power(mf,2))/Q2)*Log((-1 + Sqrt(1 + (4*Power(mf,2))/Q2))/(1 + Sqrt(1 + (4*Power(mf,2))/Q2))))/3.);
}

double NuElectronRadPXSec::Pi0FermiLoop(double mf) const{
  return ((alphas*(13/12 - Log(Power(mu,2)/Power(mf,2))))/(3.*kPi) + Log(Power(mu,2)/Power(mf,2))/3. + 
   (Power(alphas,2)*(-3847/864 - (5*Log(Power(mu,2)/Power(mf,2)))/6. - (11*Power(Log(Power(mu,2)/Power(mf,2)),2))/8. + 
        4*(361/1296 - Log(Power(mu,2)/Power(mf,2))/18. + Power(Log(Power(mu,2)/Power(mf,2)),2)/12.) + (655*std::riemann_zeta(3))/144.))/(3.*Power(kPi,2)));
}
double NuElectronRadPXSec::Pi0rFermiLoop(double mf) const{
  return (1/3 * (1 + alphas/kPi) * Log(Power(mu,2)/Power(mf,2)));
}

//soft photon correction
double NuElectronRadPXSec::deltasoft(double beta) const{ 
  return (1 + (Log((1 + beta)/(1 - beta))*(1 + Log(((1 + beta)*Sqrt(1 - Power(beta,2)))/(4.*Power(beta,2)))))/(2.*beta) - 
   (2*(beta - Log((1 + beta)/(1 - beta))/2.)*Log((2*epsilon)/lambda))/beta + (-1/6*Power(kPi,2) + Li2((1 - beta)/(1 + beta)))/beta);
}
double NuElectronRadPXSec::deltaemIR1(double beta, double omega, double omegap) const{
  return((2*(beta - Log((1 + beta)/(1 - beta))/2.)*Log((2*(1 + beta)*epsilon)/
       (beta*me*(1 + (Sqrt(1 - Power(beta,2))*(-((Power(beta,2)*Power(me,2))/(1 - Power(beta,2))) + Power(omega,2) - Power(omegap,2)))/(2.*beta*me*omegap)))))/beta);
}
double NuElectronRadPXSec::deltaemIR2(double beta, double omega, double omegap) const{
  return (-1 + Log((1 - (Sqrt(1 - Power(beta,2))*(-((Power(beta,2)*Power(me,2))/(1 - Power(beta,2))) + Power(omega,2) - Power(omegap,2)))/(2.*me*omegap))/Sqrt(1 - Power(beta,2))) + 
   (Power(kPi,2)/6. + Log((1 - beta)/(1 + beta))*(0.5 + Log((Sqrt(1 - Power(beta,2))*
             (1 + (Sqrt(1 - Power(beta,2))*(-((Power(beta,2)*Power(me,2))/(1 - Power(beta,2))) + Power(omega,2) - Power(omegap,2)))/(2.*beta*me*omegap)))/(4.*beta))) - 
      Li2((1 - beta)/(1 + beta)) - Li2((-1 + (Sqrt(1 - Power(beta,2))*(-((Power(beta,2)*Power(me,2))/(1 - Power(beta,2))) + Power(omega,2) - Power(omegap,2)))/
           (2.*beta*me*omegap))/(1 + (Sqrt(1 - Power(beta,2))*(-((Power(beta,2)*Power(me,2))/(1 - Power(beta,2))) + Power(omega,2) - Power(omegap,2)))/(2.*beta*me*omegap))) + 
      Li2(((1 + beta)*(-1 + (Sqrt(1 - Power(beta,2))*(-((Power(beta,2)*Power(me,2))/(1 - Power(beta,2))) + Power(omega,2) - Power(omegap,2)))/(2.*beta*me*omegap)))/
        ((1 - beta)*(1 + (Sqrt(1 - Power(beta,2))*(-((Power(beta,2)*Power(me,2))/(1 - Power(beta,2))) + Power(omega,2) - Power(omegap,2)))/(2.*beta*me*omegap)))))/beta);
}
//Muon differential xsec wrt outgoing neutrino energy
double NuElectronRadPXSec::dsigmadE(double X, double omega, double cLe, double cLe_mu, double cRe) const{
  double Ep = omega * (1 - X/(2*me));
  double omegap = me + omega - Ep;
  double beta = Sqrt(Power(Ep,2)-Power(me,2))/Ep;
  double c2 = Power(2*Sqrt(2)*kGF,2);

  return (c2*me/(4*kPi)*((omega*alpha*(Power(cLe_mu,2)*ILfelectronexact(beta,omega) + 
          cLe_mu*cRe*ILRfelectronexact(beta,omega) + 
          Power(cRe,2)*IRfelectronexact(beta,omega)))/Power(kPi,3) + 
     (alpha*F2(beta)*
        (Power(cLe_mu,2)*ILR0(omega,omegap) + Power(cRe,2)*ILR0(omega,omegap) + 
          2*cLe_mu*cRe*(IL0() + IR0(omega,omegap)) + 
          Power(cLe_mu + cRe,2)*IXX(omega,omegap)))/kPi + 
     (Power(cLe_mu,2)*IL0() + cLe_mu*cRe*ILR0(omega,omegap) + 
        Power(cRe,2)*IR0(omega,omegap))*(1 + 
        (alpha*(2*F1(beta) + 
             deltaemIR1(beta,omega,omegap) + 
             deltaemIR2(beta,omega,omegap) + 
             deltasoft(beta)))/kPi) + 
     (alpha*(cLe_mu*IL0() + ((cLe_mu + cRe)*ILR0(omega,omegap))/2. + 
          cRe*IR0(omega,omegap))*(Sqrt(c2)*Pi_gg*(1 - 2*fSin28w) + (cLe_mu + cRe)*Qe*PiFermiLoop(me,2*me*(-me + Ep)) + 
          (cLe + cRe)*Qe*PiFermiLoop(kMuonMass,2*me*(-me + Ep)) + (cLe_mu + cRe)*Qe*PiFermiLoop(kTauMass,2*me*(-me + Ep)) + 
          (cLu + cRu)*Nc*Qu*Pi0rFermiLoop(mcmu) + 
          2*(cLd + cRd)*Nc*Qd*Pi0rFermiLoop(2.) + (cLu + cRu)*Nc*Qu*Pi0rFermiLoop(2.)))/kPi));
}

//____________________________________________________________________________
double NuElectronRadPXSec::XSec(
                 const Interaction * interaction, KinePhaseSpace_t kps) const
{
  /*
  Two possible processes for scattering on electron with electron in final state

  Eqns. for processes in https://arxiv.org/pdf/1907.03379.pdf at LO

  nul    + e -> nul    + e (l=e,u,t)                     (14)
  nulbar + e -> nulbar + e (l=e,u,t)                     (15)
  */
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  //----- get initial state & kinematics
  const InitialState & init_state = interaction -> InitState();
  const Kinematics &   kinematics = interaction -> Kine();
  const ProcessInfo &  proc_info  = interaction -> ProcInfo();

  double Ev = init_state.ProbeE(kRfHitElRest); //Electron rest frame
  double Ev2 = Power(Ev,2);
  double y  = kinematics.y();
  double El = Ev*(1-y);
  double X = 2*me*(1-El/Ev);


  //---- integration limits over final state electron energy
  double Emin = me;
  double Emax = me + 2*Ev2/(me+2*Ev);
  if (El<Emin || El>Emax) return 0;

  double xsec = 0; // <-- dxsec/dy

  int inu = init_state.ProbePdg();
  if ((pdg::IsNuMu(inu) || pdg::IsNuTau(inu)) && proc_info.IsWeakNC()){
    xsec = NuElectronRadPXSec::dsigmadE(X, Ev, cLl, cLl_mu, cRl);
  }

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("Elastic", pDEBUG)
    << "*** dxsec(ve-)/dy [free e-](Ev="<< Ev << ", y= "<< y<< ") = "<< xsec;
#endif

  //----- The algorithm computes dxsec/dy
  //      Check whether variable tranformation is needed
  if(kps!=kPSyfE) {
    double J = utils::kinematics::Jacobian(interaction,kPSyfE,kps);
    xsec *= J;
  }

  //----- If requested return the free electron xsec even for nuclear target
  if( interaction->TestBit(kIAssumeFreeElectron) ) return xsec;

  //----- Scale for the number of scattering centers at the target
  int Ne = init_state.Tgt().Z(); // num of scattering centers
  xsec *= Ne;

  return xsec;
}
//____________________________________________________________________________
void NuElectronRadPXSec::Configure(const Registry & config)
{
  PXSecOnElectron::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void NuElectronRadPXSec::Configure(string config)
{
  PXSecOnElectron::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void NuElectronRadPXSec::LoadConfig(void)
{
  // weinberg angle
  double thw ;
  GetParam( "WeinbergAngle", thw ) ;
  fSin28w = TMath::Power(TMath::Sin(thw), 2);
  fSin48w = TMath::Power(TMath::Sin(thw), 4);
  if (!fElectronVelocity) {
    std::cerr << "Error: fElectronVelocity is not initialized correctly." << std::endl;
  }

  //Constants - Coupling constants reported at mass scale mu = 2 GeV
  GetParam("alphas",alphas); 
  GetParam("alpha",alpha);
  GetParam("mu",mu); //GeV
  GetParam("Nc",Nc);
  GetParam("lambda",lambda);
  GetParam("epsilon",epsilon);
  GetParam("lambdac",lambdac);
  GetParam("mcmu",mcmu);
  GetParam("Pi_gg",Pi_gg);
  GetParam("Qe",Qe);
  GetParam("Qd",Qd);
  GetParam("Qu",Qu);
  GetParam("cLl",cLl);
  GetParam("cLl_mu",cLl_mu);
  GetParam("cRl",cRl);
  GetParam("cLu",cLu);
  GetParam("cRu",cRu);
  GetParam("cLd",cLd);
  GetParam("cRd",cRd);
}
//____________________________________________________________________________

