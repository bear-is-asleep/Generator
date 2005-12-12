//____________________________________________________________________________
/*!

\namespace genie::pdg

\brief     Utilities for improving the code readability when using PDG codes.

           E.g. a ..............................if( pdg::IsProton(pdg_code) )
           is much easier to read/check than....if( pdgc_code == 2212 )

\author    Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
           CCLRC, Rutherford Appleton Laboratory

\created   May 06, 2004

*/
//____________________________________________________________________________

#include <cassert>

#include "BaryonResonance/BaryonResUtils.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGCodes.h"

using namespace genie;

//____________________________________________________________________________
int genie::pdg::IonPdgCodeToZ(int ion_pdgc)
{
// Decoding Z from the PDG code
// GENIE uses the MINOS PDG extension for ions: PDG Code = [1AAAZZZ000]

  int Z = (ion_pdgc/1000) - 1000*(ion_pdgc/1000000); // don't factor out
  return Z;
}
//____________________________________________________________________________
int genie::pdg::IonPdgCodeToA(int ion_pdgc)
{
// Decoding A from the PDG code
// GENIE uses the MINOS PDG extension for ions: PDG Code = [1AAAZZZ000]

  int A = (ion_pdgc/1000000) - 1000;
  return A;
}
//____________________________________________________________________________
int genie::pdg::IonPdgCode(int A, int Z)
{
  int ion_pdgc = (1000000000) + (A*1000000) + (Z*1000);
  return ion_pdgc;
}
//____________________________________________________________________________
bool genie::pdg::IsIon(int pdgc)
{
  return (pdgc > 1000000000 && pdgc < 1999999999);
}
//____________________________________________________________________________
bool genie::pdg::IsNeutrino(int pdgc)
{
  bool is_nu = (pdgc == kPdgNuE)  ||
               (pdgc == kPdgNuMu) ||
               (pdgc == kPdgNuTau);
  return is_nu;
}
//____________________________________________________________________________
bool genie::pdg::IsAntiNeutrino(int pdgc)
{
  bool is_nubar = (pdgc == kPdgNuEBar)  ||
                  (pdgc == kPdgNuMuBar) ||
                  (pdgc == kPdgNuTauBar);

  return is_nubar;
}
//____________________________________________________________________________
bool genie::pdg::IsNeutralLepton(int pdgc)
{
  bool is_neutral_lepton = IsNeutrino(pdgc) || IsAntiNeutrino(pdgc);
  return is_neutral_lepton;
}
//____________________________________________________________________________
bool genie::pdg::IsChargedLepton(int pdgc)
{
  bool is_neg_lepton = (pdgc ==  kPdgElectron) ||
                       (pdgc ==  kPdgMuon)     ||
                       (pdgc ==  kPdgTau);

  bool is_pos_lepton = (pdgc ==  kPdgPositron) ||
                       (pdgc ==  kPdgAntiMuon) ||
                       (pdgc ==  kPdgAntiTau);

  bool is_charged_lepton = is_neg_lepton || is_pos_lepton;
  return is_charged_lepton;
}
//____________________________________________________________________________
bool genie::pdg::IsNuE(int pdgc)
{
  return (pdgc == kPdgNuE);
}
//____________________________________________________________________________
bool genie::pdg::IsNuMu(int pdgc)
{
  return (pdgc == kPdgNuMu);
}
//____________________________________________________________________________
bool genie::pdg::IsNuTau(int pdgc)
{
  return (pdgc == kPdgNuTau);
}
//____________________________________________________________________________
bool genie::pdg::IsAntiNuE(int pdgc)
{
  return (pdgc == kPdgNuEBar);
}
//____________________________________________________________________________
bool genie::pdg::IsAntiNuMu(int pdgc)
{
  return (pdgc == kPdgNuMuBar);
}
//____________________________________________________________________________
bool genie::pdg::IsAntiNuTau(int pdgc)
{
  return (pdgc == kPdgNuTauBar);
}
//____________________________________________________________________________
bool genie::pdg::IsElectron(int pdgc)
{
  return (pdgc == kPdgElectron);
}
//____________________________________________________________________________
bool genie::pdg::IsPositron(int pdgc)
{
  return (pdgc == kPdgPositron);
}
//____________________________________________________________________________
bool genie::pdg::IsMuon(int pdgc)
{
  return (pdgc == kPdgMuon);
}
//____________________________________________________________________________
bool genie::pdg::IsAntiMuon(int pdgc)
{
  return (pdgc == kPdgAntiMuon);
}
//____________________________________________________________________________
bool genie::pdg::IsTau(int pdgc)
{
  return (pdgc == kPdgTau);
}
//____________________________________________________________________________
bool genie::pdg::IsAntiTau(int pdgc)
{
  return (pdgc == kPdgAntiTau);
}
//____________________________________________________________________________
int genie::pdg::Neutrino2ChargedLepton(int pdgc)
{
  switch(pdgc) {
       case (kPdgNuE)      : return kPdgElectron; break;
       case (kPdgNuEBar)   : return kPdgPositron; break;
       case (kPdgNuMu)     : return kPdgMuon;     break;
       case (kPdgNuMuBar)  : return kPdgAntiMuon; break;
       case (kPdgNuTau)    : return kPdgTau;      break;
       case (kPdgNuTauBar) : return kPdgAntiTau;  break;
  }
  return -1;
}
//____________________________________________________________________________
bool genie::pdg::IsQuark(int pdgc)
{
  return ( pdgc == kPdgDQuark || pdgc == kPdgUQuark ||
           pdgc == kPdgSQuark || pdgc == kPdgCQuark ||
           pdgc == kPdgBQuark || pdgc == kPdgTQuark
         );
}
//____________________________________________________________________________
bool genie::pdg::IsAntiQuark(int pdgc)
{
  return ( pdgc == kPdgDQuarkBar || pdgc == kPdgUQuarkBar ||
           pdgc == kPdgSQuarkBar || pdgc == kPdgCQuarkBar ||
           pdgc == kPdgBQuarkBar || pdgc == kPdgTQuarkBar
         );
}
//____________________________________________________________________________
bool genie::pdg::IsUQuark(int pdgc)
{
  return (pdgc == kPdgUQuark);
}
//____________________________________________________________________________
bool genie::pdg::IsDQuark(int pdgc)
{
  return (pdgc == kPdgDQuark);
}
//____________________________________________________________________________
bool genie::pdg::IsSQuark(int pdgc)
{
  return (pdgc == kPdgSQuark);
}
//____________________________________________________________________________
bool genie::pdg::IsCQuark(int pdgc)
{
  return (pdgc == kPdgCQuark);
}
//____________________________________________________________________________
bool genie::pdg::IsUAntiQuark(int pdgc)
{
  return (pdgc == kPdgUQuarkBar);
}
//____________________________________________________________________________
bool genie::pdg::IsDAntiQuark(int pdgc)
{
  return (pdgc == kPdgDQuarkBar);
}
//____________________________________________________________________________
bool genie::pdg::IsSAntiQuark(int pdgc)
{
  return (pdgc == kPdgSQuarkBar);
}
//____________________________________________________________________________
bool genie::pdg::IsCAntiQuark(int pdgc)
{
  return (pdgc == kPdgCQuarkBar);
}
//____________________________________________________________________________
bool genie::pdg::IsProton(int pdgc)
{
  return (pdgc == kPdgProton);
}
//____________________________________________________________________________
bool genie::pdg::IsNeutron(int pdgc)
{
  return (pdgc == kPdgNeutron);
}
//____________________________________________________________________________
bool genie::pdg::IsNeutronOrProton(int pdgc)
{
  return (pdgc == kPdgNeutron || pdgc == kPdgProton);
}
//____________________________________________________________________________
int genie::pdg::SwitchProtonNeutron(int pdgc)
{
  assert(IsProton(pdgc) || IsNeutron(pdgc));

  if (IsProton(pdgc)) return kPdgNeutron;
  else                return kPdgProton;
}
//____________________________________________________________________________
bool genie::pdg::IsBaryonResonance(int pdgc)
{
  return utils::res::IsBaryonResonance(pdgc);
}
//____________________________________________________________________________

