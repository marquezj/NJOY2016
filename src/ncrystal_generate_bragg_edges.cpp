////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "NCrystal/NCrystal.hh"
#include "NCrystal/NCDump.hh"
#include "NCrystal/NCInfo.hh"
#include <iostream>
#include <cstdlib>
#include <cstdio>

extern "C" 
{
  void generate_bragg_edges(const char* s, int* nbragg, double* data);
}

void generate_bragg_edges(const char* s, int* nbragg, double* data)
{
  double E = 0., wl=0., m_xsectfact=0., m_fdm=0.;
  const NCrystal::Info * info = NCrystal::createInfo(s);
  int i =0, n;
  NCrystal::HKLList::const_iterator itE = info->hklEnd();
  n = info->nHKL();
  *nbragg = n;
  for (NCrystal::HKLList::const_iterator it = info->hklBegin();it!=itE;++it) {
    wl = 2.0*it->dspacing;
    E = NCrystal::wl2ekin(wl);
    m_xsectfact = 0.5/(info->getStructureInfo().volume * info->getStructureInfo().n_atoms);
    m_fdm = it->fsquared_LEAPR * it->multiplicity * it->dspacing;
    data[i] = E;
    data[i+1] = E*m_fdm * m_xsectfact * wl * wl;
    i = i + 2;
  }
}

