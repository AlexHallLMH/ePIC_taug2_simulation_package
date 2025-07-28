#ifndef _EW_BORN_H_
#define _EW_BORN_H_
#include <complex>
using std::complex;

namespace TauSpinner {

struct EWborn
{
  typedef std::complex<double> complex;

  static const int FF_LEN = 7;
  static const int FS_LEN = 4;

  /* HEADER */
  double MZ;
  double MH;
  double MT;
  double SWSQ;
  double GZ;
  double MW;
  double GW;

  /* SECTION A */
  static const int NA = 101;

  double  EEa[NA];
  complex FFa[NA][FF_LEN];
  double  FSa[NA][FS_LEN];

  /* SECTION B */
  static const int NB = 121;
  static const int MB = 15;

  double  EEb[NB];
  complex FFb[NB][MB][FF_LEN];
  double  FSb[NB][FS_LEN];

  /* SECTION C */
  static const int NC = 146;
  static const int MC = 31;

  double  EEc[NC];
  double  COSc[MC];
  complex FFc[NC][MC][FF_LEN];
  double  FSc[NC][FS_LEN];

  /* SECTION D */
  static const int ND = 81;
  static const int MD = 15;

  double  EEd[ND];
  double  COSd[MD];
  complex FFd[ND][MD][FF_LEN];
  double  FSd[ND][FS_LEN];

  /* Functions */
  bool FillFromTable(const char *tableLocation);
};
} // namespace TauSpinner

#endif

