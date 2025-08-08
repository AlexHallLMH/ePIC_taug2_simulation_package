#include "TauSpinner/ew_born.h"
#include <fstream>
using namespace std;

namespace TauSpinner {

static EWborn::complex ParseComplex(std::ifstream &in)
{
	double x, y;

	in >> x >> y;

	return EWborn::complex(x, y);
}

bool EWborn::FillFromTable(const char *tableLocation)
{
  ifstream in(tableLocation);

  if (!in.is_open())
    return false;

  /* Parse HEADER */
  in >> MZ >> MH >> MT >> SWSQ >> GZ >> MW >> GW;

  if (in.fail())
    return false;

  /* Buffers */
  char buf[16];
  int idx, idx2;

  /* Parse SECTION A */
  for (int n = 0; n < NA; ++n)
  {
    in >> buf >> idx >> EEa[n];

    if (in.fail() || idx != n)
      return false;

    for (int ff = 0; ff < FF_LEN; ++ff)
    {
      FFa[n][ff] = ParseComplex(in);

      if (in.fail())
        return false;
    }

    for (int fs = 0; fs < FS_LEN; ++fs)
    {
      in >> FSa[n][fs];

      if (in.fail())
        return false;
    }
  }

  /* Parse SECTION B */
  for (int n = 0; n < NB; ++n)
  {
    for (int m = 0; m < MB; ++m)
    {
      in >> buf >> idx >> EEb[n] >> idx2;

      if (in.fail() || idx != n || idx2 != m)
        return false;

      for (int ff = 0; ff < FF_LEN; ++ff)
      {
        FFb[n][m][ff] = ParseComplex(in);

        if (in.fail())
          return false;
      }
    }

    for (int fs = 0; fs < FS_LEN; ++fs)
    {
      in >> FSb[n][fs];

      if (in.fail())
        return false;
    }
  }

  /* Parse SECTION C */
  for (int n = 0; n < NC; ++n)
  {
    for (int m = 0; m < MC; ++m)
    {
      in >> buf >> idx >> EEc[n] >> idx2 >> COSc[m];

      if (in.fail() || idx != n || idx2 != m)
        return false;

      for (int ff = 0; ff < FF_LEN; ++ff)
      {
        FFc[n][m][ff] = ParseComplex(in);

        if (in.fail())
          return false;
      }
    }

    for (int fs = 0; fs < FS_LEN; ++fs)
    {
      in >> FSc[n][fs];

      if (in.fail())
        return false;
    }
  }

  /* Parse SECTION D */
  for (int n = 0; n < ND; ++n)
  {
    for (int m = 0; m < MD; ++m)
    {
      in >> buf >> idx >> EEd[n] >> idx2 >> COSd[m];

      if (in.fail() || idx != n || idx2 != m)
        return false;

      for (int ff = 0; ff < FF_LEN; ++ff)
      {
        FFd[n][m][ff] = ParseComplex(in);

        if (in.fail())
          return false;
      }
    }

    for (int fs = 0; fs < FS_LEN; ++fs)
    {
      in >> FSd[n][fs];

      if (in.fail())
        return false;
    }
  }

  return true;
}

} // namespace TauSpinner
