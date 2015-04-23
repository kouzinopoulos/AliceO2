/// \file Accumulator.cxx
/// \brief Implementation of the Accumulator class
/// \author Anders Vestbo <mailto:vestbo@fi.uib.no>

#include "Accumulator.h"

// uncomment if you want overflow checks
//#define _IFON_

using namespace std;
using namespace AliceO2::Hough;

Accumulator::Accumulator()
{
  // Default constructor
  fNxbins = 0;
  fNybins = 0;
  fNcells = 0;
  fXmin = 0;
  fYmin = 0;
  fXmax = 0;
  fYmax = 0;
  fBinwidthX = 0;
  fBinwidthY = 0;
  fFirstXbin = 0;
  fLastXbin = 0;
  fFirstYbin = 0;
  fLastYbin = 0;
  fEntries = 0;
  fContent = 0;
  fThreshold = 0;
}

Accumulator::Accumulator(int nxbin, double xmin, double xmax, int nybin,
                         double ymin, double ymax)
{
  fNxbins = nxbin;
  fNybins = nybin;
  fNcells = (nxbin + 2) * (nybin + 2);
  fXmin = xmin;
  fYmin = ymin;
  fXmax = xmax;
  fYmax = ymax;
  fBinwidthX = (fXmax - fXmin) / fNxbins;
  fBinwidthY = (fYmax - fYmin) / fNybins;

  fEntries = 0;
  fFirstXbin = 1;
  fFirstYbin = 1;
  fLastXbin = nxbin;
  fLastYbin = nybin;
  fThreshold = 0;

  fContent = new int[fNcells];
  Reset();
}

Accumulator::~Accumulator()
{
  // Destructor
  if (fContent)
    delete[] fContent;
}

void Accumulator::Reset()
{
  // Reset histogram contents
  if (fContent)
    for (int i = 0; i < fNcells; i++)
      fContent[i] = 0;

  fEntries = 0;
}

void Accumulator::Fill(double x, double y, int weight)
{
  // Fill the weight into a bin which correspond to x and y
  int bin = FindBin(x, y);
#ifdef _IFON_
  if (bin < 0)
    return;
#endif

  AddBinContent(bin, weight);
}

void Accumulator::Fill(double x, int ybin, int weight)
{
  // Fill the weight into a bin which correspond to x and ybin
  int xbin = FindXbin(x);
  int bin = GetBin(xbin, ybin);
#ifdef _IFON_
  if (bin < 0)
    return;
#endif

  AddBinContent(bin, weight);
}

void Accumulator::Fill(int xbin, double y, int weight)
{
  // Fill the weight into a bin which correspond to xbin and y
  int ybin = FindYbin(y);
  int bin = GetBin(xbin, ybin);
#ifdef _IFON_
  if (bin < 0)
    return;
#endif

  AddBinContent(bin, weight);
}

void Accumulator::Fill(int xbin, int ybin, int weight)
{
  // Fill the weight into a bin which correspond to xbin and ybin
  int bin = GetBin(xbin, ybin);
#ifdef _IFON_
  if (bin < 0)
    return;
#endif

  AddBinContent(bin, weight);
}

int Accumulator::FindBin(double x, double y) const
{
  // Finds the bin which correspond to x and y
  int xbin = FindXbin(x);
  int ybin = FindYbin(y);
#ifdef _IFON_
  if (!xbin || !ybin)
    return -1;
#endif

  return GetBin(xbin, ybin);
}

int Accumulator::FindLabelBin(double x, double y) const
{
  // Returns the corresponding bin with the mc labels
  int xbin = FindXbin(x);
  int ybin = FindYbin(y);
#ifdef _IFON_
  if (!xbin || !ybin)
    return -1;
#endif

  return GetLabelBin(xbin, ybin);
}

int Accumulator::FindXbin(double x) const
{
  // Finds the bin which correspond to x
  if (x < fXmin || x > fXmax)
    return 0;

  return 1 + (int)(fNxbins * (x - fXmin) / (fXmax - fXmin));
}

int Accumulator::FindYbin(double y) const
{
  // Finds the bin which correspond to y
  if (y < fYmin || y > fYmax)
    return 0;

  return 1 + (int)(fNybins * (y - fYmin) / (fYmax - fYmin));
}

int Accumulator::GetBin(int xbin, int ybin) const
{
  // Returns the bin which correspond to xbin and ybin
  if (xbin < fFirstXbin || xbin > fLastXbin)
    return 0;
  if (ybin < fFirstYbin || ybin > fLastYbin)
    return 0;

  return xbin + ybin * (fNxbins + 2);
}

int Accumulator::GetLabelBin(int xbin, int ybin) const
{
  // Returns the corresponding bin with the mc labels
  if (xbin < fFirstXbin || xbin > fLastXbin)
    return -1;
  if (ybin < fFirstYbin || ybin > fLastYbin)
    return -1;

  return (int)(xbin / 2) + ((int)(ybin / 2)) * ((int)((fNxbins + 3) / 2));
}

int Accumulator::GetBinContent(int bin) const
{
  // Return the bin content
  if (bin >= fNcells) {
    std::cerr << "bin out of range " << bin << endl;
    return 0;
  }

  if (fContent[bin] < fThreshold)
    return 0;
  return fContent[bin];
}

void Accumulator::SetBinContent(int xbin, int ybin, int value)
{
  // Set bin content
  int bin = GetBin(xbin, ybin);
#ifdef _IFON_
  if (bin == 0)
    return;
#endif

  SetBinContent(bin, value);
}

void Accumulator::SetBinContent(int bin, int value)
{
  // Set bin content

  if (bin >= fNcells) {
    std::cerr << "bin out of range " << bin << endl;
    return;
  }

  if (bin == 0)
    return;
  fContent[bin] = value;
}

void Accumulator::AddBinContent(int xbin, int ybin, int weight)
{
  // Adds weight to bin content
  int bin = GetBin(xbin, ybin);
#ifdef _IFON_
  if (bin == 0)
    return;
#endif

  AddBinContent(bin, weight);
}

void Accumulator::AddBinContent(int bin, int weight)
{
  // Adds weight to bin content
  if (bin < 0 || bin > fNcells) {
    std::cerr << "bin-value out of range " << bin << endl;
    return;
  }
  if (bin == 0)
    return;
  fEntries++;
  fContent[bin] += weight;
}

void Accumulator::Add(Accumulator* h1, double /*weight*/)
{
  // Adding two histograms. Should be identical.

  if (!h1) {
    std::cerr << "Attempting to add a non-existing histogram" << endl;
    return;
  }

  if (h1->GetNbinsX() != fNxbins || h1->GetNbinsY() != fNybins) {
    std::cerr << "Mismatch in the number of bins " << endl;
    return;
  }

  if (h1->GetFirstXbin() != fFirstXbin || h1->GetLastXbin() != fLastXbin || h1->GetFirstYbin() != fFirstYbin ||
      h1->GetLastYbin() != fLastYbin) {
    std::cerr << "Mismatch in the bin numbering " << endl;
    return;
  }

  for (int bin = 0; bin < fNcells; bin++)
    fContent[bin] += h1->GetBinContent(bin);

  fEntries += h1->GetNEntries();
}

double Accumulator::GetBinCenterX(int xbin) const
{
  // Returns the position of the center of a bin
  if (xbin < fFirstXbin || xbin > fLastXbin) {
    std::cerr << "Bin-value out of range " << xbin << endl;
    return -1;
  }

  return fXmin + (xbin - 0.5) * fBinwidthX;
}

double Accumulator::GetBinCenterY(int ybin) const
{
  // Returns the position of the center of a bin
  if (ybin < fFirstYbin || ybin > fLastYbin) {
    std::cerr << "Bin-value out of range " << ybin << endl;
    return -1;
  }

  return fYmin + (ybin - 0.5) * fBinwidthY;
}

double Accumulator::GetPreciseBinCenterX(float xbin) const
{
  // Returns the position of the center of a bin using precise values inside the bin
  if (xbin < (fFirstXbin - 1.5) || xbin > (fLastXbin + 1.5)) {
    std::cerr << "Bin-value out of range " << xbin << endl;
    return -1;
  }
  //  return fXmin + (xbin-1) * fBinwidthX + 0.5*fBinwidthX;
  return fXmin + (xbin - 0.5) * fBinwidthX;
}

double Accumulator::GetPreciseBinCenterY(float ybin) const
{
  // Returns the position of the center of a bin using precise values inside the bin
  if (ybin < (fFirstYbin - 1.5) || ybin > (fLastYbin + 1.5)) {
    std::cerr << "Bin-value out of range " << ybin << endl;
    return -1;
  }
  //  return fYmin + (ybin-1) * fBinwidthY + 0.5*fBinwidthY;
  return fYmin + (ybin - 0.5) * fBinwidthY;
}
