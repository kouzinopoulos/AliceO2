/// \file Accumulator.h
/// \brief Definition of the Accumulator class
/// \author Anders Vestbo <mailto:vestbo@fi.uib.no>

#ifndef ALICEO2_HOUGH_ACCUMULATOR_H_
#define ALICEO2_HOUGH_ACCUMULATOR_H_

#include "AliHLTStandardIncludes.h"

namespace AliceO2 {
namespace Hough {

class Accumulator {
public:
  Accumulator();
  Accumulator(int nxbin, double xmin, double xmax, int nybin, double ymin,
              double ymax);
  virtual ~Accumulator();

  void Reset();
  virtual void Fill(double x, double y, int weight = 1);
  virtual void Fill(double x, int ybin, int weight = 1);
  virtual void Fill(int xbin, double y, int weight = 1);
  virtual void Fill(int xbin, int ybin, int weight = 1);
  virtual int FindBin(double x, double y) const;
  virtual int FindLabelBin(double x, double y) const;
  virtual int FindXbin(double x) const;
  virtual int FindYbin(double y) const;
  int GetBin(int xbin, int ybin) const;
  int GetLabelBin(int xbin, int ybin) const;
  int GetBinContent(int bin) const;
  void SetBinContent(int xbin, int ybin, int value);
  void SetBinContent(int bin, int value);
  void AddBinContent(int xbin, int ybin, int weight);
  void AddBinContent(int bin, int weight);
  void Add(Accumulator* h1, double weight = 1);
  void SetThreshold(int i) { fThreshold = i; }

  double GetXmin() const { return fXmin; }
  double GetXmax() const { return fXmax; }
  double GetYmin() const { return fYmin; }
  double GetYmax() const { return fYmax; }
  virtual double GetBinCenterX(int xbin) const;
  virtual double GetBinCenterY(int ybin) const;
  double GetPreciseBinCenterX(float xbin) const;
  double GetPreciseBinCenterY(float ybin) const;
  double GetBinWidthX() const { return fBinwidthX; }
  double GetBinWidthY() const { return fBinwidthY; }
  int GetFirstXbin() const { return fFirstXbin; }
  int GetLastXbin() const { return fLastXbin; }
  int GetFirstYbin() const { return fFirstYbin; }
  int GetLastYbin() const { return fLastYbin; }
  int GetNbinsX() const { return fNxbins; }
  int GetNbinsY() const { return fNybins; }
  int GetNEntries() const { return fEntries; }

  int* fContent; //!
  int* GetContentArray() const { return fContent; }

protected:
  int fNxbins;     // Number of bins in the histogram
  int fNybins;     // Number of bins in the histogram
  int fNcells;     // Overall number of bins in the histogram
  int fEntries;    // Number of entries in the histogram
  int fFirstXbin;  // First active bin
  int fFirstYbin;  // First active bin
  int fLastXbin;   // Last active bin
  int fLastYbin;   // Last active bin
  int fThreshold;  // Bin content threshold

  double fXmin; // Lower limit in X
  double fYmin; // Lower limit in Y
  double fXmax; // Upper limit in X
  double fYmax; // Upper limit in Y

private:
  double fBinwidthX; // Bin width of the Hough space
  double fBinwidthY; // Bin width of the Hough space
};
}
}
#endif
