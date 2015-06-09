#ifndef ALICEO2_HOUGH_HOUGHTRACK_H_
#define ALICEO2_HOUGH_HOUGHTRACK_H_

#include "StandardIncludes.h"
#include "TransformerRow.h"
#include "Track.h"
#include "Transform.h"

namespace AliceO2 {
namespace Hough {

class HoughTrack : public Track {

public:
  HoughTrack();
  virtual ~HoughTrack();

  virtual void Set(Track* track);
  virtual int Compare(const Track* track) const;

  bool IsHelix() const { return fIsHelix; }
  void UpdateToFirstRow();
  void SetTrackParameters(double kappa, double eangle, int weight);
  void SetTrackParametersRow(double alpha1, double alpha2, double eta, int weight);
  void SetLineParameters(double psi, double D, int weight, int* rowrange, int refrow);

  int GetWeight() const { return fWeight; }
  double GetPsiLine() const { return fPsiLine; }
  double GetDLine() const { return fDLine; }

  int GetEtaIndex() const { return fEtaIndex; }
  double GetEta() const { return fEta; }
  int GetSlice() const { return fSlice; }
  void GetLineCrossingPoint(int padrow, float* xy);

  float GetBinX() const { return fBinX; }
  float GetBinY() const { return fBinY; }
  float GetSizeX() const { return fSizeX; }
  float GetSizeY() const { return fSizeY; }

  void SetHelixTrue() { fIsHelix = true; }
  void SetSlice(int slice) { fSlice = slice; }
  void SetEta(double f);
  void SetWeight(int i, bool update = false)
  {
    if (update) {
      fWeight += i;
    }
    else {
      fWeight = i;
    }
  }
  void SetEtaIndex(int f) { fEtaIndex = f; }
  void SetBestMCid(int f, double mindist);
  void SetDLine(double f) { fDLine = f; }
  void SetPsiLine(double f) { fPsiLine = f; }

  void SetBinXY(float binx, float biny, float sizex, float sizey)
  {
    fBinX = binx;
    fBinY = biny;
    fSizeX = sizex;
    fSizeY = sizey;
  }

private:
  double fMinDist; ///< Minimum distance to a generated track while associating mc label
  int fWeight;     ///< Track weight
  int fEtaIndex;   ///< Eta slice index
  double fEta;     ///< Track Eta
  int fSlice;      ///< The slice where this track was found

  double fDLine;   ///< ??
  double fPsiLine; ///< ??

  bool fIsHelix; ///< Is the track helix or not?

  float fBinX, fBinY, fSizeX, fSizeY; ///< Size and position of the hough space peak
};

}
}

#endif
