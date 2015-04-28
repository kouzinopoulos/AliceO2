/// \file runHough.h
/// \brief Definition of a cluster loader
/// \author Charis Kouzinopoulos

#include "AliHLTTPCTrackGeometry.h"
#include "AliHLTTPCClusterDataFormat.h"
#include "AliHLTTPCSpacePointContainer.h"
#include "AliHLTComponent.h"
#include "AliHLTTPCDefinitions.h"

#include "TCanvas.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TGraph2D.h"

#include "boost/filesystem.hpp"

#include "Accumulator.h"
#include "HoughTransformerRow.h"

// FIXME: convert float into double for bigger precision
// Information for each stored cluster
struct clusterDataFormat {
  UInt_t mID;
  UInt_t mCharge;

  // The ALICE TPC detector consists out of 2 x 18 = 36 slices (18 slices at the top, 18 slices at the bottom). Each
  // slices consists out of 2 inner patches (or partitions) and 4 outer patches = 6 patches. So we have a total of 216
  // patches.
  UInt_t mTPCSlice;
  UInt_t mTPCPartition;

  Double_t mX;
  Double_t mY;
  Double_t mZ;

  Double_t mAlpha;
  Double_t mBeta;
  Double_t mEta;

  UInt_t mEtaSlice;
};

// Information for each reconstructed track
struct trackDataFormat {
  Double_t mAlpha1;
  Double_t mBeta1;
  Double_t mAlpha2;
  Double_t mBeta2;

  UInt_t mEtaSlice;
};

struct accumulatorDataFormat {
  // UShort_t bin[rResolution * thetalphaMax + thetalphaMax];
  UShort_t* bin;
};

struct etaRow {
  UChar_t fStartPad; // First pad in the cluster
  UChar_t fEndPad;   // Last pad in the cluster
  Bool_t fIsFound;   // Is the cluster already found
};

// Parameters which represent a given pad in the hough space. Used in order to avoid as much as possible floating point
// operations during the hough transform
struct houghPadParameters {
  Float_t fAlpha;      // Starting value for the hough parameter alpha1
  Float_t fDeltaAlpha; // Slope of alpha1
  Int_t fFirstBin;     // First alpha2 bin to be filled
  Int_t fLastBin;      // Last alpha2 bin to be filled
};

std::unique_ptr<AliHLTTPCSpacePointContainer> spacepoints;
std::vector<clusterDataFormat> clusterData;
std::vector<trackDataFormat> trackData;
std::vector<accumulatorDataFormat> accumulatorData;

std::vector<unsigned int> accumulator;

// TPC specification
Int_t numberOfPatches = 6;
Int_t padRowIndex[6][2] = { { 0, 29 }, { 30, 62 }, { 63, 90 }, { 91, 116 }, { 117, 139 }, { 140, 158 } };
Int_t numberOfPadRows[6] = { 30, 33, 28, 26, 23, 19 };
Int_t numberOfRows = 159;

AliceO2::Hough::Accumulator** parameterSpace;

Double_t fgX[159] = { 85.195,  85.945,  86.695,  87.445,  88.195,  88.945,  89.695,  90.445,  91.195,  91.945,  92.695,
                      93.445,  94.195,  94.945,  95.695,  96.445,  97.195,  97.945,  98.695,  99.445,  100.195, 100.945,
                      101.695, 102.445, 103.195, 103.945, 104.695, 105.445, 106.195, 106.945, 107.695, 108.445, 109.195,
                      109.945, 110.695, 111.445, 112.195, 112.945, 113.695, 114.445, 115.195, 115.945, 116.695, 117.445,
                      118.195, 118.945, 119.695, 120.445, 121.195, 121.945, 122.695, 123.445, 124.195, 124.945, 125.695,
                      126.445, 127.195, 127.945, 128.695, 129.445, 130.195, 130.945, 131.695, 135.180, 136.180, 137.180,
                      138.180, 139.180, 140.180, 141.180, 142.180, 143.180, 144.180, 145.180, 146.180, 147.180, 148.180,
                      149.180, 150.180, 151.180, 152.180, 153.180, 154.180, 155.180, 156.180, 157.180, 158.180, 159.180,
                      160.180, 161.180, 162.180, 163.180, 164.180, 165.180, 166.180, 167.180, 168.180, 169.180, 170.180,
                      171.180, 172.180, 173.180, 174.180, 175.180, 176.180, 177.180, 178.180, 179.180, 180.180, 181.180,
                      182.180, 183.180, 184.180, 185.180, 186.180, 187.180, 188.180, 189.180, 190.180, 191.180, 192.180,
                      193.180, 194.180, 195.180, 196.180, 197.180, 198.180, 199.430, 200.930, 202.430, 203.930, 205.430,
                      206.930, 208.430, 209.930, 211.430, 212.930, 214.430, 215.930, 217.430, 218.930, 220.430, 221.930,
                      223.430, 224.930, 226.430, 227.930, 229.430, 230.930, 232.430, 233.930, 235.430, 236.930, 238.430,
                      239.930, 241.430, 242.930, 244.430, 245.930 };

Int_t numberOfPads[159] = { 67,  67,  69,  69,  69,  71,  71,  71,  73,  73,  73,  75,  75,  75,  77,  77,  77,  79,
                            79,  79,  81,  81,  81,  83,  83,  83,  85,  85,  85,  87,  87,  87,  89,  89,  89,  91,
                            91,  91,  93,  93,  93,  95,  95,  95,  97,  97,  97,  99,  99,  99,  99,  101, 101, 101,
                            103, 103, 103, 105, 105, 105, 107, 107, 107, 73,  75,  75,  75,  75,  77,  77,  77,  79,
                            79,  79,  81,  81,  81,  81,  83,  83,  83,  85,  85,  85,  85,  87,  87,  87,  89,  89,
                            89,  91,  91,  91,  91,  93,  93,  93,  95,  95,  95,  95,  97,  97,  97,  99,  99,  99,
                            101, 101, 101, 101, 103, 103, 103, 105, 105, 105, 105, 107, 107, 107, 109, 109, 109, 111,
                            111, 111, 113, 113, 113, 115, 115, 117, 117, 119, 119, 121, 121, 121, 123, 123, 125, 125,
                            127, 127, 127, 129, 129, 131, 131, 133, 133, 135, 135, 135, 137, 137, 139 };

Double_t Pi() { return 3.141592653589793; }

Double_t xMin;
Double_t xMax;
Double_t yMin;
Double_t yMax;

Double_t alphaMin;
Double_t alphaMax;
Double_t betaMin;
Double_t betaMax;

float etalphaMin;
float etalphaMax;

int rMax;
int thetalphaMax;

// Project parameters
int houghThreshold = 53;
int etaResolution = 20;
int thetaResolution = 180;

Double_t padPitchWidthLow = 0.4;
Double_t padPitchWidthUp = 0.6;

// Hough transform parameters
UChar_t** fGapCount;
UChar_t** fCurrentRowCount;
UChar_t** fPreviousBin;
UChar_t** fNextBin;
UChar_t** fNextRow;

UChar_t* fTrackLastRow;

houghPadParameters** fStartPadParameters;
houghPadParameters** fEndPadParameters;

// By setting rResolution = 1, r = 18.6 and r = 19.2 for theta = 48 will both cause bin(19,48) to increase. With a
// resoution of rResolution = 10, they will cause bin(186,48) and bin(192,48) to increase respectively. Now rResolution
// = 100 means that there will be 100 bins for r: 0 - 99. Converting r -> rBin and rBin -> r is done by methods
// getRValue and getRBinValue
int rResolution = 10000;

// sin and cos expect values in radians instead of degrees
#define DEG2RAD 0.017453293f

//#FIXME: Get rid of the etaSlice paramater in the accumulator. We keep a single accumulatr per etaSlice
// Create instead a vector fParamSpace of histograms for each etaIndex (AliHLTHoughTransformRow::GetHistogram)
void setAccumulatorBin(int etaSlice, int r, int theta)
{
  accumulator[etaSlice * (rResolution * thetalphaMax) + r * thetalphaMax + theta]++;
}

int getAccumulatorBin(int etaSlice, int r, int theta)
{
  return accumulator[etaSlice * (rResolution * thetalphaMax) + r * thetalphaMax + theta];
}

// Solve getRBinValue to r
float getRValue(int rBin)
{
  int rMin = -rMax;
  return (float)((rMax - rMin) * rBin + rMin * rResolution) / rResolution;
}

int getRBinValue(float r)
{
  int rMin = -rMax;
  return ((r - rMin) * rResolution) / (rMax - rMin);
}

UInt_t getTrackEtaSlice(int trackNumber) { return trackData[trackNumber].mEtaSlice; }

Float_t getTrackAlpha1(int trackNumber) { return trackData[trackNumber].mAlpha1; }
Float_t getTrackBeta1(int trackNumber) { return trackData[trackNumber].mBeta1; }
Float_t getTrackAlpha2(int trackNumber) { return trackData[trackNumber].mAlpha2; }
Float_t getTrackBeta2(int trackNumber) { return trackData[trackNumber].mBeta2; }

UInt_t getNumberOfTracks() { return trackData.size(); }

void setTrackParameters(UInt_t etaSlice, Float_t alpha1, Float_t beta1, Float_t alpha2, Float_t beta2)
{
  trackDataFormat track;
  track.mEtaSlice = etaSlice;
  track.mAlpha1 = alpha1;
  track.mBeta1 = beta1;
  track.mAlpha2 = alpha2;
  track.mBeta2 = beta2;
  trackData.push_back(track);
}

UInt_t getClusterID(int clusterNumber) { return clusterData[clusterNumber].mID; }
UInt_t getClusterSlice(int clusterNumber) { return clusterData[clusterNumber].mTPCSlice; }
UInt_t getClusterPartition(int clusterNumber) { return clusterData[clusterNumber].mTPCPartition; }

Double_t getClusterX(int clusterNumber) { return clusterData[clusterNumber].mX; }
Double_t getClusterY(int clusterNumber) { return clusterData[clusterNumber].mY; }
Double_t getClusterZ(int clusterNumber) { return clusterData[clusterNumber].mZ; }

UInt_t getClusterCharge(int clusterNumber) { return clusterData[clusterNumber].mCharge; }

Double_t getClusterAlpha(int clusterNumber) { return clusterData[clusterNumber].mAlpha; }
Double_t getClusterBeta(int clusterNumber) { return clusterData[clusterNumber].mBeta; }
Double_t getClusterEta(int clusterNumber) { return clusterData[clusterNumber].mEta; }
Int_t getClusterEtaSlice(int clusterNumber) { return clusterData[clusterNumber].mEtaSlice; }

void setClusterAlpha(int clusterNumber, Double_t alpha) { clusterData[clusterNumber].mAlpha = alpha; }
void setClusterBeta(int clusterNumber, Double_t beta) { clusterData[clusterNumber].mBeta = beta; }
void setClusterEta(int clusterNumber, Double_t eta) { clusterData[clusterNumber].mEta = eta; }
void setClusterEtaSlice(int clusterNumber, UInt_t etaSlice) { clusterData[clusterNumber].mEtaSlice = etaSlice; }

void setClusterParameters(UInt_t clusterID, Double_t x, Double_t y, Double_t z, UInt_t charge, UInt_t tpcSlice,
                          UInt_t tpcPartition)
{
  clusterDataFormat data;
  data.mID = clusterID;
  data.mX = x;
  data.mY = y;
  data.mZ = z;
  data.mCharge = charge;
  data.mTPCSlice = tpcSlice;
  data.mTPCPartition = tpcPartition;
  clusterData.push_back(data);
}

/// Returns the first row per partition
Int_t getFirstPadRow(Int_t partition)
{
  if (partition == -1) {
    return 0;
  } else if (partition < -1 || partition >= 6) {
    std::cerr << "Wrong partition: " << partition << endl;
    std::exit(1);
  } else {
    return padRowIndex[partition][0];
  }
}

/// Returns the last row per partition
Int_t getLastPadRow(Int_t partition)
{
  if (partition == -1) {
    return padRowIndex[5][1];
  } else if (partition < -1 || partition >= 6) {
    std::cerr << "Wrong partition: " << partition << endl;
    std::exit(1);
  } else {
    return padRowIndex[partition][1];
  }
}

/// Returns the number of pads per row
Int_t getNumberOfPads(Int_t row)
{
  if (row < 0 || row >= numberOfRows) {
    std::cerr << "Wrong row " << row << endl;
    return 0;
  }
  return numberOfPads[row];
}

/// Returns the patch for a given padrow
Int_t getPatch(Int_t padrow)
{
  if (padrow < 0 || padrow >= numberOfRows) {
    std::cerr << "Wrong padrow " << padrow << endl;
    return -2;
  }
  Int_t patch = 0;
  while (patch < numberOfPatches) {
    if (padrow >= padRowIndex[patch][0] && padrow <= padRowIndex[patch][1])
      break;
    patch++;
  }
  return patch;
}

/// Returns the pad patch width for a given patch
Double_t getPadPitchWidth(Int_t patch)
{

  if (patch < 0 || patch > numberOfPatches) {
    std::cerr << "Wrong patch " << patch << endl;
    return -1;
  }
  return patch < 2 ? padPitchWidthLow : padPitchWidthUp;
}

/// Slicerow to X value (slice 0)
Double_t Row2X(Int_t slicerow)
{
  if (slicerow < 0 || slicerow >= numberOfRows) {
    std::cerr << "Wrong slicerow " << slicerow << endl;
    return 0;
  }
  return fgX[slicerow];
}

/// Returns the corresponding histogram bin
/*Int_t getHistogramBin(Int_t xbin,Int_t ybin) const
{
  if(xbin < fFirstXbin || xbin > fLastXbin)
    return 0;
  if(ybin < fFirstYbin || ybin > fLastYbin)
    return 0;

  return xbin + ybin*(fNxbins+2);
}
*/
