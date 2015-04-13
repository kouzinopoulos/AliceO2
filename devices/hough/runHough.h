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

#include <sstream>

// FIXME: convert float into double for bigger precision
struct clusterDataFormat {
  UInt_t mID;
  UInt_t mCharge;

  Float_t mX;
  Float_t mY;
  Float_t mZ;

  Float_t mAlpha;
  Float_t mBeta;
  Float_t mEta;

  UInt_t mEtaSlice;
};

struct trackDataFormat {
  Float_t mAlpha1;
  Float_t mBeta1;
  Float_t mAlpha2;
  Float_t mBeta2;

  UInt_t mEtaSlice;
};

std::unique_ptr<AliHLTTPCSpacePointContainer> spacepoints;
std::vector<clusterDataFormat> clusterData;
std::vector<trackDataFormat> trackData;

std::vector<unsigned int> accumulator;

int xMax = 0;
int yMax = 0;
int zMax = 0;
int aMin = 0;
int aMax = 0;
int bMin = 0;
int bMax = 0;

float etaMin;
float etaMax;

int rMax;
int thetaMax;

// project parameters
int houghThreshold = 53;
int etaResolution = 20;

// By setting rResolution = 1, r = 18.6 and r = 19.2 for theta = 48 will both
// cause bin(19,48) to increase. With a
// resoution of rResolution = 10, they will cause bin(186,48) and bin(192,48) to
// increase respectively.
// Now rResolution = 100 means that there will be 100 bins for r: 0 - 99.
// Converting r -> rBin and rBin -> r is done by
// methods getRValue and getRBinValue
int rResolution = 100000;

// sin and cos expect values in radians instead of degrees
#define DEG2RAD 0.017453293f

void setAccumulatorBin(int etaSlice, int r, int theta)
{
  accumulator[etaSlice * (rResolution * thetaMax) + r * thetaMax + theta]++;
}

int getAccumulatorBin(int etaSlice, int r, int theta)
{
  return accumulator[etaSlice * (rResolution * thetaMax) + r * thetaMax + theta];
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

Float_t getClusterX(int clusterNumber) { return clusterData[clusterNumber].mX; }

Float_t getClusterY(int clusterNumber) { return clusterData[clusterNumber].mY; }

Float_t getClusterZ(int clusterNumber) { return clusterData[clusterNumber].mZ; }

UInt_t getClusterCharge(int clusterNumber) { return clusterData[clusterNumber].mCharge; }

Float_t getClusterAlpha(int clusterNumber) { return clusterData[clusterNumber].mAlpha; }

Float_t getClusterBeta(int clusterNumber) { return clusterData[clusterNumber].mBeta; }

Float_t getClusterEta(int clusterNumber) { return clusterData[clusterNumber].mEta; }

Int_t getClusterEtaSlice(int clusterNumber) { return clusterData[clusterNumber].mEtaSlice; }

void setClusterAlpha(int clusterNumber, Float_t alpha) { clusterData[clusterNumber].mAlpha = alpha; }

void setClusterBeta(int clusterNumber, Float_t beta) { clusterData[clusterNumber].mBeta = beta; }

void setClusterEta(int clusterNumber, Float_t eta) { clusterData[clusterNumber].mEta = eta; }

void setClusterEtaSlice(int clusterNumber, UInt_t etaSlice) { clusterData[clusterNumber].mEtaSlice = etaSlice; }

void setClusterParameters(UInt_t clusterID, Float_t x, Float_t y, Float_t z, UInt_t charge)
{
  clusterDataFormat data;
  data.mID = clusterID;
  data.mX = x;
  data.mY = y;
  data.mZ = z;
  data.mCharge = charge;
  clusterData.push_back(data);
}
