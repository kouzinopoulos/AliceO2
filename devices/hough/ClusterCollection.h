/// \file ClusterCollection.h
/// \brief Definition of the ClusterCollection class
/// \author Charis Kouzinopoulos <charalampos.kouzinopoulos@cern.ch>

#ifndef ALICEO2_HOUGH_ClusterCollection_H_
#define ALICEO2_HOUGH_ClusterCollection_H_

#include "StandardIncludes.h"

#include "AliHLTTPCTrackGeometry.h"
#include "AliHLTTPCClusterDataFormat.h"
#include "AliHLTTPCSpacePointContainer.h"
#include "AliHLTComponent.h"
#include "AliHLTTPCDefinitions.h"

namespace AliceO2 {
namespace Hough {

struct clusterDataFormat {
  UInt_t mID;
  UInt_t mCharge;

  // The ALICE TPC detector consists out of 2 x 18 = 36 slices (or sectors?) (18 slices at the top, 18 slices at the bottom). Each
  // slice consists of 2 inner patches (or partitions) and 4 outer patches = 6 patches for a total of 216 patches.
  // Moreover, there are 15488 (5504+9984) pads and 159 (63+64+32) pad rows per sector for a total of 557568 (36x15488) pads.
  UInt_t mTPCSlice;
  UInt_t mTPCPartition;

  // There are 3 different systems: Raw: row, pad, time Local : x,y and global z Global: global x,y and global z
  UChar_t mPadRow;
  Double_t mPad;
  Double_t mTime;

  // SpacePointData retrieved data in local coordinates
  Double_t mX;
  Double_t mY;
  Double_t mZ;

  Double_t mAlpha;
  Double_t mBeta;
  Double_t mEta;

  UInt_t mEtaSlice;
};

class ClusterCollection {
public:
  ClusterCollection();
  virtual ~ClusterCollection();

  /// Load cluster information to memory from a directory
  UInt_t processData(std::string dataPath, std::string dataType, std::string dataOrigin);

  UInt_t getClusterID(int clusterNumber) { return clusterData[clusterNumber].mID; }
  UInt_t getClusterCharge(int clusterNumber) { return clusterData[clusterNumber].mCharge; }

  UInt_t getClusterSlice(int clusterNumber) { return clusterData[clusterNumber].mTPCSlice; }
  UInt_t getClusterPartition(int clusterNumber) { return clusterData[clusterNumber].mTPCPartition; }

  UChar_t getClusterPadRow(int clusterNumber) { return clusterData[clusterNumber].mPadRow; }
  Double_t getClusterPad(int clusterNumber) { return clusterData[clusterNumber].mPad; }
  Double_t getClusterTime(int clusterNumber) { return clusterData[clusterNumber].mTime; }

  Double_t getClusterX(int clusterNumber) { return clusterData[clusterNumber].mX; }
  Double_t getClusterY(int clusterNumber) { return clusterData[clusterNumber].mY; }
  Double_t getClusterZ(int clusterNumber) { return clusterData[clusterNumber].mZ; }

  Double_t getClusterAlpha(int clusterNumber) { return clusterData[clusterNumber].mAlpha; }
  Double_t getClusterBeta(int clusterNumber) { return clusterData[clusterNumber].mBeta; }
  Double_t getClusterEta(int clusterNumber) { return clusterData[clusterNumber].mEta; }
  Int_t getClusterEtaSlice(int clusterNumber) { return clusterData[clusterNumber].mEtaSlice; }

  UInt_t getNumberOfClustersPerPadRow(UInt_t padRow) { return clusterData2[padRow].size(); }

  void setClusterAlpha(int clusterNumber, Double_t alpha) { clusterData[clusterNumber].mAlpha = alpha; }
  void setClusterBeta(int clusterNumber, Double_t beta) { clusterData[clusterNumber].mBeta = beta; }
  void setClusterEta(int clusterNumber, Double_t eta) { clusterData[clusterNumber].mEta = eta; }
  void setClusterEtaSlice(int clusterNumber, UInt_t etaSlice) { clusterData[clusterNumber].mEtaSlice = etaSlice; }

  void setClusterParameters(UInt_t clusterID, Double_t x, Double_t y, Double_t z, UInt_t charge, UChar_t padRow,
                            Double_t pad, Double_t time, UInt_t tpcSlice, UInt_t tpcPartition)
  {
    clusterDataFormat data;
    data.mID = clusterID;
    data.mCharge = charge;
    data.mTPCSlice = tpcSlice;
    data.mTPCPartition = tpcPartition;
    data.mPadRow = padRow;
    data.mPad = pad;
    data.mTime = time;
    data.mX = x;
    data.mY = y;
    data.mZ = z;

    clusterData.push_back(data);
  }

  UInt_t getClusterCharge(UInt_t padRow, UInt_t clusterNumber) { return clusterData2[padRow][clusterNumber].mCharge; }

  Double_t getClusterPad(UInt_t padRow, UInt_t clusterNumber) { return clusterData2[padRow][clusterNumber].mPad; }
  Double_t getClusterTime(UInt_t padRow, UInt_t clusterNumber) { return clusterData2[padRow][clusterNumber].mTime; }

  Double_t getClusterX(UInt_t padRow, UInt_t clusterNumber) { return clusterData2[padRow][clusterNumber].mX; }
  Double_t getClusterY(UInt_t padRow, UInt_t clusterNumber) { return clusterData2[padRow][clusterNumber].mY; }
  Double_t getClusterZ(UInt_t padRow, UInt_t clusterNumber) { return clusterData2[padRow][clusterNumber].mZ; }

  Double_t getClusterAlpha(UInt_t padRow, UInt_t clusterNumber) { return clusterData2[padRow][clusterNumber].mAlpha; }
  Double_t getClusterBeta(UInt_t padRow, UInt_t clusterNumber) { return clusterData2[padRow][clusterNumber].mBeta; }

  void setClusterAlpha(UInt_t padRow, UInt_t clusterNumber, Double_t alpha) { clusterData2[padRow][clusterNumber].mAlpha = alpha; }
  void setClusterBeta(UInt_t padRow, UInt_t clusterNumber, Double_t beta) { clusterData2[padRow][clusterNumber].mBeta = beta; }

  void setClusterParameters2(UInt_t clusterID, Double_t x, Double_t y, Double_t z, UInt_t charge, UChar_t padRow,
                             Double_t pad, Double_t time, UInt_t tpcSlice, UInt_t tpcPartition)
  {
    clusterDataFormat data;
    data.mID = clusterID;
    data.mCharge = charge;
    data.mTPCSlice = tpcSlice;
    data.mTPCPartition = tpcPartition;
    data.mPadRow = padRow;
    data.mPad = pad;
    data.mTime = time;
    data.mX = x;
    data.mY = y;
    data.mZ = z;

    clusterData2[(UInt_t)padRow].push_back(data);
  }

protected:
private:
  // Information per cluster
  std::vector<clusterDataFormat> clusterData;

  // FIXME: replace clusterData with clusterData2
  std::vector<std::vector<clusterDataFormat>> clusterData2;

  std::unique_ptr<AliHLTTPCSpacePointContainer> spacepoints;
};
}
}
#endif
