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

  // Estimation of the cluster's pseudorapidity as calculated by the calculateEta and calculateEtaSlice methods
  Double_t mEta;
  UInt_t mEtaSlice;
};

class ClusterCollection {
public:
  ClusterCollection();
  virtual ~ClusterCollection();

  /// Load cluster information to memory from disk
  UInt_t readData(std::string dataPath, std::string dataType, std::string dataOrigin);

  UInt_t getNumberOfClustersPerPadRow(UInt_t padRow) { return clusterData[padRow].size(); }

  UInt_t getClusterID(UInt_t padRow, int clusterNumber) { return clusterData[padRow][clusterNumber].mID; }
  UInt_t getClusterCharge(UInt_t padRow, UInt_t clusterNumber) { return clusterData[padRow][clusterNumber].mCharge; }

  UInt_t getClusterSlice(UInt_t padRow, int clusterNumber) { return clusterData[padRow][clusterNumber].mTPCSlice; }
  UInt_t getClusterPartition(UInt_t padRow, int clusterNumber) { return clusterData[padRow][clusterNumber].mTPCPartition; }

  Double_t getClusterPad(UInt_t padRow, UInt_t clusterNumber) { return clusterData[padRow][clusterNumber].mPad; }
  Double_t getClusterTime(UInt_t padRow, UInt_t clusterNumber) { return clusterData[padRow][clusterNumber].mTime; }

  Double_t getClusterX(UInt_t padRow, UInt_t clusterNumber) { return clusterData[padRow][clusterNumber].mX; }
  Double_t getClusterY(UInt_t padRow, UInt_t clusterNumber) { return clusterData[padRow][clusterNumber].mY; }
  Double_t getClusterZ(UInt_t padRow, UInt_t clusterNumber) { return clusterData[padRow][clusterNumber].mZ; }

  Double_t getClusterEta(UInt_t padRow, UInt_t clusterNumber) { return clusterData[padRow][clusterNumber].mEta; }
  UInt_t getClusterEtaSlice(UInt_t padRow, UInt_t clusterNumber) { return clusterData[padRow][clusterNumber].mEtaSlice; }

  void setClusterEta(UInt_t padRow, UInt_t clusterNumber, Double_t eta) { clusterData[padRow][clusterNumber].mEta = eta; }
  void setClusterEtaSlice(UInt_t padRow, UInt_t clusterNumber, UInt_t etaSlice) { clusterData[padRow][clusterNumber].mEtaSlice = etaSlice; }

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

    clusterData[(UInt_t)padRow].push_back(data);
  }

protected:
private:
  // Information per cluster
  std::vector<std::vector<clusterDataFormat>> clusterData;

  std::unique_ptr<AliHLTTPCSpacePointContainer> spacepoints;
};
}
}
#endif
