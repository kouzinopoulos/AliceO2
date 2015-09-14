/// \file Draw.h
/// \brief Definition of the Draw class
/// \author Charis Kouzinopoulos <charalampos.kouzinopoulos@cern.ch>

#ifndef ALICEO2_HOUGH_DRAW_H_
#define ALICEO2_HOUGH_DRAW_H_

#include "StandardIncludes.h"

#include "TCanvas.h"
#include "TGraph.h"
#include "TGraph2D.h"

#include "ClusterCollection.h"

namespace AliceO2 {
namespace Hough {

class Draw {
public:
  Draw();

  void CartesianClusters1DEtaSlice(ClusterCollection* clusterCollection, int etaSlice);
  void CartesianClusters1DTPCSlice(ClusterCollection* clusterCollection, int TPCSlice);
  void CartesianClusters1D(ClusterCollection* clusterCollection, int TPCSlice, int etaSlice);
  void CartesianClusters2D(ClusterCollection* clusterCollection);

private:
  TCanvas* canvas[4];
  TGraph* graph[3];
  TGraph2D* graph2D[1];
};
}
}
#endif
