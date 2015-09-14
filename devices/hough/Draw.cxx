/// \file Draw.cxx
/// \brief Implementation of the Draw class
/// \author Charis Kouzinopoulos <charalampos.kouzinopoulos@cern.ch>

#include "Draw.h"
#include "TFile.h"
#include "TTree.h"

#include "TFrame.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TApplication.h"

#include "Transform.h"

using namespace std;
using namespace AliceO2::Hough;

Draw::Draw() {}

void Draw::CartesianClusters1DEtaSlice(ClusterCollection* clusterCollection, int etaSlice)
{
  canvas[0] = new TCanvas("CartesianClusters1D", "Cartesian clusters per eta slice");
  graph[0] = new TGraph();

  for (UInt_t padRow = 0; padRow < Transform::GetNRows(); padRow++) {
    for (UInt_t clusterNumber = 0; clusterNumber < clusterCollection->getNumberOfClustersPerPadRow(padRow);
         clusterNumber++) {
      if (clusterCollection->getClusterEtaSlice(padRow, clusterNumber) == etaSlice) {
        graph[0]->SetPoint(clusterNumber, clusterCollection->getClusterX(padRow, clusterNumber),
                           clusterCollection->getClusterY(padRow, clusterNumber));
      }
    }
  }

  graph[0]->SetMarkerStyle(7);
  graph[0]->Draw("AP");
  graph[0]->GetXaxis()->SetRangeUser(80, 240);
  graph[0]->GetXaxis()->SetTitle("Pad Row");

  canvas[0]->Print("cartesianClustersPerEtaSlice.pdf");
}

void Draw::CartesianClusters1DTPCSlice(ClusterCollection* clusterCollection, int TPCSlice)
{
  canvas[1] = new TCanvas("CartesianClusters1D", "Cartesian clusters per TPC slice");
  graph[1] = new TGraph();

  for (UInt_t padRow = 0; padRow < Transform::GetNRows(); padRow++) {
    for (UInt_t clusterNumber = 0; clusterNumber < clusterCollection->getNumberOfClustersPerPadRow(padRow);
         clusterNumber++) {
      if (clusterCollection->getClusterTPCSlice(padRow, clusterNumber) == TPCSlice) {
        graph[1]->SetPoint(clusterNumber, clusterCollection->getClusterX(padRow, clusterNumber),
                           clusterCollection->getClusterY(padRow, clusterNumber));
      }
    }
  }

  graph[1]->SetMarkerStyle(7);
  graph[1]->Draw("AP");
  graph[1]->GetXaxis()->SetRangeUser(80, 240);
  graph[1]->GetXaxis()->SetTitle("Pad Row");

  canvas[1]->Print("cartesianClustersPerTPCSlice.pdf");
}

void Draw::CartesianClusters1D(ClusterCollection* clusterCollection, int TPCSlice, int etaSlice)
{
  canvas[2] = new TCanvas("cartesianClusters1D", "Cartesian clusters per TPC and eta slice");
  graph[2] = new TGraph();

  for (UInt_t padRow = 0; padRow < Transform::GetNRows(); padRow++) {
    for (UInt_t clusterNumber = 0; clusterNumber < clusterCollection->getNumberOfClustersPerPadRow(padRow);
         clusterNumber++) {
      if (clusterCollection->getClusterTPCSlice(padRow, clusterNumber) == TPCSlice &&
          clusterCollection->getClusterEtaSlice(padRow, clusterNumber) == etaSlice) {
        graph[2]->SetPoint(clusterNumber, clusterCollection->getClusterX(padRow, clusterNumber),
                           clusterCollection->getClusterY(padRow, clusterNumber));
      }
    }
  }

  graph[2]->SetMarkerStyle(7);
  graph[2]->Draw("AP");
  graph[2]->GetXaxis()->SetRangeUser(80, 240);
  graph[2]->GetXaxis()->SetTitle("Pad Row");

  canvas[2]->Print("cartesianClustersPerEtaSlicePerTPCSlice.pdf");
}

void Draw::CartesianClusters2D(ClusterCollection* clusterCollection)
{
  canvas[3] = new TCanvas("cartesianClusters2D", "Cartesian clusters", 800, 600);
  graph2D[0] = new TGraph2D();

  for (UInt_t padRow = 0; padRow < Transform::GetNRows(); padRow++) {
    for (UInt_t clusterNumber = 0; clusterNumber < clusterCollection->getNumberOfClustersPerPadRow(padRow);
         clusterNumber++) {
      graph2D[0]->SetPoint(padRow, clusterCollection->getClusterX(padRow, clusterNumber),
                           clusterCollection->getClusterY(padRow, clusterNumber),
                           clusterCollection->getClusterZ(padRow, clusterNumber));
    }
  }

  // Draw with colored dots
  graph2D[0]->SetMarkerStyle(6);
  graph2D[0]->SetTitle("TPC Clusters in x, y, z");
  graph2D[0]->Draw("pcol");

  canvas[3]->Print("cartesianClusters2D.pdf");
}
