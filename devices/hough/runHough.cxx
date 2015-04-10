/// \file runHough.cxx
/// \brief Implementation of a cluster loader
/// \author Charis Kouzinopoulos

/// List of references: [1] Cheshkov, C. "Fast Hough-transform track reconstruction for the ALICE TPC." Nuclear
/// Instruments and Methods in Physics Research Section A: Accelerators, Spectrometers, Detectors and Associated
/// Equipment 566.1 (2006): 35-39.

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

std::unique_ptr<AliHLTTPCSpacePointContainer> spacepoints;
vector<float> clusterCartesianCoordinates;
vector<float> clusterConformalMappingCoordinates;
vector<unsigned int> accumulator;
vector<int> trackCoordinates;

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

// By setting rResolution = 1, r = 18.6 and r = 19.2 for theta = 48 will both cause bin(19,48) to increase. With a
// resoution of rResolution = 10, they will cause bin(186,48) and bin(192,48) to increase respectively.
// Now rResolution = 100 means that there will be 100 bins for r: 0 - 99. Converting r -> rBin and rBin -> r is done by
// methods getRValue and getRBinValue
int rResolution = 100000;

// sin and cos expect values in radians instead of degrees
#define DEG2RAD 0.017453293f

// For clusters, 4 parameters are stored: ID, x, y and z coordinates
int clusterCartesianParameters = 4;

// For clusters after conformal mapping, 5 parameters are stored: ID, a, b, η and η slice
int clusterConformalMappingParameters = 5;

// For tracks found after reconstruction, 5 parameters are stored: η slice, start coordinates (l1,m1) and end
// coordinates (l2,m2)
int trackCoordinateParameters = 5;

// FIXME: convert float into double for bigger precision

void setAccumulatorBin(int etaSlice, int r, int theta)
{
  accumulator[etaSlice * (rResolution * thetaMax) + r * thetaMax + theta]++;
}

int getAccumulatorBin(int etaSlice, int r, int theta)
{
  return accumulator[etaSlice * (rResolution * thetaMax) + r * thetaMax + theta];
}

//Solve getRBinValue to r
float getRValue(int rBin) {
  int rMin = -rMax;
  return (float)((rMax - rMin) * rBin + rMin * rResolution) / rResolution;
}

int getRBinValue(float r) {
  int rMin = -rMax;
  return ((r - rMin) * rResolution)/ (rMax - rMin);
}

void setTracks(int etaSlice, int l1, int m1, int l2, int m2)
{
  trackCoordinates.push_back(etaSlice);
  trackCoordinates.push_back(l1);
  trackCoordinates.push_back(m1);
  trackCoordinates.push_back(l2);
  trackCoordinates.push_back(m2);
}

int getTrackEtaSlice(int trackNumber)
{
  return trackCoordinates[trackNumber * trackCoordinateParameters + 0];
}

int getTrackM1(int trackNumber)
{
  return trackCoordinates[trackNumber * trackCoordinateParameters + 1];
}

int getTrackL1(int trackNumber)
{
  return trackCoordinates[trackNumber * trackCoordinateParameters + 2];
}

int getTrackM2(int trackNumber)
{
  return trackCoordinates[trackNumber * trackCoordinateParameters + 3];
}

int getTrackL2(int trackNumber)
{
  return trackCoordinates[trackNumber * trackCoordinateParameters + 4];
}

int getNumberOfTracks()
{
  return trackCoordinates.size() / trackCoordinateParameters;
}

float getClusterID(int clusterNumber)
{
  return clusterCartesianCoordinates[clusterNumber * clusterCartesianParameters + 0];
}

float getClusterCartesianX(int clusterNumber)
{
  return clusterCartesianCoordinates[clusterNumber * clusterCartesianParameters + 1];
}

float getClusterCartesianY(int clusterNumber)
{
  return clusterCartesianCoordinates[clusterNumber * clusterCartesianParameters + 2];
}

float getClusterCartesianZ(int clusterNumber)
{
  return clusterCartesianCoordinates[clusterNumber * clusterCartesianParameters + 3];
}

void setClusterCartesianParameters(AliHLTUInt32_t clusterID, float x, float y, float z)
{
  clusterCartesianCoordinates.push_back((float)clusterID);
  clusterCartesianCoordinates.push_back(x);
  clusterCartesianCoordinates.push_back(y);
  clusterCartesianCoordinates.push_back(z);
}

float getClusterConformalMappingAlpha(int clusterNumber)
{
  return clusterConformalMappingCoordinates[clusterNumber * clusterConformalMappingParameters + 1];
}

float getClusterConformalMappingBeta(int clusterNumber)
{
  return clusterConformalMappingCoordinates[clusterNumber * clusterConformalMappingParameters + 2];
}

float getClusterConformalMappingEta(int clusterNumber)
{
  return clusterConformalMappingCoordinates[clusterNumber * clusterConformalMappingParameters + 3];
}

int getClusterConformalMappingEtaSlice(int clusterNumber)
{
  return clusterConformalMappingCoordinates[clusterNumber * clusterConformalMappingParameters + 4];
}

int setClusterConformalMappingAlpha(int clusterNumber, float alpha)
{
  clusterConformalMappingCoordinates[clusterNumber * clusterConformalMappingParameters + 1] = alpha;
}

int setClusterConformalMappingBeta(int clusterNumber, float beta)
{
  clusterConformalMappingCoordinates[clusterNumber * clusterConformalMappingParameters + 2] = beta;
}

int setClusterConformalMappingEta(int clusterNumber, float eta)
{
  clusterConformalMappingCoordinates[clusterNumber * clusterConformalMappingParameters + 3] = eta;
}

int setClusterConformalMappingEtaSlice(int clusterNumber, float etaSlice)
{
  clusterConformalMappingCoordinates[clusterNumber * clusterConformalMappingParameters + 4] = etaSlice;
}

/// Calculate an approximate value for η. See [1]:p8 for more information. Values below taken from
/// AliHLTConfMapPoint.cxx
Double_t calculatePseudoRapidity(int clusterNumber)
{

  Double_t radial = sqrt(getClusterCartesianX(clusterNumber) * getClusterCartesianX(clusterNumber) +
                         getClusterCartesianY(clusterNumber) * getClusterCartesianY(clusterNumber) +
                         getClusterCartesianZ(clusterNumber) * getClusterCartesianZ(clusterNumber));
  Double_t eta =
    0.5 * log((radial + getClusterCartesianZ(clusterNumber)) / (radial - getClusterCartesianZ(clusterNumber)));

  return eta;
}

void drawTracks(int totalNumberOfClusters)
{
  TCanvas* trackCanvas = new TCanvas("trackCanvas", "Reconstructed tracks", 0, 0, 800, 600);
  TGraph* trackGraph[totalNumberOfClusters];

  cout << "Number of tracks: " << getNumberOfTracks() << endl;

  //DEBUG
  if (trackCoordinates.size() > 0) {
    cout << getTrackM1(0) << "," << getTrackL1(0) << " - " << getTrackM2(0) << "," << getTrackL2(0) << endl;
  }

  for (int i = 0; i < getNumberOfTracks(); i++) {

    //DEBUG
    if (getTrackEtaSlice(i) != 79)
      continue;

    trackGraph[i] = new TGraph();

    trackGraph[i]->SetPoint(i * 2 + 0, getTrackM1(i), getTrackL1(i));
    trackGraph[i]->SetPoint(i * 2 + 1, getTrackM2(i), getTrackL2(i));

    if (i == 0) {
      trackGraph[i]->Draw("AC");
    } else {
      trackGraph[i]->Draw("CP");
    }
  }
  trackCanvas->Print("tracks.pdf");
}

void drawAccumulatorCurves(int totalNumberOfClusters)
{
  // FIXME: The method calculates again r. Is it possible to obtain this information directly from the accumulator?
  TCanvas* c4 = new TCanvas("c4", "Accumulator curves", 0, 0, 800, 600);
  TGraph* g[totalNumberOfClusters];

  for (Int_t i = 0; i < totalNumberOfClusters; i++) {
    float x = getClusterCartesianX(i);
    float y = getClusterCartesianY(i);

    g[i] = new TGraph();

    for (Int_t theta = 0; theta < thetaMax; theta++) {
      double r = (x * cos(theta * DEG2RAD)) + (y * sin(theta * DEG2RAD));
      g[i]->SetPoint(i * thetaMax + theta, theta, r);
    }

    g[i]->SetMarkerStyle(1);

    if (i == 0) {
      g[i]->Draw("AC");
    } else {
      g[i]->Draw("CP");
    }
  }
  c4->Print("accumulatorCurves.pdf");
}

void drawAccumulatorHistogram()
{
  TCanvas* c5 = new TCanvas("c5", "Accumulator histogram", 0, 0, 800, 600);
  TH2F* h =
    new TH2F("h", "Accumulator histogram", rResolution, 0, rResolution, thetaMax, 0, thetaMax);

  h->SetFillColor(46);

  for (Int_t r = 0; r < rResolution; r++) {
    for (Int_t theta = 0; theta < thetaMax; theta++) {
      if (getAccumulatorBin(15, r, theta) > 0) {
        h->SetBinContent(r, theta, getAccumulatorBin(15, r, theta));
      }
    }
  }

  h->Draw("LEGO");
  c5->Print("accumulatorHistogram.pdf");
}

void drawCartesianClusters1D(int etaSlice, int totalNumberOfClusters)
{
  TCanvas* cartesianClustersCanvas1D = new TCanvas("cartesianClustersCanvas1D", "Cartesian clusters");
  TGraph* cartesianClustersGraph1D = new TGraph();

  for (Int_t i = 0; i < totalNumberOfClusters; i++) {
    if (getClusterConformalMappingEtaSlice(i) == etaSlice) {
      cartesianClustersGraph1D->SetPoint(i, getClusterCartesianX(i), getClusterCartesianY(i));
    }
  }

  cartesianClustersGraph1D->SetMarkerStyle(7);
  cartesianClustersGraph1D->Draw("AP");
  //cartesianClustersGraph1D->GetXaxis()->SetRangeUser(70,250);

  cartesianClustersCanvas1D->Print("cartesianClusters1D.pdf");
}

void drawConformalMappingClusters1D(int etaSlice, int totalNumberOfClusters)
{
  TCanvas* c7 = new TCanvas("c1", "Conformal mapping clusters", 0, 0, 800, 600);
  TGraph* gr7 = new TGraph();

  for (Int_t i = 0; i < totalNumberOfClusters; i++) {
    if (getClusterConformalMappingEtaSlice(i) == etaSlice) {
      gr7->SetPoint(i, getClusterConformalMappingAlpha(i), getClusterConformalMappingBeta(i));
    }
  }

  gr7->SetMarkerStyle(7);
  gr7->Draw("AP");

  c7->Print("conformalMappingClusters1D.pdf");
}

void drawCartesianClusters(int totalNumberOfClusters)
{
  TCanvas* cartesianClustersCanvas = new TCanvas("cartesianClustersCanvas", "Cartesian clusters", 800, 600);
  TGraph2D* cartesianClustersGraph2D = new TGraph2D();

  for (Int_t i = 0; i < totalNumberOfClusters; i++) {
    // cout << "graph: " << i << " " << getClusterCartesianX(i) << " " << getClusterCartesianY(i) << " " <<
    // getClusterCartesianZ(i) << " " << endl;
    cartesianClustersGraph2D->SetPoint(i, getClusterCartesianX(i), getClusterCartesianY(i), getClusterCartesianZ(i));
  }

  // Draw with colored dots
  cartesianClustersGraph2D->SetMarkerStyle(6);
  cartesianClustersGraph2D->SetTitle("TPC Clusters in x, y, z");
  cartesianClustersGraph2D->Draw("pcol");

  cartesianClustersCanvas->Print("cartesianClusters2D.pdf");
}

void drawConformalMappingClusters(int totalNumberOfClusters)
{
  TCanvas* conformalMappingClustersCanvas =
    new TCanvas("conformalMappingClustersCanvas", "Conformal Mapping clusters", 800, 600);
  TGraph2D* conformalMappingClustersGraph2D = new TGraph2D(totalNumberOfClusters);

  for (Int_t i = 0; i < totalNumberOfClusters; i++) {
    conformalMappingClustersGraph2D->SetPoint(i, getClusterConformalMappingAlpha(i), getClusterConformalMappingBeta(i),
                                              getClusterConformalMappingEta(i));
  }

  // Draw with colored dots
  conformalMappingClustersGraph2D->SetMarkerStyle(6);
  conformalMappingClustersGraph2D->SetTitle("TPC Clusters in a, b, eta");
  conformalMappingClustersGraph2D->Draw("pcol");

  conformalMappingClustersCanvas->Print("conformalMappingClusters.pdf");
}

// Determine if the current value is the maximum in an area of 9x9 points
int localAccumulatorMaxima(int etaSlice, int r, int theta)
{
  int max = getAccumulatorBin(etaSlice, r, theta);

  for (int deltaR = -4; deltaR < 5; deltaR++) {
    for (int deltaTheta = -4; deltaTheta < 5; deltaTheta++) {
      if ((deltaR + r >= 0 && deltaR + r < rMax) && (deltaTheta + theta >= 0 && deltaTheta + theta < thetaMax)) {
        if (getAccumulatorBin(etaSlice, r + deltaR, theta + deltaTheta) > max) {
          return getAccumulatorBin(etaSlice, r + deltaR, theta + deltaTheta);
        }
      }
    }
  }
  return max;
}

/// Determine the start and end coordinates of all the tracks, calculated from the local maxima of the accumulator
/// values that exceed a given threshold
void trackFinding(int etaSlice)
{
  // The track will have coordinates (l1,m1), (l2,m2)
  float l1, m1, l2, m2;

  int trackCount = 0;

  for (Int_t rBin = 0; rBin < rResolution; rBin++) {
    for (Int_t theta = 0; theta < thetaMax; theta++) {
      if (getAccumulatorBin(etaSlice, rBin, theta) >= houghThreshold) {
        if (localAccumulatorMaxima(etaSlice, rBin, theta) > getAccumulatorBin(etaSlice, rBin, theta)) {
          continue;
        }

        if (theta >= 45 && theta <= 135) {
          //y = (r - x cos(t)) / sin(t)
          l1 = 0.0;
          m1 = (getRValue(rBin) - l1 * cos(theta * DEG2RAD)) / sin(theta * DEG2RAD);
          l2 = aMax;
          m2 = (getRValue(rBin) - l2 * cos(theta * DEG2RAD)) / sin(theta * DEG2RAD);
        } else {
          //x = (r - y sin(t)) / cos(t);
          m1 = 0.0;
          l1 = (getRValue(rBin) - m1 * sin(theta * DEG2RAD)) / cos(theta * DEG2RAD);
          m2 = bMax;
          l2 = (getRValue(rBin) - m2 * sin(theta * DEG2RAD)) / cos(theta * DEG2RAD);
        }
        //cout << "Found a track at [" << etaSlice << "][" << getRValue(rBin) << "][" << theta << "] with a peak of " << getAccumulatorBin(etaSlice, rBin, theta) << " points"
        //     << endl << "Track coordinates: " << l1 << "," << m1 << " - " << l2 << "," << m2 << endl;
        //cout << "xMax: " << xMax << " yMax: " << yMax << endl;
        trackCount++;

        setTracks(etaSlice, l1, m1, l2, m2);
      }
    }
  }
  if (trackCount > 0) {
    cout << "rMax " << rMax << " rResolution " << rResolution << endl;
    cout << "Found " << trackCount << " tracks for η slice " << etaSlice << endl;
  }
}

void conformalMapping(int totalNumberOfClusters)
{
  for (int i = 0; i < totalNumberOfClusters; i++) {
    float x = getClusterCartesianX(i);
    float y = getClusterCartesianY(i);

    // Equation (2) from paper [1]
    float alpha = x / (x * x + y * y);
    float beta = y / (x * x + y * y);
    float eta = getClusterConformalMappingEta(i);
    int etaSlice;

    if (etaMax - etaMin != 0) {
      etaSlice = (etaResolution * (eta - etaMin)) / (etaMax - etaMin);
    } else {
      etaSlice = 0;
    }

    //cout << "i: " << i << " X: " << x << " Y: " << y << " A: " << alpha << " B: " << beta << " η: " << eta << " η slice: " << etaSlice
    //     << endl;

    setClusterConformalMappingAlpha(i, alpha);
    setClusterConformalMappingBeta(i, beta);
    setClusterConformalMappingEtaSlice(i, etaSlice);
  }
}

void transformCartesian(int totalNumberOfClusters)
{
  thetaMax = 180;
  // Trigonometrically, the maximum distance is designated by the square root of the summation of the squares of the x
  // and y dimensions
  rMax = ceil(sqrt(xMax * xMax + yMax * yMax));

  // The lines will have -rMax <= r <= rMax and 0 <= theta <= thetaMax. The total space needed is thus 2 * rMax *
  // thetaMax
  accumulator.resize(2 * rMax * thetaMax * rResolution, 0);

  for (Int_t i = 0; i < totalNumberOfClusters; i++) {
    float x = getClusterCartesianX(i);
    float y = getClusterCartesianY(i);

    for (Int_t theta = 0; theta < thetaMax; theta++) {
      double r = (x * cos(theta * DEG2RAD)) + (y * sin(theta * DEG2RAD));
      //      cout << "x: " << x << " y: " << y << " theta: " << theta << " r: " << r
      //           << " round r: " << round((r + rMax) * rResolution) << endl;

      // Use a discreet value for r by rounding r. Then, rMax is added to r to change its range from -rMax <= r <= rMax
      // to 0 <= r <= 2 * rMax
      setAccumulatorBin(0, round((r + rMax) * rResolution), theta);
    }
  }
}

void transformConformalMapping(int totalNumberOfClusters)
{
  thetaMax = 180;
  // Trigonometrically, the maximum distance is designated by the square root of the summation of the squares of the x
  // and y dimensions
  rMax = ceil(sqrt(aMax * aMax + bMax * bMax));

  //cout << "rMin: " << -rMax << " rMax: " << rMax << endl;

  // The lines will have -rMax <= r <= rMax and 0 <= theta <= thetaMax. The total space needed is thus 2 * rMax *
  // thetaMax
  accumulator.resize(etaResolution * rResolution * thetaMax + rResolution * thetaMax + thetaMax, 0);

  cout << "Accumulator size: " << 2 * rMax * thetaMax * rResolution * etaResolution << endl;
  cout << "Alternative accumulator size: " << etaResolution * rResolution * thetaMax + rResolution * thetaMax + thetaMax << endl;

  for (Int_t i = 0; i < totalNumberOfClusters; i++) {
    float a = getClusterConformalMappingAlpha(i);
    float b = getClusterConformalMappingBeta(i);
    float eta = getClusterConformalMappingEta(i);
    int etaSlice = getClusterConformalMappingEtaSlice(i);

    for (Int_t theta = 0; theta < thetaMax; theta++) {
      double r = (a * cos(theta * DEG2RAD)) + (b * sin(theta * DEG2RAD));
      int rBin = getRBinValue(r);

      if (etaSlice == 20 ) {
        cout << "Cluster: " << i << " a: " << a << " b: " << b << " theta: " << theta << " rMin: " << -rMax << " r: " << r
             << " rMax: " << rMax << " rBin: " << rBin << " eta: " << eta << " etaSlice: " << etaSlice << endl;
      }

      // Use a discreet value for r by rounding r. Then, rMax is added to r to change its range from -rMax <= r <= rMax
      // to 0 <= r <= 2 * rMax
      setAccumulatorBin(etaSlice, rBin, theta);
    }
  }
}

/// Determine the minimum and maximum values of the pseudorapidity (eta). That way, the TPC digits can be transformed in
/// two dimensions instead of three in slices of similar pseudorapidity
void determineMinMaxEta(int totalNumberOfClusters)
{
  for (Int_t i = 0; i < totalNumberOfClusters; i++) {

    float eta = calculatePseudoRapidity(i);

    // Store eta to be later used by the conformalMapping method
    setClusterConformalMappingEta(i, eta);

    // Set initial values to etaMin and etaMax
    if (i == 0) {
      etaMin = eta;
      etaMax = eta;
      continue;
    }

    //cout << "i: " << i << " eta: " << eta << " " << etaMin << " " << etaMax << endl;

    if (eta < etaMin) {
      etaMin = eta;
    } else if (eta > etaMax) {
      etaMax = eta;
    }
}
  cout << "Minimum eta: " << etaMin << " Maximum eta: " << etaMax << endl;
}

/// Determine the maximum ceiling to x,y,z coordinates from the clusterCartesianCoordinates vector to later allocate the
/// hough transform accumulator
void determineMaxCartesianDimensions(int totalNumberOfClusters)
{
  for (Int_t i = 0; i < totalNumberOfClusters; i++) {
    if (abs(getClusterCartesianX(i)) > xMax) {
      xMax = ceil(abs(getClusterCartesianX(i)));
    }
    if (abs(getClusterCartesianY(i)) > yMax) {
      yMax = ceil(abs(getClusterCartesianY(i)));
    }
    if (abs(getClusterCartesianZ(i)) > zMax) {
      zMax = ceil(abs(getClusterCartesianZ(i)));
    }
  }
  cout << "xMax: " << xMax << " yMax: " << yMax << " zMax: " << zMax << endl;
}

/// Determine the maximum ceiling to a and b conformal mapping coordinates from the clusterCartesianCoordinates vector
/// to later allocate the hough transform accumulator
void determineMinMaxAB(int totalNumberOfClusters)
{
  for (Int_t i = 0; i < totalNumberOfClusters; i++) {
    if (i == 0) {
      aMin = floor(getClusterConformalMappingAlpha(i));
      aMax = ceil(abs(getClusterConformalMappingAlpha(i)));
      bMin = floor(getClusterConformalMappingBeta(i));
      bMax = ceil(abs(getClusterConformalMappingBeta(i)));
      continue;
    }
    if ( getClusterConformalMappingAlpha(i) < aMin ) {
      aMin = floor(getClusterConformalMappingAlpha(i));
    }
    else if (abs(getClusterConformalMappingAlpha(i)) > aMax) {
      aMax = ceil(abs(getClusterConformalMappingAlpha(i)));
    }

    if ( getClusterConformalMappingBeta(i) < bMin ) {
      bMin = floor(getClusterConformalMappingBeta(i));
    } else if (abs(getClusterConformalMappingBeta(i)) > bMax) {
      bMax = ceil(abs(getClusterConformalMappingBeta(i)));
    }
  }
  cout << "aMax: " << aMax << " bMax: " << bMax << endl;
}

void printData(int totalNumberOfClusters)
{
  cout << "Cluster ID" << setw(13) << "X coordinate" << setw(13) << "Y coordinate" << setw(13) << "Z coordinate"
       << endl;

  for (int i = 0; i < totalNumberOfClusters; i++) {
    cout << (AliHLTUInt32_t)getClusterID(i) << setw(13) << getClusterCartesianX(i) << setw(13)
         << getClusterCartesianY(i) << setw(13) << getClusterCartesianZ(i) << endl;
  }
}

int processData(std::string dataPath, std::string dataType, std::string dataOrigin)
{
  // Open data file for reading
  std::ifstream inputData(dataPath.c_str(), std::ifstream::binary);
  if (!inputData) {
    std::cerr << "Error, cluster data file " << dataPath << " could not be accessed" << endl;
    std::exit(1);
  }

  // Get length of file
  inputData.seekg(0, inputData.end);
  int dataLength = inputData.tellg();
  inputData.seekg(0, inputData.beg);

  // Allocate memory and read file to memory
  char* inputBuffer = new char[dataLength];
  inputData.read(inputBuffer, dataLength);
  inputData.close();

  // Retrieve the TPC slice and partition from the filename
  std::string currentSliceString(dataPath, dataPath.length() - 6, 2);
  std::string currentPartitionString(dataPath, dataPath.length() - 2, 2);

  AliHLTUInt8_t currentSlice = std::stoul(currentSliceString, nullptr, 16);
  AliHLTUInt8_t currentPartition = std::stoul(currentPartitionString, nullptr, 16);

  // Initialize a cluster point collection
  spacepoints = std::unique_ptr<AliHLTTPCSpacePointContainer>(new AliHLTTPCSpacePointContainer);
  if (!spacepoints.get()) {
    std::cerr << "Error, could not create a space point collection" << endl;
    std::exit(1);
  }

  // Create an AliHLTComponentBlockData object, fill it with default values and then set its pointer to the data buffer
  AliHLTComponentBlockData bd;
  AliHLTComponent::FillBlockData(bd);
  bd.fPtr = inputBuffer;
  bd.fSize = dataLength;
  // bd.fDataType=kAliHLTVoidDataType;
  AliHLTComponent::SetDataType(bd.fDataType, dataType.c_str(), dataOrigin.c_str());
  bd.fSpecification = kAliHLTVoidDataSpec;

  // Set slice and partition
  AliHLTTPCDefinitions::EncodeDataSpecification(currentSlice, currentSlice, currentPartition, currentPartition);

  // Add the AliHLTComponentBlockData object to AliHLTTPCSpacePointContainer
  int numberOfClusters = spacepoints->AddInputBlock(&bd);

  // cout << *spacepoints << endl;

  // Retrieve the cluster information from AliHLTTPCSpacePointContainer
  std::vector<AliHLTUInt32_t> clusterIDs;
  spacepoints->GetClusterIDs(clusterIDs);

  // Append the cluster IDs and their X, Y and Z coordinates to the clusterCartesianCoordinates vector
  for (vector<AliHLTUInt32_t>::const_iterator element = clusterIDs.begin(); element != clusterIDs.end(); element++) {
    AliHLTUInt32_t clusterID = *element;

    setClusterCartesianParameters(clusterID, spacepoints->GetX(clusterID), spacepoints->GetY(clusterID),
                                  spacepoints->GetZ(clusterID));
  }

  // De-allocate memory space
  if (inputBuffer) {
    delete[] inputBuffer;
  }
  inputBuffer = NULL;

  return numberOfClusters;
}

int main(int argc, char** argv)
{
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <event number>" << endl;
    std::exit(1);
  }

  // Create data path
  std::string dataFilename = "emulated-tpc-clusters/event";
  dataFilename += argv[1];

  boost::filesystem::path dataPath(dataFilename);
  boost::filesystem::directory_iterator endIterator;

  typedef std::multimap<std::time_t, boost::filesystem::path> result_set_t;
  result_set_t result_set;

  std::string dataType = "CLUSTERS", dataOrigin = "TPC ";

  int totalNumberOfClusters = 0, totalNumberOfDataFiles = 0;

  // Traverse the filesystem and execute processData for each cluster file found
    if (boost::filesystem::exists(dataPath) && boost::filesystem::is_directory(dataPath)) {
      for (boost::filesystem::directory_iterator directoryIterator(dataPath); directoryIterator != endIterator;
           ++directoryIterator) {
        if (boost::filesystem::is_regular_file(directoryIterator->status())) {
          totalNumberOfClusters += processData(directoryIterator->path().string(), dataType, dataOrigin);
          totalNumberOfDataFiles++;
        }
      }
    } else {
      std::cerr << "Path " << dataPath.string() << "/ could not be found or does not contain any valid data files"
                << endl;
      exit(1);
    }

    cout << "Added " << totalNumberOfClusters << " clusters from " << totalNumberOfDataFiles << " data files" << endl;

  // DEBUG
/*  totalNumberOfClusters = 100;
  for (int kk = 0; kk < 6; kk++) {
    setClusterCartesianParameters(kk, (float)9.0 + kk, (float)17.0, 1);
  }

  for (int kk = 6; kk < 11; kk++) {
    setClusterCartesianParameters(kk, 9.0 + kk - 5, 17.0 + kk - 5, 1);
  }

  for (int kk = 11; kk < 150; kk++) {
    setClusterCartesianParameters(kk, 10, 19.0 + kk - 11, 1);
  }
*/
  totalNumberOfClusters = 20000;
  // printData(totalNumberOfClusters);

  // Allocate space for the conformal mapping parameter vector
  clusterConformalMappingCoordinates.resize(totalNumberOfClusters * clusterConformalMappingParameters);

  // Determine the minimum and maximum values of eta. That way the TPC digits can be grouped into pseudorapidity bins
  determineMinMaxEta(totalNumberOfClusters);

  // Convert the TPC cluster coordinates from the cartesian to the conformal mapping system
  conformalMapping(totalNumberOfClusters);

  // drawConformalMappingClusters1D(totalNumberOfClusters);

  // Determine the maximum dimensions of the clusters for the accumulator
  determineMaxCartesianDimensions(totalNumberOfClusters);

  // Perform the hough transform on the TPC clusters for the cartesian system
  // transformCartesian(totalNumberOfClusters);

  // Determine the maximum dimensions of the clusters for the accumulator
  determineMinMaxAB(totalNumberOfClusters);

  // Perform the hough transform on the TPC clusters for the conformal mapping system
  transformConformalMapping(totalNumberOfClusters);

  // Locate tracks for each η slice
/*  for (int etaSlice = 0; etaSlice <= etaResolution; etaSlice++) {
    trackFinding(etaSlice);
  }*/
  trackFinding(15);

  //drawAccumulatorHistogram();

  drawCartesianClusters(totalNumberOfClusters);
  drawCartesianClusters1D(15, totalNumberOfClusters);
  drawConformalMappingClusters(totalNumberOfClusters);
  drawConformalMappingClusters1D(15, totalNumberOfClusters);
  drawTracks(totalNumberOfClusters);

  int i;
  int array[etaResolution + 1];
  for ( i = 0; i <= etaResolution; i++) {
    array[i] = 0;
  }


  for ( i = 0; i < totalNumberOfClusters; i++) {
    array[getClusterConformalMappingEtaSlice(i)]++;
  }

  for ( i = 0; i <= etaResolution; i++) {
    cout << i << " " << array[i] << endl;
  }

/*
  for (int theta = 0; theta < thetaMax; theta++) {
    for (int rBin = 0; rBin < rResolution; rBin++) {
      if ( getAccumulatorBin(15, rBin, theta ) > 0 ) {
        cout << theta << " " << rBin << " " << getAccumulatorBin(15, rBin, theta ) << endl;
      }
    }
  }
*/

  return 0;

  // Clear the track coordinates so the vector can be used for the conformal mapping tracker
  trackCoordinates.clear();

  // Draw the accumulator data as curves
  drawAccumulatorCurves(totalNumberOfClusters);

  return 0;
}
