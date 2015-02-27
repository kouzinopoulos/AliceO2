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

int xMax = 0.0;
int yMax = 0.0;
int zMax = 0.0;

int rMax;
int rResolution;
int thetaMax;

int houghThreshold = 6;

// sin and cos expect values in radians instead of degrees
#define DEG2RAD 0.017453293f

// We store 4 parameters, ID, X, Y, Z
int clusterParameters = 4;

// FIXME: convert float into double for bigger precision

void setAccumulatorBin(int r, int theta)
{
  accumulator[r * thetaMax + theta]++;
}

int getAccumulatorBin(int r, int theta)
{
  return accumulator[r * thetaMax + theta];
}

void setTracks(int startX, int startY, int endX, int endY)
{
  trackCoordinates.push_back(startX);
  trackCoordinates.push_back(startY);
  trackCoordinates.push_back(endX);
  trackCoordinates.push_back(endY);
}

int getTrackStartX(int trackNumber)
{
  return trackCoordinates[trackNumber * 4 + 0];
}

int getTrackStartY(int trackNumber)
{
  return trackCoordinates[trackNumber * 4 + 1];
}

int getTrackEndX(int trackNumber)
{
  return trackCoordinates[trackNumber * 4 + 2];
}

int getTrackEndY(int trackNumber)
{
  return trackCoordinates[trackNumber * 4 + 3];
}

int getNumberOfTracks()
{
  return trackCoordinates.size() / 4;
}

float getClusterID(int clusterNumber)
{
  return clusterCartesianCoordinates[clusterNumber * clusterParameters + 0];
}

float getClusterCartesianX(int clusterNumber)
{
  return clusterCartesianCoordinates[clusterNumber * clusterParameters + 1];
}

float getClusterCartesianY(int clusterNumber)
{
  return clusterCartesianCoordinates[clusterNumber * clusterParameters + 2];
}

float getClusterCartesianZ(int clusterNumber)
{
  return clusterCartesianCoordinates[clusterNumber * clusterParameters + 3];
}

float getClusterConformalMappingX(int clusterNumber)
{
  return clusterConformalMappingCoordinates[clusterNumber * clusterParameters + 1];
}

float getClusterConformalMappingY(int clusterNumber)
{
  return clusterConformalMappingCoordinates[clusterNumber * clusterParameters + 2];
}

float getClusterConformalMappingZ(int clusterNumber)
{
  return clusterConformalMappingCoordinates[clusterNumber * clusterParameters + 3];
}

void drawTracks()
{
  TCanvas* c2 = new TCanvas("c2", "Reconstructed tracks", 0, 0, 800, 600);
  TGraph* gr[getNumberOfTracks()];

  cout << "Number of tracks: " << getNumberOfTracks() << endl;

  for (int i = 0; i < getNumberOfTracks(); i++) {

    gr[i] = new TGraph();

    gr[i]->SetPoint(i * 2 + 0, getTrackStartX(i), getTrackStartY(i));
    gr[i]->SetPoint(i * 2 + 1, getTrackEndX(i), getTrackEndY(i));

    if (i == 0) {
      gr[i]->Draw("AC");
    } else {
      gr[i]->Draw("CP");
    }

    gr[i]->GetXaxis()->SetLimits(130, 161);
    gr[i]->SetMinimum(0);
    gr[i]->SetMaximum(30);
  }
  c2->Print("tracks.pdf");
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

void drawAccumulatorHistogram(int totalNumberOfClusters)
{
  TCanvas* c5 = new TCanvas("c5", "Accumulator histogram", 0, 0, 800, 600);
  TH2F* h =
    new TH2F("h", "Accumulator histogram", 2 * rMax * rResolution, 0, 2 * rMax * rResolution, thetaMax, 0, thetaMax);

  h->SetFillColor(46);

  for (Int_t r = 0; r < 2 * rMax * rResolution; r++) {
    for (Int_t theta = 0; theta < thetaMax; theta++) {
      h->SetBinContent(r, theta, getAccumulatorBin(r, theta));
    }
  }

  h->Draw("LEGO");
  c5->Print("accumulatorHistogram.pdf");
}

void drawCartesianClusters1D(int totalNumberOfClusters)
{
  TCanvas* c3 = new TCanvas("c1", "Cartesian clusters", 0, 0, 800, 600);
  TGraph* gr3 = new TGraph();

  for (Int_t i = 0; i < totalNumberOfClusters; i++) {
    gr3->SetPoint(i, getClusterCartesianX(i), getClusterCartesianY(i));
  }

  // Draw with colored dots
  gr3->SetMarkerStyle(1);
  gr3->Draw("A*");

  c3->Print("cartesian_clusters.pdf");
}

void drawConformalMappingClusters1D(int totalNumberOfClusters)
{
  TCanvas* c7 = new TCanvas("c1", "Conformal mapping clusters", 0, 0, 800, 600);
  TGraph* gr7 = new TGraph();

  for (Int_t i = 0; i < totalNumberOfClusters; i++) {
    gr7->SetPoint(i, getClusterConformalMappingX(i), getClusterConformalMappingY(i));
  }

  // Draw with colored dots
  gr7->SetMarkerStyle(1);
  gr7->Draw("A*");

  c7->Print("conformal_mapping_clusters.pdf");
}

void drawClusters(int totalNumberOfClusters)
{
  TCanvas* c1 = new TCanvas("c1", "Cartesian clusters", 0, 0, 800, 600);
  TGraph2D* dt = new TGraph2D(10000);

  for (Int_t i = 0; i < totalNumberOfClusters; i++) {
    dt->SetPoint(i, getClusterCartesianX(i), getClusterCartesianY(i), getClusterCartesianZ(i));
  }

  // Draw with colored dots
  dt->SetMarkerStyle(1);
  dt->Draw("pcol");

  c1->Print("clusters.pdf");
}

// Determine if the current value is the maximum in an area of 9x9 points
int localAccumulatorMaxima(int r, int theta)
{
  int max = getAccumulatorBin(r, theta);

  for (int deltaR = -4; deltaR < 5; deltaR++) {
    for (int deltaTheta = -4; deltaTheta < 5; deltaTheta++) {
      if ((deltaR + r >= 0 && deltaR + r < rMax) && (deltaTheta + theta >= 0 && deltaTheta + theta < thetaMax)) {
        if (getAccumulatorBin(r + deltaR, theta + deltaTheta) > max) {
          return getAccumulatorBin(r + deltaR, theta + deltaTheta);
        }
      }
    }
  }
  return max;
}

void calculateTracks(int totalNumberOfClusters)
{
  /*  // DEBUG
    for (Int_t r = 0; r < 2 * rMax * rResolution; r++) {
      cout << r << ": " << setw(3);

      for (Int_t theta = 0; theta < thetaMax; theta++) {
        cout << theta << ": " << getAccumulatorBin(r, theta) << " ";
      }
      cout << endl;
    }
  */
  int startX, startY, endX, endY;

  for (Int_t r = 0; r < 2 * rMax * rResolution; r++) {
    for (Int_t theta = 0; theta < thetaMax; theta++) {
      if (getAccumulatorBin(r, theta) >= houghThreshold) {
        if (localAccumulatorMaxima(r, theta) > getAccumulatorBin(r, theta)) {
          continue;
        }

        if (theta >= 45 && theta <= 135) {
          startX = 0;
          startY = (((double)r / rResolution) - rMax - startX * cos(theta * DEG2RAD)) / sin(theta * DEG2RAD);
          endX = xMax;
          endY = (((double)r / rResolution) - rMax - endX * cos(theta * DEG2RAD)) / sin(theta * DEG2RAD);
        } else {
          startY = 0;
          startX = (((double)r / rResolution) - rMax - startY * sin(theta * DEG2RAD)) / cos(theta * DEG2RAD);
          endY = yMax;
          endX = (((double)r / rResolution) - rMax - endY * sin(theta * DEG2RAD)) / cos(theta * DEG2RAD);
        }
        cout << "Found a track at [" << ((double)r / rResolution) - rMax << "][" << theta << "]" << endl
             << "Track coordinates: " << startX << "," << startY << " - " << endX << "," << endY << endl;

        setTracks(startX, startY, endX, endY);
      }
    }
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

    // cout << "A: " << alpha << " B: " << beta << endl;

    clusterConformalMappingCoordinates.push_back((float)getClusterID(i));
    clusterConformalMappingCoordinates.push_back(alpha);
    clusterConformalMappingCoordinates.push_back(beta);
    clusterConformalMappingCoordinates.push_back(0.0);
  }
}

void transform(int totalNumberOfClusters)
{
  thetaMax = 180;
  // Trigonometrically, the maximum distance is designated by the square root of the summation of the squares of the x
  // and y dimensions
  rMax = ceil(sqrt(xMax * xMax + yMax * yMax));

  // By setting rResolution = 1, r = 18.6 and r = 19.2 for theta = 48 will both cause bin(19,48) to increase. With a
  // resoution of rResolution = 10, they will cause bin(186,48) and bin(192,48) to increase respectively
  rResolution = 100;

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
      setAccumulatorBin(round((r + rMax) * rResolution), theta);
    }
  }
}

/// Determine the maximum ceiling to X,Y,Z coordinates from the clusterCartesianCoordinates vector to later allocate the
/// hough
/// transform accumulator
void determineMaxDimensions(int totalNumberOfClusters)
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

void printData(int totalNumberOfClusters)
{
  cout << "Cluster ID" << setw(13) << "X coordinate" << setw(13) << "Y coordinate" << setw(13) << "Z coordinate"
       << endl;

  for (int i = 0; i < totalNumberOfClusters; i++) {
    cout << (AliHLTUInt32_t)getClusterID(i) << setw(13) << getClusterCartesianX(i) << setw(13)
         << getClusterCartesianY(i) << setw(13) << getClusterCartesianZ(i) << endl;
  }
}

void addDataToCoordinatesVector(AliHLTUInt32_t clusterID, float XCoordinate, float YCoordinate, float ZCoordinate)
{
  clusterCartesianCoordinates.push_back((float)clusterID);
  clusterCartesianCoordinates.push_back(XCoordinate);
  clusterCartesianCoordinates.push_back(YCoordinate);
  clusterCartesianCoordinates.push_back(ZCoordinate);
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

    addDataToCoordinatesVector(clusterID, spacepoints->GetX(clusterID), spacepoints->GetY(clusterID),
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

  /*  // DEBUG
    totalNumberOfClusters = 15;
    for (int kk = 0; kk < 6; kk++) {

      clusterCartesianCoordinates[kk * 4] = kk;
      clusterCartesianCoordinates[kk * 4 + 1] = 9.0 + kk;
      clusterCartesianCoordinates[kk * 4 + 2] = 17.0;
      clusterCartesianCoordinates[kk * 4 + 3] = 0;
    }

    for (int kk = 6; kk < 11; kk++) {

      clusterCartesianCoordinates[kk * 4] = kk;
      clusterCartesianCoordinates[kk * 4 + 1] = 9.0 + kk - 5;
      clusterCartesianCoordinates[kk * 4 + 2] = 17.0 + kk - 5;
      clusterCartesianCoordinates[kk * 4 + 3] = 0;
    }

    for (int kk = 11; kk < 15; kk++) {
      clusterCartesianCoordinates[kk * 4] = kk;
      clusterCartesianCoordinates[kk * 4 + 1] = 10.0;
      clusterCartesianCoordinates[kk * 4 + 2] = 19.0 + kk - 11;
      clusterCartesianCoordinates[kk * 4 + 3] = 0;
    }*/

  totalNumberOfClusters = 300;

  // printData(totalNumberOfClusters);

  // drawClusters(totalNumberOfClusters, dataFilename);
  drawCartesianClusters1D(totalNumberOfClusters);

  // Transform the cartesian coordinate system into the conformal mapping system
  conformalMapping(totalNumberOfClusters);

  drawConformalMappingClusters1D(totalNumberOfClusters);

  // Determine the maximum dimensions of the clusters for the accumulator
  determineMaxDimensions(totalNumberOfClusters);

  // Perform the hough transform on the clusters
  transform(totalNumberOfClusters);

  // Determine if any lines were found
  calculateTracks(totalNumberOfClusters);

  // Draw the reconstructed tracks
  drawTracks();

  // Draw the accumulator data as curves
  drawAccumulatorCurves(totalNumberOfClusters);

  // Draw the accumulator data as a histogram
  // drawAccumulatorHistogram(totalNumberOfClusters);

  return 0;
}
