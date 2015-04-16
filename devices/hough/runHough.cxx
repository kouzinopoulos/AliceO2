/// \file runHough.cxx
/// \brief Implementation of a cluster loader
/// \author Charis Kouzinopoulos

/// List of references: [1] Cheshkov, C. "Fast Hough-transform track reconstruction for the ALICE TPC." Nuclear
/// Instruments and Methods in Physics Research Section A: Accelerators, Spectrometers, Detectors and Associated
/// Equipment 566.1 (2006): 35-39.

#include "runHough.h"

void drawTracks(int totalNumberOfClusters)
{
  TCanvas* trackCanvas = new TCanvas("trackCanvas", "Reconstructed tracks", 0, 0, 800, 600);
  TGraph* trackGraph[totalNumberOfClusters];

  cout << "Number of tracks: " << getNumberOfTracks() << endl;

  for (int i = 0; i < getNumberOfTracks(); i++) {

    // DEBUG
    if (getTrackEtaSlice(i) != 79)
      continue;

    trackGraph[i] = new TGraph();

    trackGraph[i]->SetPoint(i * 2 + 0, getTrackAlpha1(i), getTrackBeta1(i));
    trackGraph[i]->SetPoint(i * 2 + 1, getTrackAlpha2(i), getTrackBeta2(i));

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
    Float_t x = getClusterX(i);
    Float_t y = getClusterY(i);

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
  TH2F* h = new TH2F("h", "Accumulator histogram", rResolution, 0, rResolution, thetaMax, 0, thetaMax);

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
    if (getClusterEtaSlice(i) == etaSlice) {
      cartesianClustersGraph1D->SetPoint(i, getClusterX(i), getClusterY(i));
    }
  }

  cartesianClustersGraph1D->SetMarkerStyle(7);
  cartesianClustersGraph1D->Draw("AP");
  cartesianClustersGraph1D->GetXaxis()->SetRangeUser(80, 240);
  cartesianClustersGraph1D->GetXaxis()->SetTitle("Pad Row");

  cartesianClustersCanvas1D->Print("cartesianClusters1D.pdf");
}

void drawConformalMappingClusters1D(int etaSlice, int totalNumberOfClusters)
{
  TCanvas* conformalMappingClustersCanvas1D =
    new TCanvas("conformalMappingClustersCanvas1D", "Conformal mapping clusters", 0, 0, 800, 600);
  TGraph* conformalMappingClustersGraph1D = new TGraph();

  for (Int_t i = 0; i < totalNumberOfClusters; i++) {
    if (getClusterEtaSlice(i) == etaSlice) {
      conformalMappingClustersGraph1D->SetPoint(i, getClusterAlpha(i), getClusterBeta(i));
    }
  }

  conformalMappingClustersGraph1D->SetMarkerStyle(7);
  conformalMappingClustersGraph1D->Draw("AP");
  conformalMappingClustersGraph1D->GetXaxis()->SetRangeUser(0.004, 0.012);

  conformalMappingClustersCanvas1D->Print("conformalMappingClusters1D.pdf");
}

void drawCartesianClusters(int totalNumberOfClusters)
{
  TCanvas* cartesianClustersCanvas = new TCanvas("cartesianClustersCanvas", "Cartesian clusters", 800, 600);
  TGraph2D* cartesianClustersGraph2D = new TGraph2D();

  for (Int_t i = 0; i < totalNumberOfClusters; i++) {
    // cout << "graph: " << i << " " << getClusterX(i) << " " << getClusterY(i) << " " <<
    // getClusterZ(i) << " " << endl;
    cartesianClustersGraph2D->SetPoint(i, getClusterX(i), getClusterY(i), getClusterZ(i));
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
    conformalMappingClustersGraph2D->SetPoint(i, getClusterAlpha(i), getClusterBeta(i), getClusterEta(i));
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
  // The track will have coordinates (a1,b1), (a2,b2)
  Float_t a1, b1, a2, b2;

  int trackCount = 0;

  for (Int_t rBin = 0; rBin < rResolution; rBin++) {
    for (Int_t theta = 0; theta < thetaMax; theta++) {
      if (getAccumulatorBin(etaSlice, rBin, theta) >= houghThreshold) {
        if (localAccumulatorMaxima(etaSlice, rBin, theta) > getAccumulatorBin(etaSlice, rBin, theta)) {
          continue;
        }

        if (theta >= 45 && theta <= 135) {
          // y = (r - x cos(t)) / sin(t)
          a1 = 0.0;
          b1 = (getRValue(rBin) - a1 * cos(theta * DEG2RAD)) / sin(theta * DEG2RAD);
          a2 = aMax;
          b2 = (getRValue(rBin) - a2 * cos(theta * DEG2RAD)) / sin(theta * DEG2RAD);
        } else {
          // x = (r - y sin(t)) / cos(t);
          b1 = 0.0;
          a1 = (getRValue(rBin) - b1 * sin(theta * DEG2RAD)) / cos(theta * DEG2RAD);
          b2 = bMax;
          a2 = (getRValue(rBin) - b2 * sin(theta * DEG2RAD)) / cos(theta * DEG2RAD);
        }
        // cout << "Found a track at [" << etaSlice << "][" << getRValue(rBin) << "][" << theta << "] with a peak of "
        // << getAccumulatorBin(etaSlice, rBin, theta) << " points"
        //     << endl << "Track coordinates: " << a1 << "," << b1 << " - " << a2 << "," << b2 << endl;
        // cout << "xMax: " << xMax << " yMax: " << yMax << endl;
        trackCount++;

        setTrackParameters(etaSlice, a1, b1, a2, b2);
      }
    }
  }
  if (trackCount > 0) {
    cout << "rMax " << rMax << " rResolution " << rResolution << endl;
    cout << "Found " << trackCount << " tracks for η slice " << etaSlice << endl;
  }
}

/// Calculate an approximate value for η. See [1]:p8 for more information. Values below taken from
/// AliHLTConfMapPoint.cxx
void calculateEta(int totalNumberOfClusters)
{
  for (int i = 0; i < totalNumberOfClusters; i++) {
    Double_t radial =
      sqrt(getClusterX(i) * getClusterX(i) + getClusterY(i) * getClusterY(i) + getClusterZ(i) * getClusterZ(i));
    Double_t eta = 0.5 * log((radial + getClusterZ(i)) / (radial - getClusterZ(i)));

    // Store eta to be later used by the conformalMapping method
    setClusterEta(i, eta);
  }
}

/// Discretize the eta values of all clusters into etaResolution bins
void calculateEtaSlice(int totalNumberOfClusters)
{
  for (int i = 0; i < totalNumberOfClusters; i++) {
    Float_t eta = getClusterEta(i);

    if (etaMax - etaMin == 0) {
      cerr << "The minimum and maximum eta value of all clusters is identical" << endl;
      exit(1);
    }

    Double_t etaSlice = (etaResolution * (eta - etaMin)) / (etaMax - etaMin);

    setClusterEtaSlice(i, (Int_t)etaSlice);
  }
}

void conformalMapping(int totalNumberOfClusters)
{
  for (int i = 0; i < totalNumberOfClusters; i++) {
    Float_t x = getClusterX(i);
    Float_t y = getClusterY(i);

    // Equation (2) from paper [1]
    Float_t alpha = x / (x * x + y * y);
    Float_t beta = y / (x * x + y * y);

    // cout << "i: " << i << " X: " << x << " Y: " << y << " A: " << alpha << " B: " << beta << " η: " << eta << " η
    // slice: " << etaSlice
    //     << endl;

    setClusterAlpha(i, alpha);
    setClusterBeta(i, beta);
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
    Float_t x = getClusterX(i);
    Float_t y = getClusterY(i);

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

  // cout << "rMin: " << -rMax << " rMax: " << rMax << endl;

  // Reserve space for the accumulator bins
  accumulator.resize(etaResolution * rResolution * thetaMax + rResolution * thetaMax + thetaMax, 0);

  cout << "Accumulator size: " << 2 * rMax* thetaMax* rResolution* etaResolution << endl;
  cout << "Alternative accumulator size: " << etaResolution* rResolution* thetaMax + rResolution* thetaMax + thetaMax
       << endl;

  for (Int_t i = 0; i < totalNumberOfClusters; i++) {
    Float_t a = getClusterAlpha(i);
    Float_t b = getClusterBeta(i);
    Float_t eta = getClusterEta(i);
    UInt_t etaSlice = getClusterEtaSlice(i);

    for (Int_t theta = 0; theta < thetaMax; theta++) {
      double r = (a * cos(theta * DEG2RAD)) + (b * sin(theta * DEG2RAD));
      int rBin = getRBinValue(r);

      if (etaSlice == 20) {
        cout << "Cluster: " << i << " a: " << a << " b: " << b << " theta: " << theta << " rMin: " << -rMax
             << " r: " << r << " rMax: " << rMax << " rBin: " << rBin << " eta: " << eta << " etaSlice: " << etaSlice
             << endl;
      }

      // Use a discreet value for r by rounding r. Then, rMax is added to r to change its range from -rMax <= r <= rMax
      // to 0 <= r <= 2 * rMax
      setAccumulatorBin(etaSlice, rBin, theta);
    }
  }
}

void houghTransform(int totalNumberOfClusters)
{
  // Iterate through all the pad rows for each slice of a TPC partition/patch
  for (UInt_t partition = 0; partition < 6; partition++) {
    for (UInt_t slice = 0; slice < 36; slice++) {
      for (UChar_t row = getFirstPadRow(partition); row <= getLastPadRow(partition); row++) {
      }
    }
  }
}

/// Determine the minimum and maximum values of the pseudorapidity (eta). That way, the TPC digits can be transformed in
/// two dimensions instead of three in slices of similar pseudorapidity
void determineMinMaxEta(int totalNumberOfClusters)
{
  for (Int_t i = 0; i < totalNumberOfClusters; i++) {

    Float_t eta = getClusterEta(i);

    // Set initial values to etaMin and etaMax
    if (i == 0) {
      etaMin = eta;
      etaMax = eta;
      continue;
    }

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
    if (abs(getClusterX(i)) > xMax) {
      xMax = ceil(abs(getClusterX(i)));
    }
    if (abs(getClusterY(i)) > yMax) {
      yMax = ceil(abs(getClusterY(i)));
    }
    if (abs(getClusterZ(i)) > zMax) {
      zMax = ceil(abs(getClusterZ(i)));
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
      aMin = floor(getClusterAlpha(i));
      aMax = ceil(abs(getClusterAlpha(i)));
      bMin = floor(getClusterBeta(i));
      bMax = ceil(abs(getClusterBeta(i)));
      continue;
    }
    if (getClusterAlpha(i) < aMin) {
      aMin = floor(getClusterAlpha(i));
    } else if (abs(getClusterAlpha(i)) > aMax) {
      aMax = ceil(abs(getClusterAlpha(i)));
    }

    if (getClusterBeta(i) < bMin) {
      bMin = floor(getClusterBeta(i));
    } else if (abs(getClusterBeta(i)) > bMax) {
      bMax = ceil(abs(getClusterBeta(i)));
    }
  }
  cout << "aMax: " << aMax << " bMax: " << bMax << endl;
}

void printData(int totalNumberOfClusters)
{
  cout << "Cluster ID" << setw(13) << "X coordinate" << setw(13) << "Y coordinate" << setw(13) << "Z coordinate"
       << endl;

  for (int i = 0; i < totalNumberOfClusters; i++) {
    cout << (AliHLTUInt32_t)getClusterID(i) << setw(13) << getClusterX(i) << setw(13) << getClusterY(i) << setw(13)
         << getClusterZ(i) << endl;
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

    setClusterParameters(clusterID, spacepoints->GetX(clusterID), spacepoints->GetY(clusterID),
                         spacepoints->GetZ(clusterID), spacepoints->GetCharge(clusterID), currentSlice,
                         currentPartition);
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

  cout << "Added " << totalNumberOfClusters << " clusters from " << totalNumberOfDataFiles
       << " data files. Cluster vector size: " << clusterData.size() << endl;

  if (totalNumberOfClusters != clusterData.size()) {
    std::cerr << "The cluster vector size does not match the reported number of clusters" << endl;
    exit(1);
  }

  // DEBUG
  calculateEta(totalNumberOfClusters);
  determineMinMaxEta(totalNumberOfClusters);
  calculateEtaSlice(totalNumberOfClusters);
  drawCartesianClusters1D(15, totalNumberOfClusters);
  conformalMapping(totalNumberOfClusters);
  drawConformalMappingClusters1D(15, totalNumberOfClusters);

  // DEBUG
  /*  totalNumberOfClusters = 100;
    for (int kk = 0; kk < 6; kk++) {
      setClusterParameters(kk, (float)9.0 + kk, (float)17.0, 1);
    }

    for (int kk = 6; kk < 11; kk++) {
      setClusterParameters(kk, 9.0 + kk - 5, 17.0 + kk - 5, 1);
    }

    for (int kk = 11; kk < 150; kk++) {
      setClusterParameters(kk, 10, 19.0 + kk - 11, 1);
    }
  */
  /*  totalNumberOfClusters = 20000;
    // printData(totalNumberOfClusters);

    // Allocate space for the conformal mapping parameter vector
    clusterConformalMappingCoordinates.resize(totalNumberOfClusters * clusterConformalMappingParameters);

    calculateEta(totalNumberOfClusters);

    // Determine the minimum and maximum values of eta. That way the TPC digits can be grouped into pseudorapidity bins
    determineMinMaxEta(totalNumberOfClusters);

    // Discretize the eta values of all clusters into etaResolution bins
    calculateEtaSlice(totalNumberOfClusters);

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
    //for (int etaSlice = 0; etaSlice <= etaResolution; etaSlice++) {
    //  trackFinding(etaSlice);
    }
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
      array[getClusterEtaSlice(i)]++;
    }

    for ( i = 0; i <= etaResolution; i++) {
      cout << i << " " << array[i] << endl;
    }

    for (int theta = 0; theta < thetaMax; theta++) {
      for (int rBin = 0; rBin < rResolution; rBin++) {
        cout << theta << " " << rBin << " " << getAccumulatorBin(15, rBin, theta ) << endl;
      }
    }

    return 0;

    // Clear the track coordinates so the vector can be used for the conformal mapping tracker
    trackCoordinates.clear();

    // Draw the accumulator data as curves
    drawAccumulatorCurves(totalNumberOfClusters);
  */
  return 0;
}
