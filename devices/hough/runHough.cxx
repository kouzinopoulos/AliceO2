/// \file runHough.cxx
/// \brief Implementation of a cluster loader
/// \author Charis Kouzinopoulos

/// List of references: [1] Cheshkov, C. "Fast Hough-transform track reconstruction for the ALICE TPC." Nuclear
/// Instruments and Methods in Physics Research Section A: Accelerators, Spectrometers, Detectors and Associated
/// Equipment 566.1 (2006): 35-39.

#include "runHough.h"

using namespace AliceO2::Hough;

// fgBeta1 and fgBeta2 are two curves which define the Hough space
Float_t fgBeta1 = 1.0 / Row2X(84);
Float_t fgBeta2 = 1.0 / (Row2X(158) * (1.0 + tan(Pi() * 10 / 180) * tan(Pi() * 10 / 180)));

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
    Float_t x = clusterCollection->getClusterX(i);
    Float_t y = clusterCollection->getClusterY(i);

    g[i] = new TGraph();

    for (Int_t theta = 0; theta < thetalphaMax; theta++) {
      double r = (x * cos(theta * DEG2RAD)) + (y * sin(theta * DEG2RAD));
      g[i]->SetPoint(i * thetalphaMax + theta, theta, r);
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
  TH2F* h = new TH2F("h", "Accumulator histogram", rResolution, 0, rResolution, thetalphaMax, 0, thetalphaMax);

  h->SetFillColor(46);

  for (Int_t r = 0; r < rResolution; r++) {
    for (Int_t theta = 0; theta < thetalphaMax; theta++) {
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
    if (clusterCollection->getClusterEtaSlice(i) == etaSlice) {
      cartesianClustersGraph1D->SetPoint(i, clusterCollection->getClusterX(i), clusterCollection->getClusterY(i));
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
    if (clusterCollection->getClusterEtaSlice(i) == etaSlice) {
      conformalMappingClustersGraph1D->SetPoint(i, clusterCollection->getClusterAlpha(i),
                                                clusterCollection->getClusterBeta(i));
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
    cartesianClustersGraph2D->SetPoint(i, clusterCollection->getClusterX(i), clusterCollection->getClusterY(i),
                                       clusterCollection->getClusterZ(i));
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
    conformalMappingClustersGraph2D->SetPoint(i, clusterCollection->getClusterAlpha(i),
                                              clusterCollection->getClusterBeta(i),
                                              clusterCollection->getClusterEta(i));
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
      if ((deltaR + r >= 0 && deltaR + r < rMax) && (deltaTheta + theta >= 0 && deltaTheta + theta < thetalphaMax)) {
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
    for (Int_t theta = 0; theta < thetalphaMax; theta++) {
      if (getAccumulatorBin(etaSlice, rBin, theta) >= houghThreshold) {
        if (localAccumulatorMaxima(etaSlice, rBin, theta) > getAccumulatorBin(etaSlice, rBin, theta)) {
          continue;
        }

        if (theta >= 45 && theta <= 135) {
          // y = (r - x cos(t)) / sin(t)
          a1 = 0.0;
          b1 = (getRValue(rBin) - a1 * cos(theta * DEG2RAD)) / sin(theta * DEG2RAD);
          a2 = alphaMax;
          b2 = (getRValue(rBin) - a2 * cos(theta * DEG2RAD)) / sin(theta * DEG2RAD);
        } else {
          // x = (r - y sin(t)) / cos(t);
          b1 = 0.0;
          a1 = (getRValue(rBin) - b1 * sin(theta * DEG2RAD)) / cos(theta * DEG2RAD);
          b2 = betaMax;
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
    Double_t radial = sqrt(clusterCollection->getClusterX(i) * clusterCollection->getClusterX(i) +
                           clusterCollection->getClusterY(i) * clusterCollection->getClusterY(i) +
                           clusterCollection->getClusterZ(i) * clusterCollection->getClusterZ(i));
    Double_t eta =
      0.5 * log((radial + clusterCollection->getClusterZ(i)) / (radial - clusterCollection->getClusterZ(i)));

    // Store eta to be later used by the conformalMapping method
    clusterCollection->setClusterEta(i, eta);
  }
}

/// Discretize the eta values of all clusters into etaResolution bins
void calculateEtaSlice(int totalNumberOfClusters)
{
  for (int i = 0; i < totalNumberOfClusters; i++) {
    Float_t eta = clusterCollection->getClusterEta(i);

    if (etalphaMax - etalphaMin == 0) {
      cerr << "The minimum and maximum eta value of all clusters is identical" << endl;
      exit(1);
    }

    Double_t etaSlice = (etaResolution * (eta - etalphaMin)) / (etalphaMax - etalphaMin);

    clusterCollection->setClusterEtaSlice(i, (Int_t)etaSlice);
  }
}

void conformalMapping(int totalNumberOfClusters)
{
  for (int i = 0; i < totalNumberOfClusters; i++) {
    Float_t x = clusterCollection->getClusterX(i);
    Float_t y = clusterCollection->getClusterY(i);

    // Equation (2) from paper [1]
    Float_t alpha = x / (x * x + y * y);
    Float_t beta = y / (x * x + y * y);

    // cout << "i: " << i << " X: " << x << " Y: " << y << " A: " << alpha << " B: " << beta << " η: " << eta << " η
    // slice: " << etaSlice
    //     << endl;

    clusterCollection->setClusterAlpha(i, alpha);
    clusterCollection->setClusterBeta(i, beta);
  }
}

void conformalMapping2()
{
  for (UInt_t padRow = 0; padRow < Transform::GetNRows(); padRow++) {
    for (UInt_t clusterNumber = 0; clusterNumber < clusterCollection->getNumberOfClustersPerPadRow(padRow);
         clusterNumber++) {
      Float_t x = clusterCollection->getClusterX(padRow, clusterNumber);
      Float_t y = clusterCollection->getClusterY(padRow, clusterNumber);

      // Equation (2) from paper [1]
      Float_t alpha = x / (x * x + y * y);
      Float_t beta = y / (x * x + y * y);

      clusterCollection->setClusterAlpha(padRow, clusterNumber, alpha);
      clusterCollection->setClusterBeta(padRow, clusterNumber, beta);
    }
  }
}

void transformCartesian(int totalNumberOfClusters)
{
  thetalphaMax = thetaResolution;
  // Trigonometrically, the maximum distance is designated by the square root of the summation of the squares of the x
  // and y dimensions
  rMax = ceil(sqrt(xMax * xMax + yMax * yMax));

  // The lines will have -rMax <= r <= rMax and 0 <= theta <= thetalphaMax. The total space needed is thus 2 * rMax *
  // thetalphaMax
  accumulator.resize(2 * rMax * thetalphaMax * rResolution, 0);

  for (Int_t i = 0; i < totalNumberOfClusters; i++) {
    Float_t x = clusterCollection->getClusterX(i);
    Float_t y = clusterCollection->getClusterY(i);

    for (Int_t theta = 0; theta < thetalphaMax; theta++) {
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
  thetalphaMax = 180;
  // Trigonometrically, the maximum distance is designated by the square root of the summation of the squares of the x
  // and y dimensions
  rMax = ceil(sqrt(alphaMax * alphaMax + betaMax * betaMax));

  // cout << "rMin: " << -rMax << " rMax: " << rMax << endl;

  // Reserve space for the accumulator bins
  accumulator.resize(etaResolution * rResolution * thetalphaMax + rResolution * thetalphaMax + thetalphaMax, 0);

  cout << "Accumulator size: " << 2 * rMax* thetalphaMax* rResolution* etaResolution << endl;
  cout << "Alternative accumulator size: "
       << etaResolution* rResolution* thetalphaMax + rResolution* thetalphaMax + thetalphaMax << endl;

  for (Int_t i = 0; i < totalNumberOfClusters; i++) {
    Float_t a = clusterCollection->getClusterAlpha(i);
    Float_t b = clusterCollection->getClusterBeta(i);
    Float_t eta = clusterCollection->getClusterEta(i);
    UInt_t etaSlice = clusterCollection->getClusterEtaSlice(i);

    for (Int_t theta = 0; theta < thetalphaMax; theta++) {
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

/// Determine the minimum and maximum values of the pseudorapidity (eta). That way, the TPC digits can be
/// transformed in
/// two dimensions instead of three in slices of similar pseudorapidity
void determineMinMaxEta(int totalNumberOfClusters)
{
  for (Int_t i = 0; i < totalNumberOfClusters; i++) {

    Float_t eta = clusterCollection->getClusterEta(i);

    // Set initial values to etalphaMin and etalphaMax
    if (i == 0) {
      etalphaMin = eta;
      etalphaMax = eta;
      continue;
    }

    if (eta < etalphaMin) {
      etalphaMin = eta;
    } else if (eta > etalphaMax) {
      etalphaMax = eta;
    }
  }
  cout << "Minimum eta: " << etalphaMin << " Maximum eta: " << etalphaMax << endl;
}

/// Determine the maximum ceiling to x,y,z coordinates from the clusterCartesianCoordinates vector to
/// later allocate the
/// hough transform accumulator
void determineMinMaxXY(int totalNumberOfClusters)
{
  for (Int_t i = 0; i < totalNumberOfClusters; i++) {
    if (i == 0) {
      xMin = clusterCollection->getClusterX(i);
      xMax = clusterCollection->getClusterX(i);
      yMin = clusterCollection->getClusterY(i);
      yMax = clusterCollection->getClusterY(i);
      continue;
    }
    if (clusterCollection->getClusterX(i) < xMin) {
      xMin = clusterCollection->getClusterX(i);
    } else if (clusterCollection->getClusterX(i) > xMax) {
      xMax = clusterCollection->getClusterX(i);
    }

    if (clusterCollection->getClusterY(i) < yMin) {
      yMin = clusterCollection->getClusterY(i);
    } else if (clusterCollection->getClusterY(i) > yMax) {
      yMax = clusterCollection->getClusterY(i);
    }
  }
  cout << "xMin: " << xMin << " xMax: " << xMax << " yMin: " << yMin << " yMax: " << yMax << endl;
}

/// Determine the maximum ceiling to a and b conformal mapping coordinates from the
/// clusterCartesianCoordinates vector
/// to later allocate the hough transform accumulator
void determineMinMaxAlphaBeta(int totalNumberOfClusters)
{
  for (Int_t i = 0; i < totalNumberOfClusters; i++) {
    if (i == 0) {
      alphaMin = clusterCollection->getClusterX(i);
      alphaMax = clusterCollection->getClusterX(i);
      betaMin = clusterCollection->getClusterY(i);
      betaMax = clusterCollection->getClusterY(i);
      continue;
    }
    if (clusterCollection->getClusterAlpha(i) < alphaMin) {
      alphaMin = clusterCollection->getClusterAlpha(i);
    } else if (clusterCollection->getClusterX(i) > alphaMax) {
      alphaMax = clusterCollection->getClusterAlpha(i);
    }

    if (clusterCollection->getClusterY(i) < betaMin) {
      betaMin = clusterCollection->getClusterBeta(i);
    } else if (clusterCollection->getClusterBeta(i) > betaMax) {
      betaMax = clusterCollection->getClusterBeta(i);
    }
  }
  cout << "alphaMin: " << alphaMin << " alphaMax: " << alphaMax << " betaMin: " << betaMin << " betaMax: " << betaMax
       << endl;
}

void printData(int totalNumberOfClusters)
{
  cout << "ID" << setw(8) << "X" << setw(8) << "Y" << setw(8) << "Z" << setw(9) << "Pad Row" << setw(10) << "Pad"
       << setw(10) << "Time" << setw(10) << "Charge" << endl;

  for (int i = 0; i < totalNumberOfClusters; i++) {
    cout << (AliHLTUInt32_t)clusterCollection->getClusterID(i) << setw(13) << clusterCollection->getClusterX(i)
         << setw(13) << clusterCollection->getClusterY(i) << setw(13) << clusterCollection->getClusterZ(i) << setw(8)
         << (UInt_t)clusterCollection->getClusterPadRow(i) << setw(8) << clusterCollection->getClusterPad(i) << setw(8)
         << clusterCollection->getClusterTime(i) << setw(7) << clusterCollection->getClusterCharge(i) << endl;
  }
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

  clusterCollection = new ClusterCollection();

  // Traverse the filesystem and execute processData for each cluster file found
  if (boost::filesystem::exists(dataPath) && boost::filesystem::is_directory(dataPath)) {
    for (boost::filesystem::directory_iterator directoryIterator(dataPath); directoryIterator != endIterator;
         ++directoryIterator) {
      if (boost::filesystem::is_regular_file(directoryIterator->status())) {
        totalNumberOfClusters +=
          clusterCollection->processData(directoryIterator->path().string(), dataType, dataOrigin);
        totalNumberOfDataFiles++;
      }
    }
  } else {
    std::cerr << "Path " << dataPath.string() << "/ could not be found or does not contain any valid data files"
              << endl;
    exit(1);
  }

  cout << "Added " << totalNumberOfClusters << " clusters from " << totalNumberOfDataFiles << " data files." << endl;

  for (int i = 0; i < Transform::GetNRows(); i++) {
    cout << clusterCollection->getNumberOfClustersPerPadRow(i) << " Clusters for PadRow " << i << endl;
  }

  // DEBUG
  determineMinMaxXY(totalNumberOfClusters);
  calculateEta(totalNumberOfClusters);
  determineMinMaxEta(totalNumberOfClusters);
  calculateEtaSlice(totalNumberOfClusters);
  drawCartesianClusters1D(15, totalNumberOfClusters);
  conformalMapping(totalNumberOfClusters);
  conformalMapping2();
  determineMinMaxAlphaBeta(totalNumberOfClusters);
  drawConformalMappingClusters1D(15, totalNumberOfClusters);
  drawCartesianClusters(totalNumberOfClusters);

  // printData(totalNumberOfClusters);

  Float_t ptmin = 0.1 * Transform::GetSolenoidField();
  // FIXME: find the correct value for zvertex
  // Float_t zvertex = GetZ();
  Float_t zvertex = 0;

  Hough* hough = new Hough();

  hough->SetThreshold(4);
  hough->CalcTransformerParams(ptmin);
  hough->SetPeakThreshold(70, -1);
  hough->Init("./", kFALSE, 100, kFALSE, 4, 0, 0, zvertex);
  hough->SetAddHistograms();

  for (Int_t slice = 0; slice <= 35; slice++) {
    //     hough->ReadData(slice,iEvent);
    hough->Transform(0, clusterCollection);
    hough->AddAllHistogramsRows();
    hough->FindTrackCandidatesRow();
    //    hough->AddTracks();
  }

  /*
    for (UInt_t patch = 0; patch < AliceO2::Hough::Transform::GetNPatches(); patch++) {

      fHoughTransformer[patch] = new AliceO2::Hough::TransformerRow(0, patch, fNEtaSegments);
      fHoughTransformer[patch]->SetLastTransformer(lasttransformer);
      lasttransformer = fHoughTransformer[patch];

      fHoughTransformer[patch]->CreateHistograms()
    }
  */

  /*
    for (UInt_t slice = 0; slice < AliceO2::Hough::Transform::GetNSlice(); slice++) {
      for (UInt_t patch = 0; patch < AliceO2::Hough::Transform::GetNPatches(); patch++) {
        AliceO2::Hough::TransformerRow row(slice, patch, etaResolution);

        // Allocate space for the accumulator and determine initial values. This must be done before executing the
    TransformCircle method
        row.CreateHistograms(rResolution, xMin, xMax, thetaResolution, yMin, yMax);

        // Perform the actual transformation
        row.TransformCircle(clusterCollection);

        cout << "Number of Eta Segments: " << row.GetNEtaSegments() << " - etaResolution: " << etaResolution << endl;

        for (Int_t segment = 0; segment < row.GetNEtaSegments(); segment++) {
          // Get a pointer for the accumulator for each eta slice
          AliceO2::Hough::Accumulator* accumulatorPointer = row.GetHistogram(segment);

          if (accumulatorPointer->GetNEntries() > 0) {
            cout << "Number of entries in the Accumulator for eta segment " << segment << ": "
                 << accumulatorPointer->GetNEntries() << endl;
            exit(0);
          }
        }
        // Resets all the histograms. Should be done when processing new slice
        row.Reset();

        exit(0);
      }
    }
  */
  /*
    AliceO2::Hough::TransformerRow row(0, 0, etaResolution);

    // Allocate space for the accumulator. This must be done before executing the TransformCircle method
    row.CreateHistograms(rResolution, xMin, xMax, thetaResolution, yMin, yMax);

  //  row.SetLowerThreshold(3);
  //  row.SetUpperThreshold(5);

    // Performs the actual transformation
    row.TransformCircle(clusterCollection);

    cout << "Number of Eta Segments: " << row.GetNEtaSegments() << " - etaResolution: " << etaResolution << endl;

    for ( Int_t segment = 0; segment < row.GetNEtaSegments(); segment++) {
      //Get a pointer for the accumulator for each eta slice
      AliceO2::Hough::Accumulator* accumulatorPointer = row.GetHistogram(segment);

      if (accumulatorPointer->GetNEntries() > 0) {
        cout << "Number of entries in the Accumulator for eta segment "<< segment << ": " <<
  accumulatorPointer->GetNEntries() << endl;
        exit(0);
      }

      //cout << "Number of entries in the Accumulator for eta segment "<< segment << ": " <<
  accumulatorPointer->GetNEntries() << endl;

      //FIXME: The accumulator array can be accessed directly through here:
      //int* GetContentArray() const { return fContent; }
    }

    // Resets all the histograms. Should be done when processing new slice
    row.Reset();
  */
  //  AliceO2::Hough::Accumulator** test = new AliceO2::Hough::Accumulator* [5];

  // AliceO2::Hough::Accumulator test2;
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

    // Determine the minimum and maximum values of eta. That way the TPC digits can be grouped into
    pseudorapidity bins
    determineMinMaxEta(totalNumberOfClusters);

    // Discretize the eta values of all clusters into etaResolution bins
    calculateEtaSlice(totalNumberOfClusters);

    // Convert the TPC cluster coordinates from the cartesian to the conformal mapping system
    conformalMapping(totalNumberOfClusters);

    // drawConformalMappingClusters1D(totalNumberOfClusters);

    // Determine the maximum dimensions of the clusters for the accumulator
    determineMinMaxXY(totalNumberOfClusters);

    // Perform the hough transform on the TPC clusters for the cartesian system
    // transformCartesian(totalNumberOfClusters);

    // Determine the maximum dimensions of the clusters for the accumulator
    determineMinMaxAlphaBeta(totalNumberOfClusters);

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

    for (int theta = 0; theta < thetalphaMax; theta++) {
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
  delete clusterCollection;

  return 0;
}
