/// \file runHough.cxx
/// \brief Implementation of a cluster loader
/// \author Charis Kouzinopoulos

/// List of references: [1] Cheshkov, C. "Fast Hough-transform track reconstruction for the ALICE TPC." Nuclear
/// Instruments and Methods in Physics Research Section A: Accelerators, Spectrometers, Detectors and Associated
/// Equipment 566.1 (2006): 35-39.

#include "runHough.h"

using namespace AliceO2::Hough;

void draw1DCartesianClustersPerEtaSlice(int etaSlice)
{
  TCanvas* cartesianClustersPerEtaSliceCanvas1D =
    new TCanvas("cartesianClustersPerEtaSliceCanvas1D", "Cartesian clusters per eta slice");
  TGraph* cartesianClustersPerEtaSliceGraph1D = new TGraph();

  for (UInt_t padRow = 0; padRow < Transform::GetNRows(); padRow++) {
    for (UInt_t clusterNumber = 0; clusterNumber < clusterCollection->getNumberOfClustersPerPadRow(padRow);
         clusterNumber++) {
      if (clusterCollection->getClusterEtaSlice(padRow, clusterNumber) == etaSlice) {
        cartesianClustersPerEtaSliceGraph1D->SetPoint(clusterNumber,
                                                      clusterCollection->getClusterX(padRow, clusterNumber),
                                                      clusterCollection->getClusterY(padRow, clusterNumber));
      }
    }
  }

  cartesianClustersPerEtaSliceGraph1D->SetMarkerStyle(7);
  cartesianClustersPerEtaSliceGraph1D->Draw("AP");
  cartesianClustersPerEtaSliceGraph1D->GetXaxis()->SetRangeUser(80, 240);
  cartesianClustersPerEtaSliceGraph1D->GetXaxis()->SetTitle("Pad Row");

  cartesianClustersPerEtaSliceCanvas1D->Print("cartesianClustersPerEtaSlice.pdf");
}

void draw1DCartesianClustersPerPadRow(int padRow)
{
  TCanvas* cartesianClustersPerPadRowCanvas1D =
    new TCanvas("cartesianClustersPerPadRow1D", "Cartesian clusters per pad row");
  TGraph* cartesianClustersPerPadRowGraph1D = new TGraph();

  for (UInt_t clusterNumber = 0; clusterNumber < clusterCollection->getNumberOfClustersPerPadRow(padRow);
       clusterNumber++) {
    cartesianClustersPerPadRowGraph1D->SetPoint(clusterNumber, clusterCollection->getClusterX(padRow, clusterNumber),
                                                clusterCollection->getClusterY(padRow, clusterNumber));
  }

  cartesianClustersPerPadRowGraph1D->SetMarkerStyle(7);
  cartesianClustersPerPadRowGraph1D->Draw("AP");
  cartesianClustersPerPadRowGraph1D->GetXaxis()->SetRangeUser(80, 240);
  cartesianClustersPerPadRowGraph1D->GetXaxis()->SetTitle("Pad Row");

  cartesianClustersPerPadRowCanvas1D->Print("cartesianClustersPerPadRow.pdf");
}

void drawCartesianClusters()
{
  TCanvas* cartesianClustersCanvas = new TCanvas("cartesianClustersCanvas", "Cartesian clusters", 800, 600);
  TGraph2D* cartesianClustersGraph2D = new TGraph2D();

  for (UInt_t padRow = 0; padRow < Transform::GetNRows(); padRow++) {
    for (UInt_t clusterNumber = 0; clusterNumber < clusterCollection->getNumberOfClustersPerPadRow(padRow);
         clusterNumber++) {
      cartesianClustersGraph2D->SetPoint(padRow, clusterCollection->getClusterX(padRow, clusterNumber),
                                         clusterCollection->getClusterY(padRow, clusterNumber),
                                         clusterCollection->getClusterZ(padRow, clusterNumber));
    }
  }

  // Draw with colored dots
  cartesianClustersGraph2D->SetMarkerStyle(6);
  cartesianClustersGraph2D->SetTitle("TPC Clusters in x, y, z");
  cartesianClustersGraph2D->Draw("pcol");

  cartesianClustersCanvas->Print("cartesianClusters2D.pdf");
}

/// Calculate an approximate value for η. See [1]:p8 for more information. Values below taken from
/// AliHLTConfMapPoint.cxx
void calculateEta(int totalNumberOfClusters)
{
  for (UInt_t padRow = 0; padRow < Transform::GetNRows(); padRow++) {
    for (UInt_t clusterNumber = 0; clusterNumber < clusterCollection->getNumberOfClustersPerPadRow(padRow);
         clusterNumber++) {
      Double_t radial = sqrt(
        clusterCollection->getClusterX(padRow, clusterNumber) * clusterCollection->getClusterX(padRow, clusterNumber) +
        clusterCollection->getClusterY(padRow, clusterNumber) * clusterCollection->getClusterY(padRow, clusterNumber) +
        clusterCollection->getClusterZ(padRow, clusterNumber) * clusterCollection->getClusterZ(padRow, clusterNumber));
      Double_t eta = 0.5 * log((radial + clusterCollection->getClusterZ(padRow, clusterNumber)) /
                               (radial - clusterCollection->getClusterZ(padRow, clusterNumber)));

      // Store eta to the cluster header
      clusterCollection->setClusterEta(padRow, clusterNumber, eta);
    }
  }
}

/// Determine the minimum and maximum values of the pseudorapidity (eta). That way, the TPC digits can be
/// transformed in two dimensions instead of three in slices of similar pseudorapidity
void determineMinMaxEta(int totalNumberOfClusters)
{
  for (UInt_t padRow = 0; padRow < Transform::GetNRows(); padRow++) {
    for (UInt_t clusterNumber = 0; clusterNumber < clusterCollection->getNumberOfClustersPerPadRow(padRow);
         clusterNumber++) {
      Double_t eta = clusterCollection->getClusterEta(padRow, clusterNumber);

      // Set initial values to etaMin and etaMax
      if (padRow == 0 && clusterNumber == 0) {
        etaMin = eta;
        etaMax = eta;
        continue;
      }

      if (eta < etaMin) {
        etaMin = eta;
      }
      if (eta > etaMax) {
        etaMax = eta;
      }
    }
  }
  cout << "Minimum eta: " << etaMin << " Maximum eta: " << etaMax << endl;
}

/// Discretize the eta values of all clusters into etaResolution bins
void calculateEtaSlice(int totalNumberOfClusters)
{
  for (UInt_t padRow = 0; padRow < Transform::GetNRows(); padRow++) {
    for (UInt_t clusterNumber = 0; clusterNumber < clusterCollection->getNumberOfClustersPerPadRow(padRow);
         clusterNumber++) {
      Double_t eta = clusterCollection->getClusterEta(padRow, clusterNumber);

      if (etaMax - etaMin == 0) {
        cerr << "The minimum and maximum eta value of all clusters are identical" << endl;
        exit(1);
      }

      Double_t etaSlice = (etaSlices * (eta - etaMin)) / (etaMax - etaMin);

      clusterCollection->setClusterEtaSlice(padRow, clusterNumber, (Int_t)etaSlice);
    }
  }
}

void printClusterInformation(int padRow)
{
  cout << "ID" << setw(16) << "X" << setw(14) << "Y" << setw(12) << "Z" << setw(9) << "Pad" << setw(10) << "Time"
       << setw(10) << "Charge" << endl;

  for (UInt_t clusterNumber = 0; clusterNumber < clusterCollection->getNumberOfClustersPerPadRow(padRow);
       clusterNumber++) {
    cout << (AliHLTUInt32_t)clusterCollection->getClusterID(padRow, clusterNumber) << setw(13)
         << clusterCollection->getClusterX(padRow, clusterNumber) << setw(13)
         << clusterCollection->getClusterY(padRow, clusterNumber) << setw(13)
         << clusterCollection->getClusterZ(padRow, clusterNumber) << setw(8)
         << clusterCollection->getClusterPad(padRow, clusterNumber) << setw(8)
         << clusterCollection->getClusterTime(padRow, clusterNumber) << setw(7)
         << clusterCollection->getClusterCharge(padRow, clusterNumber) << endl;
  }
}

int main(int argc, char** argv)
{

  boost::program_options::options_description desc("Allowed options");
  desc.add_options()("help", "Display this help and exit")(
    "event", boost::program_options::value<std::string>(),
    "Specify an event to load inside the emulated-tpc-clusters directory <mandatory>")(
    "print", boost::program_options::value<int>(), "Print cluster information for a specific pad row and exit");

  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
  boost::program_options::notify(vm);

  if (vm.count("help")) {
    cout << desc << endl;
    return 1;
  }

  if (!vm.count("event")) {
    cout << desc << endl;
    cout << "Please provide an event number!" << endl;
    return 1;
  }

  // Create data path
  std::string dataFilename = "emulated-tpc-clusters/event";
  dataFilename += vm["event"].as<std::string>();

  boost::filesystem::path dataPath(dataFilename);
  boost::filesystem::directory_iterator endIterator;

  typedef std::multimap<std::time_t, boost::filesystem::path> result_set_t;
  result_set_t result_set;

  std::string dataType = "CLUSTERS", dataOrigin = "TPC ";

  int totalNumberOfClusters = 0, totalNumberOfDataFiles = 0;

  clusterCollection = new ClusterCollection();

  // Traverse the filesystem and call readData for each cluster file found
  if (boost::filesystem::exists(dataPath) && boost::filesystem::is_directory(dataPath)) {
    for (boost::filesystem::directory_iterator directoryIterator(dataPath); directoryIterator != endIterator;
         ++directoryIterator) {
      if (boost::filesystem::is_regular_file(directoryIterator->status())) {
        totalNumberOfClusters += clusterCollection->readData(directoryIterator->path().string(), dataType, dataOrigin);
        totalNumberOfDataFiles++;
      }
    }
  } else {
    std::cerr << "Path " << dataPath.string() << "/ could not be found or does not contain any valid data files"
              << endl;
    exit(1);
  }

  cout << "Added " << totalNumberOfClusters << " clusters from " << totalNumberOfDataFiles << " data files." << endl;

  if (vm.count("print")) {
    printClusterInformation(vm["print"].as<int>());
    return 1;
  }

  // for (int i = 0; i < Transform::GetNRows(); i++) {
  //  cout << clusterCollection->getNumberOfClustersPerPadRow(i) << " Clusters for PadRow " << i << endl;
  //}

  calculateEta(totalNumberOfClusters);
  // Determine the minimum and maximum values of eta. That way the TPC digits can be grouped into pseudorapidity bins
  determineMinMaxEta(totalNumberOfClusters);
  // Discretize the eta values of all clusters into etaResolution bins
  calculateEtaSlice(totalNumberOfClusters);

//  draw1DCartesianClustersPerPadRow(25);
  draw1DCartesianClustersPerEtaSlice(25);
  drawCartesianClusters();

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
    hough->Transform(0, clusterCollection);
    hough->AddAllHistogramsRows();
    hough->FindTrackCandidatesRow();
    //    hough->AddTracks();
  }

  delete clusterCollection;

  return 0;
}
