/// \file runHough.cxx
/// \brief Implementation of a cluster loader
/// \author Charis Kouzinopoulos

/// List of references: [1] Cheshkov, C. "Fast Hough-transform track reconstruction for the ALICE TPC." Nuclear
/// Instruments and Methods in Physics Research Section A: Accelerators, Spectrometers, Detectors and Associated
/// Equipment 566.1 (2006): 35-39.

#include "runHough.h"

using namespace AliceO2::Hough;

typedef struct DeviceOptions {
  std::string event;

  int TPCSlice;
  int etaSlice;

  bool print;
  bool printAll;
} DeviceOptions_t;

void draw1DCartesianClustersPerTPCSliceAndEtaSlice(int TPCSlice, int etaSlice)
{
  TCanvas* cartesianClustersPerTPCSliceAndEtaSliceCanvas1D =
    new TCanvas("cartesianClustersPerTPCSliceAndEtaSliceCanvas1D", "Cartesian clusters per TPC and eta slice");
  TGraph* cartesianClustersPerTPCSliceAndEtaSliceGraph1D = new TGraph();

  for (UInt_t padRow = 0; padRow < Transform::GetNRows(); padRow++) {
    for (UInt_t clusterNumber = 0; clusterNumber < clusterCollection->getNumberOfClustersPerPadRow(padRow);
         clusterNumber++) {
      if (clusterCollection->getClusterTPCSlice(padRow, clusterNumber) == TPCSlice &&
          clusterCollection->getClusterEtaSlice(padRow, clusterNumber) == etaSlice) {
        cartesianClustersPerTPCSliceAndEtaSliceGraph1D->SetPoint(clusterNumber,
                                                                 clusterCollection->getClusterX(padRow, clusterNumber),
                                                                 clusterCollection->getClusterY(padRow, clusterNumber));
      }
    }
  }

  cartesianClustersPerTPCSliceAndEtaSliceGraph1D->SetMarkerStyle(7);
  cartesianClustersPerTPCSliceAndEtaSliceGraph1D->Draw("AP");
  cartesianClustersPerTPCSliceAndEtaSliceGraph1D->GetXaxis()->SetRangeUser(80, 240);
  cartesianClustersPerTPCSliceAndEtaSliceGraph1D->GetXaxis()->SetTitle("Pad Row");

  cartesianClustersPerTPCSliceAndEtaSliceCanvas1D->Print("cartesianClustersPerEtaSlicePerTPCSlice.pdf");
}

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

void printClusterInformation()
{
  cout << "ID" << setw(13) << "X" << setw(10) << "Y" << setw(10) << "Z" << setw(9) << "Pad" << setw(9) << "Time"
       << setw(10) << "Charge" << setw(10) << "Sector" << setw(10) << "Eta" << endl;

  for (UInt_t padRow = 0; padRow < Transform::GetNRows(); padRow++) {
    for (UInt_t clusterNumber = 0; clusterNumber < clusterCollection->getNumberOfClustersPerPadRow(padRow);
         clusterNumber++) {
      cout << (AliHLTUInt32_t)clusterCollection->getClusterID(padRow, clusterNumber) << setw(9)
           << clusterCollection->getClusterX(padRow, clusterNumber) << setw(10)
           << clusterCollection->getClusterY(padRow, clusterNumber) << setw(10)
           << clusterCollection->getClusterZ(padRow, clusterNumber) << setw(8)
           << clusterCollection->getClusterPad(padRow, clusterNumber) << setw(7)
           << clusterCollection->getClusterTime(padRow, clusterNumber) << setw(5)
           << clusterCollection->getClusterCharge(padRow, clusterNumber) << setw(5)
           << clusterCollection->getClusterTPCSlice(padRow, clusterNumber) << setw(5)
           << clusterCollection->getClusterEtaSlice(padRow, clusterNumber) << endl;
    }
  }
}

void printClusterInformation(int TPCSlice, int etaSlice)
{
  cout << "ID" << setw(13) << "X" << setw(10) << "Y" << setw(10) << "Z" << setw(9) << "Pad" << setw(9) << "Time"
       << setw(10) << "Charge" << endl;

  for (UInt_t padRow = 0; padRow < Transform::GetNRows(); padRow++) {
    for (UInt_t clusterNumber = 0; clusterNumber < clusterCollection->getNumberOfClustersPerPadRow(padRow);
         clusterNumber++) {

      if (clusterCollection->getClusterTPCSlice(padRow, clusterNumber) != TPCSlice ||
          clusterCollection->getClusterEtaSlice(padRow, clusterNumber) != etaSlice) {
        continue;
      }

      cout << (AliHLTUInt32_t)clusterCollection->getClusterID(padRow, clusterNumber) << setw(9)
           << clusterCollection->getClusterX(padRow, clusterNumber) << setw(10)
           << clusterCollection->getClusterY(padRow, clusterNumber) << setw(10)
           << clusterCollection->getClusterZ(padRow, clusterNumber) << setw(8)
           << clusterCollection->getClusterPad(padRow, clusterNumber) << setw(9)
           << clusterCollection->getClusterTime(padRow, clusterNumber) << setw(5)
           << clusterCollection->getClusterCharge(padRow, clusterNumber) << endl;
    }
  }
}

inline bool parseCommandLine(int _argc, char* _argv[], DeviceOptions* _options)
{
  if (_options == NULL) {
    throw runtime_error("Internal error: options' container is empty.");
  }

  namespace bpo = boost::program_options;
  bpo::options_description desc("Options");
  desc.add_options()("event", bpo::value<std::string>()->required(), "Specify an event to load from disk <mandatory>")(
    "TPCSlice", bpo::value<int>()->default_value(25), "Specify a TPC slice to draw clusters from")(
    "etaSlice", bpo::value<int>()->default_value(29), "Specify an eta slice to draw clusters from")(
    "print", "Print information on clusters for the given TPC and eta slices and exit")(
    "print-all", "Print information on all clusters and exit")("help", "Print this help message");

  bpo::variables_map vm;
  bpo::store(bpo::parse_command_line(_argc, _argv, desc), vm);

  if (vm.count("help")) {
    cout << desc << endl;
    return false;
  }

  bpo::notify(vm);

  if (vm.count("event")) {
    _options->event = vm["event"].as<std::string>();
  }
  if (vm.count("TPCSlice")) {
    _options->TPCSlice = vm["TPCSlice"].as<int>();
  }
  if (vm.count("etaSlice")) {
    _options->etaSlice = vm["etaSlice"].as<int>();
  }
  if (vm.count("print")) {
    _options->print = true;
  } else {
    _options->print = false;
  }
  if (vm.count("print-all")) {
    _options->printAll = true;
  } else {
    _options->printAll = false;
  }

  return true;
}

int main(int argc, char** argv)
{

  DeviceOptions_t options;
  try {
    if (!parseCommandLine(argc, argv, &options)) {
      return 0;
    }
  }
  catch (const exception& e) {
    cerr << e.what() << endl;
    return 1;
  }

  // Create data path
  std::string dataFilename = "emulated-tpc-clusters/event";
  dataFilename += options.event;

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

  calculateEta(totalNumberOfClusters);
  // Determine the minimum and maximum values of eta. That way the TPC digits can be grouped into pseudorapidity bins
  determineMinMaxEta(totalNumberOfClusters);
  // Discretize the eta values of all clusters into etaResolution bins
  calculateEtaSlice(totalNumberOfClusters);

  //  draw1DCartesianClustersPerPadRow(25);
  draw1DCartesianClustersPerEtaSlice(options.etaSlice);
  draw1DCartesianClustersPerTPCSliceAndEtaSlice(options.TPCSlice, options.etaSlice);
  drawCartesianClusters();

  if (options.print || options.printAll) {
    for (int i = 0; i < Transform::GetNRows(); i++) {
      cout << clusterCollection->getNumberOfClustersPerPadRow(i) << " Clusters for PadRow " << i << endl;
    }
  }

  if (options.print) {
    printClusterInformation(options.TPCSlice, options.etaSlice);
    return 1;
  }

  if (options.printAll) {
    printClusterInformation();
    return 1;
  }

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
