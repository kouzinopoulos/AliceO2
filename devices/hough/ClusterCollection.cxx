/// \file ClusterCollection.cxx
/// \brief Implementation of the ClusterCollection class
/// \author Charis Kouzinopoulos <charalampos.kouzinopoulos@cern.ch>

#include "Transform.h"
#include "ClusterCollection.h"

using namespace std;
using namespace AliceO2::Hough;

// One clusterDataFormat vector is allocated per pad row. Total vectors are Transform::GetNRows() (159)
ClusterCollection::ClusterCollection() : clusterData(Transform::GetNRows()) {}

ClusterCollection::~ClusterCollection() {}

UInt_t ClusterCollection::readData(std::string dataPath, std::string dataType, std::string dataOrigin, int TPCSlice)
{
  // Retrieve the TPC slice and partition from the filename
  std::string currentSliceString(dataPath, dataPath.length() - 6, 2);
  std::string currentPartitionString(dataPath, dataPath.length() - 2, 2);

  AliHLTUInt8_t currentSlice = std::stoul(currentSliceString, nullptr, 16);
  AliHLTUInt8_t currentPartition = std::stoul(currentPartitionString, nullptr, 16);

  //If the file does not correspond to the given TPC slice provided by the user, return
  if (currentSlice != TPCSlice) {
    return 0;
  }

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

  // Initialize a cluster point collection
  spacepoints = std::unique_ptr<AliHLTTPCSpacePointContainer>(new AliHLTTPCSpacePointContainer);
  if (!spacepoints.get()) {
    std::cerr << "Error, could not create a space point collection" << endl;
    std::exit(1);
  }

  // Create an AliHLTComponentBlockData object, fill it with default values and then set its pointer to
  // the data buffer
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

  // cout << "ID" << setw(8) << "fX" << setw(8) << "fY" << setw(8) << "fZ" << setw(9) << "Pad Row" << setw(10)
  //     << "SigmaY2" << setw(10) << "SigmaZ2" << setw(7) << "Charge" << setw(7) << "QMax" << setw(10) << "Track Id"
  //     << setw(7) << "MC Id" << setw(7) << "Used" << endl;

  // cout << *spacepoints << endl;

  // Retrieve the cluster information from AliHLTTPCSpacePointContainer
  std::vector<AliHLTUInt32_t> clusterIDs;
  spacepoints->GetClusterIDs(clusterIDs);

  // Append the cluster IDs and their X, Y and Z coordinates to the clusterCartesianCoordinates vector
  for (vector<AliHLTUInt32_t>::const_iterator element = clusterIDs.begin(); element != clusterIDs.end(); element++) {
    AliHLTUInt32_t clusterID = *element;

    // Convert local coordinates to raw coordinates and store them
    Float_t coordinates[3] = { spacepoints->GetX(clusterID), spacepoints->GetY(clusterID),
                               spacepoints->GetZ(clusterID) };
    Transform::LocHLT2Raw(coordinates, currentSlice, spacepoints->GetPadRow(clusterID));

    setClusterParameters(clusterID, spacepoints->GetX(clusterID), spacepoints->GetY(clusterID),
                          spacepoints->GetZ(clusterID), spacepoints->GetCharge(clusterID),
                          spacepoints->GetPadRow(clusterID), coordinates[1], coordinates[2], currentSlice,
                          currentPartition);
  }

  // De-allocate memory space
  if (inputBuffer) {
    delete[] inputBuffer;
  }
  inputBuffer = NULL;

  return numberOfClusters;
}
