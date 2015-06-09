/// \file BaseTransformer.cxx
/// \brief Implementation of the BaseTransformer class
/// \author Anders Vestbo <vestbo@fi.uib.no>

#include "BaseTransformer.h"

using namespace AliceO2::Hough;

BaseTransformer::BaseTransformer()
  : fLastTransformer(0x0),
    fSlice(0),
    fPatch(0),
    fLastPatch(-1),
    fNEtaSegments(0),
    fEtaMin(0),
    fEtaMax(0),
    fLowerThreshold(0),
    fUpperThreshold(1023),
    //fDigitRowData(0x0),
    fZVertex(0.0)

{
}

BaseTransformer::BaseTransformer(Int_t slice, Int_t patch, Int_t netasegments, Float_t zvertex)
  : fLastTransformer(0x0),
    fSlice(0),
    fPatch(0),
    fLastPatch(-1),
    fNEtaSegments(0),
    fEtaMin(0),
    fEtaMax(0),
    fLowerThreshold(3),
    fUpperThreshold(1023),
    /*fDigitRowData(0x0),*/
    fZVertex(zvertex)
{
  Init(slice, patch, netasegments);
}

BaseTransformer::~BaseTransformer()
{
}

void BaseTransformer::Init(Int_t slice, Int_t patch, Int_t netasegments, Int_t /*n_seqs*/)
{
  fSlice = slice;
  fPatch = patch;
  fLastPatch = -1;
  fNEtaSegments = netasegments;
  fEtaMin = 0;
  fEtaMax = fSlice < 18 ? 1. : -1.;
}
