// @(#) $Id: HoughBaseTransformer.cxx 15995 2006-11-30 17:45:45Z hristov $

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group
//-------------------------------------------------------------------------
//          Implementation of the HoughBaseTransformer class
//  that is the base class for AliHLTHoughTransformer,
//  AliHLTHoughTransformerVhdl, AliHLTHoughTransformerGlobal,
//  AliHLTHoughTransformerRow
//-------------------------------------------------------------------------

#include "HoughBaseTransformer.h"

/** \class HoughBaseTransformer
<pre>
//_____________________________________________________________
// HoughBaseTransformer
//
// The base class for implementations of Hough Transform on ALICE TPC data.
//
// This is an abstract class, and is only meant to provide the interface
// to the different implementations.
//
</pre>
*/

using namespace AliceO2::Hough;

HoughBaseTransformer::HoughBaseTransformer()
  : fLastTransformer(0x0),
    fSlice(0),
    fPatch(0),
    fLastPatch(-1),
    fNEtaSegments(0),
    fEtaMin(0),
    fEtaMax(0),
    fLowerThreshold(0),
    fUpperThreshold(1023),
    /*fDigitRowData(0x0),*/
    fZVertex(0.0)

{
  // Default constructor
}

HoughBaseTransformer::HoughBaseTransformer(Int_t slice, Int_t patch, Int_t netasegments, Float_t zvertex)
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
  // normal ctor

  Init(slice, patch, netasegments);
}

HoughBaseTransformer::~HoughBaseTransformer()
{
  // dtor
}

void HoughBaseTransformer::Init(Int_t slice, Int_t patch, Int_t netasegments, Int_t /*n_seqs*/)
{
  // Transformer init
  fSlice = slice;
  fPatch = patch;
  fLastPatch = -1;
  fNEtaSegments = netasegments;
  fEtaMin = 0;
  fEtaMax = fSlice < 18 ? 1. : -1.;
}
