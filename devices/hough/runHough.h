/// \file runHough.h
/// \brief Definition of a cluster loader
/// \author Charis Kouzinopoulos

#include "TCanvas.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TGraph2D.h"

#include "boost/filesystem.hpp"
#include "boost/format.hpp"
#include "boost/program_options.hpp"

#include "Accumulator.h"
#include "BaseTransformer.h"
#include "ClusterCollection.h"
#include "Hough.h"
#include "TransformerRow.h"
#include "Transform.h"

AliceO2::Hough::ClusterCollection* clusterCollection;

Double_t etaMin;
Double_t etaMax;

// Project parameters
int etaSlices = 100;
