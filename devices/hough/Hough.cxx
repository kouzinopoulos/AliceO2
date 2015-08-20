/// \file Hough.cxx
/// \brief Implementation of the Hough class
/// \author Anders Vestbo, Cvetan Cheshkov

#include <sys/time.h>

//#include "AliHLTTPCLogging.h"

#ifdef HAVE_ALIHLTHOUGHMERGER
#include "AliHLTHoughMerger.h"
#endif // HAVE_ALIHLTHOUGHMERGER

#ifdef HAVE_ALIHLTHOUGHINTMERGER
#include "AliHLTHoughIntMerger.h"
#endif // HAVE_ALIHLTHOUGHINTMERGER

#ifdef HAVE_ALIHLTHOUGHGLOBALMERGER
#include "AliHLTHoughGlobalMerger.h"
#endif // HAVE_ALIHLTHOUGHGLOBALMERGER

#include "Accumulator.h"
#include "BaseTransformer.h"
#include "Hough.h"
#include "HoughTrack.h"
#include "MaxFinder.h"
#include "StandardIncludes.h"
#include "Transform.h"
#include "TransformerRow.h"

//#include "AliHLTDigitData.h"
//#include "HoughEval.h"
//#include "AliHLTTPCTransform.h"
//#include "AliHLTTPCTrackArray.h"

#ifdef HAVE_ALIHLTDDLDATAFILEHANDLER
#include "AliHLTDDLDataFileHandler.h"
#endif // HAVE_ALIHLTDDLDATAFILEHANDLER

//#include "HoughKalmanTrack.h"

#ifdef HAVE_THREAD
#include "TThread.h"
#endif // HAVE_THREAD

#if __GNUC__ >= 3
using namespace std;
#endif

using namespace AliceO2::Hough;

Hough::Hough()
{
  // Constructor

  fBinary = false;
  fAddHistograms = false;
  fDoIterative = false;
  fWriteDigits = false;
  fUse8bits = false;

  //  fMemHandler = 0;
  fHoughTransformer = 0;
  //  fEval = 0;
  fPeakFinder = 0;
  //  fTracks = 0;
  //  fGlobalTracks = 0;
  //  fMerger = 0;
  //  fInterMerger = 0;
  //  fGlobalMerger = 0;
  //  fBenchmark = 0;

  fNEtaSegments = 0;
  fNPatches = 0;
  fLastPatch = -1;
  //  fVersion = 0;
  fCurrentSlice = 0;
  fEvent = 0;

  fKappaSpread = 6;
  fPeakRatio = 0.5;
  fInputFile = 0;
  fInputPtr = 0;
  //  fRawEvent = 0;

  SetTransformerParams();
  SetThreshold();
  SetNSaveIterations();
  SetPeakThreshold();
  // just be sure that index is empty for new event
  /*  AliHLTTPCFileHandler::CleanStaticIndex();
    fRunLoader = 0;
  #ifdef HAVE_THREAD
    fThread = 0;
  #endif // HAVE_THREAD
  */
}

Hough::Hough(Char_t* path, Bool_t binary, Int_t netasegments, Bool_t bit8, Int_t tv, Char_t* infile, Char_t* ptr)
{
  // Normal constructor
  fBinary = binary;
  strcpy(fPath, path);
  fNEtaSegments = netasegments;
  fAddHistograms = false;
  fDoIterative = false;
  fWriteDigits = false;
  fUse8bits = bit8;
  //  fVersion = tv;
  fKappaSpread = 6;
  fPeakRatio = 0.5;
  if (!fBinary) {
    if (infile) {
      fInputFile = infile;
      fInputPtr = 0;
    } else {
      fInputFile = 0;
      fInputPtr = ptr;
    }
  } else {
    fInputFile = 0;
    fInputPtr = 0;
  }
// just be sure that index is empty for new event
//  AliHLTTPCFileHandler::CleanStaticIndex();
//  fRunLoader = 0;
#ifdef HAVE_THREAD
  fThread = 0;
#endif // HAVE_THREAD
}

Hough::~Hough()
{
  // dtor

  CleanUp();
  /*
  #ifdef HAVE_ALIHLTHOUGHMERGER
    if (fMerger)
      delete fMerger;
  #endif // HAVE_ALIHLTHOUGHMERGER
  // cout << "Cleaned class merger " << endl;
  #ifdef HAVE_ALIHLTHOUGHINTMERGER
    if (fInterMerger)
      delete fInterMerger;
  #endif // HAVE_ALIHLTHOUGHINTMERGER
    // cout << "Cleaned class inter " << endl;
    if (fPeakFinder)
      delete fPeakFinder;
  // cout << "Cleaned class peak " << endl;
  #ifdef HAVE_ALIHLTHOUGHGLOBALMERGER
    if (fGlobalMerger)
      delete fGlobalMerger;
  #endif // HAVE_ALIHLTHOUGHGLOBALMERGER
    // cout << "Cleaned class global " << endl;
    if (fBenchmark)
      delete fBenchmark;
    // cout << "Cleaned class bench " << endl;
    if (fGlobalTracks)
      delete fGlobalTracks;
  // cout << "Cleaned class globaltracks " << endl;
  #ifdef HAVE_THREAD
    if (fThread) {
      //    fThread->Delete();
      delete fThread;
      fThread = 0;
    }
  #endif // HAVE_THREAD*/
}

void Hough::CleanUp()
{
  // Cleanup memory

  for (Int_t i = 0; i < fNPatches; i++) {
    /*
        if (fTracks[i])
          delete fTracks[i];
        // cout << "Cleaned tracks " << i << endl;
        if (fEval[i])
          delete fEval[i];
        // cout << "Cleaned eval " << i << endl;
    */

    if (fHoughTransformer[i])
      delete fHoughTransformer[i];
    // cout << "Cleaned traf " << i << endl;
    /*
        if (fMemHandler[i])
          delete fMemHandler[i];
        // cout << "Cleaned mem " << i << endl;
    */
  }
  /*
    if (fTracks)
      delete[] fTracks;
    // cout << "Cleaned class tracks " << endl;
    if (fEval)
      delete[] fEval;
    // cout << "Cleaned class eval " << endl;
    */
  if (fHoughTransformer)
    delete[] fHoughTransformer;
  // cout << "Cleaned cleass trafo " << endl;
  /*  if (fMemHandler)
      delete[] fMemHandler;
      */
  // cout << "Cleaned class mem " << endl;
}
/*
void Hough::Init(Int_t netasegments, Int_t tv, AliRawEvent* rawevent, Float_t zvertex)
{
  // Normal constructor
  fNEtaSegments = netasegments;
  fVersion = tv;
  fRawEvent = rawevent;
  fZVertex = zvertex;

  Init();
}
*/
void Hough::Init(Char_t* path, Bool_t binary, Int_t netasegments, Bool_t bit8, Int_t tv, Char_t* infile, Char_t* ptr,
                 Float_t zvertex)
{
  // Normal init of the Hough
  fBinary = binary;
  strcpy(fPath, path);
  fNEtaSegments = netasegments;
  fWriteDigits = false;
  fUse8bits = bit8;
  //  fVersion = tv;
  if (!fBinary) {
    if (infile) {
      fInputFile = infile;
      fInputPtr = 0;
    } else {
      fInputFile = 0;
      fInputPtr = ptr;
    }
  } else {
    fInputFile = 0;
    fInputPtr = 0;
  }
  fZVertex = zvertex;

  Init(); // do the rest
}

void Hough::Init(Bool_t doit, Bool_t addhists)
{
  // Init
  fDoIterative = doit;
  fAddHistograms = addhists;

  fNPatches = Transform::GetNPatches();
  fHoughTransformer = new TransformerRow* [fNPatches];
  /*  fMemHandler = new AliHLTTPCMemHandler* [fNPatches];


    fTracks = new AliHLTTPCTrackArray* [fNPatches];
    fEval = new HoughEval* [fNPatches];

    fGlobalTracks = new AliHLTTPCTrackArray("HoughTrack");
  */
  BaseTransformer* lasttransformer = 0;

  for (Int_t i = 0; i < fNPatches; i++) {
    fHoughTransformer[i] = new TransformerRow(0, i, fNEtaSegments, false, fZVertex);
    fHoughTransformer[i]->SetLastTransformer(lasttransformer);
    lasttransformer = fHoughTransformer[i];
    //      fHoughTransformer[i]->CreateHistograms(fNBinX[i],fLowPt[i],fNBinY[i],-fPhi[i],fPhi[i]);
    fHoughTransformer[i]->CreateHistograms(fNBinX[i], -fLowPt[i], fLowPt[i], fNBinY[i], -fPhi[i], fPhi[i]);
    // fHoughTransformer[i]->CreateHistograms(fLowPt[i],fUpperPt[i],fPtRes[i],fNBinY[i],fPhi[i]);

    fHoughTransformer[i]->SetLowerThreshold(fThreshold[i]);
    fHoughTransformer[i]->SetUpperThreshold(100);

    cout << "Initializing Hough transformer" << endl;
  }

  fPeakFinder = new MaxFinder("KappaPhi", 50000);
  //  fMerger = 0;
  //  fInterMerger = 0;

  //  fGlobalMerger = 0;
  //    fBenchmark = new AliHLTTPCBenchmark();
}

void Hough::SetTransformerParams(Float_t ptres, Float_t ptmin, Float_t ptmax, Int_t ny, Int_t patch)
{
  // Setup the parameters for the Hough Transformer
  // This includes the bin size and limits for
  // the parameter space histograms

  Int_t mrow;
  Float_t psi = 0;
  if (patch == -1)
    mrow = 80;
  else
    mrow = Transform::GetLastRow(patch);
  if (ptmin) {
    Double_t lineradius = sqrt(pow(Transform::Row2X(mrow), 2) + pow(Transform::GetMaxY(mrow), 2));
    Double_t kappa = -1 * Transform::GetBField() * Transform::GetBFact() / ptmin;
    psi = Transform::Deg2Rad(10) - asin(lineradius * kappa / 2);
    cout << "Calculated psi range " << psi << " in patch " << patch << endl;
  }

  if (patch == -1) {
    Int_t i = 0;
    while (i < 6) {
      fPtRes[i] = ptres;
      fLowPt[i] = ptmin;
      fUpperPt[i] = ptmax;
      fNBinY[i] = ny;
      fPhi[i] = psi;
      fNBinX[i] = 0;
      i++;
    }
    return;
  }

  fPtRes[patch] = ptres;
  fLowPt[patch] = ptmin;
  fUpperPt[patch] = ptmax;
  fNBinY[patch] = ny;
  fPhi[patch] = psi;
}
/*
void Hough::SetTransformerParams(Int_t nx,Int_t ny,Float_t ptmin,Int_t patch)
{
  // Setup the parameters for the Hough Transformer

  Int_t mrow=80;
  Double_t lineradius = sqrt(pow(Transform::Row2X(mrow),2) + pow(Transform::GetMaxY(mrow),2));
  Double_t kappa = -1*Transform::GetBField()*Transform::GetBFact()/ptmin;
  Double_t psi = Transform::Deg2Rad(10) - asin(lineradius*kappa/2);
  cout<<"Calculated psi range "<<psi<<" in patch "<<patch<<endl;

  Int_t i=0;
  while(i < 6)
    {
      fLowPt[i] = ptmin;
      fNBinY[i] = ny;
      fNBinX[i] = nx;
      fPhi[i] = psi;
      i++;
    }
}
*/
void Hough::SetTransformerParams(Int_t nx, Int_t ny, Float_t ptmin, Int_t /*patch*/)
{
  // Setup the parameters for the Hough Transformer

  Double_t lineradius =
    1.0 / (TransformerRow::GetBeta1() * sqrt(1.0 + tan(Transform::Pi() * 10 / 180) * tan(Transform::Pi() * 10 / 180)));
  Double_t alpha1 = TransformerRow::GetBeta1() * tan(Transform::Pi() * 10 / 180);
  Double_t kappa = 1 * Transform::GetBField() * Transform::GetBFact() / (ptmin * 0.9);
  Double_t psi = Transform::Deg2Rad(10) - asin(lineradius * kappa / 2);
  //  cout<<"Calculated psi range "<<psi<<" in patch "<<patch<<endl;
  Double_t alpha2 = alpha1 - (TransformerRow::GetBeta1() - TransformerRow::GetBeta2()) * tan(psi);
  //  cout<<"Calculated alphas range "<<alpha1<<" "<<alpha2<<" in patch "<<patch<<endl;

  Int_t i = 0;
  while (i < 6) {
    fLowPt[i] = 1.1 * alpha1;
    fNBinY[i] = ny;
    fNBinX[i] = nx;
    fPhi[i] = alpha2;
    i++;
  }
}

void Hough::CalcTransformerParams(Float_t ptmin)
{
  // Setup the parameters for the Row Hough Transformer
  // Automatically adjusts the number of bins in X and Y in a way
  // that the size of the hough bin is 2x (in X) and 2.5 (in Y) the
  // size of the tpc pads

  Double_t lineradius =
    1.0 / (TransformerRow::GetBeta1() * sqrt(1.0 + tan(Transform::Pi() * 10 / 180) * tan(Transform::Pi() * 10 / 180)));
  Double_t alpha1 = TransformerRow::GetBeta1() * tan(Transform::Pi() * 10 / 180);
  Double_t kappa = 1 * Transform::GetBField() * Transform::GetBFact() / (ptmin * 0.9);
  Double_t psi = Transform::Deg2Rad(10) - asin(lineradius * kappa / 2);
  //  cout<<"Calculated psi range "<<psi<<endl;
  Double_t alpha2 = alpha1 - (TransformerRow::GetBeta1() - TransformerRow::GetBeta2()) * tan(psi);
  alpha1 *= 1.1;
  //  cout<<"Calculated alphas range "<<alpha1<<" "<<alpha2<<endl;

  Double_t sizex = 2.0 * Transform::GetPadPitchWidthLow() * TransformerRow::GetBeta1() * TransformerRow::GetBeta1();
  Double_t sizey = 2.5 * Transform::GetPadPitchWidthUp() * TransformerRow::GetBeta2() * TransformerRow::GetBeta2();

  Int_t nx = 2 * (Int_t)(alpha1 / sizex) + 1;
  Int_t ny = 2 * (Int_t)(alpha2 / sizey) + 1;
  //  cout<<"Calculated number of bins "<<nx<<" "<<ny<<endl;

  Int_t i = 0;
  while (i < 6) {
    fLowPt[i] = alpha1;
    fNBinY[i] = ny;
    fNBinX[i] = nx;
    fPhi[i] = alpha2;
    i++;
  }
}

void Hough::SetTransformerParams(Int_t nx, Int_t ny, Float_t lpt, Float_t phi)
{
  // SetTransformerParams

  Int_t i = 0;
  while (i < 6) {
    fLowPt[i] = lpt;
    fNBinY[i] = ny;
    fNBinX[i] = nx;
    fPhi[i] = phi;
    i++;
  }
}

void Hough::SetThreshold(Int_t t3, Int_t patch)
{
  // Set digits threshold
  if (patch == -1) {
    Int_t i = 0;
    while (i < 6)
      fThreshold[i++] = t3;
    return;
  }
  fThreshold[patch] = t3;
}

void Hough::SetPeakThreshold(Int_t threshold, Int_t patch)
{
  // Set Peak Finder threshold
  if (patch == -1) {
    Int_t i = 0;
    while (i < 6)
      fPeakThreshold[i++] = threshold;
    return;
  }
  fPeakThreshold[patch] = threshold;
}

// void Hough::DoBench(Char_t* name) { fBenchmark->Analyze(name); }

void Hough::Transform(Int_t* rowrange, ClusterCollection* clusterCollection)
{
  // Transform all data given to the transformer within the given slice
  //(after ReadData(slice))

  Int_t patchorder[6] = { 5, 2, 0, 1, 3, 4 }; // The order in which patches are processed
  //  Int_t patchorder[6] = {0,1,2,3,4,5}; //The order in which patches are processed
  //  Int_t patchorder[6] = {5,4,3,2,1,0}; //The order in which patches are processed
  //  Int_t patchorder[6] = {5,2,4,3,1,0}; //The order in which patches are processed
  fLastPatch = -1;
  for (Int_t i = 0; i < fNPatches; i++) {
    // In case of Row transformer reset the arrays only once
    PrepareForNextPatch(patchorder[i]);
    if (!rowrange) {
      fHoughTransformer[patchorder[i]]->SetLastPatch(fLastPatch);
      fHoughTransformer[patchorder[i]]->TransformCircle(clusterCollection);
    } else {
      fHoughTransformer[i]->TransformCircleC(rowrange, 1);
    }
    fLastPatch = patchorder[i];
  }

  // FIXME: Add a real timer here
  cout << "Transform done in average per patch of " << 0 / fNPatches << " ms" << endl;
}

void Hough::AddAllHistogramsRows()
{
  // Add the histograms within one etaslice.
  // Resulting histogram are in patch=0.

  UChar_t lastpatchlastrow = Transform::GetLastRowOnDDL(fLastPatch) + 1;

  UChar_t* tracklastrow = ((TransformerRow*)fHoughTransformer[0])->GetTrackLastRow();

  for (Int_t i = 0; i < fNEtaSegments; i++) {
    UChar_t* gapcount = ((TransformerRow*)fHoughTransformer[0])->GetGapCount(i);
    UChar_t* currentrowcount = ((TransformerRow*)fHoughTransformer[0])->GetCurrentRowCount(i);

    Accumulator* hist = fHoughTransformer[0]->GetHistogram(i);
    Int_t xmin = hist->GetFirstXbin();
    Int_t xmax = hist->GetLastXbin();
    Int_t ymin = hist->GetFirstYbin();
    Int_t ymax = hist->GetLastYbin();
    Int_t nxbins = hist->GetNbinsX() + 2;

    for (Int_t ybin = ymin; ybin <= ymax; ybin++) {
      for (Int_t xbin = xmin; xbin <= xmax; xbin++) {
        Int_t bin = xbin + ybin * nxbins; // Int_t bin = hist->GetBin(xbin,ybin);
        if (gapcount[bin] < MAX_N_GAPS) {
          if (tracklastrow[bin] > lastpatchlastrow) {
            if (lastpatchlastrow > currentrowcount[bin])
              gapcount[bin] += (lastpatchlastrow - currentrowcount[bin] - 1);
          } else {
            if (tracklastrow[bin] > currentrowcount[bin])
              gapcount[bin] += (tracklastrow[bin] - currentrowcount[bin] - 1);
          }
          if (gapcount[bin] < MAX_N_GAPS)
            hist->AddBinContent(bin, (159 - gapcount[bin]));
        }
      }
    }
  }

  fAddHistograms = kTRUE;
  cout << "Adding histograms in " << 0 << " ms" << endl;
}

void Hough::PrepareForNextPatch(Int_t nextpatch)
{
  // Prepare the parameter space for the processing of
  // the next read patch. According to the already
  // accumulated number of gaps in parameter space
  // bins, the routine updates the dynamic
  // pointers used in order to jump rapidly during the
  // filling of the parameter space.

  UChar_t lastpatchlastrow;
  if (fLastPatch == -1) {
    lastpatchlastrow = 0;
  } else {
    lastpatchlastrow = Transform::GetLastRowOnDDL(fLastPatch) + 1;
  }

  UChar_t nextpatchfirstrow;

  if (nextpatch == 0) {
    nextpatchfirstrow = 0;
  } else {
    nextpatchfirstrow = Transform::GetFirstRowOnDDL(nextpatch) - 1;
  }

  UChar_t* trackfirstrow = ((TransformerRow*)fHoughTransformer[0])->GetTrackFirstRow();
  UChar_t* tracklastrow = ((TransformerRow*)fHoughTransformer[0])->GetTrackLastRow();

  for (Int_t i = 0; i < fNEtaSegments; i++) {
    UChar_t* gapcount = ((TransformerRow*)fHoughTransformer[0])->GetGapCount(i);
    UChar_t* currentrowcount = ((TransformerRow*)fHoughTransformer[0])->GetCurrentRowCount(i);
    UChar_t* prevbin = ((TransformerRow*)fHoughTransformer[0])->GetPrevBin(i);
    UChar_t* nextbin = ((TransformerRow*)fHoughTransformer[0])->GetNextBin(i);
    UChar_t* nextrow = ((TransformerRow*)fHoughTransformer[0])->GetNextRow(i);

    Accumulator* hist = fHoughTransformer[0]->GetHistogram(i);
    Int_t xmin = hist->GetFirstXbin();
    Int_t xmax = hist->GetLastXbin();
    Int_t ymin = hist->GetFirstYbin();
    Int_t ymax = hist->GetLastYbin();
    Int_t nxbins = hist->GetNbinsX() + 2;

    if (fLastPatch != -1) {
      UChar_t lastyvalue = 0;
      Int_t endybin = ymin - 1;
      for (Int_t ybin = nextrow[ymin]; ybin <= ymax; ybin = nextrow[++ybin]) {
        UChar_t lastxvalue = 0;
        UChar_t maxvalue = 0;
        Int_t endxbin = xmin - 1;
        for (Int_t xbin = xmin; xbin <= xmax; xbin++) {
          Int_t bin = xbin + ybin * nxbins;
          UChar_t value = 0;
          if (gapcount[bin] < MAX_N_GAPS) {
            if (tracklastrow[bin] > lastpatchlastrow) {
              if (lastpatchlastrow > currentrowcount[bin])
                gapcount[bin] += (lastpatchlastrow - currentrowcount[bin] - 1);
            } else {
              if (tracklastrow[bin] > currentrowcount[bin])
                gapcount[bin] += (tracklastrow[bin] - currentrowcount[bin] - 1);
            }
            if (gapcount[bin] < MAX_N_GAPS) {
              value = 1;
              maxvalue = 1;
              if (trackfirstrow[bin] < nextpatchfirstrow) {
                currentrowcount[bin] = nextpatchfirstrow;
              } else {
                currentrowcount[bin] = trackfirstrow[bin];
              }
            }
          }
          if (value > 0) {
            nextbin[xbin + ybin * nxbins] = (UChar_t)xbin;
            prevbin[xbin + ybin * nxbins] = (UChar_t)xbin;
            if (value > lastxvalue) {
              UChar_t* tempnextbin = nextbin + endxbin + 1 + ybin * nxbins;
              memset(tempnextbin, (UChar_t)xbin, xbin - endxbin - 1);
            }
            endxbin = xbin;
          } else {
            prevbin[xbin + ybin * nxbins] = (UChar_t)endxbin;
          }
          lastxvalue = value;
        }
        UChar_t* tempnextbin = nextbin + endxbin + 1 + ybin * nxbins;
        memset(tempnextbin, (UChar_t)(xmax + 1), xmax - endxbin);
        if (maxvalue > 0) {
          nextrow[ybin] = (UChar_t)ybin;
          if (maxvalue > lastyvalue) {
            UChar_t* tempnextrow = nextrow + endybin + 1;
            memset(tempnextrow, (UChar_t)ybin, ybin - endybin - 1);
          }
          endybin = ybin;
        }
        lastyvalue = maxvalue;
      }
      UChar_t* tempnextrow = nextrow + endybin + 1;
      memset(tempnextrow, (UChar_t)(ymax + 1), ymax - endybin + 1);
    } else {
      UChar_t lastyvalue = 0;
      Int_t endybin = ymin - 1;
      for (Int_t ybin = ymin; ybin <= ymax; ybin++) {
        UChar_t maxvalue = 0;
        for (Int_t xbin = xmin; xbin <= xmax; xbin++) {
          Int_t bin = xbin + ybin * nxbins;
          if (gapcount[bin] < MAX_N_GAPS) {
            maxvalue = 1;
            if (trackfirstrow[bin] < nextpatchfirstrow) {
              currentrowcount[bin] = nextpatchfirstrow;
            } else {
              currentrowcount[bin] = trackfirstrow[bin];
            }
          }
        }
        if (maxvalue > 0) {
          nextrow[ybin] = (UChar_t)ybin;
          if (maxvalue > lastyvalue) {
            UChar_t* tempnextrow = nextrow + endybin + 1;
            memset(tempnextrow, (UChar_t)ybin, ybin - endybin - 1);
          }
          endybin = ybin;
        }
        lastyvalue = maxvalue;
      }
      UChar_t* tempnextrow = nextrow + endybin + 1;
      memset(tempnextrow, (UChar_t)(ymax + 1), ymax - endybin + 1);
    }
  }
}

void Hough::FindTrackCandidatesRow()
{
  // Find TransformerRow track candidates
  // Look for peaks in histograms, and find the track candidates
  Int_t npatches;
  if (fAddHistograms)
    npatches = 1; // Histograms have been added.
  else
    npatches = fNPatches;

  for (Int_t i = 0; i < npatches; i++) {
    TransformerRow* tr = fHoughTransformer[i];
    Accumulator* h = tr->GetHistogram(0);
    Float_t deltax = h->GetBinWidthX() * TransformerRow::GetDAlpha();
    Float_t deltay = h->GetBinWidthY() * TransformerRow::GetDAlpha();
    Float_t deltaeta = (tr->GetEtaMax() - tr->GetEtaMin()) / tr->GetNEtaSegments() * TransformerRow::GetDEta();
    Float_t zvertex = tr->GetZVertex();
    //    fTracks[i]->Reset();
    fPeakFinder->Reset();

    for (Int_t j = 0; j < fNEtaSegments; j++) {
      Accumulator* hist = tr->GetHistogram(j);
      if (hist->GetNEntries() == 0)
        continue;
      fPeakFinder->SetHistogram(hist);
      fPeakFinder->SetEtaSlice(j);
      fPeakFinder->SetTrackLUTs(((TransformerRow*)tr)->GetTrackNRows(), ((TransformerRow*)tr)->GetTrackFirstRow(),
                                ((TransformerRow*)tr)->GetTrackLastRow(), ((TransformerRow*)tr)->GetNextRow(j));
#ifdef do_mc
      cout << "Starting " << j << " etaslice" << endl;
#endif
      fPeakFinder->SetThreshold(fPeakThreshold[i]);
      fPeakFinder->FindAdaptedRowPeaks(1, 0, 0); // Maxima finder for HoughTransformerRow
    }

    for (Int_t k = 0; k < fPeakFinder->GetEntries(); k++) {
      //    if(fPeakFinder->GetWeight(k) < 0) continue;
      HoughTrack* track = (HoughTrack*)fTracks[i]->NextTrack();
      Double_t starteta = tr->GetEta(fPeakFinder->GetStartEta(k), fCurrentSlice);
      Double_t endeta = tr->GetEta(fPeakFinder->GetEndEta(k), fCurrentSlice);
      Double_t eta = (starteta + endeta) / 2.0;
      track->SetTrackParametersRow(fPeakFinder->GetXPeak(k), fPeakFinder->GetYPeak(k), eta, fPeakFinder->GetWeight(k));
      track->SetPterr(deltax);
      track->SetPsierr(deltay);
      track->SetTglerr(deltaeta);
      track->SetBinXY(fPeakFinder->GetXPeak(k), fPeakFinder->GetYPeak(k), fPeakFinder->GetXPeakSize(k),
                      fPeakFinder->GetYPeakSize(k));
      track->SetZ0(zvertex);
      Int_t etaindex = (fPeakFinder->GetStartEta(k) + fPeakFinder->GetEndEta(k)) / 2;
      track->SetEtaIndex(etaindex);
      Int_t rows[2];
      ((TransformerRow*)tr)->GetTrackLength(fPeakFinder->GetXPeak(k), fPeakFinder->GetYPeak(k), rows);
      track->SetRowRange(rows[0], rows[1]);
      track->SetSector(fCurrentSlice);
      track->SetSlice(fCurrentSlice);
#ifdef do_mc
      Int_t label = tr->GetTrackID(etaindex, fPeakFinder->GetXPeak(k), fPeakFinder->GetYPeak(k));
      track->SetMCid(label);
#endif
    }

    // cout << "Found " << fTracks[i]->GetNTracks() << " tracks in slice " << fCurrentSlice << endl;
    // fTracks[i]->QSort();
  }
}
