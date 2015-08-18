/// \file TransformerRow.cxx
/// \brief Implementation of the TransformerRow class
/// \author Cvetan Cheshkov <cvetan.cheshkov@cern.ch>

#include "Transform.h"
#include "HoughTrack.h"
#include "TransformerRow.h"

using namespace std;
using namespace AliceO2::Hough;

Float_t TransformerRow::fgBeta1 = 1.0 / Transform::Row2X(84);
Float_t TransformerRow::fgBeta2 =
  1.0 / (Transform::Row2X(158) * (1.0 + tan(Transform::Pi() * 10 / 180) * tan(Transform::Pi() * 10 / 180)));
Float_t TransformerRow::fgDAlpha = 0.22;
Float_t TransformerRow::fgDEta = 0.40;
Double_t TransformerRow::fgEtaCalcParam1 = 1.0289;
Double_t TransformerRow::fgEtaCalcParam2 = 0.15192;
Double_t TransformerRow::fgEtaCalcParam3 = 1. / (32. * 600. * 600.);

TransformerRow::TransformerRow()
{
  // Default constructor
  fParamSpace = 0;

  fGapCount = 0;
  fCurrentRowCount = 0;
#ifdef do_mc
  fTrackID = 0;
#endif
  fTrackNRows = 0;
  fTrackFirstRow = 0;
  fTrackLastRow = 0;
  fInitialGapCount = 0;

  fPrevBin = 0;
  fNextBin = 0;
  fNextRow = 0;

  fStartPadParams = 0;
  fEndPadParams = 0;
  fLUTr = 0;
  fLUTforwardZ = 0;
  fLUTbackwardZ = 0;
}

TransformerRow::TransformerRow(Int_t slice, Int_t patch, Int_t netasegments, Bool_t /*DoMC*/, Float_t zvertex)
  : BaseTransformer(slice, patch, netasegments, zvertex)
{
  // Normal constructor
  fParamSpace = 0;

  fGapCount = 0;
  fCurrentRowCount = 0;
#ifdef do_mc
  fTrackID = 0;
#endif

  fTrackNRows = 0;
  fTrackFirstRow = 0;
  fTrackLastRow = 0;
  fInitialGapCount = 0;

  fPrevBin = 0;
  fNextBin = 0;
  fNextRow = 0;

  fStartPadParams = 0;
  fEndPadParams = 0;
  fLUTr = 0;
  fLUTforwardZ = 0;
  fLUTbackwardZ = 0;
}

TransformerRow::~TransformerRow()
{
  // Destructor
  if (fLastTransformer) {
    return;
  }
  DeleteHistograms();
#ifdef do_mc
  if (fTrackID) {
    for (Int_t i = 0; i < GetNEtaSegments(); i++) {
      if (!fTrackID[i]) {
        continue;
      }
      delete fTrackID[i];
    }
    delete[] fTrackID;
    fTrackID = 0;
  }
#endif

  if (fGapCount) {
    for (Int_t i = 0; i < GetNEtaSegments(); i++) {
      if (!fGapCount[i]) {
        continue;
      }
      delete[] fGapCount[i];
    }
    delete[] fGapCount;
    fGapCount = 0;
  }
  if (fCurrentRowCount) {
    for (Int_t i = 0; i < GetNEtaSegments(); i++) {
      if (fCurrentRowCount[i]) {
        delete[] fCurrentRowCount[i];
      }
    }
    delete[] fCurrentRowCount;
    fCurrentRowCount = 0;
  }
  if (fTrackNRows) {
    delete[] fTrackNRows;
    fTrackNRows = 0;
  }
  if (fTrackFirstRow) {
    delete[] fTrackFirstRow;
    fTrackFirstRow = 0;
  }
  if (fTrackLastRow) {
    delete[] fTrackLastRow;
    fTrackLastRow = 0;
  }
  if (fInitialGapCount) {
    delete[] fInitialGapCount;
    fInitialGapCount = 0;
  }
  if (fPrevBin) {
    for (Int_t i = 0; i < GetNEtaSegments(); i++) {
      if (!fPrevBin[i]) {
        continue;
      }
      delete[] fPrevBin[i];
    }
    delete[] fPrevBin;
    fPrevBin = 0;
  }
  if (fNextBin) {
    for (Int_t i = 0; i < GetNEtaSegments(); i++) {
      if (!fNextBin[i]) {
        continue;
      }
      delete[] fNextBin[i];
    }
    delete[] fNextBin;
    fNextBin = 0;
  }
  if (fNextRow) {
    for (Int_t i = 0; i < GetNEtaSegments(); i++) {
      if (!fNextRow[i]) {
        continue;
      }
      delete[] fNextRow[i];
    }
    delete[] fNextRow;
    fNextRow = 0;
  }
  if (fStartPadParams) {
    for (Int_t i = Transform::GetFirstRow(0); i <= Transform::GetLastRow(5); i++) {
      if (!fStartPadParams[i]) {
        continue;
      }
      delete[] fStartPadParams[i];
    }
    delete[] fStartPadParams;
    fStartPadParams = 0;
  }
  if (fEndPadParams) {
    for (Int_t i = Transform::GetFirstRow(0); i <= Transform::GetLastRow(5); i++) {
      if (!fEndPadParams[i]) {
        continue;
      }
      delete[] fEndPadParams[i];
    }
    delete[] fEndPadParams;
    fEndPadParams = 0;
  }
  if (fLUTr) {
    for (Int_t i = Transform::GetFirstRow(0); i <= Transform::GetLastRow(5); i++) {
      if (!fLUTr[i]) {
        continue;
      }
      delete[] fLUTr[i];
    }
    delete[] fLUTr;
    fLUTr = 0;
  }
  if (fLUTforwardZ) {
    delete[] fLUTforwardZ;
    fLUTforwardZ = 0;
  }
  if (fLUTbackwardZ) {
    delete[] fLUTbackwardZ;
    fLUTbackwardZ = 0;
  }
}

void TransformerRow::DeleteHistograms()
{
  // Clean up
  if (!fParamSpace) {
    return;
  }
  for (Int_t i = 0; i < GetNEtaSegments(); i++) {
    if (!fParamSpace[i]) {
      continue;
    }
    delete fParamSpace[i];
  }
  delete[] fParamSpace;
}

struct AliHLTTrackLength {
  // Structure is used for temporarely storage of the LUT
  // which contains the track lengths associated to each hough
  // space bin
  Bool_t fIsFilled; // Is bin already filled?
  UInt_t fFirstRow; // First TPC row crossed by the track
  UInt_t fLastRow;  // Last TPC row crossed by the track
  Float_t fTrackPt; // Pt of the track
};

/// Create the histograms (parameter space)
/// nxbin = #bins in X
/// nybin = #bins in Y
/// xmin xmax ymin ymax = histogram limits in X and Y
void TransformerRow::CreateHistograms(Int_t nxbin, Float_t xmin, Float_t xmax, Int_t nybin, Float_t ymin, Float_t ymax)
{
  cout << "fLastTransformer: " << fLastTransformer << endl;
  if (fLastTransformer) {
    SetTransformerArrays((TransformerRow*)fLastTransformer);
    return;
  }
  fParamSpace = new Accumulator* [GetNEtaSegments()];

  cout << "GetNEtaSegments(): " << GetNEtaSegments() << endl;

  for (Int_t i = 0; i < GetNEtaSegments(); i++) {
    fParamSpace[i] = new Accumulator(nxbin, xmin, xmax, nybin, ymin, ymax);
  }
#ifdef do_mc
  {
    Accumulator* hist = fParamSpace[0];
    Int_t ncellsx = (hist->GetNbinsX() + 3) / 2;
    Int_t ncellsy = (hist->GetNbinsY() + 3) / 2;
    Int_t ncells = ncellsx * ncellsy;
    if (!fTrackID) {
      cout << "Transformer: Allocating " << GetNEtaSegments() * ncells * sizeof(AliHLTTrackIndex)
           << " bytes to fTrackID" << endl;
      fTrackID = new AliHLTTrackIndex* [GetNEtaSegments()];
      for (Int_t i = 0; i < GetNEtaSegments(); i++) {
        fTrackID[i] = new AliHLTTrackIndex[ncells];
      }
    }
  }
#endif
  Accumulator* hist = fParamSpace[0];


  Int_t ncells = (hist->GetNbinsX() + 2) * (hist->GetNbinsY() + 2);
  cout << "ncells: " << ncells << endl;

  if (!fGapCount) {
    cout << "Transformer: Allocating " << GetNEtaSegments() * ncells * sizeof(UChar_t) << " bytes to fGapCount" << endl;
    fGapCount = new UChar_t* [GetNEtaSegments()];
    for (Int_t i = 0; i < GetNEtaSegments(); i++) {
      fGapCount[i] = new UChar_t[ncells];
    }
  }

  if (!fCurrentRowCount) {
    cout << "Transformer: Allocating " << GetNEtaSegments() * ncells * sizeof(UChar_t) << " bytes to fCurrentRowCount"
         << endl;
    fCurrentRowCount = new UChar_t* [GetNEtaSegments()];
    for (Int_t i = 0; i < GetNEtaSegments(); i++) {
      fCurrentRowCount[i] = new UChar_t[ncells];
    }
  }
  if (!fPrevBin) {
    cout << "Transformer: Allocating " << GetNEtaSegments() * ncells * sizeof(UChar_t) << " bytes to fPrevBin" << endl;
    fPrevBin = new UChar_t* [GetNEtaSegments()];
    for (Int_t i = 0; i < GetNEtaSegments(); i++) {
      fPrevBin[i] = new UChar_t[ncells];
    }
  }
  if (!fNextBin) {
    cout << "Transformer: Allocating " << GetNEtaSegments() * ncells * sizeof(UChar_t) << " bytes to fNextBin" << endl;
    fNextBin = new UChar_t* [GetNEtaSegments()];
    for (Int_t i = 0; i < GetNEtaSegments(); i++) {
      fNextBin[i] = new UChar_t[ncells];
    }
  }
  Int_t ncellsy = hist->GetNbinsY() + 2;
  cout << "ncellsy: " << ncellsy << endl;

  if (!fNextRow) {
    cout << "Transformer: Allocating " << GetNEtaSegments() * ncellsy * sizeof(UChar_t) << " bytes to fNextRow" << endl;
    fNextRow = new UChar_t* [GetNEtaSegments()];
    for (Int_t i = 0; i < GetNEtaSegments(); i++) {
      fNextRow[i] = new UChar_t[ncellsy];
    }
  }

  if (!fTrackNRows) {
    cout << "Transformer: Allocating " << ncells * sizeof(UChar_t) << " bytes to fTrackNRows" << endl;
    fTrackNRows = new UChar_t[ncells];
    cout << "Transformer: Allocating " << ncells * sizeof(UChar_t) << " bytes to fTrackFirstRow" << endl;
    fTrackFirstRow = new UChar_t[ncells];
    cout << "Transformer: Allocating " << ncells * sizeof(UChar_t) << " bytes to fTrackLastRow" << endl;
    fTrackLastRow = new UChar_t[ncells];
    cout << "Transformer: Allocating " << ncells * sizeof(UChar_t) << " bytes to fInitialGapCount" << endl;
    fInitialGapCount = new UChar_t[ncells];

    HoughTrack track;
    Int_t xmin = hist->GetFirstXbin();
    Int_t xmax = hist->GetLastXbin();
    Int_t xmiddle = (hist->GetNbinsX() + 1) / 2;
    Int_t ymin = hist->GetFirstYbin();
    Int_t ymax = hist->GetLastYbin();
    Int_t nxbins = hist->GetNbinsX() + 2;
    Int_t nxgrid = (hist->GetNbinsX() + 3) / 2 + 1;
    Int_t nygrid = hist->GetNbinsY() + 3;

    cout << "xmin: " << xmin << " xmax: " << xmax << " xmiddle: " << xmiddle << " ymin: " << ymin << " ymax: " << ymax << " nxbins: " << nxbins << " nxgrid: " << nxgrid << " nygrid: " << nygrid << endl;
    AliHLTTrackLength* tracklength = new AliHLTTrackLength[nxgrid * nygrid];
    memset(tracklength, 0, nxgrid * nygrid * sizeof(AliHLTTrackLength));

    /*cout << "fTrackNRows[0 + 0 * 10002] = 255" << endl;
    cout << "fTrackNRows[1 + 0 * 10002] = 255" << endl;
    cout << "..." << endl;
    cout << "fTrackNRows[5000 + 0 * 10002] = 255" << endl;
    cout << "fTrackNRows[0 + 1 * 10002] = 255" << endl;
    cout << "..." << endl;
    cout << "fTrackNRows[5000 + 181 * 10002] = 255" << endl;
*/
    for (Int_t ybin = ymin - 1; ybin <= (ymax + 1); ybin++) {
      for (Int_t xbin = xmin - 1; xbin <= xmiddle; xbin++) {

        cout << "========== " << "ybin: " << ybin << " " << " xbin: " << xbin << " ==========" << endl;
        cout << " fTrackNRows[" << xbin << " + " << ybin << " * " << nxbins << "] = 255" << endl;
        fTrackNRows[xbin + ybin * nxbins] = 255;
        for (Int_t deltay = 0; deltay <= 1; deltay++) {
          for (Int_t deltax = 0; deltax <= 1; deltax++) {

            cout << "====== " << "deltay: " <<  deltay << " " << " deltax: " << deltax << " ======" << endl;

            AliHLTTrackLength* curtracklength = &tracklength[(xbin + deltax) + (ybin + deltay) * nxgrid];
            UInt_t maxfirstrow = 0;
            UInt_t maxlastrow = 0;
            Float_t maxtrackpt = 0;

            cout << "maxfirstrow: " << maxfirstrow << " maxlastrow: " << maxlastrow << " maxtrackpt: " << maxtrackpt << endl;
            if (curtracklength->fIsFilled) {
              maxfirstrow = curtracklength->fFirstRow;
              maxlastrow = curtracklength->fLastRow;
              maxtrackpt = curtracklength->fTrackPt;
              cout << "=== Current track is filled! ===" << endl;
              cout << "maxfirstrow: " << maxfirstrow << " maxlastrow: " << maxlastrow << " maxtrackpt: " << maxtrackpt << endl;
            } else {
              Float_t xtrack = hist->GetPreciseBinCenterX((Float_t)xbin + 0.5 * (Float_t)(2 * deltax - 1));
              Float_t ytrack = hist->GetPreciseBinCenterY((Float_t)ybin + 0.5 * (Float_t)(2 * deltay - 1));

              Float_t psi = atan((xtrack - ytrack) / (fgBeta1 - fgBeta2));
              Float_t kappa = 2.0 * (xtrack * cos(psi) - fgBeta1 * sin(psi));
              track.SetTrackParameters(kappa, psi, 1);
              maxtrackpt = track.GetPt();

              cout << "=== Current track is not filled! ===" << endl;
              cout << "xtrack: " << xtrack << " ytrack: " << ytrack << " psi: " << psi << " kappa: " << kappa << " maxtrackpt: " << maxtrackpt << " SolenoidField: " << 0.9 * 0.1 * Transform::GetSolenoidField() << endl;

              if (maxtrackpt < 0.9 * 0.1 * Transform::GetSolenoidField()) {
                maxfirstrow = maxlastrow = 0;
                curtracklength->fIsFilled = true;
                curtracklength->fFirstRow = maxfirstrow;
                curtracklength->fLastRow = maxlastrow;
                curtracklength->fTrackPt = maxtrackpt;
              } else {

                Bool_t firstrow = false;
                UInt_t curfirstrow = 0;
                UInt_t curlastrow = 0;

                Double_t centerx = track.GetCenterX();
                Double_t centery = track.GetCenterY();
                Double_t radius = track.GetRadius();

                cout << "centerx: " << centerx << " centery: " << centery << " radius: " << radius << endl;

                for (Int_t j = Transform::GetFirstRow(0); j <= Transform::GetLastRow(5); j++) {
                  Float_t hit[3];
                  //		      if(!track.GetCrossingPoint(j,hit)) continue;
                  hit[0] = Transform::Row2X(j);
                  Double_t aa = (hit[0] - centerx) * (hit[0] - centerx);
                  Double_t r2 = radius * radius;

                  if (aa > r2) {
                    continue;
                  }

                  Double_t aa2 = sqrt(r2 - aa);
                  Double_t y1 = centery + aa2;
                  Double_t y2 = centery - aa2;
                  hit[1] = y1;

                  cout << "aa2: " << aa2 << " y1: " << y1 << " y2: " << y2 << " hit[1]: " << hit[1] << " maxtrackpt: " << maxtrackpt << endl;

                  if (fabs(y2) < fabs(y1)) {
                    hit[1] = y2;
                  }

                  hit[2] = 0;

                  Transform::LocHLT2Raw(hit, 0, j);
                  hit[1] += 0.5;
                  if (hit[1] >= 0 && hit[1] < Transform::GetNPads(j)) {
                    if (!firstrow) {
                      curfirstrow = j;
                      firstrow = true;
                    }
                    curlastrow = j;
                  } else {
                    if (firstrow) {
                      firstrow = false;
                      if ((curlastrow - curfirstrow) >= (maxlastrow - maxfirstrow)) {
                        maxfirstrow = curfirstrow;
                        maxlastrow = curlastrow;
                      }
                    }
                  }
                }
                if ((curlastrow - curfirstrow) >= (maxlastrow - maxfirstrow)) {
                  maxfirstrow = curfirstrow;
                  maxlastrow = curlastrow;
                }

                curtracklength->fIsFilled = true;
                curtracklength->fFirstRow = maxfirstrow;
                curtracklength->fLastRow = maxlastrow;
                curtracklength->fTrackPt = maxtrackpt;
              }
            }
            if ((maxlastrow - maxfirstrow) < fTrackNRows[xbin + ybin * nxbins]) {
              fTrackNRows[xbin + ybin * nxbins] = maxlastrow - maxfirstrow;
              fInitialGapCount[xbin + ybin * nxbins] = 1;

              if ((maxlastrow - maxfirstrow + 1) <= MIN_TRACK_LENGTH) {
                fInitialGapCount[xbin + ybin * nxbins] = MAX_N_GAPS + 1;
              }

              if (maxtrackpt < 0.9 * 0.1 * Transform::GetSolenoidField()) {
                fInitialGapCount[xbin + ybin * nxbins] = MAX_N_GAPS;
              }

              fTrackFirstRow[xbin + ybin * nxbins] = maxfirstrow;
              fTrackLastRow[xbin + ybin * nxbins] = maxlastrow;

              Int_t xmirror = xmax - xbin + 1;
              Int_t ymirror = ymax - ybin + 1;
              fTrackNRows[xmirror + ymirror * nxbins] = fTrackNRows[xbin + ybin * nxbins];
              fInitialGapCount[xmirror + ymirror * nxbins] = fInitialGapCount[xbin + ybin * nxbins];
              fTrackFirstRow[xmirror + ymirror * nxbins] = fTrackFirstRow[xbin + ybin * nxbins];
              fTrackLastRow[xmirror + ymirror * nxbins] = fTrackLastRow[xbin + ybin * nxbins];
            }
          }
        }
        //	      cout<<" fTrackNRows "<<(Int_t)fInitialGapCount[xbin + ybin*nxbins]<<" "<<xbin<<" "<<ybin<<"
        //"<<(Int_t)fTrackNRows[xbin + ybin*nxbins]<<" "<<(Int_t)fTrackFirstRow[xbin + ybin*nxbins]<<"
        //"<<(Int_t)fTrackLastRow[xbin + ybin*nxbins]<<" "<<endl;
      }
    }
    delete[] tracklength;
  }

  if (!fStartPadParams) {
    Int_t nrows = Transform::GetLastRow(5) - Transform::GetFirstRow(0) + 1;
    cout << "Transformer: Allocating about " << nrows * 100 * sizeof(AliHLTPadHoughParams)
         << " bytes to fStartPadParams" << endl;
    fStartPadParams = new AliHLTPadHoughParams* [nrows];
    cout << "Transformer: Allocating about " << nrows * 100 * sizeof(AliHLTPadHoughParams) << " bytes to fEndPadParams"
         << endl;
    fEndPadParams = new AliHLTPadHoughParams* [nrows];
    cout << "Transformer: Allocating about " << nrows * 100 * sizeof(Float_t) << " bytes to fLUTr" << endl;
    fLUTr = new Float_t* [nrows];

    Float_t beta1 = fgBeta1;
    Float_t beta2 = fgBeta2;
    Float_t beta1minusbeta2 = fgBeta1 - fgBeta2;
    Float_t ymin = hist->GetYmin();
    Float_t histbin = hist->GetBinWidthY();
    Float_t xmin = hist->GetXmin();
    Float_t xmax = hist->GetXmax();
    Float_t xbin = (xmax - xmin) / hist->GetNbinsX();
    Int_t firstbinx = hist->GetFirstXbin();
    Int_t lastbinx = hist->GetLastXbin();
    Int_t nbinx = hist->GetNbinsX() + 2;
    Int_t firstbin = hist->GetFirstYbin();
    Int_t lastbin = hist->GetLastYbin();
    for (Int_t i = Transform::GetFirstRow(0); i <= Transform::GetLastRow(5); i++) {
      cout << "========== " << "Pad row: " << i << " ==========" << endl;
      Int_t npads = Transform::GetNPads(i);
      Int_t ipatch = Transform::GetPatch(i);
      Double_t padpitch = Transform::GetPadPitchWidth(ipatch);
      Float_t x = Transform::Row2X(i);
      Float_t x2 = x * x;

      fStartPadParams[i] = new AliHLTPadHoughParams[npads];
      fEndPadParams[i] = new AliHLTPadHoughParams[npads];
      fLUTr[i] = new Float_t[npads];
      for (Int_t pad = 0; pad < npads; pad++) {

        cout << "=== " << "Pad: " << pad << " ===" << endl;

        Float_t y = (pad - 0.5 * (npads - 1)) * padpitch;
        fLUTr[i][pad] = sqrt(x2 + y * y);
        Float_t starty = (pad - 0.5 * npads) * padpitch;
        Float_t r1 = x2 + starty * starty;
        Float_t xoverr1 = x / r1;
        Float_t startyoverr1 = starty / r1;
        Float_t endy = (pad - 0.5 * (npads - 2)) * padpitch;
        Float_t r2 = x2 + endy * endy;
        Float_t xoverr2 = x / r2;
        Float_t endyoverr2 = endy / r2;
        Float_t a1 = beta1minusbeta2 / (xoverr1 - beta2);
        Float_t b1 = (xoverr1 - beta1) / (xoverr1 - beta2);
        Float_t a2 = beta1minusbeta2 / (xoverr2 - beta2);
        Float_t b2 = (xoverr2 - beta1) / (xoverr2 - beta2);

        Float_t alpha1 = (a1 * startyoverr1 + b1 * ymin - xmin) / xbin;





        Float_t deltaalpha1 = b1 * histbin / xbin;
        if (b1 < 0) {
          alpha1 += deltaalpha1;
        }

        cout << "beta1: " << beta1 << " beta2: " << beta2 << " beta1minusbeta2: " << beta1minusbeta2 << " ymin: " << ymin << " histbin: " << histbin
             << "xmin: " << xmin << " xmax: " << xmax << " xbin: " << xbin << " firstbinx: " << firstbinx << " lastbinx: " << lastbinx
             << "nbinx: " << nbinx << " firstbin: " << firstbin << " lastbin: " << lastbin << " npads: " << npads << " ipatch: " << ipatch
             << "padpitch: " << padpitch << " x: " << x << " x2: " << x2 << " y: " << y << " fLUTr[i][pad]: " << fLUTr[i][pad]
             << "starty: " << starty << " r1: " << r1 << " xoverr1: " << xoverr1 << " startyoverr1: " << startyoverr1 << " endy: " << endy
             << "r2: " << r2 << " xoverr2: " << xoverr2 << " endyoverr2: " << endyoverr2 << " a1: " << a1 << " b1: " << b1
             << "a2: " << a2 << " b2: " << b2 << " alpha1: " << alpha1 << " deltaalpha1: " << deltaalpha1 << endl;

        Float_t alpha2 = (a2 * endyoverr2 + b2 * ymin - xmin) / xbin;
        Float_t deltaalpha2 = b2 * histbin / xbin;
        if (b2 >= 0) {
          alpha2 += deltaalpha2;
        }

        fStartPadParams[i][pad].fAlpha = alpha1;
        fStartPadParams[i][pad].fDeltaAlpha = deltaalpha1;
        fEndPadParams[i][pad].fAlpha = alpha2;
        fEndPadParams[i][pad].fDeltaAlpha = deltaalpha2;

        // Find the first and last bin rows to be filled
        Bool_t binfound1 = false;
        Bool_t binfound2 = false;
        Int_t firstbin1 = lastbin;
        Int_t firstbin2 = lastbin;
        Int_t lastbin1 = firstbin;
        Int_t lastbin2 = firstbin;
        for (Int_t b = firstbin; b <= lastbin; b++, alpha1 += deltaalpha1, alpha2 += deltaalpha2) {
          Int_t binx1 = 1 + (Int_t)alpha1;
//          cout << "binx1: " << binx1 << " binx2: " << 1 + (Int_t)alpha2 << " firstbinx: " << firstbinx
//               << " lastbinx: " << lastbinx << endl;
          if (binx1 <= lastbinx) {
            UChar_t initialgapcount;
            if (binx1 >= firstbinx) {
              initialgapcount = fInitialGapCount[binx1 + b * nbinx];
            } else {
              initialgapcount = fInitialGapCount[firstbinx + b * nbinx];
            }
            if (initialgapcount != MAX_N_GAPS) {
              if (!binfound1) {
                firstbin1 = b;
                binfound1 = true;
              }
              lastbin1 = b;
            }
          }
          Int_t binx2 = 1 + (Int_t)alpha2;
          if (binx2 >= firstbin) {
            UChar_t initialgapcount;
            if (binx2 <= lastbinx) {
              initialgapcount = fInitialGapCount[binx2 + b * nbinx];
            } else {
              initialgapcount = fInitialGapCount[lastbinx + b * nbinx];
            }
            if (initialgapcount != MAX_N_GAPS) {
              if (!binfound2) {
                firstbin2 = b;
                binfound2 = true;
              }
              lastbin2 = b;
            }
          }
        }
        fStartPadParams[i][pad].fFirstBin = firstbin1;
        fStartPadParams[i][pad].fLastBin = lastbin1;
        fEndPadParams[i][pad].fFirstBin = firstbin2;
        fEndPadParams[i][pad].fLastBin = lastbin2;
      }
    }
  }

  // create lookup table for z of the digits
  if (!fLUTforwardZ) {
    Int_t ntimebins = Transform::GetNTimeBins();
    cout << "Transformer: Allocating " << ntimebins * sizeof(Float_t) << " bytes to fLUTforwardZ" << endl;
    fLUTforwardZ = new Float_t[ntimebins];
    cout << "Transformer: Allocating " << ntimebins * sizeof(Float_t) << " bytes to fLUTbackwardZ" << endl;
    fLUTbackwardZ = new Float_t[ntimebins];
    for (Int_t i = 0; i < ntimebins; i++) {
      Float_t z;
      z = Transform::GetZFast(0, i, GetZVertex());
      fLUTforwardZ[i] = z;
      z = Transform::GetZFast(18, i, GetZVertex());
      fLUTbackwardZ[i] = z;
    }
  }
}

void TransformerRow::Reset()
{
  // Reset all the histograms. Should be done when processing new slice
  if (fLastTransformer) {
    return;
  }

  if (!fParamSpace) {
    cout << "No histograms to reset" << endl;
    return;
  }

  for (Int_t i = 0; i < GetNEtaSegments(); i++) {
    fParamSpace[i]->Reset();
  }

#ifdef do_mc
  {
    Accumulator* hist = fParamSpace[0];
    Int_t ncellsx = (hist->GetNbinsX() + 3) / 2;
    Int_t ncellsy = (hist->GetNbinsY() + 3) / 2;
    Int_t ncells = ncellsx * ncellsy;
    for (Int_t i = 0; i < GetNEtaSegments(); i++)
      memset(fTrackID[i], 0, ncells * sizeof(AliHLTTrackIndex));
  }
#endif
  Accumulator* hist = fParamSpace[0];
  Int_t ncells = (hist->GetNbinsX() + 2) * (hist->GetNbinsY() + 2);
  for (Int_t i = 0; i < GetNEtaSegments(); i++) {
    memcpy(fGapCount[i], fInitialGapCount, ncells * sizeof(UChar_t));
    memcpy(fCurrentRowCount[i], fTrackFirstRow, ncells * sizeof(UChar_t));
  }
}

Int_t TransformerRow::GetEtaIndex(Double_t eta) const
{
  Double_t etaslice = (GetEtaMax() - GetEtaMin()) / GetNEtaSegments();
  Double_t index = (eta - GetEtaMin()) / etaslice;
  return (Int_t)index;
}

inline Accumulator* TransformerRow::GetHistogram(Int_t etaindex)
{
  if (!fParamSpace || etaindex >= GetNEtaSegments() || etaindex < 0) {
    return 0;
  }

  if (!fParamSpace[etaindex]) {
    return 0;
  }

  return fParamSpace[etaindex];
}

Double_t TransformerRow::GetEta(Int_t etaindex, Int_t /*slice*/) const
{
  Double_t etaslice = (GetEtaMax() - GetEtaMin()) / GetNEtaSegments();
  Double_t eta = 0;
  eta = (Double_t)((etaindex + 0.5) * etaslice);
  return eta;
}

void TransformerRow::TransformCircle(ClusterCollection* clusterCollection)
{
  TransformCircleFromDigitArray(clusterCollection);
}

// Do the Hough Transform reading the data from the HLT binary file
void TransformerRow::TransformCircleFromDigitArray(ClusterCollection* clusterCollection)
{
  // Load the parameters used by the fast calculation of eta
  Double_t etaparam1 = GetEtaCalcParam1();
  Double_t etaparam2 = GetEtaCalcParam2();
  Double_t etaparam3 = GetEtaCalcParam3();

  Int_t netasegments = GetNEtaSegments();
  Double_t etamin = GetEtaMin();
  Double_t etaslice = (GetEtaMax() - etamin) / netasegments;

  Int_t lowerthreshold = GetLowerThreshold();

  // Assumes that all the etaslice histos are the same!
  Accumulator* h = fParamSpace[0];
  Int_t firstbiny = h->GetFirstYbin();
  Int_t firstbinx = h->GetFirstXbin();
  Int_t lastbinx = h->GetLastXbin();
  Int_t nbinx = h->GetNbinsX() + 2;

  UChar_t lastpad;
  Int_t lastetaindex = -1;
  AliHLTEtaRow* etaclust = new AliHLTEtaRow[netasegments];
  /*
    AliHLTDigitRowData* tempPt = GetDataPointer();
    if (!tempPt) {
      LOG(AliHLTLog::kError, "AliHLTHoughTransformer::TransformCircle", "Data") << "No input data " << endl;
      return;
    }*/

  Int_t ipatch = GetPatch();
  Int_t ilastpatch = GetLastPatch();
  Int_t islice = GetSlice();
  Float_t* lutz;

  if (islice < 18) {
    lutz = fLUTforwardZ;
  } else {
    lutz = fLUTbackwardZ;
  }

  // Loop over the padrows:
  for (UChar_t row = Transform::GetFirstRow(ipatch); row <= Transform::GetLastRow(ipatch); row++) {
    lastpad = 255;
    // Flush eta clusters array
    memset(etaclust, 0, netasegments * sizeof(AliHLTEtaRow));

    Float_t radius = 0;

    // Loop over the data on this padrow
    for (UInt_t cluster = 0; cluster < clusterCollection->getNumberOfClustersPerPadRow(row); cluster++) {
      UShort_t charge = clusterCollection->getClusterCharge(row, cluster);

      if ((Int_t)charge <= lowerthreshold) {
        continue;
      }

      UChar_t pad = clusterCollection->getClusterPad(row, cluster);
      UShort_t time = clusterCollection->getClusterTime(row, cluster);

      if (pad != lastpad) {
        radius = fLUTr[(Int_t)row][(Int_t)pad];
        lastetaindex = -1;
      }

      Float_t z = lutz[(Int_t)time];
      Double_t radiuscorr = radius * (1. + etaparam3 * radius * radius);
      Double_t zovr = z / radiuscorr;
      Double_t eta = (etaparam1 - etaparam2 * fabs(zovr)) * zovr;
      // Get the corresponding index, which determines which histogram to fill:
      Int_t etaindex = (Int_t)((eta - etamin) / etaslice);

#ifndef do_mc
      if (etaindex == lastetaindex) {
        continue;
      }
#endif

      /*      cout << "Cluster at patch " << ipatch << " row " << (Int_t)row << " pad " << (Int_t)pad << " time " <<
         time
                 << " x: " << clusterCollection->getClusterX(row, cluster)
                 << " y: " << clusterCollection->getClusterY(row, cluster) << " z: " << z
                 << " a: " << clusterCollection->getClusterAlpha(row, cluster)
                 << " b: " << clusterCollection->getClusterBeta(row, cluster)
                 << " etaindex " << etaindex
                 << endl; */

      if (etaindex < 0 || etaindex >= netasegments) {
        continue;
      }

      if (!etaclust[etaindex].fIsFound) {
        etaclust[etaindex].fStartPad = pad;
        etaclust[etaindex].fEndPad = pad;
        etaclust[etaindex].fIsFound = 1;
        continue;
      } else {
        if (pad <= (etaclust[etaindex].fEndPad + 1)) {
          etaclust[etaindex].fEndPad = pad;
        } else {
          FillCluster(row, etaindex, etaclust, ilastpatch, firstbinx, lastbinx, nbinx, firstbiny);

          etaclust[etaindex].fStartPad = pad;
          etaclust[etaindex].fEndPad = pad;
        }
      }
      lastpad = pad;
      lastetaindex = etaindex;
    }
    // Write remaining clusters
    for (Int_t etaindex = 0; etaindex < netasegments; etaindex++) {
      // Check for empty row
      if ((etaclust[etaindex].fStartPad == 0) && (etaclust[etaindex].fEndPad == 0)) {
        continue;
      }

      FillCluster(row, etaindex, etaclust, ilastpatch, firstbinx, lastbinx, nbinx, firstbiny);
    }
  }

  delete[] etaclust;
}

/*
void TransformerRow::TransformCircleFromRawStream()
{
  // Do the Hough Transform

  // Load the parameters used by the fast calculation of eta
  Double_t etaparam1 = GetEtaCalcParam1();
  Double_t etaparam2 = GetEtaCalcParam2();
  Double_t etaparam3 = GetEtaCalcParam3();

  Int_t netasegments = GetNEtaSegments();
  Double_t etamin = GetEtaMin();
  Double_t etaslice = (GetEtaMax() - etamin) / netasegments;

  Int_t lowerthreshold = GetLowerThreshold();

  // Assumes that all the etaslice histos are the same!
  Accumulator* h = fParamSpace[0];
  Int_t firstbiny = h->GetFirstYbin();
  Int_t firstbinx = h->GetFirstXbin();
  Int_t lastbinx = h->GetLastXbin();
  Int_t nbinx = h->GetNbinsX() + 2;

  Int_t lastetaindex = -1;
  AliHLTEtaRow* etaclust = new AliHLTEtaRow[netasegments];

  if (!fTPCRawStream) {
    cerr << "No input data " << endl;
    return;
  }

  Int_t ipatch = GetPatch();
  UChar_t rowmin = AliceO2::Hough::Transform::GetFirstRowOnDDL(ipatch);
  UChar_t rowmax = AliceO2::Hough::Transform::GetLastRowOnDDL(ipatch);
  //  Int_t ntimebins = AliHLTTransform::GetNTimeBins();
  Int_t ilastpatch = GetLastPatch();
  Int_t islice = GetSlice();
  Float_t* lutz;
  if (islice < 18)
    lutz = fLUTforwardZ;
  else
    lutz = fLUTbackwardZ;

  // Flush eta clusters array
  memset(etaclust, 0, netasegments * sizeof(AliHLTEtaRow));

  UChar_t i = 0;
  Int_t npads = 0;
  Float_t radius = 0;
  UChar_t pad = 0;

  // Loop over the rawdata:
  while (fTPCRawStream->Next()) {

    if (fTPCRawStream->IsNewSector() || fTPCRawStream->IsNewRow()) {

      // Write remaining clusters
      for (Int_t etaindex = 0; etaindex < netasegments; etaindex++) {
        // Check for empty row
        if ((etaclust[etaindex].fStartPad == 0) && (etaclust[etaindex].fEndPad == 0))
          continue;

        FillCluster(i, etaindex, etaclust, ilastpatch, firstbinx, lastbinx, nbinx, firstbiny);
      }

      Int_t sector = fTPCRawStream->GetSector();
      Int_t row = fTPCRawStream->GetRow();
      Int_t slice, srow;
      AliceO2::Hough::Transform::Sector2Slice(slice, srow, sector, row);
      if (slice != islice) {
        cerr << "Found slice " << slice << ", expected " << islice << endl;
        continue;
      }

      i = (UChar_t)srow;
      npads = AliceO2::Hough::Transform::GetNPads(srow) - 1;

      // Flush eta clusters array
      memset(etaclust, 0, netasegments * sizeof(AliHLTEtaRow));

      radius = 0;
    }

    if ((i < rowmin) || (i > rowmax))
      continue;

    //    cout<<" Starting row "<<(UInt_t)i<<endl;
    // Loop over the data on this padrow:
    if (fTPCRawStream->IsNewRow() || fTPCRawStream->IsNewPad()) {
      pad = fTPCRawStream->GetPad();

      radius = fLUTr[(Int_t)i][(Int_t)pad];
      lastetaindex = -1;
    }

    UShort_t time = fTPCRawStream->GetTime();

    if (fTPCRawStream->GetSignal() <= lowerthreshold)
      continue;

    Float_t z = lutz[(Int_t)time];
    Double_t radiuscorr = radius * (1. + etaparam3 * radius * radius);
    Double_t zovr = z / radiuscorr;
    Double_t eta = (etaparam1 - etaparam2 * fabs(zovr)) * zovr;
    // Get the corresponding index, which determines which histogram to fill:
    Int_t etaindex = (Int_t)((eta - etamin) / etaslice);

#ifndef do_mc
    if (etaindex == lastetaindex)
      continue;
#endif
    //    cout<<" Digit at patch "<<ipatch<<" row "<<i<<" pad "<<(Int_t)pad<<" time "<<time<<" etaslice
    //    "<<etaindex<<endl;

    if (etaindex < 0 || etaindex >= netasegments)
      continue;

    if (!etaclust[etaindex].fIsFound) {
      etaclust[etaindex].fStartPad = pad;
      etaclust[etaindex].fEndPad = pad;
      etaclust[etaindex].fIsFound = 1;
      continue;
    } else {
      if (pad <= (etaclust[etaindex].fEndPad + 1)) {
        etaclust[etaindex].fEndPad = pad;
      } else {
        FillCluster(i, etaindex, etaclust, ilastpatch, firstbinx, lastbinx, nbinx, firstbiny);

        etaclust[etaindex].fStartPad = pad;
        etaclust[etaindex].fEndPad = pad;
      }
    }
    lastetaindex = etaindex;
  }

  // Write remaining clusters
  for (Int_t etaindex = 0; etaindex < netasegments; etaindex++) {
    // Check for empty row
    if ((etaclust[etaindex].fStartPad == 0) && (etaclust[etaindex].fEndPad == 0))
      continue;

    FillCluster(i, etaindex, etaclust, ilastpatch, firstbinx, lastbinx, nbinx, firstbiny);
  }

  delete[] etaclust;
}
*/

#ifndef do_mc
Int_t TransformerRow::GetTrackID(Int_t /*etaindex*/, Double_t /*alpha1*/, Double_t /*alpha2*/) const
{
  // Does nothing if do_mc undefined
  cout << "Flag switched off" << endl;
  return -1;
#else
Int_t TransformerRow::GetTrackID(Int_t etaindex, Double_t alpha1, Double_t alpha2) const
{
  // Returns the MC label for a given peak found in the Hough space
  if (etaindex < 0 || etaindex > GetNEtaSegments()) {
    cout << "Wrong etaindex " << etaindex << endl;
    return -1;
  }
  Accumulator* hist = fParamSpace[etaindex];
  Int_t bin = hist->FindLabelBin(alpha1, alpha2);
  if (bin == -1) {
    cout << "Track candidate outside Hough space boundaries: " << alpha1 << " " << alpha2 << endl;
    return -1;
  }
  Int_t label = -1;
  Int_t max = 0;
  for (UInt_t i = 0; i < (MaxTrack - 1); i++) {
    Int_t nhits = fTrackID[etaindex][bin].fNHits[i];
    if (nhits == 0) {
      break;
    }
    if (nhits > max) {
      max = nhits;
      label = fTrackID[etaindex][bin].fLabel[i];
    }
  }
  Int_t label2 = -1;
  Int_t max2 = 0;
  for (UInt_t i = 0; i < (MaxTrack - 1); i++) {
    Int_t nhits = fTrackID[etaindex][bin].fNHits[i];
    if (nhits == 0) {
      break;
    }
    if (nhits > max2) {
      if (fTrackID[etaindex][bin].fLabel[i] != label) {
        max2 = nhits;
        label2 = fTrackID[etaindex][bin].fLabel[i];
      }
    }
  }
  if (max2 != 0) {
    cout << " TrackID"
         << " label " << label << " max " << max << " label2 " << label2 << " max2 " << max2 << " "
         << (Float_t)max2 / (Float_t)max << " " << fTrackID[etaindex][bin].fLabel[MaxTrack - 1] << " "
         << (Int_t)fTrackID[etaindex][bin].fNHits[MaxTrack - 1] << endl;
  }
  return label;
#endif
}

Int_t TransformerRow::GetTrackLength(Double_t alpha1, Double_t alpha2, Int_t* rows) const
{
  // Returns the track length for a given peak found in the Hough space

  Accumulator* hist = fParamSpace[0];
  Int_t bin = hist->FindBin(alpha1, alpha2);
  if (bin == -1) {
    cout << "Track candidate outside Hough space boundaries: " << alpha1 << " " << alpha2 << endl;
    return -1;
  }
  rows[0] = fTrackFirstRow[bin];
  rows[1] = fTrackLastRow[bin];

  return 0;
}

inline void TransformerRow::FillClusterRow(UChar_t i, Int_t binx1, Int_t binx2, UChar_t* ngaps2, UChar_t* currentrow2,
                                           UChar_t* lastrow2
#ifdef do_mc
                                           ,
                                           AliHLTEtaRow etaclust, AliHLTTrackIndex* trackid
#endif
                                           )
{
  for (Int_t bin = binx1; bin <= binx2; bin++) {
    if (ngaps2[bin] < MAX_N_GAPS) {
      if (i < lastrow2[bin] && i > currentrow2[bin]) {
        ngaps2[bin] += (i - currentrow2[bin] - 1);
        currentrow2[bin] = i;
      }
#ifdef do_mc
      if (i < lastrow2[bin] && i >= currentrow2[bin]) {
        for (UInt_t t = 0; t < (MaxTrack - 1); t++) {
          Int_t label = etaclust.fMcLabels[t];
          if (label == 0) {
            break;
          }
          UInt_t c;
          Int_t tempbin2 = (Int_t)(bin / 2);
          for (c = 0; c < (MaxTrack - 1); c++) {
            if (trackid[tempbin2].fLabel[c] == label || trackid[tempbin2].fNHits[c] == 0) {
              break;
            }
          }
          trackid[tempbin2].fLabel[c] = label;
          if (trackid[tempbin2].fCurrentRow[c] != i) {
            trackid[tempbin2].fNHits[c]++;
            trackid[tempbin2].fCurrentRow[c] = i;
          }
        }
      }
#endif
    }
  }
}

inline void TransformerRow::FillCluster(UChar_t i, Int_t etaindex, AliHLTEtaRow* etaclust, Int_t ilastpatch,
                                        Int_t firstbinx, Int_t lastbinx, Int_t nbinx, Int_t firstbiny)
{
  UChar_t* ngaps = fGapCount[etaindex];
  UChar_t* currentrow = fCurrentRowCount[etaindex];
  UChar_t* lastrow = fTrackLastRow;
  UChar_t* prevbin = fPrevBin[etaindex];
  UChar_t* nextbin = fNextBin[etaindex];
  UChar_t* nextrow = fNextRow[etaindex];
#ifdef do_mc
  AliHLTTrackIndex* trackid = fTrackID[etaindex];
#endif

  // Do the transformation:
  AliHLTPadHoughParams* startparams = &fStartPadParams[(Int_t)i][etaclust[etaindex].fStartPad];
  AliHLTPadHoughParams* endparams = &fEndPadParams[(Int_t)i][etaclust[etaindex].fEndPad];

  Float_t alpha1 = startparams->fAlpha;
  Float_t deltaalpha1 = startparams->fDeltaAlpha;
  Float_t alpha2 = endparams->fAlpha;
  Float_t deltaalpha2 = endparams->fDeltaAlpha;

  // cout << "row: " << (int)i <<" etaIndex: " << etaindex << " etaclustPad: " << (int)etaclust->fStartPad << " alpha1:
  // " << alpha1 << " alpha2: " << alpha2 << " deltaalpha1: " << deltaalpha1 << " deltaalpha2: " << deltaalpha1 << endl;

  Int_t firstbin1 = startparams->fFirstBin;
  Int_t firstbin2 = endparams->fFirstBin;
  Int_t firstbin = firstbin1;
  if (firstbin > firstbin2)
    firstbin = firstbin2;

  Int_t lastbin1 = startparams->fLastBin;
  Int_t lastbin2 = endparams->fLastBin;
  Int_t lastbin = lastbin1;
  if (lastbin < lastbin2) {
    lastbin = lastbin2;
  }

  alpha1 += (firstbin - firstbiny) * deltaalpha1;
  alpha2 += (firstbin - firstbiny) * deltaalpha2;

//  cout << "firstbin: " << firstbin << " firstbin1: " << firstbin1 << " firstbin2: " << firstbin2
//       << " lastbin: " << lastbin << " lastbin1: " << lastbin1 << " lastbin2: " << lastbin2 << endl;

  // Fill the histogram along the alpha2 range
  if (ilastpatch == -1) {
    for (Int_t b = firstbin; b <= lastbin; b++, alpha1 += deltaalpha1, alpha2 += deltaalpha2) {
      Int_t binx1 = 1 + (Int_t)alpha1;
      if (binx1 > lastbinx) {
        continue;
      }
      if (binx1 < firstbinx) {
        binx1 = firstbinx;
      }
      Int_t binx2 = 1 + (Int_t)alpha2;
      if (binx2 < firstbinx) {
        continue;
      }
      if (binx2 > lastbinx) {
        binx2 = lastbinx;
      }
#ifdef do_mc
      if (binx2 < binx1) {
        cout << "Wrong filling " << binx1 << " " << binx2 << " " << i << " " << etaclust[etaindex].fStartPad << " "
             << etaclust[etaindex].fEndPad << endl;
      }
#endif
      Int_t tempbin = b * nbinx;
      UChar_t* ngaps2 = ngaps + tempbin;
      UChar_t* currentrow2 = currentrow + tempbin;
      UChar_t* lastrow2 = lastrow + tempbin;
#ifdef do_mc
      Int_t tempbin2 = ((Int_t)(b / 2)) * ((Int_t)((nbinx + 1) / 2));
      AliHLTTrackIndex* trackid2 = trackid + tempbin2;
#endif
      FillClusterRow(i, binx1, binx2, ngaps2, currentrow2, lastrow2
#ifdef do_mc
                     ,
                     etaclust[etaindex], trackid2
#endif
                     );
    }
  } else {
    for (Int_t b = firstbin; b <= lastbin; b++, alpha1 += deltaalpha1, alpha2 += deltaalpha2) {
      Int_t binx1 = 1 + (Int_t)alpha1;
      if (binx1 > lastbinx) {
        continue;
      }
      if (binx1 < firstbinx) {
        binx1 = firstbinx;
      }
      Int_t binx2 = 1 + (Int_t)alpha2;
      if (binx2 < firstbinx) {
        continue;
      }
      if (binx2 > lastbinx) {
        binx2 = lastbinx;
      }
#ifdef do_mc
      if (binx2 < binx1) {
        cout << "Wrong filling " << binx1 << " " << binx2 << " " << i << " " << etaclust[etaindex].fStartPad << " "
             << etaclust[etaindex].fEndPad << endl;
      }
#endif
      if (nextrow[b] > b) {
        Int_t deltab = (nextrow[b] - b - 1);
        b += deltab;
        alpha1 += deltaalpha1 * deltab;
        alpha2 += deltaalpha2 * deltab;
        continue;
      }
      Int_t tempbin = b * nbinx;
      binx1 = (UInt_t)nextbin[binx1 + tempbin];
      binx2 = (UInt_t)prevbin[binx2 + tempbin];
      if (binx2 < binx1) {
        continue;
      }
      UChar_t* ngaps2 = ngaps + tempbin;
      UChar_t* currentrow2 = currentrow + tempbin;
      UChar_t* lastrow2 = lastrow + tempbin;
#ifdef do_mc
      Int_t tempbin2 = ((Int_t)(b / 2)) * ((Int_t)((nbinx + 1) / 2));
      AliHLTTrackIndex* trackid2 = trackid + tempbin2;
#endif
      FillClusterRow(i, binx1, binx2, ngaps2, currentrow2, lastrow2
#ifdef do_mc
                     ,
                     etaclust[etaindex], trackid2
#endif
                     );
    }
  }
}

#ifdef do_mc
inline void TransformerRow::FillClusterMCLabels(AliHLTDigitData digpt, AliHLTEtaRow* etaclust)
{
  for (Int_t t = 0; t < 3; t++) {
    Int_t label = digpt.fTrackID[t];
    if (label < 0) {
      break;
    }
    UInt_t c;
    for (c = 0; c < (MaxTrack - 1); c++) {
      if (etaclust->fMcLabels[c] == label || etaclust->fMcLabels[c] == 0) {
        break;
      }
    }

    etaclust->fMcLabels[c] = label;
  }
}
#endif

void TransformerRow::SetTransformerArrays(TransformerRow* tr)
{
  fGapCount = tr->fGapCount;
  fCurrentRowCount = tr->fCurrentRowCount;
#ifdef do_mc
  fTrackID = tr->fTrackID;
#endif
  fTrackNRows = tr->fTrackNRows;
  fTrackFirstRow = tr->fTrackFirstRow;
  fTrackLastRow = tr->fTrackLastRow;
  fInitialGapCount = tr->fInitialGapCount;

  fPrevBin = tr->fPrevBin;
  fNextBin = tr->fNextBin;
  fNextRow = tr->fNextRow;

  fStartPadParams = tr->fStartPadParams;
  fEndPadParams = tr->fEndPadParams;
  fLUTr = tr->fLUTr;
  fLUTforwardZ = tr->fLUTforwardZ;
  fLUTbackwardZ = tr->fLUTbackwardZ;

  fParamSpace = tr->fParamSpace;

  return;
}
