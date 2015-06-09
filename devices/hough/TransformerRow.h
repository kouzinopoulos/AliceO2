/// \file TransformerRow.h
/// \brief Definition of the TransformerRow class
/// \author Anders Vestbo <vestbo@fi.uib.no>

#ifndef ALICEO2_HOUGH_TRANSFORMERROW_H_
#define ALICEO2_HOUGH_TRANSFORMERROW_H_

#include "Accumulator.h"
#include "BaseTransformer.h"
#include "ClusterCollection.h"
#include "StandardIncludes.h"

#define MAX_N_GAPS 5
#define MIN_TRACK_LENGTH 70

namespace AliceO2 {
namespace Hough {

struct AliHLTEtaRow {
  UChar_t fStartPad; ///< First pad in the cluster
  UChar_t fEndPad;   ///< Last pad in the cluster
  Bool_t fIsFound;   ///< Is the cluster already found
#ifdef do_mc
  Int_t fMcLabels[MaxTrack]; ///< Array to store mc labels inside cluster
#endif
};

/// Parameters which represent given pad in the hough space. Used in order to avoid as much as possible floating point
/// operations during the hough transform
struct AliHLTPadHoughParams {
  Float_t fAlpha;      ///< Starting value for the hough parameter alpha1
  Float_t fDeltaAlpha; ///< Slope of alpha1
  Int_t fFirstBin;     ///< First alpha2 bin to be filled
  Int_t fLastBin;      ///< Last alpha2 bin to be filled
};

class Accumulator;

/// Transforms the TPC data into the hough space and counts the missed TPC rows corresponding to each track candidate -
/// hough space bin
class TransformerRow : public BaseTransformer {

public:
  TransformerRow();
  TransformerRow(Int_t slice, Int_t patch, Int_t netasegments, Bool_t DoMC = false, Float_t zvertex = 0.0);
  virtual ~TransformerRow();

  void CreateHistograms(Float_t ptmin, Float_t ptmax, Float_t pres, Int_t nybin, Float_t psi)
  {
    BaseTransformer::CreateHistograms(ptmin, ptmax, pres, nybin, psi);
  }
  void CreateHistograms(Int_t /*nxbin*/, Float_t /*ptmin*/, Int_t /*nybin*/, Float_t /*phimin*/, Float_t /*phimax*/)
  {
    STDCERR << "This method for creation of parameter space histograms is not supported for this Transformer!"
            << STDENDL;
  }
  void CreateHistograms(Int_t nxbin, Float_t xmin, Float_t xmax, Int_t nybin, Float_t ymin, Float_t ymax);
  void Reset();

  // Contains the hough transformation. It reads as an input the preloaded array with digits
  void TransformCircle(ClusterCollection* clusterCollection);

  void TransformCircle(Int_t* row_range, Int_t every) { BaseTransformer::TransformCircle(row_range, every); }

  /// Returns the histogram index of the corresponding eta.
  Int_t GetEtaIndex(Double_t eta) const;

  /// Returns a pointer to the histogram which contains etaindex eta slice
  Accumulator* GetHistogram(Int_t etaindex);

  /// Returns eta calculated in the middle of the eta slice
  Double_t GetEta(Int_t etaindex, Int_t slice) const;
  Int_t GetTrackID(Int_t etaindex, Double_t alpha1, Double_t alpha2) const;
  Int_t GetTrackLength(Double_t alpha1, Double_t alpha2, Int_t* rows) const;
  UChar_t* GetGapCount(Int_t etaindex) const { return fGapCount[etaindex]; }
  UChar_t* GetCurrentRowCount(Int_t etaindex) const { return fCurrentRowCount[etaindex]; }
  UChar_t* GetPrevBin(Int_t etaindex) const { return fPrevBin[etaindex]; }
  UChar_t* GetNextBin(Int_t etaindex) const { return fNextBin[etaindex]; }
  UChar_t* GetNextRow(Int_t etaindex) const { return fNextRow[etaindex]; }
  UChar_t* GetTrackNRows() const { return fTrackNRows; }
  UChar_t* GetTrackFirstRow() const { return fTrackFirstRow; }
  UChar_t* GetTrackLastRow() const { return fTrackLastRow; }
  static Float_t GetBeta1() { return fgBeta1; }
  static Float_t GetBeta2() { return fgBeta2; }
  static Float_t GetDAlpha() { return fgDAlpha; }
  static Float_t GetDEta() { return fgDEta; }
  static Double_t GetEtaCalcParam1() { return fgEtaCalcParam1; }
  static Double_t GetEtaCalcParam2() { return fgEtaCalcParam2; }
  static Double_t GetEtaCalcParam3() { return fgEtaCalcParam3; }
  /*
    void SetTPCRawStream(AliTPCRawStream* rawstream)
    {
      fTPCRawStream = rawstream;
    }
  */
private:
  UChar_t** fGapCount;        //!
  UChar_t** fCurrentRowCount; //!
#ifdef do_mc
  AliHLTTrackIndex** fTrackID; //!
#endif

  UChar_t* fTrackNRows;      //!
  UChar_t* fTrackFirstRow;   //!
  UChar_t* fTrackLastRow;    //!
  UChar_t* fInitialGapCount; //!

  UChar_t** fPrevBin; //!
  UChar_t** fNextBin; //!
  UChar_t** fNextRow; //!

  AliHLTPadHoughParams** fStartPadParams; //!
  AliHLTPadHoughParams** fEndPadParams;   //!
  Float_t** fLUTr;                        //!

  Float_t* fLUTforwardZ;  //!
  Float_t* fLUTbackwardZ; //!

  Accumulator** fParamSpace; //!

  void TransformCircleFromDigitArray(ClusterCollection* clusterCollection);
  void TransformCircleFromRawStream();

  void DeleteHistograms(); // Method to clean up the histograms containing Hough space

  /// Part of the fast hough transform. It fills one row of the hough space. It is called by the FillCluster() method
  /// inside the loop over alpha2 bins
  inline void FillClusterRow(UChar_t i, Int_t binx1, Int_t binx2, UChar_t* ngaps2, UChar_t* currentrow2,
                             UChar_t* lastrow2
#ifdef do_mc
                             ,
                             AliHLTEtaRow etaclust, AliHLTTrackIndex* trackid
#endif
                             );

  /// Part of the fast hough transform. It fills a TPC cluster into the hough space.
  inline void FillCluster(UChar_t i, Int_t etaindex, AliHLTEtaRow* etaclust, Int_t ilastpatch, Int_t firstbinx,
                          Int_t lastbinx, Int_t nbinx, Int_t firstbiny);
#ifdef do_mc
  /// Part of the fast hough transform. It fills the MC labels of a TPC cluster into a special hough space array.
  inline void FillClusterMCLabels(AliHLTDigitData digpt, AliHLTEtaRow* etaclust);
#endif

  // In case of sequential filling of the hough space, it is used to transmit the pointers to the hough arrays from one
  // transformer to the following one.
  void SetTransformerArrays(TransformerRow* tr);

  static Float_t fgBeta1, fgBeta2; ///< Two curves which define the Hough space
  static Float_t fgDAlpha, fgDEta; ///< Correlation factor between Hough space bin size and resolution
  static Double_t fgEtaCalcParam1,
    fgEtaCalcParam2;               ///< Parameters used for fast calculation of eta during the binning of Hough space
  static Double_t fgEtaCalcParam3; ///< Parameter used during the eta binning of the Hough Space in order to account for
  /// finite track radii

  //  AliTPCRawStream* fTPCRawStream; ///< Pointer to the raw stream in case of fast reading of the raw data (fast_raw
  //  flag)
};
}
}

#endif
