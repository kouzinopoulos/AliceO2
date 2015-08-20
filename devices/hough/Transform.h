/// \file Transform.h
/// \brief Definition of the Transform class
/// \author Anders Vestbo <vestbo@fi.uib.no>, Uli Frankenfeld <franken@fi.uib.no>, Constantin Loizides
/// <loizides@ikf.uni-frankfurt.de>

#ifndef ALICEO2_HOUGH_TRANSFORM_H_
#define ALICEO2_HOUGH_TRANSFORM_H_

namespace AliceO2 {
namespace Hough {

/// Transformation class for ALICE TPC.
/// Class which contains all detector specific parameters for the TPC and different useful functions for coordinate
/// transforms. The class is completely static, which means that no object needs to be instantiated. Function calls
/// should then be done like, e.g.: double eta = Transform::GetEta(xyz);
/// IMPORTANT: If used as is, default detector parameters will be used and you really have to make sure that these
/// correspond to the AliROOT version you are currently working on. You should therefore always initialize the
/// parameters by
///
/// Transform::Init(path);
///
/// where path is a char*, giving the path to where file containing the detector parameter is located. This file should
/// be called "l3transform.config", and can be created with the function MakeInitFile.
/// You can also force reading the parameters from a AliTPCParam object by setting the flag;
///
/// Transform::Init(path,true);
///
/// where path is a char* either providing the rootfile name containing the geometry or the path to the rootfile which
/// should then be called alirunfile.root. Note that for both of these cases you have to compile with USEPACKAGE=ALIROOT
/// set (see level3code/Makefile.conf).
/// Currently, there are 4 versions of the Transformer:
///             fVersion==kValiroot: ALIROOT-head compatible
///             fVersion==kVcosmics: Cosmics data run (2003) compatible
///             fVersion==kVdefault: means no config file has been loaded
///             fVersion==kVdeprecated: dont use old (before July 2003) style of transformer

class Transform {

public:
  enum VersionType { kVdefault = 0, kVdeprecated = 1, kValiroot = 10, kVcosmics = 100 };

private:
  static const double fgkBFACT;            ///< bfield
  static const double fgkPi;               ///< pi
  static const double fgkPi2;              ///< 2pi
  static const double fgk2Pi;              ///< pi/2
  static const double fgkAnodeWireSpacing; ///< anode wire spacing
  static const double fgkToDeg;            ///< rad to deg

  static int fgNPatches;   ///< 6 (dont change this)
  static int fgRows[6][2]; ///< rows per patch
  static int fgNRows[6];   ///< rows per patch

  static double fgBField;           ///< field
  static double fgBFieldFactor;     ///< field
  static double fgSolenoidBField;   ///< field
  static int fgNTimeBins;           ///< ntimebins
  static int fgNRowLow;             ///< nrows
  static int fgNRowUp;              ///< nrows
  static int fgNRowUp1;             ///< nrows
  static int fgNRowUp2;             ///< nrows
  static int fgNSectorLow;          ///< nsector
  static int fgNSectorUp;           ///< nsector
  static int fgSlice2Sector[36][2]; ///< nslice
  static int fgSector2Slice[72];    ///< nslice
  static int fgSectorLow[72];       ///< nsector
  static double fgPadPitchWidthLow; ///< pad pitch
  static double fgPadPitchWidthUp;  ///< pad pitch
  static double fgZWidth;           ///< width
  static double fgZSigma;           ///< sigma
  static double fgZLength;          ///< length
  static double fgZOffset;          ///< offset
  static int fgNSector;             ///< 72  (dont change this)
  static int fgNSlice;              ///< 36  (dont change this)
  static int fgNRow;                ///< 159 (dont change this)
  static double fgNRotShift;        ///< Rotation shift (eg. 0.5 for 10 degrees)
  static int fgNPads[159];          ///< fill this following Init and fVersion
  static double fgX[159];           ///< X position in local coordinates
  static int fgVersion;             ///< flags the version
  static double fgDiffT;            ///< Transversal diffusion constant
  static double fgDiffL;            ///< Longitudinal diffusion constant
  static double fgOmegaTau;         ///< ExB effects
  static double fgInnerPadLength;   ///< innner pad length
  static double fgOuter1PadLength;  ///< outer pad length
  static double fgOuter2PadLength;  ///< outer pad length
  static double fgInnerPRFSigma;    ///< inner pad response function
  static double fgOuter1PRFSigma;   ///< outer pad response function
  static double fgOuter2PRFSigma;   ///< outer pad response function
  static double fgTimeSigma;        ///< Minimal longitudinal width
  static int fgADCSat;              ///< ADC Saturation (1024 = 10 bit)
  static int fgZeroSup;             ///< Zero suppression threshold
  static double fgCos[36];          ///< stores the cos value for local to global rotations
  static double fgSin[36];          ///< stores the sin value for local to global rotations

public:
  virtual ~Transform() {}

  // setters
  static void SetNPatches(int i) { fgNPatches = i; }
  static void SetNRows(int s[6])
  {
    for (int i = 0; i < fgNPatches; i++)
      fgNRows[i] = s[i];
  }
  static void SetRows(int s[6][2])
  {
    for (int i = 0; i < fgNPatches; i++) {
      fgRows[i][0] = s[i][0];
      fgRows[i][1] = s[i][1];
    }
  }
  static void SetBField(double f) { fgBField = f; } // careful, these 3 are not independent!
  static void SetBFieldFactor(double f)
  {
    fgBFieldFactor = f;
    fgBField = fgBFieldFactor * fgSolenoidBField * 0.1;
  }
  static void SetSolenoidBField(double f)
  {
    fgSolenoidBField = f;
    fgBField = fgBFieldFactor * fgSolenoidBField * 0.1;
  }
  static void SetNTimeBins(int i) { fgNTimeBins = i; }
  static void SetNRowLow(int i) { fgNRowLow = i; }
  static void SetNRowUp(int i) { fgNRowUp = i; }
  static void SetNRowUp1(int i) { fgNRowUp1 = i; }
  static void SetNRowUp2(int i) { fgNRowUp2 = i; }
  static void SetSlice2Sector(int s[36][2])
  {
    for (int i = 0; i < fgNSlice; i++) {
      fgSlice2Sector[i][0] = s[i][0];
      fgSlice2Sector[i][1] = s[i][1];
    }
  }
  static void SetSector2Slice(int s[72])
  {
    for (int i = 0; i < fgNSector; i++)
      fgSector2Slice[i] = s[i];
  }
  static void SetSectorLow(int s[72])
  {
    for (int i = 0; i < fgNSector; i++)
      fgSectorLow[i] = s[i];
  }
  static void SetNSectorLow(int i) { fgNSectorLow = i; }
  static void SetNSectorUp(int i) { fgNSectorUp = i; }
  static void SetPadPitchWidthLow(double f) { fgPadPitchWidthLow = f; }
  static void SetPadPitchWidthUp(double f) { fgPadPitchWidthUp = f; }
  static void SetZWidth(double f) { fgZWidth = f; }
  static void SetZSigma(double f) { fgZSigma = f; }
  static void SetZLength(double f) { fgZLength = f; }
  static void SetZOffset(double f) { fgZOffset = f; }
  static void SetNSector(int i) { fgNSector = i; }
  static void SetNSlice(int i) { fgNSlice = i; }
  static void SetNRow(int i) { fgNRow = i; }
  static void SetNRotShift(double f) { fgNRotShift = f; }
  static void SetNPads(int pads[159])
  {
    for (int i = 0; i < fgNRow; i++)
      fgNPads[i] = pads[i];
  }
  static void SetX(double xs[159])
  {
    for (int i = 0; i < fgNRow; i++)
      fgX[i] = xs[i];
  }
  static void SetVersion(int i) { fgVersion = i; }
  static void SetDiffT(double f) { fgDiffT = f; }
  static void SetDiffL(double f) { fgDiffL = f; }
  static void SetOmegaTau(double f) { fgOmegaTau = f; }
  static void SetInnerPadLength(double f) { fgInnerPadLength = f; }
  static void SetOuter1PadLength(double f) { fgOuter1PadLength = f; }
  static void SetOuter2PadLength(double f) { fgOuter2PadLength = f; }
  static void SetInnerPRFSigma(double f) { fgInnerPRFSigma = f; }
  static void SetOuter1PRFSigma(double f) { fgOuter1PRFSigma = f; }
  static void SetOuter2PRFSigma(double f) { fgOuter2PRFSigma = f; }
  static void SetTimeSigma(double f) { fgTimeSigma = f; }
  static void SetADCSat(int i) { fgADCSat = i; }
  static void SetZeroSup(int i) { fgZeroSup = i; }

  // getters
  static const char* GetParamName() { return "75x40_100x60_150x60"; }
  static double Pi() { return fgkPi; }
  static double PiHalf() { return fgkPi2; }
  static double TwoPi() { return fgk2Pi; }
  static double GetAnodeWireSpacing() { return fgkAnodeWireSpacing; }
  static double GetBFact() { return fgkBFACT; }
  static double ToRad() { return 1. / fgkToDeg; }
  static double ToDeg() { return fgkToDeg; }

  /// Returns the first row per patch
  static int GetFirstRow(int patch);

  /// Returns the last row per patch
  static int GetLastRow(int patch);

  /// Returns the first row per patch
  static int GetFirstRowOnDDL(int patch);

  /// Returns the last row per patch
  static int GetLastRowOnDDL(int patch);
  static int GetNRows(int patch);

  /// Returns patch for the given padrow
  static int GetPatch(int padrow);

  /// Returns the number of rows per patch
  static int GetNRows() { return fgNRow; }
  static int GetNRowLow() { return fgNRowLow; }
  static int GetNRowUp1() { return fgNRowUp1; }
  static int GetNRowUp2() { return fgNRowUp2; }

  /// Returns the padrow number corresponding to cartesian _local_ x value
  static int GetPadRow(float x);

  static int GetNPatches() { return fgNPatches; }

  /// Returns the number of pads per row
  static int GetNPads(int row);
  static int GetNTimeBins() { return fgNTimeBins; }
  static double GetBField() { return fgBField; }
  static double GetSolenoidField() { return fgSolenoidBField; }
  static double GetBFactFactor() { return fgBFieldFactor; }
  static double GetBFieldValue() { return (fgBField * fgkBFACT); }
  static float Deg2Rad(float angle) { return angle / fgkToDeg; }
  static float Rad2Deg(float angle) { return angle * fgkToDeg; }
  static int GetVersion() { return fgVersion; }
  static double GetPadPitchWidthLow() { return fgPadPitchWidthLow; }
  static double GetPadPitchWidthUp() { return fgPadPitchWidthUp; }

  /// Returns the pad patch width for the given patch
  static double GetPadPitchWidth(int patch);
  static double GetZWidth() { return fgZWidth; }
  static double GetZLength() { return fgZLength; }
  static double GetZOffset() { return fgZOffset; }
  static double GetDiffT() { return fgDiffT; }
  static double GetDiffL() { return fgDiffL; }

  /// Calculates the expected transverse cluster width as a function of drift distance and crossing angle.
  /// \param z local z-coordinate of cluster
  /// \param angle track crossing angle with normal to padrow plane
  /// return value = sigma^2 (cartesian coordinates)
  static double GetParSigmaY2(int padrow, float z, float angle);

  /// Calculates the expected longitudinal cluster width as a function of drift distance and track crossing angle.
  /// \param z local z-coordinate of cluster
  /// \param tgl tan(dipangle)
  /// return value = sigma^2 (cartesian coordinates)
  static double GetParSigmaZ2(int padrow, float z, float tgl);
  static double GetOmegaTau() { return fgOmegaTau; }

  /// Returns the pad length for the given padrow
  static double GetPadLength(int padrow);

  /// Returns sigma of pad response function for the given padrow
  static double GetPRFSigma(int padrow);
  static double GetTimeSigma() { return fgTimeSigma; }
  static double GetZSigma() { return fgZSigma; }
  static int GetADCSat() { return fgADCSat; }
  static int GetZeroSup() { return fgZeroSup; }
  static int GetNSlice() { return fgNSlice; }
  static int GetNSector() { return fgNSector; }
  static int GetNSectorLow() { return fgNSectorLow; }
  static int GetNSectorUp() { return fgNSectorUp; }

  // Slice to sector number
  static bool Slice2Sector(int slice, int slicerow, int& sector, int& row);

  /// Sector number to slice
  static bool Sector2Slice(int& slice, int sector);

  /// Sector number to slice
  static bool Sector2Slice(int& slice, int& slicerow, int sector, int row);

  /// Slicerow to X value (slice 0)
  static double Row2X(int slicerow);

  /// Returns the maximum y value (for slice 0)
  static double GetMaxY(int slicerow);

  /// Returns Eta
  static double GetEta(float* xyz);

  /// Returns Eta
  static double GetEta(int slice, int padrow, int pad, int time);

  /// Returns phi
  static double GetPhi(float* xyz);

  /// Returns z
  static double GetZFast(int slice, int time, float vertex = 0.);

  /// Transforms xyz into rpe
  static void XYZtoRPhiEta(float* rpe, float* xyz);

  /// Transformation to global coordinate system
  static void Local2Global(float* xyz, int slice);

  /// Returns angle global
  static void Local2GlobalAngle(float* angle, int slice);
  static void Global2LocalAngle(float* angle, int slice);

  static void Raw2Local(float* xyz, int sector, int row, float pad, float time);
  static void RawHLT2Local(float* xyz, int slice, int slicerow, float pad, float time);
  static void Raw2Local(float* xyz, int sector, int row, int pad, int time);
  static void RawHLT2Local(float* xyz, int slice, int slicerow, int pad, int time);
  static void Local2Global(float* xyz, int sector, int row);
  static void LocHLT2Global(float* xyz, int slice, int slicerow);
  static void Global2Local(float* xyz, int sector);
  static void Global2LocHLT(float* xyz, int slice);
  static void Raw2Global(float* xyz, int sector, int row, float pad, float time);
  static void RawHLT2Global(float* xyz, int slice, int slicerow, float pad, float time);
  static void Raw2Global(float* xyz, int sector, int row, int pad, int time);
  static void RawHLT2Global(float* xyz, int slice, int slicerow, int pad, int time);
  static void Local2Raw(float* xyz, int sector, int row);
  static void LocHLT2Raw(float* xyz, int slice, int slicerow);
  static void Global2Raw(float* xyz, int sector, int row);
  static void Global2HLT(float* xyz, int slice, int slicerow);
};
}
}

#endif
