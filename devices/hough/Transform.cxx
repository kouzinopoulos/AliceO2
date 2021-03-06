/// \file Transform.cxx
/// \brief Implementation of the Transform class
/// \author Anders Vestbo <vestbo@fi.uib.no>, Uli Frankenfeld <franken@fi.uib.no>, Constantin Loizides
/// <loizides@ikf.uni-frankfurt.de>

#include "StandardIncludes.h"
#include "Transform.h"

using namespace AliceO2::Hough;

const double Transform::fgkAnodeWireSpacing = 0.25; // Taken from the TDR
const double Transform::fgkBFACT = 0.0029980;       // Conversion Factor
const double Transform::fgkPi = 3.141592653589793;
const double Transform::fgk2Pi = 2 * 3.141592653589793;
const double Transform::fgkPi2 = 0.5 * 3.141592653589793;
const double Transform::fgkToDeg = 180 / 3.141592653589793;

// Defined by HLT and GSI
int Transform::fgNPatches = 6;
int Transform::fgRows[6][2] = { { 0, 29 }, { 30, 62 }, { 63, 90 }, { 91, 116 }, { 117, 139 }, { 140, 158 } };
int Transform::fgNRows[6] = { 30, 33, 28, 26, 23, 19 };

// The following definition is generated by MakeInitFile function
double Transform::fgBField = 0.2;
double Transform::fgSolenoidBField = 2;
double Transform::fgBFieldFactor = 1;
int Transform::fgVersion = kVdefault;
int Transform::fgNTimeBins = 446;
int Transform::fgNRowLow = 63;
int Transform::fgNRowUp = 96;
int Transform::fgNRowUp1 = 64;
int Transform::fgNRowUp2 = 32;
int Transform::fgNSectorLow = 36;
int Transform::fgNSectorUp = 36;
int Transform::fgNSector = 72;
double Transform::fgPadPitchWidthLow = 0.4;
double Transform::fgPadPitchWidthUp = 0.6;
double Transform::fgZWidth = 0.5660;
double Transform::fgZSigma = 0.2288;
double Transform::fgZLength = 250.0000;
double Transform::fgZOffset = 0.6864;
double Transform::fgDiffT = 0.0220;
double Transform::fgDiffL = 0.0220;
double Transform::fgOmegaTau = 0.1450;
double Transform::fgInnerPadLength = 0.75;
double Transform::fgOuter1PadLength = 1.00;
double Transform::fgOuter2PadLength = 1.50;
double Transform::fgInnerPRFSigma = 0.203811;
double Transform::fgOuter1PRFSigma = 0.299325;
double Transform::fgOuter2PRFSigma = 0.299323;
double Transform::fgTimeSigma = 0.228809;
int Transform::fgADCSat = 1024;
int Transform::fgZeroSup = 0;
int Transform::fgNSlice = 36;
int Transform::fgNRow = 159;
double Transform::fgNRotShift = 0.5;
int Transform::fgSlice2Sector[36][2] = { { 0, 36 },
                                         { 1, 37 },
                                         { 2, 38 },
                                         { 3, 39 },
                                         { 4, 40 },
                                         { 5, 41 },
                                         { 6, 42 },
                                         { 7, 43 },
                                         { 8, 44 },
                                         { 9, 45 },
                                         { 10, 46 },
                                         { 11, 47 },
                                         { 12, 48 },
                                         { 13, 49 },
                                         { 14, 50 },
                                         { 15, 51 },
                                         { 16, 52 },
                                         { 17, 53 },
                                         { 18, 54 },
                                         { 19, 55 },
                                         { 20, 56 },
                                         { 21, 57 },
                                         { 22, 58 },
                                         { 23, 59 },
                                         { 24, 60 },
                                         { 25, 61 },
                                         { 26, 62 },
                                         { 27, 63 },
                                         { 28, 64 },
                                         { 29, 65 },
                                         { 30, 66 },
                                         { 31, 67 },
                                         { 32, 68 },
                                         { 33, 69 },
                                         { 34, 70 },
                                         { 35, 71 } };

int Transform::fgSector2Slice[72] = { 0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15, 16, 17,
                                      18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35,
                                      0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15, 16, 17,
                                      18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35 };

int Transform::fgSectorLow[72] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

double Transform::fgX[159] = {
  85.195,  85.945,  86.695,  87.445,  88.195,  88.945,  89.695,  90.445,  91.195,  91.945,  92.695,  93.445,  94.195,
  94.945,  95.695,  96.445,  97.195,  97.945,  98.695,  99.445,  100.195, 100.945, 101.695, 102.445, 103.195, 103.945,
  104.695, 105.445, 106.195, 106.945, 107.695, 108.445, 109.195, 109.945, 110.695, 111.445, 112.195, 112.945, 113.695,
  114.445, 115.195, 115.945, 116.695, 117.445, 118.195, 118.945, 119.695, 120.445, 121.195, 121.945, 122.695, 123.445,
  124.195, 124.945, 125.695, 126.445, 127.195, 127.945, 128.695, 129.445, 130.195, 130.945, 131.695, 135.180, 136.180,
  137.180, 138.180, 139.180, 140.180, 141.180, 142.180, 143.180, 144.180, 145.180, 146.180, 147.180, 148.180, 149.180,
  150.180, 151.180, 152.180, 153.180, 154.180, 155.180, 156.180, 157.180, 158.180, 159.180, 160.180, 161.180, 162.180,
  163.180, 164.180, 165.180, 166.180, 167.180, 168.180, 169.180, 170.180, 171.180, 172.180, 173.180, 174.180, 175.180,
  176.180, 177.180, 178.180, 179.180, 180.180, 181.180, 182.180, 183.180, 184.180, 185.180, 186.180, 187.180, 188.180,
  189.180, 190.180, 191.180, 192.180, 193.180, 194.180, 195.180, 196.180, 197.180, 198.180, 199.430, 200.930, 202.430,
  203.930, 205.430, 206.930, 208.430, 209.930, 211.430, 212.930, 214.430, 215.930, 217.430, 218.930, 220.430, 221.930,
  223.430, 224.930, 226.430, 227.930, 229.430, 230.930, 232.430, 233.930, 235.430, 236.930, 238.430, 239.930, 241.430,
  242.930, 244.430, 245.930
};

int Transform::fgNPads[159] = {
  67,  67,  69,  69,  69,  71,  71,  71,  73,  73,  73,  75,  75,  75,  77,  77,  77,  79,  79,  79,  81,  81,  81,
  83,  83,  83,  85,  85,  85,  87,  87,  87,  89,  89,  89,  91,  91,  91,  93,  93,  93,  95,  95,  95,  97,  97,
  97,  99,  99,  99,  99,  101, 101, 101, 103, 103, 103, 105, 105, 105, 107, 107, 107, 73,  75,  75,  75,  75,  77,
  77,  77,  79,  79,  79,  81,  81,  81,  81,  83,  83,  83,  85,  85,  85,  85,  87,  87,  87,  89,  89,  89,  91,
  91,  91,  91,  93,  93,  93,  95,  95,  95,  95,  97,  97,  97,  99,  99,  99,  101, 101, 101, 101, 103, 103, 103,
  105, 105, 105, 105, 107, 107, 107, 109, 109, 109, 111, 111, 111, 113, 113, 113, 115, 115, 117, 117, 119, 119, 121,
  121, 121, 123, 123, 125, 125, 127, 127, 127, 129, 129, 131, 131, 133, 133, 135, 135, 135, 137, 137, 139
};

double Transform::fgCos[36] = { 0.9848077297,  0.8660253882,  0.6427876353,  0.3420201540,  0.0000000000,
                                -0.3420201540, -0.6427876353, -0.8660253882, -0.9848077297, -0.9848077297,
                                -0.8660253882, -0.6427876353, -0.3420201540, -0.0000000000, 0.3420201540,
                                0.6427876353,  0.8660253882,  0.9848077297,  0.9848077297,  0.8660253882,
                                0.6427876353,  0.3420201540,  0.0000000000,  -0.3420201540, -0.6427876353,
                                -0.8660253882, -0.9848077297, -0.9848077297, -0.8660253882, -0.6427876353,
                                -0.3420201540, -0.0000000000, 0.3420201540,  0.6427876353,  0.8660253882,
                                0.9848077297 };

double Transform::fgSin[36] = { 0.1736481786,  0.5000000000,  0.7660444379,  0.9396926165,  1.0000000000,
                                0.9396926165,  0.7660444379,  0.5000000000,  0.1736481786,  -0.1736481786,
                                -0.5000000000, -0.7660444379, -0.9396926165, -1.0000000000, -0.9396926165,
                                -0.7660444379, -0.5000000000, -0.1736481786, 0.1736481786,  0.5000000000,
                                0.7660444379,  0.9396926165,  1.0000000000,  0.9396926165,  0.7660444379,
                                0.5000000000,  0.1736481786,  -0.1736481786, -0.5000000000, -0.7660444379,
                                -0.9396926165, -1.0000000000, -0.9396926165, -0.7660444379, -0.5000000000,
                                -0.1736481786 };

int Transform::GetNPads(int row)
{
  if (row < 0 || row >= fgNRow) {
    cout << "Wrong row " << row << endl;
    return 0;
  }
  return fgNPads[row];
}

int Transform::GetFirstRow(int patch)
{
  if (patch == -1) {
    return 0;
  } else if (patch < -1 || patch >= 6) {
    cout << "Wrong patch " << patch << endl;
    return 0;
  } else {
    return fgRows[patch][0];
  }
}

int Transform::GetLastRow(int patch)
{
  if (patch == -1) {
    return fgRows[5][1];
  } else if (patch < -1 || patch >= 6) {
    cout << "Wrong patch " << patch << endl;
    return 0;
  } else {
    return fgRows[patch][1];
  }
}

int Transform::GetFirstRowOnDDL(int patch)
{
  if (patch == -1) {
    return 0;
  } else if (patch < -1 || patch >= 6) {
    cout << "Wrong patch " << patch << endl;
    return 0;
  } else {
    if (patch == 1) {
      return fgRows[patch][0] + 1;
    }
    return fgRows[patch][0];
  }
}

int Transform::GetLastRowOnDDL(int patch)
{
  if (patch == -1) {
    return fgRows[5][1];
  } else if (patch < -1 || patch >= 6) {
    cerr << "Wrong patch " << patch << endl;
    return 0;
  } else {
    if (patch == 2 || patch == 4) {
      return fgRows[patch][1] - 1;
    }
    return fgRows[patch][1];
  }
}

int Transform::GetNRows(int patch)
{
  if (patch == -1) {
    return fgNRow;
  } else if (patch < -1 || patch >= 6) {
    cerr << "Wrong patch " << patch << endl;
    return 0;
  } else {
    return fgNRows[patch];
  }
}

int Transform::GetPadRow(float xvalue)
{
  if (xvalue < 0 || xvalue > 250) {
    cerr << "Suspicious x-value, make sure it is in local coordinate! " << xvalue << endl;
    return -1;
  }

  int x = (int)rint(xvalue * 10);
  if (x < (int)rint(fgX[1] * 10)) {
    return 0;
  } else if (x > (int)rint(fgX[fgNRow - 2] * 10)) {
    return fgNRow - 1;
  } else {
    int padrow = 1; // Of course, a more clever algorithm could help here
    while (padrow < fgNRow - 2) {
      if (x > (int)rint(fgX[padrow - 1] * 10) && x < (int)rint(fgX[padrow + 1] * 10)) {
        break;
      }
      padrow++;
    }
    return padrow;
  }
}

int Transform::GetPatch(int padrow)
{
  if (padrow < 0 || padrow >= fgNRow) {
    cerr << "Wrong padrow " << padrow << endl;
    return -2;
  }
  int patch = 0;
  while (patch < fgNPatches) {
    if (padrow >= fgRows[patch][0] && padrow <= fgRows[patch][1]) {
      break;
    }
    patch++;
  }
  return patch;
}

double Transform::GetPadLength(int padrow)
{
  if (padrow >= fgNRow) {
    cerr << "Wrong padrow " << padrow << endl;
    return 0;
  }

  if (padrow < fgNRowLow) {
    return fgInnerPadLength;
  }
  if (padrow >= fgNRowLow && padrow < fgNRowLow + fgNRowUp1 - 1) {
    return fgOuter1PadLength;
  }
  if (padrow >= fgNRowLow + fgNRowUp1 - 1) {
    return fgOuter2PadLength;
  }

  // should never happen
  cerr << "Wrong padrow " << padrow << endl;
  return -1.0;
}

double Transform::GetPadPitchWidth(int patch)
{
  if (patch < 0 || patch > fgNPatches) {
    cerr << "Wrong patch " << patch << endl;
    return -1;
  }
  return patch < 2 ? fgPadPitchWidthLow : fgPadPitchWidthUp;
}

double Transform::GetParSigmaY2(int padrow, float z, float angle)
{
  double drift;
  if (z > 0) {
    drift = fgZLength - z;
  } else {
    drift = fgZLength + z;
  }

  double t1 = GetPRFSigma(padrow) * GetPRFSigma(padrow);
  double t2 = fgDiffT * fgDiffT * drift;
  double t3 = GetPadLength(padrow) * GetPadLength(padrow) * tan(angle) * tan(angle) / 12;
  double t4 = fgkAnodeWireSpacing * fgkAnodeWireSpacing * (tan(angle) - fgOmegaTau) * (tan(angle) - fgOmegaTau) / 12;

  return (t1 + t2 + t3 + t4);
}

double Transform::GetParSigmaZ2(int padrow, float z, float tgl)
{
  double drift;
  if (z > 0) {
    drift = Transform::GetZLength() - 0.275 - z;
  } else {
    drift = Transform::GetZLength() - 0.302 + z;
  }

  double t1 = fgZSigma * fgZSigma;
  double t2 = fgDiffL * fgDiffL * drift;
  double t3 = GetPadLength(padrow) * GetPadLength(padrow) * tgl * tgl / 12;

  return (t1 + t2 + t3);
}

double Transform::GetPRFSigma(int padrow)
{
  if (padrow >= fgNRow) {
    cerr << "Wrong padrow " << padrow << endl;
    return 0;
  }
  if (padrow < fgNRowLow) {
    return fgInnerPRFSigma;
  }
  if (padrow >= fgNRowLow && padrow < fgNRowLow + fgNRowUp1 - 1) {
    return fgOuter1PRFSigma;
  }
  if (padrow >= fgNRowLow + fgNRowUp1 - 1) {
    return fgOuter2PRFSigma;
  }

  // should never happen
  cerr << "Wrong padrow " << padrow << endl;
  return -1.;
}

double Transform::GetEta(float* xyz)
{
  double r3 = sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2]);
  double eta = 0.5 * log((r3 + xyz[2]) / (r3 - xyz[2]));
  return eta;
}

void Transform::XYZtoRPhiEta(float* rpe, float* xyz)
{
  rpe[0] = sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2]);
  rpe[1] = atan2(xyz[1], xyz[0]);
  rpe[2] = 0.5 * log((rpe[0] + xyz[2]) / (rpe[0] - xyz[2]));
}

double Transform::GetEta(int slice, int padrow, int pad, int time)
{
  float xyz[3];
  int sector, row;
  Slice2Sector(slice, padrow, sector, row);
  Raw2Local(xyz, sector, row, pad, time);

  return GetEta(xyz);
}

double Transform::GetPhi(float* xyz)
{
  double phi = atan2(xyz[1], xyz[0]);
  return phi;
}

bool Transform::Slice2Sector(int slice, int slicerow, int& sector, int& row)
{
  if (slicerow < 0 && slicerow >= fgNRow) {
    cout << "Wrong slicerow " << slicerow << endl;
    return false;
  }
  if (slice < 0 || slice >= fgNSlice) {
    cout << "Wrong slice " << slice << endl;
    return false;
  }

  if (slicerow < fgNRowLow) {
    sector = fgSlice2Sector[slice][0];
    row = slicerow;
  } else {
    sector = fgSlice2Sector[slice][1];
    row = slicerow - fgNRowLow;
  }

  return true;
}

bool Transform::Sector2Slice(int& slice, int sector)
{
  if (sector < 0 || sector >= fgNSector) {
    cout << "Wrong sector " << sector << endl;
    return false;
  }

  slice = fgSector2Slice[sector];

  return true;
}

bool Transform::Sector2Slice(int& slice, int& slicerow, int sector, int row)
{
  if (sector < 0 || sector >= fgNSector) {
    cout << "Wrong sector " << sector << endl;
    return false;
  }
  if (row < 0) {
    cout << "Wrong row " << row << endl;
    return false;
  }

  if (fgSectorLow[sector]) {
    if (row >= fgNRowLow) {
      cout << "Wrong row " << row << endl;
      return false;
    }
    slice = fgSector2Slice[sector];
    slicerow = row;
  } else {
    if (row >= fgNRowUp) {
      cout << "Wrong row " << row << endl;
      return false;
    }
    slice = fgSector2Slice[sector];
    slicerow = row + fgNRowLow;
  }

  return true;
}

double Transform::GetMaxY(int slicerow)
{
  if (slicerow < fgNRowLow) {
    return fgPadPitchWidthLow * fgNPads[slicerow] / 2;
  } else {
    return fgPadPitchWidthUp * fgNPads[slicerow] / 2;
  }
}

double Transform::Row2X(int slicerow)
{
  if (slicerow < 0 || slicerow >= fgNRow) {
    cout << "Wrong slicerow " << slicerow << endl;
    return 0;
  }
  return fgX[slicerow];
}

double Transform::GetZFast(int slice, int time, float vertex)
{
  double z = fgZWidth * time - fgZOffset;
  if (slice < 18) {
    z = fgZLength - z - vertex;
  } else {
    z = z - fgZLength - vertex;
  }
  return z;
}

void Transform::Local2Global(float* xyz, int slice)
{
  float x0 = xyz[0];
  float y0 = xyz[1];

  xyz[0] = x0 * fgCos[slice] - y0 * fgSin[slice];
  xyz[1] = x0 * fgSin[slice] + y0 * fgCos[slice];
  xyz[2] = xyz[2]; // global z=local z
}

void Transform::Local2GlobalAngle(float* angle, int slice)
{
  angle[0] = fmod(angle[0] + (slice + fgNRotShift) * (2 * fgkPi / 18), 2 * fgkPi);
}

void Transform::Global2LocalAngle(float* angle, int slice)
{
  // get angle local
  angle[0] = angle[0] - (slice + fgNRotShift) * (2 * fgkPi / 18);
  if (angle[0] < 0)
    angle[0] += 2 * fgkPi;
}

void Transform::Raw2Local(float* xyz, int sector, int row, float pad, float time)
{
  // Transformation from rawdata to local coordinate system

  int slice, slicerow;
  Sector2Slice(slice, slicerow, sector, row);

  // X-Value
  xyz[0] = Row2X(slicerow);

  // Y-Value
  int npads = fgNPads[slicerow];

  if (fgSectorLow[sector])
    xyz[1] = (pad - 0.5 * (npads - 1)) * fgPadPitchWidthLow;
  else
    xyz[1] = (pad - 0.5 * (npads - 1)) * fgPadPitchWidthUp;

  // Z-Value (remember PULSA Delay)
  if (slice < 18)
    xyz[2] = fgZLength - fgZWidth * time + fgZOffset;
  else
    xyz[2] = fgZWidth * time - fgZOffset - fgZLength;
}

void Transform::Raw2Local(float* xyz, int sector, int row, int pad, int time)
{
  // Transformation from rawdata to local coordinate system

  int slice, slicerow;
  Sector2Slice(slice, slicerow, sector, row);

  // X-Value
  xyz[0] = Row2X(slicerow);

  // Y-Value
  int npads = fgNPads[slicerow];

  if (fgSectorLow[sector])
    xyz[1] = (pad - 0.5 * (npads - 1)) * fgPadPitchWidthLow;
  else
    xyz[1] = (pad - 0.5 * (npads - 1)) * fgPadPitchWidthUp;

  // Z-Value (remember PULSA Delay)
  if (slice < 18)
    xyz[2] = fgZLength - fgZWidth * time + fgZOffset;
  else
    xyz[2] = fgZWidth * time - fgZOffset - fgZLength;
}

void Transform::RawHLT2Local(float* xyz, int slice, int slicerow, float pad, float time)
{
  // Transformation from HLT rawdata to local coordinate system

  // X-Value
  xyz[0] = Row2X(slicerow);

  // Y-Value
  int npads = fgNPads[slicerow];
  if (slicerow < fgNRowLow)
    xyz[1] = (pad - 0.5 * (npads - 1)) * fgPadPitchWidthLow;
  else
    xyz[1] = (pad - 0.5 * (npads - 1)) * fgPadPitchWidthUp;

  // Z-Value
  if (slice < 18)
    xyz[2] = fgZLength - fgZWidth * time + fgZOffset;
  else
    xyz[2] = fgZWidth * time - fgZOffset - fgZLength;
}

void Transform::RawHLT2Local(float* xyz, int slice, int slicerow, int pad, int time)
{
  // Transformation from HLT rawdata to local coordinate system

  // X-Value
  xyz[0] = Row2X(slicerow);

  // Y-Value
  int npads = fgNPads[slicerow];
  if (slicerow < fgNRowLow)
    xyz[1] = (pad - 0.5 * (npads - 1)) * fgPadPitchWidthLow;
  else
    xyz[1] = (pad - 0.5 * (npads - 1)) * fgPadPitchWidthUp;

  // Z-Value
  if (slice < 18)
    xyz[2] = fgZLength - fgZWidth * time + fgZOffset;
  else
    xyz[2] = fgZWidth * time - fgZOffset - fgZLength;
}

void Transform::Local2Global(float* xyz, int sector, int row)
{
  // Transformation to global coordinate system
  int slice, slicerow;
  Sector2Slice(slice, slicerow, sector, row);
  float r = Row2X(slicerow); // have to get x value first

  xyz[0] = r * fgCos[slice] - xyz[1] * fgSin[slice];
  xyz[1] = r * fgSin[slice] + xyz[1] * fgCos[slice];
  xyz[2] = xyz[2]; // global z=local z
}

void Transform::LocHLT2Global(float* xyz, int slice, int slicerow)
{
  // Transformation from HLT to global coordinate system
  float r = Row2X(slicerow); // have to get x value first

  xyz[0] = r * fgCos[slice] - xyz[1] * fgSin[slice];
  xyz[1] = r * fgSin[slice] + xyz[1] * fgCos[slice];
  xyz[2] = xyz[2]; // global z=local z
}

void Transform::Global2Local(float* xyz, int sector)
{ // check code
  int slice;
  Sector2Slice(slice, sector);

  float x1 = xyz[0] * fgCos[slice] + xyz[1] * fgSin[slice];
  float y1 = -xyz[0] * fgSin[slice] + xyz[1] * fgCos[slice];
  xyz[0] = x1;
  xyz[1] = y1;
}

void Transform::Global2LocHLT(float* xyz, int slice)
{
  float x1 = xyz[0] * fgCos[slice] + xyz[1] * fgSin[slice];
  float y1 = -xyz[0] * fgSin[slice] + xyz[1] * fgCos[slice];
  xyz[0] = x1;
  xyz[1] = y1;
}

void Transform::Raw2Global(float* xyz, int sector, int row, float pad, float time)
{
  // Transformation from raw to global coordinates

  Raw2Local(xyz, sector, row, pad, time);
  Local2Global(xyz, sector, row);
}

void Transform::Raw2Global(float* xyz, int sector, int row, int pad, int time)
{
  // Transformation from raw to global coordinates

  Raw2Local(xyz, sector, row, pad, time);
  Local2Global(xyz, sector, row);
}

void Transform::RawHLT2Global(float* xyz, int slice, int slicerow, float pad, float time)
{
  // Transformation from raw to global coordinates

  RawHLT2Local(xyz, slice, slicerow, pad, time);
  LocHLT2Global(xyz, slice, slicerow);
}

void Transform::RawHLT2Global(float* xyz, int slice, int slicerow, int pad, int time)
{
  // Transformation from raw to global coordinates

  RawHLT2Local(xyz, slice, slicerow, pad, time);
  LocHLT2Global(xyz, slice, slicerow);
}

void Transform::Local2Raw(float* xyz, int sector, int row)
{
  // Transformation from local coordinates to raw

  int slice, slicerow;
  Sector2Slice(slice, slicerow, sector, row);

  xyz[0] = slicerow;

  if (fgSectorLow[sector])
    xyz[1] = xyz[1] / fgPadPitchWidthLow + 0.5 * (fgNPads[slicerow] - 1);
  else
    xyz[1] = xyz[1] / fgPadPitchWidthUp + 0.5 * (fgNPads[slicerow] - 1);

  if (slice < 18)
    xyz[2] = (fgZLength - xyz[2] + fgZOffset) / fgZWidth;
  else
    xyz[2] = (fgZLength + xyz[2] + fgZOffset) / fgZWidth;
}

void Transform::LocHLT2Raw(float* xyz, int slice, int slicerow)
{
  // Transformation from local coordinates to raw

  xyz[0] = slicerow;

  if (slicerow < fgNRowLow)
    xyz[1] = xyz[1] / fgPadPitchWidthLow + 0.5 * (fgNPads[slicerow] - 1);
  else
    xyz[1] = xyz[1] / fgPadPitchWidthUp + 0.5 * (fgNPads[slicerow] - 1);

  if (slice < 18)
    xyz[2] = (fgZLength - xyz[2] + fgZOffset) / fgZWidth;
  else
    xyz[2] = (fgZLength + xyz[2] + fgZOffset) / fgZWidth;
}

void Transform::Global2Raw(float* xyz, int sector, int row)
{
  // Transformation from global coordinates to raw.

  Global2Local(xyz, sector);
  Local2Raw(xyz, sector, row);
}

void Transform::Global2HLT(float* xyz, int slice, int slicerow)
{
  // Transformation from global coordinates to raw.

  Global2LocHLT(xyz, slice);
  LocHLT2Raw(xyz, slice, slicerow);
}
