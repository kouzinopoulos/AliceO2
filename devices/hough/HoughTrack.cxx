// @(#) $Id: HoughTrack.cxx 15995 2006-11-30 17:45:45Z hristov $

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "HoughTrack.h"

using namespace AliceO2::Hough;

HoughTrack::HoughTrack() : Track()
{
  fWeight = 0;
  fMinDist = 0;
  fDLine = 0;
  fPsiLine = 0;
  fIsHelix = true;
  fEtaIndex = -1;
  fEta = 0;
  ComesFromMainVertex(true);
}

HoughTrack::~HoughTrack()
{
}

void HoughTrack::Set(Track* track)
{
  // Basically copy constructor
  HoughTrack* tpt = (HoughTrack*)track;
  SetTrackParameters(tpt->GetKappa(), tpt->GetPsi(), tpt->GetWeight());
  SetEtaIndex(tpt->GetEtaIndex());
  SetEta(tpt->GetEta());
  SetTgl(tpt->GetTgl());
  SetPsi(tpt->GetPsi());
  SetPterr(tpt->GetPterr());
  SetTglerr(tpt->GetTglerr());
  SetPsierr(tpt->GetPsierr());
  SetCenterX(tpt->GetCenterX());
  SetCenterY(tpt->GetCenterY());
  SetFirstPoint(tpt->GetFirstPointX(), tpt->GetFirstPointY(), tpt->GetFirstPointZ());
  SetLastPoint(tpt->GetLastPointX(), tpt->GetLastPointY(), tpt->GetLastPointZ());
  SetCharge(tpt->GetCharge());
  SetRowRange(tpt->GetFirstRow(), tpt->GetLastRow());
  SetSlice(tpt->GetSlice());
  SetHits(tpt->GetNHits(), (unsigned int*)tpt->GetHitNumbers());
  SetMCid(tpt->GetMCid());
  SetBinXY(tpt->GetBinX(), tpt->GetBinY(), tpt->GetSizeX(), tpt->GetSizeY());
  SetSector(tpt->GetSector());
  return;

  //    fWeight = tpt->GetWeight();
  //    fDLine = tpt->GetDLine();
  //    fPsiLine = tpt->GetPsiLine();
  //    SetNHits(tpt->GetWeight());
  //    SetRowRange(tpt->GetFirstRow(),tpt->GetLastRow());
  //    fIsHelix = false;
}

int HoughTrack::Compare(const Track* tpt) const
{
  // Compare 2 hough tracks according to their weight
  HoughTrack* track = (HoughTrack*)tpt;
  if (track->GetWeight() < GetWeight())
    return 1;
  if (track->GetWeight() > GetWeight())
    return -1;
  return 0;
}

void HoughTrack::SetEta(double f)
{
  // Set eta, and calculate fTanl, which is the tan of dipangle

  fEta = f;
  double theta = 2 * atan(exp(-1. * fEta));
  double dipangle = AliceO2::Hough::Transform::PiHalf() - theta;
  double tgl = tan(dipangle);
  SetTgl(tgl);
}

void HoughTrack::UpdateToFirstRow()
{
  // Update the track parameters to the point where track cross
  // its first padrow.`

  // Get the crossing point with the first padrow:
  float xyz[3];
  if (!GetCrossingPoint(GetFirstRow(), xyz))
    cout << "Track does not cross padrow " << GetFirstRow() << " centerx " << GetCenterX() << " centery "
         << GetCenterY() << " Radius " << GetRadius() << " tgl " << GetTgl() << endl;

  // printf("Track with eta %f tgl %f crosses at x %f y %f z %f on padrow
  // %d\n",GetEta(),GetTgl(),xyz[0],xyz[1],xyz[2],GetFirstRow());
  // printf("Before: first %f %f %f tgl %f center %f %f charge
  // %d\n",GetFirstPointX(),GetFirstPointY(),GetFirstPointZ(),GetTgl(),GetCenterX(),GetCenterY(),GetCharge());

  double radius = sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1]);

  // Get the track parameters

  /*
    double x0    = GetR0() * cos(GetPhi0()) ;
    double y0    = GetR0() * sin(GetPhi0()) ;
  */
  double rc = GetRadius(); // fabs(GetPt()) / AliceO2::Hough::Transform::GetBFieldValue();
  double tPhi0 = GetPsi() + GetCharge() * AliceO2::Hough::Transform::PiHalf() / abs(GetCharge());
  double xc = GetCenterX(); // x0 - rc * cos(tPhi0) ;
  double yc = GetCenterY(); // y0 - rc * sin(tPhi0) ;

  // Check helix and cylinder intersect
  double fac1 = xc * xc + yc * yc;
  double sfac = sqrt(fac1);

  if (fabs(sfac - rc) > radius || fabs(sfac + rc) < radius) {
    cerr << "Track does not intersect" << endl;
    return;
  }

  // Find intersection
  double fac2 = (radius * radius + fac1 - rc * rc) / (2.00 * radius * sfac);
  double phi = atan2(yc, xc) + GetCharge() * acos(fac2);
  double td = atan2(radius * sin(phi) - yc, radius * cos(phi) - xc);

  // Intersection in z
  if (td < 0)
    td = td + AliceO2::Hough::Transform::TwoPi();
  double deltat = fmod((-GetCharge() * td + GetCharge() * tPhi0), AliceO2::Hough::Transform::TwoPi());
  if (deltat < 0.)
    deltat += AliceO2::Hough::Transform::TwoPi();
  else if (deltat > AliceO2::Hough::Transform::TwoPi())
    deltat -= AliceO2::Hough::Transform::TwoPi();
  double z = GetZ0() + rc * GetTgl() * deltat;

  double xExtra = radius * cos(phi);
  double yExtra = radius * sin(phi);

  double tPhi = atan2(yExtra - yc, xExtra - xc);

  // if ( tPhi < 0 ) tPhi += 2. * M_PI ;
  double tPsi = tPhi - GetCharge() * AliceO2::Hough::Transform::PiHalf() / abs(GetCharge());
  if (tPsi > AliceO2::Hough::Transform::TwoPi())
    tPsi -= AliceO2::Hough::Transform::TwoPi();
  else if (tPsi < 0.)
    tPsi += AliceO2::Hough::Transform::TwoPi();

  // And finally, update the track parameters
  SetR0(radius);
  SetPhi0(phi);
  SetZ0(z);
  SetPsi(tPsi);
  SetFirstPoint(xyz[0], xyz[1], z);
  // printf("After: first %f %f %f tgl %f center %f %f charge
  // %d\n",GetFirstPointX(),GetFirstPointY(),GetFirstPointZ(),GetTgl(),GetCenterX(),GetCenterY(),GetCharge());

  // printf("First point set %f %f %f\n",xyz[0],xyz[1],z);

  // Also, set the coordinates of the point where track crosses last padrow:
  GetCrossingPoint(GetLastRow(), xyz);
  SetLastPoint(xyz[0], xyz[1], xyz[2]);
  // printf("last point %f %f %f\n",xyz[0],xyz[1],xyz[2]);
}

void HoughTrack::SetTrackParameters(double kappa, double eangle, int weight)
{
  // Set track parameters - sort of ctor
  fWeight = weight;
  fMinDist = 100000;
  SetKappa(kappa);
  double pt = fabs(AliceO2::Hough::Transform::GetBFieldValue() / kappa);
  SetPt(pt);
  double radius = 1 / fabs(kappa);
  SetRadius(radius);
  SetFirstPoint(0, 0, 0);
  SetPsi(eangle); // Psi = emission angle when first point is vertex
  SetPhi0(0);     // not defined for vertex reference point
  SetR0(0);
  double charge = -1. * kappa;
  SetCharge((int)copysign(1., charge));
  double trackPhi0 = GetPsi() + charge * 0.5 * AliceO2::Hough::Transform::Pi() / fabs(charge);
  double xc = GetFirstPointX() - GetRadius() * cos(trackPhi0);
  double yc = GetFirstPointY() - GetRadius() * sin(trackPhi0);
  SetCenterX(xc);
  SetCenterY(yc);
  SetNHits(1); // just for the trackarray IO
  fIsHelix = true;
}

void HoughTrack::SetTrackParametersRow(double alpha1, double alpha2, double eta, int weight)
{
  // Set track parameters for TransformerRow
  // This includes curvature,emission angle and eta
  double psi = atan((alpha1 - alpha2) / (TransformerRow::GetBeta1() - TransformerRow::GetBeta2()));
  double kappa = 2.0 * (alpha1 * cos(psi) - TransformerRow::GetBeta1() * sin(psi));
  SetTrackParameters(kappa, psi, weight);

  double zovr;
  double etaparam1 = TransformerRow::GetEtaCalcParam1();
  double etaparam2 = TransformerRow::GetEtaCalcParam2();
  if (eta > 0)
    zovr = (etaparam1 - sqrt(etaparam1 * etaparam1 - 4. * etaparam2 * eta)) / (2. * etaparam2);
  else
    zovr = -1. * (etaparam1 - sqrt(etaparam1 * etaparam1 + 4. * etaparam2 * eta)) / (2. * etaparam2);
  double r = sqrt(1. + zovr * zovr);
  double exacteta = 0.5 * log((1 + zovr / r) / (1 - zovr / r));
  SetEta(exacteta);
}

void HoughTrack::SetLineParameters(double psi, double D, int weight, int* rowrange, int /*ref_row*/)
{
  // Initialize a track piece, not yet a track
  // Used in case of straight line transformation

  // Transform line parameters to coordinate system of slice:

  //  D = D + fTransform->Row2X(ref_row)*cos(psi);

  fDLine = D;
  fPsiLine = psi;
  fWeight = weight;
  SetNHits(1);
  SetRowRange(rowrange[0], rowrange[1]);
  fIsHelix = false;
}

void HoughTrack::SetBestMCid(int mcid, double mindist)
{
  // Finds and set the closest mc label
  if (mindist < fMinDist) {
    fMinDist = mindist;
    SetMCid(mcid);
  }
}

void HoughTrack::GetLineCrossingPoint(int padrow, float* xy)
{
  // Returns the crossing point of the track with a given padrow
  if (fIsHelix) {
    printf("HoughTrack::GetLineCrossingPoint : Track is not a line\n");
    return;
  }

  float xhit = AliceO2::Hough::Transform::Row2X(padrow) - AliceO2::Hough::Transform::Row2X(GetFirstRow());
  float a = -1 / tan(fPsiLine);
  float b = fDLine / sin(fPsiLine);
  float yhit = a * xhit + b;
  xy[0] = xhit;
  xy[1] = yhit;
}
