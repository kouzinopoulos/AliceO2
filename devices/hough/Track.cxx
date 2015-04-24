// @(#) $Id: Track.cxx 15995 2006-11-30 17:45:45Z hristov $

// Author: Anders Vestbo <mailto:vestbo$fi.uib.no>, Uli Frankenfeld <mailto:franken@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliHLTStandardIncludes.h"
#include "Track.h"
#include "Transform.h"
/*#include "AliHLTVertex.h"
#include "AliHLTSpacePointData.h"
*/

#if __GNUC__ >= 3
using namespace std;
#endif

Track::Track()
{
  // Constructor
  fNHits = 0;
  fMCid = -1;
  fKappa = 0;
  fRadius = 0;
  fCenterX = 0;
  fCenterY = 0;
  ComesFromMainVertex(false);
  fQ = 0;
  fPhi0 = 0;
  fPsi = 0;
  fR0 = 0;
  fTanl = 0;
  fZ0 = 0;
  fPt = 0;
  fLength = 0;
  fIsLocal = true;
  fRowRange[0] = 0;
  fRowRange[1] = 0;
  SetFirstPoint(0, 0, 0);
  SetLastPoint(0, 0, 0);
  memset(fHitNumbers, 0, 159 * sizeof(unsigned int));
  fPID = 0;

  fSector = 0;
  fPterr = 0;
  fPsierr = 0;
  fZ0err = 0;
  fTanlerr = 0;
  fPoint[0] = fPoint[1] = fPoint[2] = 0;
  fPointPsi = 0;
}

void Track::Set(Track* tpt)
{
  // setter
  SetRowRange(tpt->GetFirstRow(), tpt->GetLastRow());
  SetPhi0(tpt->GetPhi0());
  SetKappa(tpt->GetKappa());
  SetNHits(tpt->GetNHits());
  SetFirstPoint(tpt->GetFirstPointX(), tpt->GetFirstPointY(), tpt->GetFirstPointZ());
  SetLastPoint(tpt->GetLastPointX(), tpt->GetLastPointY(), tpt->GetLastPointZ());
  SetPt(tpt->GetPt());
  SetPsi(tpt->GetPsi());
  SetTgl(tpt->GetTgl());
  SetPterr(tpt->GetPterr());
  SetPsierr(tpt->GetPsierr());
  SetTglerr(tpt->GetTglerr());
  SetCharge(tpt->GetCharge());
  SetHits(tpt->GetNHits(), (unsigned int*)tpt->GetHitNumbers());
#ifdef do_mc
  SetMCid(tpt->GetMCid());
#endif
  SetPID(tpt->GetPID());
  SetSector(tpt->GetSector());
}

int Track::Compare(const Track* track) const
{
  // compare tracks
  if (track->GetNHits() < GetNHits())
    return 1;
  if (track->GetNHits() > GetNHits())
    return -1;
  return 0;
}

Track::~Track()
{
  // Nothing to do
}

double Track::GetP() const
{
  // Returns total momentum.
  return fabs(GetPt()) * sqrt(1. + GetTgl() * GetTgl());
}

double Track::GetPseudoRapidity() const { /* get pseudo rap return 0.5 * log((GetP() + GetPz()) / (GetP() - GetPz()));*/ }

/*
double Track::GetEta() const
{
  return GetPseudoRapidity();
}
*/

double Track::GetRapidity() const
{
  // get rap
  const double kmpi = 0.13957;
  return 0.5 * log((kmpi + GetPz()) / (kmpi - GetPz()));
}

void Track::Rotate(int slice, bool tolocal)
{
  // Rotate track to global parameters
  // If flag tolocal is set, the track is rotated
  // to local coordinates.

  float psi[1] = { GetPsi() };
  if (!tolocal)
    AliceO2::Hough::Transform::Local2GlobalAngle(psi, slice);
  else
    AliceO2::Hough::Transform::Global2LocalAngle(psi, slice);
  SetPsi(psi[0]);
  float first[3];
  first[0] = GetFirstPointX();
  first[1] = GetFirstPointY();
  first[2] = GetFirstPointZ();
  if (!tolocal)
    AliceO2::Hough::Transform::Local2Global(first, slice);
  else
    AliceO2::Hough::Transform::Global2LocHLT(first, slice);
  // AliceO2::Hough::Transform::Global2Local(first,slice,true);

  SetFirstPoint(first[0], first[1], first[2]);
  float last[3];
  last[0] = GetLastPointX();
  last[1] = GetLastPointY();
  last[2] = GetLastPointZ();
  if (!tolocal)
    AliceO2::Hough::Transform::Local2Global(last, slice);
  else
    AliceO2::Hough::Transform::Global2LocHLT(last, slice);
  // AliceO2::Hough::Transform::Global2Local(last,slice,true);
  SetLastPoint(last[0], last[1], last[2]);

  float center[3] = { GetCenterX(), GetCenterY(), 0 };
  if (!tolocal)
    AliceO2::Hough::Transform::Local2Global(center, slice);
  else
    AliceO2::Hough::Transform::Global2LocHLT(center, slice);
  // AliceO2::Hough::Transform::Global2Local(center,slice,true);
  SetCenterX(center[0]);
  SetCenterY(center[1]);

  SetPhi0(atan2(fFirstPoint[1], fFirstPoint[0]));
  SetR0(sqrt(fFirstPoint[0] * fFirstPoint[0] + fFirstPoint[1] * fFirstPoint[1]));

  if (!tolocal)
    fIsLocal = false;
  else
    fIsLocal = true;
}

void Track::CalculateHelix()
{
  // Calculate Radius, CenterX and CenterY from Psi, X0, Y0
  fRadius = fPt / (AliceO2::Hough::Transform::GetBFieldValue());
  if (fRadius)
    fKappa = -fQ * 1. / fRadius;
  else
    fRadius = 999999; // just zero
  double trackPhi0 = fPsi + fQ * AliceO2::Hough::Transform::PiHalf();

  fCenterX = fFirstPoint[0] - fRadius * cos(trackPhi0);
  fCenterY = fFirstPoint[1] - fRadius * sin(trackPhi0);

  SetPhi0(atan2(fFirstPoint[1], fFirstPoint[0]));
  SetR0(sqrt(fFirstPoint[0] * fFirstPoint[0] + fFirstPoint[1] * fFirstPoint[1]));
}

double Track::GetCrossingAngle(int padrow, int slice)
{
  // Calculate the crossing angle between track and given padrow.
  // Take the dot product of the tangent vector of the track, and
  // vector perpendicular to the padrow.
  // In order to do this, we need the tangent vector to the track at the
  // point. This is done by rotating the radius vector by 90 degrees;
  // rotation matrix: (  0  1 )
  //                 ( -1  0 )

  float angle = 0; // Angle perpendicular to the padrow in local coordinates
  if (slice >= 0) // Global coordinates
  {
    AliceO2::Hough::Transform::Local2GlobalAngle(&angle, slice);
    if (!CalculateReferencePoint(angle, AliceO2::Hough::Transform::Row2X(padrow)))
      cerr << "Track::GetCrossingAngle : Track does not cross line in slice " << slice << " row " << padrow << endl;
  } else // should be in local coordinates
  {
    float xyz[3];
    GetCrossingPoint(padrow, xyz);
    fPoint[0] = xyz[0];
    fPoint[1] = xyz[1];
    fPoint[2] = xyz[2];
  }

  double tangent[2];

  tangent[0] = (fPoint[1] - GetCenterY()) / GetRadius();
  tangent[1] = -1. * (fPoint[0] - GetCenterX()) / GetRadius();

  double perppadrow[2] = { cos(angle), sin(angle) };
  double cosbeta = fabs(tangent[0] * perppadrow[0] + tangent[1] * perppadrow[1]);
  if (cosbeta > 1)
    cosbeta = 1;
  return acos(cosbeta);
}

bool Track::GetCrossingPoint(int padrow, float* xyz)
{
  // Assumes the track is given in local coordinates

  if (!IsLocal()) {
    cerr << "GetCrossingPoint: Track is given on global coordinates" << endl;
    return false;
  }

  double xHit = AliceO2::Hough::Transform::Row2X(padrow);

  xyz[0] = xHit;
  double aa = (xHit - GetCenterX()) * (xHit - GetCenterX());
  double r2 = GetRadius() * GetRadius();
  if (aa > r2)
    return false;

  double aa2 = sqrt(r2 - aa);
  double y1 = GetCenterY() + aa2;
  double y2 = GetCenterY() - aa2;
  xyz[1] = y1;
  if (fabs(y2) < fabs(y1))
    xyz[1] = y2;

  double yHit = xyz[1];
  double angle1 = atan2((yHit - GetCenterY()), (xHit - GetCenterX()));
  if (angle1 < 0)
    angle1 += 2. * AliceO2::Hough::Transform::Pi();
  double angle2 = atan2((GetFirstPointY() - GetCenterY()), (GetFirstPointX() - GetCenterX()));
  if (angle2 < 0)
    angle2 += AliceO2::Hough::Transform::TwoPi();
  double diffangle = angle1 - angle2;
  diffangle = fmod(diffangle, AliceO2::Hough::Transform::TwoPi());
  if ((GetCharge() * diffangle) > 0)
    diffangle = diffangle - GetCharge() * AliceO2::Hough::Transform::TwoPi();
  double stot = fabs(diffangle) * GetRadius();
  double zHit = GetFirstPointZ() + stot * GetTgl();
  xyz[2] = zHit;

  return true;
}

bool Track::CalculateReferencePoint(double angle, double radius)
{
  // Global coordinate: crossing point with y = ax+ b;
  // a=tan(angle-AliceO2::Hough::Transform::PiHalf());
  //
  const double krr = radius; // position of reference plane
  const double kxr = cos(angle) * krr;
  const double kyr = sin(angle) * krr;

  double a = tan(angle - AliceO2::Hough::Transform::PiHalf());
  double b = kyr - a * kxr;

  double pp = (fCenterX + a * fCenterY - a * b) / (1 + pow(a, 2));
  double qq = (pow(fCenterX, 2) + pow(fCenterY, 2) - 2 * fCenterY * b + pow(b, 2) - pow(fRadius, 2)) / (1 + pow(a, 2));

  double racine = pp * pp - qq;
  if (racine < 0)
    return IsPoint(false); // no Point

  double rootRacine = sqrt(racine);
  double x0 = pp + rootRacine;
  double x1 = pp - rootRacine;
  double y0 = a * x0 + b;
  double y1 = a * x1 + b;

  double diff0 = sqrt(pow(x0 - kxr, 2) + pow(y0 - kyr, 2));
  double diff1 = sqrt(pow(x1 - kxr, 2) + pow(y1 - kyr, 2));

  if (diff0 < diff1) {
    fPoint[0] = x0;
    fPoint[1] = y0;
  } else {
    fPoint[0] = x1;
    fPoint[1] = y1;
  }

  double pointPhi0 = atan2(fPoint[1] - fCenterY, fPoint[0] - fCenterX);
  double trackPhi0 = atan2(fFirstPoint[1] - fCenterY, fFirstPoint[0] - fCenterX);
  if (fabs(trackPhi0 - pointPhi0) > AliceO2::Hough::Transform::Pi()) {
    if (trackPhi0 < pointPhi0)
      trackPhi0 += AliceO2::Hough::Transform::TwoPi();
    else
      pointPhi0 += AliceO2::Hough::Transform::TwoPi();
  }
  double stot = -fQ * (pointPhi0 - trackPhi0) * fRadius;
  fPoint[2] = fFirstPoint[2] + stot * fTanl;

  fPointPsi = pointPhi0 - fQ * AliceO2::Hough::Transform::PiHalf();
  if (fPointPsi < 0.)
    fPointPsi += AliceO2::Hough::Transform::TwoPi();
  fPointPsi = fmod(fPointPsi, AliceO2::Hough::Transform::TwoPi());

  return IsPoint(true);
}

bool Track::CalculateEdgePoint(double angle)
{
  // Global coordinate: crossing point with y = ax; a=tan(angle);
  //
  double rmin = AliceO2::Hough::Transform::Row2X(AliceO2::Hough::Transform::GetFirstRow(-1)); // min Radius of TPC
  double rmax = AliceO2::Hough::Transform::Row2X(AliceO2::Hough::Transform::GetLastRow(-1)); // max Radius of TPC

  double a = tan(angle);
  double pp = (fCenterX + a * fCenterY) / (1 + pow(a, 2));
  double qq = (pow(fCenterX, 2) + pow(fCenterY, 2) - pow(fRadius, 2)) / (1 + pow(a, 2));
  double racine = pp * pp - qq;
  if (racine < 0)
    return IsPoint(false); // no Point
  double rootRacine = sqrt(racine);
  double x0 = pp + rootRacine;
  double x1 = pp - rootRacine;
  double y0 = a * x0;
  double y1 = a * x1;

  double r0 = sqrt(pow(x0, 2) + pow(y0, 2));
  double r1 = sqrt(pow(x1, 2) + pow(y1, 2));
  // find the right crossing point:
  // inside the TPC modules
  bool ok0 = false;
  bool ok1 = false;

  if (r0 > rmin && r0 < rmax) {
    double da = atan2(y0, x0);
    if (da < 0)
      da += AliceO2::Hough::Transform::TwoPi();
    if (fabs(da - angle) < 0.5)
      ok0 = true;
  }
  if (r1 > rmin && r1 < rmax) {
    double da = atan2(y1, x1);
    if (da < 0)
      da += AliceO2::Hough::Transform::TwoPi();
    if (fabs(da - angle) < 0.5)
      ok1 = true;
  }
  if (!(ok0 || ok1))
    return IsPoint(false); // no Point

  if (ok0 && ok1) {
    double diff0 = sqrt(pow(fFirstPoint[0] - x0, 2) + pow(fFirstPoint[1] - y0, 2));
    double diff1 = sqrt(pow(fFirstPoint[0] - x1, 2) + pow(fFirstPoint[1] - y1, 2));
    if (diff0 < diff1)
      ok1 = false; // use ok0
    else
      ok0 = false; // use ok1
  }
  if (ok0) {
    fPoint[0] = x0;
    fPoint[1] = y0;
  } else {
    fPoint[0] = x1;
    fPoint[1] = y1;
  }

  double pointPhi0 = atan2(fPoint[1] - fCenterY, fPoint[0] - fCenterX);
  double trackPhi0 = atan2(fFirstPoint[1] - fCenterY, fFirstPoint[0] - fCenterX);
  if (fabs(trackPhi0 - pointPhi0) > AliceO2::Hough::Transform::Pi()) {
    if (trackPhi0 < pointPhi0)
      trackPhi0 += AliceO2::Hough::Transform::TwoPi();
    else
      pointPhi0 += AliceO2::Hough::Transform::TwoPi();
  }
  double stot = -fQ * (pointPhi0 - trackPhi0) * fRadius;
  fPoint[2] = fFirstPoint[2] + stot * fTanl;

  fPointPsi = pointPhi0 - fQ * AliceO2::Hough::Transform::PiHalf();
  if (fPointPsi < 0.)
    fPointPsi += AliceO2::Hough::Transform::TwoPi();
  fPointPsi = fmod(fPointPsi, AliceO2::Hough::Transform::TwoPi());

  return IsPoint(true);
}

bool Track::CalculatePoint(double xplane)
{
  // Local coordinate: crossing point with x plane
  //
  double racine = pow(fRadius, 2) - pow(xplane - fCenterX, 2);
  if (racine < 0)
    return IsPoint(false);
  double rootRacine = sqrt(racine);

  double y0 = fCenterY + rootRacine;
  double y1 = fCenterY - rootRacine;
  // double diff0 = sqrt(pow(fFirstPoint[0]-xplane)+pow(fFirstPoint[1]-y0));
  // double diff1 = sqrt(pow(fFirstPoint[0]-xplane)+pow(fFirstPoint[1]-y1));
  double diff0 = fabs(y0 - fFirstPoint[1]);
  double diff1 = fabs(y1 - fFirstPoint[1]);

  fPoint[0] = xplane;
  if (diff0 < diff1)
    fPoint[1] = y0;
  else
    fPoint[1] = y1;

  double pointPhi0 = atan2(fPoint[1] - fCenterY, fPoint[0] - fCenterX);
  double trackPhi0 = atan2(fFirstPoint[1] - fCenterY, fFirstPoint[0] - fCenterX);
  if (fabs(trackPhi0 - pointPhi0) > AliceO2::Hough::Transform::Pi()) {
    if (trackPhi0 < pointPhi0)
      trackPhi0 += AliceO2::Hough::Transform::TwoPi();
    else
      pointPhi0 += AliceO2::Hough::Transform::TwoPi();
  }
  double stot = -fQ * (pointPhi0 - trackPhi0) * fRadius;
  fPoint[2] = fFirstPoint[2] + stot * fTanl;

  fPointPsi = pointPhi0 - fQ * AliceO2::Hough::Transform::PiHalf();
  if (fPointPsi < 0.)
    fPointPsi += AliceO2::Hough::Transform::TwoPi();
  fPointPsi = fmod(fPointPsi, AliceO2::Hough::Transform::TwoPi());

  return IsPoint(true);
}

void Track::UpdateToFirstPoint()
{
  // Update track parameters to the innermost point on the track.
  // This means that the parameters of the track will be given in the point
  // of closest approach to the first innermost point, i.e. the point
  // lying on the track fit (and not the coordinates of the innermost point itself).
  // This function assumes that fFirstPoint is already set to the coordinates of the innermost
  // assigned cluster.
  //
  // During the helix-fit, the first point on the track is set to the coordinates
  // of the innermost assigned cluster. This may be ok, if you just want a fast
  // estimate of the "global" track parameters; such as the momentum etc.
  // However, if you later on want to do more precise local calculations, such
  // as impact parameter, residuals etc, you need to give the track parameters
  // according to the actual fit.

  double xc = GetCenterX() - GetFirstPointX();
  double yc = GetCenterY() - GetFirstPointY();

  double distx1 = xc * (1 + GetRadius() / sqrt(xc * xc + yc * yc));
  double disty1 = yc * (1 + GetRadius() / sqrt(xc * xc + yc * yc));
  double distance1 = sqrt(distx1 * distx1 + disty1 * disty1);

  double distx2 = xc * (1 - GetRadius() / sqrt(xc * xc + yc * yc));
  double disty2 = yc * (1 - GetRadius() / sqrt(xc * xc + yc * yc));
  double distance2 = sqrt(distx2 * distx2 + disty2 * disty2);

  // Choose the closest:
  double point[2];
  if (distance1 < distance2) {
    point[0] = distx1 + GetFirstPointX();
    point[1] = disty1 + GetFirstPointY();
  } else {
    point[0] = distx2 + GetFirstPointX();
    point[1] = disty2 + GetFirstPointY();
  }

  double pointpsi = atan2(point[1] - GetCenterY(), point[0] - GetCenterX());
  pointpsi -= GetCharge() * AliceO2::Hough::Transform::PiHalf();
  if (pointpsi < 0)
    pointpsi += AliceO2::Hough::Transform::TwoPi();

  // Update the track parameters
  SetR0(sqrt(point[0] * point[0] + point[1] * point[1]));
  SetPhi0(atan2(point[1], point[0]));
  SetFirstPoint(point[0], point[1], GetZ0());
  SetPsi(pointpsi);
}
/*
void Track::GetClosestPoint(AliHLTVertex *vertex,double &closestx,double &closesty,double &closestz)
{
  //Calculate the point of closest approach to the vertex
  //This function calculates the minimum distance from the helix to the vertex, and choose
  //the corresponding point lying on the helix as the point of closest approach.

  double xc = GetCenterX() - vertex->GetX();
  double yc = GetCenterY() - vertex->GetY();

  double distx1 = xc*(1 + GetRadius()/sqrt(xc*xc + yc*yc));
  double disty1 = yc*(1 + GetRadius()/sqrt(xc*xc + yc*yc));
  double distance1 = sqrt(distx1*distx1 + disty1*disty1);

  double distx2 = xc*(1 - GetRadius()/sqrt(xc*xc + yc*yc));
  double disty2 = yc*(1 - GetRadius()/sqrt(xc*xc + yc*yc));
  double distance2 = sqrt(distx2*distx2 + disty2*disty2);

  //Choose the closest:
  if(distance1 < distance2)
    {
      closestx = distx1 + vertex->GetX();
      closesty = disty1 + vertex->GetY();
    }
  else
    {
      closestx = distx2 + vertex->GetX();
      closesty = disty2 + vertex->GetY();
    }

  //Get the z coordinate:
  double angle1 = atan2((closesty-GetCenterY()),(closestx-GetCenterX()));
  if(angle1 < 0) angle1 = angle1 + AliceO2::Hough::Transform::TwoPi();

  double angle2 = atan2((GetFirstPointY()-GetCenterY()),(GetFirstPointX()-GetCenterX()));
  if(angle2 < 0) angle2 = angle2 + AliceO2::Hough::Transform::TwoPi();

  double diff_angle = angle1 - angle2;
  diff_angle = fmod(diff_angle,AliceO2::Hough::Transform::TwoPi());

  if((GetCharge()*diff_angle) < 0) diff_angle = diff_angle + GetCharge()*AliceO2::Hough::Transform::TwoPi();
  double stot = fabs(diff_angle)*GetRadius();
  closestz = GetFirstPointZ() - stot*GetTgl();
}
*/
void Track::Print() const
{ // print out parameters of track
  cout << fNHits << " " << fMCid << " " << fKappa << " " << fRadius << " " << fCenterX << " " << fCenterY << " "
       << fFromMainVertex << " " << fRowRange[0] << " " << fRowRange[1] << " " << fSector << " " << fQ << " " << fTanl
       << " " << fPsi << " " << fPt << " " << fLength << " " << fPterr << " " << fPsierr << " " << fZ0err << " "
       << fTanlerr << " " << fPhi0 << " " << fR0 << " " << fZ0 << " " << fFirstPoint[0] << " " << fFirstPoint[1] << " "
       << fFirstPoint[2] << " " << fLastPoint[0] << " " << fLastPoint[1] << " " << fLastPoint[2] << " " << fPoint[0]
       << " " << fPoint[1] << " " << fPoint[2] << " " << fPointPsi << " " << fIsPoint << " " << fIsLocal << " " << fPID
       << endl;
}
