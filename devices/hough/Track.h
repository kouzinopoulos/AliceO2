// @(#) $Id: Track.h 15995 2006-11-30 17:45:45Z hristov $

#ifndef ALIL3TRACK_H
#define ALIL3TRACK_H

class AliHLTVertex;
class AliHLTSpacePointData;

class Track {

public:
  Track();
  virtual ~Track();

  virtual void Set(Track* track);
  virtual int Compare(const Track* track) const;
  virtual void CalculateHelix();

  bool CalculateReferencePoint(double angle, double radius = 132); // Calculate Reference Point
  bool CalculateEdgePoint(double angle); // Calculate crossing point with line
  bool CalculatePoint(double xplane); // Calculate crossing point with X-plane
  bool IsPoint() { return fIsPoint; }
  double GetCrossingAngle(int padrow, int slice = -1);
  bool GetCrossingPoint(int padrow, float* xyz);
  double GetDistance(double /*x0*/, double /*x1*/) { return 0; }
  void UpdateToFirstPoint();
  // void GetClosestPoint(AliHLTVertex *vertex,double &closest_x,double &closest_y,double &closest_z);
  void Rotate(int slice, bool tolocal = false);
  bool IsLocal() const { return fIsLocal; }
  void Print() const;

  // getter
  double GetFirstPointX() const { return fFirstPoint[0]; }
  double GetFirstPointY() const { return fFirstPoint[1]; }
  double GetFirstPointZ() const { return fFirstPoint[2]; }
  double GetLastPointX() const { return fLastPoint[0]; }
  double GetLastPointY() const { return fLastPoint[1]; }
  double GetLastPointZ() const { return fLastPoint[2]; }

  double GetPointPsi() const { return fPointPsi; }
  double GetPointX() const { return fPoint[0]; }
  double GetPointY() const { return fPoint[1]; }
  double GetPointZ() const { return fPoint[2]; }

  double GetPt() const { return fPt; }
  double GetTgl() const { return fTanl; }
  double GetPsi() const { return fPsi; }
  double GetPhi0() const { return fPhi0; }
  double GetR0() const { return fR0; }
  double GetZ0() const { return fFirstPoint[2]; }
  float GetPID() const { return fPID; }

  double GetPterr() const { return fPterr; }
  double GetPsierr() const { return fPsierr; }
  double GetTglerr() const { return fTanlerr; }

  double GetKappa() const { return fKappa; }
  double GetRadius() const { return fRadius; }
  double GetCenterX() const { return fCenterX; }
  double GetCenterY() const { return fCenterY; }

  int GetNHits() const { return fNHits; }
  int GetNumberOfPoints() const { return fNHits; }
  bool ComesFromMainVertex() const { return fFromMainVertex; }

  double GetPx() const { return fPt * cos(fPsi); }
  double GetPy() const { return fPt * sin(fPsi); }
  double GetPz() const { return fPt * fTanl; }

  double GetP() const;
  double GetPseudoRapidity() const;
  double GetRapidity() const;

  int GetCharge() const { return fQ; }
  int GetMCid() const { return fMCid; }
  double GetLength() const { return fLength; }

  int GetFirstRow() const { return fRowRange[0]; }
  int GetLastRow() const { return fRowRange[1]; }
  int GetSector() const { return fSector; }

  unsigned int* GetHitNumbers() { return fHitNumbers; }

  // setter
  void SetPID(float pid) { fPID = pid; }
  void SetMCid(int f) { fMCid = f; }
  void SetFirstPoint(double f, double g, double h)
  {
    fFirstPoint[0] = f;
    fFirstPoint[1] = g;
    fFirstPoint[2] = h;
  }
  void SetLastPoint(double f, double g, double h)
  {
    fLastPoint[0] = f;
    fLastPoint[1] = g;
    fLastPoint[2] = h;
  }
  void SetHits(int nhits, unsigned int* hits) { memcpy(fHitNumbers, hits, nhits * sizeof(unsigned int)); }
  void SetPhi0(double f) { fPhi0 = f; }
  void SetPsi(double f) { fPsi = f; }
  void SetR0(double f) { fR0 = f; }
  void SetTgl(double f) { fTanl = f; }
  void SetZ0(double f) { fFirstPoint[2] = f; }
  void SetPt(double f) { fPt = f; }
  void SetLength(double f) { fLength = f; }
  void SetPterr(double f) { fPterr = f; }
  void SetPsierr(double f) { fPsierr = f; }
  void SetZ0err(double f) { fZ0err = f; }
  void SetTglerr(double f) { fTanlerr = f; }
  void SetKappa(double f) { fKappa = f; }
  void SetNHits(int f) { fNHits = f; }
  void SetRowRange(int f, int g)
  {
    fRowRange[0] = f;
    fRowRange[1] = g;
  }
  void SetSector(int f) { fSector = f; }
  void SetRadius(double f) { fRadius = f; }
  void SetCenterX(double f) { fCenterX = f; }
  void SetCenterY(double f) { fCenterY = f; }
  void SetCharge(int f) { fQ = f; }
  void ComesFromMainVertex(bool f) { fFromMainVertex = f; }

private:
  int fNHits; // Number of hits
  int fMCid; // Assigned id from MC data.

  double fKappa;        // Signed curvature (projected to a circle)
  double fRadius;       // Radius of the helix (projected to a circle)
  double fCenterX;      // x coordinate of the center of the helix (projected to a circle)
  double fCenterY;      // y coordinate of the center of the helix (projected to a circle)
  bool fFromMainVertex; // true if tracks origin is the main vertex, otherwise false

  int fRowRange[2]; // Subsector where this track was build
  int fSector; // Sector # where  this track was build

  // data from momentum fit
  int fQ; // charge measured fit

  // track parameters:
  double fTanl; // tan of dipangle
  double fPsi; // azimuthal angle of the momentum
  double fPt; // transverse momentum
  double fLength; // length of track (s)

  double fPterr; // error in pt
  double fPsierr; // error in psi
  double fZ0err; // error in first point
  double fTanlerr; // error in tanl

  double fPhi0; // azimuthal angle of the first point
  double fR0; // radius of the first point
  double fZ0; // z coordinate of the first point (fFirstPoint[2])

  double fFirstPoint[3]; // first point
  double fLastPoint[3]; // last point
  double fPoint[3]; // point
  double fPointPsi; // azimuthal angle of the momentum at Point

  bool fIsPoint; // Helix crosses the X-plane
  bool fIsLocal; // Track given in local coordinates.

  float fPID; // pid
  unsigned int fHitNumbers[159]; // Array of hit numbers for this track

  bool IsPoint(bool ispoint)
  {
    fIsPoint = ispoint;
    return fIsPoint;
  }
};

typedef Track AliL3Track; // for backward compatibility

#endif
