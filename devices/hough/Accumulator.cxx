/// \file Accumulator.cxx
/// \brief Implementation of the Accumulator class
/// \author Anders Vestbo <mailto:vestbo@fi.uib.no>

#include "Accumulator.h"

//uncomment if you want overflow checks
//#define _IFON_

using namespace std;

using namespace AliceO2::Hough;

ClassImp(AliceO2::Hough::Accumulator)

Accumulator::Accumulator()
{
  // Default constructor
  fNxbins = 0;
  fNybins = 0;
  fNcells = 0;
  fXmin = 0;
  fYmin = 0;
  fXmax = 0;
  fYmax = 0;
  fBinwidthX = 0;
  fBinwidthY = 0;
  fFirstXbin = 0;
  fLastXbin = 0;
  fFirstYbin = 0;
  fLastYbin = 0;
  fEntries = 0;
  fContent = 0;
  fThreshold = 0;
#ifdef use_root
  fRootHisto = 0;
#endif
}

Accumulator::Accumulator(const Char_t *name,const Char_t */*id*/,
			       Int_t nxbin,Double_t xmin,Double_t xmax,
			       Int_t nybin,Double_t ymin,Double_t ymax)
{
  // Normal constructor
  strcpy(fName,name);

  fNxbins = nxbin;
  fNybins = nybin;
  fNcells = (nxbin+2)*(nybin+2);
  fXmin = xmin;
  fYmin = ymin;
  fXmax = xmax;
  fYmax = ymax;
  fBinwidthX = (fXmax - fXmin) / fNxbins;
  fBinwidthY = (fYmax - fYmin) / fNybins;

  fEntries = 0;
  fFirstXbin = 1;
  fFirstYbin = 1;
  fLastXbin = nxbin;
  fLastYbin = nybin;
#ifdef use_root
  fRootHisto = 0;
#endif
  fThreshold = 0;

  fContent = new Int_t[fNcells];
  Reset();
}

Accumulator::~Accumulator()
{
  //Destructor
  if(fContent)
    delete [] fContent;
#ifdef use_root
  if(fRootHisto)
    delete fRootHisto;
#endif
}

void Accumulator::Reset()
{
  // Reset histogram contents
  if(fContent)
    for(Int_t i=0; i<fNcells; i++) fContent[i] = 0;

  fEntries=0;
}

void Accumulator::Fill(Double_t x,Double_t y,Int_t weight)
{
  // Fill the weight into a bin which correspond to x and y
  Int_t bin = FindBin(x,y);
#ifdef _IFON_
  if(bin < 0)
    return;
#endif

  AddBinContent(bin,weight);
}

void Accumulator::Fill(Double_t x,Int_t ybin,Int_t weight)
{
  // Fill the weight into a bin which correspond to x and ybin
  Int_t xbin = FindXbin(x);
  Int_t bin = GetBin(xbin,ybin);
#ifdef _IFON_
  if(bin < 0)
    return;
#endif

  AddBinContent(bin,weight);
}

void Accumulator::Fill(Int_t xbin,Double_t y,Int_t weight)
{
  // Fill the weight into a bin which correspond to xbin and y
  Int_t ybin = FindYbin(y);
  Int_t bin = GetBin(xbin,ybin);
#ifdef _IFON_
  if(bin < 0)
    return;
#endif

  AddBinContent(bin,weight);
}

void Accumulator::Fill(Int_t xbin,Int_t ybin,Int_t weight)
{
  // Fill the weight into a bin which correspond to xbin and ybin
  Int_t bin = GetBin(xbin,ybin);
#ifdef _IFON_
  if(bin < 0)
    return;
#endif

  AddBinContent(bin,weight);
}

Int_t Accumulator::FindBin(Double_t x,Double_t y) const
{
  // Finds the bin which correspond to x and y
  Int_t xbin = FindXbin(x);
  Int_t ybin = FindYbin(y);
#ifdef _IFON_
  if(!xbin || !ybin)
    return -1;
#endif

  return GetBin(xbin,ybin);
}

Int_t Accumulator::FindLabelBin(Double_t x,Double_t y) const
{
  // Returns the corresponding bin with the mc labels
  Int_t xbin = FindXbin(x);
  Int_t ybin = FindYbin(y);
#ifdef _IFON_
  if(!xbin || !ybin)
    return -1;
#endif

  return GetLabelBin(xbin,ybin);
}

Int_t Accumulator::FindXbin(Double_t x) const
{
  // Finds the bin which correspond to x
  if(x < fXmin || x > fXmax)
    return 0;

  return 1 + (Int_t)(fNxbins*(x-fXmin)/(fXmax-fXmin));
}

Int_t Accumulator::FindYbin(Double_t y) const
{
  // Finds the bin which correspond to y
  if(y < fYmin || y > fYmax)
    return 0;

  return 1 + (Int_t)(fNybins*(y-fYmin)/(fYmax-fYmin));
}

Int_t Accumulator::GetBin(Int_t xbin,Int_t ybin) const
{
  // Returns the bin which correspond to xbin and ybin
  if(xbin < fFirstXbin || xbin > fLastXbin)
    return 0;
  if(ybin < fFirstYbin || ybin > fLastYbin)
    return 0;

  return xbin + ybin*(fNxbins+2);
}

Int_t Accumulator::GetLabelBin(Int_t xbin,Int_t ybin) const
{
  // Returns the corresponding bin with the mc labels
  if(xbin < fFirstXbin || xbin > fLastXbin)
    return -1;
  if(ybin < fFirstYbin || ybin > fLastYbin)
    return -1;

  return (Int_t)(xbin/2) + ((Int_t)(ybin/2))*((Int_t)((fNxbins+3)/2));
}

Int_t Accumulator::GetBinContent(Int_t bin) const
{
  // Return the bin content
  if(bin >= fNcells)
    {
      LOG(AliHLTLog::kError,"Accumulator::GetBinContent","array")<<AliHLTLog::kDec<<
	"bin out of range "<<bin<<ENDLOG;
      return 0;
    }

  if(fContent[bin] < fThreshold)
    return 0;
  return fContent[bin];
}

void Accumulator::SetBinContent(Int_t xbin,Int_t ybin,Int_t value)
{
  // Set bin content
  Int_t bin = GetBin(xbin,ybin);
#ifdef _IFON_
  if(bin == 0)
    return;
#endif

  SetBinContent(bin,value);
}

void Accumulator::SetBinContent(Int_t bin,Int_t value)
{
  // Set bin content

  if(bin >= fNcells)
    {
      LOG(AliHLTLog::kError,"Accumulator::SetBinContent","array")<<AliHLTLog::kDec<<
	"bin out of range "<<bin<<ENDLOG;
      return;
    }

  if(bin == 0)
    return;
  fContent[bin]=value;
}

void Accumulator::AddBinContent(Int_t xbin,Int_t ybin,Int_t weight)
{
  // Adds weight to bin content
  Int_t bin = GetBin(xbin,ybin);
#ifdef _IFON_
  if(bin == 0)
    return;
#endif

  AddBinContent(bin,weight);
}

void Accumulator::AddBinContent(Int_t bin,Int_t weight)
{
  // Adds weight to bin content
  if(bin < 0 || bin > fNcells)
    {
      LOG(AliHLTLog::kError,"Accumulator::AddBinContent","array")<<AliHLTLog::kDec<<
	"bin-value out of range "<<bin<<ENDLOG;
      return;
    }
  if(bin == 0)
    return;
  fEntries++;
  fContent[bin] += weight;
}

void Accumulator::Add(Accumulator *h1,Double_t /*weight*/)
{
  //Adding two histograms. Should be identical.

  if(!h1)
    {
      LOG(AliHLTLog::kError,"Accumulator::Add","Pointer")<<
	"Attempting to add a non-existing histogram"<<ENDLOG;
      return;
    }

  if(h1->GetNbinsX()!=fNxbins || h1->GetNbinsY()!=fNybins)
    {
      LOG(AliHLTLog::kError,"Accumulator::Add","array")<<
	"Mismatch in the number of bins "<<ENDLOG;
      return;
    }

  if(h1->GetFirstXbin()!=fFirstXbin || h1->GetLastXbin()!=fLastXbin ||
     h1->GetFirstYbin()!=fFirstYbin || h1->GetLastYbin()!=fLastYbin)
    {
      LOG(AliHLTLog::kError,"Accumulator::Add","array")<<
	"Mismatch in the bin numbering "<<ENDLOG;
      return;
    }

  for(Int_t bin=0; bin<fNcells; bin++)
    fContent[bin] += h1->GetBinContent(bin);

  fEntries += h1->GetNEntries();
}

Double_t Accumulator::GetBinCenterX(Int_t xbin) const
{
  // Returns the position of the center of a bin
  if(xbin < fFirstXbin || xbin > fLastXbin)
    {
      LOG(AliHLTLog::kError,"Accumulator::GetBinCenterX","xbin")
	<<"Bin-value out of range "<<xbin<<ENDLOG;
      return -1;
    }

  return fXmin + (xbin-0.5) * fBinwidthX;
}

Double_t Accumulator::GetBinCenterY(Int_t ybin) const
{
  // Returns the position of the center of a bin
  if(ybin < fFirstYbin || ybin > fLastYbin)
    {
      LOG(AliHLTLog::kError,"Accumulator::GetBinCenterY","ybin")
	<<"Bin-value out of range "<<ybin<<ENDLOG;
      return -1;
    }

  return fYmin + (ybin-0.5) * fBinwidthY;
}

Double_t Accumulator::GetPreciseBinCenterX(Float_t xbin) const
{
  // Returns the position of the center of a bin using precise values inside the bin
  if(xbin < (fFirstXbin-1.5) || xbin > (fLastXbin+1.5))
    {
      LOG(AliHLTLog::kError,"Accumulator::GetBinCenterX","xbin")
	<<"Bin-value out of range "<<xbin<<ENDLOG;
      return -1;
    }
  //  return fXmin + (xbin-1) * fBinwidthX + 0.5*fBinwidthX;
  return fXmin + (xbin-0.5) * fBinwidthX;
}

Double_t Accumulator::GetPreciseBinCenterY(Float_t ybin) const
{
  // Returns the position of the center of a bin using precise values inside the bin
  if(ybin < (fFirstYbin-1.5) || ybin > (fLastYbin+1.5))
    {
      LOG(AliHLTLog::kError,"Accumulator::GetBinCenterY","ybin")
	<<"Bin-value out of range "<<ybin<<ENDLOG;
      return -1;
    }
  //  return fYmin + (ybin-1) * fBinwidthY + 0.5*fBinwidthY;
  return fYmin + (ybin-0.5) * fBinwidthY;
}

void Accumulator::Draw(const Char_t *option)
{
  // Fill the contents of the corresponding ROOT histogram and draws it
#ifdef use_root
  if(!fRootHisto)
    CreateRootHisto();

  for(Int_t xbin=GetFirstXbin(); xbin<=GetLastXbin(); xbin++)
    {
      for(Int_t ybin=GetFirstYbin(); ybin<=GetLastYbin(); ybin++)
	{
	  Int_t bin = GetBin(xbin,ybin);
	  fRootHisto->Fill(GetBinCenterX(xbin),GetBinCenterY(ybin),GetBinContent(bin));
	}
    }

  //fRootHisto->SetStats(kFALSE);
  fRootHisto->Draw(option);
  return;
#else
  cerr<<"Accumulator::Draw : You need to compile with ROOT in order to draw histogram"<<endl;
#endif
}

void Accumulator::CreateRootHisto()
{
  // Create ROOT histogram out of Accumulator
#ifdef use_root
  fRootHisto = new TH2F(fName,"",fNxbins,fXmin,fXmax,fNybins,fYmin,fYmax);
  return;
#else
  cerr<<"Accumulator::CreateRootHisto : You need to compile with ROOT in order to create ROOT histogram"<<endl;
#endif
}

ofstream& operator<<(ofstream &o, const Accumulator &h)
{
  for(Int_t xbin=h.GetFirstXbin(); xbin<=h.GetLastXbin(); xbin++)
    {
      for(Int_t ybin=h.GetFirstYbin(); ybin<=h.GetLastYbin(); ybin++)
	{
	  Int_t bin = h.GetBin(xbin,ybin);
	  o << h.GetBinContent(bin) << " ";
	}
      o << endl;
    }
  return o;
}
