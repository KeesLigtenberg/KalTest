//*************************************************************************
//* =====================
//*  TKalTrackSite Class
//* =====================
//*
//* (Description)
//*   Track measurement site class used by Kalman filter.
//* (Requires)
//*     TKalTrackState
//* (Provides)
//*     class TKalTrackSite
//* (Update Recored)
//*   2003/09/30  Y.Nakashima       Original version.
//*   2010/04/06  K.Fujii           Modified a c-tor to allow a 1-dim hit,
//*                                 for which pivot is at the xpected hit.
//*                                 Modified IsAccepted() to allow user-
//*                                 defined filter conditions.
//*
//*************************************************************************

#include "TKalTrackSite.h"    // from KalTrackLib
#include "TKalTrackState.h"   // from KalTrackLib
#include "TVTrackHit.h"       // from KalTrackLib
#include "TKalFilterCond.h"   // from KalTrackLib
#include "TVSurface.h"        // from GeomLib
#include "TBField.h"          // from Bfield

#include <iostream>           // from STL
#include <memory>             // from STL

using namespace std;

//_________________________________________________________________________
//  ----------------------------------
//   Class for measurement site
//  ----------------------------------
//
ClassImp(TKalTrackSite)

//_________________________________________________________________________
//  ----------------------------------
//  Ctors and Dtor
//  ----------------------------------
TKalTrackSite::TKalTrackSite(Int_t m, Int_t p)
              : TVKalSite(m,p), fHitPtr(0), fX0(), fIsHitOwner(kFALSE),
                fCondPtr(0)
{
	if(!TBField::IsUsingUniformBfield()) {
	   fBfield = TBField::GetGlobalBfield(fX0).Mag(); // wrong for 1-dim hit.
	}
}

TKalTrackSite::TKalTrackSite(const TVTrackHit &ht,
                                   Int_t       p)
              : TVKalSite(ht.GetDimension(),p),
                fHitPtr(static_cast<const TVTrackHit *>(&ht)), 
                fX0(),
                fIsHitOwner(kFALSE),
                fCondPtr(0)
{
   for (Int_t i=0; i<ht.GetDimension(); i++) {
      GetMeasVec     ()(i,0) = ht.GetX(i);
      GetMeasNoiseMat()(i,i) = TMath::Power(ht.GetDX(i),2);
   }
   // Leave the pivot at the origin for a 1-dim hit
   if (ht.GetDimension() > 1) fX0 = ht.GetMeasLayer().HitToXv(ht);

   // B field is unkonwn for a 1-dim hit untill prediction is made.
   // Temporarily set it to that at the origin.
   if(!TBField::IsUsingUniformBfield()) { 
	fBfield = TBField::GetGlobalBfield(fX0).Mag();
   }
}

TKalTrackSite::~TKalTrackSite()
{
   if (IsHitOwner() && fHitPtr) delete fHitPtr;
}

//_________________________________________________________________________
//  ----------------------------------
//  Implementation of public methods
//  ----------------------------------

TVKalState & TKalTrackSite::CreateState(const TKalMatrix &sv, Int_t type)
{
   SetOwner();
   return *(new TKalTrackState(sv,*this,type));
}

TVKalState & TKalTrackSite::CreateState(const TKalMatrix &sv,
                                        const TKalMatrix &c,
                                              Int_t       type)
{
   SetOwner();
   return *(new TKalTrackState(sv,c,*this,type));
}

Int_t TKalTrackSite::CalcXexp(const TVKalState &a, 
                                    TVector3   &xx, 
                                    Double_t   &phi) const
{
   std::auto_ptr<TVTrack> hel(&static_cast<const TKalTrackState &>(a).CreateTrack());

   const TVSurface &ms = dynamic_cast<const TVSurface &>(GetHit().GetMeasLayer());

   if(!TBField::IsUsingUniformBfield()) {
	   const double eps = 1.e-5; 
	   Int_t intersect=ms.CalcXingPointWith(*hel,xx,phi,eps);
	   if(!intersect) //use closest point if no intersection
		   return ms.calcClosestPointWith(*hel,xx,phi,eps);
	   return intersect;
   }
   else {
	   Int_t intersect = ms.CalcXingPointWith(*hel,xx,phi);
	   if(!intersect)
		   return ms.calcClosestPointWith(*hel,xx,phi);
	   return intersect;
   }

}



Int_t TKalTrackSite::CalcExpectedMeasVec(const TVKalState &a, TKalMatrix &h)
{
   Double_t phi = 0.;
   TVector3 xxv;
   if (!CalcXexp(a,xxv,phi)) return 0;	// no hit

   if (a.GetNrows() == 6) h = GetHit().XvToMv(xxv,a(5,0),a);
   else                   h = GetHit().XvToMv(xxv,0., a);

   return 1;
}

Int_t TKalTrackSite::CalcMeasVecDerivative(const TVKalState &a,
                                                 TKalMatrix &H)
{
   // Calculate 
   //    H = (@h/@a) = (@d/@a, @z/@a)^t
   // where
   //        h(a) = (d, z)^t: expected meas vector
   //        a = (drho, phi0, kappa, dz, tanl, t0)
   //

   TVector3      xxv;
   Double_t      phi = 0.;

   int nCrossings=CalcXexp(a,xxv,phi);
   if (!nCrossings) return 0; // hit on S(x) = 0

   const TVSurface &ms = dynamic_cast<const TVSurface &>(GetHit().GetMeasLayer());
   TKalMatrix    dsdx(ms.CalcDSDx(xxv));   // (@S(x)/@x)

   if(nCrossings==-2) {//-2 indicates xxv was a closest point, not a crossing!
	   auto x=dsdx(0,0), y=dsdx(0,1);
	   dsdx(0,1)=x; dsdx(0,0)=-y; //use vector perpendicular in xy plane instead.
	   //fixme: If I understand correctly it does not matter how this vector is changed as long as it is not perpendicular to the difference. Or does this work only for cylinders along z?
   }

   std::auto_ptr<TVTrack> hel(&static_cast<const TKalTrackState &>(a).CreateTrack());

   TKalMatrix    dxda   = hel->CalcDxDa(phi);  	     // (@x(phi,a)/@a)
   TKalMatrix    dxdphi = hel->CalcDxDphi(phi);	     // (@x(phi,a)/@phi)

   TKalMatrix dphida = dsdx * dxda;
   TKalMatrix dsdphi = dsdx * dxdphi;
   Double_t denom = -dsdphi(0,0);
   dphida *= 1/denom;

   TKalMatrix dxphiada = dxdphi * dphida + dxda; // (@x(phi(a),a)/@a)

   GetHit().GetMeasLayer().CalcDhDa(GetHit(),a, xxv, dxphiada, H); // H = (@h/@a)

   return 1;
}

TVector3 TKalTrackSite::GetLocalPivot() const
{	
	if(!TBField::IsUsingUniformBfield()) {
		//get the local pviot in a non-uniform magnetic field
  		TVector3 localPivot = fFrame.Transform(fX0, TTrackFrame::kGlobalToLocal);
  	
  		return localPivot;
	}
	else {
		//return global pivot if magetic field is uniform
		return fX0;
	}
}

Bool_t TKalTrackSite::IsAccepted()
{
   if (fCondPtr) return fCondPtr->IsAccepted(*this);
   else          return kTRUE;
}

void TKalTrackSite::DebugPrint() const
{
   cout << " dchi2 = " << GetDeltaChi2()   << endl;
   cout << " res_d = " << (*(TKalTrackSite *)this).GetResVec()(0,0) << endl;
   cout << " res_z = " << (*(TKalTrackSite *)this).GetResVec()(1,0) << endl;
   fHitPtr->DebugPrint();
}

