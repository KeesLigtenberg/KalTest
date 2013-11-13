//*************************************************************************
//* =====================
//*  TKalDetCradle Class
//* =====================
//*
//* (Description)
//*   A singleton class to hold information of detector system
//*   used in Kalman filter classes.
//* (Requires)
//*     TObjArray
//*     TVKalDetector
//* (Provides)
//*     class TKalDetCradle
//* (Update Recored)
//*   2005/02/23  A.Yamaguchi  	Original version.
//*   2005/08/14  K.Fujii        Removed CalcTable(), GetMeasLayerTable(),
//*                              GetPhiTable(), and GetDir() and added
//*                              Transport() to do their functions.
//*   2010/04/06  K.Fujii        Modified Transport() to allow a 1-dim hit,
//*                              for which pivot is at the expected hit.
//*   2012/11/29  K.Fujii        Moved GetEnergyLoss and CalcQms from
//*                              TKalDetCradle.
//*
//*************************************************************************

#include "TKalDetCradle.h"   // from KalTrackLib
#include "TVMeasLayer.h"     // from KalTrackLib
#include "TVKalDetector.h"   // from KalTrackLib
#include "TKalTrackSite.h"   // from KalTrackLib
#include "TKalTrackState.h"  // from KalTrackLib
#include "TKalTrack.h"       // from KalTrackLib
#include "TVSurface.h"       // from GeomLib
#include "TVTrack.h"         // from GeomLib
#include "J4Timer.h"
#include <memory>            // from STL
#include <iostream>          // from STL

ClassImp(TKalDetCradle)

//_________________________________________________________________________
//  ----------------------------------
//   Ctors and Dtor
//  ----------------------------------

TKalDetCradle::TKalDetCradle(Int_t n)
             : TObjArray(n), fIsMSON(kTRUE), fIsDEDXON(kTRUE),
               fDone(kFALSE), fIsClosed(kFALSE)
{
}

TKalDetCradle::~TKalDetCradle()
{
}

//_________________________________________________________________________
//  ----------------------------------
//   Utility Methods
//  ----------------------------------
//_________________________________________________________________________
// -----------------
//  Install
// -----------------
//    installs a sub-detector into this cradle.
//
void TKalDetCradle::Install(TVKalDetector &det)
{
    if (IsClosed()) {
        std::cerr << ">>>> Error!! >>>> TKalDetCradle::Install" << std::endl
        << "      Cradle already closed. Abort!!"     << std::endl;
        abort();
    }
    TIter next(&det);
    TObject *mlp = 0;  // measment layer pointer
    while ((mlp = next())) {
        Add(mlp);
        dynamic_cast<TAttElement *>(mlp)->SetParentPtr(&det);
        det.SetParentPtr(this);
    }
    fDone = kFALSE;
}

void TKalDetCradle::Transport(const TKalTrackSite  &from,   // site from
                                    TKalTrackSite  &to,     // site to
                                    TKalMatrix     &sv,     // state vector
                                    TKalMatrix     &F,      // propagator matrix
                                    TKalMatrix     &Q)      // process noise matrix
{
    // ---------------------------------------------------------------------
    //  Sort measurement layers in this cradle if not
    // ---------------------------------------------------------------------
#ifdef __J4TIMER__
   static int timerid = -1;
   J4Timer timer(timerid, "TKalDetCradle", "Transport");
   timer.Start();
#endif

    if (!fDone) Update();

    // ---------------------------------------------------------------------
    //  Move to site "to"
    // ---------------------------------------------------------------------

    std::auto_ptr<TVTrack> help(&static_cast<TKalTrackState &>
                                (from.GetCurState()).CreateTrack()); // tmp track
    
    const TVMeasLayer& ml_to = to.GetHit().GetMeasLayer() ;
    
    TVector3  x0; // local pivot at the "to" site
    this->Transport(from, ml_to, x0, sv, F, Q, help);

    TVTrack &hel = *help;
    
    // ---------------------------------------------------------------------
    //  Move pivot from last expected hit to actural hit at site "to"
    // ---------------------------------------------------------------------
    
    if (to.GetDimension() > 1) {
        
        double fid = 0.;
        Int_t sdim = sv.GetNrows();              // number of track parameters
        TKalMatrix DF(sdim, sdim);               // propagator matrix segment
        
        hel.MoveTo(to.GetGlobalPivot(), fid, &DF, 0, kFALSE);     // move pivot to actual hit (to)

        F = DF * F;                              // update F accordingly
        hel.PutInto(sv);                         // save updated hel to sv
        
    } else {
        TVector3 x0g = (hel.GetFrame()).Transform(x0, TTrackFrame::kLocalToGlobal);
        to.SetPivot(x0g);                         // if it is a 1-dim hit
    }

    to.SetFrame(hel.GetFrame());     // set frame
    to.SetBfield(hel.GetMagField()); // set B field

#ifdef __J4TIMER__
	timer.Stop();
#endif

#ifdef __DEBUG__
    std::cout << "b field of site to is: " << hel.GetMagField() << ", "
                                           << to.GetBfield() << std::endl;
    std::cerr << "local pivot of site to is: ";
    std::cerr << "x0 = (" << to.GetLocalPivot().X() << ", "
                          << to.GetLocalPivot().Y() << ", "
                          << to.GetLocalPivot().Z() << ")"
    << std::endl;
#endif
}

int TKalDetCradle::Transport(const TKalTrackSite  &from,  // site from
                             const TVMeasLayer    &ml_to, // layer to reach
                             TVector3       &x0,    // pivot for sv
                             TKalMatrix     &sv,    // state vector
                             TKalMatrix     &F,     // propagator matrix
                             TKalMatrix     &Q)     // process noise matrix
{
    if (!fDone) Update();
    
    std::auto_ptr<TVTrack> help(&static_cast<TKalTrackState &>
                                (from.GetCurState()).CreateTrack()); // tmp track
    return Transport(from, ml_to, x0, sv, F, Q, help);
}

//
//
//
//_________________________________________________________________________
// -----------------
//  Transport
// -----------------
//    transports state (sv) from site (from) to layer (ml_to), taking into
//    account multiple scattering and energy loss and updates state (sv),
//    fills pivot in x0, propagator matrix (F), and process noise matrix (Q).
//

int TKalDetCradle::Transport(const TKalTrackSite  &from,  // site from
                             const TVMeasLayer    &ml_to, // layer to reach
                                   TVector3       &x0,    // pivot for sv
                                   TKalMatrix     &sv,    // state vector
                                   TKalMatrix     &F,     // propagator matrix
                                   TKalMatrix     &Q,     // process noise matrix
                           std::auto_ptr<TVTrack> &help)  // pointer to update track object
{
    // ---------------------------------------------------------------------
    //  Sort measurement layers in this cradle if not
    // ---------------------------------------------------------------------
    if (!fDone) Update();
    
    
    // ---------------------------------------------------------------------
    //  Decide the mass of the current track
    // ---------------------------------------------------------------------
    static const Double_t kMpi = 0.13957018; // pion mass [GeV]
	const double eps = 1.e-5;
    TKalTrack *ktp  = static_cast<TKalTrack *>(from.GetParentPtr());
    Double_t   mass = ktp ? ktp->GetMass() : kMpi;
    
    // ---------------------------------------------------------------------
    //  Locate sites from and to in this cradle
    // ---------------------------------------------------------------------
    Int_t  fridx = from.GetHit().GetMeasLayer().GetIndex();      // index of site from
    Int_t  toidx = ml_to.GetIndex();                             // index of layer to
    Int_t  di    = fridx > toidx ? -1 : 1;                       // layer increment
    TVTrack &hel = *help;
    
    //=====================
    // FIXME
    //=====================
    TVector3 xfrom = from.GetGlobalPivot();         // get the referenece point
    TVector3 xto;                                   // reference point at destination to be returned by CalcXingPointWith
    Double_t fito = 0;                              // deflection angle to destination to be returned by CalcXingPointWith
    
    const TVSurface *sfp = dynamic_cast<const TVSurface *>(&ml_to);// surface at destination
    
    sfp->CalcXingPointWith(hel, xto, fito, 0, eps);                    // the default tolerance is used
    // as mode is 0 here the closest point crossing point is taken
    // this means that if we are at the top of a looping track
    // and the point to which we want to move is on the other side of
    // the loop but has a lower radius the transport will move down
    // through all layers and segfault on reaching index -1
    
    //   if( does_cross < 1 ) return does_cross ;
    
    TMatrixD dxdphi = hel.CalcDxDphi(fito);                       // tangent vector at destination surface
    TVector3 dxdphiv(dxdphi(0,0),dxdphi(1,0),dxdphi(2,0));        // convert matirix diagonal to vector
    //  Double_t cpa = hel.GetKappa();                                // get pt
    
    Bool_t isout = -fito*dxdphiv.Dot(sfp->GetOutwardNormal(xto)) < 0 ? kTRUE : kFALSE;  // out-going or in-coming at the destination surface
    //=====================
    // ENDFIXME
    //=====================
    
    TVector3 xx;                               // expected hit position vector
    Double_t fid     = 0.;                     // deflection angle from the last hit
    
    Int_t sdim = sv.GetNrows();                // number of track parameters
    F.UnitMatrix();                            // set the propagator matrix to the unit matrix
    Q.Zero();                                  // zero the noise matrix
    
    TKalMatrix DF(sdim, sdim);                 // propagator matrix segment
    
    // ---------------------------------------------------------------------
    //  Loop over layers and transport sv, F, and Q step by step
    // ---------------------------------------------------------------------
    Int_t ifr = fridx; // set index to the index of the intitial starting layer
    
    // here we make first make sure that the helix is at the crossing point of the current surface.
    // this is necessary to ensure that the material is only accounted for between fridx and toidx
    // otherwise it is possible to have inconsistencies with material treatment.
    // loop until we reach the index toidx, which is the surface we need to reach
    for (Int_t ito=fridx; (di>0 && ito<=toidx)||(di<0 && ito>=toidx); ito += di) {
        
        Double_t fid_temp = fid; // deflection angle from the last layer crossing
        
        int mode = ito!=fridx ? di : 0; // need to move to the from site as the helix may not be on the crossing point yet, meaning that the eloss and ms will be incorrectely attributed ...
        
        if (dynamic_cast<TVSurface *>(At(ito))->CalcXingPointWith(hel, xx, fid, mode,eps)) { // if we have a crossing point at this surface, note di specifies if we are moving forwards or backwards

            //=====================
            // FIXME
            //=====================
            static const Double_t kMergin = 1.0;
            // if the distance from the current crossing point to the starting point - kMergin(1mm) is greater than the distance from the destination to the starting point
            // this is needed to skip crossing points which come from the far side of the IP, for a cylinder this would not be a problem
            // but for the bounded planes it is perfectly posible due to the sorting in R
            // reset the deflection angle and skip this layer
            // this would at stop layers being added which are too far away but I am not sure how this will work with the problem described above.
            if( (xx-xfrom).Mag() - kMergin > (xto-xfrom).Mag() ){
                fid = fid_temp;
                continue ;
            }
            //=====================
            // ENDFIXME
            //=====================

            const TVMeasLayer   &ml  = *dynamic_cast<TVMeasLayer *>(At(ifr)); // get the last layer
            const TMaterial     &mat = ml.GetMaterial(isout);
            
            TKalMatrix Qms(sdim, sdim);
            if (IsMSOn()&& ito!=fridx ){
                
                CalcQms(mat, hel, fid, Qms, mass); // Qms for this step, using the fact that the material was found to be outgoing or incomming above, and the distance from the last layer
            }
            
            hel.MoveTo(xx, fid, &DF);         // move the helix to the present crossing point, DF will simply have its values overwritten so it could be explicitly set to unity here
            if (sdim == 6) DF(5, 5) = 1.;     // t0 stays the same
            F = DF * F;                       // update F
            TKalMatrix DFt  = TKalMatrix(TMatrixD::kTransposed, DF);
            
            Q = DF * (Q + Qms) * DFt;         // transport Q to the present crossing point
            
            if (IsDEDXOn() && ito!=fridx) {
                hel.PutInto(sv);                               // copy hel to sv
                // whether the helix is moving forwards or backwards is calculated using the sign of the charge and the sign of the deflection angle
                // Bool_t isfwd = ((cpa > 0 && df < 0) || (cpa <= 0 && df > 0)) ? kForward : kBackward;  // taken from TVMeasurmentLayer::GetEnergyLoss  not df = fid
                sv(2,0) += GetEnergyLoss(mat, hel, fid, mass); // correct for dE/dx, returns delta kappa i.e. the change in pt
                hel.SetTo(sv, hel.GetPivot());                 // save sv back to hel
            }
            ifr = ito; // for the next iteration set the "previous" layer to the current layer moved to
        } else {
            fid = fid_temp;
        }
    } // end of loop over surfaces
    
    //   // ---------------------------------------------------------------------
    //   //  Move pivot to crossing point with layer to move to
    //   // ---------------------------------------------------------------------
    //   dynamic_cast<const TVSurface *>(&ml_to)->CalcXingPointWith(hel, xx, fid);
    //   hel.MoveTo(xx, fid, &DF); // move pivot to expected hit, DF will simply have its values overwritten so it could be explicitly set to unity here
    //   F = DF * F;                          // update F accordingly
    
    x0 = hel.GetPivot();                 // local pivot corresponding to sv
    hel.PutInto(sv);                     // save updated hel to sv
    
    return 0;
    
}



//_________________________________________________________________________
// -----------------
//  Update
// -----------------
//    sorts meaurement layers according to layer's sorting policy
//    and puts index to layers from inside to outside.
//
void TKalDetCradle::Update()
{
    fDone = kTRUE;
    
    UnSort();   // unsort
    Sort();     // sort layers according to sorting policy
    
    TIter next(this);
    TVMeasLayer *mlp = 0;
    Int_t i = 0;
    
    while ((mlp = dynamic_cast<TVMeasLayer *>(next()))) {
        mlp->SetIndex(i++);
    }
    
}

//_________________________________________________________________________
// -----------------
//  GetEnergyLoss
// -----------------
//    returns energy loss.
//
Double_t TKalDetCradle::GetEnergyLoss(const TMaterial &mat,
                                      const TVTrack   &hel,
                                      Double_t   df,
                                      Double_t   mass) const
{
    Double_t cpa    = hel.GetKappa();
    Double_t tnl    = hel.GetTanLambda();
    Double_t tnl2   = tnl * tnl;
    Double_t tnl21  = 1. + tnl2;
    Double_t cslinv = TMath::Sqrt(tnl21);
    Double_t mom2   = tnl21 / (cpa * cpa);
    
    // -----------------------------------------
    // Bethe-Bloch eq. (Physical Review D P195.)
    // -----------------------------------------
    static const Double_t kK   = 0.307075e-3;     // [GeV*cm^2]
    static const Double_t kMe  = 0.510998902e-3;  // electron mass [GeV]
    
    Double_t dnsty = mat.GetDensity();		// density
    Double_t A     = mat.GetA();                 // atomic mass
    Double_t Z     = mat.GetZ();                 // atomic number
    //Double_t I    = Z * 1.e-8;			// mean excitation energy [GeV]
    //Double_t I    = (2.4 +Z) * 1.e-8;		// mean excitation energy [GeV]
    Double_t I    = (9.76 * Z + 58.8 * TMath::Power(Z, -0.19)) * 1.e-9;
    Double_t hwp  = 28.816 * TMath::Sqrt(dnsty * Z/A) * 1.e-9;
    Double_t bg2  = mom2 / (mass * mass);
    Double_t gm2  = 1. + bg2;
    Double_t meM  = kMe / mass;
    Double_t x    = log10(TMath::Sqrt(bg2));
    Double_t C0   = - (2. * log(I/hwp) + 1.);
    Double_t a    = -C0/27.;
    Double_t del;
    if (x >= 3.)            del = 4.606 * x + C0;
    else if (0.<=x && x<3.) del = 4.606 * x + C0 + a * TMath::Power(3.-x, 3.);
    else                    del = 0.;
    Double_t tmax = 2.*kMe*bg2 / (1. + meM*(2.*TMath::Sqrt(gm2) + meM));
    Double_t dedx = kK * Z/A * gm2/bg2 * (0.5*log(2.*kMe*bg2*tmax / (I*I))
                                          - bg2/gm2 - del);
    
    Double_t path = hel.IsInB()
    ? TMath::Abs(hel.GetRho()*df)*cslinv
    : TMath::Abs(df)*cslinv;
    
    //fg: switched from using cm to mm in KalTest - material (density) and energy still in GeV and cm
    path /= 10. ;
    
    Double_t edep = dedx * dnsty * path;
    
    
    Double_t cpaa = TMath::Sqrt(tnl21 / (mom2 + edep
                                         * (edep + 2. * TMath::Sqrt(mom2 + mass * mass))));
    Double_t dcpa = TMath::Abs(cpa) - cpaa;
    
    static const Bool_t kForward  = kTRUE;
    static const Bool_t kBackward = kFALSE;
    Bool_t isfwd = ((cpa > 0 && df < 0) || (cpa <= 0 && df > 0)) ? kForward : kBackward;
    return isfwd ? (cpa > 0 ? dcpa : -dcpa) : (cpa > 0 ? -dcpa : dcpa);
}

//_________________________________________________________________________
// -----------------
//  CalQms
// -----------------
//    calculates process noise matrix for multiple scattering with
//    thin layer approximation.
//
void TKalDetCradle::CalcQms(const TMaterial  &mat,
                            const TVTrack    &hel,
                            Double_t    df,
                            TKalMatrix &Qms,
                            Double_t    mass) const
{
    Double_t cpa    = hel.GetKappa();
    Double_t tnl    = hel.GetTanLambda();
    Double_t tnl2   = tnl * tnl;
    Double_t tnl21  = 1. + tnl2;
    Double_t cpatnl = cpa * tnl;
    Double_t cslinv = TMath::Sqrt(tnl21);
    Double_t mom    = TMath::Abs(1. / cpa) * cslinv;
    Double_t   beta = mom / TMath::Sqrt(mom * mom + mass * mass);
    
    Double_t x0inv = 1. / mat.GetRadLength();  // radiation length inverse
    
    // *Calculate sigma_ms0 =============================================
    static const Double_t kMS1  = 0.0136;
    static const Double_t kMS12 = kMS1 * kMS1;
    static const Double_t kMS2  = 0.038;
    
    Double_t path = hel.IsInB()
    ? TMath::Abs(hel.GetRho()*df)*cslinv
    : TMath::Abs(df)*cslinv;
    
    //fg: switched from using cm to mm in KalTest - material (density) and energy still in GeV and cm
    path /= 10. ;
    
    Double_t xl   = path * x0inv;
    // ------------------------------------------------------------------
    // Very Crude Treatment!!
    Double_t tmp = 1. + kMS2 * TMath::Log(TMath::Max(1.e-4, xl));
    tmp /= (mom * beta);
    Double_t sgms2 = kMS12 * xl * tmp * tmp;
    // ------------------------------------------------------------------
    
    Qms(1,1) = sgms2 * tnl21;
    Qms(2,2) = sgms2 * cpatnl * cpatnl;
    Qms(2,4) = sgms2 * cpatnl * tnl21;
    Qms(4,2) = sgms2 * cpatnl * tnl21;
    Qms(4,4) = sgms2 * tnl21  * tnl21;
}
