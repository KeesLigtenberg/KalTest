#ifndef __EXTPCDETECTOR__
#define __EXTPCDETECTOR__

#include "EXVKalDetector.h"

class EXTPCKalDetector : public EXVKalDetector {
public:
   EXTPCKalDetector(Int_t m = 100);
   ~EXTPCKalDetector();

   void ProcessHit(const TVector3    &xx,
                   const TVMeasLayer &ms,
                         TObjArray   &hits);

private:
   static Double_t fgVdrift;   // drift velocity

   ClassDef(EXTPCKalDetector,1)   // Sample hit class
};

#endif
