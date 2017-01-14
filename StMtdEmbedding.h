#ifndef STMTDEMBEDDING_HH
#define STMTDEMBEDDING_HH

#include <vector>
#include "StMaker.h"
#ifndef ST_NO_NAMESPACES
using std::vector;
#endif

class TH1F;
class TH2F;
class TH3F;
class TFile;
class TTree;
class THnSparse;
class TRandom3;

class StEvent;
class StTrack;
class StMcEvent;
class StMcTrack;
class StMtdHit;
class StMtdPidTraits;
class StPrimaryVertex;
class StAssociationMaker;
class StRefMultCorrVPDMBZDCNoVtx;

#include "TFile.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "StThreeVectorD.hh"
#include "StEnumerations.h"
#include "StMtdUtil/StMtdGeometry.h"
#include "StAssociationMaker/StAssociationMaker.h"

#if !defined(ST_NO_TEMPLATE_DEF_ARGS) || defined(__CINT__)
typedef vector<Int_t> IntVec;
typedef vector<TLorentzVector> LorentzVec;
typedef vector<Double_t> DoubleVec;
typedef vector<StThreeVectorD > PointVec;
#else
typedef vector<Int_t, allocator<Int_t>> IntVec;
typedef vector<TLorentzVector, allocator<TLorentzVector>> LorentzVec;
typedef vector<Double_t, allocator<Double_t>> DoubleVec;
typedef vector<StThreeVectorD, allocator<StThreeVectorD>> PointVec;
#endif

class StMtdEmbedding : public StMaker {
    public:
        StMtdEmbedding(const Char_t *name = "MtdEmbedding");
        ~StMtdEmbedding();

        Int_t    Init();
        //Int_t    InitRun(const Int_t runNumber);
        Int_t    Make();
        Int_t    Finish();
        void     setOutFileName(const TString name);

    protected:
        Int_t    processStEvent();
        Bool_t   passTrigger(StEvent *event);
        Int_t    getTofMult(StEvent *event);
        Int_t    selectPrimaryVertex(StEvent *event);
        Bool_t   passTrack(StTrack *track, StVertex *vtx) const;
        Bool_t   passEvent(StEvent *event);
        Int_t    gRefMult(StEvent *event, StThreeVectorF vtxPos);
        Bool_t   isGoodVertex(const Double_t vpd_vz, const Double_t tpc_vz, const Double_t tpc_vr);
        double   getDca(StGlobalTrack *globalTrack, StThreeVectorF vtxPos) const;
        Int_t    findMatchedRcTrack(StMcTrack *track);
        StMtdPidTraits *getMtdPidTraits(StTrack *track);
        Bool_t   isMcMtdHit(StMtdHit *hit);
        Bool_t   isMcMuon(StMcTrack *track);
        Double_t rotatePhi(Double_t phi) const;
        Double_t getNSigmaPi(StTrack *track);
        Bool_t   checkMuCandidate(Float_t pt, Float_t nSigmaPion, Float_t dy, Float_t dz, Float_t dtof);
        double   weightFunc(double pT);
        void     calPolarization(TLorentzVector iVec,TLorentzVector vec, TH3F*h,THnSparse *hn , THnSparse *hnCS, double w );
        bool     mtdTriggerEff(double pT);
        bool     mtdResponseEff(int bkl, int mod, double pT);
        void     bookHistos();
        void     printConfig();

        struct StMtdEmbedData
        {
            //---event level
            int    runID;
            int    tofMult;
            float mBField;
            float zdcRate;// * 1e-3
            float bbcRate;
            float tpcVx;
            float tpcVy;
            float tpcVz;
            float vpdVz;
            //---track level
            int    nMcTracks;
            int    nMcMuon;
            int    nMcJpsi;
            int    nRcJpsi;
            //mc muon tracks
            int    mcPkey[500];//MC track parent key
            int    mcGeantId[500];//muon Id = 5 or 6
            double mcpt[500];
            double mcphi[500];
            double mceta[500];
            int    mccharge[500];
            int    rcPkey[500];//RC track parent key
            int    rcNHitsFit[500];
            int    rcNHitsPoss[500];
            int    rcNHitsDedx[500];
            double rcDca[500];
            double rcpt[500];
            double rcphi[500];
            double rceta[500];
            double rcNSigmaPi[500];
            int    rcCharge[500];
            // both rc track match with mtd without track quality cut 
            int rcBackleg[500];
            int rcModule[500];
            double rcDz[500];
            double rcDy[500];
            double rcDtof[500];
            int passTrkCut[500];
            int passMuonCut[500];
        };

    private:
        StEvent            *mStEvent;
        StMcEvent          *mMcEvent;
        StAssociationMaker *mAssoMaker;                
        rcTrackMapType     *mRcTrackMap;                              
        mcTrackMapType     *mMcTrackMap;                          
        StPrimaryVertex    *mPriVtx;
        StRefMultCorrVPDMBZDCNoVtx      *mRefMultCorr;    
        StMtdGeometry      *mMtdGeom; 

        Bool_t             mDebug;
        Bool_t             mPrintConfig;
        Bool_t             mUseVpdVtx;

        // event level
        Int_t              mRunId; 
        Int_t              mTofMult;
        Double_t           mZdcRate;
        Double_t           mBbcRate;

        // vertex cuts
        Bool_t           mRequireVtxRanking;                               
        Bool_t           mUsePrimVtxClosestToVpd;  
        Bool_t           mRequireMtdHitForPrimVtx;
        Double_t         mMaxVtxZ;
        Double_t         mMaxVtxR;
        Double_t         mMaxDiffz;

        // track cuts
        Int_t            mTrackType;                                 // 0 - primary tracks; 1 - global tracks
        Double_t         mMinTrkPt;                                  // Minimum track pt
        Double_t         mMaxTrkPt;                                  // Maximum track pt
        Double_t         mMinTrkPhi;                                 // Minimum track phi
        Double_t         mMaxTrkPhi;                                 // Maximum track phi
        Double_t         mMinTrkEta;                                 // Minimum track eta
        Double_t         mMaxTrkEta;                                 // Maximum track eta
        Int_t            mMinNHitsFit;                               // Minimum number of hits used for track fit
        Int_t            mMinNHitsDedx;                              // Minimum number of hits used for de/dx
        Double_t         mMinFitHitsFaction;                         // Minimum fraction of NHitsFit/NHitsPoss
        Double_t         mMaxDca;                                    // Maximum track dc


        // muon PID
        Double_t         mMinNSigmaPiCut;
        Double_t         mMaxNSigmaPiCut;
        Double_t         mMuDYcmCut;
        Double_t         mMuDYSigCut;
        Double_t         mMuDZcmCut;
        Double_t         mMuDZSigCut;
        Double_t         mMuMindTofCut;
        Double_t         mMuMaxdTofCut;

        // maps
        map<Int_t, Int_t> mRcIndices;
        map<Int_t, Int_t> mMcIndices;

        TTree            *mOutTree;
        TFile            *mOutTreeFile;
        TString           mOutFileName;      // Name of the output file 
        StMtdEmbedData    mEmbedData;

        // list of histograms
        TH1F             *mhEventStat;
        TH2F             *mhTofMultVsZdcRate;
        TH2F             *mhTofMultVzBbcRate;
        TH2F             *mhTpcZvsVpdZ;
        TH1F             *mhDefaultVtxZ;
        TH1F             *mhVzVsOrigin;
        TH3F             *mhMcMPtCos;
        TH3F             *mhRcMPtCos;
        THnSparse        *mhMcTrkInfo;
        THnSparse        *mhMcMPtCosPhi;
        THnSparse        *mhMcMPtCosPhiCS;
        THnSparse        *mhRcMPtCosPhi;
        THnSparse        *mhRcMPtCosPhiCS;
        // for QA
        TH2F             *mhNmcJpsiVsNrcJpsi;
        TH1F             *mhRcJpsiParentPt; 
        THnSparse        *mhMcJpsiQA;
        THnSparse        *mhRcJpsiQA;
        THnSparse        *mhMcJpsiW;
        THnSparse        *mhRcJpsiW;

        ClassDef(StMtdEmbedding, 1)
};

inline void StMtdEmbedding::setOutFileName(const TString name) { mOutFileName = name; }
#endif
