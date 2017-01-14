#include "StMtdEmbedding.h"
#include "headers.h"

const Double_t muMass = 0.10566;

ClassImp(StMtdEmbedding)

    //_____________________________________________________________________________
StMtdEmbedding::StMtdEmbedding(const Char_t *name) : StMaker(name)
                                                     ,mOutFileName("")
{
    // default constructor
    mStEvent                 = NULL;
    mMcEvent                 = NULL;
    mAssoMaker               = NULL;
    mRcTrackMap              = NULL;
    mMcTrackMap              = NULL;
    mMtdGeom                 = NULL;
    mDebug                   = kFALSE;

    mPrintConfig             = kTRUE;

    mRunId                   = -1;
    mTofMult                 = -1;
    mZdcRate                 = 0;
    mBbcRate                 = 0;


    // vertex cuts
    mRequireVtxRanking       = kFALSE;
    mUsePrimVtxClosestToVpd  = kFALSE;
    mRequireMtdHitForPrimVtx = kFALSE;
    mMaxVtxZ                 = 100;
    mMaxVtxR                 = 2;
    mMaxDiffz                = 3;

    // track cuts
    mTrackType               = primary;  
    mMinTrkPt                = 1.;
    mMaxTrkPt                = 1e4;
    mMinNHitsFit             = 15;
    mMinNHitsDedx            = 10;
    mMinFitHitsFaction       = 0.52;
    mMaxDca                  = 3;
    mMinTrkPhi               = 0;
    mMaxTrkPhi               = 2*pi;
    mMinTrkEta               = -0.8;
    mMaxTrkEta               = 0.8;

    // muon PID
    mMinNSigmaPiCut          =-2;
    mMaxNSigmaPiCut          =3;
    mMuDYcmCut               =20.;
    mMuDYSigCut              =3.;
    mMuDZcmCut               =20.;
    mMuDZSigCut              =3.; 
    mMuMindTofCut            =-1e+4;
    mMuMaxdTofCut            =1.0;

    mOutTree                 = NULL;
    mOutTreeFile             = NULL;
    mUseVpdVtx               = false;

}

//_____________________________________________________________________________
StMtdEmbedding::~StMtdEmbedding()
{
}

//_____________________________________________________________________________
Int_t StMtdEmbedding::Init()
{
    if(mPrintConfig) printConfig();

    if(!mOutFileName.Length()){
        LOG_ERROR << "StMiniTreeMaker:: no output file specified for tree and histograms." << endm;
        return kStERR;
    }
    mOutTreeFile = new TFile(mOutFileName.Data(),"recreate");
    LOG_INFO << "StMiniTreeMaker:: create the output file to store the tree and histograms: " << mOutFileName.Data() << endm;

    bookHistos();
    return kStOK;
}

//_____________________________________________________________________________
Int_t StMtdEmbedding::Finish()
{
    mOutTreeFile->cd();
    mOutTreeFile->Write();

    // ThnSparse
    mhMcTrkInfo->Write();
    mhRcJpsiQA->Write();
    mhMcJpsiQA->Write();
    mhRcJpsiW->Write();
    mhMcJpsiW->Write();

    mhMcMPtCosPhi->Write();
    mhMcMPtCosPhiCS->Write();
    mhRcMPtCosPhi->Write();
    mhRcMPtCosPhiCS->Write();

    mOutTreeFile->Close();

    LOG_INFO << "StMiniTreeMaker::Finish() -> write out tree in " << mOutFileName.Data() << endm;

    return kStOK;
}

//_____________________________________________________________________________
Int_t StMtdEmbedding::Make()
{
    Int_t iret = -1;

    memset(&mEmbedData, 0, sizeof(mEmbedData));
    mRcIndices.clear();
    mMcIndices.clear();

    mStEvent=(StEvent *) GetInputDS("StEvent");
    if(mStEvent)
    {
        mMcEvent = (StMcEvent*) GetDataSet("StMcEvent");
        if(!mMcEvent)
        {
            LOG_WARN << "No McEvent is available " << endm;
            return kStWarn;
        }

        mAssoMaker = (StAssociationMaker*) GetMaker("StAssociationMaker");
        if (mAssoMaker)
        {
            mRcTrackMap = mAssoMaker->rcTrackMap();
            mMcTrackMap = mAssoMaker->mcTrackMap();
        }
        else
        {
            LOG_WARN << "No association map" << endm;
            return kStWarn;
        }
        iret = processStEvent();
    }
    mOutTree->Fill();
    return 0;
}

//_____________________________________________________________________________
Int_t StMtdEmbedding::processStEvent()
{
    mhEventStat->Fill(0.5);
    mRunId = mStEvent->runInfo()->runId();
    if(mDebug) cout<<"mRunId:"<<mRunId<<endl;
    mEmbedData.runID = mRunId;

    // dimuon triggered event
    if(!passTrigger(mStEvent)) return kStOK; 
    mhEventStat->Fill(1.5);
    if(!passEvent(mStEvent)) return kStOK;
    mhEventStat->Fill(2.5);

    // choose vertex 
    int rc_vtx_index = selectPrimaryVertex(mStEvent);
    if(mDebug) cout<<"rc_vtx_index: "<<rc_vtx_index<<endl;
    mPriVtx = mStEvent->primaryVertex(rc_vtx_index);
    StThreeVectorF mVerPosition = mPriVtx->position();
    double tpcVz = mVerPosition.z();
    if(mDebug) cout<<"choose TPC vertex Vz: "<<tpcVz<<endl;
    // all the particles are embedded to the default vertex
    mhDefaultVtxZ->Fill(tpcVz);
    mEmbedData.tpcVx = mVerPosition.x();
    mEmbedData.tpcVy = mVerPosition.y();
    mEmbedData.tpcVz = mVerPosition.z();

    // multiplicity
    mTofMult = getTofMult(mStEvent);
    if(mDebug) cout<<"mTofMult: "<<mTofMult<<endl;
    mZdcRate = mStEvent->runInfo()->zdcCoincidenceRate() * 1e-3;
    if(mDebug) cout<<"mZdcRate: "<<mZdcRate<<endl;
    mBbcRate = mStEvent->runInfo()->bbcCoincidenceRate() * 1e-3;
    if(mDebug) cout<<"mBbcRate: "<<mBbcRate<<endl;
    mhTofMultVsZdcRate->Fill(mZdcRate,mTofMult);
    mhTofMultVzBbcRate->Fill(mBbcRate,mTofMult);
    mEmbedData.tofMult = mTofMult;
    mEmbedData.bbcRate = mBbcRate;
    mEmbedData.zdcRate = mZdcRate;
    mEmbedData.mBField = mStEvent->runInfo()->magneticField();

    double vtx_x = -999;
    double vtx_y = -999;
    double vtx_z = -999;

    if(mUseVpdVtx)//false
    {
        vtx_x = mStEvent->primaryVertex()->position().x();
        vtx_y = mStEvent->primaryVertex()->position().y();
        StBTofCollection * tofCollection = mStEvent->btofCollection(); 
        if(tofCollection)
        {
            StBTofHeader *tofHeader = tofCollection->tofHeader();
            if(tofHeader)
                vtx_z = tofHeader->vpdVz();
        }
    }
    else// true & use default vertex (rc_vtx_index = 0)  
    {
        if(rc_vtx_index>-1 && passEvent(mStEvent))
        {
            StPrimaryVertex *vtx = mStEvent->primaryVertex(rc_vtx_index);
            vtx_x = vtx->position().x();
            vtx_y = vtx->position().y();
            vtx_z = vtx->position().z();
        }
    }
    if(mDebug) cout<<"TPC Vz: "<<vtx_z<<endl; 
    if(fabs(vtx_z)< mMaxVtxZ)
    {
        StThreeVectorF vtxPos(vtx_x,vtx_y,vtx_z);
        StBTofCollection * tofCollection = mStEvent->btofCollection(); 
        if(tofCollection)
        {
            StBTofHeader *tofHeader = tofCollection->tofHeader();
            if(tofHeader){
                double vpdVz = tofHeader->vpdVz();
                mEmbedData.vpdVz = vpdVz;
                mhTpcZvsVpdZ->Fill(vtx_z,vpdVz);
            }
        }
    }

    //--- build map between mc and reconstructed tracks ------//
    if(mDebug) cout<<"/****** build map betweet MC & RC ******/"<<endl;
    StSPtrVecMcTrack mctracks = mMcEvent->tracks();
    StSPtrVecTrackNode &trackNodes = mStEvent->trackNodes();
    Int_t nMcTracks = mctracks.size();
    for(Int_t i=0; i<nMcTracks; i++){
        StMcTrack *mctrack = dynamic_cast<StMcTrack *>(mctracks[i]);
        if(!mctrack) continue;
        Int_t rcIndex = findMatchedRcTrack(mctrack);
        mMcIndices[i] = rcIndex;
        if(rcIndex>=0){
            mRcIndices[rcIndex] = i;
            StTrack *rcTrack = trackNodes[rcIndex]->track(primary);
            StMtdPidTraits *mtdPid = getMtdPidTraits(rcTrack);
            if(mtdPid && isMcMtdHit(mtdPid->mtdHit())){
                pair<Double_t,Double_t> leading = mtdPid->mtdHit()->leadingEdgeTime();
                double tof_mc = (leading.first+leading.second)/2 ;   // smear the timing in MTD
                mtdPid->setTimeOfFlight(tof_mc);
            }
        }
    }
    mEmbedData.nMcTracks = nMcTracks;
    if(mDebug) cout<< "# of mc tracks: " << nMcTracks<<endl;

    // Reset the dy in the PidTraits for MC
    //Int_t nNodes = trackNodes.size();
    //const double cell_width = 3.8; //cm
    //for(Int_t i=0; i<nNodes; i++)
    //{
    //    StTrack *pTrack = trackNodes[i]->track(primary);
    //    if(!pTrack) continue;
    //    StMtdPidTraits *mtdPid = getMtdPidTraits(pTrack);
    //    if(!mtdPid) continue;
    //    StMtdHit* hit = mtdPid->mtdHit();
    //    if(!isMcMtdHit(hit)) continue;
    //    int backleg = hit->backleg();
    //    double dy = mtdPid->deltaY();
    //    if(backleg==8)  dy -= 3 * cell_width;
    //    if(backleg==24) dy += 2 * cell_width;
    //    mtdPid->setDeltaY(dy);
    //}

    // MC truth
    if(mDebug) cout<<"/****** MC Truth ******/"<<endl;
    Int_t nMcMuon = 0;
    Int_t nEmbedJpsi = 0;
    Int_t nRcJpsi = 0;
    for(Int_t i=0; i<nMcTracks; i++){
        StMcTrack *mctrack = dynamic_cast<StMcTrack *>(mctracks[i]);
        if(!mctrack) continue;
        if(!isMcMuon(mctrack)) continue;
        nMcMuon++;

        double pt = mctrack->pt();
        double eta = mctrack->pseudoRapidity();
        double phi = rotatePhi(mctrack->momentum().phi());
        double charge = mctrack->particleDefinition()->charge();
        double mc_fill[] = {pt, eta, phi, charge};
        mhMcTrkInfo->Fill(mc_fill);

        if(mDebug)cout<<"geantId: "<<mctrack->geantId()<<endl;
        if(mDebug)cout<<"parent geantId: "<<mctrack->parent()->geantId()<<endl;

        mEmbedData.mcPkey[i]     = mctrack->parent()->key();
        mEmbedData.mcGeantId[i]     = mctrack->geantId();
        mEmbedData.mcpt[i]     = mctrack->pt();
        mEmbedData.mcphi[i]    = mctrack->momentum().phi();
        mEmbedData.mceta[i]    = eta = mctrack->pseudoRapidity();
        mEmbedData.mccharge[i]    = mctrack->particleDefinition()->charge();

        const StThreeVectorF mcMom = mctrack->momentum();
        TLorentzVector mcMuon1(mcMom.x(),mcMom.y(),mcMom.z(),mctrack->energy());

        Int_t gq = mctrack->particleDefinition()->charge();
        for(Int_t j=i+1; j<nMcTracks; j++){
            StMcTrack *mctrack2 = dynamic_cast<StMcTrack *>(mctracks[j]);
            if(!mctrack2 || !isMcMuon(mctrack2)) continue;
            if(!mctrack->parent() || !mctrack2->parent() ||
                    mctrack->parent()->key() != mctrack2->parent()->key()) continue;
            const StThreeVectorF mcMom2 = mctrack2->momentum();
            TLorentzVector mcMuon2(mcMom2.x(),mcMom2.y(),mcMom2.z(),mctrack2->energy());
            Int_t gq2 = mctrack2->particleDefinition()->charge();
            if(gq*gq2>0) continue;
            // Mc Jpsi QA
            TLorentzVector mcParent = mcMuon1 + mcMuon2;
            //if(mcParent.Rapidity()<-0.5 || mcParent.Rapidity()>0.5) continue;
            nEmbedJpsi++;
            Double_t mcJpsi_m  = mcParent.M();
            Double_t mcJpsi_pt = mcParent.Pt();
            Double_t mcJpsi_y  = mcParent.Rapidity();
            Double_t mcJpsi_phi= mcParent.Phi();
            Double_t mcJpsiFill[] = {mcJpsi_m,mcJpsi_pt,mcJpsi_y,mcJpsi_phi};
            mhMcJpsiQA->Fill(mcJpsiFill);
            mhMcJpsiW->Fill(mcJpsiFill,weightFunc(mcJpsi_pt));

            if(gq>0) calPolarization(mcMuon1,mcParent,mhMcMPtCos,mhMcMPtCosPhi, mhMcMPtCosPhiCS,weightFunc(mcJpsi_pt));
            else calPolarization(mcMuon2,mcParent,mhMcMPtCos,mhMcMPtCosPhi, mhMcMPtCosPhiCS, weightFunc(mcJpsi_pt));
        }
    }

    if(mDebug) cout<<"/********** RC **************/"<<endl;
    for(Int_t i=0; i<nMcTracks; i++){
        StMcTrack *mctrack = dynamic_cast<StMcTrack *>(mctracks[i]);
        if(!mctrack) continue;
        if(!isMcMuon(mctrack)) continue;
        if(mMcIndices[i]<0) continue;

        StTrack *rcTrack1 = trackNodes[mMcIndices[i]]->track(primary);
        const StThreeVectorF imom = rcTrack1->geometry()->momentum();
        TLorentzVector muoni;
        muoni.SetXYZM(imom.x(),imom.y(),imom.z(),muMass);

        if(rcTrack1){
            mEmbedData.rcPkey[i] = rcTrack1 -> key();
            mEmbedData.rcNHitsFit[i] = rcTrack1 -> fitTraits().numberOfFitPoints(kTpcId);
            mEmbedData.rcNHitsPoss[i] = rcTrack1 -> numberOfPossiblePoints(kTpcId); 
            StTpcDedxPidAlgorithm pidAlgorithm;
            const StParticleDefinition *pd = rcTrack1 -> pidTraits(pidAlgorithm);
            if( pd && pidAlgorithm.traits() ){
                mEmbedData.rcNHitsDedx[i] = pidAlgorithm.traits()->numberOfPoints();
            }
            StGlobalTrack *globalTrack = dynamic_cast<StGlobalTrack*>(rcTrack1->node()->track(global));
            mEmbedData.rcDca[i] = getDca(globalTrack,mPriVtx->position());
            mEmbedData.rcpt[i] = muoni.Pt();
            mEmbedData.rcphi[i] = muoni.Phi();
            mEmbedData.rceta[i] = muoni.PseudoRapidity();
            mEmbedData.rcNSigmaPi[i] = getNSigmaPi(rcTrack1);
            mEmbedData.rcCharge[i] = rcTrack1->geometry()->charge();
        } 
        Int_t gq = rcTrack1->geometry()->charge();
        StMtdPidTraits *mtdPid1 = getMtdPidTraits(rcTrack1);
        if(!mtdPid1) continue;
        StMtdHit* hit1 = mtdPid1->mtdHit();if(!isMcMtdHit(hit1)) continue;
        int backleg1 = hit1->backleg();
        int module1 = hit1->module();
        mEmbedData.rcBackleg[i] = backleg1;
        mEmbedData.rcModule[i] = module1;

        double rcPt1 = muoni.Pt();
        double rcNSigmaPi1 = getNSigmaPi(rcTrack1);
        double dy1 = mtdPid1->deltaY();
        double dz1 = mtdPid1->deltaZ();
        double dtof1 = mtdPid1->timeOfFlight()-mtdPid1->expTimeOfFlight();
        mEmbedData.rcDy[i] = dy1;
        mEmbedData.rcDz[i] = dz1;
        mEmbedData.rcDtof[i] = dtof1;

        for(Int_t j=i+1; j<nMcTracks; j++){
            StMcTrack *mctrack2 = dynamic_cast<StMcTrack *>(mctracks[j]);
            if(!mctrack2 || !isMcMuon(mctrack2)) continue;
            if(mMcIndices[j]<0) continue;

            StTrack *rcTrack2 = trackNodes[mMcIndices[j]]->track(primary);
            const StThreeVectorF jmom = rcTrack2->geometry()->momentum();
            TLorentzVector muonj;
            muonj.SetXYZM(jmom.x(),jmom.y(),jmom.z(),muMass);
            //if(!passTrack(rcTrack1,mPriVtx) || !passTrack(rcTrack2,mPriVtx)) continue;

            Int_t gq2 = rcTrack2->geometry()->charge();

            StMtdPidTraits *mtdPid2 = getMtdPidTraits(rcTrack2);
            if(!mtdPid2) continue;

            StMtdHit* hit2 = mtdPid2->mtdHit();if(!isMcMtdHit(hit2)) continue;
            int backleg2 = hit2->backleg();
            int module2 = hit2->module();

            // check difference between default vertex and primary vertex first point
            StThreeVectorF origin1 = rcTrack1->geometry()->origin();
            StThreeVectorF diff1 = origin1 - mVerPosition;
            mhVzVsOrigin->Fill(diff1.mag());
            StThreeVectorF origin2 = rcTrack2->geometry()->origin();
            StThreeVectorF diff2 = origin2 - mVerPosition;
            mhVzVsOrigin->Fill(diff2.mag());

            // Mc parent pair and Jpsi pt
            if(!mctrack->parent() || !mctrack2->parent() ||
                    mctrack->parent()->key() != mctrack2->parent()->key()) continue;
            const StThreeVectorF mcMom = mctrack->momentum();
            TLorentzVector mcMuon1(mcMom.x(),mcMom.y(),mcMom.z(),mctrack->energy());
            const StThreeVectorF mcMom2 = mctrack2->momentum();
            TLorentzVector mcMuon2(mcMom2.x(),mcMom2.y(),mcMom2.z(),mctrack2->energy());
            TLorentzVector mcParent = mcMuon1 + mcMuon2;
            Double_t mcJpsi_pt = mcParent.Pt();
            mhRcJpsiParentPt->Fill(mcJpsi_pt);

            mEmbedData.passTrkCut[i]=-999;
            mEmbedData.passMuonCut[i]=-999;

            if(!passTrack(rcTrack1,mPriVtx) || !passTrack(rcTrack2,mPriVtx)) continue;
            mEmbedData.passTrkCut[i]=1;
            if(!checkMuCandidate(rcPt1,rcNSigmaPi1,dy1,dz1,dtof1)) continue;
            mEmbedData.passMuonCut[i]=1;
            if(!mtdTriggerEff(rcPt1)) continue;
            if(!mtdResponseEff(backleg1,module1,rcPt1)) continue;

            double rcPt2 = muonj.Pt();
            if(mDebug) cout<<"rcPt2: "<< rcPt2 <<endl;
            double rcNSigmaPi2 = getNSigmaPi(rcTrack2);
            double dy2 = mtdPid2->deltaY();
            double dz2 = mtdPid2->deltaZ();
            double dtof2 = mtdPid2->timeOfFlight()-mtdPid2->expTimeOfFlight();
            if(!checkMuCandidate(rcPt2,rcNSigmaPi2,dy2,dz2,dtof2)) continue;
            if(!mtdTriggerEff(rcPt2)) continue;
            if(!mtdResponseEff(backleg2,module2,rcPt2)) continue;

            // rc Jpsi QA
            TLorentzVector rcParent = muoni + muonj;
            //if(rcParent.Rapidity()<-0.5 || rcParent.Rapidity()>0.5) continue;
            nRcJpsi++;
            Double_t rcJpsi_m  = rcParent.M();
            Double_t rcJpsi_pt = rcParent.Pt();
            Double_t rcJpsi_y  = rcParent.Rapidity();
            Double_t rcJpsi_phi= rcParent.Phi();
            Double_t rcJpsiFill[] = {rcJpsi_m,rcJpsi_pt,rcJpsi_y,rcJpsi_phi};
            mhRcJpsiQA->Fill(rcJpsiFill);
            mhRcJpsiW->Fill(rcJpsiFill,weightFunc(mcJpsi_pt));

            if(gq>0) calPolarization(muoni,rcParent,mhRcMPtCos,mhRcMPtCosPhi, mhRcMPtCosPhiCS, weightFunc(mcJpsi_pt));
            else calPolarization(muonj,rcParent,mhRcMPtCos,mhRcMPtCosPhi, mhRcMPtCosPhiCS, weightFunc(mcJpsi_pt));
        }
    }

    if(mDebug) cout<< "# of mc muon tracks: " << nMcMuon<<endl;
    mEmbedData.nMcMuon = nMcMuon;
    if(mDebug) cout<<" nEmbedJpsi: "<<nEmbedJpsi<<endl;
    if(mDebug) cout<<" nRcJpsi: "<<nRcJpsi<<endl;
    mEmbedData.nMcJpsi = nEmbedJpsi;
    mEmbedData.nRcJpsi = nRcJpsi;
    mhNmcJpsiVsNrcJpsi ->Fill(nEmbedJpsi,nRcJpsi);
    return kStOK;
}

//_____________________________________________________________________________
Bool_t StMtdEmbedding::passTrigger(StEvent *event)
{
    // run14 AuAu 200 GeV TriggerId
    Int_t di_muon[3] = {470602, 480602, 490602};
    Int_t trigger = -1;
    for(Int_t i =0; i<3; i++){
        if(event->triggerIdCollection()->nominal()->isTrigger(di_muon[i]))
        {
            trigger = i;
            break;
        }
    }
    if(trigger==-1) return kFALSE;
    else return kTRUE;
}
//_____________________________________________________________________________
Int_t StMtdEmbedding::getTofMult(StEvent *event)
{
    int nTofMatch = 0; 
    StPrimaryVertex *vertex = event->primaryVertex();
    int nTracks = vertex->numberOfDaughters();
    for(int i=0; i<nTracks; i++) 
    {    
        StTrack *pTrack = vertex->daughter(i);
        if(pTrack->type()!=primary) continue;
        if(pTrack->fitTraits().numberOfFitPoints(kTpcId)<10) continue;
        StGlobalTrack *globalTrack = dynamic_cast<StGlobalTrack*>(pTrack->node()->track(global));
        StDcaGeometry* trDcaGeom = globalTrack->dcaGeometry();
        if(!trDcaGeom) continue;
        StPhysicalHelixD dcahh = trDcaGeom->helix();
        double dca = dcahh.distance(vertex->position(),kFALSE);
        if(dca>3) continue;
        if(pTrack->isBToFMatched()) nTofMatch++;
    }    
    return nTofMatch;
}
//_____________________________________________________________________________
Int_t StMtdEmbedding::selectPrimaryVertex(StEvent *event)
{
    Int_t nPrimVtx = event->numberOfPrimaryVertices();

    Int_t vtx_index = 0;//default vertex

    if(mUsePrimVtxClosestToVpd)//kTRUE_Run14 kFALSE_Run15  
    { 
        // Get VPD vz
        Double_t mVpdVz = 0.;
        StBTofCollection * tofCollection = event->btofCollection();
        if(!tofCollection) return vtx_index;
        StBTofHeader *tofHeader = tofCollection->tofHeader();
        if(!tofHeader)     return vtx_index;
        mVpdVz  = tofHeader->vpdVz();
        if(mDebug) cout<<"selectPrimaryVertex: mVpdVz:"<<mVpdVz<<endl;

        // find the primary vertex closest to VPD vz
        for(Int_t i=0; i<nPrimVtx; i++)
        { 
            StPrimaryVertex *vertex = event->primaryVertex(i);
            Double_t dz = vertex->position().z() - mVpdVz;
            if(abs(mVpdVz)<200. && fabs(dz)<3.)
            { 
                vtx_index = i;
                break;
            }
        }
        return vtx_index;
    }

    if(mRequireMtdHitForPrimVtx)//kFALSE
    {
        // Find the primary vertex matched with MTD hits
        vtx_index = -1;
        for(Int_t i=0; i<nPrimVtx; i++)
        {
            StSPtrVecPrimaryTrack &  tracks = event->primaryVertex(i)->daughters();
            Int_t nPrimTrk = tracks.size(); 
            Int_t nMtdTrk = 0;
            for(Int_t j=0; j<nPrimTrk; j++)
            {
                StTrack* pTrack = tracks[j];
                if(!pTrack || !passTrack(pTrack,event->primaryVertex(i))) continue;
                StSPtrVecTrackPidTraits& traits = pTrack->pidTraits();
                StMtdPidTraits* mtdpid = NULL;
                for(UInt_t it=0; it<traits.size(); it++)
                {
                    if (traits[it]->detector() == kMtdId)
                    {
                        mtdpid = dynamic_cast<StMtdPidTraits*>(traits[it]);
                        break; 
                    } 
                }   
                if(mtdpid)
                {
                    StMtdHit* hit = mtdpid->mtdHit();
                    if(!hit) continue;
                    nMtdTrk ++;
                } 
            }   
            if(nMtdTrk>=2)
            {
                vtx_index = i;
                break;
            } 
        }   
        return vtx_index;
    }

}

//_____________________________________________________________________________
Bool_t StMtdEmbedding::passEvent(StEvent *event)
{
    Int_t vtx_index = selectPrimaryVertex(event);
    if(vtx_index<0) return kFALSE;
    StPrimaryVertex* priVertex = event->primaryVertex(vtx_index);
    if(!priVertex) return kFALSE;
    StThreeVectorF verPos = priVertex->position();

    StBTofCollection * tofCollection = event->btofCollection();
    if(!tofCollection) return kFALSE;
    StBTofHeader *tofHeader = tofCollection->tofHeader();
    if(!tofHeader)     return kFALSE;
    if(mRequireVtxRanking && priVertex->ranking()<0) return kFALSE;//mRequireVtxRanking = kFALSE 
    if(! (TMath::Abs(verPos.x())>0 || TMath::Abs(verPos.y())>0 || TMath::Abs(verPos.z())>0) ) return kFALSE;
    double vr = TMath::Sqrt(verPos.x()*verPos.x() + verPos.y()*verPos.y());
    if(!isGoodVertex(tofHeader->vpdVz(), verPos.z(), vr)) return kFALSE;

    return kTRUE;
}

//_____________________________________________________________________________
Bool_t StMtdEmbedding::isGoodVertex(const Double_t vpd_vz, const Double_t tpc_vz, const Double_t tpc_vr)
{
    if(mMaxVtxZ<1e4 && TMath::Abs(tpc_vz)>=mMaxVtxZ)          return kFALSE;
    //if(mMaxDiffz<1e4 && TMath::Abs(tpc_vz-vpd_vz)>=mMaxDiffz) return kFALSE;
    if(mMaxVtxR<1e4 && tpc_vr>=mMaxVtxR)                      return kFALSE;
    return kTRUE;
}
//_____________________________________________________________________________
Bool_t StMtdEmbedding::passTrack(StTrack *track, StVertex *vtx) const 
{
    if(!track) {
        if(mDebug) cout<<"passTrack::!track"<<endl;
        return kFALSE;}
    StThreeVectorF mom = track->geometry()->momentum();
    Float_t pt = mom.perp();
    Float_t eta = mom.pseudoRapidity();
    Float_t phi = mom.phi();
    if(phi<0) phi += 2*pi;
    Int_t nHitsFit = track->fitTraits().numberOfFitPoints(kTpcId);
    Int_t nHitsPoss = track->numberOfPossiblePoints(kTpcId);

    if(pt < mMinTrkPt   || pt > mMaxTrkPt)       {if(mDebug) cout<<"passTrack::!Pt"<<endl;return kFALSE;}
    //if(eta < mMinTrkEta || eta > mMaxTrkEta)     {if(mDebug) cout<<"passTrack::!eta"<<endl;return kFALSE;}
    //if(phi < mMinTrkPhi || phi > mMaxTrkPhi)     {if(mDebug) cout<<"passTrack::!phi"<<endl;return kFALSE;}
    if(nHitsFit<mMinNHitsFit)                    {if(mDebug) cout<<"passTrack::!nHitFit"<<endl;return kFALSE;}
    if(1.*nHitsFit/nHitsPoss<mMinFitHitsFaction) {if(mDebug) cout<<"passTrack::!ratio"<<endl;return kFALSE;}

    StTpcDedxPidAlgorithm pidAlgorithm;
    const StParticleDefinition *pd = track->pidTraits(pidAlgorithm);
    if(!pd || !pidAlgorithm.traits()) return kFALSE;
    if(pidAlgorithm.traits()->numberOfPoints()<mMinNHitsDedx) {if(mDebug) cout<<"passTrack::!nHitdEdx"<<endl;return kFALSE;}

    StGlobalTrack *globalTrack = 0x0; 
    if(mTrackType==primary) globalTrack = dynamic_cast<StGlobalTrack*>(track->node()->track(global));
    else              globalTrack = dynamic_cast<StGlobalTrack*>(track);
    double dca = getDca(globalTrack,vtx->position());
    if(mMaxDca<1e4 && dca>mMaxDca) {if(mDebug) cout<<"passTrack::!dca"<<endl;return kFALSE;}
    return kTRUE;
}
//_____________________________________________________________________________
Int_t StMtdEmbedding::gRefMult(StEvent *event, StThreeVectorF vtxPos)
{
    int grefmult = 0;
    if(!event) return grefmult;

    StSPtrVecTrackNode &trackNodes = event->trackNodes();
    int nNodes = trackNodes.size();
    for(int i=0; i<nNodes; i++)
    {
        StGlobalTrack *pTrack = dynamic_cast<StGlobalTrack*>(trackNodes[i]->track(global));
        if(!pTrack) continue;
        if(pTrack->idTruth()<200) continue;

        StThreeVectorF momentum = pTrack->geometry()->momentum();
        if (fabs(momentum.pseudoRapidity()) <  0.5 &&
                fabs(getDca(pTrack,vtxPos)) < 3 &&
                pTrack->fitTraits().numberOfFitPoints(kTpcId) >= 10)
            grefmult++;
    }
    return grefmult;
}
//_____________________________________________________________________________
double StMtdEmbedding::getDca(StGlobalTrack *globalTrack, StThreeVectorF vtxPos) const
{
    if(!globalTrack) return 999;
    StDcaGeometry* trDcaGeom = globalTrack->dcaGeometry();
    if(!trDcaGeom) return 999;
    StPhysicalHelixD dcahh = trDcaGeom->helix();
    return dcahh.distance(vtxPos,kFALSE);
}
//_____________________________________________________________________________
Int_t StMtdEmbedding::findMatchedRcTrack(StMcTrack *track)
{
    pair<mcTrackMapIter,mcTrackMapIter> mcBounds = mMcTrackMap->equal_range(track);
    Int_t maxCommonHits = 0;
    const StGlobalTrack *rcCandTrack = 0;
    for(mcTrackMapIter mcMapIter = mcBounds.first; mcMapIter != mcBounds.second; mcMapIter ++)
    {
        StTrackPairInfo *pair = mcMapIter->second;
        const StGlobalTrack *rcTrack = pair->partnerTrack();
        Int_t commonHits = pair->commonTpcHits();
        if(commonHits > maxCommonHits)
        {
            maxCommonHits = commonHits;
            rcCandTrack = rcTrack;
        }
    }

    Int_t rcIndex = -1;
    if(maxCommonHits>=10)
    {
        StSPtrVecTrackNode &trackNodes = mStEvent->trackNodes();
        Int_t nNodes = trackNodes.size();
        for(Int_t i=0; i<nNodes; i++)
        {
            StTrack *rcTrack = trackNodes[i]->track(primary);
            if(!rcTrack) continue;
            if(rcTrack->key()==rcCandTrack->key())
            {
                rcIndex = i;
                break;
            }
        }
    }
    return rcIndex;
}
//_____________________________________________________________________________
Bool_t StMtdEmbedding::isMcMtdHit(StMtdHit *hit)
{
    if(mDebug) cout<<"isMcMtdHit"<<endl;
    return (hit && hit->idTruth()>0);
}
//_____________________________________________________________________________
StMtdPidTraits *StMtdEmbedding::getMtdPidTraits(StTrack *track)
{
    StSPtrVecTrackPidTraits &traits = track->pidTraits();
    StMtdPidTraits *mtdPid = 0;
    for(UInt_t itrait=0; itrait<traits.size(); itrait++)
    {
        if(traits[itrait]->detector() == kMtdId)
        {
            mtdPid = dynamic_cast<StMtdPidTraits*>(traits[itrait]);
            break;
        }
    }
    return mtdPid;
}
//_____________________________________________________________________________
Bool_t StMtdEmbedding::isMcMuon(StMcTrack *track)
{
    Int_t geantId = track->geantId();
    if(geantId==5 || geantId==6) return kTRUE;
    else return kFALSE;
}
//_____________________________________________________________________________
Double_t StMtdEmbedding::rotatePhi(Double_t phi) const
{
    Double_t outPhi = phi;
    while(outPhi<0) outPhi += 2*pi;
    while(outPhi>2*pi) outPhi -= 2*pi;
    return outPhi;
}

//_____________________________________________________________________________
Double_t StMtdEmbedding::getNSigmaPi(StTrack *track)
{
    if(!track) return kFALSE;
    double nSigmaPi = -999.;
    StTpcDedxPidAlgorithm pidAlgorithm;
    const StParticleDefinition *pd = track->pidTraits(pidAlgorithm);
    if(pd && pidAlgorithm.traits()){
        static StPionPlus* Pion = StPionPlus::instance();
        nSigmaPi = pidAlgorithm.numberOfSigma(Pion);
    }
    return nSigmaPi;
}
//_____________________________________________________________________________
Bool_t StMtdEmbedding::checkMuCandidate(Float_t pt, Float_t nSigmaPion, Float_t dy, Float_t dz, Float_t dtof)
{//no dy dz dtof cut
    bool passNSigmaPi = false;
    bool passMuDeltaY = true;
    bool passMuDeltaZ = true;
    bool passMuDeltaTof = true;
    if(mMinNSigmaPiCut < nSigmaPion &&  nSigmaPion < mMaxNSigmaPiCut)  passNSigmaPi = true;

    double marg   = 0.;
    if(pt>3.) marg = 0.5;
    // dy pt dependent cut
    double sigy   =  -17.6867 + 18.4528*exp(0.637142/pt);
    if( fabs(dy)  <= (mMuDYSigCut+marg)*sigy ) passMuDeltaY = true;
    //if( fabs(dy)  <= mMuDYcmCut ) passMuDeltaY = true;
    // dz pt dependent cut
    double sigz   =  -32.6793 + 32.6034*exp(0.444217/pt);
    if( fabs(dz)  <= (mMuDZSigCut+marg)*sigz ) passMuDeltaZ = true;
    //if( fabs(dz)  <= mMuDZcmCut ) passMuDeltaZ = true;
    // Pt dependent cut
    double tofcut = 0.0817528 + 0.0169419*exp(4.34897/pt);
    if( mMuMindTofCut <= dtof && dtof  <= mMuMaxdTofCut ) passMuDeltaTof = true;

    if( passNSigmaPi && passMuDeltaY && passMuDeltaZ && passMuDeltaTof ) return true;
    else return false;
}

//_____________________________________________________________________________
double StMtdEmbedding::weightFunc(double pT)
{
    pT = (3.4948*TMath::Power(TMath::Exp((-0.395305)*pT)+pT/2.91793,(-8.46161)))*4*TMath::Pi()*pT;
    return pT;
}
//_____________________________________________________________________________
void StMtdEmbedding::calPolarization(TLorentzVector iVec,TLorentzVector vec, TH3F*h,THnSparse *hn, THnSparse *hnCS, double w )
{
    //iVec:positron TLorentzVector  vec:JPSI LorentzVector
    TLorentzVector mu(iVec);//positron
    TLorentzVector Proton1(0.,0.,100.,100.),Proton2(0.,0.,-100.,100.);//pp200GeV

    TVector3 XXHX,YYHX,ZZHX;//hilicity frame XYZaxis
    ZZHX = vec.Vect();
    YYHX = vec.Vect().Cross(Proton1.Vect());
    XXHX = YYHX.Cross(ZZHX); 

    TVector3 jpsi = vec.BoostVector();

    mu.Boost(-1*jpsi);

    Float_t theta = mu.Angle(ZZHX);
    Float_t phi = TMath::ATan2((mu.Vect().Dot(YYHX.Unit())),(mu.Vect().Dot(XXHX.Unit())));
    Float_t cosTheta = TMath::Cos(theta);
    h->Fill(vec.M(),vec.Pt(),cosTheta);
    Double_t fill[]={vec.M(),vec.Pt(),cosTheta,phi};
    hn->Sumw2();
    hn->Fill(fill,w);

    Proton1.Boost(-1*jpsi);
    Proton2.Boost(-1*jpsi);
    TVector3 XX,YY,ZZ;//Collins-Soper frame
    ZZ = Proton1.Vect()*(1/Proton1.Vect().Mag())-Proton2.Vect()*(1/Proton2.Vect().Mag());
    YY = Proton1.Vect().Cross(Proton2.Vect());
    XX = YY.Cross(ZZ);
    Float_t thetaCS = mu.Angle(ZZ);
    Float_t phiCS = TMath::ATan2((mu.Vect().Dot(YY.Unit())),(mu.Vect().Dot(XX.Unit())));
    Float_t cosThetaCS = TMath::Cos(thetaCS);
    Double_t fillCS[]={vec.M(),vec.Pt(),cosThetaCS,phiCS};
    hnCS->Sumw2();
    hnCS->Fill(fillCS,w);

}

//_____________________________________________________________________________
bool StMtdEmbedding::mtdTriggerEff(double pT)
{
    TFile *f = TFile::Open("~/Run15/Polarizatoin/EmbeddingTest/StRoot/StMtdEmbedding/Run13MTDTriggerEfficiency.root","read");
    TF1 *fEff = (TF1*) f -> Get("fittrigerr");
    double pEff = fEff -> Eval(pT);
    double p = gRandom -> Uniform(1);
    f->Close();
    if(p<pEff) return true;
    else {
        if(mDebug) cout<<"mtdTriggerEff: random: "<<p<<" ;pt: "<<pT<<" ;eff: "<<pEff<<endl; 
        return false;
    }
}

//_____________________________________________________________________________
bool StMtdEmbedding::mtdResponseEff(int bkl, int mod, double pT)
{
    TFile *f = TFile::Open("~/Run15/Polarizatoin/EmbeddingTest/StRoot/StMtdEmbedding/Run15ResponseEffViaPtTemplate.root","read");
    TF1 *fResLow = (TF1*) f -> Get(Form("fPtMtdEffBkl%d_Mod%d",bkl-1,mod-1));
    TF1 *fResHigh = (TF1*) f -> Get(Form("fSclPtTmpBkl%d_Mod%d",bkl-1,mod-1));
    double pEff=1;
    if(pT<6)pEff = fResLow->Eval(pT);
    else pEff = fResHigh->Eval(pT);
    double p = gRandom -> Uniform(1);
    f->Close();
    if(p<pEff) return true;
    else {
        if(mDebug) cout<<"mtdRespose: random: "<<p<<" ;pt: "<<pT<<" ;eff: "<<pEff<<endl; 
        return false;
    }
}
//_____________________________________________________________________________
void StMtdEmbedding::bookHistos()
{
    //mOutTreeFile = new TFile("test.root","recreate");
    mOutTree = new TTree("EmbedTree","Embedding tree");
    mOutTree->SetAutoSave(100000);
    //--event level
    mOutTree->Branch("runID",   &mEmbedData.runID,   "runID/I");
    mOutTree->Branch("tofMult",   &mEmbedData.tofMult,   "tofMult/I");
    mOutTree->Branch("mBField",   &mEmbedData.mBField,   "mBField/F");
    mOutTree->Branch("bbcRate",   &mEmbedData.bbcRate,   "bbcRate/F");
    mOutTree->Branch("zdcRate",   &mEmbedData.zdcRate,   "zdcRate/F");
    mOutTree->Branch("tpcVx",   &mEmbedData.tpcVx,   "tpcVx/F");
    mOutTree->Branch("tpcVy",   &mEmbedData.tpcVy,   "tpcVy/F");
    mOutTree->Branch("tpcVz",   &mEmbedData.tpcVz,   "tpcVz/F");
    mOutTree->Branch("vpdVz",   &mEmbedData.vpdVz,   "vpdVz/F");
    //--track level
    mOutTree->Branch("nMcTracks",   &mEmbedData.nMcTracks,   "nMcTracks/I");
    mOutTree->Branch("nMcMuon",   &mEmbedData.nMcMuon,   "nMcMuon/I");
    mOutTree->Branch("nMcJpsi",   &mEmbedData.nMcJpsi,   "nMcJpsi/I");
    mOutTree->Branch("nRcJpsi",   &mEmbedData.nRcJpsi,   "nRcJpsi/I");
    mOutTree->Branch("mcPkey",   mEmbedData.mcPkey,   "mcPkey[nMcTracks]/I");
    mOutTree->Branch("mcGeantId",   mEmbedData.mcGeantId,   "mcGeantId[nMcTracks]/I");
    mOutTree->Branch("mcpt",      mEmbedData.mcpt,      "mcpt[nMcTracks]/D");
    mOutTree->Branch("mcphi",     mEmbedData.mcphi,     "mcphi[nMcTracks]/D");
    mOutTree->Branch("mceta",     mEmbedData.mceta,     "mceta[nMcTracks]/D");
    mOutTree->Branch("mccharge",     mEmbedData.mccharge,     "mccharge[nMcTracks]/I");
    mOutTree->Branch("rcPkey",  mEmbedData.rcPkey,   "rcPkey[nMcTracks]/I");
    mOutTree->Branch("rcNHitsFit",   mEmbedData.rcNHitsFit,   "rcNHitsFit[nMcTracks]/I");
    mOutTree->Branch("rcNHitsPoss",   mEmbedData.rcNHitsPoss,   "rcNHitsPoss[nMcTracks]/I");
    mOutTree->Branch("rcNHitsDedx",   mEmbedData.rcNHitsDedx,   "rcNHitsDedx[nMcTracks]/I");
    mOutTree->Branch("rcDca",   mEmbedData.rcDca,   "rcDca[nMcTracks]/D");
    mOutTree->Branch("rcpt",      mEmbedData.rcpt,      "rcpt[nMcTracks]/D");
    mOutTree->Branch("rcphi",     mEmbedData.rcphi,     "rcphi[nMcTracks]/D");
    mOutTree->Branch("rceta",     mEmbedData.rceta,     "rceta[nMcTracks]/D");
    mOutTree->Branch("rcNSigmaPi",  mEmbedData.rcNSigmaPi,  "rcNSigmaPi[nMcTracks]/D");
    mOutTree->Branch("rcCharge",     mEmbedData.rcCharge,     "rcCharge[nMcTracks]/I");
    mOutTree->Branch("rcBackleg",        mEmbedData.rcBackleg,        "rcBackleg[nMcTracks]/I");
    mOutTree->Branch("rcModule",        mEmbedData.rcModule,        "rcModule[nMcTracks]/I");
    mOutTree->Branch("rcDz",        mEmbedData.rcDz,        "rcDz[nMcTracks]/D");
    mOutTree->Branch("rcDy",        mEmbedData.rcDy,        "rcDy[nMcTracks]/D");
    mOutTree->Branch("rcDtof",      mEmbedData.rcDtof,      "rcDtof[nMcTracks]/D");
    mOutTree->Branch("passTrkCut",      mEmbedData.passTrkCut,      "passTrkCut[nMcTracks]/I");
    mOutTree->Branch("passMuonCut",      mEmbedData.passMuonCut,      "passMuonCut[nMcTracks]/I");


    // event histograms
    mhEventStat = new TH1F("hEventStat","Event statistics",10,0,10);
    mhEventStat->GetXaxis()->SetBinLabel(1,"All events");
    mhEventStat->GetXaxis()->SetBinLabel(2,"Good dimuon trigger");
    mhEventStat->GetXaxis()->SetBinLabel(3,"Good vertex");

    mhTofMultVsZdcRate = new TH2F(Form("mhTofMultVsZdcRate"),Form("dimuon: mTofMult vs zdc rate;zdc rate (kHz);mTofMult"),1000,0,1000,50,0,50);
    mhTofMultVzBbcRate = new TH2F(Form("mhTofMultVzBbcRate"),Form("dimuon: mTofMult vs bbc rate;bbc rate (kHz);mTofMult"),1000,0,10000,50,0,50);
    mhTpcZvsVpdZ  = new TH2F(Form("mhTpcZvsVpdZ"),Form("z distribution of default primary vertex vs Vpd vertex;TPC vz (cm);VPD vz (cm)"),300,-150,150,300,-150,150);
    mhDefaultVtxZ  = new TH1F(Form("mhDefaultVtxZ"),Form("z distribution of default primary vertex;vz (cm)"),300,-150,150);
    mhVzVsOrigin  = new TH1F("mhVzVsOrigin","difference between Tpc default vertex and primary track origin point",300,-150,150);

    mhNmcJpsiVsNrcJpsi = new TH2F("mhNmcJpsiVsNrcJpsi","MC input VS RC J/#psi per event;MC input # of J/#psi;Reconstructed # of J/#psi",50,0,50,50,0,50);

    mhRcJpsiParentPt = new TH1F("mhRcJpsiParentPt","used in RC loop",400,0,20);

    const int nTrkPtBin = 200;
    const double lowTrkPtBin = 0, hiTrkPtBin = 20.;
    const int dimTrkInfo = 4;
    const int nBinsTrkInfo[dimTrkInfo] = {nTrkPtBin, 40, 360, 3};
    const double lowBinTrkInfo[dimTrkInfo] = {lowTrkPtBin, -2, 0, -1};
    const double upBinTrkInfo[dimTrkInfo] = {hiTrkPtBin, 2, 2*pi, 2};
    // MC truth
    mhMcTrkInfo = new THnSparseF(Form("mhMcTrkInfo"),Form("dimuon: p_{T} vs #eta vs #varphi vs charge  of MC tracks;p_{T}^{mc} (GeV/c);#eta_{mc};#varphi_{mc}"),dimTrkInfo,nBinsTrkInfo,lowBinTrkInfo,upBinTrkInfo);

    // for QA
    const Int_t dim = 4;
    const Int_t nBins[dim]={200,400,40,96};
    const Double_t lowBins[dim]={2,0,-1,-TMath::Pi()};
    const Double_t highBins[dim]={4,20,1,TMath::Pi()};
    mhMcJpsiQA = new THnSparseF("mhMcJpsiQA","MC input: invariant mass distribution;M_{#mu#mu} (GeV/c^{2});p_{T} (GeV/c); y;#varphi",dim,nBins,lowBins,highBins);
    mhRcJpsiQA = new THnSparseF("mhRcJpsiQA","RC: invariant mass distribution;M_{#mu#mu} (GeV/c^{2});p_{T} (GeV/c); y;#varphi",dim,nBins,lowBins,highBins);
    mhMcJpsiW = new THnSparseF("mhMcJpsiW","MC input: after weight ;M_{#mu#mu} (GeV/c^{2});p_{T} (GeV/c); y;#varphi",dim,nBins,lowBins,highBins);
    mhRcJpsiW = new THnSparseF("mhRcJpsiW","RC: after weight ;M_{#mu#mu} (GeV/c^{2});p_{T} (GeV/c); y;#varphi",dim,nBins,lowBins,highBins);

    // for calPolarization 
    mhMcMPtCos = new TH3F("mhMcMPtCos","Mc unLikeSignPair;m_{#mu^{+}#mu^{-}};p_{T} (GeV/c);cos#theta",200,2.,4., 400,0.,20.,20,-1.,1.);
    mhRcMPtCos = new TH3F("mhRcMPtCos","Rc unLikeSignPair;m_{#mu^{+}#mu^{-}};p_{T} (GeV/c);cos#theta",200,2.,4., 400,0.,20.,20,-1.,1.);

    mhMcMPtCosPhi = new THnSparseF("mhMcMPtCosPhi","Mc HF UnlikeSignPair;m_{#mu#mu} (GeV/c^{2});p_{T} (GeV/c);cos#theta;#phi",dim,nBins,lowBins,highBins);
    mhMcMPtCosPhiCS = new THnSparseF("mhMcMPtCosPhiCS","Mc CS UnlikeSignPair;m_{#mu#mu} (GeV/c^{2});p_{T} (GeV/c);cos#theta;#phi",dim,nBins,lowBins,highBins);
    mhRcMPtCosPhi = new THnSparseF("mhRcMPtCosPhi","Rc HF UnlikeSignPair;m_{#mu#mu} (GeV/c^{2});p_{T} (GeV/c);cos#theta;#phi",dim,nBins,lowBins,highBins);
    mhRcMPtCosPhiCS = new THnSparseF("mhRcMPtCosPhiCS","Rc CS UnlikeSignPair;m_{#mu#mu} (GeV/c^{2});p_{T} (GeV/c);cos#theta;#phi",dim,nBins,lowBins,highBins);

}

//_____________________________________________________________________________
void StMtdEmbedding::printConfig()
{
    printf("\n=== Configuration for StMtdEmbedding ===\n");
    if(mRequireVtxRanking)       printf("Require vertex rannking\n");
    if(mUsePrimVtxClosestToVpd)  printf("Use the primary vtx closest to VPD\n");
    if(mRequireMtdHitForPrimVtx) printf("Use the primary vtx with MTD hits\n");
    if(mTrackType==primary) printf("Use primary tracks\n");
    if(mTrackType==global) printf("Use global tracks\n");
    printf("Maximum vertex r: %1.1f cm\n",mMaxVtxR);
    printf("Maximum vertex z: %1.1f cm\n",mMaxVtxZ);
    printf("|TPC-VPD| < %1.1f cm\n",mMaxDiffz);
    printf("Track pt  range: [%1.2f, %1.2f]\n",mMinTrkPt,mMaxTrkPt);
    printf("Track phi range: [%1.2f, %1.2f]\n",mMinTrkPhi,mMaxTrkPhi);
    printf("Track eta range: [%1.2f, %1.2f]\n",mMinTrkEta,mMaxTrkEta);
    printf("Minimum number of fit hits: %d\n",mMinNHitsFit);
    printf("Minimum number of dedx hits: %d\n",mMinNHitsDedx);
    printf("Minimum fraction of fit hits: %4.2f\n",mMinFitHitsFaction);
    printf("Maximum dca: %1.1f cm\n",mMaxDca);
    printf("=======================================\n\n");
}

