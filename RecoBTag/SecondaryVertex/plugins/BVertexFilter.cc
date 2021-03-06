// -*- C++ -*-
//
// Package:    BVertexFilter
// Class:      BVertexFilter
// 
/**\class BVertexFilter BVertexFilter.cc DPGAnalysis/BVertexFilter/src/BVertexFilter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Andrea RIZZI
//         Created:  Mon Dec  7 18:02:10 CET 2009
// $Id: BVertexFilter.cc,v 1.4 2010/02/28 20:10:01 wmtan Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "RecoBTag/SecondaryVertex/interface/SecondaryVertex.h"
//
// class declaration
//
#include "RecoBTag/SecondaryVertex/interface/VertexFilter.h"
class BVertexFilter : public edm::EDFilter {
   public:
      explicit BVertexFilter(const edm::ParameterSet&);
      ~BVertexFilter();

   private:
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      edm::EDGetTokenT<reco::VertexCollection> token_primaryVertex;
      edm::EDGetTokenT<reco::VertexCollection> token_secondaryVertex;
      reco::VertexFilter                      svFilter;
      bool                                    useVertexKinematicAsJetAxis;
      int                                     minVertices;
};


BVertexFilter::BVertexFilter(const edm::ParameterSet& params):
      svFilter(params.getParameter<edm::ParameterSet>("vertexFilter")),
      useVertexKinematicAsJetAxis(params.getParameter<bool>("useVertexKinematicAsJetAxis")),
      minVertices(params.getParameter<int>("minVertices"))

{
	token_primaryVertex = consumes<reco::VertexCollection>(params.getParameter<edm::InputTag>("primaryVertices"));
	token_secondaryVertex = consumes<reco::VertexCollection>(params.getParameter<edm::InputTag>("secondaryVertices"));
        produces<reco::VertexCollection>();

}


BVertexFilter::~BVertexFilter()
{
}

bool
BVertexFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
 int count = 0; 
 edm::Handle<reco::VertexCollection> pvHandle; 
 iEvent.getByToken(token_primaryVertex, pvHandle);
 edm::Handle<reco::VertexCollection> svHandle; 
 iEvent.getByToken(token_secondaryVertex, svHandle);
 const reco::Vertex & primary = (*pvHandle.product())[0];
 const reco::VertexCollection & vertices = *svHandle.product();
 std::auto_ptr<reco::VertexCollection> recoVertices(new reco::VertexCollection);

 if(! primary.isFake()) 
 {
   for(reco::VertexCollection::const_iterator it=vertices.begin() ; it!=vertices.end() ; ++it)
    {
          GlobalVector axis(0,0,0);
          if(useVertexKinematicAsJetAxis) axis = GlobalVector(it->p4().X(),it->p4().Y(),it->p4().Z());
          if(svFilter(primary,reco::SecondaryVertex(primary,*it,axis,true),axis))  {
                count++;
                recoVertices->push_back(*it);
           }
   }
 }
   iEvent.put(recoVertices);

   return(count >= minVertices);
}


//define this as a plug-in
DEFINE_FWK_MODULE(BVertexFilter);
