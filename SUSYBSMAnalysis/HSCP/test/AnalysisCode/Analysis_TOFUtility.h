// Original Author:  Loic Quertenmont
 

#include <exception>
#include <vector>
#include <unordered_map>

class moduleGeom{
   public:
      float TP; TVector3 pos; TVector3 width; TVector3 length; TVector3 thick;
      static std::unordered_map<unsigned int, moduleGeom*> static_geomMap;
   public :
        moduleGeom(){}
        ~moduleGeom(){}

        TVector3 toGlobal(TVector3 local ){ return (pos + local.x()*width.Unit() + local.y()*length.Unit() + local.z()*thick.Unit());}
        TVector3 toLocal (TVector3 global){ TVector3 o = global-pos;  return TVector3((o*width.Unit()), (o*length.Unit()), (o*thick.Unit()));}

        static moduleGeom* get(unsigned int detId){return static_geomMap[detId]; }

        static void loadGeometry(string path){
            moduleGeom::static_geomMap.clear(); //reset the geometry map

            TChain* t = new TChain("GeomDumper/geom");
            t->Add(path.c_str());

            unsigned int rawId;                 t->SetBranchAddress("rawId", &rawId);
            float trapezeParam;                 t->SetBranchAddress("trapezeParam", &trapezeParam);
            TVector3* posV    = new TVector3(); t->SetBranchAddress("pos",    &posV);
            TVector3* widthV  = new TVector3(); t->SetBranchAddress("width",  &widthV);
            TVector3* lengthV = new TVector3(); t->SetBranchAddress("length", &lengthV);
            TVector3* thickV  = new TVector3(); t->SetBranchAddress("thick",  &thickV);

            for (unsigned int ientry = 0; ientry < t->GetEntries(); ientry++) {
                t->GetEntry(ientry);
                
                moduleGeom* mod = new moduleGeom();
                mod->TP     = trapezeParam;
                mod->pos    = *posV;
                mod->width  = *widthV;
                mod->length = *lengthV;
                mod->thick  = *thickV;
                static_geomMap[rawId] = mod;
            }
            delete t;
        }



        static string getStationName(unsigned int detId){
            char stationName[255]; 
            DetId    geomDetId(detId);
            if(geomDetId.subdetId()==1){
               DTChamberId dtId(detId);
               sprintf(stationName,"DT%i%i",dtId.wheel(),dtId.station() );
            }else if(geomDetId.subdetId()==2){
               CSCDetId cscId(detId);
               CSCDetId cscchamberId = cscId.chamberId();
               sprintf(stationName,"ME%c%i%i",cscchamberId.zendcap()>0?'+':'-',cscchamberId.station(), cscchamberId.ring() );
            }
            return stationName;
        }

        static string getChamberName(unsigned int detId){
            char chamberName[255]; 
            DetId    geomDetId(detId);
            if(geomDetId.subdetId()==1){
               DTChamberId dtId(detId);
               sprintf(chamberName,"%s_%02i" ,getStationName(detId).c_str(), dtId.sector() );
            }else if(geomDetId.subdetId()==2){
               CSCDetId cscId(detId);
               CSCDetId cscchamberId = cscId.chamberId();
               sprintf(chamberName,"%s_%02i" ,getStationName(detId).c_str(), cscchamberId.chamber() );
            }
            return chamberName;
        }

};
std::unordered_map<unsigned int, moduleGeom*> moduleGeom::static_geomMap; //need to define this here to reference the object


class muonTimingCalculator{
   private:
      std::unordered_map<string, float>* t0OffsetMap; 
      std::map<unsigned int, std::unordered_map<string, float> > t0OffsetMapPerRuns;

   public:
      struct TimeMeasurementCSC {float distIP;float timeCorr; int station; float weightVertex; float weightInvbeta; };
      struct TimeMeasurementDT  {float distIP;float timeCorr; int station; bool isLeft; bool isPhi; float posInLayer; DetId driftCell;};
      struct TimeMeasurement    {double distIP; double localt0; double weightVertex; double weightInvBeta;
         TimeMeasurement(double distIP_, double localt0_, double weightVertex_, double weightInvBeta_){distIP = distIP_; localt0=localt0_; weightVertex=weightVertex_; weightInvBeta=weightInvBeta_;}
      };
      
      vector<const CSCSegment*> cscSegs;
      vector<const DTRecSegment4D*> dtSegs;
      vector<TimeMeasurement> tmSeq;

      //////////Global variables for Segment matcher
      double DTradius;
      bool   TightMatchDT;
      bool   cscTightMatch;

      /////////Global Variables for the CSC Time Extractor
      double thePruneCut_;
      double theStripTimeOffset_;
      double theWireTimeOffset_;
      double theStripError_;
      double theWireError_;
      double UseWireTime_;
      double UseStripTime_;

      /////////Global Variables for the DT Time Extractor
      int    theHitsMin_;
      double thePruneCutDT_;
      double theTimeOffset_;
      double theError_;
      bool   useSegmentT0_;
      bool   doWireCorr_;
      bool   dropTheta_;
      bool   requireBothProjections_;

   public:
      muonTimingCalculator(){
         //////////Global variables for Segment matcher
         DTradius      = 0.01;
         TightMatchDT  = false;
         cscTightMatch = true;

         /////////Global Variables for the CSC Time Extractor
         thePruneCut_ = 9.0;
         theStripTimeOffset_ = 0.0;
         theWireTimeOffset_ = 0.0;
         theStripError_ = 7.0;
         theWireError_ = 8.6;
         UseWireTime_ = true;
         UseStripTime_ = true;

         /////////Global Variables for the DT Time Extractor
         theHitsMin_ = 3;
         thePruneCutDT_ = 10000.0;
         theTimeOffset_ = 0.0;
         theError_ = 6.0;
         useSegmentT0_ = false;
         doWireCorr_ = true;
         dropTheta_ = true;
         requireBothProjections_ = false;
      }
     ~muonTimingCalculator(){}

      ////////////////////////////////////////////
      //all code related to t0 corrections
      ////////////////////////////////////////////
     
      void loadTimeOffset(string path){
         t0OffsetMap = NULL;
         t0OffsetMapPerRuns.clear();

         std::vector<string> chambers;

         FILE* pFile = fopen(path.c_str(), "r");
         char line[16384];
         while(fgets(line, 16384, pFile)){
            unsigned int run; char chamber[128];  float correction;
            char* pch=strtok(line,",");
            if(string(pch).find("runs")!=std::string::npos){  //get the list of runs
               while((pch=strtok(NULL,","))){sscanf(pch, "%u", &run); t0OffsetMapPerRuns[run] = std::unordered_map<string, float>(); }         
            }else if(string(pch).find("chambers")!=std::string::npos){  //get the list of chambers
               while((pch=strtok(NULL,","))){sscanf(pch, "%s", chamber); chambers.push_back(chamber);}         
            }else if(string(pch).find("run")!=std::string::npos){
               sscanf(pch, "run %u", &run);
               std::unordered_map<string, float>& t0OffsetMap =  t0OffsetMapPerRuns[run];
               int Index=0;
               while((pch=strtok(NULL,","))){sscanf(pch, "%f", &correction); t0OffsetMap[chambers[Index]] = correction; Index++;}
            }else{
               printf("Unknown type of line");
            }
         }
         fclose(pFile);
      }

      void setRun(unsigned int currentRun){
         std::map<unsigned int, std::unordered_map<string, float> >::iterator it, itPrev=t0OffsetMapPerRuns.begin();
         for(it=t0OffsetMapPerRuns.begin(); it!=t0OffsetMapPerRuns.end(); it++){
            if(it->first>currentRun){t0OffsetMap = &(itPrev->second); return;}//runs are ordered, so the previous iterator correspond to our run
            itPrev=it;
         }
         t0OffsetMap = &(it->second); //just in case we go beyond the list of run for which we have a correciton
      }

      double t0Offset(unsigned int detId){ return (*t0OffsetMap)[moduleGeom::getStationName(detId)]; } //stupid to save the information in a map of string key
      double t0OffsetChamber(unsigned int detId){  return (*t0OffsetMap)[moduleGeom::getChamberName(detId)]; } //stupid to save the information in a map of string key


      ////////////////////////////////////////////
      //all code related to muon segment matching
      ////////////////////////////////////////////

         //Fwlite version of the code found in https://raw.githubusercontent.com/cms-sw/cmssw/CMSSW_7_6_X/RecoMuon/TrackingTools/src/MuonSegmentMatcher.cc
         vector<const CSCSegment*>& matchCSC(const reco::Track& muon, const CSCSegmentCollection& CSCSegmentColl)
         {
           cscSegs.clear();

           double matchRatioCSC=0;
           double CSCXCut = 0.001;
           double CSCYCut = 0.001;
           double countMuonCSCHits = 0;

           for(CSCSegmentCollection::const_iterator segmentCSC = CSCSegmentColl.begin(); segmentCSC != CSCSegmentColl.end(); segmentCSC++) {
             double CSCcountAgreeingHits=0;
             if ( !segmentCSC->isValid()) continue; 

             const vector<CSCRecHit2D>& CSCRechits2D = segmentCSC->specificRecHits();
             countMuonCSCHits = 0;
             CSCDetId myChamber((*segmentCSC).geographicalId().rawId());

             bool segments = false;

             for(trackingRecHit_iterator hitC = muon.recHitsBegin(); hitC != muon.recHitsEnd(); ++hitC) {
               if (!(*hitC)->isValid()) continue; 
               if ( (*hitC)->geographicalId().det() != DetId::Muon ) continue; 
               if ( (*hitC)->geographicalId().subdetId() != MuonSubdetId::CSC ) continue;
               if (!(*hitC)->isValid()) continue;
               if ( (*hitC)->recHits().size()>1) segments = true;

               //DETECTOR CONSTRUCTION
               DetId id = (*hitC)->geographicalId();
               CSCDetId cscDetIdHit(id.rawId());

               if (segments) {
                 if(!(myChamber.rawId()==cscDetIdHit.rawId())) continue; 

                 // and compare the local positions
                 LocalPoint positionLocalCSC = (*hitC)->localPosition();
                 LocalPoint segLocalCSC = segmentCSC->localPosition();
                 if ((fabs(positionLocalCSC.x()-segLocalCSC.x())<CSCXCut) && 
                     (fabs(positionLocalCSC.y()-segLocalCSC.y())<CSCYCut)) 
                   cscSegs.push_back(&(*segmentCSC)); 
                 continue;
               }

               if(!(cscDetIdHit.ring()==myChamber.ring())) continue;
               if(!(cscDetIdHit.station()==myChamber.station())) continue;
               if(!(cscDetIdHit.endcap()==myChamber.endcap())) continue;
               if(!(cscDetIdHit.chamber()==myChamber.chamber())) continue;

               countMuonCSCHits++;

               LocalPoint positionLocalCSC = (*hitC)->localPosition();
                 
               for (vector<CSCRecHit2D>::const_iterator hiti=CSCRechits2D.begin(); hiti!=CSCRechits2D.end(); hiti++) {

                 if ( !hiti->isValid()) continue; 
                 CSCDetId cscDetId((hiti->geographicalId()).rawId());
                         
                 if ((*hitC)->geographicalId().rawId()!=(hiti->geographicalId()).rawId()) continue;

                 LocalPoint segLocalCSC = hiti->localPosition();
                 //		cout<<"Layer Id (MuonHit) =  "<<cscDetIdHit<<" Muon Local Position (det frame) "<<positionLocalCSC <<endl;
                 //		cout<<"Layer Id  (CSCHit) =  "<<cscDetId<<"  Hit Local Position (det frame) "<<segLocalCSC <<endl;
                 if((fabs(positionLocalCSC.x()-segLocalCSC.x())<CSCXCut) && 
                    (fabs(positionLocalCSC.y()-segLocalCSC.y())<CSCYCut)) {
                   CSCcountAgreeingHits++;
                   //		  cout << "   Matched." << endl;
                 }  
               }//End 2D rechit iteration
             }//End muon hit iteration
             
             matchRatioCSC = countMuonCSCHits == 0 ? 0 : CSCcountAgreeingHits/countMuonCSCHits;
                         
             if ((matchRatioCSC>0.9) && ((countMuonCSCHits>1) || !cscTightMatch)) cscSegs.push_back(&(*segmentCSC));

           } //End CSC Segment Iteration 
           return cscSegs;
         }

         //Fwlite version of the code found in https://raw.githubusercontent.com/cms-sw/cmssw/CMSSW_7_6_X/RecoMuon/TrackingTools/src/MuonSegmentMatcher.cc
         vector<const DTRecSegment4D*>& matchDT(const reco::Track &muon, const DTRecSegment4DCollection& DTSegmentColl)
         {
           dtSegs.clear();

           std::vector<TrackingRecHit const *> dtHits;

           bool segments = false;

           // Loop and select DT recHits
           for(trackingRecHit_iterator hit = muon.recHitsBegin(); hit != muon.recHitsEnd(); ++hit) {
             if (!(*hit)->isValid()) continue; 
             if ( (*hit)->geographicalId().det() != DetId::Muon ) continue; 
             if ( (*hit)->geographicalId().subdetId() != MuonSubdetId::DT ) continue;
             if ((*hit)->recHits().size()) 
               if ((*(*hit)->recHits().begin())->recHits().size()>1) segments = true;
             dtHits.push_back(*hit);
           }
           
           //  cout << "Muon DT hits found: " << dtHits.size() << " segments " << segments << endl;
           
           double PhiCutParameter=DTradius;
           double ZCutParameter=DTradius;
           double matchRatioZ=0;
           double matchRatioPhi=0;

           for (DTRecSegment4DCollection::const_iterator rechit = DTSegmentColl.begin(); rechit!=DTSegmentColl.end();++rechit) {
           
             LocalPoint pointLocal = rechit->localPosition();

             if (segments) {
               // Loop over muon recHits
               for(auto hit = dtHits.begin(); hit != dtHits.end(); ++hit) {
                                                 
                 // Pick the one in the same DT Chamber as the muon
                 DetId idT = (*hit)->geographicalId();
                 if(!(rechit->geographicalId().rawId()==idT.rawId())) continue; 

                 // and compare the local positions
                 LocalPoint segLocal = (*hit)->localPosition();
                 if ((fabs(pointLocal.x()-segLocal.x())<ZCutParameter) && 
                     (fabs(pointLocal.y()-segLocal.y())<ZCutParameter)) 
                   dtSegs.push_back(&(*rechit));
               }
               continue;
             }

             double nhitsPhi = 0;
             double nhitsZ = 0;
                   
             if(rechit->hasZed()) {
               double countMuonDTHits = 0;
               double countAgreeingHits=0;

               const DTRecSegment2D* segmZ;
               segmZ = dynamic_cast<const DTRecSegment2D*>(rechit->zSegment());
               nhitsZ = segmZ->recHits().size(); 
                         
               const vector<DTRecHit1D> hits1d = segmZ->specificRecHits();
               DTChamberId chamberSegIdT((segmZ->geographicalId()).rawId());
                         
               // Loop over muon recHits
               for(auto hit = dtHits.begin(); hit != dtHits.end(); ++hit) {

                 if ( !(*hit)->isValid()) continue; 
                 
                 DetId idT = (*hit)->geographicalId();
                 DTChamberId dtDetIdHitT(idT.rawId());
                 DTSuperLayerId dtDetLayerIdHitT(idT.rawId());

                 LocalPoint  pointLocal = (*hit)->localPosition();
                                                 
                 if ((chamberSegIdT==dtDetIdHitT) && (dtDetLayerIdHitT.superlayer()==2)) countMuonDTHits++;

                 for (vector<DTRecHit1D>::const_iterator hiti=hits1d.begin(); hiti!=hits1d.end(); hiti++) {

                   if ( !hiti->isValid()) continue; 

                   // Pick the one in the same DT Layer as the 1D hit
                   if(!(hiti->geographicalId().rawId()==idT.rawId())) continue; 

                   // and compare the local positions
                   LocalPoint segLocal = hiti->localPosition();
         //        cout << "Zed Segment Point = "<<pointLocal<<"    Muon Point = "<<segLocal<<"  Dist:  "
         //             << (fabs(pointLocal.x()-segLocal.x()))+(fabs(pointLocal.y()-segLocal.y()))<< endl;
                   if ((fabs(pointLocal.x()-segLocal.x())<ZCutParameter) && 
                       (fabs(pointLocal.y()-segLocal.y())<ZCutParameter)) 
                     countAgreeingHits++;
                 } //End Segment Hit Iteration
               } //End Muon Hit Iteration
                         
               matchRatioZ = countMuonDTHits == 0 ? 0 : countAgreeingHits/countMuonDTHits;
               if (nhitsZ)
                 if (countAgreeingHits/nhitsZ>matchRatioZ) matchRatioZ=countAgreeingHits/nhitsZ;
             } //End HasZed Check
                                 
             if(rechit->hasPhi()) {
               double countMuonDTHits = 0;
               double countAgreeingHits=0;

               //PREPARE PARAMETERS FOR SEGMENT DETECTOR GEOMETRY
               const DTRecSegment2D* segmPhi;
               segmPhi = dynamic_cast<const DTRecSegment2D*>(rechit->phiSegment());
               nhitsPhi = segmPhi->recHits().size();
                         
               const vector<DTRecHit1D> hits1d = segmPhi->specificRecHits();
               DTChamberId chamberSegIdT((segmPhi->geographicalId()).rawId());
                         
               // Loop over muon recHits
               for(auto hit = dtHits.begin(); hit != dtHits.end(); ++hit) {

                 if ( !(*hit)->isValid()) continue; 

                 DetId idT = (*hit)->geographicalId();
                 DTChamberId dtDetIdHitT(idT.rawId());
                 DTSuperLayerId dtDetLayerIdHitT(idT.rawId());

                 LocalPoint pointLocal = (*hit)->localPosition(); //Localposition is in DTLayer http://cmslxr.fnal.gov/lxr/source/DataFormats/DTRecHit/interface/DTRecHit1D.h

                 if ((chamberSegIdT==dtDetIdHitT)&&((dtDetLayerIdHitT.superlayer()==1)||(dtDetLayerIdHitT.superlayer()==3))) 
                   countMuonDTHits++;

                 for (vector<DTRecHit1D>::const_iterator hiti=hits1d.begin(); hiti!=hits1d.end(); hiti++) {

                   if ( !hiti->isValid()) continue; 

                   // Pick the one in the same DT Layer as the 1D hit
                   if(!(hiti->geographicalId().rawId()==idT.rawId())) continue; 
                                                  
                   // and compare the local positions
                   LocalPoint segLocal = hiti->localPosition();
         //        cout << "     Phi Segment Point = "<<pointLocal<<"    Muon Point = "<<segLocal<<"  Dist:   " 
         //             << (fabs(pointLocal.x()-segLocal.x()))+(fabs(pointLocal.y()-segLocal.y()))<< endl;

                   if ((fabs(pointLocal.x()-segLocal.x())<PhiCutParameter) && 
                       (fabs(pointLocal.y()-segLocal.y())<PhiCutParameter))
                     countAgreeingHits++; 
                 } // End Segment Hit Iteration
               } // End Muon Hit Iteration

               matchRatioPhi = countMuonDTHits != 0 ? countAgreeingHits/countMuonDTHits : 0;
               if (nhitsPhi)
                 if (countAgreeingHits/nhitsPhi>matchRatioPhi) matchRatioPhi=countAgreeingHits/nhitsPhi;
             } // End HasPhi Check
         //    DTChamberId chamberSegId2((rechit->geographicalId()).rawId());
             if (TightMatchDT && nhitsPhi && nhitsZ) {
               if((matchRatioPhi>0.9)&&(matchRatioZ>0.9)) {
         //      cout<<"Making a tight match in Chamber "<<chamberSegId2<<endl;
                 dtSegs.push_back(&(*rechit));
               }
             } else {
               if((matchRatioPhi>0.9 && nhitsPhi)||(matchRatioZ>0.9 && nhitsZ)) {
         //      cout<<"Making a loose match in Chamber "<<chamberSegId2<<endl;
                 dtSegs.push_back(&(*rechit));
               }
             }
             
           } //End DT Segment Iteration
           return dtSegs;
         }


      ////////////////////////////////////////////
      //all code related to TOF computation
      ////////////////////////////////////////////


         double fitT0(double &a, double &b, const std::vector<double>& xl, const std::vector<double>& yl, const std::vector<double>& xr, const std::vector<double>& yr ) {

           double sx=0,sy=0,sxy=0,sxx=0,ssx=0,ssy=0,s=0,ss=0;

           for (unsigned int i=0; i<xl.size(); i++) {
             sx+=xl[i];
             sy+=yl[i];
             sxy+=xl[i]*yl[i];
             sxx+=xl[i]*xl[i];
             s++;
             ssx+=xl[i];
             ssy+=yl[i];
             ss++;
           } 

           for (unsigned int i=0; i<xr.size(); i++) {
             sx+=xr[i];
             sy+=yr[i];
             sxy+=xr[i]*yr[i];
             sxx+=xr[i]*xr[i];
             s++;
             ssx-=xr[i];
             ssy-=yr[i];
             ss--;
           } 

           double delta = ss*ss*sxx+s*sx*sx+s*ssx*ssx-s*s*sxx-2*ss*sx*ssx;
           
           double t0_corr=0.;

           if (delta) {
             a=(ssy*s*ssx+sxy*ss*ss+sy*sx*s-sy*ss*ssx-ssy*sx*ss-sxy*s*s)/delta;
             b=(ssx*sy*ssx+sxx*ssy*ss+sx*sxy*s-sxx*sy*s-ssx*sxy*ss-sx*ssy*ssx)/delta;
             t0_corr=(ssx*s*sxy+sxx*ss*sy+sx*sx*ssy-sxx*s*ssy-sx*ss*sxy-ssx*sx*sy)/delta;
           }

           // convert drift distance to time
           t0_corr/=-0.00543;

           return t0_corr;
         }



         double iBetaFromCSC(reco::MuonRef& muon, int CORRECTION_LEVEL){
                        //try to compute 1/beta for csc
                        std::vector<muonTimingCalculator::TimeMeasurementCSC> tms;
                        for(unsigned int ic=0;ic<cscSegs.size();ic++){ 
                           if(cscSegs[ic]->specificRecHits().size()<=0)continue;  

                           const std::vector<CSCRecHit2D> hits2d = cscSegs[ic]->specificRecHits();
                           for (std::vector<CSCRecHit2D>::const_iterator hiti=hits2d.begin(); hiti!=hits2d.end(); hiti++) {
                               muonTimingCalculator::TimeMeasurementCSC thisHit;

            //                  std::pair< TrajectoryStateOnSurface, double> tsos;
            //                  tsos=propag->propagateWithPath(muonFTS,cscDet->surface());

            //                  if (tsos.first.isValid()) thisHit.distIP = tsos.second+posp.mag(); 
            //                    else thisHit.distIP = cscDet->toGlobal(hiti->localPosition()).mag();
                               thisHit.distIP = moduleGeom::get(hiti->geographicalId())->pos.Mag();

                              if(UseStripTime_){
                                  thisHit.weightInvbeta = thisHit.distIP*thisHit.distIP/(theStripError_*theStripError_*30.*30.);
                                  thisHit.weightVertex = 1./(theStripError_*theStripError_);
                                  thisHit.timeCorr = hiti->tpeak();
                                  if(CORRECTION_LEVEL==1)thisHit.timeCorr -= t0Offset       (cscSegs[ic]->cscDetId().rawId());
                                  if(CORRECTION_LEVEL==2)thisHit.timeCorr -= t0OffsetChamber(cscSegs[ic]->cscDetId().rawId());

                                  tms.push_back(thisHit);
                              }

                              if(UseWireTime_){
                                  thisHit.weightInvbeta = thisHit.distIP*thisHit.distIP/(theWireError_*theWireError_*30.*30.);
                                  thisHit.weightVertex = 1./(theWireError_*theWireError_);
                                  thisHit.timeCorr = hiti->wireTime();
                                  if(CORRECTION_LEVEL==1)thisHit.timeCorr -= t0Offset       (cscSegs[ic]->cscDetId().rawId());
                                  if(CORRECTION_LEVEL==2)thisHit.timeCorr -= t0OffsetChamber(cscSegs[ic]->cscDetId().rawId());
                                  tms.push_back(thisHit);
                              }
                          } // rechit
                        }//csc segment



                          // Now loop over the measurements, calculate 1/beta and cut away outliers
                          double invbeta=0;
                          double invbetaerr=0;
                          double totalWeightInvbeta=0;
                          double totalWeightVertex=0;
                          bool modified = false;
                          std::vector <double> dstnc, dsegm, dtraj, hitWeightInvbeta, hitWeightVertex;
                          do {    
                            modified = false;
                            dstnc.clear();
                            dsegm.clear();
                            dtraj.clear();
                            hitWeightInvbeta.clear();
                            hitWeightVertex.clear();
                              
                            totalWeightInvbeta=0;
                            totalWeightVertex=0;
                              
                                for (std::vector<muonTimingCalculator::TimeMeasurementCSC>::iterator tm=tms.begin(); tm!=tms.end(); ++tm) {
                                  dstnc.push_back(tm->distIP);
                                  dsegm.push_back(tm->timeCorr);
                                  hitWeightInvbeta.push_back(tm->weightInvbeta);
                                  hitWeightVertex.push_back(tm->weightVertex);
                                  totalWeightInvbeta+=tm->weightInvbeta;
                                  totalWeightVertex+=tm->weightVertex;
                                }
                                  
                            if (totalWeightInvbeta==0) break;        

                            // calculate the value and error of 1/beta from the complete set of 1D hits
                            // inverse beta - weighted average of the contributions from individual hits
                            invbeta=0;
                            for (unsigned int i=0;i<dstnc.size();i++) 
                              invbeta+=(1.+dsegm.at(i)/dstnc.at(i)*30.)*hitWeightInvbeta.at(i)/totalWeightInvbeta;

                            double chimax=0.;
                            std::vector<muonTimingCalculator::TimeMeasurementCSC>::iterator tmmax;
                            
                            // the dispersion of inverse beta
                            double diff;
                            for (unsigned int i=0;i<dstnc.size();i++) {
                              diff=(1.+dsegm.at(i)/dstnc.at(i)*30.)-invbeta;
                              diff=diff*diff*hitWeightInvbeta.at(i);
                              invbetaerr+=diff;
                              if (diff>chimax) { 
                                tmmax=tms.begin()+i;
                                chimax=diff;
                              }
                            }
                            
                            invbetaerr=sqrt(invbetaerr/totalWeightInvbeta); 
                         
                            // cut away the outliers
                            if (chimax>thePruneCut_) {
                              tms.erase(tmmax);
                              modified=true;
                            }    

                          }while (modified);

                    

                          // std::cout << " *** FINAL Measured 1/beta: " << invbeta << " +/- " << invbetaerr << std::endl;

                          //save hit info to the keeping vector
                          for (unsigned int i=0;i<dstnc.size();i++) {
                             tmSeq.push_back(muonTimingCalculator::TimeMeasurement(dstnc[i], dsegm[i], hitWeightInvbeta[i], hitWeightVertex[i]));
                          }

                          return invbeta;
         }


         double iBetaFromDT(reco::MuonRef& muon, int CORRECTION_LEVEL){
         //                 printf("\n#Muon pT= %f Eta=%f Phi = %f\n", muon->pt(), muon->eta(), muon->phi());

                          std::vector<muonTimingCalculator::TimeMeasurementDT> tms;
                          // create a collection on TimeMeasurements for the track        
                          for (std::vector<const DTRecSegment4D*>::iterator rechit = dtSegs.begin(); rechit!=dtSegs.end();++rechit) {

                            // Create the ChamberId
                            DetId id = (*rechit)->geographicalId();
                            DTChamberId chamberId(id.rawId());
                            int station = chamberId.station();
         //                   if (debug) std::cout << "Matched DT segment in station " << station << std::endl;

                            // use only segments with both phi and theta projections present (optional)
                            bool bothProjections = ( ((*rechit)->hasPhi()) && ((*rechit)->hasZed()) );
                            
                            if (requireBothProjections_ && !bothProjections) continue;

                            // loop over (theta, phi) segments
                            for (int phi=0; phi<2; phi++) {

                              if (dropTheta_ && !phi) continue;

                              const DTRecSegment2D* segm;
                              if (phi) segm = dynamic_cast<const DTRecSegment2D*>((*rechit)->phiSegment()); 
                                  else segm = dynamic_cast<const DTRecSegment2D*>((*rechit)->zSegment());

                              if(segm == 0) continue;
                              if (!segm->specificRecHits().size()) continue;

                              moduleGeom* geomDet = moduleGeom::get(segm->geographicalId());
                              const std::vector<DTRecHit1D> hits1d = segm->specificRecHits();


                              // store all the hits from the segment
                              for (std::vector<DTRecHit1D>::const_iterator hiti=hits1d.begin(); hiti!=hits1d.end(); hiti++) {

                                moduleGeom* dtcell = moduleGeom::get(hiti->geographicalId());
                                muonTimingCalculator::TimeMeasurementDT thisHit;

         //                       std::pair< TrajectoryStateOnSurface, double> tsos;
         //                       tsos=propag->propagateWithPath(muonFTS,dtcell->surface());

                                double dist;            
                                double dist_straight = dtcell->toGlobal(TVector3(hiti->localPosition().x(), hiti->localPosition().y(), hiti->localPosition().z()) ).Mag(); 
         //                       if (tsos.first.isValid()) { 
         //                         dist = tsos.second+posp.mag(); 
         //               //        std::cout << "Propagate distance: " << dist << " ( innermost: " << posp.mag() << ")" << std::endl; 
         //                       } else { 
                                  dist = dist_straight;
         //               //        std::cout << "Geom distance: " << dist << std::endl; 
         //                       }

                                thisHit.driftCell = hiti->geographicalId();
                                if (hiti->lrSide()==DTEnums::Left) thisHit.isLeft=true; else thisHit.isLeft=false;
                                thisHit.isPhi = phi;
                                thisHit.posInLayer = geomDet->toLocal(dtcell->toGlobal(TVector3(hiti->localPosition().x(), hiti->localPosition().y(), hiti->localPosition().z()) )).x();

                                
                                thisHit.distIP = dist;
                                thisHit.station = station;
                                if (useSegmentT0_ && segm->ist0Valid()) thisHit.timeCorr=segm->t0();
                                else thisHit.timeCorr=0.;
         //                       thisHit.timeCorr += theTimeOffset_;
                                if(CORRECTION_LEVEL==1)thisHit.timeCorr -= t0Offset       (hiti->geographicalId());
                                if(CORRECTION_LEVEL==2)thisHit.timeCorr -= t0OffsetChamber(hiti->geographicalId());

                                  
                                // signal propagation along the wire correction for unmached theta or phi segment hits
         //                       if (doWireCorr_ && !bothProjections && tsos.first.isValid()) {
         //                         const DTLayer* layer = theDTGeom->layer(hiti->wireId());
         //                         float propgL = layer->toLocal( tsos.first.globalPosition() ).y();
         //                         float wirePropCorr = propgL/24.4*0.00543; // signal propagation speed along the wire
         //                         if (thisHit.isLeft) wirePropCorr=-wirePropCorr;
         //                         thisHit.posInLayer += wirePropCorr;
         //                         const DTSuperLayer *sl = layer->superLayer();
         //                         float tofCorr = sl->surface().position().mag()-tsos.first.globalPosition().mag();
         //                         tofCorr = (tofCorr/29.979)*0.00543;
         //                         if (thisHit.isLeft) tofCorr=-tofCorr;
         //                         thisHit.posInLayer += tofCorr;
         //                       } else {
                                  // non straight-line path correction
                                  float slCorr = (dist_straight-dist)/29.979*0.00543;
                                  if (thisHit.isLeft) slCorr=-slCorr;
                                  thisHit.posInLayer += slCorr;
         //                       }

         //                         printf("##Hit dist=%f   time=%f  posInLayer=%f\n", thisHit.distIP, thisHit.timeCorr, thisHit.posInLayer);


         //                       if (debug) std::cout << " dist: " << dist << "  t0: " << thisHit.posInLayer << std::endl;
                         
                                tms.push_back(thisHit);
                              }
                            } // phi = (0,1)            
                          } // rechit

                             
                          double invbeta=0;
                          double invbetaerr=0;
                          double totalWeightInvbeta=0;
                          double totalWeightVertex=0;
                          bool modified = false;
                          std::vector <double> dstnc, dsegm, dtraj, hitWeightVertex, hitWeightInvbeta, left;
                            
                          // Now loop over the measurements, calculate 1/beta and cut away outliers
                          do {    

                            modified = false;
                            dstnc.clear();
                            dsegm.clear();
                            dtraj.clear();
                            hitWeightVertex.clear();
                            hitWeightInvbeta.clear();
                            left.clear();
                              
                            std::vector <int> hit_idx;
                            totalWeightInvbeta=0;
                            totalWeightVertex=0;
                              
                            // Rebuild segments
                            for (int sta=1;sta<5;sta++)
                              for (int phi=0;phi<2;phi++) {
                                std::vector <muonTimingCalculator::TimeMeasurementDT> seg;
                                std::vector <int> seg_idx;
                                int tmpos=0;
                                for (std::vector<muonTimingCalculator::TimeMeasurementDT>::iterator tm=tms.begin(); tm!=tms.end(); ++tm) {
                                  if ((tm->station==sta) && (tm->isPhi==phi)) {
                                    seg.push_back(*tm);
                                    seg_idx.push_back(tmpos);
                                  }
                                  tmpos++;  
                                }

                                unsigned int segsize = seg.size();
                                if (segsize<theHitsMin_) continue;

                                double a=0, b=0;
                                std::vector <double> hitxl,hitxr,hityl,hityr;

                                for (std::vector<muonTimingCalculator::TimeMeasurementDT>::iterator tm=seg.begin(); tm!=seg.end(); ++tm) {

                                  DetId id = tm->driftCell;
                                  moduleGeom* dtcell = moduleGeom::get(id.rawId());
         //                         const GeomDet* dtcell = theTrackingGeometry->idToDet(id);
                                  DTChamberId chamberId(id.rawId());
         //                         const GeomDet* dtcham = theTrackingGeometry->idToDet(chamberId);
                                  moduleGeom* dtcham = moduleGeom::get(chamberId.rawId());

                                  double celly=dtcham->toLocal(dtcell->pos).z();
                                    
                                  if (tm->isLeft) {
                                    hitxl.push_back(celly);
                                    hityl.push_back(tm->posInLayer);
                                  } else {
                                    hitxr.push_back(celly);
                                    hityr.push_back(tm->posInLayer);
                                  }    
                                }

                                if (!fitT0(a,b,hitxl,hityl,hitxr,hityr)) {
                                  continue;
                                }

                                // a segment must have at least one left and one right hit
                                if ((!hitxl.size()) || (!hityl.size())) continue;

                                int segidx=0;
                                for (std::vector<muonTimingCalculator::TimeMeasurementDT>::const_iterator tm=seg.begin(); tm!=seg.end(); ++tm) {
                                  DetId id = tm->driftCell;
                                  moduleGeom* dtcell = moduleGeom::get(id.rawId());
         //                         const GeomDet* dtcell = theTrackingGeometry->idToDet(id);
                                  DTChamberId chamberId(id.rawId());
         //                         const GeomDet* dtcham = theTrackingGeometry->idToDet(chamberId);
                                  moduleGeom* dtcham = moduleGeom::get(chamberId.rawId());                       

                                  double layerZ  = dtcham->toLocal(dtcell->pos).z();
                                  double segmLocalPos = b+layerZ*a;
                                  double hitLocalPos = tm->posInLayer;
                                  int hitSide = -tm->isLeft*2+1;
                                  double t0_segm = (-(hitSide*segmLocalPos)+(hitSide*hitLocalPos))/0.00543+tm->timeCorr;
                                  
         //                         if (debug) std::cout << "   Segm hit.  dstnc: " << tm->distIP << "   t0: " << t0_segm << std::endl;
                                    
                                  dstnc.push_back(tm->distIP);
                                  dsegm.push_back(t0_segm);
                                  left.push_back(hitSide);
                                  hitWeightInvbeta.push_back(((double)seg.size()-2.)*tm->distIP*tm->distIP/((double)seg.size()*30.*30.*theError_*theError_));
                                  hitWeightVertex.push_back(((double)seg.size()-2.)/((double)seg.size()*theError_*theError_));
                                  hit_idx.push_back(seg_idx.at(segidx));
                                  segidx++;
                                  totalWeightInvbeta+=((double)seg.size()-2.)*tm->distIP*tm->distIP/((double)seg.size()*30.*30.*theError_*theError_);
                                  totalWeightVertex+=((double)seg.size()-2.)/((double)seg.size()*theError_*theError_);
                                }
                              }

                            if (totalWeightInvbeta==0) break;        

                            // calculate the value and error of 1/beta from the complete set of 1D hits
         //                   if (debug)
         //                     std::cout << " Points for global fit: " << dstnc.size() << std::endl;

                            // inverse beta - weighted average of the contributions from individual hits
                            invbeta=0;
                            for (unsigned int i=0;i<dstnc.size();i++){
                              invbeta+=(1.+dsegm.at(i)/dstnc.at(i)*30.)*hitWeightInvbeta.at(i)/totalWeightInvbeta;
                            }

                            double chimax=0.;
                            std::vector<muonTimingCalculator::TimeMeasurementDT>::iterator tmmax;
                            
                            // the dispersion of inverse beta
                            double diff;
                            for (unsigned int i=0;i<dstnc.size();i++) {
                              diff=(1.+dsegm.at(i)/dstnc.at(i)*30.)-invbeta;
                              diff=diff*diff*hitWeightInvbeta.at(i);
                              invbetaerr+=diff;
                              if (diff>chimax) { 
                                tmmax=tms.begin()+hit_idx.at(i);
                                chimax=diff;
                              }
                            }
                            
                            invbetaerr=sqrt(invbetaerr/totalWeightInvbeta); 
                         
                            // cut away the outliers
                            if (chimax>thePruneCutDT_) {
                              tms.erase(tmmax);
                              modified=true;
                            }    

         //                   if (debug)
         //                     std::cout << " Measured 1/beta: " << invbeta << " +/- " << invbetaerr << std::endl;

                          } while (modified);

         //                 printf("iBeta %f (fly) vs %f (aod)\n", invbeta, dttof->inverseBeta());


                          //save hit info to the keeping vector
                          for (unsigned int i=0;i<dstnc.size();i++) {
                             tmSeq.push_back(muonTimingCalculator::TimeMeasurement(dstnc[i], dsegm[i], hitWeightInvbeta[i], hitWeightVertex[i]));
                          }

                          return invbeta;
         }



         void rawFit(double &a, double &da, double &b, double &db, const std::vector<double>& hitsx, const std::vector<double>& hitsy) {
           double s=0,sx=0,sy=0,x,y;
           double sxx=0,sxy=0;

           a=b=0;
           if (hitsx.size()==0) return;
           if (hitsx.size()==1) {
             b=hitsy[0];
           } else {
             for (unsigned int i =0; i<hitsx.size(); i++) {
               x=hitsx[i];
               y=hitsy[i];
               sy += y;
               sxy+= x*y;
               s += 1.;
               sx += x;
               sxx += x*x;
             }

             double d = s*sxx - sx*sx;
             b = (sxx*sy- sx*sxy)/ d;
             a = (s*sxy - sx*sy) / d;
             da = sqrt(sxx/d);
             db = sqrt(s/d);
           }
         }



         reco::MuonTimeExtra fillTimeFromMeasurements(){
              reco::MuonTimeExtra muTime;

              std::vector <double> x,y;
              double invbeta=0, invbetaerr=0;
              double vertexTime=0, vertexTimeErr=0, vertexTimeR=0, vertexTimeRErr=0;
              double freeBeta, freeBetaErr, freeTime, freeTimeErr;

              if(tmSeq.size()>1){
                 float totalWeightInvBeta = 0;
                 float totalWeightVertex  = 0;
                 for (unsigned int i=0;i<tmSeq.size();i++){
                    totalWeightInvBeta += tmSeq[i].weightInvBeta;
                    totalWeightVertex += tmSeq[i].weightVertex;
                 }

                 for (unsigned int i=0;i<tmSeq.size();i++) {
                   invbeta+=(1.+tmSeq[i].localt0/tmSeq[i].distIP*30.)*tmSeq[i].weightInvBeta/totalWeightInvBeta;
                   x.push_back(tmSeq[i].distIP/30.);
                   y.push_back(tmSeq[i].localt0+tmSeq[i].distIP/30.);
                   vertexTime+=tmSeq[i].localt0*tmSeq[i].weightVertex/totalWeightVertex;
                   vertexTimeR+=(tmSeq[i].localt0+2*tmSeq[i].distIP/30.)*tmSeq[i].weightVertex/totalWeightVertex;
                 }

                 double diff;
                 for (unsigned int i=0;i<tmSeq.size();i++) {
                   diff=(1.+tmSeq[i].localt0/tmSeq[i].distIP*30.)-invbeta;
                   invbetaerr+=diff*diff*tmSeq[i].weightInvBeta;
                   diff=tmSeq[i].localt0-vertexTime;
                   vertexTimeErr+=diff*diff*tmSeq[i].weightVertex;
                   diff=tmSeq[i].localt0+2*tmSeq[i].distIP/30.-vertexTimeR;
                   vertexTimeRErr+=diff*diff*tmSeq[i].weightVertex;
                 }
                 
                 double cf = 1./(tmSeq.size()-1);
                 invbetaerr=sqrt(invbetaerr/totalWeightInvBeta*cf);
                 vertexTimeErr=sqrt(vertexTimeErr/totalWeightVertex*cf);
                 vertexTimeRErr=sqrt(vertexTimeRErr/totalWeightVertex*cf);
              }

              muTime.setInverseBeta(invbeta);
              muTime.setInverseBetaErr(invbetaerr);
              muTime.setTimeAtIpInOut(vertexTime);
              muTime.setTimeAtIpInOutErr(vertexTimeErr);
              muTime.setTimeAtIpOutIn(vertexTimeR);
              muTime.setTimeAtIpOutInErr(vertexTimeRErr);
                  
              rawFit(freeBeta, freeBetaErr, freeTime, freeTimeErr, x, y);

              muTime.setFreeInverseBeta(freeBeta);
              muTime.setFreeInverseBetaErr(freeBetaErr);
                
              muTime.setNDof(tmSeq.size());

              return muTime;
            }


};

