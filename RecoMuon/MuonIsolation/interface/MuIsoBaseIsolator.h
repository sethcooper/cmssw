#ifndef MuonIsolation_MuIsoBaseIsolator_H
#define MuonIsolation_MuIsoBaseIsolator_H

#include <vector>
#include "DataFormats/MuonReco/interface/MuIsoDeposit.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Candidate/interface/Candidate.h"

namespace muonisolation {
  class MuIsoBaseIsolator {

  public:
    typedef reco::MuIsoDeposit::Veto Veto;
    typedef reco::MuIsoDeposit::Vetos Vetos;

    struct DepositAndVetos {
      DepositAndVetos(): dep(0), vetos(0) {}
      DepositAndVetos(const reco::MuIsoDeposit* depA, const Vetos* vetosA = 0):
	dep(depA), vetos(vetosA) {}
      const reco::MuIsoDeposit* dep;
      const Vetos* vetos;
    };
    typedef std::vector<DepositAndVetos> DepositContainer;
  
    enum ResultType {
      ISOL_INT_TYPE = 0,
      ISOL_FLOAT_TYPE,
      ISOL_BOOL_TYPE,
      ISOL_INVALID_TYPE
    };

    class Result {
    public:
      Result() : valInt(-999), valFloat(-999), valBool(false), typeF_(ISOL_INVALID_TYPE) {}
	Result(ResultType typ) : valInt(-999), valFloat(-999.), valBool(false), typeF_(typ) {}
	  
	  template <typename T> T val() const;
	  
	  int valInt;
	  float valFloat;
	  bool valBool;
	  ResultType typeF() const {return typeF_;}

    protected:
	  ResultType typeF_;
    };
    

    virtual ~MuIsoBaseIsolator(){}

    //! Compute and return the isolation variable
    virtual Result result(DepositContainer deposits) const = 0;
    //! Compute and return the isolation variable, with vetoes and the muon
    virtual Result result(DepositContainer deposits, const reco::Candidate& muon) const {
      return result(deposits);
    }
    //! Compute and return the isolation variable, with vetoes and the muon
    virtual Result result(DepositContainer deposits, const reco::Track& muon) const {
      return result(deposits);
    }

    virtual ResultType resultType() const = 0;

  };

  template<> inline int MuIsoBaseIsolator::Result::val<int>() const { return valInt;}
  template<> inline float MuIsoBaseIsolator::Result::val<float>() const { return valFloat;}
  template<> inline bool MuIsoBaseIsolator::Result::val<bool>() const { return valBool;}
  
}
#endif

