#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include <iostream>

class FloatVectorTableProducer : public edm::stream::EDProducer<> {
public:
  explicit FloatVectorTableProducer(const edm::ParameterSet& cfg)
  : src_(consumes<std::vector<float>>(cfg.getParameter<edm::InputTag>("src"))),
    name_(cfg.getParameter<std::string>("name")),
    doc_(cfg.getParameter<std::string>("doc")),
    extension_(cfg.getParameter<bool>("extension")) {
    produces<nanoaod::FlatTable>();
    std::cout<<"Producing:"<<std::endl;
  }
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src", edm::InputTag("SSLPuppiProducer:massdiff"));
    desc.add<std::string>("name", "");
    desc.add<std::string>("doc", "");
    desc.add<bool>("extension", true);
    descriptions.add("FloatVectorTableProducer",desc);
  }
  void produce(edm::Event& iEvent, const edm::EventSetup&) override {
    edm::Handle<std::vector<float>> vecHandle;
    iEvent.getByToken(src_, vecHandle);
    const std::vector<float>& vec = *vecHandle;
    std::cout<<"Producing:"<<vec[0]<<std::endl;
    nanoaod::FlatTable tab(1, name_, /*singleton=*/false, extension_);
    tab.addColumnValue<float>("values", vec[0], doc_, 1);
    std::unique_ptr<nanoaod::FlatTable> out(new nanoaod::FlatTable(tab));
    iEvent.put(std::move(out));
  }

private:
  edm::EDGetTokenT<std::vector<float>> src_;
  std::string name_, doc_;
  bool extension_;
};


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(FloatVectorTableProducer);
