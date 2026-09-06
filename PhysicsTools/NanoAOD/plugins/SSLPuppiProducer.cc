// Original Author: Yuji Li
// Created: Wed, 13 Mar 2024 09:19:58 GMT

#include <cmath>
#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "HeterogeneousCore/SonicTriton/interface/TritonEDProducer.h"

#include "TLorentzVector.h"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"

namespace {
  constexpr double kPi = 3.14159265358979323846;
  constexpr float kDeltaRCut = 0.4F;
  constexpr float kJetRadius = 0.4F;
  constexpr float kMaximumAbsEta = 2.5F;
  constexpr float kMinimumParticlePt = 0.5F;
  constexpr float kMinimumJetPt = 30.F;
  constexpr double kJetMatchingRadius = 0.1;
  constexpr float kOutputPaddingValue = -100.F;
  constexpr std::size_t kMinimumOutputSize = 40;
}  // namespace

class SSLPuppiProducer : public TritonEDProducer<> {
public:
  explicit SSLPuppiProducer(const edm::ParameterSet&);

  void acquire(edm::Event const& event, edm::EventSetup const&, Input& input) override;
  void produce(edm::Event& event, edm::EventSetup const&, Output const& output) override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  const edm::EDGetTokenT<std::vector<pat::PackedGenParticle>> genParticleToken_;
  const edm::EDGetTokenT<std::vector<pat::PackedCandidate>> pfToken_;
  const unsigned int maxNPF_;
};

SSLPuppiProducer::SSLPuppiProducer(const edm::ParameterSet& config)
    : TritonEDProducer<>(config),
      genParticleToken_(
          mayConsume<std::vector<pat::PackedGenParticle>>(config.getParameter<edm::InputTag>("genParticleSrc"))),
      pfToken_(consumes<std::vector<pat::PackedCandidate>>(config.getParameter<edm::InputTag>("pf_src"))),
      maxNPF_(config.getParameter<unsigned int>("max_n_pf")) {
  produces<std::vector<float>>("ptdiff");
  produces<std::vector<float>>("ptpfdiff");
  produces<std::vector<float>>("ptchsdiff");
  produces<std::vector<float>>("ptlvdiff");
  produces<std::vector<float>>("ptpuppidiff");
  produces<std::vector<float>>("massdiff");
  produces<std::vector<float>>("masspfdiff");
  produces<std::vector<float>>("masschsdiff");
  produces<std::vector<float>>("masslvdiff");
  produces<std::vector<float>>("masspuppidiff");
  produces<std::vector<int>>("massdiffgenidx");
  produces<std::vector<int>>("masspfdiffgenidx");
  produces<std::vector<int>>("masschsdiffgenidx");
  produces<std::vector<int>>("masslvdiffgenidx");
  produces<std::vector<int>>("masspuppidiffgenidx");
  produces<std::vector<float>>("GenJetMass");
  produces<std::vector<float>>("GenJetPt");
  produces<std::vector<float>>("PFJetMass");
  produces<std::vector<float>>("PFJetPt");
  produces<std::vector<float>>("CHSJetMass");
  produces<std::vector<float>>("CHSJetPt");
  produces<std::vector<float>>("LVJetMass");
  produces<std::vector<float>>("LVJetPt");
  produces<std::vector<float>>("PUPPIJetMass");
  produces<std::vector<float>>("PUPPIJetPt");
  produces<std::vector<float>>("SSLJetMass");
  produces<std::vector<float>>("SSLJetPt");
  produces<std::vector<float>>("GenJetEta");
  produces<std::vector<float>>("GenJetPhi");
  produces<std::vector<float>>("PFJetEta");
  produces<std::vector<float>>("PFJetPhi");
  produces<std::vector<float>>("CHSJetEta");
  produces<std::vector<float>>("CHSJetPhi");
  produces<std::vector<float>>("LVJetEta");
  produces<std::vector<float>>("LVJetPhi");
  produces<std::vector<float>>("PUPPIJetEta");
  produces<std::vector<float>>("PUPPIJetPhi");
  produces<std::vector<float>>("SSLJetEta");
  produces<std::vector<float>>("SSLJetPhi");
}

void SSLPuppiProducer::acquire(edm::Event const& event, edm::EventSetup const&, Input& input) {
  client_->setBatchSize(1);

  const auto& particleFlowCandidates = event.get(pfToken_);
  auto& nodeInput = input.at("INPUT0");
  auto& edgeInput = input.at("INPUT1");

  int numberOfNodes = 0;
  int numberOfEdges = 0;
  std::vector<float> particleEta;
  std::vector<float> particlePhi;
  std::size_t particleIndex = 0;

  for (const auto& candidate : particleFlowCandidates) {
    if (std::abs(candidate.eta()) > kMaximumAbsEta) {
      continue;
    }
    ++numberOfNodes;
  }

  nodeInput.setShape(0, numberOfNodes);
  auto nodeData = nodeInput.allocate<float>();
  auto& nodeValues = (*nodeData)[0];
  nodeValues.clear();

  for (const auto& candidate : particleFlowCandidates) {
    if (std::abs(candidate.eta()) > kMaximumAbsEta) {
      continue;
    }

    nodeValues.push_back(candidate.eta());
    nodeValues.push_back(candidate.phi());
    nodeValues.push_back(candidate.pt());
    particleEta.push_back(candidate.eta());
    particlePhi.push_back(candidate.phi());

    if (candidate.charge() != 0) {
      nodeValues.push_back(1);
      nodeValues.push_back(0);
      nodeValues.push_back(0);
    }
    if (candidate.charge() == 0) {
      if (candidate.pdgId() == 22) {
        nodeValues.push_back(0);
        nodeValues.push_back(1);
        nodeValues.push_back(0);
      } else {
        nodeValues.push_back(0);
        nodeValues.push_back(0);
        nodeValues.push_back(1);
      }
    }

    if (candidate.charge() != 0) {
      if (candidate.puppiWeight() == 1) {
        nodeValues.push_back(0);
        nodeValues.push_back(1);
        nodeValues.push_back(0);
      } else {
        nodeValues.push_back(1);
        nodeValues.push_back(0);
        nodeValues.push_back(0);
      }
    }
    if (candidate.charge() == 0) {
      nodeValues.push_back(0);
      nodeValues.push_back(0);
      nodeValues.push_back(1);
    }

    nodeValues.push_back(0);
    ++particleIndex;
    if (particleIndex == maxNPF_) {
      break;
    }
  }

  if (particleEta.empty()) {
    return;
  }

  for (std::size_t first = 0; first < particleEta.size(); ++first) {
    for (std::size_t second = 0; second < particleEta.size(); ++second) {
      float deltaPhi = std::fabs(particlePhi[first] - particlePhi[second]);
      const float deltaEta = std::fabs(particleEta[first] - particleEta[second]);
      if (deltaPhi > kPi) {
        const float offset = std::ceil((deltaPhi - kPi) / (2 * kPi)) * (2 * kPi);
        deltaPhi = deltaPhi - offset;
      }
      const float deltaR = std::sqrt(deltaPhi * deltaPhi + deltaEta * deltaEta);
      if (first == second) {
        continue;
      }
      if (deltaR < kDeltaRCut) {
        ++numberOfEdges;
      }
    }
  }

  edgeInput.setShape(1, numberOfEdges);
  auto edgeData = edgeInput.allocate<long>();
  auto& edgeValues = (*edgeData)[0];
  edgeValues.clear();
  std::vector<long> edgeSources;
  std::vector<long> edgeTargets;

  for (std::size_t source = 0; source < particleEta.size(); ++source) {
    for (std::size_t target = 0; target < particleEta.size(); ++target) {
      float deltaPhi = std::fabs(particlePhi[source] - particlePhi[target]);
      const float deltaEta = std::fabs(particleEta[source] - particleEta[target]);
      if (deltaPhi > kPi) {
        const float offset = std::ceil((deltaPhi - kPi) / (2 * kPi)) * (2 * kPi);
        deltaPhi = deltaPhi - offset;
      }
      const float deltaR = std::sqrt(deltaPhi * deltaPhi + deltaEta * deltaEta);

      if (deltaR < kDeltaRCut) {
        if (source == target) {
          continue;
        }
        edgeSources.push_back(static_cast<long>(source));
        edgeTargets.push_back(static_cast<long>(target));
      }
    }
  }

  for (const long source : edgeSources) {
    edgeValues.push_back(source);
  }
  for (const long target : edgeTargets) {
    edgeValues.push_back(target);
  }

  nodeValues.resize(8 * maxNPF_);
  nodeInput.toServer(nodeData);
  edgeInput.toServer(edgeData);
}

void SSLPuppiProducer::produce(edm::Event& event, edm::EventSetup const&, Output const& output) {
  const auto& outputTensor = output.at("OUTPUT0");
  const auto& inferenceOutput = outputTensor.fromServer<float>();
  const auto& sslScores = inferenceOutput[0];

  const auto& particleFlowCandidates = event.get(pfToken_);
  edm::Handle<std::vector<pat::PackedGenParticle>> genParticles;
  event.getByToken(genParticleToken_, genParticles);

  std::vector<fastjet::PseudoJet> pfJetInputs;
  std::vector<fastjet::PseudoJet> puppiJetInputs;
  std::vector<fastjet::PseudoJet> sslJetInputs;
  std::vector<fastjet::PseudoJet> chsJetInputs;
  std::vector<fastjet::PseudoJet> lvJetInputs;

  int numberOfLVParticles = 0;
  int numberOfPUParticles = 0;
  std::size_t scoreIndex = 0;

  for (const auto& candidate : particleFlowCandidates) {
    if (std::abs(candidate.eta()) > kMaximumAbsEta) {
      continue;
    }

    const float sslScore = sslScores[scoreIndex];
    if (candidate.charge() != 0 && candidate.pt() > kMinimumParticlePt && candidate.fromPV() > 2 &&
        candidate.puppiWeight() > 0.99) {
      ++numberOfLVParticles;
    }
    if (candidate.charge() != 0 && candidate.pt() > kMinimumParticlePt && candidate.fromPV() < 1 &&
        candidate.puppiWeight() < 0.01) {
      ++numberOfPUParticles;
    }

    if (candidate.pt() > kMinimumParticlePt) {
      TLorentzVector pfP4;
      TLorentzVector puppiP4;
      TLorentzVector sslP4;
      TLorentzVector chsP4;
      TLorentzVector lvP4;

      pfP4.SetPtEtaPhiM(candidate.pt(), candidate.eta(), candidate.phi(), 0);
      if (candidate.charge() == 0) {
        sslP4.SetPtEtaPhiM(candidate.pt() * sslScore, candidate.eta(), candidate.phi(), 0);
      } else {
        sslP4.SetPtEtaPhiM(candidate.pt() * candidate.puppiWeight(), candidate.eta(), candidate.phi(), 0);
      }
      if (candidate.charge() == 0) {
        chsP4.SetPtEtaPhiM(candidate.pt(), candidate.eta(), candidate.phi(), 0);
      } else {
        chsP4.SetPtEtaPhiM(candidate.pt() * candidate.puppiWeight(), candidate.eta(), candidate.phi(), 0);
      }
      puppiP4.SetPtEtaPhiM(candidate.pt() * candidate.puppiWeight(), candidate.eta(), candidate.phi(), 0);
      if (candidate.puppiWeight() > 0.99) {
        lvP4.SetPtEtaPhiM(candidate.pt(), candidate.eta(), candidate.phi(), 0);
      }

      pfJetInputs.emplace_back(pfP4.Px(), pfP4.Py(), pfP4.Pz(), pfP4.E());
      sslJetInputs.emplace_back(sslP4.Px(), sslP4.Py(), sslP4.Pz(), sslP4.E());
      puppiJetInputs.emplace_back(puppiP4.Px(), puppiP4.Py(), puppiP4.Pz(), puppiP4.E());
      chsJetInputs.emplace_back(chsP4.Px(), chsP4.Py(), chsP4.Pz(), chsP4.E());
      lvJetInputs.emplace_back(lvP4.Px(), lvP4.Py(), lvP4.Pz(), lvP4.E());
    }
    ++scoreIndex;
  }

  std::vector<fastjet::PseudoJet> genJetInputs;
  for (const auto& particle : *genParticles) {
    if (particle.status() != 1) {
      continue;
    }
    if (std::abs(particle.pdgId()) == 12) {
      continue;
    }
    if (std::abs(particle.pdgId()) == 14) {
      continue;
    }
    if (std::abs(particle.pdgId()) == 16) {
      continue;
    }
    if (std::abs(particle.eta()) > kMaximumAbsEta) {
      continue;
    }

    if (particle.pt() > kMinimumParticlePt) {
      TLorentzVector genParticleP4;
      genParticleP4.SetPtEtaPhiM(particle.pt(), particle.eta(), particle.phi(), 0);
      genJetInputs.emplace_back(genParticleP4.Px(), genParticleP4.Py(), genParticleP4.Pz(), genParticleP4.E());
    }
  }

  const fastjet::JetDefinition jetDefinition(fastjet::antikt_algorithm, kJetRadius);
  const fastjet::ClusterSequence pfClusterSequence(pfJetInputs, jetDefinition);
  const fastjet::ClusterSequence sslClusterSequence(sslJetInputs, jetDefinition);
  const fastjet::ClusterSequence puppiClusterSequence(puppiJetInputs, jetDefinition);
  const fastjet::ClusterSequence chsClusterSequence(chsJetInputs, jetDefinition);
  const fastjet::ClusterSequence lvClusterSequence(lvJetInputs, jetDefinition);
  const fastjet::ClusterSequence genClusterSequence(genJetInputs, jetDefinition);

  const std::vector<fastjet::PseudoJet> pfJets = fastjet::sorted_by_pt(pfClusterSequence.inclusive_jets(kMinimumJetPt));
  const std::vector<fastjet::PseudoJet> puppiJets =
      fastjet::sorted_by_pt(puppiClusterSequence.inclusive_jets(kMinimumJetPt));
  const std::vector<fastjet::PseudoJet> sslJets =
      fastjet::sorted_by_pt(sslClusterSequence.inclusive_jets(kMinimumJetPt));
  const std::vector<fastjet::PseudoJet> genJets =
      fastjet::sorted_by_pt(genClusterSequence.inclusive_jets(kMinimumJetPt));
  const std::vector<fastjet::PseudoJet> chsJets =
      fastjet::sorted_by_pt(chsClusterSequence.inclusive_jets(kMinimumJetPt));
  const std::vector<fastjet::PseudoJet> lvJets = fastjet::sorted_by_pt(lvClusterSequence.inclusive_jets(kMinimumJetPt));

  std::vector<float> ptDiff;
  std::vector<float> ptPFDiff;
  std::vector<float> ptCHSDiff;
  std::vector<float> ptLVDiff;
  std::vector<float> ptPUPPIDiff;
  std::vector<float> massDiff;
  std::vector<float> massPFDiff;
  std::vector<float> massCHSDiff;
  std::vector<float> massLVDiff;
  std::vector<float> massPUPPIDiff;
  std::vector<int> massDiffGenIndex;
  std::vector<int> massPFDiffGenIndex;
  std::vector<int> massCHSDiffGenIndex;
  std::vector<int> massLVDiffGenIndex;
  std::vector<int> massPUPPIDiffGenIndex;
  std::vector<float> genJetMass;
  std::vector<float> genJetPt;
  std::vector<float> genJetEta;
  std::vector<float> genJetPhi;
  std::vector<float> sslJetMass;
  std::vector<float> sslJetPt;
  std::vector<float> sslJetEta;
  std::vector<float> sslJetPhi;
  std::vector<float> chsJetMass;
  std::vector<float> chsJetPt;
  std::vector<float> chsJetEta;
  std::vector<float> chsJetPhi;
  std::vector<float> lvJetMass;
  std::vector<float> lvJetPt;
  std::vector<float> lvJetEta;
  std::vector<float> lvJetPhi;
  std::vector<float> pfJetMass;
  std::vector<float> pfJetPt;
  std::vector<float> pfJetEta;
  std::vector<float> pfJetPhi;
  std::vector<float> puppiJetMass;
  std::vector<float> puppiJetPt;
  std::vector<float> puppiJetEta;
  std::vector<float> puppiJetPhi;

  int genIndex = 0;
  for (const auto& genJet : genJets) {
    if (numberOfLVParticles < 5 || numberOfPUParticles < 50) {
      continue;
    }

    TLorentzVector genP4;
    genP4.SetPxPyPzE(genJet.px(), genJet.py(), genJet.pz(), genJet.e());
    genJetEta.push_back(genP4.Eta());
    genJetPhi.push_back(genP4.Phi());
    genJetPt.push_back(genP4.Pt());
    genJetMass.push_back(genP4.M());

    for (const auto& pfJet : pfJets) {
      TLorentzVector pfP4;
      pfP4.SetPxPyPzE(pfJet.px(), pfJet.py(), pfJet.pz(), pfJet.e());
      if (pfP4.DeltaR(genP4) < kJetMatchingRadius) {
        massPFDiff.push_back((pfP4.M() - genP4.M()) / genP4.M());
        ptPFDiff.push_back((pfP4.Pt() - genP4.Pt()) / genP4.Pt());
        massPFDiffGenIndex.push_back(genIndex);
        pfJetMass.push_back(pfP4.M());
        pfJetPt.push_back(pfP4.Pt());
        pfJetEta.push_back(pfP4.Eta());
        pfJetPhi.push_back(pfP4.Phi());
      }
    }

    for (const auto& puppiJet : puppiJets) {
      TLorentzVector puppiP4;
      puppiP4.SetPxPyPzE(puppiJet.px(), puppiJet.py(), puppiJet.pz(), puppiJet.e());
      if (puppiP4.DeltaR(genP4) < kJetMatchingRadius) {
        massPUPPIDiff.push_back((puppiP4.M() - genP4.M()) / genP4.M());
        ptPUPPIDiff.push_back((puppiP4.Pt() - genP4.Pt()) / genP4.Pt());
        massPUPPIDiffGenIndex.push_back(genIndex);
        puppiJetMass.push_back(puppiP4.M());
        puppiJetPt.push_back(puppiP4.Pt());
        puppiJetEta.push_back(puppiP4.Eta());
        puppiJetPhi.push_back(puppiP4.Phi());
      }
    }

    for (const auto& chsJet : chsJets) {
      TLorentzVector chsP4;
      chsP4.SetPxPyPzE(chsJet.px(), chsJet.py(), chsJet.pz(), chsJet.e());
      if (chsP4.DeltaR(genP4) < kJetMatchingRadius) {
        massCHSDiff.push_back((chsP4.M() - genP4.M()) / genP4.M());
        ptCHSDiff.push_back((chsP4.Pt() - genP4.Pt()) / genP4.Pt());
        massCHSDiffGenIndex.push_back(genIndex);
        chsJetMass.push_back(chsP4.M());
        chsJetPt.push_back(chsP4.Pt());
        chsJetEta.push_back(chsP4.Eta());
        chsJetPhi.push_back(chsP4.Phi());
      }
    }

    for (const auto& lvJet : lvJets) {
      TLorentzVector lvP4;
      lvP4.SetPxPyPzE(lvJet.px(), lvJet.py(), lvJet.pz(), lvJet.e());
      if (lvP4.DeltaR(genP4) < kJetMatchingRadius) {
        massLVDiff.push_back((lvP4.M() - genP4.M()) / genP4.M());
        ptLVDiff.push_back((lvP4.Pt() - genP4.Pt()) / genP4.Pt());
        massLVDiffGenIndex.push_back(genIndex);
        lvJetMass.push_back(lvP4.M());
        lvJetPt.push_back(lvP4.Pt());
        lvJetEta.push_back(lvP4.Eta());
        lvJetPhi.push_back(lvP4.Phi());
      }
    }

    for (const auto& sslJet : sslJets) {
      TLorentzVector sslP4;
      sslP4.SetPxPyPzE(sslJet.px(), sslJet.py(), sslJet.pz(), sslJet.e());
      if (sslP4.DeltaR(genP4) < kJetMatchingRadius) {
        massDiff.push_back((sslP4.M() - genP4.M()) / genP4.M());
        ptDiff.push_back((sslP4.Pt() - genP4.Pt()) / genP4.Pt());
        massDiffGenIndex.push_back(genIndex);
        sslJetMass.push_back(sslP4.M());
        sslJetPt.push_back(sslP4.Pt());
        sslJetEta.push_back(sslP4.Eta());
        sslJetPhi.push_back(sslP4.Phi());
      }
    }

    ++genIndex;
  }

  for (std::size_t index = massDiff.size(); index < kMinimumOutputSize; ++index) {
    massDiff.push_back(kOutputPaddingValue);
    ptDiff.push_back(kOutputPaddingValue);
    massDiffGenIndex.push_back(-100);
    sslJetMass.push_back(kOutputPaddingValue);
    sslJetPt.push_back(kOutputPaddingValue);
    sslJetEta.push_back(kOutputPaddingValue);
    sslJetPhi.push_back(kOutputPaddingValue);
  }
  for (std::size_t index = massPFDiff.size(); index < kMinimumOutputSize; ++index) {
    massPFDiff.push_back(kOutputPaddingValue);
    ptPFDiff.push_back(kOutputPaddingValue);
    massPFDiffGenIndex.push_back(-100);
    pfJetMass.push_back(kOutputPaddingValue);
    pfJetPt.push_back(kOutputPaddingValue);
    pfJetEta.push_back(kOutputPaddingValue);
    pfJetPhi.push_back(kOutputPaddingValue);
  }
  for (std::size_t index = massPUPPIDiff.size(); index < kMinimumOutputSize; ++index) {
    massPUPPIDiff.push_back(kOutputPaddingValue);
    ptPUPPIDiff.push_back(kOutputPaddingValue);
    massPUPPIDiffGenIndex.push_back(-100);
    puppiJetMass.push_back(kOutputPaddingValue);
    puppiJetPt.push_back(kOutputPaddingValue);
    puppiJetEta.push_back(kOutputPaddingValue);
    puppiJetPhi.push_back(kOutputPaddingValue);
  }
  for (std::size_t index = massCHSDiff.size(); index < kMinimumOutputSize; ++index) {
    massCHSDiff.push_back(kOutputPaddingValue);
    ptCHSDiff.push_back(kOutputPaddingValue);
    massCHSDiffGenIndex.push_back(-100);
    chsJetMass.push_back(kOutputPaddingValue);
    chsJetPt.push_back(kOutputPaddingValue);
    chsJetEta.push_back(kOutputPaddingValue);
    chsJetPhi.push_back(kOutputPaddingValue);
  }
  for (std::size_t index = massLVDiff.size(); index < kMinimumOutputSize; ++index) {
    massLVDiff.push_back(kOutputPaddingValue);
    ptLVDiff.push_back(kOutputPaddingValue);
    massLVDiffGenIndex.push_back(-100);
    lvJetMass.push_back(kOutputPaddingValue);
    lvJetPt.push_back(kOutputPaddingValue);
    lvJetEta.push_back(kOutputPaddingValue);
    lvJetPhi.push_back(kOutputPaddingValue);
  }
  for (std::size_t index = genJetPt.size(); index < kMinimumOutputSize; ++index) {
    genJetPt.push_back(kOutputPaddingValue);
    genJetPhi.push_back(kOutputPaddingValue);
    genJetEta.push_back(kOutputPaddingValue);
    genJetMass.push_back(kOutputPaddingValue);
  }

  event.put(std::make_unique<std::vector<float>>(std::move(ptDiff)), "ptdiff");
  event.put(std::make_unique<std::vector<float>>(std::move(ptPFDiff)), "ptpfdiff");
  event.put(std::make_unique<std::vector<float>>(std::move(ptCHSDiff)), "ptchsdiff");
  event.put(std::make_unique<std::vector<float>>(std::move(ptLVDiff)), "ptlvdiff");
  event.put(std::make_unique<std::vector<float>>(std::move(ptPUPPIDiff)), "ptpuppidiff");
  event.put(std::make_unique<std::vector<float>>(std::move(massDiff)), "massdiff");
  event.put(std::make_unique<std::vector<float>>(std::move(massPFDiff)), "masspfdiff");
  event.put(std::make_unique<std::vector<float>>(std::move(massCHSDiff)), "masschsdiff");
  event.put(std::make_unique<std::vector<float>>(std::move(massLVDiff)), "masslvdiff");
  event.put(std::make_unique<std::vector<float>>(std::move(massPUPPIDiff)), "masspuppidiff");
  event.put(std::make_unique<std::vector<int>>(std::move(massDiffGenIndex)), "massdiffgenidx");
  event.put(std::make_unique<std::vector<int>>(std::move(massPFDiffGenIndex)), "masspfdiffgenidx");
  event.put(std::make_unique<std::vector<int>>(std::move(massCHSDiffGenIndex)), "masschsdiffgenidx");
  event.put(std::make_unique<std::vector<int>>(std::move(massLVDiffGenIndex)), "masslvdiffgenidx");
  event.put(std::make_unique<std::vector<int>>(std::move(massPUPPIDiffGenIndex)), "masspuppidiffgenidx");
  event.put(std::make_unique<std::vector<float>>(std::move(genJetEta)), "GenJetEta");
  event.put(std::make_unique<std::vector<float>>(std::move(genJetPt)), "GenJetPt");
  event.put(std::make_unique<std::vector<float>>(std::move(genJetPhi)), "GenJetPhi");
  event.put(std::make_unique<std::vector<float>>(std::move(genJetMass)), "GenJetMass");
  event.put(std::make_unique<std::vector<float>>(std::move(puppiJetEta)), "PUPPIJetEta");
  event.put(std::make_unique<std::vector<float>>(std::move(puppiJetPhi)), "PUPPIJetPhi");
  event.put(std::make_unique<std::vector<float>>(std::move(puppiJetPt)), "PUPPIJetPt");
  event.put(std::make_unique<std::vector<float>>(std::move(puppiJetMass)), "PUPPIJetMass");
  event.put(std::make_unique<std::vector<float>>(std::move(sslJetEta)), "SSLJetEta");
  event.put(std::make_unique<std::vector<float>>(std::move(sslJetPhi)), "SSLJetPhi");
  event.put(std::make_unique<std::vector<float>>(std::move(sslJetPt)), "SSLJetPt");
  event.put(std::make_unique<std::vector<float>>(std::move(sslJetMass)), "SSLJetMass");
  event.put(std::make_unique<std::vector<float>>(std::move(chsJetEta)), "CHSJetEta");
  event.put(std::make_unique<std::vector<float>>(std::move(chsJetPhi)), "CHSJetPhi");
  event.put(std::make_unique<std::vector<float>>(std::move(chsJetPt)), "CHSJetPt");
  event.put(std::make_unique<std::vector<float>>(std::move(chsJetMass)), "CHSJetMass");
  event.put(std::make_unique<std::vector<float>>(std::move(lvJetEta)), "LVJetEta");
  event.put(std::make_unique<std::vector<float>>(std::move(lvJetPhi)), "LVJetPhi");
  event.put(std::make_unique<std::vector<float>>(std::move(lvJetPt)), "LVJetPt");
  event.put(std::make_unique<std::vector<float>>(std::move(lvJetMass)), "LVJetMass");
  event.put(std::make_unique<std::vector<float>>(std::move(pfJetEta)), "PFJetEta");
  event.put(std::make_unique<std::vector<float>>(std::move(pfJetPhi)), "PFJetPhi");
  event.put(std::make_unique<std::vector<float>>(std::move(pfJetPt)), "PFJetPt");
  event.put(std::make_unique<std::vector<float>>(std::move(pfJetMass)), "PFJetMass");
}

void SSLPuppiProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription description;
  TritonClient::fillPSetDescription(description);
  description.add<edm::InputTag>("genParticleSrc", edm::InputTag("packedGenParticles"));
  description.add<edm::InputTag>("pf_src", edm::InputTag("packedPFCandidates"));

  // Keep these legacy parameters so configurations made before histogram removal remain valid.
  description.add<edm::InputTag>("JetsAK8", edm::InputTag("slimmedJetsAK8"));
  description.add<edm::InputTag>("GENJetsAK8", edm::InputTag("slimmedGenJetsAK8"));

  description.add<unsigned int>("max_n_pf", 9000);
  descriptions.add("SSLPuppiProducer", description);
}

DEFINE_FWK_MODULE(SSLPuppiProducer);
