import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.nano_eras_cff import *
from PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.simpleCandidateFlatTableProducer_cfi import simpleCandidateFlatTableProducer
from PhysicsTools.NanoAOD.simpleSingletonCandidateFlatTableProducer_cfi import simpleSingletonCandidateFlatTableProducer
from PhysicsTools.NanoAOD.globalVariablesTableProducer_cfi import globalVariablesTableProducer
from PhysicsTools.NanoAOD.SSLPuppiProducer_cfi import SSLPuppiProducer as _SSLPuppiProducer

SSLPuppiProducer = _SSLPuppiProducer.clone(
    Client = dict(
        timeout = 300,
        mode = "Async",
        modelName = "PuppiGNN",
        modelConfigPath = "/afs/cern.ch/user/y/yujil/Test/CMSSW_13_3_1/src/TritonDemo/models/PuppiGNN/config.pbtxt",
        # version "1" is the resolutionTune
        # version "2" is the responeTune
        modelVersion = "1",
    ),
    pf_src = "packedPFCandidates",
    
)

SSLTable = globalVariablesTableProducer.clone(
    #src = cms.InputTag("SSLPuppiProducer","SSL"),
    name = cms.string("massdiff"),
    extension = cms.bool(True),
    variables = cms.PSet(
       massdiff = ExtVar(cms.InputTag("SSLPuppiProducer:massdiff"), "float", doc="massdiff from SSL",precision=10),
    ),
)
SSLTask = cms.Task(SSLPuppiProducer)
SSLTablesTask = cms.Task(SSLTable)