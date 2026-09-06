import FWCore.ParameterSet.Config as cms

from PhysicsTools.NanoAOD.SSLPuppiProducer_cfi import SSLPuppiProducer as _sslPuppiProducer
from PhysicsTools.NanoAOD.common_cff import ExtVar as _ExtVar
from PhysicsTools.NanoAOD.fastVariablesTableProducer_cfi import (
    fastVariablesTableProducer as _fastVariablesTableProducer,
)
from PhysicsTools.NanoAOD.fastintTableProducer_cfi import fastintTableProducer as _fastIntTableProducer
from PhysicsTools.NanoAOD.globalVariablesTableProducer_cfi import (
    globalVariablesTableProducer as _globalVariablesTableProducer,
)


SSLPuppiProducer = _sslPuppiProducer.clone(
    Client=dict(
        timeout=300,
        mode="Async",
        modelName="PuppiGNN",
        modelConfigPath=cms.FileInPath("TritonDemo/models/PuppiGNN/config.pbtxt"),
        # Version 1 is resolution-tuned; version 2 is response-tuned.
        modelVersion="1",
    ),
    pf_src=cms.InputTag("packedPFCandidates"),
)


def _make_table(producer, source_product, value_type, *, output_name=None):
    """Build a one-column table with one canonical name for all output metadata."""
    output_name = source_product if output_name is None else output_name
    return producer.clone(
        name=cms.string(output_name),
        extension=cms.bool(False),
        variables=cms.PSet(
            **{
                output_name: _ExtVar(
                    cms.InputTag("SSLPuppiProducer", source_product),
                    value_type,
                    doc=f"{output_name} from SSL",
                    precision=10,
                )
            }
        ),
    )


def _make_float_table(source_product, *, output_name=None):
    return _make_table(
        _fastVariablesTableProducer,
        source_product,
        "float",
        output_name=output_name,
    )


def _make_int_table(source_product, *, output_name=None):
    return _make_table(
        _fastIntTableProducer,
        source_product,
        "int",
        output_name=output_name,
    )


# This legacy global table is intentionally not scheduled by the standard SSL
# table task.
sslLegacyJetRelativeMassResidualGlobalTable = _make_table(
    _globalVariablesTableProducer,
    "massdiff",
    "float",
)

# Mass residuals.
sslJetRelativeMassResidualTable = _make_float_table("massdiff")
sslPFJetRelativeMassResidualTable = _make_float_table("masspfdiff")
sslCHSJetRelativeMassResidualTable = _make_float_table("masschsdiff")
sslPuppiJetRelativeMassResidualTable = _make_float_table("masspuppidiff")
sslLeadingVertexJetRelativeMassResidualTable = _make_float_table("masslvdiff")

# Generator-level jet kinematics.
sslGenJetEtaTable = _make_float_table("GenJetEta")
sslGenJetPtTable = _make_float_table("GenJetPt")
sslGenJetPhiTable = _make_float_table("GenJetPhi")
sslGenJetMassTable = _make_float_table("GenJetMass")

# Generator-match indices.
sslJetGenMatchIndexTable = _make_int_table("massdiffgenidx")
sslPFJetGenMatchIndexTable = _make_int_table("masspfdiffgenidx")
sslCHSJetGenMatchIndexTable = _make_int_table("masschsdiffgenidx")
sslPuppiJetGenMatchIndexTable = _make_int_table("masspuppidiffgenidx")
sslLeadingVertexJetGenMatchIndexTable = _make_int_table("masslvdiffgenidx")

# Reconstructed jet transverse momenta and masses.
sslPFJetPtTable = _make_float_table("PFJetPt")
sslPFJetMassTable = _make_float_table("PFJetMass")
sslCHSJetPtTable = _make_float_table("CHSJetPt")
sslCHSJetMassTable = _make_float_table("CHSJetMass")
# PUPPIJet* are legacy EDM product labels; use Puppi consistently in NanoAOD.
sslPuppiJetPtTable = _make_float_table(
    source_product="PUPPIJetPt",
    output_name="PuppiJetPt",
)
sslPuppiJetMassTable = _make_float_table(
    source_product="PUPPIJetMass",
    output_name="PuppiJetMass",
)
sslJetPtTable = _make_float_table("SSLJetPt")
sslJetMassTable = _make_float_table("SSLJetMass")
sslLeadingVertexJetPtTable = _make_float_table("LVJetPt")
sslLeadingVertexJetMassTable = _make_float_table("LVJetMass")

# Transverse-momentum residuals.
sslJetRelativePtResidualTable = _make_float_table("ptdiff")
sslPFJetRelativePtResidualTable = _make_float_table("ptpfdiff")
sslCHSJetRelativePtResidualTable = _make_float_table("ptchsdiff")
sslPuppiJetRelativePtResidualTable = _make_float_table("ptpuppidiff")
sslLeadingVertexJetRelativePtResidualTable = _make_float_table("ptlvdiff")

# Reconstructed jet directions.
sslPFJetEtaTable = _make_float_table("PFJetEta")
sslPFJetPhiTable = _make_float_table("PFJetPhi")
sslLeadingVertexJetEtaTable = _make_float_table("LVJetEta")
sslLeadingVertexJetPhiTable = _make_float_table("LVJetPhi")
sslPuppiJetEtaTable = _make_float_table(
    source_product="PUPPIJetEta",
    output_name="PuppiJetEta",
)
sslPuppiJetPhiTable = _make_float_table(
    source_product="PUPPIJetPhi",
    output_name="PuppiJetPhi",
)
sslCHSJetEtaTable = _make_float_table("CHSJetEta")
sslCHSJetPhiTable = _make_float_table("CHSJetPhi")
sslJetEtaTable = _make_float_table("SSLJetEta")
sslJetPhiTable = _make_float_table("SSLJetPhi")


sslPuppiProducerTask = cms.Task(SSLPuppiProducer)

sslPuppiTablesTask = cms.Task(
    sslJetRelativeMassResidualTable,
    sslPFJetRelativeMassResidualTable,
    sslCHSJetRelativeMassResidualTable,
    sslPuppiJetRelativeMassResidualTable,
    sslLeadingVertexJetRelativeMassResidualTable,
    sslGenJetEtaTable,
    sslGenJetPtTable,
    sslGenJetPhiTable,
    sslGenJetMassTable,
    sslJetGenMatchIndexTable,
    sslPFJetGenMatchIndexTable,
    sslCHSJetGenMatchIndexTable,
    sslPuppiJetGenMatchIndexTable,
    sslLeadingVertexJetGenMatchIndexTable,
    sslPFJetPtTable,
    sslPFJetMassTable,
    sslCHSJetPtTable,
    sslCHSJetMassTable,
    sslPuppiJetPtTable,
    sslPuppiJetMassTable,
    sslJetPtTable,
    sslJetMassTable,
    sslLeadingVertexJetPtTable,
    sslLeadingVertexJetMassTable,
    sslJetRelativePtResidualTable,
    sslPFJetRelativePtResidualTable,
    sslCHSJetRelativePtResidualTable,
    sslPuppiJetRelativePtResidualTable,
    sslLeadingVertexJetRelativePtResidualTable,
    sslPFJetEtaTable,
    sslPFJetPhiTable,
    sslLeadingVertexJetEtaTable,
    sslLeadingVertexJetPhiTable,
    sslPuppiJetEtaTable,
    sslPuppiJetPhiTable,
    sslCHSJetEtaTable,
    sslCHSJetPhiTable,
    sslJetEtaTable,
    sslJetPhiTable,
)
