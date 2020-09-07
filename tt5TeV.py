from bamboo.analysismodules import NanoAODModule, NanoAODHistoModule, NanoAODSkimmerModule
import os.path

from bamboo.treedecorators import NanoAODDescription, nanoPUWeightVar, ReadJetMETVar
description_nanov4_5TeV_data = NanoAODDescription.get("v5", year="2018", isMC=False,
        systVariations=[ ReadJetMETVar("Jet", "MET", jetsNomName="nom", metNomName="nom") ],
        removeGroups=["ChsMET_", "L1_", "METFixEE2017_"],
        removeCollections=["nCorrT1METJet"]
        )
description_nanov4_5TeV_MC = NanoAODDescription.get("v5", year="2018", isMC=False,
        systVariations=[ nanoPUWeightVar, ReadJetMETVar("Jet", "MET", jetsExclVars=["raw", "jesHF", "jesHF_2018", "jesHF_2017", "jesHF_2016"], metNomName="nom", metExclVars=["raw", "nom"], bTaggers=["csvv2", "deepcsv", "deepjet", "cmva"], bTagWPs=["L", "M", "T", "shape"]) ],
        removeGroups=["ChsMET_", "L1_", "METFixEE2017_", "CaloMET_", "PuppiMET_", "TkMET_", "HTXS_"],
        removeCollections=["nCorrT1METJet", "nFatJet", "nOtherPV", "nPhoton", "nSV", "nSoftActivityJet", "nSubJet", "nTau", "nTrigObj"]
        )

class Nano5TeVBase(NanoAODModule):
    """ Base module for custom 5TeV NanoAODv4 samples """
    def mergeCounters(self, outF, infileNames, sample=None):
        ## merge all histograms in the input files
        from bamboo.root import gbl
        aFile = gbl.TFile.Open(infileNames[0])
        outF.cd()
        histos = { ky.GetName(): aFile.Get(ky.GetName()).Clone() for ky in aFile.GetListOfKeys() if ky.GetClassName() == "TH1F" }
        for iFN in infileNames[1:]:
            for i in range(5):
                oFile = gbl.TFile.Open(iFN)
                if oFile:
                    break
                else:
                    print(f"Could not open {iFN}, sleeping a bit and then retrying...")
                    time.sleep(10)
            if not oFile:
                raise RuntimeError(f"Could not open {iFN} after five attempts")
            for hName, h in histos.items():
                hO = oFile.Get(hName)
                if hO:
                    h.Add(oFile.Get(hName))
        outF.cd()
        for h in histos.values():
            h.Write()
    def readCounters(self, resultsFile):
        dCount = {}
        for ky in resultsFile.GetListOfKeys():
            if ky.GetClassName() == "TH1F":
                obj = resultsFile.Get(ky.GetName())
                if obj.GetNbinsX() == 1:
                    dCount[ky.GetName()] = obj.GetBinContent(1)
                else:
                    for i in range(obj.GetNbinsX()):
                        dCount[f"{ky.GetName()}{i:d}"] = obj.GetBinContent(i+1)
        return dCount
    def prepareTree(self, tree, sample=None, sampleCfg=None):
        ## Decorate the tree
        from bamboo.treedecorators import NanoAODDescription, nanoPUWeightVar, ReadJetMETVar
        tree,noSel,be,lumiArgs = super(Nano5TeVBase, self).prepareTree(tree, sample=sample, sampleCfg=sampleCfg,
                description=(description_nanov4_5TeV_MC if self.isMC(sample) else description_nanov4_5TeV_data))
        return tree,noSel,be,lumiArgs

class Nano5TeVHistoModule(Nano5TeVBase, NanoAODHistoModule):
    pass

class HelloWorld(Nano5TeVHistoModule):
    def definePlots(self, t, noSel, sample=None, sampleCfg=None):
        from bamboo.plots import Plot, CutFlowReport, SummedPlot
        from bamboo.plots import EquidistantBinning as EqB
        from bamboo import treefunctions as op

        isMC = self.isMC(sample)
        if isMC:
            noSel = noSel.refine("mcWeight", weight=t.genWeight)
        noSel = noSel.refine("trig", cut=op.OR(t.HLT.HIL3DoubleMu0, t.HLT.HIEle20_Ele12_CaloIdL_TrackIdL_IsoVL_DZ))

        plots = []

        muons = op.select(t.Muon, lambda mu : mu.pt > 20.)
        twoMuSel = noSel.refine("twoMuons", cut=[ op.rng_len(muons) > 1 ])
        plots.append(Plot.make1D("dimu_M",
            op.invariant_mass(muons[0].p4, muons[1].p4), twoMuSel, EqB(100, 20., 120.),
            title="Dimuon invariant mass", plotopts={"show-overflow":False,
            "legend-position": [0.2, 0.6, 0.5, 0.9]}))

        return plots
