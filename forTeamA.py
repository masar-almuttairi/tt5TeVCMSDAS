from itertools import chain, product, repeat
from functools import partial 

def isGoodElectron(el, ptCut=10.):
    from bamboo import treefunctions as op
    return op.AND(
            el.pt > ptCut,
            op.abs(el.eta) < 2.5,
            ## TODO add the other cuts
            op.NOT(op.AND(el.jet.isValid, op.OR(el.jet.btagDeepB > .1522, el.jet.btagDeepB <= -999.)))
            ## ^^ this one is tricky, the rest are straightforward
            )
def isGoodMuon(mu, ptCut=10.):
    from bamboo import treefunctions as op
    return op.AND(
            mu.pt > ptCut,
            op.abs(mu.eta) < 2.4,
            ## TODO add the other cuts
            )

class DiLepton_A(Nano5TeVHistoModule):
    def definePlots(self, t, noSel, sample=None, sampleCfg=None):
        from bamboo.plots import Plot, CutFlowReport, SummedPlot
        from bamboo.plots import EquidistantBinning as EqB
        from bamboo import treefunctions as op

        isMC = self.isMC(sample)
        trigCut, trigWeight = None, None
        if isMC:
            noSel = noSel.refine("mcWeight", weight=[ t.genWeight, t.puWeight, t.PrefireWeight ])
            trigCut = op.OR(t.HLT.HIEle20_Ele12_CaloIdL_TrackIdL_IsoVL_DZ, t.HLT.HIL3DoubleMu0, t.HLT.HIL3Mu20, t.HLT.HIEle20_WPLoose_Gsf)
            ## TODO add a correction for prescaled triggers
        else:
            ## suggested trigger order: dielectron, dimuon or single muon, single electron (to minimise loss due to prescales). Electron triggered-events should be taken from the HighEGJet primary datasets, muon-triggered events from the SingleMuon primary datset
            pd = sample.split("_")[0]
            if pd == "SingleMuon":
                ## TODO fill trigger cut
            elif pd == "HighEGJet":
                ## TODO fill trigger cut
        noSel = noSel.refine("trig", cut=trigCut, weight=trigWeight)

        plots = []

        goodLeptons = {
            "el" : op.select(t.Electron, partial(isGoodElectron, ptCut=15.)),
            "mu" : op.select(t.Muon, partial(isGoodMuon, ptCut=15.))
            }
        plots += [
            Plot.make1D("trig_nLeptons15", op.rng_len(goodLeptons["el"])+op.rng_len(goodLeptons["mu"]), noSel, EqB(15, 0., 15.)),
            Plot.make1D("trig_nEl15", op.rng_len(goodLeptons["el"]), noSel, EqB(15, 0., 15.)),
            Plot.make1D("trig_nMu15", op.rng_len(goodLeptons["mu"]), noSel, EqB(15, 0., 15.)) 
            ]
        from bamboo.scalefactors import get_scalefactor
        sf_loose = {
            "mu": get_scalefactor("lepton", "Muon_RecoToLoose", sfLib=scalefactors_lepMVA, paramDefs=binningVariables_nano_noScaleSyst, systName="muLoose"),
            "el": get_scalefactor("lepton", "Electron_RecoToLoose", sfLib=scalefactors_lepMVA, paramDefs=binningVariables_nano_noScaleSyst, systName="elLoose")
            }
        sf_tight = {
            "mu": get_scalefactor("lepton", "Muon_LooseToTight", sfLib=scalefactors_lepMVA, paramDefs=binningVariables_nano_noScaleSyst, systName="muTight"),
            "el": get_scalefactor("lepton", "Electron_LooseToTight", sfLib=scalefactors_lepMVA, paramDefs=binningVariables_nano_noScaleSyst, systName="elTight")
            }

        nGoodLeptons = op.rng_len(goodLeptons["el"])+op.rng_len(goodLeptons["mu"])
        hasTwoGoodLeptons = noSel.refine("has2Lep", cut=(nGoodLeptons > 1)) # avoid overlap with 1l
        jets = op.sort(op.select(t.Jet, lambda j : op.AND(
            j.pt > 25., ## you decide...
            op.abs(j.eta) < 2.4,
            j.jetId & 0x2, ## tight JetID
            op.AND( ## lepton-jet cross-cleaning
                op.NOT(op.rng_any(goodLeptons["el"], lambda l : op.deltaR(l.p4, j.p4) < 0.4)),
                op.NOT(op.rng_any(goodLeptons["mu"], lambda l : op.deltaR(l.p4, j.p4) < 0.4)))
            )), lambda j : -j.pt)
        for fl1,fl2 in product(*repeat(goodLeptons.keys(), 2)):
            dilepSel = lambda l1,l2 : op.AND(
                    l1.charge != l2.charge,
                    (l1.p4+l2.p4).M() > 12.
                    )
            if fl1 == fl2:
                lGood = op.sort(goodLeptons[fl1], lambda l : -l.pt)
                dilep = op.combine(lGood, N=2, pred=dilepSel)
            else:
                l1Good = op.sort(goodLeptons[fl1], lambda l : -l.pt)
                l2Good = op.sort(goodLeptons[fl2], lambda l : -l.pt)
                dilep = op.combine((l1Good, l2Good), pred=dilepSel)
            ll = dilep[0]
            hasDilep = hasTwoGoodLeptons.refine(f"hasDilep{fl1}{fl2}", cut=(op.rng_len(dilep) > 0, ll[0].pt > 25.),
                    weight=([ sf_loose[fl1](ll[0]), sf_loose[fl2](ll[1]), sf_tight[fl1](ll[0]), sf_tight[fl2](ll[1]) ] if isMC else None))
            plots += [
                Plot.make1D(f"dilepton_{fl1}{fl2}_Mll", (ll[0].p4+ll[1].p4).M(), hasDilep, EqB(50, 70, 120.), title="Dilepton mass"),
                ]
            for il,ifl in enumerate((fl1, fl2)):
                plots += [
                    Plot.make1D(f"dilepton_{fl1}{fl2}_L{il:d}PT", ll[il].pt, hasDilep, EqB(50, 0., 100.), title=f"Lepton {il:d} PT"),
                    Plot.make1D(f"dilepton_{fl1}{fl2}_L{il:d}ETA", ll[il].eta, hasDilep, EqB(50, -2.5, 2.5), title=f"Lepton {il:d} ETA"),
                    ]
            plots += [
                Plot.make1D(f"dilepton_{fl1}{fl2}_nJets", op.rng_len(jets), hasDilep, EqB(15, 0, 15.), title="Jet multiplicity"),
                ]

        return plots
