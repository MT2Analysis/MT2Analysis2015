# Script to write the cutflow for several regions of the MT2 analysis
# specify the sample name of the babytree, the year of data-taking, and the region where you want the cutflow

# setup:
# a recent root version and a python version >=2.7

# For the moment implemented only for SR and CR1lep
# TODO: QCD region (and Photon CR)
# TODO: missing id cut on electrons for zll CR
# Note: the so-called 'cleanings' are commented out because take forever to run (why?); anyhow the impact is very small


from ROOT import TChain, TH1F

cuts={}
#cut[<year>][<cutName>] = <cutExpr>

cuts['2016'] = {}
# cuts of the signal region and single-lepton CR
cuts['2016']['noCut'] =             '(1)'
cuts['2016']['isGolden'] =          '(isGolden)'
cuts['2016']['goodVertex'] =        '(nVert > 0)'
cuts['2016']['filters'] =           '(Flag_goodVertices>0 && Flag_HBHENoiseFilter>0 && Flag_HBHENoiseIsoFilter>0 && Flag_globalTightHalo2016Filter>0 && Flag_EcalDeadCellTriggerPrimitiveFilter>0 && Flag_eeBadScFilter>0 && Flag_badMuonFilterV2>0 && Flag_badChargedHadronFilterV2>0)'
cuts['2016']['cleanings'] =         '(nJet200MuFrac50DphiMet == 0 && met_miniaodPt/met_caloPt<5.0 && jet_pt[0] < 13000)'
cuts['2016']['triggers'] =          '(HLT_PFMET120_PFMHT120 || HLT_PFHT900 || HLT_PFHT300_PFMET110 || HLT_PFJet450 || HLT_PFMETNoMu120_PFMHTNoMu120)'
cuts['2016']['HT_Etmiss'] =         '(( nJet30>1 && ht<1000. && met_pt>250.) || ( nJet30>1 && ht>=1000. && met_pt>30.) || (nJet30==1 && met_pt>250.) )'
cuts['2016']['nJets>0'] =           '(nJet30>=1 && nJet30FailId == 0)'
cuts['2016']['deltaPhi>0.3'] =      '(deltaPhiMin > 0.3)'
cuts['2016']['HTmiss-Etmiss'] =     '(diffMetMht < 0.5*met_pt)'
cuts['2016']['MT2'] =               '((((mt2>200 && ht<1500) || (mt2>400&&ht>=1500)) &&(nJet30>1) )  ||  ( ht>200 && nJet30==1 )  )'
cuts['2016']['leptonVeto'] =        '(nMuons10==0 && nElectrons10==0)'
cuts['2016']['isoTrackVeto'] =      '(nPFLep5LowMT==0 && nPFHad10LowMT==0)'
cuts['2016']['1lepton'] =           '(nLepLowMT==1)' # should be equivalent to '(nMuons10==1 || nElectrons10==1 || nPFLep5LowMT==1 || nPFHad10LowMT==1 )'
cuts['2016']['incl1leptonCR'] =     '(ht>250. && met_pt>250  && nJet30>1 && mt2>200.)'
# cuts that are dilepton CR specific
cuts['2016']['zll_SF_triggers'] =   '( lep_pdgId[0] == -lep_pdgId[1] && (HLT_DoubleMu || HLT_DoubleMu_NonIso || HLT_SingleMu_NonIso || HLT_DoubleEl || HLT_DoubleEl33 || HLT_Photon165_HE10))' 
cuts['2016']['zll_OF_triggers'] =   '( lep_pdgId[0] != -lep_pdgId[1] && ( HLT_MuX_Ele12 || HLT_Mu8_EleX || HLT_Mu33_Ele33_NonIso || HLT_Mu30_Ele30_NonIso || HLT_Photon165_HE10 || HLT_SingleMu_NonIso ))'
cuts['2016']['zll_HT_Etmiss'] =     '((nJet30>1 && zll_ht<1000. && zll_met_pt>250.) || ( nJet30>1 && zll_ht>=1000. && zll_met_pt>30.) || (nJet30==1 && zll_met_pt>250.) )'
cuts['2016']['zll_deltaPhi>0.3'] =  '(zll_deltaPhiMin > 0.3)'
cuts['2016']['zll_HTmiss-Etmiss'] = '(zll_diffMetMht < 0.5*zll_met_pt)'
cuts['2016']['zll_MT2'] =           '((((zll_mt2>200 && zll_ht<1500) || (zll_mt2>400 && zll_ht>=1500)) &&(nJet30>1) )  ||  ( zll_ht>200 && nJet30==1 )  )'
cuts['2016']['ZmassPt'] =           '(zll_pt > 200. && fabs(zll_mass-91.19)<20)'
cuts['2016']['2leptons'] =          '(nLep==2 && lep_pdgId[0]*lep_pdgId[1])>0 && lep_pt[0]>100 && lep_pt[1]>30 &&) '
# abs(myTree.lep_pdgId[0])==11 && myTree.lep_tightId[0]< 0.5
# abs(myTree.lep_pdgId[1])==11 && myTree.lep_tightId[1]< 0.5 

cuts['2017'] = {}
cuts['2017']['noCut'] =             cuts['2016']['noCut']
cuts['2017']['isGolden'] =          cuts['2016']['isGolden']
cuts['2017']['goodVertex'] =        cuts['2016']['goodVertex']
cuts['2017']['filters'] =           '(Flag_goodVertices>0 && Flag_globalTightHalo2016Filter>0 && Flag_HBHENoiseFilter>0 && Flag_HBHENoiseIsoFilter>0 && Flag_EcalDeadCellTriggerPrimitiveFilter>0 && Flag_BadPFMuonFilter>0 && Flag_eeBadScFilter>0 && Flag_ecalBadCalibFilter>0)'
cuts['2017']['cleanings'] =         cuts['2016']['cleanings']
cuts['2017']['triggers'] =          '(HLT_PFMET120_PFMHT120 || HLT_PFHT1050 || HLT_PFHT500_PFMET100_PFMHT100 || HLT_PFJet500 || HLT_PFMETNoMu120_PFMHTNoMu120 || HLT_PFMETNoMu120_PFMHTNoMu120_PFHT60 )'
cuts['2017']['HT_Etmiss'] =         '(( nJet30>1 && ht<1200. && met_pt>250.) || ( nJet30>1 && ht>=1200. && met_pt>30.) || (nJet30==1 && met_pt>250.) )'
cuts['2017']['nJets>0'] =           cuts['2016']['nJets>0']
cuts['2017']['deltaPhi>0.3'] =      cuts['2016']['deltaPhi>0.3']
cuts['2017']['HTmiss-Etmiss'] =     cuts['2016']['HTmiss-Etmiss']
cuts['2017']['leptonVeto'] =        cuts['2016']['leptonVeto']
cuts['2017']['isoTrackVeto'] =      cuts['2016']['isoTrackVeto']
cuts['2017']['MT2'] =               cuts['2016']['MT2']
cuts['2017']['1lepton'] =           cuts['2016']['1lepton']
cuts['2017']['incl1leptonCR'] =     cuts['2016']['incl1leptonCR']
# cuts that are dilepton CR specific
cuts['2017']['zll_SF_triggers'] =   cuts['2016']['zll_SF_triggers']
cuts['2017']['zll_OF_triggers'] =   cuts['2016']['zll_OF_triggers']
cuts['2017']['zll_HT_Etmiss'] =     '((nJet30>1 && zll_ht<1000. && zll_met_pt>250.) || ( nJet30>1 && zll_ht>=1000. && zll_met_pt>30.) || (nJet30==1 && zll_met_pt>250.) )'
cuts['2017']['zll_deltaPhi>0.3'] =  cuts['2016']['zll_deltaPhi>0.3']
cuts['2017']['zll_HTmiss-Etmiss'] = cuts['2016']['zll_HTmiss-Etmiss']
cuts['2017']['zll_MT2'] =           cuts['2016']['zll_MT2']
cuts['2017']['ZmassPt'] =           cuts['2016']['ZmassPt']
cuts['2017']['2leptons'] =          cuts['2016']['2leptons']

ordered_cutNames = {}
ordered_cutNames['SR'] =     ['noCut' ,'isGolden', 'goodVertex', 'filters', 'triggers', 'HT_Etmiss', 'nJets>0', 'deltaPhi>0.3', 'HTmiss-Etmiss', 'MT2', 'leptonVeto','isoTrackVeto', ]
ordered_cutNames['CR1lep'] = ['noCut' ,'isGolden', 'goodVertex', 'filters', 'triggers', 'HT_Etmiss', 'nJets>0', 'deltaPhi>0.3', 'HTmiss-Etmiss', 'MT2', '1lepton', 'incl1leptonCR']
ordered_cutNames['CR2lep'] = ['noCut' ,'isGolden', 'goodVertex', 'filters', 'zll_SF_triggers', 'zll_HT_Etmiss', 'nJets>0', 'zll_deltaPhi>0.3', 'zll_HTmiss-Etmiss', 'zll_MT2', 'ZmassPt', '2leptons']
#ordered_cutNames['CR1ph'] =
#ordered_cutNames['SR'] = ['noCut', 'cleanings']

class SampleCutflow(object):
  def __init__(self, sampleName, fileNames,treeName, year, region):
    self.sampleName=sampleName
    self.fileNames=fileNames
    self.treeName=treeName
    self.year=year
    self.region=region
    cutflow={}
    for cutName in ordered_cutNames:
      cutflow[cutName] = -9999
    self.cutflow=cutflow
    efficiency={}
    efficiency['rel']={}
    efficiency['abs']={}
    for cutName in ordered_cutNames[self.region]:
      efficiency['abs'][cutName]= -9999
      efficiency['rel'][cutName]= -9999
    self.efficiency=efficiency

  def fill(self):
    chain = TChain(self.treeName)
    for fileName in self.fileNames:
      chain.Add(fileName)
    selection = '(1)'
    for cutName in ordered_cutNames[self.region]:
      selection = selection + '&&' + cuts[year][cutName]
      histo = TH1F('histo', 'histo', 10, 0., 14000)
      print 'going to use this selection', selection
      chain.Project(histo.GetName(), 'met_pt', selection)
      self.cutflow[cutName] = float(histo.Integral())
      del histo

  def fillEfficiency(self):
    for i,cutName in enumerate(ordered_cutNames[self.region]):
      self.efficiency['abs'][cutName]=self.cutflow[cutName]/self.cutflow[ordered_cutNames[self.region][0]]*100
      if i==0:
        self.efficiency['rel'][cutName]=self.efficiency['abs'][cutName]
      else:
        self.efficiency['rel'][cutName]=self.cutflow[cutName]/self.cutflow[ordered_cutNames[self.region][i-1]]*100

  def stamp(self):
    print '\n\n\n'
    print 'Cutflow for sample={s}, year={y}, region={r}'.format(s=self.sampleName, y=self.year, r=self.region)
    print 'File names {n}'.format(n='\t\n'.join(self.fileNames))
    print ''

    print '{:>15} {:>15} {:>15} {:>15}'.format('Cut', 'Cutflow', 'Rel eff (%)', 'Abs eff (%)')
    for cutName in ordered_cutNames[self.region]:
      print '{:>15} {:>15} {:>15} {:>15}'.format(cutName, self.cutflow[cutName], '{:.2f}'.format(self.efficiency['rel'][cutName]), '{:.2f}'.format(self.efficiency['abs'][cutName]))
    print ''

if __name__ == "__main__":


  fileNames = ['root://t3dcachedb.psi.ch//pnfs/psi.ch/cms/trivcat/store/user/mratti/MT2production/80X/PostProcessed/TEST0_postProc_data2017_TEST0_v3/merged_MET_HTMHT_JetHT_SingleElectron_DoubleEG_DoubleMuon_MuonEG.root']
  sampleName = 'Data 2017, periods=BCDEF, Streams=MET, JetHT, HTMHT, SingleElectron, DoubleMuon, DoubleEG, MuonEG'
  year = '2017'
  treeName = 'mt2'
  regions = ['CR1lep'] #'SR'


  #fileNames = ['root://t3dcachedb.psi.ch//pnfs/psi.ch/cms/trivcat/store/user/mschoene/MT2production/80X/PostProcessed/reMini2017_runB_Feb28_postProc_Mar04/skimAndPrune/merged_DoubleEG_DoubleMuon_HTMHT_JetHT_MET_MuonEG_SingleElectron_SingleMuon_SinglePhoton.root',
#'root://t3dcachedb.psi.ch//pnfs/psi.ch/cms/trivcat/store/user/mschoene/MT2production/80X/PostProcessed/reMini2017_runC_Feb28_postProc_Mar04/skimAndPrune/merged_DoubleEG_DoubleMuon_HTMHT_JetHT_MET_MuonEG_SingleElectron_SingleMuon_SinglePhoton.root',
#'root://t3dcachedb.psi.ch//pnfs/psi.ch/cms/trivcat/store/user/mschoene/MT2production/80X/PostProcessed/reMini2017_runD_Feb28_postProc_Mar04/skimAndPrune/merged_DoubleEG_DoubleMuon_HTMHT_JetHT_MET_MuonEG_SingleElectron_SingleMuon_SinglePhoton.root',
#'root://t3dcachedb.psi.ch//pnfs/psi.ch/cms/trivcat/store/user/mschoene/MT2production/80X/PostProcessed/reMini2017_runE_Feb28_postProc_Mar04/skimAndPrune/merged_DoubleEG_DoubleMuon_HTMHT_JetHT_MET_MuonEG_SingleElectron_SingleMuon_SinglePhoton.root',
#'root://t3dcachedb.psi.ch//pnfs/psi.ch/cms/trivcat/store/user/mschoene/MT2production/80X/PostProcessed/reMini2017_runF_Feb28_postProc_Mar06/skimAndPrune/merged_DoubleEG_DoubleMuon_HTMHT_JetHT_MET_MuonEG_SingleElectron_SingleMuon_SinglePhoton.root',
#'root://t3dcachedb.psi.ch//pnfs/psi.ch/cms/trivcat/store/user/mschoene/MT2production/80X/PostProcessed/reMini2017_runG_Feb28_postProc_Mar06/skimAndPrune/merged_DoubleEG_DoubleMuon_HTMHT_JetHT_MET_MuonEG_SingleElectron_SingleMuon_SinglePhoton.root',
#'root://t3dcachedb.psi.ch//pnfs/psi.ch/cms/trivcat/store/user/mschoene/MT2production/80X/PostProcessed/reMini2017_runH_Feb28_postProc_Mar04/skimAndPrune/merged_DoubleEG_DoubleMuon_HTMHT_JetHT_MET_MuonEG_SingleElectron_SingleMuon_SinglePhoton.root'
#              ]
#  sampleName = 'Data 2016, periods=CDEFGH, Streams=MET, JetHT, HTMHT, SingleElectron, DoubleMuon, DoubleEG, MuonEG, SinglePhoton, SingleMuon'
#  year = '2016'
  treeName = 'mt2'

  for region in regions:
    myCutflow = SampleCutflow(sampleName, fileNames, treeName, year, region)
    myCutflow.fill()
    myCutflow.fillEfficiency()
    myCutflow.stamp()
