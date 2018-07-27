# Script to write the cutflow for several regions of the MT2 analysis
# specify the sample name of the babytree, the year of data-taking, and the region where you want the cutflow

# setup:
# a recent root version and a python version >=2.7

# example:
# python cutflow.py -f /shome/mratti/heppy_workarea/CMSSW_9_4_1/src/CMGTools/TTHAnalysis/cfg/testmc_94X_std_5K/WJetsToLNu_HT-600To800/mt2.root -s Wlv -t mt2 -y 2017 -g SR_kinematic -d mini

# For the moment implemented only for SR and CR1lep
# TODO: QCD region (and Photon CR?)
# TODO: missing id cut on electrons for zll CR, add nlep==2 requirement, same signa requirement yielding 0 evts ? lep_pdgId[0]*lep_pdgId[1]>0 && 
# TODO: add 2016 for nano
# TODO: eventually converge on nano and remove mini
# NOTE: the so-called 'cleanings' are commented out because take forever to run (why?); anyhow the impact is very small


from ROOT import TChain, TH1F
import os.path

cuts={}
#cut[<year>][<cutName>] = <cutExpr>
cuts['mini']={}
cuts['mini']['2016'] = {}
# cuts of the signal region and single-lepton CR
cuts['mini']['2016']['noCut'] =             '(1)'
cuts['mini']['2016']['isGolden'] =          '(isGolden)'
cuts['mini']['2016']['goodVertex'] =        '(nVert > 0)'
cuts['mini']['2016']['filters'] =           '(Flag_goodVertices>0 && Flag_HBHENoiseFilter>0 && Flag_HBHENoiseIsoFilter>0 && Flag_globalTightHalo2016Filter>0 && Flag_EcalDeadCellTriggerPrimitiveFilter>0 && Flag_eeBadScFilter>0 && Flag_badMuonFilterV2>0 && Flag_badChargedHadronFilterV2>0)'
cuts['mini']['2016']['cleanings'] =         '(nJet200MuFrac50DphiMet == 0 && met_miniaodPt/met_caloPt<5.0 && jet_pt[0] < 13000)'
cuts['mini']['2016']['triggers'] =          '(HLT_PFMET120_PFMHT120 || HLT_PFHT900 || HLT_PFHT300_PFMET110 || HLT_PFJet450 || HLT_PFMETNoMu120_PFMHTNoMu120)'
cuts['mini']['2016']['HT_Etmiss'] =         '(( nJet30>1 && ht<1000. && met_pt>250.) || ( nJet30>1 && ht>=1000. && met_pt>30.) || (nJet30==1 && met_pt>250.) )'
cuts['mini']['2016']['nJets>0'] =           '(nJet30>=1 && nJet30FailId == 0)'
cuts['mini']['2016']['deltaPhi>0.3'] =      '(deltaPhiMin > 0.3)'
cuts['mini']['2016']['HTmiss-Etmiss'] =     '(diffMetMht < 0.5*met_pt)'
cuts['mini']['2016']['MT2'] =               '((((mt2>200 && ht<1500) || (mt2>400&&ht>=1500)) &&(nJet30>1) )  ||  ( ht>200 && nJet30==1 )  )'
cuts['mini']['2016']['leptonVeto'] =        '(nMuons10==0 && nElectrons10==0)'
cuts['mini']['2016']['isoTrackVeto'] =      '(nPFLep5LowMT==0 && nPFHad10LowMT==0)'
cuts['mini']['2016']['isoTrackVetoLep'] =    '(nPFLep5LowMT==0)'
cuts['mini']['2016']['isoTrackVetoHad'] =    '(nPFHad10LowMT==0)'
cuts['mini']['2016']['1lepton'] =           '(nLepLowMT==1)' # should be equivalent to '(nMuons10==1 || nElectrons10==1 || nPFLep5LowMT==1 || nPFHad10LowMT==1 )'
cuts['mini']['2016']['incl1leptonCR'] =     '(ht>250. && met_pt>250  && nJet30>1 && mt2>200.)'
# cuts['mini'] that are dilepton CR specific
cuts['mini']['2016']['zll_SF_triggers'] =   '( lep_pdgId[0] == -lep_pdgId[1] && (HLT_DoubleMu || HLT_DoubleMu_NonIso || HLT_SingleMu_NonIso || HLT_DoubleEl || HLT_DoubleEl33 || HLT_Photon165_HE10))'
cuts['mini']['2016']['zll_OF_triggers'] =   '( lep_pdgId[0] != -lep_pdgId[1] && ( HLT_MuX_Ele12 || HLT_Mu8_EleX || HLT_Mu33_Ele33_NonIso || HLT_Mu30_Ele30_NonIso || HLT_Photon165_HE10 || HLT_SingleMu_NonIso ))'
cuts['mini']['2016']['zll_HT_Etmiss'] =     '((nJet30>1 && zll_ht<1000. && zll_met_pt>250.) || ( nJet30>1 && zll_ht>=1000. && zll_met_pt>30.) || (nJet30==1 && zll_met_pt>250.) )'
cuts['mini']['2016']['zll_deltaPhi>0.3'] =  '(zll_deltaPhiMin > 0.3)'
cuts['mini']['2016']['zll_HTmiss-Etmiss'] = '(zll_diffMetMht < 0.5*zll_met_pt)'
cuts['mini']['2016']['zll_MT2'] =           '((((zll_mt2>200 && zll_ht<1500) || (zll_mt2>400 && zll_ht>=1500)) &&(nJet30>1) )  ||  ( zll_ht>200 && nJet30==1 )  )'
cuts['mini']['2016']['ZmassPt'] =           '(zll_pt > 200. && fabs(zll_mass-91.19)<20)'
cuts['mini']['2016']['2leptons'] =          '(lep_pt[0]>100 && lep_pt[1]>30)'
# abs(myTree.lep_pdgId[0])==11 && myTree.lep_tightId[0]< 0.5
# abs(myTree.lep_pdgId[1])==11 && myTree.lep_tightId[1]< 0.5

cuts['mini']['2017'] = {}
cuts['mini']['2017']['noCut'] =             cuts['mini']['2016']['noCut']
cuts['mini']['2017']['isGolden'] =          cuts['mini']['2016']['isGolden']
cuts['mini']['2017']['goodVertex'] =        cuts['mini']['2016']['goodVertex']
cuts['mini']['2017']['filters'] =           '(Flag_goodVertices>0 && Flag_globalTightHalo2016Filter>0 && Flag_HBHENoiseFilter>0 && Flag_HBHENoiseIsoFilter>0 && Flag_EcalDeadCellTriggerPrimitiveFilter>0 && Flag_BadPFMuonFilter>0 && Flag_eeBadScFilter>0 && Flag_ecalBadCalibFilter>0)'
cuts['mini']['2017']['cleanings'] =         cuts['mini']['2016']['cleanings']
cuts['mini']['2017']['triggers'] =          '(HLT_PFMET120_PFMHT120 || HLT_PFHT1050 || HLT_PFHT500_PFMET100_PFMHT100 || HLT_PFJet500 || HLT_PFMETNoMu120_PFMHTNoMu120 || HLT_PFMETNoMu120_PFMHTNoMu120_PFHT60 )'
cuts['mini']['2017']['HT_Etmiss'] =         '(( nJet30>1 && ht<1200. && met_pt>250.) || ( nJet30>1 && ht>=1200. && met_pt>30.) || (nJet30==1 && met_pt>250.) )'
cuts['mini']['2017']['nJets>0'] =           cuts['mini']['2016']['nJets>0']
cuts['mini']['2017']['deltaPhi>0.3'] =      cuts['mini']['2016']['deltaPhi>0.3']
cuts['mini']['2017']['HTmiss-Etmiss'] =     cuts['mini']['2016']['HTmiss-Etmiss']
cuts['mini']['2017']['leptonVeto'] =        cuts['mini']['2016']['leptonVeto']
cuts['mini']['2017']['isoTrackVeto'] =      cuts['mini']['2016']['isoTrackVeto']
cuts['mini']['2017']['MT2'] =               cuts['mini']['2016']['MT2']
cuts['mini']['2017']['isoTrackVetoLep'] =   cuts['mini']['2016']['isoTrackVetoLep']
cuts['mini']['2017']['isoTrackVetoHad'] =   cuts['mini']['2016']['isoTrackVetoHad']
cuts['mini']['2017']['1lepton'] =           cuts['mini']['2016']['1lepton']
cuts['mini']['2017']['incl1leptonCR'] =     cuts['mini']['2016']['incl1leptonCR']
# cuts['mini'] that are dilepton CR specific
cuts['mini']['2017']['zll_SF_triggers'] =   cuts['mini']['2016']['zll_SF_triggers']
cuts['mini']['2017']['zll_OF_triggers'] =   cuts['mini']['2016']['zll_OF_triggers']
cuts['mini']['2017']['zll_HT_Etmiss'] =     '((nJet30>1 && zll_ht<1000. && zll_met_pt>250.) || ( nJet30>1 && zll_ht>=1000. && zll_met_pt>30.) || (nJet30==1 && zll_met_pt>250.) )'
cuts['mini']['2017']['zll_deltaPhi>0.3'] =  cuts['mini']['2016']['zll_deltaPhi>0.3']
cuts['mini']['2017']['zll_HTmiss-Etmiss'] = cuts['mini']['2016']['zll_HTmiss-Etmiss']
cuts['mini']['2017']['zll_MT2'] =           cuts['mini']['2016']['zll_MT2']
cuts['mini']['2017']['ZmassPt'] =           cuts['mini']['2016']['ZmassPt']
cuts['mini']['2017']['2leptons'] =          cuts['mini']['2016']['2leptons']

cuts['nano']={}
cuts['nano']['2017'] = {}
for key in cuts['mini']['2017'].keys():
  if 'trigger' not in key:
    cuts['nano']['2017'][key] = cuts['mini']['2017'][key]
  else:
    cuts['nano']['2017']['triggers'] = '(1)'
    cuts['nano']['2017']['zll_SF_triggers'] = '(1)'
    cuts['nano']['2017']['zll_OF_triggers'] = '(1)'


ordered_cutNames = {}
ordered_cutNames['SR'] =     ['noCut' ,'isGolden', 'goodVertex', 'filters', 'triggers', 'HT_Etmiss', 'nJets>0', 'deltaPhi>0.3', 'HTmiss-Etmiss', 'MT2', 'leptonVeto','isoTrackVetoLep', 'isoTrackVetoHad' ]
ordered_cutNames['CR1lep'] = ['noCut' ,'isGolden', 'goodVertex', 'filters', 'triggers', 'HT_Etmiss', 'nJets>0', 'deltaPhi>0.3', 'HTmiss-Etmiss', 'MT2', '1lepton', 'incl1leptonCR']
ordered_cutNames['CR2lep'] = ['noCut' ,'isGolden', 'goodVertex', 'filters', 'zll_SF_triggers', 'zll_HT_Etmiss', 'nJets>0', 'zll_deltaPhi>0.3', 'zll_HTmiss-Etmiss', 'zll_MT2', 'ZmassPt', '2leptons']
#ordered_cutNames['CR1ph'] =
#ordered_cutNames['SR'] = ['noCut', 'cleanings']
ordered_cutNames['SR_kinematic'] = ['noCut' , 'HT_Etmiss',  'deltaPhi>0.3', 'HTmiss-Etmiss', 'MT2', 'leptonVeto','isoTrackVetoLep', 'isoTrackVetoHad' ]
ordered_cutNames['CR1lep_kinematic'] = ['noCut', 'HT_Etmiss',  'deltaPhi>0.3', 'HTmiss-Etmiss', 'MT2', '1lepton', 'incl1leptonCR']
ordered_cutNames['CR2lep_kinematic'] = ['noCut', 'zll_HT_Etmiss',  'zll_deltaPhi>0.3', 'zll_HTmiss-Etmiss', 'zll_MT2', 'ZmassPt', '2leptons']


class SampleCutflow(object):
  def __init__(self, sampleName, fileNames,treeName, year, region, dataFormat):
    self.sampleName=sampleName
    self.fileNames=fileNames
    self.treeName=treeName
    self.year=year
    self.region=region
    self.dataFormat=dataFormat
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
      if not os.path.isfile(fileName): raise RuntimeError('{} not a file'.format(fileName))
      chain.Add(fileName)
    selection = '(1)'
    for cutName in ordered_cutNames[self.region]:
      selection = selection + '&&' + cuts[self.dataFormat][self.year][cutName]
      histo = TH1F('histo', 'histo', 10, 0., 14000)
      #print 'going to use this selection', selection
      chain.Project(histo.GetName(), 'met_pt', selection)
      self.cutflow[cutName] = float(histo.Integral())
      del histo

  def fillEfficiency(self):
    for i,cutName in enumerate(ordered_cutNames[self.region]):
      self.efficiency['abs'][cutName]=self.cutflow[cutName]/self.cutflow[ordered_cutNames[self.region][0]]*100
      if i==0:
        self.efficiency['rel'][cutName]=self.efficiency['abs'][cutName]
      else:
        den = self.cutflow[ordered_cutNames[self.region][i-1]]*100
        if den!=0:
          self.efficiency['rel'][cutName]=self.cutflow[cutName]/self.cutflow[ordered_cutNames[self.region][i-1]]*100
        else:
          self.efficiency['rel'][cutName]=0

  def stamp(self):
    print '\n\n\n'
    print 'Cutflow for sample={s}, year={y}, region={r}'.format(s=self.sampleName, y=self.year, r=self.region)
    print 'File name {n}'.format(n='\t\n'.join(self.fileNames))
    print ''

    print '{:>20} {:>20} {:>20} {:>20}'.format('Cut', 'Cutflow', 'Rel eff (%)', 'Abs eff (%)')
    for cutName in ordered_cutNames[self.region]:
      print '{:>20} {:>20} {:>20} {:>20}'.format(cutName, self.cutflow[cutName], '{:.2f}'.format(self.efficiency['rel'][cutName]), '{:.2f}'.format(self.efficiency['abs'][cutName]))
    print ''

if __name__ == "__main__":

  from argparse import ArgumentParser
  import os
  parser = ArgumentParser(description='', add_help=True)

  allowed_regionGroups = ordered_cutNames.keys()
  allowed_formats = ['mini', 'nano']

  parser.add_argument('-f', '--fileName', type=str, dest='fileName', help='full path to filename')
  parser.add_argument('-s', '--sampleName', type=str, dest='sampleName', help='name of sample, e.g. Wlv miniAOD')
  parser.add_argument('-t', '--treeName', type=str, dest='treeName', help='name of the tree', default='Events')
  parser.add_argument('-y', '--year', type=str, dest='year', help='year of data taking', default='2017')
  parser.add_argument('-g'  '--regionGroup', type=str, dest='regionGroup', help='cut series', default='SR', choices=allowed_regionGroups)
  parser.add_argument('-d'  '--dataFormat', type=str, dest='dataFormat', help='mini or nano', default='nano', choices=allowed_formats)

  options = parser.parse_args()

  myCutflow = SampleCutflow(sampleName=options.sampleName, fileNames=[options.fileName], treeName=options.treeName, year=options.year, region=options.regionGroup, dataFormat=options.dataFormat )
  myCutflow.fill()
  myCutflow.fillEfficiency()
  myCutflow.stamp()
