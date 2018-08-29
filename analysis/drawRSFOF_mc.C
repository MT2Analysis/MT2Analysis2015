{

//  TFile *_file0 = TFile::Open("/mnt/t3nfs01/data01/shome/mschoene/8_0_12_analysisPlayArea/src/mschoene_newBinning/analysis/EventYields_data_Run2016_12p9ifb_zurich2016/zllControlRegion/ZllPurityTrees.root");
//  TTree* SF_top = (TTree*) _file0->Get("Top/HT250toInf_j1toInf_b0toInf/tree_Top_HT250toInf_j1toInf_b0toInf")->Clone("SF_top");
//  TTree* SF_dy = (TTree*) _file0->Get("DYJets/HT250toInf_j1toInf_b0toInf/tree_DYJets_HT250toInf_j1toInf_b0toInf")->Clone("SF_dy");
  TChain* SF=new TChain("SF");
  //  SF->Add("/mnt/t3nfs01/data01/shome/mschoene//8_0_12_analysisPlayArea/src/mschoene_newBinning/analysis/EventYields_data_2016_36p2ifb/zllControlRegion/ZllPurityTrees.root/Top/HT250toInf_j1toInf_b0toInf/tree_Top_HT250toInf_j1toInf_b0toInf");
//  SF->Add("/mnt/t3nfs01/data01/shome/mschoene/8_0_12_analysisPlayArea/src/mschoene_newBinning/analysis/EventYields_data_2016_36p2ifb/zllControlRegion/ZllPurityTrees.root/Top/HT250toInf_j1toInf_b0toInf/tree_Top_HT250toInf_j1toInf_b0toInf");
//  SF->Add("/mnt/t3nfs01/data01/shome/mschoene/8_0_12_analysisPlayArea/src/mschoene_newBinning/analysis/EventYields_data_2016_36p2ifb/zllControlRegion/ZllPurityTrees.root/DYJets/HT250toInf_j1toInf_b0toInf/tree_DYJets_HT250toInf_j1toInf_b0toInf");
  SF->Add("/mnt/t3nfs01/data01/shome/mmasciov/8_0_12/src/myMT2Analysis/analysis/EventYields_data_2016_SnTMC_362ifb_noBtagSF//zllControlRegion/ZllPurityTrees.root/Top/HT250toInf_j1toInf_b0toInf/tree_Top_HT250toInf_j1toInf_b0toInf");
  SF->Add("/mnt/t3nfs01/data01/shome/mmasciov/8_0_12/src/myMT2Analysis/analysis/EventYields_data_2016_SnTMC_362ifb_noBtagSF//zllControlRegion/ZllPurityTrees.root/DYJets/HT250toInf_j1toInf_b0toInf/tree_DYJets_HT250toInf_j1toInf_b0toInf");

  //  SF->Add("/mnt/t3nfs01/data01/shome/mschoene//8_0_12_analysisPlayArea/src/mschoene_newBinning/analysis/EventYields_data_2016_36p2ifb/zllControlRegion/ZllPurityTrees.root/Top/HT250toInf_j1toInf_b0toInf/tree_Top_HT250toInf_j1toInf_b0toInf");

//  SF->Add("/mnt/t3nfs01/data01/shome/mschoene/8_0_12_analysisPlayArea/src/mschoene_newBinning/analysis/EventYields_data_Run2016_12p9ifb_zurich2016/zllControlRegion/ZllPurityTrees.root/Top/HT250toInf_j1toInf_b0toInf/tree_Top_HT250toInf_j1toInf_b0toInf");
//  SF->Add("/mnt/t3nfs01/data01/shome/mschoene/8_0_12_analysisPlayArea/src/mschoene_newBinning/analysis/EventYields_data_Run2016_12p9ifb_zurich2016/zllControlRegion/ZllPurityTrees.root/DYJets/HT250toInf_j1toInf_b0toInf/tree_DYJets_HT250toInf_j1toInf_b0toInf");

//  TFile *_file1 = TFile::Open("/mnt/t3nfs01/data01/shome/mschoene/8_0_12_analysisPlayArea/src/mschoene_newBinning/analysis/EventYields_data_Run2016_12p9ifb_zurich2016/zllControlRegion/ZllPurityTrees_of.root");
//  TTree* OF = (TTree*) _file1->Get("Top/HT250toInf_j1toInf_b0toInf/tree_Top_HT250toInf_j1toInf_b0toInf")->Clone("OF");

  TChain* OF=new TChain("OF");
  //  OF->Add("/mnt/t3nfs01/data01/shome/mschoene//8_0_12_analysisPlayArea/src/mschoene_newBinning/analysis/EventYields_data_2016_36p2ifb/zllControlRegion/ZllPurityTrees_of.root/Top/HT250toInf_j1toInf_b0toInf/tree_Top_HT250toInf_j1toInf_b0toInf");
  OF->Add("/mnt/t3nfs01/data01/shome/mmasciov/8_0_12/src/myMT2Analysis/analysis/EventYields_data_2016_SnTMC_362ifb_noBtagSF//zllControlRegion/ZllPurityTrees_of.root/Top/HT250toInf_j1toInf_b0toInf/tree_Top_HT250toInf_j1toInf_b0toInf");
  OF->Add("/mnt/t3nfs01/data01/shome/mmasciov/8_0_12/src/myMT2Analysis/analysis/EventYields_data_2016_SnTMC_362ifb_noBtagSF//zllControlRegion/ZllPurityTrees_of.root/DYJets/HT250toInf_j1toInf_b0toInf/tree_DYJets_HT250toInf_j1toInf_b0toInf");

//  OF->Add("/mnt/t3nfs01/data01/shome/mschoene/8_0_12_analysisPlayArea/src/mschoene_newBinning/analysis/EventYields_data_Run2016_12p9ifb_zurich2016/zllControlRegion/ZllPurityTrees_of.root/Top/HT250toInf_j1toInf_b0toInf/tree_Top_HT250toInf_j1toInf_b0toInf");
//  OF->Add("/mnt/t3nfs01/data01/shome/mschoene/8_0_12_analysisPlayArea/src/mschoene_newBinning/analysis/EventYields_data_Run2016_12p9ifb_zurich2016/zllControlRegion/ZllPurityTrees_of.root/DYJets/HT250toInf_j1toInf_b0toInf/tree_DYJets_HT250toInf_j1toInf_b0toInf");


  gDirectory->cd("");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  
  cout << SF->GetEntries() << endl;
//  cout << SF_top->GetEntries() << endl;
//  cout << SF_dy->GetEntries() << endl;
  cout << OF->GetEntries() << endl;

  float bins_mll[]={50,60,71.20,80,90,100,111.20,120,130,140,150,160,170,180,200};
  TH1D* mllSF = new TH1D("mllSF", "", 14, bins_mll);
  mllSF->Sumw2();

  TH1D* mllOF = new TH1D("mllOF", "", 14, bins_mll);
  mllOF->Sumw2();

  mllSF->SetLineColor(4);
  mllSF->SetMarkerColor(4);
  mllSF->SetMarkerStyle(4);
  mllOF->SetLineColor(2);
  mllOF->SetMarkerColor(2);
  mllOF->SetMarkerStyle(20);

  mllSF->GetXaxis()->SetTitle("M_{ll} [GeV]");
  mllOF->GetXaxis()->SetTitle("M_{ll} [GeV]");
  mllSF->GetYaxis()->SetTitle("Entries");
  mllOF->GetYaxis()->SetTitle("Entries");


  SF->Draw("Z_mass>>mllSF", "(Z_pt<=200 && Z_mass>=50 && (Z_mass<=71.19 || Z_mass>111.19) && ht>=250 && mt2>200 && lep_pt0>100. && lep_pt1>30.)*weight/HLT_weight", "goff");
  OF->Draw("Z_mass>>mllOF", "(Z_pt<=200 && Z_mass>=50 && (Z_mass<=71.19 || Z_mass>111.19) && ht>=250 && mt2>200 && lep_pt0>100. && lep_pt1>30.)*weight", "goff");

  mllSF->Scale(36.459);
  mllOF->Scale(36.459);
  
  double integral_sf, error_sf;
  double integral_of, error_of;
  
  integral_sf = mllSF->IntegralAndError(1,-1, error_sf);
  integral_of = mllOF->IntegralAndError(1,-1, error_of);
//  error_sf = TMath::Sqrt(integral_sf);
//  error_of = TMath::Sqrt(integral_of);
  
  cout << "R(SF/OF) = ("<<integral_sf<<" +- "<<error_sf<<")/("<<integral_of<<" +- "<<error_of<<")"<<endl;
  double RSFOF = integral_sf/integral_of;
  double err_R = TMath::Sqrt( (error_sf/integral_of)*(error_sf/integral_of) + (error_of*integral_sf/(integral_of*integral_of))*(error_of*integral_sf/(integral_of*integral_of)) );
  cout << "R(SF/OF) = "<< RSFOF << "+-" << err_R << endl;


  TCanvas* c1=new TCanvas("c1", "", 600, 600);
  c1->cd();
  mllSF->Draw("PE");
  mllOF->Draw("PE,same");
  

  TH1D* mllRSFOF = (TH1D*) mllSF->Clone("mllRSFOF");
  mllRSFOF->Divide(mllOF);

  mllRSFOF->SetMinimum(0.);
  mllRSFOF->SetMaximum(2.);
  
  mllRSFOF->GetYaxis()->SetRangeUser(0.,2.);
  mllRSFOF->GetYaxis()->SetTitle("R(SF/OF)");

  TCanvas* c1r=new TCanvas("c1r", "", 600, 600);
  c1r->cd();

  mllRSFOF->Draw("PE");
  //  mllRSFOF->Fit("pol0");
  
  TLine* lr = new TLine(50, (float) RSFOF, 200, (float)RSFOF);
  TLine* lrup = new TLine(50, (float)(RSFOF+err_R), 200, (float)(RSFOF+err_R));
  TLine* lrdn = new TLine(50, (float)(RSFOF-err_R), 200, (float)(RSFOF-err_R));

  lr->SetLineColor(1);
  lrup->SetLineColor(1);
  lrdn->SetLineColor(1);
  lrup->SetLineStyle(2);
  lrdn->SetLineStyle(2);
  
  lr->Draw("same");
  lrup->Draw("same");
  lrdn->Draw("same");

  //HT  
  float bins_ht[]={250,450,575,1000,2000};
  
  TH1D* htSF = new TH1D("htSF", "", 4, bins_ht);
  htSF->Sumw2();

  TH1D* htOF = new TH1D("htOF", "", 4, bins_ht);
  htOF->Sumw2();

  htSF->SetLineColor(4);
  htSF->SetMarkerColor(4);
  htSF->SetMarkerStyle(4);
  htOF->SetLineColor(2);
  htOF->SetMarkerColor(2);
  htOF->SetMarkerStyle(20);

  htSF->GetXaxis()->SetTitle("H_{T} [GeV]");
  htOF->GetXaxis()->SetTitle("H_{T} [GeV]");
  htSF->GetYaxis()->SetTitle("Entries");
  htOF->GetYaxis()->SetTitle("Entries");


  SF->Draw("ht>>htSF", "(Z_pt<=200 && Z_mass>=50 && (Z_mass<=71.19 || Z_mass>111.19) && ht>=250 && mt2>200 && lep_pt0>100. )*weight/HLT_weight", "goff");
  OF->Draw("ht>>htOF", "(Z_pt<=200 && Z_mass>=50 && (Z_mass<=71.19 || Z_mass>111.19) && ht>=250 && mt2>200 && lep_pt0>100. )*weight", "goff");
  
  htSF->Scale(36.459);
  htOF->Scale(36.459);
 

//  integral_sf = htSF->IntegralAndError(1,-1, error_sf);
//  integral_of = htOF->IntegralAndError(1,-1, error_of);
//  
//  cout << "R(SF/OF) = ("<<integral_sf<<" +- "<<error_sf<<")/("<<integral_of<<" +- "<<error_of<<")"<<endl;
//  RSFOF = integral_sf/integral_of;
//  err_R = TMath::Sqrt( (error_sf/integral_of)*(error_sf/integral_of) + (error_of*integral_sf/(integral_of*integral_of))*(error_of*integral_sf/(integral_of*integral_of)) );
//  cout << "R(SF/OF) = "<< RSFOF << "+-" << err_R << endl;

  TCanvas* c2=new TCanvas("c2", "", 600, 600);
  c2->cd();
  htSF->Draw("PE");
  htOF->Draw("PE,same");
  
  TH1D* htRSFOF = (TH1D*) htSF->Clone("htRSFOF");
  htRSFOF->Divide(htOF);

  htRSFOF->SetMinimum(0.);
  htRSFOF->SetMaximum(2.);
  
  htRSFOF->GetYaxis()->SetRangeUser(0.,2.);
  htRSFOF->GetYaxis()->SetTitle("R(SF/OF)");

  TCanvas* c2r=new TCanvas("c2r", "", 600, 600);
  c2r->cd();

  htRSFOF->Draw("PE");
  //  htRSFOF->Fit("pol0");

  TLine* lrht = new TLine(250, (float) RSFOF, 2000, (float)RSFOF);
  TLine* lrhtup = new TLine(250, (float)(RSFOF+err_R), 2000, (float)(RSFOF+err_R));
  TLine* lrhtdn = new TLine(250, (float)(RSFOF-err_R), 2000, (float)(RSFOF-err_R));

  lrht->SetLineColor(1);
  lrhtup->SetLineColor(1);
  lrhtdn->SetLineColor(1);
  lrhtup->SetLineStyle(2);
  lrhtdn->SetLineStyle(2);
  
  lrht->Draw("same");
  lrhtup->Draw("same");
  lrhtdn->Draw("same");


  //NJets
  float bins_nj[]={1,2,4,7,10};
  
  TH1D* njSF = new TH1D("njSF", "", 4, bins_nj);
  njSF->Sumw2();

  TH1D* njOF = new TH1D("njOF", "", 4, bins_nj);
  njOF->Sumw2();

  njSF->SetLineColor(4);
  njSF->SetMarkerColor(4);
  njSF->SetMarkerStyle(4);
  njOF->SetLineColor(2);
  njOF->SetMarkerColor(2);
  njOF->SetMarkerStyle(20);

  njSF->GetXaxis()->SetTitle("N_{j}");
  njOF->GetXaxis()->SetTitle("N_{j}");
  njSF->GetYaxis()->SetTitle("Entries");
  njOF->GetYaxis()->SetTitle("Entries");


  SF->Draw("nJets>>njSF", "(Z_pt<=200 && Z_mass>=50 && (Z_mass<=71.19 || Z_mass>111.19) && ht>=250 && mt2>200 && lep_pt0>100. )*weight/HLT_weight", "goff");
  OF->Draw("nJets>>njOF", "(Z_pt<=200 && Z_mass>=50 && (Z_mass<=71.19 || Z_mass>111.19) && ht>=250 && mt2>200 && lep_pt0>100. )*weight", "goff");
  
  njSF->Scale(36.459);
  njOF->Scale(36.459);
  

//  integral_sf = njSF->IntegralAndError(1,-1, error_sf);
//  integral_of = njOF->IntegralAndError(1,-1, error_of);
//  
//  cout << "R(SF/OF) = ("<<integral_sf<<" +- "<<error_sf<<")/("<<integral_of<<" +- "<<error_of<<")"<<endl;
//  RSFOF = integral_sf/integral_of;
//  err_R = TMath::Sqrt( (error_sf/integral_of)*(error_sf/integral_of) + (error_of*integral_sf/(integral_of*integral_of))*(error_of*integral_sf/(integral_of*integral_of)) );
//  cout << "R(SF/OF) = "<< RSFOF << "+-" << err_R << endl;

  TCanvas* c3=new TCanvas("c3", "", 600, 600);
  c3->cd();
  njSF->Draw("PE");
  njOF->Draw("PE,same");
  
  TH1D* njRSFOF = (TH1D*) njSF->Clone("njRSFOF");
  njRSFOF->Divide(njOF);

  njRSFOF->SetMinimum(0.);
  njRSFOF->SetMaximum(2.);
  
  njRSFOF->GetYaxis()->SetRangeUser(0.,2.);
  njRSFOF->GetYaxis()->SetTitle("R(SF/OF)");

  TCanvas* c3r=new TCanvas("c3r", "", 600, 600);
  c3r->cd();

  njRSFOF->Draw("PE");
  //  njRSFOF->Fit("pol0");

  TLine* lrnj = new TLine(1, (float) RSFOF, 10, (float)RSFOF);
  TLine* lrnjup = new TLine(1, (float)(RSFOF+err_R), 10, (float)(RSFOF+err_R));
  TLine* lrnjdn = new TLine(1, (float)(RSFOF-err_R), 10, (float)(RSFOF-err_R));

  lrnj->SetLineColor(1);
  lrnjup->SetLineColor(1);
  lrnjdn->SetLineColor(1);
  lrnjup->SetLineStyle(2);
  lrnjdn->SetLineStyle(2);
  
  lrnj->Draw("same");
  lrnjup->Draw("same");
  lrnjdn->Draw("same");


  //NB
  float bins_nb[]={0,1,2,3,4};
  
  TH1D* nbSF = new TH1D("nbSF", "", 4, bins_nb);
  nbSF->Sumw2();

  TH1D* nbOF = new TH1D("nbOF", "", 4, bins_nb);
  nbOF->Sumw2();

  nbSF->SetLineColor(4);
  nbSF->SetMarkerColor(4);
  nbSF->SetMarkerStyle(4);
  nbOF->SetLineColor(2);
  nbOF->SetMarkerColor(2);
  nbOF->SetMarkerStyle(20);

  nbSF->GetXaxis()->SetTitle("N_{b}");
  nbOF->GetXaxis()->SetTitle("N_{b}");
  nbSF->GetYaxis()->SetTitle("Entries");
  nbOF->GetYaxis()->SetTitle("Entries");


  SF->Draw("nBJets>>nbSF", "(Z_pt<=200 && Z_mass>=50 && (Z_mass<=71.19 || Z_mass>111.19) && ht>=250 && mt2>200 && lep_pt0>100. && lep_pt1>30.)*weight/HLT_weight", "goff");
  OF->Draw("nBJets>>nbOF", "(Z_pt<=200 && Z_mass>=50 && (Z_mass<=71.19 || Z_mass>111.19) && ht>=250 && mt2>200 && lep_pt0>100. && lep_pt1>30.)*weight", "goff");
  
  nbSF->Scale(36.459);
  nbOF->Scale(36.459);
  

//  integral_sf = nbSF->IntegralAndError(1,-1, error_sf);
//  integral_of = nbOF->IntegralAndError(1,-1, error_of);
//  
//  cout << "R(SF/OF) = ("<<integral_sf<<" +- "<<error_sf<<")/("<<integral_of<<" +- "<<error_of<<")"<<endl;
//  RSFOF = integral_sf/integral_of;
//  err_R = TMath::Sqrt( (error_sf/integral_of)*(error_sf/integral_of) + (error_of*integral_sf/(integral_of*integral_of))*(error_of*integral_sf/(integral_of*integral_of)) );
//  cout << "R(SF/OF) = "<< RSFOF << "+-" << err_R << endl;

  TCanvas* c4=new TCanvas("c4", "", 600, 600);
  c4->cd();
  nbSF->Draw("PE");
  nbOF->Draw("PE,same");
  
  TH1D* nbRSFOF = (TH1D*) nbSF->Clone("nbRSFOF");
  nbRSFOF->Divide(nbOF);

  nbRSFOF->SetMinimum(0.);
  nbRSFOF->SetMaximum(2.);
  
  nbRSFOF->GetYaxis()->SetRangeUser(0.,2.);
  nbRSFOF->GetYaxis()->SetTitle("R(SF/OF)");

  TCanvas* c4r=new TCanvas("c4r", "", 600, 600);
  c4r->cd();

  nbRSFOF->Draw("PE");
  //  nbRSFOF->Fit("pol0");

  TLine* lrnb   = new TLine(0, (float) RSFOF, 4, (float)RSFOF);
//  TLine* lrnbup = new TLine(0, (float)(RSFOF+err_R), 4, (float)(RSFOF+err_R));
//  TLine* lrnbdn = new TLine(0, (float)(RSFOF-err_R), 4, (float)(RSFOF-err_R));

  TLine* lrnbup = new TLine(0, (float)(RSFOF+0.15), 4, (float)(RSFOF+0.15));
  TLine* lrnbdn = new TLine(0, (float)(RSFOF-0.15), 4, (float)(RSFOF-0.15));

  lrnb->SetLineColor(1);
  lrnbup->SetLineColor(2);
  lrnbdn->SetLineColor(2);
  lrnbup->SetLineStyle(2);
  lrnbdn->SetLineStyle(2);
  
  lrnb->Draw("same");
  lrnbup->Draw("same");
  lrnbdn->Draw("same");


  //MT2
  float bins_mt2[]={200,300,400,600,1500};
  
  TH1D* mt2SF = new TH1D("mt2SF", "", 4, bins_mt2);
  mt2SF->Sumw2();

  TH1D* mt2OF = new TH1D("mt2OF", "", 4, bins_mt2);
  mt2OF->Sumw2();

  mt2SF->SetLineColor(4);
  mt2SF->SetMarkerColor(4);
  mt2SF->SetMarkerStyle(4);
  mt2OF->SetLineColor(2);
  mt2OF->SetMarkerColor(2);
  mt2OF->SetMarkerStyle(20);

  mt2SF->GetXaxis()->SetTitle("M_{T2}");
  mt2OF->GetXaxis()->SetTitle("M_{T2}");
  mt2SF->GetYaxis()->SetTitle("Entries");
  mt2OF->GetYaxis()->SetTitle("Entries");


  SF->Draw("mt2>>mt2SF", "(Z_pt<=200 && Z_mass>=50 && (Z_mass<=71.19 || Z_mass>111.19) && ht>=250 && mt2>200 && lep_pt0>100. )*weight/HLT_weight", "goff");
  OF->Draw("mt2>>mt2OF", "(Z_pt<=200 && Z_mass>=50 && (Z_mass<=71.19 || Z_mass>111.19) && ht>=250 && mt2>200 && lep_pt0>100. )*weight", "goff");
  
  mt2SF->Scale(36.459);
  mt2OF->Scale(36.459);
  

//  integral_sf = mt2SF->IntegralAndError(1,-1, error_sf);
//  integral_of = mt2OF->IntegralAndError(1,-1, error_of);
//  
//  cout << "R(SF/OF) = ("<<integral_sf<<" +- "<<error_sf<<")/("<<integral_of<<" +- "<<error_of<<")"<<endl;
//  RSFOF = integral_sf/integral_of;
//  err_R = TMath::Sqrt( (error_sf/integral_of)*(error_sf/integral_of) + (error_of*integral_sf/(integral_of*integral_of))*(error_of*integral_sf/(integral_of*integral_of)) );
//  cout << "R(SF/OF) = "<< RSFOF << "+-" << err_R << endl;

  TCanvas* c5=new TCanvas("c5", "", 600, 600);
  c5->cd();
  mt2SF->Draw("PE");
  mt2OF->Draw("PE,same");
  
  TH1D* mt2RSFOF = (TH1D*) mt2SF->Clone("mt2RSFOF");
  mt2RSFOF->Divide(mt2OF);

  mt2RSFOF->SetMinimum(0.);
  mt2RSFOF->SetMaximum(2.);
  
  mt2RSFOF->GetYaxis()->SetRangeUser(0.,2.);
  mt2RSFOF->GetYaxis()->SetTitle("R(SF/OF)");

  TCanvas* c5r=new TCanvas("c5r", "", 600, 600);
  c5r->cd();

  mt2RSFOF->Draw("PE");
  //  mt2RSFOF->Fit("pol0");

  TLine* lrmt2   = new TLine(200, (float) RSFOF,        1500., (float)RSFOF);
  TLine* lrmt2up = new TLine(200, (float)(RSFOF+err_R), 1500., (float)(RSFOF+err_R));
  TLine* lrmt2dn = new TLine(200, (float)(RSFOF-err_R), 1500., (float)(RSFOF-err_R));

  lrmt2->SetLineColor(1);
  lrmt2up->SetLineColor(1);
  lrmt2dn->SetLineColor(1);
  lrmt2up->SetLineStyle(2);
  lrmt2dn->SetLineStyle(2);
  
  lrmt2->Draw("same");
  lrmt2up->Draw("same");
  lrmt2dn->Draw("same");


  // dR(l1,l2)
  TH1D* drllSF = new TH1D("drllSF", "", 25,0,5);
  drllSF->Sumw2();

  TH1D* drllOF = new TH1D("drllOF", "", 25,0,5);
  drllOF->Sumw2();

  drllSF->SetLineColor(4);
  drllSF->SetMarkerColor(4);
  drllSF->SetMarkerStyle(4);
  drllOF->SetLineColor(2);
  drllOF->SetMarkerColor(2);
  drllOF->SetMarkerStyle(20);

  drllSF->GetXaxis()->SetTitle("#Delta R(l_{1}, l_{2})");
  drllOF->GetXaxis()->SetTitle("#Delta R(l_{1}, l_{2})");
  drllSF->GetYaxis()->SetTitle("Entries");
  drllOF->GetYaxis()->SetTitle("Entries");

  //  SF->Draw("TMath::Sqrt( (lep_eta0-lep_eta1)*(lep_eta0-lep_eta1) + ((lep_phi0-lep_phi1)>TMath::Pi() ? lep_phi0-lep_phi1-2*TMath::Pi() : ((lep_phi0-lep_phi1)<=-TMath::Pi() ? lep_phi0-lep_phi1+2*TMath::Pi() : lep_phi0-lep_phi1) )*((lep_phi0-lep_phi1)>TMath::Pi() ? lep_phi0-lep_phi1-2*TMath::Pi() : ((lep_phi0-lep_phi1)<=-TMath::Pi() ? lep_phi0-lep_phi1+2*TMath::Pi() : lep_phi0-lep_phi1) ) )>>drllSF", "(Z_pt<=200 && Z_mass>=50 && (Z_mass<=71.19 || Z_mass>111.19) && ht>=250 && mt2>200 && lep_pt0>100.  && (lep_eta0<1.4 || lep_eta0>1.6) && (lep_eta1<1.4 || lep_eta1>1.6))*weight/HLT_weight", "goff");
  //  OF->Draw("TMath::Sqrt( (lep_eta0-lep_eta1)*(lep_eta0-lep_eta1) + ((lep_phi0-lep_phi1)>TMath::Pi() ? lep_phi0-lep_phi1-2*TMath::Pi() : ((lep_phi0-lep_phi1)<=-TMath::Pi() ? lep_phi0-lep_phi1+2*TMath::Pi() : lep_phi0-lep_phi1) )*((lep_phi0-lep_phi1)>TMath::Pi() ? lep_phi0-lep_phi1-2*TMath::Pi() : ((lep_phi0-lep_phi1)<=-TMath::Pi() ? lep_phi0-lep_phi1+2*TMath::Pi() : lep_phi0-lep_phi1) ) )>>drllOF", "(Z_pt<=200 && Z_mass>=50 && (Z_mass<=71.19 || Z_mass>111.19) && ht>=250 && mt2>200 && lep_pt0>100.  && (lep_eta0<1.4 || lep_eta0>1.6) && (lep_eta1<1.4 || lep_eta1>1.6))*weight", "goff");

  SF->Draw("TMath::Sqrt( (lep_eta0-lep_eta1)*(lep_eta0-lep_eta1) + ((lep_phi0-lep_phi1)>TMath::Pi() ? lep_phi0-lep_phi1-2*TMath::Pi() : ((lep_phi0-lep_phi1)<=-TMath::Pi() ? lep_phi0-lep_phi1+2*TMath::Pi() : lep_phi0-lep_phi1) )*((lep_phi0-lep_phi1)>TMath::Pi() ? lep_phi0-lep_phi1-2*TMath::Pi() : ((lep_phi0-lep_phi1)<=-TMath::Pi() ? lep_phi0-lep_phi1+2*TMath::Pi() : lep_phi0-lep_phi1) ) )>>drllSF", "(Z_mass>50. && (Z_mass<=71.19 || Z_mass>111.19) && Z_pt<=200)*weight/HLT_weight", "goff");
  OF->Draw("TMath::Sqrt( (lep_eta0-lep_eta1)*(lep_eta0-lep_eta1) + ((lep_phi0-lep_phi1)>TMath::Pi() ? lep_phi0-lep_phi1-2*TMath::Pi() : ((lep_phi0-lep_phi1)<=-TMath::Pi() ? lep_phi0-lep_phi1+2*TMath::Pi() : lep_phi0-lep_phi1) )*((lep_phi0-lep_phi1)>TMath::Pi() ? lep_phi0-lep_phi1-2*TMath::Pi() : ((lep_phi0-lep_phi1)<=-TMath::Pi() ? lep_phi0-lep_phi1+2*TMath::Pi() : lep_phi0-lep_phi1) ) )>>drllOF", "(Z_mass>50. && (Z_mass<=71.19 || Z_mass>111.19) && Z_pt<=200)*weight", "goff");

  drllSF->Scale(36.459);
  drllOF->Scale(36.459);


  TCanvas* c6=new TCanvas("c6", "", 600, 600);
  c6->cd();
  drllSF->Draw("PE");
  drllOF->Draw("PE,same");

  TH1D* drllRSFOF = (TH1D*) drllSF->Clone("drllRSFOF");
  drllRSFOF->Divide(drllOF);

  drllRSFOF->SetMinimum(0.);
  drllRSFOF->SetMaximum(2.);

  drllRSFOF->GetYaxis()->SetRangeUser(0.,2.);
  drllRSFOF->GetYaxis()->SetTitle("R(SF/OF)");

  TCanvas* c6r=new TCanvas("c6r", "", 600, 600);
  c6r->cd();

  drllRSFOF->Draw("PE");

  TLine* lrdrll   = new TLine(0, (float) RSFOF,        5, (float)RSFOF);
  TLine* lrdrllup = new TLine(0, (float)(RSFOF+err_R), 5, (float)(RSFOF+err_R));
  TLine* lrdrlldn = new TLine(0, (float)(RSFOF-err_R), 5, (float)(RSFOF-err_R));

  lrdrll->SetLineColor(1);
  lrdrllup->SetLineColor(1);
  lrdrlldn->SetLineColor(1);
  lrdrllup->SetLineStyle(2);
  lrdrlldn->SetLineStyle(2);

  lrdrll->Draw("same");
  lrdrllup->Draw("same");
  lrdrlldn->Draw("same");


  //ZpT
  //  float bins_zpt[]={0,100,200,300,400,1500};
  float bins_zpt[]={0,200,1500};
  
  TH1D* zptSF = new TH1D("zptSF", "", 2, bins_zpt);
  zptSF->Sumw2();

  TH1D* zptOF = new TH1D("zptOF", "", 2, bins_zpt);
  zptOF->Sumw2();

  zptSF->SetLineColor(4);
  zptSF->SetMarkerColor(4);
  zptSF->SetMarkerStyle(4);
  zptOF->SetLineColor(2);
  zptOF->SetMarkerColor(2);
  zptOF->SetMarkerStyle(20);

  zptSF->GetXaxis()->SetTitle("p_{T}(ll)");
  zptOF->GetXaxis()->SetTitle("p_{T}(ll)");
  zptSF->GetYaxis()->SetTitle("Entries");
  zptOF->GetYaxis()->SetTitle("Entries");


  SF->Draw("Z_pt>>zptSF", "(Z_mass>=50 && (Z_mass<=71.19 || Z_mass>111.19) && ht>=250 && mt2>200 && lep_pt0>100. && lep_pt1>30.)*weight/HLT_weight", "goff");
  OF->Draw("Z_pt>>zptOF", "(Z_mass>=50 && (Z_mass<=71.19 || Z_mass>111.19) && ht>=250 && mt2>200 && lep_pt0>100. && lep_pt1>30.)*weight", "goff");
  
  zptSF->Scale(36.459);
  zptOF->Scale(36.459);
  

  integral_sf = zptSF->IntegralAndError(1,-1, error_sf);
  integral_of = zptOF->IntegralAndError(1,-1, error_of);
  
  cout << "R(SF/OF) = ("<<integral_sf<<" +- "<<error_sf<<")/("<<integral_of<<" +- "<<error_of<<")"<<endl;
  RSFOF = integral_sf/integral_of;
  err_R = TMath::Sqrt( (error_sf/integral_of)*(error_sf/integral_of) + (error_of*integral_sf/(integral_of*integral_of))*(error_of*integral_sf/(integral_of*integral_of)) );
  cout << "R(SF/OF) = "<< RSFOF << "+-" << err_R << endl;

  TCanvas* c7=new TCanvas("c7", "", 600, 600);
  c7->cd();
  zptSF->Draw("PE");
  zptOF->Draw("PE,same");
  
  TH1D* zptRSFOF = (TH1D*) zptSF->Clone("zptRSFOF");
  zptRSFOF->Divide(zptOF);

  zptRSFOF->SetMinimum(0.);
  zptRSFOF->SetMaximum(2.);
  
  zptRSFOF->GetYaxis()->SetRangeUser(0.,2.);
  zptRSFOF->GetYaxis()->SetTitle("R(SF/OF)");

  TCanvas* c7r=new TCanvas("c7r", "", 600, 600);
  c7r->cd();

  zptRSFOF->Draw("PE");
  //  zptRSFOF->Fit("pol0");

  TLine* lrzpt   = new TLine(0, (float) RSFOF,        1500., (float)RSFOF);
//  TLine* lrzptup = new TLine(0, (float)(RSFOF+err_R), 1500., (float)(RSFOF+err_R));
//  TLine* lrzptdn = new TLine(0, (float)(RSFOF-err_R), 1500., (float)(RSFOF-err_R));

  TLine* lrzptup = new TLine(0, (float)(RSFOF+0.15), 1500., (float)(RSFOF+0.15));
  TLine* lrzptdn = new TLine(0, (float)(RSFOF-0.15), 1500., (float)(RSFOF-0.15));

  lrzpt->SetLineColor(1);
//  lrzptup->SetLineColor(1);
//  lrzptdn->SetLineColor(1);
  lrzptup->SetLineColor(2);
  lrzptdn->SetLineColor(2);
  lrzptup->SetLineStyle(2);
  lrzptdn->SetLineStyle(2);
  
  lrzpt->Draw("same");
  lrzptup->Draw("same");
  lrzptdn->Draw("same");



  // Save As:

  c1 ->SaveAs("RSFOF_mc_noZwindow_36p2ifb_onlyTop_l1pt100/SFnOF_mll.pdf");
  c1r->SaveAs("RSFOF_mc_noZwindow_36p2ifb_onlyTop_l1pt100/RSFOF_mll.pdf");

  c2 ->SaveAs("RSFOF_mc_noZwindow_36p2ifb_onlyTop_l1pt100/SFnOF_ht.pdf");
  c2r->SaveAs("RSFOF_mc_noZwindow_36p2ifb_onlyTop_l1pt100/RSFOF_ht.pdf");

  c3 ->SaveAs("RSFOF_mc_noZwindow_36p2ifb_onlyTop_l1pt100/SFnOF_nj.pdf");
  c3r->SaveAs("RSFOF_mc_noZwindow_36p2ifb_onlyTop_l1pt100/RSFOF_nj.pdf");

  c4 ->SaveAs("RSFOF_mc_noZwindow_36p2ifb_onlyTop_l1pt100/SFnOF_nb.pdf");
  c4r->SaveAs("RSFOF_mc_noZwindow_36p2ifb_onlyTop_l1pt100/RSFOF_nb.pdf");

  c5 ->SaveAs("RSFOF_mc_noZwindow_36p2ifb_onlyTop_l1pt100/SFnOF_mt2.pdf");
  c5r->SaveAs("RSFOF_mc_noZwindow_36p2ifb_onlyTop_l1pt100/RSFOF_mt2.pdf");

  c6 ->SaveAs("RSFOF_mc_noZwindow_36p2ifb_onlyTop_l1pt100/SFnOF_drll.pdf");
  c6r->SaveAs("RSFOF_mc_noZwindow_36p2ifb_onlyTop_l1pt100/RSFOF_drll.pdf");

  c7 ->SaveAs("RSFOF_mc_noZwindow_36p2ifb_onlyTop_l1pt100/SFnOF_zpt.pdf");
  c7r->SaveAs("RSFOF_mc_noZwindow_36p2ifb_onlyTop_l1pt100/RSFOF_zpt.pdf");


}
