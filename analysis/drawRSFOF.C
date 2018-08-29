{

  // set the TStyle
  TStyle* style = new TStyle("DrawBaseStyle", "");
  gStyle->SetCanvasColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetOptStat(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->cd();
  // For the canvas:
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetCanvasDefH(600); //Height of canvas
  gStyle->SetCanvasDefW(600); //Width of canvas
  gStyle->SetCanvasDefX(0); //POsition on screen
  gStyle->SetCanvasDefY(0);
  // For the Pad:
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(kWhite);
  gStyle->SetPadGridX(false);
  gStyle->SetPadGridY(false);
  gStyle->SetGridColor(0);
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);
  // For the frame:
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(1);
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameFillStyle(0);
  gStyle->SetFrameLineColor(1);
  gStyle->SetFrameLineStyle(1);
  gStyle->SetFrameLineWidth(1);
  // Margins:
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.15);//0.13);
  gStyle->SetPadLeftMargin(0.15);//0.16);
  gStyle->SetPadRightMargin(0.05);//0.02);
  // For the Global title:
  gStyle->SetOptTitle(0);
  gStyle->SetTitleFont(42);
  gStyle->SetTitleColor(1);
  gStyle->SetTitleTextColor(1);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleFontSize(0.05);
  // For the axis titles:
  gStyle->SetTitleColor(1, "XYZ");
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetTitleSize(0.05, "XYZ");
  gStyle->SetTitleXOffset(1.15);//0.9);
  gStyle->SetTitleYOffset(1.5); // => 1.15 if exponents
  // For the axis labels:
  gStyle->SetLabelColor(1, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetLabelOffset(0.007, "XYZ");
  gStyle->SetLabelSize(0.045, "XYZ");
  // For the axis:
  gStyle->SetAxisColor(1, "XYZ");
  gStyle->SetStripDecimals(kTRUE);
  gStyle->SetTickLength(0.03, "XYZ");
  gStyle->SetNdivisions(510, "XYZ");
  gStyle->SetPadTickX(1); // To get tick marks on the opposite side of the frame
  gStyle->SetPadTickY(1);
  // for histograms:
  gStyle->SetHistLineColor(1);
  
  //  TFile( *_file0 = TFile::Open("/mnt/t3nfs01/data01/shome/mmasciov/8_0_12/src/myMT2Analysis/analysis/EventYields_data_2016_SnTMC_362ifb_noBtagSF/zllControlRegion/data.root");
  TFile *_file0 = TFile::Open("/mnt/t3nfs01/data01/shome/mschoene/8_0_12_analysisPlayArea/src/mschoene_newBinning/analysis/EventYields_dataETH_SnTMC_35p9ifb/zllControlRegion/data.root");
  TTree* SF = (TTree*) _file0->Get("data/HT250toInf_j1toInf_b0toInf/tree_data_HT250toInf_j1toInf_b0toInf")->Clone("SF");

  //   TFile *_file1 = TFile::Open("/mnt/t3nfs01/data01/shome/mmasciov/8_0_12/src/myMT2Analysis/analysis/EventYields_data_2016_SnTMC_362ifb_noBtagSF/zllControlRegion/data_of.root");
  TFile *_file1 = TFile::Open("/mnt/t3nfs01/data01/shome/mschoene/8_0_12_analysisPlayArea/src/mschoene_newBinning/analysis/EventYields_dataETH_SnTMC_35p9ifb/zllControlRegion//data_of.root");
  TTree* OF = (TTree*) _file1->Get("data_of/HT250toInf_j1toInf_b0toInf/tree_data_of_HT250toInf_j1toInf_b0toInf")->Clone("OF");

  gDirectory->cd("");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  cout << SF->GetEntries() << endl;
  cout << OF->GetEntries() << endl;

  float bins_mll[]={50,60,71.20,80,90,100,111.20,120,130,140,150,160,170,180,200};
  TH1D* mllSF = new TH1D("mllSF", "", 14, bins_mll);
  mllSF->Sumw2();

  TH1D* mllOF = new TH1D("mllOF", "", 14, bins_mll);
  mllOF->Sumw2();

  mllSF->SetLineColor(1);
  mllSF->SetMarkerColor(1);
  mllSF->SetMarkerStyle(20);
  mllOF->SetLineColor(2);
  mllOF->SetMarkerColor(2);
  mllOF->SetMarkerStyle(20);


  TH1D* hstat = new TH1D("hstat", "", 1, 0,1);
  hstat->SetLineColor(1);
  hstat->SetLineStyle(2);

  TH1D* hsyst = new TH1D("hyst", "", 1, 0,1);
  hsyst->SetLineColor(2);
  hsyst->SetLineStyle(2);

  TLegend* leg=new TLegend(0.6, 0.7+0.05, 0.89, 0.85+0.05);
  leg->AddEntry(mllSF, "Data", "PL");
  leg->AddEntry(hstat, "Stat. unc.", "L");
  leg->AddEntry(hsyst, "Syst. unc.", "L");
  //  leg->AddEntry(mllOF, "Syst. uncertainty", "L");
  
  leg->SetLineColor(0);
  leg->SetFillColor(0);

  TPaveText* label_top = new TPaveText(0.4,0.959,0.975,0.963, "brNDC");
  //  TPaveText* label_top = new TPaveText(0.4,0.953,0.975,0.975, "brNDC");
  label_top->SetBorderSize(0);
  label_top->SetFillColor(kWhite);
  label_top->SetTextSize(0.038);
  label_top->SetTextAlign(31); // align right
  label_top->SetTextFont(42);  // label_top->SetTextFont(62);
  label_top->AddText("35.9 fb^{-1} (13 TeV)");

//  TPaveText* label_top = new TPaveText(0.8+0.04,0.91+0.035,0.9+0.04,0.94+0.035, "brNDC");
//  label_top->SetBorderSize(0);
//  label_top->SetFillColor(kWhite);
//  label_top->SetTextSize(0.038);
//  label_top->SetTextAlign(31); // align right
//  label_top->SetTextFont(62);
//  label_top->AddText("35.9 fb^{-1} (13 TeV)");

  TPaveText* label_cms = new TPaveText(0.143,0.96,0.27,0.965, "brNDC");
  label_cms->SetBorderSize(0);
  label_cms->SetFillColor(kWhite);
  label_cms->SetTextSize(0.042);
  label_cms->SetTextAlign(11); // align left
  label_cms->SetTextFont(62);
  label_cms->AddText("CMS");

//   TPaveText* label_cms = new TPaveText(0.1+0.05,0.91+0.035,0.4+0.04,0.94+0.035, "brNDC");
//  label_cms->SetBorderSize(0);
//  label_cms->SetFillColor(kWhite);
//  label_cms->SetTextSize(0.038);
//  label_cms->SetTextAlign(11); // align left
//  label_cms->SetTextFont(62);
//  label_cms->AddText("CMS Preliminary");


  TPaveText* label_nj = new TPaveText(0.1+0.06,0.7,0.3+0.06,0.89, "brNDC");
  label_nj->SetBorderSize(0);
  label_nj->SetFillColor(kWhite);
  label_nj->SetTextSize(0.038);
  label_nj->SetTextAlign(31);                                                                                                                                                       
  label_nj->SetTextFont(62);
  label_nj->AddText("N_{j} #geq 2");

  mllSF->GetXaxis()->SetTitle("M_{ll} [GeV]");
  mllOF->GetXaxis()->SetTitle("M_{ll} [GeV]");
  mllSF->GetYaxis()->SetTitle("Entries");
  mllOF->GetYaxis()->SetTitle("Entries");

  SF->Draw("Z_mass>>mllSF", "Z_pt<=200 && Z_mass>=50 && (Z_mass<=71.19 || Z_mass>111.19) && ht>=250 && mt2>200. && lep_pt0>100. && lep_pt1>30. && (nJets>1 || (mt2>250.))", "goff");
  OF->Draw("Z_mass>>mllOF", "Z_pt<=200 && Z_mass>=50 && (Z_mass<=71.19 || Z_mass>111.19) && ht>=250 && mt2>200. && lep_pt0>100. && lep_pt1>30. && (nJets>1 || (mt2>250.))", "goff");
  
  double integral_sf, error_sf;
  double integral_of, error_of;
  
  integral_sf = mllSF->IntegralAndError(1,-1, error_sf);
  integral_of = mllOF->IntegralAndError(1,-1, error_of);
  
  cout << "R^{SF/OF} = ("<<integral_sf<<" +- "<<error_sf<<")/("<<integral_of<<" +- "<<error_of<<")"<<endl;
  double RSFOF = integral_sf/integral_of;
  double err_R = TMath::Sqrt( (error_sf/integral_of)*(error_sf/integral_of) + (error_of*integral_sf/(integral_of*integral_of))*(error_of*integral_sf/(integral_of*integral_of)) );
  cout << "R^{SF/OF} = "<< RSFOF << "+-" << err_R << endl;


  SF->Draw("Z_mass>>mllSF", "Z_pt<=200 && Z_mass>=50 && (Z_mass<=71.19 || Z_mass>111.19) && ht>=250 && mt2>200. && lep_pt0>100. && lep_pt1>30. && (nJets>1 || (mt2>250.)) &&  diffMetMht/met<0.5 && deltaPhiMin >0.3", "goff");
  OF->Draw("Z_mass>>mllOF", "Z_pt<=200 && Z_mass>=50 && (Z_mass<=71.19 || Z_mass>111.19) && ht>=250 && mt2>200. && lep_pt0>100. && lep_pt1>30. && (nJets>1 || (mt2>250.)) &&  diffMetMht/met<0.5 && deltaPhiMin >0.3", "goff");

  integral_sf = mllSF->IntegralAndError(1,-1, error_sf);
  integral_of = mllOF->IntegralAndError(1,-1, error_of);
  
  cout << "R^{SF/OF} = ("<<integral_sf<<" +- "<<error_sf<<")/("<<integral_of<<" +- "<<error_of<<")"<<endl;
  RSFOF = integral_sf/integral_of;
  err_R = TMath::Sqrt( (error_sf/integral_of)*(error_sf/integral_of) + (error_of*integral_sf/(integral_of*integral_of))*(error_of*integral_sf/(integral_of*integral_of)) );
  cout << "R^{SF/OF} = "<< RSFOF << "+-" << err_R << endl;



  TCanvas* c1=new TCanvas("c1", "", 600, 600);
  c1->cd();
  mllSF->Draw("PE");
  mllOF->Draw("PE,same");
  

  TH1D* mllRSFOF = (TH1D*) mllSF->Clone("mllRSFOF");
  mllRSFOF->Divide(mllOF);

  mllRSFOF->SetMinimum(0.);
  mllRSFOF->SetMaximum(2.);
  
  mllRSFOF->GetYaxis()->SetRangeUser(0.,2.);
  mllRSFOF->GetYaxis()->SetTitle("R^{SF/OF}");
  mllRSFOF->GetYaxis()->SetTitleOffset(1.4);

//  TCanvas* c1r=new TCanvas("c1r", "", 600, 600);
//  c1r->cd();
//
//  mllRSFOF->Draw("PE");
//  //  mllRSFOF->Fit("pol0");
  
  TLine* lr = new TLine(50, (float) RSFOF, 200, (float)RSFOF);
  TLine* lrup = new TLine(50, (float)(RSFOF+err_R), 200, (float)(RSFOF+err_R));
  TLine* lrdn = new TLine(50, (float)(RSFOF-err_R), 200, (float)(RSFOF-err_R));
  TLine* lrup15 = new TLine(50, (float)(RSFOF+0.15), 200, (float)(RSFOF+0.15));
  TLine* lrdn15 = new TLine(50, (float)(RSFOF-0.15), 200, (float)(RSFOF-0.15));

//  lr->SetLineColor(1);
//  lrup->SetLineColor(1);
//  lrdn->SetLineColor(1);
//  lrup->SetLineStyle(2);
//  lrdn->SetLineStyle(2);
//  
//  lrup15->SetLineColor(2);
//  lrdn15->SetLineColor(2);
//  lrup15->SetLineStyle(2);
//  lrdn15->SetLineStyle(2);
//  
//  lr->Draw("same");
//  lrup->Draw("same");
//  lrdn->Draw("same");
//  lrup15->Draw("same");
//  lrdn15->Draw("same");

  //HT  
  float bins_ht[]={250,450,575,1000,2000};
  
  TH1D* htSF = new TH1D("htSF", "", 4, bins_ht);
  htSF->Sumw2();

  TH1D* htOF = new TH1D("htOF", "", 4, bins_ht);
  htOF->Sumw2();

  htSF->SetLineColor(1);
  htSF->SetMarkerColor(1);
  htSF->SetMarkerStyle(20);
  htOF->SetLineColor(2);
  htOF->SetMarkerColor(2);
  htOF->SetMarkerStyle(20);

  htSF->GetXaxis()->SetTitle("H_{T} [GeV]");
  htOF->GetXaxis()->SetTitle("H_{T} [GeV]");
  htSF->GetYaxis()->SetTitle("Entries");
  htOF->GetYaxis()->SetTitle("Entries");


  SF->Draw("ht>>htSF", "Z_pt<=200 && Z_mass>=50 && (Z_mass<=71.19 || Z_mass>111.19) && ht>=250 && mt2>200. && lep_pt0>100. && lep_pt1>30.  && (nJets>1 || (mt2>250.))", "goff");
  OF->Draw("ht>>htOF", "Z_pt<=200 && Z_mass>=50 && (Z_mass<=71.19 || Z_mass>111.19) && ht>=250 && mt2>200. && lep_pt0>100. && lep_pt1>30.  && (nJets>1 || (mt2>250.))", "goff");
  
  integral_sf = htSF->IntegralAndError(1,-1, error_sf);
  integral_of = htOF->IntegralAndError(1,-1, error_of);
  
  cout << "R^{SF/OF} = ("<<integral_sf<<" +- "<<error_sf<<")/("<<integral_of<<" +- "<<error_of<<")"<<endl;
  RSFOF = integral_sf/integral_of;
  err_R = TMath::Sqrt( (error_sf/integral_of)*(error_sf/integral_of) + (error_of*integral_sf/(integral_of*integral_of))*(error_of*integral_sf/(integral_of*integral_of)) );
  cout << "R^{SF/OF} = "<< RSFOF << "+-" << err_R << endl;

  TCanvas* c2=new TCanvas("c2", "", 600, 600);
  c2->cd();
  htSF->Draw("PE");
  htOF->Draw("PE,same");
  
  TH1D* htRSFOF = (TH1D*) htSF->Clone("htRSFOF");
  htRSFOF->Divide(htOF);

  htRSFOF->SetMinimum(0.);
  htRSFOF->SetMaximum(2.);
  
  htRSFOF->GetYaxis()->SetRangeUser(0.,2.);
  htRSFOF->GetYaxis()->SetTitle("R^{SF/OF}");
  htRSFOF->GetYaxis()->SetTitleOffset(1.4);

  TCanvas* c2r=new TCanvas("c2r", "", 600, 600);
  c2r->cd();

  htRSFOF->Draw("PE");
  //  htRSFOF->Fit("pol0");

  TLine* lrht = new TLine(250, (float) RSFOF, 2000, (float)RSFOF);
  TLine* lrhtup = new TLine(250, (float)(RSFOF+err_R), 2000, (float)(RSFOF+err_R));
  TLine* lrhtdn = new TLine(250, (float)(RSFOF-err_R), 2000, (float)(RSFOF-err_R));
  TLine* lrhtup15 = new TLine(250, (float)(RSFOF+0.15), 2000, (float)(RSFOF+0.15));
  TLine* lrhtdn15 = new TLine(250, (float)(RSFOF-0.15), 2000, (float)(RSFOF-0.15));

  lrht->SetLineColor(1);
  lrhtup->SetLineColor(1);
  lrhtdn->SetLineColor(1);
  lrhtup->SetLineStyle(2);
  lrhtdn->SetLineStyle(2);
  
  lrhtup15->SetLineColor(2);
  lrhtdn15->SetLineColor(2);
  lrhtup15->SetLineStyle(2);
  lrhtdn15->SetLineStyle(2);
  
  lrht->Draw("same");
  lrhtup->Draw("same");
  lrhtdn->Draw("same");

  lrhtup15->Draw("same");
  lrhtdn15->Draw("same");

  leg->Draw("same");
  label_top->Draw("same");
  label_cms->Draw("same");
  gPad->RedrawAxis();

  //NJets
  float bins_nj[]={1,2,4,7,10};
  
  TH1D* njSF = new TH1D("njSF", "", 4, bins_nj);
  njSF->Sumw2();

  TH1D* njOF = new TH1D("njOF", "", 4, bins_nj);
  njOF->Sumw2();

  njSF->SetLineColor(1);
  njSF->SetMarkerColor(1);
  njSF->SetMarkerStyle(20);
  njOF->SetLineColor(2);
  njOF->SetMarkerColor(2);
  njOF->SetMarkerStyle(20);

  njSF->GetXaxis()->SetTitle("N_{j}");
  njOF->GetXaxis()->SetTitle("N_{j}");
  njSF->GetYaxis()->SetTitle("Entries");
  njOF->GetYaxis()->SetTitle("Entries");


  SF->Draw("nJets>>njSF", "Z_pt<=200 && Z_mass>=50 && (Z_mass<=71.19 || Z_mass>111.19) && ht>=250 && mt2>200. && lep_pt0>100. && lep_pt1>30. && (nJets>1 || mt2>250)", "goff");
  OF->Draw("nJets>>njOF", "Z_pt<=200 && Z_mass>=50 && (Z_mass<=71.19 || Z_mass>111.19) && ht>=250 && mt2>200. && lep_pt0>100. && lep_pt1>30. && (nJets>1 || mt2>250)", "goff");
  
//  integral_sf = njSF->IntegralAndError(1,-1, error_sf);
//  integral_of = njOF->IntegralAndError(1,-1, error_of);
//  
//  cout << "R^{SF/OF} = ("<<integral_sf<<" +- "<<error_sf<<")/("<<integral_of<<" +- "<<error_of<<")"<<endl;
//  RSFOF = integral_sf/integral_of;
//  err_R = TMath::Sqrt( (error_sf/integral_of)*(error_sf/integral_of) + (error_of*integral_sf/(integral_of*integral_of))*(error_of*integral_sf/(integral_of*integral_of)) );
//  cout << "R^{SF/OF} = "<< RSFOF << "+-" << err_R << endl;

  TCanvas* c3=new TCanvas("c3", "", 600, 600);
  c3->cd();
  njSF->Draw("PE");
  njOF->Draw("PE,same");
  
  TH1D* njRSFOF = (TH1D*) njSF->Clone("njRSFOF");
  njRSFOF->Divide(njOF);

  njRSFOF->SetMinimum(0.);
  njRSFOF->SetMaximum(2.);
  
  njRSFOF->GetYaxis()->SetRangeUser(0.,2.);
  njRSFOF->GetYaxis()->SetTitle("R^{SF/OF}");
  njRSFOF->GetYaxis()->SetTitleOffset(1.4);

  TCanvas* c3r=new TCanvas("c3r", "", 600, 600);
  c3r->cd();

  njRSFOF->Draw("PE");
  //  njRSFOF->Fit("pol0");

  TLine* lrnj = new TLine(1, (float) RSFOF, 10, (float)RSFOF);
  TLine* lrnjup = new TLine(1, (float)(RSFOF+err_R), 10, (float)(RSFOF+err_R));
  TLine* lrnjdn = new TLine(1, (float)(RSFOF-err_R), 10, (float)(RSFOF-err_R));

  TLine* lrnjup15 = new TLine(1, (float)(RSFOF+0.15), 10, (float)(RSFOF+0.15));
  TLine* lrnjdn15 = new TLine(1, (float)(RSFOF-0.15), 10, (float)(RSFOF-0.15));

  lrnj->SetLineColor(1);
  lrnjup->SetLineColor(1);
  lrnjdn->SetLineColor(1);
  lrnjup->SetLineStyle(2);
  lrnjdn->SetLineStyle(2);
  
  lrnjup15->SetLineColor(2);
  lrnjdn15->SetLineColor(2);
  lrnjup15->SetLineStyle(2);
  lrnjdn15->SetLineStyle(2);
  
  lrnj->Draw("same");
  lrnjup->Draw("same");
  lrnjdn->Draw("same");
  lrnjup15->Draw("same");
  lrnjdn15->Draw("same");

  leg->Draw("same");
  label_top->Draw("same");
  label_cms->Draw("same");
  gPad->RedrawAxis();


  //NB
  float bins_nb[]={0,1,2,3,4};
  
  TH1D* nbSF = new TH1D("nbSF", "", 4, bins_nb);
  nbSF->Sumw2();

  TH1D* nbOF = new TH1D("nbOF", "", 4, bins_nb);
  nbOF->Sumw2();

  nbSF->SetLineColor(1);
  nbSF->SetMarkerColor(1);
  nbSF->SetMarkerStyle(20);
  nbOF->SetLineColor(2);
  nbOF->SetMarkerColor(2);
  nbOF->SetMarkerStyle(20);

  nbSF->GetXaxis()->SetTitle("N_{b}");
  nbOF->GetXaxis()->SetTitle("N_{b}");
  nbSF->GetYaxis()->SetTitle("Entries");
  nbOF->GetYaxis()->SetTitle("Entries");


  SF->Draw("nBJets>>nbSF", "Z_pt<=200 && Z_mass>=50 && (Z_mass<=71.19 || Z_mass>111.19) && ht>=250 && mt2>200. && lep_pt0>100. && lep_pt1>30. && (nJets>1 || (mt2>250.)) ", "goff");
  OF->Draw("nBJets>>nbOF", "Z_pt<=200 && Z_mass>=50 && (Z_mass<=71.19 || Z_mass>111.19) && ht>=250 && mt2>200. && lep_pt0>100. && lep_pt1>30. && (nJets>1 || (mt2>250.)) ", "goff");
  
//  integral_sf = nbSF->IntegralAndError(1,-1, error_sf);
//  integral_of = nbOF->IntegralAndError(1,-1, error_of);
//  
//  cout << "R^{SF/OF} = ("<<integral_sf<<" +- "<<error_sf<<")/("<<integral_of<<" +- "<<error_of<<")"<<endl;
//  RSFOF = integral_sf/integral_of;
//  err_R = TMath::Sqrt( (error_sf/integral_of)*(error_sf/integral_of) + (error_of*integral_sf/(integral_of*integral_of))*(error_of*integral_sf/(integral_of*integral_of)) );
//  cout << "R^{SF/OF} = "<< RSFOF << "+-" << err_R << endl;

  TCanvas* c4=new TCanvas("c4", "", 600, 600);
  c4->cd();
  nbSF->Draw("PE");
  nbOF->Draw("PE,same");
  
  TH1D* nbRSFOF = (TH1D*) nbSF->Clone("nbRSFOF");
  nbRSFOF->Divide(nbOF);

  nbRSFOF->SetMinimum(0.);
  nbRSFOF->SetMaximum(2.);
  
  nbRSFOF->GetYaxis()->SetRangeUser(0.,2.);
  nbRSFOF->GetYaxis()->SetTitle("R^{SF/OF}");
  nbRSFOF->GetYaxis()->SetTitleOffset(1.4);

  TCanvas* c4r=new TCanvas("c4r", "", 600, 600);
  c4r->cd();

  nbRSFOF->Draw("PE");
  //  nbRSFOF->Fit("pol0");

  TLine* lrnb   = new TLine(0, (float) RSFOF, 4, (float)RSFOF);
  TLine* lrnbup = new TLine(0, (float)(RSFOF+err_R), 4, (float)(RSFOF+err_R));
  TLine* lrnbdn = new TLine(0, (float)(RSFOF-err_R), 4, (float)(RSFOF-err_R));

  TLine* lrnbup15 = new TLine(0, (float)(RSFOF+0.15), 4, (float)(RSFOF+0.15));
  TLine* lrnbdn15 = new TLine(0, (float)(RSFOF-0.15), 4, (float)(RSFOF-0.15));

  lrnb->SetLineColor(1);
  lrnbup->SetLineColor(1);
  lrnbdn->SetLineColor(1);
  lrnbup->SetLineStyle(2);
  lrnbdn->SetLineStyle(2);
  
  lrnbup15->SetLineColor(2);
  lrnbdn15->SetLineColor(2);
  lrnbup15->SetLineStyle(2);
  lrnbdn15->SetLineStyle(2);
  
  lrnb->Draw("same");
  lrnbup->Draw("same");
  lrnbdn->Draw("same");
  lrnbup15->Draw("same");
  lrnbdn15->Draw("same");

  leg->Draw("same");
  label_top->Draw("same");
  label_cms->Draw("same");
  gPad->RedrawAxis();

  //MT2
  float bins_mt2[]={200,300,400,1500};
  
  TH1D* mt2SF = new TH1D("mt2SF", "", 3, bins_mt2);
  mt2SF->Sumw2();

  TH1D* mt2OF = new TH1D("mt2OF", "", 3, bins_mt2);
  mt2OF->Sumw2();

  mt2SF->SetLineColor(1);
  mt2SF->SetMarkerColor(1);
  mt2SF->SetMarkerStyle(20);
  mt2OF->SetLineColor(2);
  mt2OF->SetMarkerColor(2);
  mt2OF->SetMarkerStyle(20);

  mt2SF->GetXaxis()->SetTitle("M_{T2} [GeV]");
  mt2OF->GetXaxis()->SetTitle("M_{T2} [GeV]");
  mt2SF->GetYaxis()->SetTitle("Entries");
  mt2OF->GetYaxis()->SetTitle("Entries");


  SF->Draw("mt2>>mt2SF", "Z_pt<=200 && Z_mass>=50 && (Z_mass<=71.19 || Z_mass>111.19) && ht>=250 && mt2>200. && lep_pt0>100. && lep_pt1>30. && (nJets>1 || (mt2>250.)) && nJets>1", "goff");
  OF->Draw("mt2>>mt2OF", "Z_pt<=200 && Z_mass>=50 && (Z_mass<=71.19 || Z_mass>111.19) && ht>=250 && mt2>200. && lep_pt0>100. && lep_pt1>30. && (nJets>1 || (mt2>250.)) && nJets>1", "goff");
  
//  integral_sf = mt2SF->IntegralAndError(1,-1, error_sf);
//  integral_of = mt2OF->IntegralAndError(1,-1, error_of);
//  
//  cout << "R^{SF/OF} = ("<<integral_sf<<" +- "<<error_sf<<")/("<<integral_of<<" +- "<<error_of<<")"<<endl;
//  RSFOF = integral_sf/integral_of;
//  err_R = TMath::Sqrt( (error_sf/integral_of)*(error_sf/integral_of) + (error_of*integral_sf/(integral_of*integral_of))*(error_of*integral_sf/(integral_of*integral_of)) );
//  cout << "R^{SF/OF} = "<< RSFOF << "+-" << err_R << endl;

  TCanvas* c5=new TCanvas("c5", "", 600, 600);
  c5->cd();
  mt2SF->Draw("PE");
  mt2OF->Draw("PE,same");
  
  TH1D* mt2RSFOF = (TH1D*) mt2SF->Clone("mt2RSFOF");
  mt2RSFOF->Divide(mt2OF);

  mt2RSFOF->SetMinimum(0.);
  mt2RSFOF->SetMaximum(2.);
  
  mt2RSFOF->GetYaxis()->SetRangeUser(0.,2.);
  mt2RSFOF->GetYaxis()->SetTitle("R^{SF/OF}");
  mt2RSFOF->GetYaxis()->SetTitleOffset(1.4);

  TCanvas* c5r=new TCanvas("c5r", "", 600, 600);
  c5r->cd();

  mt2RSFOF->Draw("PE");
  //  mt2RSFOF->Fit("pol0");

  TLine* lrmt2   = new TLine(200, (float) RSFOF,        1500., (float)RSFOF);
  TLine* lrmt2up = new TLine(200, (float)(RSFOF+err_R), 1500., (float)(RSFOF+err_R));
  TLine* lrmt2dn = new TLine(200, (float)(RSFOF-err_R), 1500., (float)(RSFOF-err_R));

  TLine* lrmt2up15 = new TLine(200, (float)(RSFOF+0.15), 1500., (float)(RSFOF+0.15));
  TLine* lrmt2dn15 = new TLine(200, (float)(RSFOF-0.15), 1500., (float)(RSFOF-0.15));

  lrmt2->SetLineColor(1);
  lrmt2up->SetLineColor(1);
  lrmt2dn->SetLineColor(1);
  lrmt2up->SetLineStyle(2);
  lrmt2dn->SetLineStyle(2);
  
  lrmt2up15->SetLineColor(2);
  lrmt2dn15->SetLineColor(2);
  lrmt2up15->SetLineStyle(2);
  lrmt2dn15->SetLineStyle(2);
  
  lrmt2->Draw("same");
  lrmt2up->Draw("same");
  lrmt2dn->Draw("same");
  lrmt2up15->Draw("same");
  lrmt2dn15->Draw("same");

  leg->Draw("same");
  label_top->Draw("same");
  label_cms->Draw("same");
  label_nj->Draw("same");
  gPad->RedrawAxis();

//  // dR(l1,l2)
//  TH1D* drllSF = new TH1D("drllSF", "", 25,0,5);
//  drllSF->Sumw2();
//
//  TH1D* drllOF = new TH1D("drllOF", "", 25,0,5);
//  drllOF->Sumw2();
//
//  drllSF->SetLineColor(1);
//  drllSF->SetMarkerColor(1);
//  drllSF->SetMarkerStyle(20);
//  drllOF->SetLineColor(2);
//  drllOF->SetMarkerColor(2);
//  drllOF->SetMarkerStyle(20);
//
//  drllSF->GetXaxis()->SetTitle("#Delta R(l_{1}, l_{2})");
//  drllOF->GetXaxis()->SetTitle("#Delta R(l_{1}, l_{2})");
//  drllSF->GetYaxis()->SetTitle("Entries");
//  drllOF->GetYaxis()->SetTitle("Entries");
//
//  SF->Draw("TMath::Sqrt( (lep_eta0-lep_eta1)*(lep_eta0-lep_eta1) + ((lep_phi0-lep_phi1)>TMath::Pi() ? lep_phi0-lep_phi1-2*TMath::Pi() : ((lep_phi0-lep_phi1)<=-TMath::Pi() ? lep_phi0-lep_phi1+2*TMath::Pi() : lep_phi0-lep_phi1) )*((lep_phi0-lep_phi1)>TMath::Pi() ? lep_phi0-lep_phi1-2*TMath::Pi() : ((lep_phi0-lep_phi1)<=-TMath::Pi() ? lep_phi0-lep_phi1+2*TMath::Pi() : lep_phi0-lep_phi1) ) )>>drllSF", "Z_pt<=200 && Z_mass>=50 && (Z_mass<=71.19 || Z_mass>111.19) && ht>=250 && mt2>200. && lep_pt0>100. && lep_pt1>30. && (nJets>1 || (mt2>250.)) ", "goff");
//  OF->Draw("TMath::Sqrt( (lep_eta0-lep_eta1)*(lep_eta0-lep_eta1) + ((lep_phi0-lep_phi1)>TMath::Pi() ? lep_phi0-lep_phi1-2*TMath::Pi() : ((lep_phi0-lep_phi1)<=-TMath::Pi() ? lep_phi0-lep_phi1+2*TMath::Pi() : lep_phi0-lep_phi1) )*((lep_phi0-lep_phi1)>TMath::Pi() ? lep_phi0-lep_phi1-2*TMath::Pi() : ((lep_phi0-lep_phi1)<=-TMath::Pi() ? lep_phi0-lep_phi1+2*TMath::Pi() : lep_phi0-lep_phi1) ) )>>drllOF", "Z_pt<=200 && Z_mass>=50 && (Z_mass<=71.19 || Z_mass>111.19) && ht>=250 && mt2>200. && lep_pt0>100. && lep_pt1>30. && (nJets>1 || (mt2>250.)) ", "goff");
//  
////  integral_sf = drllSF->IntegralAndError(1,-1, error_sf);
////  integral_of = drllOF->IntegralAndError(1,-1, error_of);
////  
////  cout << "R^{SF/OF} = ("<<integral_sf<<" +- "<<error_sf<<")/("<<integral_of<<" +- "<<error_of<<")"<<endl;
////  RSFOF = integral_sf/integral_of;
////  err_R = TMath::Sqrt( (error_sf/integral_of)*(error_sf/integral_of) + (error_of*integral_sf/(integral_of*integral_of))*(error_of*integral_sf/(integral_of*integral_of)) );
////  cout << "R^{SF/OF} = "<< RSFOF << "+-" << err_R << endl;
//
//  TCanvas* c6=new TCanvas("c6", "", 600, 600);
//  c6->cd();
//  drllSF->Draw("PE");
//  drllOF->Draw("PE,same");
//  
//  TH1D* drllRSFOF = (TH1D*) drllSF->Clone("drllRSFOF");
//  drllRSFOF->Divide(drllOF);
//
//  drllRSFOF->SetMinimum(0.);
//  drllRSFOF->SetMaximum(2.);
//  
//  drllRSFOF->GetYaxis()->SetRangeUser(0.,2.);
//  drllRSFOF->GetYaxis()->SetTitle("R^{SF/OF}");
//  drllRSFOF->GetYaxis()->SetTitleOffset(1.4);
//
//  TCanvas* c6r=new TCanvas("c6r", "", 600, 600);
//  c6r->cd();
//
//  drllRSFOF->Draw("PE");
//  //  drllRSFOF->Fit("pol0");
//
//  TLine* lrdrll   = new TLine(0, (float) RSFOF,        5, (float)RSFOF);
//  TLine* lrdrllup = new TLine(0, (float)(RSFOF+err_R), 5, (float)(RSFOF+err_R));
//  TLine* lrdrlldn = new TLine(0, (float)(RSFOF-err_R), 5, (float)(RSFOF-err_R));
//
//  TLine* lrdrllup15 = new TLine(0, (float)(RSFOF+0.15), 5, (float)(RSFOF+0.15));
//  TLine* lrdrlldn15 = new TLine(0, (float)(RSFOF-0.15), 5, (float)(RSFOF-0.15));
//
//  lrdrll->SetLineColor(1);
//  lrdrllup->SetLineColor(1);
//  lrdrlldn->SetLineColor(1);
//  lrdrllup->SetLineStyle(2);
//  lrdrlldn->SetLineStyle(2);
//  
//  lrdrllup15->SetLineColor(2);
//  lrdrlldn15->SetLineColor(2);
//  lrdrllup15->SetLineStyle(2);
//  lrdrlldn15->SetLineStyle(2);
//  
//  lrdrll->Draw("same");
//  lrdrllup->Draw("same");
//  lrdrlldn->Draw("same");
//
//  lrdrllup15->Draw("same");
//  lrdrlldn15->Draw("same");
//
//  leg->Draw("same");
//  label_top->Draw("same");
//  label_cms->Draw("same");
//  gPad->RedrawAxis();
//
//  // MET (untransformed)
//  TH1D* rawmetSF = new TH1D("rawmetSF", "", 20,0,1000);
//  rawmetSF->Sumw2();
//
//  TH1D* rawmetOF = new TH1D("rawmetOF", "", 20,0,1000);
//  rawmetOF->Sumw2();
//
//  rawmetSF->SetLineColor(1);
//  rawmetSF->SetMarkerColor(1);
//  rawmetSF->SetMarkerStyle(20);
//  rawmetOF->SetLineColor(2);
//  rawmetOF->SetMarkerColor(2);
//  rawmetOF->SetMarkerStyle(20);
//
//  rawmetSF->GetXaxis()->SetTitle("E_{T}^{miss} (un-transformed) [GeV]");
//  rawmetOF->GetXaxis()->SetTitle("E_{T}^{miss} (un-transformed) [GeV]");
//  rawmetSF->GetYaxis()->SetTitle("Entries");
//  rawmetOF->GetYaxis()->SetTitle("Entries");
//
//  SF->Draw("raw_met>>rawmetSF", "Z_pt<=200 && Z_mass>=50 && (Z_mass<=71.19 || Z_mass>111.19) && ht>=250 && mt2>200. && lep_pt0>100. && lep_pt1>30. && (nJets>1 || (mt2>250.)) ", "goff");
//  OF->Draw("raw_met>>rawmetOF", "Z_pt<=200 && Z_mass>=50 && (Z_mass<=71.19 || Z_mass>111.19) && ht>=250 && mt2>200. && lep_pt0>100. && lep_pt1>30. && (nJets>1 || (mt2>250.)) ", "goff");
//  
////  integral_sf = rawmetSF->IntegralAndError(1,-1, error_sf);
////  integral_of = rawmetOF->IntegralAndError(1,-1, error_of);
////  
////  cout << "R^{SF/OF} = ("<<integral_sf<<" +- "<<error_sf<<")/("<<integral_of<<" +- "<<error_of<<")"<<endl;
////  RSFOF = integral_sf/integral_of;
////  err_R = TMath::Sqrt( (error_sf/integral_of)*(error_sf/integral_of) + (error_of*integral_sf/(integral_of*integral_of))*(error_of*integral_sf/(integral_of*integral_of)) );
////  cout << "R^{SF/OF} = "<< RSFOF << "+-" << err_R << endl;
//
//  TCanvas* c7=new TCanvas("c7", "", 600, 600);
//  c7->cd();
//  rawmetSF->Draw("PE");
//  rawmetOF->Draw("PE,same");
//  
//  TH1D* rawmetRSFOF = (TH1D*) rawmetSF->Clone("rawmetRSFOF");
//  rawmetRSFOF->Divide(rawmetOF);
//
//  rawmetRSFOF->SetMinimum(0.);
//  rawmetRSFOF->SetMaximum(2.);
//  
//  rawmetRSFOF->GetYaxis()->SetRangeUser(0.,2.);
//  rawmetRSFOF->GetYaxis()->SetTitle("R^{SF/OF}");
//  rawmetRSFOF->GetYaxis()->SetTitleOffset(1.4);
//
//  TCanvas* c7r=new TCanvas("c7r", "", 600, 600);
//  c7r->cd();
//
//  rawmetRSFOF->Draw("PE");
//  //  rawmetRSFOF->Fit("pol0");
//
//  TLine* lrrawmet   = new TLine(0, (float) RSFOF,        1000, (float)RSFOF);
//  TLine* lrrawmetup = new TLine(0, (float)(RSFOF+err_R), 1000, (float)(RSFOF+err_R));
//  TLine* lrrawmetdn = new TLine(0, (float)(RSFOF-err_R), 1000, (float)(RSFOF-err_R));
//
//  TLine* lrrawmetup15 = new TLine(0, (float)(RSFOF+0.15), 1000, (float)(RSFOF+0.15));
//  TLine* lrrawmetdn15 = new TLine(0, (float)(RSFOF-0.15), 1000, (float)(RSFOF-0.15));
//
//  lrrawmet->SetLineColor(1);
//  lrrawmetup->SetLineColor(1);
//  lrrawmetdn->SetLineColor(1);
//  lrrawmetup->SetLineStyle(2);
//  lrrawmetdn->SetLineStyle(2);
//  
//  lrrawmetup15->SetLineColor(2);
//  lrrawmetdn15->SetLineColor(2);
//  lrrawmetup15->SetLineStyle(2);
//  lrrawmetdn15->SetLineStyle(2);
//  
//  lrrawmet->Draw("same");
//  lrrawmetup->Draw("same");
//  lrrawmetdn->Draw("same");
//
//  lrrawmetup15->Draw("same");
//  lrrawmetdn15->Draw("same");
//
//  leg->Draw("same");
//  label_top->Draw("same");
//  label_cms->Draw("same");
//  gPad->RedrawAxis();
//
//  // Leading lepton pT
//  float bins_lep1pt[]={0,25,50,75,100,125,150,175,200,300};
//  TH1D* lep1ptSF = new TH1D("lep1ptSF", "", 9, bins_lep1pt);
//  lep1ptSF->Sumw2();
//
//  TH1D* lep1ptOF = new TH1D("lep1ptOF", "", 9, bins_lep1pt);
//  lep1ptOF->Sumw2();
//
//  lep1ptSF->SetLineColor(1);
//  lep1ptSF->SetMarkerColor(1);
//  lep1ptSF->SetMarkerStyle(20);
//  lep1ptOF->SetLineColor(2);
//  lep1ptOF->SetMarkerColor(2);
//  lep1ptOF->SetMarkerStyle(20);
//
//  lep1ptSF->GetXaxis()->SetTitle("p_{T}(leading lepton) [GeV]");
//  lep1ptOF->GetXaxis()->SetTitle("p_{T}(leading lepton) [GeV]");
//  lep1ptSF->GetYaxis()->SetTitle("Entries");
//  lep1ptOF->GetYaxis()->SetTitle("Entries");
//
//  SF->Draw("lep_pt0>>lep1ptSF", "Z_pt<=200 && Z_mass>=50 && (Z_mass<=71.19 || Z_mass>111.19) && ht>=250 && mt2>200. && lep_pt1>30. && (nJets>1 || (mt2>250.))", "goff");
//  OF->Draw("lep_pt0>>lep1ptOF", "Z_pt<=200 && Z_mass>=50 && (Z_mass<=71.19 || Z_mass>111.19) && ht>=250 && mt2>200. && lep_pt1>30. && (nJets>1 || (mt2>250.))", "goff");
//  
////  integral_sf = lep1ptSF->IntegralAndError(1,-1, error_sf);
////  integral_of = lep1ptOF->IntegralAndError(1,-1, error_of);
////  
////  cout << "R^{SF/OF} = ("<<integral_sf<<" +- "<<error_sf<<")/("<<integral_of<<" +- "<<error_of<<")"<<endl;
////  RSFOF = integral_sf/integral_of;
////  err_R = TMath::Sqrt( (error_sf/integral_of)*(error_sf/integral_of) + (error_of*integral_sf/(integral_of*integral_of))*(error_of*integral_sf/(integral_of*integral_of)) );
////  cout << "R^{SF/OF} = "<< RSFOF << "+-" << err_R << endl;
//
//  TCanvas* c8=new TCanvas("c8", "", 600, 600);
//  c8->cd();
//  lep1ptSF->Draw("PE");
//  lep1ptOF->Draw("PE,same");
//  
//  TH1D* lep1ptRSFOF = (TH1D*) lep1ptSF->Clone("lep1ptRSFOF");
//  lep1ptRSFOF->Divide(lep1ptOF);
//
//  lep1ptRSFOF->SetMinimum(0.);
//  lep1ptRSFOF->SetMaximum(2.);
//  
//  lep1ptRSFOF->GetYaxis()->SetRangeUser(0.,2.);
//  lep1ptRSFOF->GetYaxis()->SetTitle("R^{SF/OF}");
//
//  TCanvas* c8r=new TCanvas("c8r", "", 600, 600);
//  c8r->cd();
//
//  lep1ptRSFOF->Draw("PE");
//  //  lep1ptRSFOF->Fit("pol0");
//
//  TLine* lrlep1pt   = new TLine(25, (float) RSFOF,        300, (float)RSFOF);
//  TLine* lrlep1ptup = new TLine(25, (float)(RSFOF+err_R), 300, (float)(RSFOF+err_R));
//  TLine* lrlep1ptdn = new TLine(25, (float)(RSFOF-err_R), 300, (float)(RSFOF-err_R));
//
//  TLine* lrlep1ptup15 = new TLine(25, (float)(RSFOF+0.15), 300, (float)(RSFOF+0.15));
//  TLine* lrlep1ptdn15 = new TLine(25, (float)(RSFOF-0.15), 300, (float)(RSFOF-0.15));
//
//  lrlep1pt->SetLineColor(1);
//  lrlep1ptup->SetLineColor(1);
//  lrlep1ptdn->SetLineColor(1);
//  lrlep1ptup->SetLineStyle(2);
//  lrlep1ptdn->SetLineStyle(2);
//  
//  lrlep1ptup15->SetLineColor(2);
//  lrlep1ptdn15->SetLineColor(2);
//  lrlep1ptup15->SetLineStyle(2);
//  lrlep1ptdn15->SetLineStyle(2);
//  
//  lrlep1pt->Draw("same");
//  lrlep1ptup->Draw("same");
//  lrlep1ptdn->Draw("same");
//
//  lrlep1ptup15->Draw("same");
//  lrlep1ptdn15->Draw("same");
//
//  leg->Draw("same");
//  label_top->Draw("same");
//  label_cms->Draw("same");
//  gPad->RedrawAxis();
//
//
//  // Sub-leading lepton pT
//  float bins_lep2pt[]={0,20,30,40,60,80,100,150,300};
//  TH1D* lep2ptSF = new TH1D("lep2ptSF", "", 8, bins_lep2pt);
//  lep2ptSF->Sumw2();
//
//  TH1D* lep2ptOF = new TH1D("lep2ptOF", "", 8, bins_lep2pt);
//  lep2ptOF->Sumw2();
//
//  lep2ptSF->SetLineColor(1);
//  lep2ptSF->SetMarkerColor(1);
//  lep2ptSF->SetMarkerStyle(20);
//  lep2ptOF->SetLineColor(2);
//  lep2ptOF->SetMarkerColor(2);
//  lep2ptOF->SetMarkerStyle(20);
//
//  lep2ptSF->GetXaxis()->SetTitle("p_{T}(sub-leading lepton) [GeV]");
//  lep2ptOF->GetXaxis()->SetTitle("p_{T}(sub-leading lepton) [GeV]");
//  lep2ptSF->GetYaxis()->SetTitle("Entries");
//  lep2ptOF->GetYaxis()->SetTitle("Entries");
//
//  SF->Draw("lep_pt1>>lep2ptSF", "Z_pt<=200 && Z_mass>=50 && (Z_mass<=71.19 || Z_mass>111.19) && ht>=250 && mt2>200. && lep_pt0>100.", "goff");
//  OF->Draw("lep_pt1>>lep2ptOF", "Z_pt<=200 && Z_mass>=50 && (Z_mass<=71.19 || Z_mass>111.19) && ht>=250 && mt2>200. && lep_pt0>100.", "goff");
//  
////  integral_sf = lep2ptSF->IntegralAndError(1,-1, error_sf);
////  integral_of = lep2ptOF->IntegralAndError(1,-1, error_of);
////  
////  cout << "R^{SF/OF} = ("<<integral_sf<<" +- "<<error_sf<<")/("<<integral_of<<" +- "<<error_of<<")"<<endl;
////  RSFOF = integral_sf/integral_of;
////  err_R = TMath::Sqrt( (error_sf/integral_of)*(error_sf/integral_of) + (error_of*integral_sf/(integral_of*integral_of))*(error_of*integral_sf/(integral_of*integral_of)) );
////  cout << "R^{SF/OF} = "<< RSFOF << "+-" << err_R << endl;
//
//  TCanvas* c9=new TCanvas("c9", "", 600, 600);
//  c9->cd();
//  lep2ptSF->Draw("PE");
//  lep2ptOF->Draw("PE,same");
//  
//  TH1D* lep2ptRSFOF = (TH1D*) lep2ptSF->Clone("lep2ptRSFOF");
//  lep2ptRSFOF->Divide(lep2ptOF);
//
//  lep2ptRSFOF->SetMinimum(0.);
//  lep2ptRSFOF->SetMaximum(2.);
//  
//  lep2ptRSFOF->GetYaxis()->SetRangeUser(0.,2.);
//  lep2ptRSFOF->GetYaxis()->SetTitle("R^{SF/OF}");
//
//  TCanvas* c9r=new TCanvas("c9r", "", 600, 600);
//  c9r->cd();
//
//  lep2ptRSFOF->Draw("PE");
//  //  lep2ptRSFOF->Fit("pol0");
//
//  TLine* lrlep2pt   = new TLine(20, (float) RSFOF,        300, (float)RSFOF);
//  TLine* lrlep2ptup = new TLine(20, (float)(RSFOF+err_R), 300, (float)(RSFOF+err_R));
//  TLine* lrlep2ptdn = new TLine(20, (float)(RSFOF-err_R), 300, (float)(RSFOF-err_R));
//
//  TLine* lrlep2ptup15 = new TLine(20, (float)(RSFOF+0.15), 300, (float)(RSFOF+0.15));
//  TLine* lrlep2ptdn15 = new TLine(20, (float)(RSFOF-0.15), 300, (float)(RSFOF-0.15));
//
//  lrlep2pt->SetLineColor(1);
//  lrlep2ptup->SetLineColor(1);
//  lrlep2ptdn->SetLineColor(1);
//  lrlep2ptup->SetLineStyle(2);
//  lrlep2ptdn->SetLineStyle(2);
//  
//  lrlep2ptup15->SetLineColor(2);
//  lrlep2ptdn15->SetLineColor(2);
//  lrlep2ptup15->SetLineStyle(2);
//  lrlep2ptdn15->SetLineStyle(2);
//  
//  lrlep2pt->Draw("same");
//  lrlep2ptup->Draw("same");
//  lrlep2ptdn->Draw("same");
//
//  lrlep2ptup15->Draw("same");
//  lrlep2ptdn15->Draw("same");
//
//  leg->Draw("same");
//  label_top->Draw("same");
//  label_cms->Draw("same");
//  gPad->RedrawAxis();
//
//  // Leptons pT
//  TH2D* lepptSF_2d = new TH2D("lepptSF_2d", "", 9, bins_lep1pt, 8, bins_lep2pt);
//  lepptSF_2d->Sumw2();
//
//  TH2D* lepptOF_2d = new TH2D("lepptOF_2d", "", 9, bins_lep1pt, 8, bins_lep2pt);
//  lepptOF_2d->Sumw2();
//
//  lepptSF_2d->SetLineColor(1);
//  lepptSF_2d->SetMarkerColor(1);
//  lepptSF_2d->SetMarkerStyle(20);
//  lepptOF_2d->SetLineColor(2);
//  lepptOF_2d->SetMarkerColor(2);
//  lepptOF_2d->SetMarkerStyle(20);
//
//  lepptSF_2d->GetXaxis()->SetTitle("p_{T}(leading lepton) [GeV]");
//  lepptSF_2d->GetYaxis()->SetTitle("p_{T}(sub-leading lepton) [GeV]");
//  lepptOF_2d->GetXaxis()->SetTitle("p_{T}(leading lepton) [GeV]");
//  lepptOF_2d->GetYaxis()->SetTitle("p_{T}(sub-leading lepton) [GeV]");
////  lepptSF_2d->GetZaxis()->SetTitle("Entries");
////  lepptOF_2d->GetZaxis()->SetTitle("Entries");
//
//  SF->Draw("lep_pt1:lep_pt0>>lepptSF_2d", "Z_pt<=200 && Z_mass>=50 && (Z_mass<=71.19 || Z_mass>111.19) && ht>=250 && mt2>200.", "goff");
//  OF->Draw("lep_pt1:lep_pt0>>lepptOF_2d", "Z_pt<=200 && Z_mass>=50 && (Z_mass<=71.19 || Z_mass>111.19) && ht>=250 && mt2>200.", "goff");
//  
//
//  TCanvas* c10=new TCanvas("c10", "", 600, 600);
//  c10->cd();
//  lepptSF_2d->Draw("colz, text");
//
//  TCanvas* c10_2=new TCanvas("c10_2", "", 600, 600);
//  c10_2->cd();
//  lepptOF_2d->Draw("colz, text");
//
//  gStyle->SetPaintTextFormat("2.1f");
//  
//  TH2D* lepptRSFOF_2d = (TH2D*) lepptSF_2d->Clone("lepptRSFOF_2d");
//  lepptRSFOF_2d->Divide(lepptOF_2d);
//
//  lepptRSFOF_2d->SetMarkerSize(0.6);
//  lepptRSFOF_2d->GetYaxis()->SetTitleOffset(1.3);
//  lepptRSFOF_2d->GetZaxis()->SetRangeUser(0.,2);
//  //  lepptRSFOF_2d->GetZaxis()->SetTitle("R^{SF/OF}");
//
//  TCanvas* c10r=new TCanvas("c10r", "", 900, 900);
//  c10r->cd();
//
//  lepptRSFOF_2d->Draw("colz,texte");
//  //  lepptRSFOF->Fit("pol0");
//
////  TLine* lrleppt   = new TLine(20, (float) RSFOF,        300, (float)RSFOF);
////  TLine* lrlepptup = new TLine(20, (float)(RSFOF+err_R), 300, (float)(RSFOF+err_R));
////  TLine* lrlepptdn = new TLine(20, (float)(RSFOF-err_R), 300, (float)(RSFOF-err_R));
////
////  TLine* lrlepptup15 = new TLine(20, (float)(RSFOF+0.15), 300, (float)(RSFOF+0.15));
////  TLine* lrlepptdn15 = new TLine(20, (float)(RSFOF-0.15), 300, (float)(RSFOF-0.15));
////
////  lrleppt->SetLineColor(1);
////  lrlepptup->SetLineColor(1);
////  lrlepptdn->SetLineColor(1);
////  lrlepptup->SetLineStyle(2);
////  lrlepptdn->SetLineStyle(2);
////  
////  lrlepptup15->SetLineColor(2);
////  lrlepptdn15->SetLineColor(2);
////  lrlepptup15->SetLineStyle(2);
////  lrlepptdn15->SetLineStyle(2);
////  
////  lrleppt->Draw("same");
////  lrlepptup->Draw("same");
////  lrlepptdn->Draw("same");
////
////  lrlepptup15->Draw("same");
////  lrlepptdn15->Draw("same");



  //SAVE as:
//  c1 ->SaveAs("RSFOF_noZwindow_2016re-reco_full_l1pt100_l2pt30/SFnOF_mll.pdf");
//  c1r->SaveAs("RSFOF_noZwindow_2016re-reco_full_l1pt100_l2pt30/RSFOF_mll.pdf");

  c2 ->SaveAs("RSFOF_noZwindow_2016re-reco_full_l1pt100_l2pt30/SFnOF_ht.pdf");
  c2r->SaveAs("RSFOF_noZwindow_2016re-reco_full_l1pt100_l2pt30/RSFOF_ht.pdf");

  c3 ->SaveAs("RSFOF_noZwindow_2016re-reco_full_l1pt100_l2pt30/SFnOF_nj.pdf");
  c3r->SaveAs("RSFOF_noZwindow_2016re-reco_full_l1pt100_l2pt30/RSFOF_nj.pdf");

  c4 ->SaveAs("RSFOF_noZwindow_2016re-reco_full_l1pt100_l2pt30/SFnOF_nb.pdf");
  c4r->SaveAs("RSFOF_noZwindow_2016re-reco_full_l1pt100_l2pt30/RSFOF_nb.pdf");

  c5 ->SaveAs("RSFOF_noZwindow_2016re-reco_full_l1pt100_l2pt30/SFnOF_mt2.pdf");
  c5r->SaveAs("RSFOF_noZwindow_2016re-reco_full_l1pt100_l2pt30/RSFOF_mt2.pdf");

//  c6 ->SaveAs("RSFOF_noZwindow_2016re-reco_full_l1pt100_l2pt30/SFnOF_drll.pdf");
//  c6r->SaveAs("RSFOF_noZwindow_2016re-reco_full_l1pt100_l2pt30/RSFOF_drll.pdf");
//
//  c7 ->SaveAs("RSFOF_noZwindow_2016re-reco_full_l1pt100_l2pt30/SFnOF_rawmet.pdf");
//  c7r->SaveAs("RSFOF_noZwindow_2016re-reco_full_l1pt100_l2pt30/RSFOF_rawmet.pdf");
//
//  c8 ->SaveAs("RSFOF_noZwindow_2016re-reco_full_l1pt100_l2pt30/SFnOF_lep1pt_lep2pt30.pdf");
//  c8r->SaveAs("RSFOF_noZwindow_2016re-reco_full_l1pt100_l2pt30/RSFOF_lep1pt_lep2pt30.pdf");
//
//  c9 ->SaveAs("RSFOF_noZwindow_2016re-reco_full_l1pt100_l2pt30/SFnOF_lep2pt_lep1pt100.pdf");
//  c9r->SaveAs("RSFOF_noZwindow_2016re-reco_full_l1pt100_l2pt30/RSFOF_lep2pt_lep1pt100.pdf");
//
//  c10 ->SaveAs("RSFOF_noZwindow_2016re-reco_full_l1pt100_l2pt30/SF_leppt2d.pdf");
//  c10_2 ->SaveAs("RSFOF_noZwindow_2016re-reco_full_l1pt100_l2pt30/OF_leppt2d.pdf");
//  c10r->SaveAs("RSFOF_noZwindow_2016re-reco_full_l1pt100_l2pt30/RSFOF_leppt2d.pdf");


}

//  LocalWords:  TFile
