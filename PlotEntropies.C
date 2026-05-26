#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TMath.h>

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

const int NENERGIES = 4;
const double lgE[NENERGIES] = { 17.5, 18.0, 18.5, 19.0 };
const double ScalcNentryFit[NENERGIES] = {7.21, 7.26, 7.42, 7.46};
const double dScalcNentryFit[NENERGIES] = {0.02, 0.02, 0.02, 0.03};

void PlotEntropies()
{
	TFile * infile_entropies = TFile::Open("combined_entropy_per_energy.root", "READ");
	
	// TGraphErrors * comb_gr_Scalc    = new TGraphErrors();
	// TGraphErrors * comb_gr_avgScalc = new TGraphErrors();
	// TGraphErrors * comb_gr_lnnavg   = new TGraphErrors();
	TGraphErrors * comb_gr_Scalc    = (TGraphErrors*)infile_entropies->Get("Scalc");
	TGraphErrors * comb_gr_avgScalc = (TGraphErrors*)infile_entropies->Get("avgScalc");
	TGraphErrors * comb_gr_lnnavg   = (TGraphErrors*)infile_entropies->Get("lnnavg");
	TGraphAsymmErrors * syst_gr_Scalc = (TGraphAsymmErrors*)infile_entropies->Get("syst_gr_Scalc");
	TGraphErrors * S_lnnavg				    = new TGraphErrors(NENERGIES);
	TGraphAsymmErrors * syst_S_lnnavg = new TGraphAsymmErrors(NENERGIES);
	TGraphErrors * S_lnnavg_infty	    = new TGraphErrors(NENERGIES);
	TGraphErrors * delta_S 				    = new TGraphErrors(NENERGIES);

	
	// comb_gr_Scalc = 
	// comb_gr_avgScalc = 
	// comb_gr_lnnavg = 
	
	for( int ien = 0 ; ien < NENERGIES; ien++ )
	{
		double _lgE = 0.;
		double _Scalc = 0.;
		double _lnnavg = 0.;
		double _avgScalc = 0.;
		double _systScalc = 0.;
		
    comb_gr_Scalc->GetPoint(ien,_lgE,_Scalc);
    comb_gr_lnnavg->GetPoint(ien,_lgE,_lnnavg);
    comb_gr_avgScalc->GetPoint(ien,_lgE,_avgScalc);
    syst_gr_Scalc->GetPoint(ien,_lgE,_systScalc);
		double d_lnnavg = comb_gr_lnnavg->GetErrorY(ien);
		double d_Scalc = comb_gr_Scalc->GetErrorY(ien);
		
		double syst_S_dn = syst_gr_Scalc->GetErrorYlow(ien);
		double syst_S_up = syst_gr_Scalc->GetErrorYhigh(ien);
		
		delta_S->SetPoint(ien,_lgE,(_Scalc-_avgScalc)/_Scalc*100.);
		S_lnnavg->SetPoint(ien,_lnnavg,_Scalc);
		S_lnnavg->SetPointError(ien,d_lnnavg,d_Scalc);
		
		S_lnnavg_infty->SetPoint(ien,_lnnavg,ScalcNentryFit[ien]);
		S_lnnavg_infty->SetPointError(ien,d_lnnavg,dScalcNentryFit[ien]);
		
		syst_S_lnnavg->SetPoint(ien,_lnnavg,_Scalc);
		syst_S_lnnavg->SetPointError(ien,0.03,0.03,syst_S_dn,syst_S_up);
		
	}
	
	double grmin = 5.;
	double grmax = 8.49;

	TCanvas* c1 = new TCanvas("c1", "c1", 900, 700);
	c1->SetMargin(0.12, 0.04, 0.12, 0.06);
	c1->SetLogy(0);
	c1->SetLogx(0);

	TLegend * leg = new TLegend(7.0,grmin,grmax,5.75,"","false");
	
	S_lnnavg_infty->SetTitle(";ln#LTn#GT; S_{h}");
	S_lnnavg_infty->SetMarkerStyle(24);
	S_lnnavg_infty->SetMarkerColor(kGray+1);
	S_lnnavg_infty->SetLineColor(kGray+1);
	S_lnnavg_infty->SetMarkerSize(1.50);
	S_lnnavg_infty->SetLineWidth(1);
	S_lnnavg_infty->GetXaxis()->SetLimits(grmin, grmax);
	S_lnnavg_infty->SetMinimum(grmin);
	S_lnnavg_infty->SetMaximum(grmax);
	S_lnnavg_infty->Draw("APE");
	
	leg->AddEntry(S_lnnavg_infty,"MC extrapolated");
	
	S_lnnavg->SetMarkerStyle(20);
	S_lnnavg->SetMarkerColor(2);
	S_lnnavg->SetLineColor(2);
	S_lnnavg->SetMarkerSize(1.0);
	S_lnnavg->SetLineWidth(1);
	S_lnnavg->GetXaxis()->SetLimits(grmin, grmax);
	S_lnnavg->SetMinimum(grmin);
	S_lnnavg->SetMaximum(grmax);
	S_lnnavg->Draw("PE SAME");
	
	S_lnnavg->Print();
	
	leg->AddEntry(S_lnnavg,"Auger based MC");
	
	syst_S_lnnavg->Print();
	syst_S_lnnavg->SetLineWidth(1);
	syst_S_lnnavg->SetLineColor(2);
	syst_S_lnnavg->SetFillColorAlpha(2, 0.3);
	syst_S_lnnavg->SetLineStyle(2);
	syst_S_lnnavg->Draw("E2 SAME");
	
	syst_S_lnnavg->Print();
	
	TF1* universal2 = new TF1("ln_navg","x+1", grmin, grmax);
	universal2->SetLineColor(kBlack); // Customize appearance
	universal2->SetLineStyle(2); // Customize appearance
	universal2->SetLineWidth(3);
	universal2->Draw("SAME");
	leg->AddEntry(universal2,"ln#LTn#GT+1","l");
	
	leg->Draw("SAME");
	
		
	c1->Update();
	c1->Print("figs/S_lnnavg.png");
	
	
	
	TCanvas* c2 = new TCanvas("c2", "c2", 900, 700);
	c2->SetMargin(0.12, 0.04, 0.12, 0.06);
	c2->SetLogy(0);
	c2->SetLogx(0);

	TLegend * leg2 = new TLegend(18.8,6.0,19.5,7.,"","false");
	
	delta_S->SetTitle(";ln#LTn#GT; #delta S [%]");
	delta_S->SetMarkerStyle(20);
	delta_S->SetMarkerSize(1.0);
	delta_S->SetLineWidth(1);
	delta_S->GetXaxis()->SetLimits(17., 19.5);
	delta_S->SetMinimum(0.);
	delta_S->SetMaximum(7.);
	delta_S->Draw("APE");
	
	leg2->AddEntry(delta_S,"Auger based MC");
		
	leg2->Draw("SAME");
	
	c2->Update();
	c2->Print("figs/delta_S.png");
	
	
}
