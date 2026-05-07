#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TMath.h>
#include "TRandom3.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

const int NFILES = 4;
const int NENERGIES = 4;
const int NTOY = 1000;
const double lgE[NENERGIES] = { 17.5, 18.0, 18.5, 19.0 };
const char * name_primaries[NFILES] = { "H", "He", "N", "Fe" };

const double eposlhc_fractions[NENERGIES][NFILES] = {
																										 { 0.232517 , 0.129347 , 0.480098 , 0.158038 },
																										 { 0.398609 , 0.032087 , 0.414095 , 0.155210 },
																										 { 0.258519 , 0.212748 , 0.468386 , 0.060348 },
																										 { 0.140560 , 0.000140 , 0.843293 , 0.016006 }
																										};

const double eposlhc_fractions_uncert[NENERGIES][NFILES][2] = { //   H donw			H up			 He down   He up			 N down	  	N up			 Fe down		Fe up
																																{ { 0.086101, 0.073556}, {0.080740, 0.087469}, {0.182908, 0.065389}, {0.116334, 0.195281} },
																																{ { 0.124844, 0.069415}, {0.039534, 0.152586}, {0.167176, 0.042384}, {0.120518, 0.184949} },
																																{ { 0.139268, 0.076605}, {0.087848, 0.119351}, {0.089119, 0.058500}, {0.059245, 0.143624} },
																																{ { 0.064260, 0.054008}, {0.090952, 0.299223}, {0.282756, 0.067671}, {0.000199, 0.112543} }
																															};






// =====================================================
// Negative binomial PMF with parameters mean and k
// =====================================================
double NBDpmf(int n, double mean, double k)
{
  if (n < 0 || mean <= 0.0 || k <= 0.0) {
    return 0.0;
  }

  const double logP =
      TMath::LnGamma(n + k)
    - TMath::LnGamma(k)
    - TMath::LnGamma(n + 1.0)
    + n * log(mean / (mean + k))
    + k * log(k / (mean + k));

  return exp(logP);
}

double TripleNBD(double* x, double* p)
{
  // p[0] = norm
  // p[1] = alpha1
  // p[2] = alpha2
  // p[3], p[4] = mean1, k1
  // p[5], p[6] = mean2, k2
  // p[7], p[8] = mean3, k3

  double n = x[0];
  double norm   = p[0];
  double alpha1 = p[1];
  double alpha2 = p[2];
  double alpha3 = 1.0 - alpha1 - alpha2;

  if (alpha1 < 0.0 || alpha1 > 1.0) return 1e30;
  if (alpha2 < 0.0 || alpha2 > 1.0) return 1e30;
  if (alpha3 < 0.0 || alpha3 > 1.0) return 1e30;

  double f1 = NBDpmf((int)TMath::Nint(n), p[3], p[4]);
  double f2 = NBDpmf((int)TMath::Nint(n), p[5], p[6]);
  double f3 = NBDpmf((int)TMath::Nint(n), p[7], p[8]);

  return norm * (alpha1 * f1 + alpha2 * f2 + alpha3 * f3);
}

double QuadrupleNBD(double* x, double* p)
{
  // p[0]  = norm
  // p[1]  = alpha1
  // p[2]  = alpha2
  // p[3]  = alpha3
  // p[4]  = mean1
  // p[5]  = k1
  // p[6]  = mean2
  // p[7]  = k2
  // p[8]  = mean3
  // p[9]  = k3
  // p[10] = mean4
  // p[11] = k4

  const double xx = x[0];
  const int n = static_cast<int>(TMath::Nint(xx));

  if (fabs(xx - n) > 0.5001) {
    return 0.0;
  }

  const double norm   = p[0];
  const double alpha1 = p[1];
  const double alpha2 = p[2];
  const double alpha3 = p[3];
  const double alpha4 = 1.0 - alpha1 - alpha2 - alpha3;

  // Enforce physical weights
  if (alpha1 < 0.0 || alpha1 > 1.0) return 1e30;
  if (alpha2 < 0.0 || alpha2 > 1.0) return 1e30;
  if (alpha3 < 0.0 || alpha3 > 1.0) return 1e30;
  if (alpha4 < 0.0 || alpha4 > 1.0) return 1e30;

  const double mean1 = p[4];
  const double k1    = p[5];
  const double mean2 = p[6];
  const double k2    = p[7];
  const double mean3 = p[8];
  const double k3    = p[9];
  const double mean4 = p[10];
  const double k4    = p[11];

  if (mean1 <= 0.0 || k1 <= 0.0) return 1e30;
  if (mean2 <= 0.0 || k2 <= 0.0) return 1e30;
  if (mean3 <= 0.0 || k3 <= 0.0) return 1e30;
  if (mean4 <= 0.0 || k4 <= 0.0) return 1e30;

  const double f1 = NBDpmf(n, mean1, k1);
  const double f2 = NBDpmf(n, mean2, k2);
  const double f3 = NBDpmf(n, mean3, k3);
  const double f4 = NBDpmf(n, mean4, k4);

  return norm * (alpha1 * f1 + alpha2 * f2 + alpha3 * f3 + alpha4 * f4);
}

double ComputeMeanFromHistogram(const TH1* h)
{
  if (!h) {
    return -1.0;
  }

  const int nBins = h->GetNbinsX();

  double sumCounts = 0.0;
  double sumWeighted = 0.0;

  for (int iBin = 1; iBin <= nBins; ++iBin) {
    const double n = h->GetBinCenter(iBin);
    const double c = h->GetBinContent(iBin);

    if (c > 0.0) {
      sumCounts   += c;
      sumWeighted += n * c;
    }
  }

  if (sumCounts <= 0.0) {
    return -1.0;
  }

  return sumWeighted / sumCounts;
}

double ComputeMeanUncertaintyFromHistogram(const TH1* h)
{
  if (!h) {
    return -1.0;
  }

  const int nBins = h->GetNbinsX();

  double sumCounts = 0.0;
  for (int iBin = 1; iBin <= nBins; ++iBin) {
    const double c = h->GetBinContent(iBin);
    if (c > 0.0) {
      sumCounts += c;
    }
  }

  if (sumCounts <= 0.0) {
    return -1.0;
  }

  double error2 = 0.0;

  for (int iBin = 1; iBin <= nBins; ++iBin) {
    const double n  = h->GetBinCenter(iBin);
    const double c  = h->GetBinContent(iBin);
    const double dc = h->GetBinError(iBin);

    if (c <= 0.0 || dc <= 0.0) {
      continue;
    }

    const double dp = dc / sumCounts;

    error2 += (n * dp) * (n * dp);
  }

  return sqrt(error2);
}

double ComputeLogMeanUncertaintyFromHistogram(const TH1* h)
{
  const double mean = ComputeMeanFromHistogram(h);
  const double dmean = ComputeMeanUncertaintyFromHistogram(h);

  if (mean <= 0.0 || dmean < 0.0) {
    return -1.0;
  }

  return dmean / mean;
}

// =====================================================
// Entropy from fitted NBD
// S = -sum P(n) ln P(n)
// =====================================================
double ComputeEntropyFromNBD(double mean, double k, int nMax)
{
  if (mean <= 0.0 || k <= 0.0 || nMax < 0) {
    return -1.0;
  }

  vector<double> probs;
  probs.reserve(nMax + 1);

  double sumP = 0.0;
  for (int n = 0; n <= nMax; ++n) {
    const double p = NBDpmf(n, mean, k);
    probs.push_back(p);
    sumP += p;
  }

  if (sumP <= 0.0) {
    return -1.0;
  }

  double S = 0.0;
  for (double p : probs) {
    const double pn = p / sumP;
    if (pn > 0.0) {
      S -= pn * log(pn);
    }
  }

  return S;
}

double ComputeEntropyFromQuadNBD(double alpha1, double alpha2, double alpha3,
                                 double mean1,  double k1,
                                 double mean2,  double k2,
                                 double mean3,  double k3,
                                 double mean4,  double k4,
                                 int nMax)
{
  const double alpha4 = 1.0 - alpha1 - alpha2 - alpha3;

  if (nMax < 0) {
    return -1.0;
  }

  if (alpha1 < 0.0 || alpha2 < 0.0 || alpha3 < 0.0 || alpha4 < 0.0) {
    return -1.0;
  }

  if (mean1 <= 0.0 || k1 <= 0.0 ||
      mean2 <= 0.0 || k2 <= 0.0 ||
      mean3 <= 0.0 || k3 <= 0.0 ||
      mean4 <= 0.0 || k4 <= 0.0) {
    return -1.0;
  }

  vector<double> probs;
  probs.reserve(nMax + 1);

  double sumP = 0.0;

  for (int n = 0; n <= nMax; ++n) {
    const double p1 = NBDpmf(n, mean1, k1);
    const double p2 = NBDpmf(n, mean2, k2);
    const double p3 = NBDpmf(n, mean3, k3);
    const double p4 = NBDpmf(n, mean4, k4);

    const double p =
        alpha1 * p1 +
        alpha2 * p2 +
        alpha3 * p3 +
        alpha4 * p4;

    probs.push_back(p);
    sumP += p;
  }

  if (sumP <= 0.0) {
    return -1.0;
  }

  double S = 0.0;
  for (double p : probs) {
    const double pn = p / sumP;
    if (pn > 0.0) {
      S -= pn * log(pn);
    }
  }

  return S;
}

double ComputeEntropyFromHistogram(const TH1* h)
{
  if (!h) {
    return -1.0;
  }

  const int nBins = h->GetNbinsX();
  double sumCounts = 0.0;

  for (int iBin = 1; iBin <= nBins; ++iBin) {
    const double c = h->GetBinContent(iBin);
    if (c > 0.0) {
      sumCounts += c;
    }
  }

  if (sumCounts <= 0.0) {
    return -1.0;
  }

  double S = 0.0;
  for (int iBin = 1; iBin <= nBins; ++iBin) {
    const double c = h->GetBinContent(iBin);
    if (c > 0.0) {
      const double p = c / sumCounts;
      S -= p * log(p);
    }
  }

  return S;
}

double ComputeEntropyUncertaintyFromHistogram(const TH1* h)
{
  if (!h) {
    return -1.0;
  }

  const int nBins = h->GetNbinsX();

  double sumCounts = 0.0;
  for (int iBin = 1; iBin <= nBins; ++iBin) {
    const double c = h->GetBinContent(iBin);
    if (c > 0.0) {
      sumCounts += c;
    }
  }

  if (sumCounts <= 0.0) {
    return -1.0;
  }

  double entropy_error2 = 0.0;

  for (int iBin = 1; iBin <= nBins; ++iBin) {
    const double c  = h->GetBinContent(iBin);
    const double dc = h->GetBinError(iBin);

    if (c <= 0.0 || dc <= 0.0) {
      continue;
    }

    const double p  = c / sumCounts;
    const double dp = dc / sumCounts;

    if (p > 0.0 && dp > 0.0) {
      const double dSdp = -log(p) - 1.0;
      entropy_error2 += (dSdp * dp) * (dSdp * dp);
    }
  }

  return sqrt(entropy_error2);
}

void CombinedEntropy()
{
	TRandom3 rng(0);
	
	TFile * infile_mults = TFile::Open("multiplicities_per_primaries.root", "READ");
	TH1D  * hMult[NFILES][NENERGIES];
	
	TFile * infile_entropies = TFile::Open("entropy_per_primaries.root", "READ");
	TGraphErrors * gr_Sfit[NFILES];
	TGraphErrors * gr_Scalc[NFILES];
	TGraphErrors * comb_gr_avgScalc = new TGraphErrors(NENERGIES);
	TGraphErrors * comb_gr_avgSfit  = new TGraphErrors(NENERGIES);
	TGraphErrors * comb_gr_Sfit     = new TGraphErrors(NENERGIES);
	TGraphErrors * comb_gr_Scalc    = new TGraphErrors(NENERGIES);
	TGraphErrors * comb_gr_lnnavg   = new TGraphErrors(NENERGIES);
	TGraphAsymmErrors * syst_gr_Scalc = new TGraphAsymmErrors(NENERGIES);
	
	for( int ifile = 0 ; ifile < NFILES ; ifile++)
	{
		cerr << "Processing primary " << name_primaries[ifile] << endl;
		
		gr_Sfit[ifile] = (TGraphErrors*)infile_entropies->Get(Form("Sfit_primary_%s",name_primaries[ifile]));
		gr_Scalc[ifile] = (TGraphErrors*)infile_entropies->Get(Form("Scalc_primary_%s",name_primaries[ifile]));
	}
	
	for ( int ien = 0 ; ien < NENERGIES ; ien++)
	{
		double avgScalc = 0.;
		double avgSfit = 0.;
			
		for( int ifile = 0 ; ifile < NFILES ; ifile++)
		{
			double lgE = 0.;
			double Sfit = 0.;
			double Scalc = 0.;
			gr_Sfit[ifile]->GetPoint(ifile,lgE,Sfit);
			gr_Scalc[ifile]->GetPoint(ifile,lgE,Scalc);
			
			avgScalc += eposlhc_fractions[ien][ifile] * Scalc;
			avgSfit += eposlhc_fractions[ien][ifile] * Sfit;
		}
		
		comb_gr_avgScalc->SetName("avgScalc");
		comb_gr_avgScalc->SetPoint(ien,lgE[ien],avgScalc);
		comb_gr_avgSfit->SetName("avgSfit");
		comb_gr_avgSfit->SetPoint(ien,lgE[ien],avgSfit);
	}	
	
	TH1D* hMix[NENERGIES];
	
	for( int ien = 0 ; ien < NENERGIES ; ien++)
	{
		hMult[0][ien] = (TH1D*)infile_mults->Get(Form("hMult_en%i_prim%i", ien, 0));
		hMix[ien] = (TH1D*)hMult[0][ien]->Clone(Form("hPnMix_lgE%.1f", lgE[ien]));
		hMix[ien]->Reset();
		hMix[ien]->SetDirectory(nullptr);
		
		for( int ifile = 0 ; ifile < NFILES ; ifile++)
		{
			hMult[ifile][ien] = (TH1D*)infile_mults->Get(Form("hMult_en%i_prim%i",ien,ifile));
			hMix[ien]->Add(hMult[ifile][ien], eposlhc_fractions[ien][ifile]);
		}
		
		vector<double> entropyValues_dn;
		vector<double> entropyValues_up;
		
		for (int itoy = 0; itoy < NTOY; ++itoy)
		{
			double gup[NFILES];
			double gdn[NFILES];
			bool bad = false;

			for (int ifile = 0; ifile < NFILES; ++ifile)
			{
				gdn[ifile] = eposlhc_fractions[ien][ifile] - rng.Gaus(0.0, eposlhc_fractions_uncert[ien][ifile][0]);
				gup[ifile] = eposlhc_fractions[ien][ifile] + rng.Gaus(0.0, eposlhc_fractions_uncert[ien][ifile][1]);
				if (gdn[ifile] < 0.0 || gup[ifile] < 0.0) bad = true;
			}

			if (bad)
			{
				--itoy;
				continue;
			}

			double sumgdn = gdn[0] + gdn[1] + gdn[2] + gdn[3];
			double sumgup = gup[0] + gup[1] + gup[2] + gup[3];

			double fdn[NFILES];
			double fup[NFILES];
			for (int ifile = 0; ifile < NFILES; ++ifile) 
			{
				fdn[ifile] = gdn[ifile] / sumgdn;
				fup[ifile] = gup[ifile] / sumgup;
			}

			TH1D* hMix_uncert_dn = (TH1D*)hMix[ien]->Clone("hMix_uncert_dn");
			TH1D* hMix_uncert_up = (TH1D*)hMix[ien]->Clone("hMix_uncert_up");
			hMix_uncert_dn->Reset();
			hMix_uncert_up->Reset();

			for (int ifile = 0; ifile < NFILES; ++ifile)
			{
				hMix_uncert_dn->Add(hMult[ifile][ien], fdn[ifile]);
				hMix_uncert_up->Add(hMult[ifile][ien], fup[ifile]);
			}
			
			entropyValues_dn.push_back(ComputeEntropyFromHistogram(hMix_uncert_dn));
			entropyValues_up.push_back(ComputeEntropyFromHistogram(hMix_uncert_up));

			delete hMix_uncert_dn;
			delete hMix_uncert_up;
		}
		
		for(int i = 0 ; i < entropyValues_dn.size() ; i++)
			if(i%1000==0) cerr << i << " " << entropyValues_dn.at(i) << " " << entropyValues_up.at(i) << endl;
		
		cerr << "Mix integral at lgE = " << lgE[ien] << " : " << hMix[ien]->Integral() << endl;
			
	
		// -----------------------------------------
		// Fit range
		// -----------------------------------------
		/*
		const int firstBin = hMix[ien]->FindFirstBinAbove(0.0);
		const int lastBin  = hMix[ien]->FindLastBinAbove(0.0);

		if (firstBin < 1 || lastBin < firstBin) {
			cerr << "ERROR: no non-empty bins found." << endl;
			return;
		}

		const double xMin = hMix[ien]->GetBinCenter(firstBin);
		const double xMax = hMix[ien]->GetBinCenter(lastBin);

		// -----------------------------------------
		// Initial parameters
		// -----------------------------------------
		const double normInit = hMix[ien]->Integral(firstBin, lastBin);
		const double meanInit = hMix[ien]->GetMean();
		const double rmsInit  = hMix[ien]->GetRMS();
		const double varInit  = rmsInit * rmsInit;

		double kInit = 2.0;
		if (varInit > meanInit && meanInit > 0.0) {
			kInit = meanInit * meanInit / (varInit - meanInit);
			if (!std::isfinite(kInit) || kInit <= 0.0) {
				kInit = 2.0;
			}
		}
		
		TF1* fNBD = new TF1("fNBD", QuadrupleNBD, 1, xMax, 12);
		// fNBD->SetParNames("Norm", "alpha_{1}", "alpha_{2}", "<n>_{1}", "k_{1}", "<n>_{2}", "k_{2}", "<n>_{3}", "k_{3}");
		// TF1* f4NBD = new TF1("f4NBD", QuadrupleNBD, xMin, xMax, 12);
		fNBD->SetParName(0,  "Norm");
		fNBD->SetParName(1,  "alpha1");
		fNBD->SetParName(2,  "alpha2");
		fNBD->SetParName(3,  "alpha3");
		fNBD->SetParName(4,  "mean1");
		fNBD->SetParName(5,  "k1");
		fNBD->SetParName(6,  "mean2");
		fNBD->SetParName(7,  "k2");
		fNBD->SetParName(8,  "mean3");
		fNBD->SetParName(9,  "k3");
		fNBD->SetParName(10, "mean4");
		fNBD->SetParName(11, "k4");

		double pars[12] = {
			1, 0.4, 0.3, 0.1,
			1, 2,
			60, 40,
			160, 10,
			1000, 2
		};
		fNBD->SetParameters(pars);

		fNBD->SetParLimits(0, 1.0, 1e6);   // norm
		fNBD->SetParLimits(1, 0.0, 1.0);   // alpha1
		fNBD->SetParLimits(2, 0.0, 1.0);   // alpha2
		fNBD->SetParLimits(3, 0.0, 1.0);   // alpha3
		fNBD->SetParLimits(4, 0.1, 5e1);   // mean1
		fNBD->SetParLimits(5, 0.1, 1e2);   // k1
		fNBD->SetParLimits(6, 0.1, 1e3);   // mean2
		fNBD->SetParLimits(7, 0.1, 1e3);   // k2
		fNBD->SetParLimits(8, 1e2, 2e2);   // mean3
		fNBD->SetParLimits(9, 0.1, 1e2);   // k3
		fNBD->SetParLimits(10, 2e2, 1e4);  // mean4
		fNBD->SetParLimits(11, 0.1, 1e2);  // k4

		hMix[ien]->Fit(fNBD, "R0");

		const double alpha1 = fNBD->GetParameter(1);
		const double alpha2 = fNBD->GetParameter(2);
		const double alpha3 = fNBD->GetParameter(3);

		const double mean1 = fNBD->GetParameter(4);
		const double k1    = fNBD->GetParameter(5);
		const double mean2 = fNBD->GetParameter(6);
		const double k2    = fNBD->GetParameter(7);
		const double mean3 = fNBD->GetParameter(8);
		const double k3    = fNBD->GetParameter(9);
		const double mean4 = fNBD->GetParameter(10);
		const double k4    = fNBD->GetParameter(11);
		
		const double chi2 = fNBD->GetChisquare();
		const int ndf     = fNBD->GetNDF();
		const double chi2ndf = (ndf > 0) ? chi2 / ndf : -1.0;
		const double conflev = TMath::Prob(chi2,ndf);
		const double Sfit = ComputeEntropyFromQuadNBD(alpha1, alpha2, alpha3, mean1, k1, mean2, k2, mean3, k3, mean4, k4, xMax);
		
		
		// -----------------------------------------
		// Draw
		// -----------------------------------------
		TCanvas* c = new TCanvas("cMultNBD", "LeadingInteractions mult NBD fit", 900, 700);
		c->SetMargin(0.12, 0.04, 0.12, 0.06);
		c->SetLogy(1);
		c->SetLogx(1);

		hMix[ien]->SetMarkerStyle(20);
		hMix[ien]->SetMarkerSize(0.9);
		hMix[ien]->SetLineWidth(2);
		hMix[ien]->GetXaxis()->SetRangeUser(0, 2000);
		hMix[ien]->Draw("HIST");

		fNBD->SetLineWidth(2);
		fNBD->Draw("SAME");

		TLatex lat;
		lat.SetNDC();
		lat.SetTextFont(42);
		lat.SetTextSize(0.035);

		const double xText = 0.60;
		const double yText = 0.86;
		const double dy    = 0.05;

		lat.DrawLatex(xText, yText,            Form("#chi^{2}/NDF = %.3f, CL=%.1f%%", chi2ndf, conflev*100));
		// lat.DrawLatex(xText, yText - 1.0 * dy, Form("k = %.4f", k));
		lat.DrawLatex(xText, yText - 1.0 * dy, Form("S_{fit} = %.4f", Sfit));
		lat.DrawLatex(xText, yText - 2.0 * dy, Form("S_{hist} = %.4f", Shist));
		lat.DrawLatex(xText, yText - 3.0 * dy, Form("ln<n_{hist}> = %.4f #pm %.4f", log(navghist), navgherr));

		c->Update();
		c->Print(Form("figs/combined_fitmult_en%i.png",ien));
		*/
		
		const double Shist = ComputeEntropyFromHistogram(hMix[ien]);
		const double Sherr = ComputeEntropyUncertaintyFromHistogram(hMix[ien]);
		const double navghist = ComputeMeanFromHistogram(hMix[ien]);
		const double navgherr = ComputeLogMeanUncertaintyFromHistogram(hMix[ien]);
		
		double syst_dn = Shist - *min_element(entropyValues_dn.begin(), entropyValues_dn.end());
		double syst_up = *max_element(entropyValues_up.begin(), entropyValues_up.end())-Shist;
		
		cout << "========================================" << endl;
		cout << "Entries   : " << hMix[ien]->GetEntries() << endl;
		// cout << "chi2/NDF  : " << chi2ndf << endl;
		// cout << "k         : " << k << endl;
		cout << "<n>       : " << navghist << endl;
		cout << "1+ln<n>   : " << 1.+log(navghist) << "#pm" << navgherr << endl;
		// cout << "Sfit      : " << Sfit << endl;
		cout << "Shist     : " << Shist << "#pm" << Sherr << "-" << syst_dn << "+" << syst_up << endl;
		cout << "========================================" << endl;
		
		// comb_gr_Sfit->SetPoint(ien,lgE[ien],Sfit);
		// comb_gr_Sfit->SetPointError(ien,0.,0.);
		// comb_gr_Sfit->SetName("Sfit");
		
		comb_gr_Scalc->SetPoint(ien,lgE[ien],Shist);
		comb_gr_Scalc->SetPointError(ien,0.,Sherr);
		comb_gr_Scalc->SetName("Scalc");
		
		comb_gr_lnnavg->SetPoint(ien,lgE[ien],log(navghist));
		comb_gr_lnnavg->SetPointError(ien,0.,navgherr/navghist);
		comb_gr_lnnavg->SetName("lnnavg");
		
		syst_gr_Scalc->SetPoint(ien,lgE[ien],Shist);
		syst_gr_Scalc->SetPointError(ien,0.1,0.1,syst_dn,syst_up);
		syst_gr_Scalc->SetName("syst_gr_Scalc");
	}
	
	cerr << "Writing out entropies ... " << endl;
	TFile * outfile = new TFile(Form("combined_entropy_per_energy.root"), "recreate");
	outfile->cd();
	// comb_gr_Sfit->Write();
	comb_gr_Scalc->Write();
	syst_gr_Scalc->Write();
	comb_gr_lnnavg->Write();
	// comb_gr_avgSfit->Write();
	comb_gr_avgScalc->Write();
}