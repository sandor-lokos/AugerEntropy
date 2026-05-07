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

const int NFILES = 4;
const int NENERGIES = 4;
const int NENTRY = 20;
const int lim_entries[NENTRY] = {  100,  200,  300,  400,  500,
																	1000, 1200, 1500, 1750, 2000,
																	2250, 2500, 2750, 3000, 3250,
																	3500, 3750, 4000, 4500, 5000};

const char* filenames[NFILES][NENERGIES] = {
																						{ "mcfiles/conex_100_17.5.root",
																							"mcfiles/conex_100_18.0.root",
																							"mcfiles/conex_100_18.5.root",
																							"mcfiles/conex_100_19.0.root"},
																						{	"mcfiles/conex_400_17.5.root",
																							"mcfiles/conex_400_18.0.root",
																							"mcfiles/conex_400_18.5.root",
																							"mcfiles/conex_400_19.0.root"},
																						{	"mcfiles/conex_1600_17.5.root",
																							"mcfiles/conex_1600_18.0.root",
																							"mcfiles/conex_1600_18.5.root",
																							"mcfiles/conex_1600_19.0.root"},
																						{	"mcfiles/conex_5600_17.5.root",
																							"mcfiles/conex_5600_18.0.root",
																							"mcfiles/conex_5600_18.5.root",
																							"mcfiles/conex_5600_19.0.root"},
																					 };

// =====================================================
// Main macro
// =====================================================
void SaveMult()
{
	TH1D  * hMult[NFILES][NENERGIES];
  TFile * infile[NFILES][NENERGIES];
	
	for( int ien = 0 ; ien < NENERGIES ; ien++)
	{
		for( int ifile = 0 ; ifile < NFILES ; ifile++)
		{
			infile[ifile][ien] = TFile::Open(filenames[ifile][ien], "READ");
			if (!infile[ifile][ien] || infile[ifile][ien]->IsZombie()) {
				cerr << "ERROR: cannot open file " << filenames[ifile][ien] << endl;
				return;
			}
			cerr << "Input file: " << filenames[ifile][ien] << endl;
			
			TTree* tree = dynamic_cast<TTree*>(infile[ifile][ien]->Get("LeadingInteractions"));
			if (!tree) {
				cerr << "ERROR: tree 'LeadingInteractions' not found." << endl;
				infile[ifile][ien]->Close();
				return;
			}
			
			Int_t nInt = 0;
			Int_t mult[10000] = {0};

			tree->SetBranchAddress("nInt", &nInt);
			tree->SetBranchAddress("mult", mult);

			const Long64_t nEntries = tree->GetEntries();
			if (nEntries <= 0) {
				cerr << "ERROR: tree has no entries." << endl;
				infile[ifile][ien]->Close();
				return;
			}
			cerr << "NUMBER OF ENTRIES: " << nEntries << endl;			
			cerr << "Initializing the histogram..." << endl;
			hMult[ifile][ien] = new TH1D(Form("hMult_en%i_prim%i",ien,ifile), "N;Counts", 10000, -0.5, 10000.5 );
			hMult[ifile][ien]->SetDirectory(nullptr);
			
			for (int iEntry = 0; iEntry < nEntries; ++iEntry)
			{
				tree->GetEntry(iEntry);
				if (nInt > 0)
					hMult[ifile][ien]->Fill(mult[0]);
			}
			
			if (hMult[ifile][ien]->GetEntries() == 0)
			{
				cerr << Form("ERROR: histogram for energy %i, primary %i is empty.", ien, ifile) << endl;
				return;
			}
			
			// infile[ifile][ien]->Close();
		} // end of file loop
	} // end of energy loop
	
	cerr << "Writing out the normalized multiplicities ... " << endl;
	TFile * outfile = new TFile(Form("multiplicities_per_primaries.root"), "recreate");
	outfile->cd();
	for( int ien = 0 ; ien < NENERGIES ; ien++)
	{
		for( int ifile = 0 ; ifile < NFILES ; ifile++)
		{
			cerr << "\t file: " << ifile << ", energy: " << ien << endl;
			hMult[ifile][ien]->Scale(1./hMult[ifile][ien]->Integral());	
			hMult[ifile][ien]->Write();
			cerr << "\t Sanity check int(h): " << hMult[ifile][ien]->Integral() << endl;
		} // end of file loop
	} // end of energy loop
	
	outfile->Close();
	
	
	
	
	
  // TFile * infile[NFILES][NENERGIES];
	TH1D  * hMultEntryDep[NENTRY][NFILES][NENERGIES];
	
	for( int iNentr = 0 ; iNentr < NENTRY ; iNentr++)
	{
		for( int ien = 0 ; ien < NENERGIES ; ien++)
		{
			for( int ifile = 0 ; ifile < NFILES ; ifile++)
			{
				infile[ifile][ien] = TFile::Open(filenames[ifile][ien], "READ");
				if (!infile[ifile][ien] || infile[ifile][ien]->IsZombie()) {
					cerr << "ERROR: cannot open file " << filenames[ifile][ien] << endl;
					return;
				}
				cerr << "Input file: " << filenames[ifile][ien] << endl;
				
				TTree* tree = dynamic_cast<TTree*>(infile[ifile][ien]->Get("LeadingInteractions"));
				if (!tree) {
					cerr << "ERROR: tree 'LeadingInteractions' not found." << endl;
					infile[ifile][ien]->Close();
					return;
				}
				
				Int_t nInt = 0;
				Int_t mult[10000] = {0};

				tree->SetBranchAddress("nInt", &nInt);
				tree->SetBranchAddress("mult", mult);

				const Long64_t nEntries = lim_entries[iNentr]; //tree->GetEntries();
				if (nEntries <= 0) {
					cerr << "ERROR: tree has no entries." << endl;
					infile[ifile][ien]->Close();
					return;
				}
				cerr << "NUMBER OF ENTRIES: " << nEntries << endl;			
				cerr << "Initializing the histogram..." << endl;
				hMultEntryDep[iNentr][ifile][ien] = new TH1D(Form("hMult_en%i_prim%i_entry%i",ien,ifile,iNentr), "N;Counts", 10000, -0.5, 10000.5 );
				hMultEntryDep[iNentr][ifile][ien]->SetDirectory(nullptr);
				
				for (int iEntry = 0; iEntry < nEntries; ++iEntry)
				{
					tree->GetEntry(iEntry);
					if (nInt > 0)
						hMultEntryDep[iNentr][ifile][ien]->Fill(mult[0]);
				}
				
				if (hMultEntryDep[iNentr][ifile][ien]->GetEntries() == 0)
				{
					cerr << Form("ERROR: histogram for energy %i, primary %i is empty.", ien, ifile) << endl;
					return;
				}
				
				infile[ifile][ien]->Close();
			} // end of file loop
		} // end of energy loop
	} // end of entries loop
	cerr << "Writing out the normalized multiplicities ... " << endl;
	outfile->Close();
	
	TFile * outfile_entry = new TFile(Form("entrydep_multiplicities_per_primaries.root"), "recreate");
	outfile_entry->cd();
	for( int iNentr = 0 ; iNentr < NENTRY ; iNentr++)
	{
		for( int ien = 0 ; ien < NENERGIES ; ien++)
		{
			for( int ifile = 0 ; ifile < NFILES ; ifile++)
			{
				cerr << "\t file: " << ifile << ", energy: " << ien << endl;
				hMultEntryDep[iNentr][ifile][ien]->Scale(1./hMultEntryDep[iNentr][ifile][ien]->Integral());	
				hMultEntryDep[iNentr][ifile][ien]->Write();
				cerr << "\t Sanity check int(h): " << hMultEntryDep[iNentr][ifile][ien]->Integral() << endl;
			} // end of file loop
		} // end of energy loop
	} // end of entries loop
	cerr << "Writing out the normalized entry dependent multiplicities ... " << endl;
	outfile_entry->Close();
	
}