#include "ext_data_clnt.hh"
#include "ext_data_struct_info.hh"
#include "ext_data_client.h"
#include "ext_h101_unpack.h"
#include "ext_h101_foot.h"
#include"libs.hh"
#define LENGTH(x) (sizeof x / sizeof *x)
#define NDETS LENGTH(foot_id)
using namespace std;
typedef struct EXT_STR_h101_t
{
	EXT_STR_h101_unpack_t unpack;
	EXT_STR_h101_FOOT_onion foot;
	//EXT_STR_h101_FOOT foot;
} EXT_STR_h101;
namespace{//don't change this!!!
	int foot_id[] = {1, 2, 13, 4, 11, 6, 7, 12};
}
Double_t pedestal[16][640];
Double_t sigma[16][640];
TH2F * h2_peds[16];
TH2F * h2_calib[16];
TH2F * h2_calib_fine[16];
TH1D * h1_peds[16];
TH1D * h1_sigma[16];
TH1D * h1_sigma_fine[16];
auto pedfilename = "/u/land/r3broot/202402_s091_s118/R3BParams_s091_s118/macros/exp/online/foot/event_display_dev/pedestal.dat";

bool is_good_strip(UInt_t det, UInt_t strip)
{
	switch(det){
		case 13:
			if(
					strip == 203 ||	strip==204 || strip==205
			  ) return false;
	}
	//if((strip%64)>61 || (strip%64)<4) return false;//ASIC edge
	if((strip%64)<3) return false;//ASIC edge
	return true;
}
void Draw_Asic_Borders(double ymin, double ymax)
{
	for(int i_asic=1; i_asic<10; i_asic++)
	{
		TLine* l = new TLine(64*i_asic,ymin,64*i_asic,ymax);
		l->SetLineStyle(7);
		l->SetLineWidth(1);
		l->SetLineColor(13);
		l->Draw("same");
	}
}
TH2F Get_Clean_TH2F(TH2F* h2_in, int thr)
{
	TString hname = TString(h2_in->GetTitle()) + "_clean";
	TH2F h2_out_clean(hname.Data(),hname.Data(),
			h2_in->GetNbinsX(),
			h2_in->GetXaxis()->GetXmin(), h2_in->GetXaxis()->GetXmax(),
			h2_in->GetNbinsY(),
			h2_in->GetYaxis()->GetXmin(), h2_in->GetYaxis()->GetXmax());
	for(int xbin=1; xbin<=h2_in->GetNbinsX(); ++xbin){//clean up low stat bins
		for(int ybin=1; ybin<=h2_in->GetNbinsY(); ++ybin){
			double binContent = h2_in->GetBinContent(xbin, ybin);
			h2_out_clean.SetBinContent(xbin, ybin,(binContent < thr) ? 0.0 : binContent); 
		}
	}
	return h2_out_clean;
}
void Calculate_Pedestals(int threshold)
{
	cout << "\n-- Calculating pedestals";
	std::ofstream fout(pedfilename);
	if (!fout.is_open()){
		std::cerr << "Error opening file: " << pedfilename << std::endl;
		return;
	}
	TF1 foo("foo","gaus",50,700);
	for(int i=0; i<16; i++)
	{
		if(h2_peds[i]->Integral()<1) continue;//do not remove this cond
		auto h2_peds_clean = Get_Clean_TH2F(h2_peds[i], threshold);
		gErrorIgnoreLevel = kFatal;//to suppress printout info on the terminal
		h2_peds_clean.FitSlicesY(&foo,1,h2_peds_clean.GetNbinsX(),0,"QNR",0);
		gErrorIgnoreLevel = kPrint; // Reset to default level
		h1_peds[i] =  dynamic_cast<TH1D*> (gDirectory->Get(TString(h2_peds_clean.GetTitle()) +"_1"));
		h1_sigma[i] = dynamic_cast<TH1D*> (gDirectory->Get(TString(h2_peds_clean.GetTitle()) +"_2"));
		fout << i << std::endl;
		for(int j=0; j<640; j++){
			pedestal[i][j] = h1_peds[i]->GetBinContent(j+1);
			sigma[i][j]    = h1_sigma[i]->GetBinContent(j+1);
			fout << j+1 << "  " <<  pedestal[i][j] << "  " <<  sigma[i][j] << std::endl;
		}
	}
	std::cout << "\n-- Writing pedestals and sigmas into file " << pedfilename << std::endl;
	fout.close();
}
void Calculate_Fine_Sigmas(int threshold)
{
	cout << "\n-- Calculating fine sigmas";
	TF1 foo("foo","gaus",-200,400);
	for(int i=0; i<16; i++)
	{
		if(h2_calib_fine[i]->Integral()<1) continue;//do not remove this cond
		auto h2_clean = Get_Clean_TH2F(h2_calib_fine[i], threshold);
		gErrorIgnoreLevel = kFatal;//to suppress printout info on the terminal
		h2_clean.FitSlicesY(&foo,1,h2_clean.GetNbinsX(),0,"QNR",0);
		gErrorIgnoreLevel = kPrint; // Reset to default level
		h1_sigma_fine[i] = dynamic_cast<TH1D*> (gDirectory->Get(TString(h2_clean.GetTitle()) +"_2"));
	}
}
void Draw_Pedestals(TCanvas * can, Bool_t is_fine_sigma)
{
	std::cout << "\n-- Plotting pedestals for all FOOTs\n";
	can->Divide(4,4);
	for(int i=0;i<NDETS;++i)
	{
		int d = foot_id[i]-1;
		if(h2_peds[d]->Integral()<1) continue;//don't remove this cond!
		can->cd(i*2+1);
		h2_peds[d]->Draw("colz");
		Draw_Asic_Borders(h2_peds[d]->GetYaxis()->GetXmin(), 
				h2_peds[d]->GetYaxis()->GetXmax());
		h1_peds[d]->SetMarkerStyle(kFullCircle);
		h1_peds[d]->SetMarkerSize(0.1);
		h1_peds[d]->SetMarkerColor(kRed);
		h1_peds[d]->SetLineColor(kRed);
		//h1_peds[d]->Draw("same");		
		can->cd(i*2+2);
		h1_sigma[d]->SetMarkerStyle(kFullCircle);
		h1_sigma[d]->SetMarkerSize(0.3);
		h1_sigma[d]->SetMarkerColor(kBlue);
		h1_sigma[d]->SetLineColor(kBlue);
		h1_sigma[d]->GetYaxis()->SetRangeUser(0,15);
		h1_sigma[d]->GetXaxis()->SetTitle("Strip No.");
		h1_sigma[d]->GetYaxis()->SetTitle("ADC sigma");
		TString htitle;
		htitle.Form("Sigmas for FOOT%d",foot_id[i]);
		h1_sigma[d]->SetTitle(htitle.Data());
		h1_sigma[d]->Draw();
		if(is_fine_sigma)
		{	
			h1_sigma_fine[d]->SetMarkerStyle(kFullCircle);
			h1_sigma_fine[d]->SetMarkerSize(0.3);
			h1_sigma_fine[d]->SetMarkerColor(kRed);
			h1_sigma_fine[d]->SetLineColor(kRed);
			h1_sigma_fine[d]->Draw("same");
		}
		Draw_Asic_Borders(0, 15);
	}
}
bool Read_Pedestals()
{
	ifstream pedfile(pedfilename,ifstream::in);
	if ( !pedfile.is_open() ){
		cout << "\n\nERROR: Cannot open pedestal file! Try to run first with the option --pedestal \n\n";
		return false;
	}
	std::cout << "\n-- Reading pedestals from the file " << pedfilename;
	int Nlines = 0;
	int strip_no, det_no;
	std::string line;
	while( std::getline(pedfile, line) ) Nlines++;
	cout << "\n-- Number of lines in the pedestal file: " << Nlines << "\n\n";
	pedfile.clear();
	pedfile.seekg( 0, std::ios::beg );
	for(int d=0; d<Nlines/641; d++){
		pedfile >>det_no;
		for(Int_t i=0 ; i<640 ; i++){
			pedfile >> strip_no >>  pedestal[det_no][i] >> sigma[det_no][i];
		}
	}
	pedfile.close();
	return true;
}
void Make_Histograms()
{
	for(int i=0; i<16; ++i){
		TString hname = TString::Format("FOOT%d raw",i+1);
		h2_peds[i] = new TH2F(hname.Data(),hname.Data(),640,1,641,1000,0,1000);
		hname = TString::Format("FOOT%d calib",i+1);
		h2_calib[i] = new TH2F(hname.Data(),hname.Data(),640,1,641,200,-100,100);
		hname = TString::Format("FOOT%d calib fine",i+1);
		h2_calib_fine[i] = new TH2F(hname.Data(),hname.Data(),640,1,641,200,-100,100);
		h2_peds[i]->SetDirectory(0);//do not store histo for http server
		h2_calib[i]->SetDirectory(0);
		h2_calib_fine[i]->SetDirectory(0);
	}
}
int main(Int_t argc, Char_t* argv[])
{
	Bool_t NeedHelp = kFALSE;
	Bool_t MakePedestals = kFALSE;
	Bool_t MakeFineSigma = kFALSE;
	Bool_t RunOnline = kFALSE;
	int max_events = -1; //default stat
	for (Int_t i = 0; i < argc; i++){
		if (strncmp(argv[i],"--max-events=",13) == 0){//if not set, unlimited
			max_events = atoi(argv[i]+13);
		}
		else if (strncmp(argv[i],"--pedestals",11) == 0){//determine/plot pedestals
			MakePedestals = kTRUE;
		}
		else if (strncmp(argv[i],"--fine-sigma",12) == 0){//pedestals with fine sigmas
			MakeFineSigma = kTRUE;
		}
		else if (strncmp(argv[i],"--online",8) == 0){//send data to http
			RunOnline = kTRUE;
			gROOT->SetBatch(kTRUE);
		}
		else if (strncmp(argv[i],"--help",6) == 0){//print usage info
			NeedHelp = kTRUE;
		}
	}
	if (NeedHelp || ((MakePedestals || MakeFineSigma) && RunOnline)){
		std::cout << "\nThe following options are supported:" << "\n\n";
		std::cout << " --max-events=NUMBER : desired statistics (default NUMBER = -1 (i.e. unlimited))\n" << std::endl;
		std::cout << " --pedestals : to flag pedestal run \n" << std::endl;
		std::cout << " --finei-sigma : caclulate fine sigmas (first run --pedestals) \n" << std::endl;
		std::cout << " --online : to send online data to http port \n" << std::endl;
		std::cout << " --help :  to print usage info \n" << std::endl;
		return 0;
	}
	std::cout << "\n-- Requested statistics: " << max_events << " events";
	std::cout << "\n-- Collecting data for ";
	if (MakePedestals) std::cout << "pedestals determination";
	else if(RunOnline) std::cout << "online monitor (http)";
	else if(MakeFineSigma) std::cout << "pedestals with fine sigmas";
	else cout << "the analysis";
	gROOT->Macro("rootlogon.C");
	gStyle->SetPalette(kBird);

	EXT_STR_h101 ucesb_struct;
	ext_data_clnt fClient;
	Bool_t status;
	std::ostringstream command;

	command << "$UCESB_DIR/../upexps/202205_s522/202205_s522 --stream=krpc003:8003 --allow-errors --input-buffer=300Mi --ntuple=RAW,FOOT,STRUCT,-";

	/* Fork off ucesb (calls fork() and pipe()) */
	FILE* fFd = popen(command.str().c_str(), "r");
	if (nullptr == fFd){
		cout << "\npopen() failed\n";
		return 0;
	}
	/* Connect to forked instance */
	status = fClient.connect(fileno(fFd));
	if (kFALSE == status){
		cout << "\next_data_clnt::connect() failed";
		cout << "\nucesb error: " << fClient.last_error();
		return 0;
	}
	uint32_t struct_map_success = 0;
	auto fStructInfo = new ext_data_struct_info();
	int ok;
	size_t fOffset = offsetof(EXT_STR_h101, foot);
	EXT_STR_h101_FOOT_ITEMS_INFO(ok, *fStructInfo, fOffset, EXT_STR_h101_FOOT, 0);

	if (!ok){
		perror("ext_data_struct_info_item");
		fprintf (stderr,"Failed to setup structure information.\n");
		exit(1);
	}
	status = fClient.setup(NULL, 0, fStructInfo, &struct_map_success, sizeof(ucesb_struct));
	if (status != 0){
		cout << "\next_data_clnt::setup() failed";
		cout << "\nUCESB error: " << fClient.last_error();
	}
	Make_Histograms();//creating empty histograms
	if(!MakePedestals){
		if(!Read_Pedestals()){
			std::cerr << "ERROR!  Cannot read pedestals! Try to run with the option '--pedestals'" <<  std::endl;
			return 0;
		}
	}
	TApplication* theApp = new TApplication("App", 0, 0);
	TCanvas * c_fine;
	TCanvas * c_raw;
	TCanvas * c_noped;
	if(RunOnline)
	{
		auto serv = new THttpServer("http:8891");//setup web display
		gROOT->SetWebDisplay("browser"); // or "cef" for embedded browser
		c_raw = new TCanvas("Raw signals","Raw signals",1000,500);
		c_raw->Divide(4,2);
		c_noped = new TCanvas("Subtracted pedestals","Subtracted pedestals",1000,500);
		c_noped->Divide(4,2);
		c_fine = new TCanvas("Fine baseline","Fine baseline",1000,500);
		c_fine->Divide(4,2);
		serv->Register("/", c_raw);
		serv->Register("/", c_noped);
		serv->Register("/", c_fine);
	}
	int ret;
	const void* raw;
	ssize_t raw_words;
	int Nevents=0;
	int asic_pitch = 32;//chunk of strips, 64 in 1 asic
	const int Nasic = 640/asic_pitch;
	double  asic_offset[Nasic];
	double  signal = 0;
	int counter_asic =0;
	int stat=0;
	double signal_raw;
	int strip;
	while(1)//Eventloop
	{
		ret = fClient.fetch_event(&ucesb_struct, sizeof(ucesb_struct));
		if (-1 == ret){
			cout << "\next_data_clnt::fetch_event() failed! \nUCESB error: " << fClient.last_error() << endl;
			return 0;
		}
		ret = fClient.get_raw_data(&raw, &raw_words);
		if (0 != ret){
			cout << "\nFailed to get raw data.\n";
			return 0;
		}
		for(int d=0;d<16;d++)
		{
			if(ucesb_struct.foot.FOOT[d]._ !=640) continue;
			stat=0; counter_asic=0;
			for(int i=0; i<Nasic; i++){  asic_offset[i]=0.; }//reset asic baselines
			for(int  j=0 ; j<ucesb_struct.foot.FOOT[d]._; j++)
			{
				strip = ucesb_struct.foot.FOOT[d].I[j];
				signal_raw = ucesb_struct.foot.FOOT[d].E[j];
				h2_peds[d]->Fill(strip, signal_raw); 
				if(MakePedestals && !MakeFineSigma) continue;
				//calculate baseline correction	
				signal = signal_raw - pedestal[d][j];
				h2_calib[d]->Fill(strip, signal);
				if(fabs(signal) < (6 * sigma[d][j]) && 
						is_good_strip(foot_id[d],strip)){
					stat++;
					asic_offset[counter_asic] += signal;
				}
				if((strip % asic_pitch)==0 && strip>1){//switch to next asic
					if(stat>1) asic_offset[counter_asic] /= stat;
					counter_asic++;  stat=0;
				}
			}
			//----- apply fine baseline correction
			counter_asic=0;  stat=0;
			for(int  j=0 ; j<ucesb_struct.foot.FOOT[d]._; j++)
			{
				strip = ucesb_struct.foot.FOOT[d].I[j];
				signal_raw = ucesb_struct.foot.FOOT[d].E[j];
				if((strip % asic_pitch) == 1 && strip>1) counter_asic++; 
				signal = signal_raw - pedestal[d][j] - asic_offset[counter_asic];
				h2_calib_fine[d]->Fill(strip, signal);
			}
		}
		Nevents++;
		if(max_events>0 &&  Nevents>=max_events) break;
		else if(RunOnline && (Nevents % 5000)==0)//update every N events
		{
			for(int i=0;i<NDETS;++i)
			{
				int d = foot_id[i]-1;
				c_raw->cd(i+1);
				h2_peds[d]->Draw("colz");
				Draw_Asic_Borders(h2_peds[d]->GetYaxis()->GetXmin(),
						h2_peds[d]->GetYaxis()->GetXmax());
				gPad->Update();
				c_raw->Update();
				c_noped->cd(i+1);
				h2_calib[d]->Draw("colz");
				Draw_Asic_Borders(h2_calib[d]->GetYaxis()->GetXmin(),
						h2_calib[d]->GetYaxis()->GetXmax());
				gPad->Update();
				c_noped->Update();
				c_fine->cd(i+1);
				h2_calib_fine[d]->Draw("colz");
				gPad->Update();
				c_fine->Update();
				Draw_Asic_Borders(h2_calib_fine[d]->GetYaxis()->GetXmin(),
						h2_calib_fine[d]->GetYaxis()->GetXmax());
				gPad->Update();
				c_fine->Update();
				gSystem->ProcessEvents();
			}
		}
	}
	if(MakePedestals || MakeFineSigma){
		auto * c1 = new TCanvas("c1","c1",1000,1000);
		Calculate_Pedestals(10);
		if(MakeFineSigma) Calculate_Fine_Sigmas(10);
		Draw_Pedestals(c1, MakeFineSigma);
	}
	else //display calib spectra online (canvas or browser)
	{
		auto * c2 = new TCanvas("c2","c2",1000,500);
		c2->Divide(4,2);
		for(int i=0;i<NDETS;++i)
		{
			int d = foot_id[i]-1;
			if(h2_calib_fine[d]->Integral()<1) continue;//don't remove this cond!
			if(RunOnline) c_fine->cd(i+1);//displayed in browser
			else c2->cd(i+1);//dispalyed on PC
			h2_calib_fine[d]->Draw("colz");
			Draw_Asic_Borders(h2_calib_fine[d]->GetYaxis()->GetXmin(),
					h2_calib_fine[d]->GetYaxis()->GetXmax());
		}
	}
	gPad->Update();
	theApp->Run(kTRUE);
	return 0;
}
