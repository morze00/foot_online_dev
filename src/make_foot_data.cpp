#include"libs.hh"
#include <TMathBase.h>
#define LENGTH(x) (sizeof x / sizeof *x)
#define NDETS LENGTH(foot_id)
using namespace std;
namespace{//don't change this!!!
    int foot_id[] = {15, 16, 1, 2, 13, 4, 11, 6, 7, 12, 9, 10};
}
//Histogram binning and ranges
int Nbins = 1000; double ymin = -500; double ymax = 500;
const double    FOOT_LENGTH   = 96.;//mm
const double    NSIGMA        = 4;
const int MAX_STAT_PEDS = 3e4;
const int MAX_STAT_ANA = -1;

//Clustering data structures
typedef std::pair<int, double> strip_data;
typedef std::vector<strip_data> cluster; 
typedef std::vector<cluster> foot_data; 

Double_t pedestal[NDETS][640];
Double_t sigma[NDETS][640];
Double_t sigma_fine[NDETS][640];

class AnaWrapper : public TObject
{
    public:
        AnaWrapper();
        virtual ~AnaWrapper();
        void Init();
        bool is_good_strip(UInt_t det, UInt_t strip);
        double get_cog(cluster c);//calculate center of gravity of a cluster
        double get_esum(cluster c);//calculate cluster sum
        void Check_Strip(UInt_t number, double energy, foot_data& fdata);
        void analyse(int firstEvent, int max_events, TChain * ch);
        void Make_FineSigmas(int threshold);
        void Make_Pedestals(int threshold);
        void Draw_Everything();

        //---------- Definition of main  histograms ----------------
        TH2F * h2_peds_raw[NDETS];
        TH2F * h2_peds_raw_clean[NDETS];
        TH2F * h2_cal_fine[NDETS];
        TH2F * h2_baseline[NDETS];
        TH2F * h2_baseline_clean[NDETS];
        TH1D * h1_baseline[NDETS];
        TH1D * h1_peds[NDETS];
        TH1D * h1_sigma_raw[NDETS];
        TH1D * h1_sigma_fine[NDETS];
        TH1D * h1_cluster_e[NDETS];
        TH1I * h1_cluster_size[NDETS];
        TH2D * h2_cluster_e_vs_cog[NDETS];

        //---------- Definition of canvases
        TCanvas * canvas_raw_peds;
        TCanvas * canvas_raw_sigma;
        TCanvas * canvas_baseline;
        TCanvas * canvas_fine;
        TCanvas * canvas_cluster_energy;
        TCanvas * canvas_cluster_size;
        TCanvas * canvas_cluster_e_vs_cog;
        bool is_plot;
        bool is_good;
};

AnaWrapper::AnaWrapper()
    : is_plot(false)
      , is_good(false)
{
}

AnaWrapper::~AnaWrapper() {}

void AnaWrapper::Init()
{
    TString hname;
    for(int i=0; i<NDETS; i++)
    {
        hname.Form("hraw%d",foot_id[i]);
        h2_peds_raw[i]     = new TH2F(hname.Data(),hname.Data(),640,1,641,1000,0,1000);

        hname.Form("h%d",foot_id[i]);
        h2_peds_raw_clean[i]  = new TH2F(hname.Data(),hname.Data(),640,1,641,1000,0,1000);


        hname.Form("h2_cal_fine_FOOT%d",foot_id[i]);
        h2_cal_fine[i]     = new TH2F(hname.Data(),hname.Data(),640,1,641,Nbins,ymin,ymax);

        hname.Form("h2_baseline_FOOT%d",foot_id[i]);
        h2_baseline[i]     = new TH2F(hname.Data(),hname.Data(), 640,1,641,200,-50,50);

        hname.Form("h2_baseline_clean_FOOT%d",foot_id[i]);
        h2_baseline_clean[i] = new TH2F(hname.Data(),hname.Data(), 640,1,641,500,-100,100);

        hname.Form("h1_cluster_E_FOOT%d",foot_id[i]);
        h1_cluster_e[i]    = new TH1D(hname.Data(),hname.Data(), 2500, -500, 2000);

        hname.Form("h1_cluster_size_FOOT%d",foot_id[i]);
        h1_cluster_size[i] = new TH1I(hname.Data(),hname.Data(), 10, 0, 10);

        hname.Form("h1_cluster_e_vs_cog_FOOT%d",foot_id[i]);
        h2_cluster_e_vs_cog[i] = new TH2D(hname.Data(),hname.Data(),640,1,641,200,-10,100);
    }
    canvas_raw_peds = new TCanvas("canvas_raw_peds","canvas_raw_peds",1800,1200);
    canvas_raw_peds->Divide(4,3);

    canvas_raw_sigma = new TCanvas("canvas_raw_sigma","canvas_raw_sigma",1800,1200);
    canvas_raw_sigma->Divide(4,3);

    canvas_baseline = new TCanvas("canvas_baseline","canvas_baseline",1800,1200);
    canvas_baseline->Divide(4,3);


    canvas_fine   = new TCanvas("fine","fine",1800,1200);
    canvas_fine->Divide(4,3);

    canvas_cluster_energy   = new TCanvas("cluster_energy","cluster_energy",1800,1200);
    canvas_cluster_energy->Divide(4,3);

    canvas_cluster_size   = new TCanvas("cluster_size","cluster_size",1800,1200);
    canvas_cluster_size->Divide(4,3);

    canvas_cluster_e_vs_cog   = new TCanvas("cluster_e_vs_cog","cluster_e_vs_cog",1800,1200);
    canvas_cluster_e_vs_cog->Divide(4,3);
}

//list of bad or dead strips
bool AnaWrapper::is_good_strip(UInt_t det, UInt_t strip)
{
    switch(det){
        case 1:
            if(
                    strip==94 || (strip>506 && strip<513)
                    || (strip >528 && strip<532)
              ) return false;
        case 2:
            if(
                    strip==229 
              ) return false;
        case 7:
            if(
                    strip==344 || strip==597 || strip==598
              ) return false;

        case 12:
            if(
                    strip==232
              ) return false;
        case 15:
            if(
                    strip==204 || strip==203
              ) return false;
        case 16:
            if(
                    strip==338
              ) return false;
    }
    //if((strip%64)>62 || (strip%64)<3) return false;
    return true;
}

double AnaWrapper::get_cog(cluster c)//calculate center of gravity of a cluster
{
    double value = 0;  double esum = 0.;
    for(auto & s: c){
        value  += (s.first * s.second);
        esum += s.second;
    }
    value /= esum;//center of gravity
    return value;
}

double AnaWrapper::get_esum(cluster c)//calculate cluster sum
{
    double esum = 0.;
    for(auto & s: c){
        esum += s.second;
    }
    return esum;
}

void AnaWrapper::Check_Strip(UInt_t number, double energy, foot_data& fdata)
{
    //cout << "\n\n---------------------------------------";
    //cout << "\nChecking strip: " << number << "\t Energy: " << energy << endl;
    strip_data strip = std::make_pair(number,energy);
    cluster clust;
    if(fdata.size()==0)//no cluster yet, create new
    {
        clust.push_back(strip);
        fdata.push_back(clust);
        //cout << "\n\t New cluster is created for this strip";
        return;
    }
    cluster    this_clust = fdata.back();
    strip_data this_strip = this_clust.back();
    if(abs(strip.first-this_strip.first)<2)//neighbour found 
    {
        //cout << "\n\tStrip belong to exisitng cluster! Adding it...";
        fdata.back().push_back(strip);
        return;
    }
    else
    {
        //cout << "\n\tStrip is a new cluster! Making it...";
        clust.clear();
        clust.push_back(strip);
        fdata.push_back(clust);
        return;
    }
}

void AnaWrapper::Make_FineSigmas(int threshold)
{
    TF1* foo = new TF1("foo","gaus",-100,100);
    for(int i=0; i<NDETS; i++)
    {
        for(int xbin=1; xbin<=h2_baseline[i]->GetNbinsX(); ++xbin)
        {
            for(int ybin=1; ybin<=h2_baseline[i]->GetNbinsY(); ++ybin){
                double binContent = h2_baseline[i]->GetBinContent(xbin, ybin);
                if(binContent<threshold){
                    h2_baseline_clean[i]->SetBinContent(xbin, ybin,0.0); 
                }
                else{
                    h2_baseline_clean[i]->SetBinContent(xbin, ybin, binContent) ;
                }
            }
        }
        h2_baseline_clean[i]->FitSlicesY(foo,1,640,0,"QNR",0);
        h1_baseline[i] = (TH1D*)gDirectory->Get(Form("h2_baseline_clean_FOOT%d_1",foot_id[i]))->Clone(Form("h1_baseline_%d",foot_id[i]));
        h1_sigma_fine[i]  = (TH1D*)gDirectory->Get(Form("h2_baseline_clean_FOOT%d_2",foot_id[i]))->Clone(Form("h1_sigma_fine_%d",foot_id[i]));
        for(int j=0; j<640; j++){
            sigma_fine[i][j]    = h1_sigma_fine[i]->GetBinContent(j+1);
            cout << "\nFOOT: " << foot_id[i] << "\tStrip: " << j  << "\tFine sigma: " << sigma_fine[i][j];
        }
    }
    return;
}

void AnaWrapper::Make_Pedestals(int threshold)
{
    TF1* foo = new TF1("foo","gaus",50,700);
    for(int i=0; i<NDETS; i++)
    {
        for(int xbin=1; xbin<=h2_peds_raw[i]->GetNbinsX(); ++xbin){
            for(int ybin=1; ybin<=h2_peds_raw[i]->GetNbinsY(); ++ybin){
                double binContent = h2_peds_raw[i]->GetBinContent(xbin, ybin);
                if(binContent<threshold){
                    h2_peds_raw_clean[i]->SetBinContent(xbin, ybin,0.0); 
                }
                else{
                    h2_peds_raw_clean[i]->SetBinContent(xbin, ybin, binContent) ;
                }
            }
        }
        h2_peds_raw_clean[i]->FitSlicesY(foo,1,640,0,"QNR",0);
        h1_peds[i]      = (TH1D*)gDirectory->Get(Form("h%d_1",foot_id[i]))->Clone(Form("h1_peds_%d",  foot_id[i]));
        h1_sigma_raw[i] = (TH1D*)gDirectory->Get(Form("h%d_2",foot_id[i]))->Clone(Form("h1_sigma_raw_%d",foot_id[i]));
        for(int j=0; j<640; j++){
            pedestal[i][j] = h1_peds[i]->GetBinContent(j+1);
            sigma[i][j]    = h1_sigma_raw[i]->GetBinContent(j+1);
            cout << "\nFOOT: " << foot_id[i] << "\tStrip: " << j 
                << "\t Ped: " << pedestal[i][j] << "\tSig: " << sigma[i][j];
        }
    }
    return;
}

void AnaWrapper::Draw_Everything()
{
    //------- Drawing everything-------
    for(int i=0; i<NDETS; i++)
    {
        canvas_raw_peds->cd(i+1);
        h2_peds_raw[i]->Draw("colz");
        //h2_peds_raw_clean[i]->Draw("colz");
        h1_peds[i]->SetMarkerStyle(kFullCircle);
        h1_peds[i]->SetMarkerSize(0.1);
        h1_peds[i]->SetMarkerColor(kRed);
        h1_peds[i]->SetLineColor(kRed);
        h1_peds[i]->Draw("same");

        canvas_baseline->cd(i+1);
        gPad->SetLogz();
        h2_baseline[i]->Draw("colz");
        h1_baseline[i]->SetMarkerStyle(kFullCircle);
        h1_baseline[i]->SetMarkerSize(0.1);
        h1_baseline[i]->SetMarkerColor(kBlue);
        h1_baseline[i]->SetLineColor(kBlue);
        h1_baseline[i]->Draw("same");


        canvas_raw_sigma->cd(i+1);
        h1_sigma_raw[i]->SetMarkerStyle(kFullCircle);
        h1_sigma_raw[i]->SetMarkerSize(0.2);
        h1_sigma_raw[i]->SetMarkerColor(kBlue);
        h1_sigma_raw[i]->SetLineColor(kBlue);
        h1_sigma_raw[i]->GetYaxis()->SetRangeUser(-2,10);
        h1_sigma_raw[i]->GetXaxis()->SetTitle("Strip No.");
        h1_sigma_raw[i]->GetYaxis()->SetTitle("ADC sigma");
        //h1_sigma_raw[i]->SetTitle("Sigmas before baseline correction");
        h1_sigma_raw[i]->Draw();

        h1_sigma_fine[i]->SetMarkerStyle(kFullCircle);
        h1_sigma_fine[i]->SetMarkerSize(0.2);
        h1_sigma_fine[i]->SetMarkerColor(kRed);
        h1_sigma_fine[i]->SetLineColor(kRed);
        h1_sigma_fine[i]->Draw("same");

    }
    return;
}

void AnaWrapper::analyse(int firstEvent, int max_events, TChain * ch)
{
    TApplication* theApp = new TApplication("App", 0, 0);
    Init();
    //------ Define input tree data for all detectors ----------
    UInt_t TRIGGER;
    UInt_t FOOT[NDETS];//multiplicity of foot detectors 
    UInt_t FOOTE[NDETS][640];
    UInt_t FOOTI[NDETS][640];
    UInt_t TPATv[1];
    ch->SetBranchAddress("TRIGGER",&TRIGGER);
    //ch->SetBranchAddress("TPATv",TPATv);
    for(int i=0; i<NDETS; i++){
        TString bname   = Form("FOOT%d",  foot_id[i]);
        TString bname_E = Form("FOOT%dE", foot_id[i]);
        TString bname_I = Form("FOOT%dI", foot_id[i]);
        ch->SetBranchAddress(bname.Data(),&FOOT[i]);
        ch->SetBranchAddress(bname_E.Data(),FOOTE[i]);
        ch->SetBranchAddress(bname_I.Data(),FOOTI[i]);
    }
    //----- Collect pedestals data ------
    int Nevents = ch->GetEntries();
    if(max_events>0) Nevents = max_events; 
    int stat=0;
    int det_mul=0;
    cout << "\n-- Analysing raw pedestals " << endl;
    cout << "\n-- Max Stat: " << MAX_STAT_PEDS << "\n\n";
    for(int ev=0; ev<Nevents; ev++){
        cout << "\r-- Event # : " << ev << flush;
        ch->GetEntry(ev);
        int det_mul=0;
        for(int i=0; i<NDETS; i++){
            if(FOOT[i]>600) det_mul++;
        }
        //cout << "Seen coincident detectors: " << det_mul << "\n";
        if(det_mul!=NDETS) continue;//Use correct time stitching!!!
        for(int f=0; f<NDETS; f++){
            for(int  j=0 ; j<640 ; j++){
                h2_peds_raw[f]->Fill(FOOTI[f][j],FOOTE[f][j]); 
            }
        }
        stat++;
        if(stat==MAX_STAT_PEDS) break;
    }
    Make_Pedestals(50);

    //-------- Make baseline histo to get fine sigmas -------
    double  asic_offset[10];
    double  signal = 0;
    int     counter_asic =0;
    int     stat_baseline =0;
    cout << "\n-- Analysing fine baseline " << endl;
    for(int ev=0; ev<Nevents; ev++)
    {
        ch->GetEntry(ev);
        for(int f=0; f<NDETS; f++){//loop FOOTs
            cout << "\r-- Event # : " << ev << flush;
            //------------ Calculating correction for individual asics ---------
            stat=0; counter_asic=0;
            for(int i=0; i<10; i++){  asic_offset[i]=0.; }//reset asic baselines
            for(int i=0; i<640; i++){
                signal = FOOTE[f][i] - pedestal[f][i];
                if(fabs(signal) < (6 * sigma[f][i]) &&
                        is_good_strip(foot_id[f],FOOTI[f][i]) && //no bad strips
                        ((FOOTI[f][i])%64)<63 && ((FOOTI[f][i])%64)>1)
                {
                    stat++;
                    asic_offset[counter_asic] += signal;
                }
                if((FOOTI[f][i]%64)==0){//switch to next asic
                    asic_offset[counter_asic] /= stat;
                    counter_asic++;  stat=0;
                }
            }
            //-------- Applying fine correction and fill baseline histograms ----------
            counter_asic=0;
            stat=0;
            for(int i=0; i<640; i++)
            {
                if((FOOTI[f][i]%64) == 1 && FOOTI[f][i]>1) counter_asic++; 
                signal = FOOTE[f][i] - pedestal[f][i] - asic_offset[counter_asic];
                //if(fabs(signal) < (6 * sigma[f][i])){
                h2_baseline[f]->Fill(FOOTI[f][i],signal);
                //}
            }
        }//end detector loop
        stat_baseline++; if(stat_baseline==MAX_STAT_PEDS) break;
    }//end evetloop
    Make_FineSigmas(50);

    //--------- Final analysis 
    TFile* outfile = new TFile("output.root","Recreate","Write");
    TTree *tree = new TTree("tree","Tree with vectors of clusters");
    std::vector<int> id_foot;
    std::vector<int> size_foot;
    std::vector<double> e_foot;
    std::vector<double> cog_foot;
    std::vector<double> eta_foot;

    tree->Branch("id_foot",&id_foot);
    tree->Branch("e_foot",&e_foot);
    tree->Branch("cog_foot",&cog_foot);
    tree->Branch("size_foot",&size_foot);
    tree->Branch("eta_foot",&eta_foot);

    foot_data fdata[NDETS];
    cout << "\n\n-- Final analysis \n\n";
    if(MAX_STAT_ANA<0) Nevents=ch->GetEntries();
    else Nevents = MAX_STAT_ANA;
    for(int ev=0; ev<Nevents; ev++)
    {
        cout << "\r-- Event # : " << ev << flush;
        ch->GetEntry(ev);

        id_foot.clear();
        size_foot.clear();
        e_foot.clear();
        cog_foot.clear();
        eta_foot.clear();

        for(int f=0; f<NDETS; f++)//loop over all foots
        {
            fdata[f].clear(); 
            stat=0; counter_asic=0;
            for(int i=0; i<10; i++){  asic_offset[i]=0.; }//reset asic baselines
            for(int i=0; i<640; i++){
                signal = FOOTE[f][i] - pedestal[f][i];
                if(fabs(signal) < (6 * sigma[f][i]) &&
                        is_good_strip(foot_id[f],FOOTI[f][i]) && //no bad strips
                        ((FOOTI[f][i])%64)<63 && ((FOOTI[f][i])%64)>1)
                {
                    stat++;
                    asic_offset[counter_asic] += signal;
                }
                if((FOOTI[f][i]%64)==0){//switch to next asic
                    asic_offset[counter_asic] /= stat;
                    counter_asic++;  stat=0;
                }
            }
            //--- Fill cluster data using accurate sigmas 
            counter_asic=0;
            for(int i=0; i<640; i++)
            {
                if((FOOTI[f][i]%64) == 1 && FOOTI[f][i]>1) counter_asic++;
                signal = FOOTE[f][i] - pedestal[f][i] - asic_offset[counter_asic];
                if(signal>(NSIGMA * sigma_fine[f][i]) && 
                        is_good_strip(foot_id[f],FOOTI[f][i]))
                {
                    Check_Strip(FOOTI[f][i], signal, fdata[f]);
                }
            }
        }//end loop detectors
        //Filling output tree
        for(int f=0; f<NDETS; f++)
        {
            for(auto & clust: fdata[f])
            {
                id_foot.push_back(foot_id[f]);
                e_foot.push_back(get_esum(clust));
                cog_foot.push_back(get_cog(clust));
                size_foot.push_back(clust.size());
            }
        }
        tree->Fill();
    }//end of eventloop
    tree->AutoSave();
    Draw_Everything();

    outfile->Close();
    theApp->Run();
    return;
}

int main(Int_t argc, Char_t* argv[])
{
    gRandom = new TRandom3();
    gRandom->SetSeed(0);
    gROOT->Macro("rootlogon.C");
    gStyle->SetPalette(kRainBow);

    AnaWrapper ana;

    TChain * ch = new TChain("h101");
    //ch->Add("../roots_foots/main0131_0041.root");
    ch->Add("../roots_foots/run131_0001_stitched.root");
    //ch->Add("/Users/vpanin/Desktop/GSI/Experiments/S522/analysis/Tracking/data_unpacked/");
    ana.analyse(0,-1,ch);
    return 0;
}
