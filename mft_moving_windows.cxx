// download_qc_objects_mw.cxx
// David Grund, 2024

// cpp headers
#include <iostream>
#include <string>
#include <ctime>
#include <vector>
// root headers
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TH1.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
// o2 headers
#include "CCDB/CcdbApi.h"

o2::ccdb::CcdbApi api_ccdb;
o2::ccdb::CcdbApi api_qcdb;

std::string str_sor = "STF";
std::string str_eor = "EOR";
std::string mw_path = "qc_async/MFT/MO/Tracks/mw/";
bool rewrite_root = false;

long get_timestamp (int run, std::string str) 
{
  long timestamp = -1;
  std::map<std::string, std::string> hdRCT = api_ccdb.retrieveHeaders("RCT/Info/RunInformation", map<std::string, std::string>(), run);
  const auto startRCT = hdRCT.find(str);
  if (startRCT != hdRCT.end()) {
    timestamp = stol(startRCT->second);
    std::cout << str << " found, timestamp: " << timestamp << "\n";
  } else {
    std::cout << str << " not found in headers!" << "\n";
  }
  return timestamp;
}

// https://stackoverflow.com/questions/63419684/convert-epoch-milliseconds-to-utc-date-and-time-using-c
std::string timestamp_to_date_string (long ts_ms, bool hours_only = false, bool verbose = false)
{
  long ts_sec = ts_ms / 1000; // from miliseconds (e.g. 1719668288021) to seconds
  std::time_t ts_as_time_t = ts_sec; // convert from long to time_t
  auto ts_as_tm = std::localtime(&ts_as_time_t);
  char buff[80];
  if (hours_only) std::strftime(buff, sizeof(buff), "%H:%M:%S", ts_as_tm);
  else            std::strftime(buff, sizeof(buff), "%m/%d/%Y %H:%M:%S", ts_as_tm);
  std::string date(buff);
  if(verbose) {
    std::cout << "\ntimestamp to date string:\n"
      << "ts (ms): " << ts_ms << " (input)\n"
      << "ts (s): " << ts_sec << "\n"
      << "as time_t: " << ts_as_time_t << "\n"
      << "as tm: " << ts_as_tm << "\n";
    if(hours_only) std::cout << "as hh:mm:ss: " << date << "\n\n";
    else           std::cout << "as mm/dd/yyyy hh:mm:ss: " << date << "\n\n";
  }
      
  return date;
}

std::vector<long> get_validity_from_name (std::string s, std::string delimiter = "_", bool verbose = false) 
{
  std::vector<std::string> substrings;
  size_t pos = 0;
  std::string substr;
  while ((pos = s.find(delimiter)) != std::string::npos) 
  {
    substr = s.substr(0, pos);
    substrings.push_back(substr);
    s.erase(0, pos + delimiter.length());
  }
  substrings.push_back(s);
  if(verbose) for(auto v : substrings) std::cout << v << "\n";
  std::vector<long> validity;
  validity.push_back(std::stol(substrings[0]));
  validity.push_back(std::stol(substrings[1]));
  return validity;
}

template <typename TH>
TH* download_histo(int run, std::string pass, std::string hname, long ts, std::vector<long>* val = NULL, bool verbose = false)
{
  std::string s_run = std::to_string(run);
  // create metadata:
  std::map<std::string, std::string> metadata;
  metadata["RunNumber"] = s_run;
  metadata["PassName"] = pass;
  // print the used metadata:
  std::string date = timestamp_to_date_string(ts);
  std::cout << "\nDownloading for: " << date << " (" << ts << ")\n";
  if(verbose) {
    std::cout << "metadata requirements: \n -> ";
    std::map<std::string,std::string>::iterator it;
    for (it = metadata.begin(); it != metadata.end(); it++) std::cout << it->first << "=" << it->second << ", ";
    std::cout << "\n";
  }
  // retrieve the histogram
  TH* h = api_qcdb.retrieveFromTFileAny<TH>(hname,metadata,ts);
  // retrieve its headers
  if(h) {
    std::map<std::string, std::string> headers = api_qcdb.retrieveHeaders(hname,metadata,ts);
    // print all the headers
    if(verbose) {
      std::map<std::string, std::string>::iterator it;
      for (it = headers.begin(); it != headers.end(); it++) std::cout << it->first << "\t" << it->second << "\n";
    }
    // extract the validities
    if(val) {
      long valid_from = std::stol(headers.at("Valid-From"));
      long valid_until = std::stol(headers.at("Valid-Until"));
      val->clear();
      val->push_back(valid_from);
      val->push_back(valid_until);
    }
  } else {
    if(val) {
      val->clear();
      val->push_back(0);
      val->push_back(0);
    }
  }
  return h;
}

void mft_moving_windows (int run, std::string pass, std::string hname, bool verbose = false)
{
  gSystem->Exec("mkdir -p root_files/");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  // connect to ccdb, qcdb, qcdbmc    
  api_ccdb.init("alice-ccdb.cern.ch");
  api_qcdb.init("ali-qcdb-gpn.cern.ch:8083");

  long ts_sor = get_timestamp(run, str_sor);
  long ts_eor = get_timestamp(run, str_eor);

  std::string fname = Form("root_files/%i_%s.root", run, pass.data());

  // download the objects
  bool fexists = !gSystem->AccessPathName(fname.data());
  if(fexists && !rewrite_root) {
    std::cout << fname << " already downloaded -> skipping\n";
  } else {
    TFile* f = new TFile(fname.data(), "recreate");
    TObjArray* arr = new TObjArray();
    long current_ts = ts_sor;
    bool next_obj_exists = true;
    while (next_obj_exists)
    {
      std::vector<long> val;
      TH1F* h = download_histo<TH1F>(run, pass, mw_path+hname, current_ts, &val);
      if(h) {
        int n_bins = h->GetNbinsX();
        float x_low = h->GetBinLowEdge(1);
        float x_upp = h->GetBinLowEdge(n_bins+1);
        if(verbose) std::cout << "n_bins: " << n_bins << ", x_low: " << x_low << ", x_upp: " << x_upp << "\n";
        long valid_from = val[0];
        long valid_until = val[1];
        TH1F *h_new = new TH1F(Form("%li_%li", valid_from, valid_until), h->GetTitle(), n_bins, x_low, x_upp);
        h_new->GetXaxis()->SetTitle(h->GetXaxis()->GetTitle());
        h_new->GetYaxis()->SetTitle(h->GetYaxis()->GetTitle());
        for(int i = 1; i <= n_bins; i++) h_new->SetBinContent(i, h->GetBinContent(i));
        arr->Add(h_new);
        current_ts = valid_until + 1;
        if (current_ts >= ts_eor) next_obj_exists = false;
      } else {
        std::cout << "The object does not exist!\n";
        next_obj_exists = false;
      }
    }
    f->cd();
    arr->Write(hname.data(), TObject::kSingleKey);
    f->Write("",TObject::kWriteDelete);
    f->Close();
  }

  // plot the objects
  gSystem->Exec(Form("mkdir -p plots/%i_%s/", run, pass.data()));
  TFile* f = TFile::Open(fname.data());
  TObjArray *arr = (TObjArray*) f->Get(hname.data());

  // get the maximum
  float y_max = 0;
  for (int i = 0; i < arr->GetEntries(); i++)
  {
    TH1F* h_curr = (TH1F*)arr->At(i);
    float curr_max = h_curr->GetMaximum() * 1.05;
    if(curr_max > y_max) y_max = curr_max;
  }

  for (int i = 0; i < arr->GetEntries(); i++)
  {
    TH1F* h_curr = (TH1F*)arr->At(i);

    std::vector<long> val_curr = get_validity_from_name(h_curr->GetName());
    if(verbose) {
      std::cout << "Histogram validity: \n"
        << " * from: " << timestamp_to_date_string(val_curr.at(0)) << " (" << val_curr.at(0) << ")\n"
        << " * until: " << timestamp_to_date_string(val_curr.at(1)) << " (" << val_curr.at(1) << ")\n";
    }
    h_curr->SetTitle(Form("%s (%lis window)", h_curr->GetTitle(), (val_curr.at(1)-val_curr.at(0)) / 1000));

    TCanvas* c = new TCanvas("c", "", 900, 700);
    c->SetTopMargin(0.11);
    c->SetBottomMargin(0.1);
    c->SetLeftMargin(0.1);
    c->SetRightMargin(0.02);
    c->cd();
    h_curr->SetLineColor(kBlue);
    h_curr->SetLineStyle(1);
    h_curr->SetLineWidth(2);
    h_curr->GetYaxis()->SetRangeUser(0, y_max);
    h_curr->GetYaxis()->SetLabelSize(0.036);
    h_curr->GetXaxis()->SetLabelSize(0.036);
    h_curr->Draw();

    TH1F* h_next = NULL;
    std::vector<long> val_next;
    int n_rows = 1;
    if (i < arr->GetEntries()-1)
    {
      n_rows++;
      h_next = (TH1F*)arr->At(i+1);
      val_next = get_validity_from_name(h_next->GetName());
      c->cd();
      h_next->SetLineColor(kRed);
      h_next->SetLineStyle(2);
      h_next->SetLineWidth(2);
      h_next->Draw("same");
    }

    TLegend l(0.66, 0.99-n_rows*0.05, 0.98, 0.99);
    l.AddEntry(h_curr,Form("%s#minus%s", 
      timestamp_to_date_string(val_curr.at(0), true, false).data(), 
      timestamp_to_date_string(val_curr.at(1), true, false).data()),"L");
    if(h_next) {
      l.AddEntry(h_next,Form("%s#minus%s", 
        timestamp_to_date_string(val_next.at(0), true, false).data(), 
        timestamp_to_date_string(val_next.at(1), true, false).data()),"L");
    }
    l.SetTextSize(0.036);
    l.SetBorderSize(0);
    l.SetFillStyle(0);
    l.SetMargin(0.30);
    l.Draw();

    TLatex ltx;
    ltx.SetTextSize(0.04);
    ltx.SetTextAlign(12);
    ltx.SetNDC();
    ltx.DrawLatex(0.1,0.94,Form("%s", h_curr->GetTitle()));

    c->Print(Form("plots/%i_%s/%s_%li.pdf", run, pass.data(), hname.data(), val_curr.at(0)));
    delete c;
  }

  f->Close();

  std::cout << "\nDone\n";
  return;
}
  