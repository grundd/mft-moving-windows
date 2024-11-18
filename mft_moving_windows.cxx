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
// custom headers
#include "binning_rof.h"

o2::ccdb::CcdbApi api_ccdb;
o2::ccdb::CcdbApi api_qcdb;

int _run;
std::string _pass;
std::string _path;
std::string _hname;
std::string _title;
std::string _opt_hist;
std::string _opt_plot;
std::string _str_sor;
std::string _str_eor;
bool _plot_next;
bool _rewrite_root;
int _aggr_histo;
long increment = 60000; // 1 min = 60 * 1000 ms

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
TH* download_histo(long ts, std::vector<long>* val = NULL, bool verbose = false)
{
  std::string s_run = std::to_string(_run);
  // create metadata:
  std::map<std::string, std::string> metadata;
  metadata["RunNumber"] = s_run;
  if(_pass != "online") metadata["PassName"] = _pass;
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
  TH* h = api_qcdb.retrieveFromTFileAny<TH>(_path+_hname,metadata,ts);
  // retrieve its headers
  if(h) {
    std::map<std::string, std::string> headers = api_qcdb.retrieveHeaders(_path+_hname,metadata,ts);
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

TCanvas* plot_histos (TH1F* h_curr, std::vector<long>* val_curr, float y_max = 0, 
  TH1F* h_next = NULL, std::vector<long>* val_next = NULL)
{
  int n_rows = 1;
  if(h_next) n_rows++;
  if(y_max == 0) y_max = h_curr->GetMaximum() * 1.05;

  TCanvas* c = new TCanvas("c", "", 900, 700);
  c->SetTopMargin(0.11);
  c->SetBottomMargin(0.1);
  c->SetLeftMargin(0.1);
  c->SetRightMargin(0.02);
  c->cd();

  h_curr->SetLineColor(kBlue);
  h_curr->SetLineStyle(1);
  h_curr->SetLineWidth(2);

  if(_opt_hist.find("logx") != string::npos) c->SetLogx();
  if(_opt_hist.find("logy") != string::npos) {
    float y_min = h_curr->GetMinimum(0);
    h_curr->GetYaxis()->SetRangeUser(y_min, y_max);
    c->SetLogy(); 
  } else {
    h_curr->GetYaxis()->SetRangeUser(0, y_max);
  }

  h_curr->GetYaxis()->SetLabelSize(0.036);
  h_curr->GetXaxis()->SetLabelSize(0.036);
  h_curr->Draw(_opt_plot.data());

  if (h_next) {
    h_next->SetLineColor(kRed);
    h_next->SetLineStyle(2);
    h_next->SetLineWidth(2);
    h_next->Draw(Form("%s same", _opt_plot.data()));
  }

  TLegend* l = new TLegend(0.66, 0.99-n_rows*0.05, 0.98, 0.99);
  l->AddEntry(h_curr,Form("%s#minus%s", 
    timestamp_to_date_string(val_curr->at(0), true, false).data(), 
    timestamp_to_date_string(val_curr->at(1), true, false).data()),"L");
  if(h_next) {
    l->AddEntry(h_next,Form("%s#minus%s", 
      timestamp_to_date_string(val_next->at(0), true, false).data(), 
      timestamp_to_date_string(val_next->at(1), true, false).data()),"L");
  }
  l->SetTextSize(0.036);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->SetMargin(0.30);
  l->Draw();

  TLatex* ltx = new TLatex();
  ltx->SetTextSize(0.04);
  ltx->SetTextAlign(12);
  ltx->SetNDC();
  ltx->DrawLatex(0.1,0.94,Form("%s", h_curr->GetTitle()));

  return c;
}

void run_moving_windows (bool verbose = false)
{
  gSystem->Exec("mkdir -p root_files/");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  // connect to ccdb, qcdb, qcdbmc    
  api_ccdb.init("alice-ccdb.cern.ch");
  api_qcdb.init("ali-qcdb-gpn.cern.ch:8083");

  long ts_sor = get_timestamp(_run, _str_sor);
  long ts_eor = get_timestamp(_run, _str_eor);

  std::string fname = Form("root_files/%i_%s.root", _run, _pass.data());

  // download the objects
  bool fexists = !gSystem->AccessPathName(fname.data());
  if(fexists && !_rewrite_root) {
    std::cout << fname << " already downloaded -> skipping\n";
  } else {
    TFile* f = new TFile(fname.data(), "recreate");
    TObjArray* arr = new TObjArray();
    long current_ts = ts_sor;
    bool next_obj_exists = true;
    while (next_obj_exists)
    {
      std::vector<long> val;
      TH1F* h = download_histo<TH1F>(current_ts, &val);
      if (h) {
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
        if(_opt_hist.find("rebinROF") != string::npos) h_new = rebin_rof(h_new);
        arr->Add(h_new);
        current_ts = valid_until + 1;
        if (current_ts >= ts_eor) next_obj_exists = false;
      } else if (current_ts < ts_eor) {
        std::cout << "The object for " << timestamp_to_date_string(current_ts) << " not found!\n";
        current_ts += increment;
        std::cout << " -> increasing the timestamp to " << timestamp_to_date_string(current_ts) << "\n";
      } else {
        std::cout << "EOR timestamp exceeded and the object does not exist!\n";
        next_obj_exists = false;
      }
    }
    f->cd();
    arr->Write(_hname.data(), TObject::kSingleKey);
    f->Write("",TObject::kWriteDelete);
    f->Close();
  }

  // plot the objects:

  gSystem->Exec(Form("mkdir -p plots/%i_%s/", _run, _pass.data()));
  TFile* f = TFile::Open(fname.data());
  TObjArray *arr = (TObjArray*) f->Get(_hname.data());

  if (arr->IsEmpty()) {
    std::cout << "Empty object array, no histograms to plot!\n";
    return;
  } 

  TH1F* h_total;
  std::vector<long> val_total;
  TH1F* h_first = (TH1F*)arr->At(0);
  if(h_first) {
    int n_bins = h_first->GetNbinsX();
    float x_low = h_first->GetBinLowEdge(1);
    float x_upp = h_first->GetBinLowEdge(n_bins+1);
    h_total = new TH1F("total", _title.data(), n_bins, x_low, x_upp);
    if(_opt_hist.find("rebinROF") != string::npos) h_total = rebin_rof(h_total);
  }

  TH1F* h_aggr;
  std::vector<long> val_aggr;
  if(_aggr_histo > 0) h_aggr = (TH1F*)h_total->Clone();

  // get the maximum
  float y_max = 0;
  for (int i = 0; i < arr->GetEntries(); i++)
  {
    TH1F* h_curr = (TH1F*)arr->At(i);
    h_total->Add(h_curr);
    float curr_max = h_curr->GetMaximum() * 1.05;
    if(curr_max > y_max) y_max = curr_max;
  }
  h_total->Scale(1./arr->GetEntries()); // an approximation!

  int n_aggr = 0;
  for (int i = 0; i < arr->GetEntries(); i++)
  {
    TH1F* h_curr = (TH1F*)arr->At(i);

    std::vector<long> val_curr = get_validity_from_name(h_curr->GetName());
    if(i == 0) val_total.push_back(val_curr.at(0));
    if(i == arr->GetEntries()-1) val_total.push_back(val_curr.at(0));
    if(verbose) {
      std::cout << "Histogram validity: \n"
        << " * from: " << timestamp_to_date_string(val_curr.at(0)) << " (" << val_curr.at(0) << ")\n"
        << " * until: " << timestamp_to_date_string(val_curr.at(1)) << " (" << val_curr.at(1) << ")\n";
    }
    std::string htitle = h_curr->GetTitle();
    if(!_title.empty()) htitle = _title;
    h_curr->SetTitle(Form("%s (%lis window)", htitle.data(), (val_curr.at(1)-val_curr.at(0)) / 1000));

    TH1F* h_next = NULL;
    std::vector<long> val_next;
    if (i < arr->GetEntries()-1 && _plot_next) {
      h_next = (TH1F*)arr->At(i+1);
      val_next = get_validity_from_name(h_next->GetName());
    }

    TCanvas* c = plot_histos(h_curr, &val_curr, y_max, h_next, &val_next);
    c->Print(Form("plots/%i_%s/%s_%li.pdf", _run, _pass.data(), _hname.data(), val_curr.at(0)));
    delete c;

    if (_aggr_histo > 0) 
    {
      std::vector<long> val_new = get_validity_from_name(h_curr->GetName());
      if(h_aggr->GetEntries() == 0) {
        val_aggr.push_back(val_new.at(0));
        val_aggr.push_back(val_new.at(1));
      } else {
        val_aggr.at(1) = val_new.at(1);
      }
      h_aggr->Add(h_curr);


      if (((i+1) % _aggr_histo == 0) || i == arr->GetEntries()-1) 
      {
        h_aggr->Scale(1./_aggr_histo);

        n_aggr++;
        TCanvas* c_aggr = plot_histos(h_aggr, &val_aggr, y_max);
        c_aggr->Print(Form("plots/%i_%s/%s_%02i.pdf", _run, _pass.data(), _hname.data(), n_aggr));
        delete c_aggr;

        h_aggr->Reset();
        val_aggr.clear();
      }
    }
  }

  TCanvas* c_tot = plot_histos(h_total, &val_total);
  c_tot->Print(Form("plots/%i_%s/%s_total.pdf", _run, _pass.data(), _hname.data()));
  delete c_tot;

  f->Close();

  std::cout << "\nDone\n";
  return;
}

void mft_moving_windows (
  int run, 
  std::string pass, 
  std::string path, 
  std::string hname, 
  std::string title, 
  std::string opt_hist,
  std::string opt_plot,
  std::string str_sor,
  std::string str_eor,
  bool plot_next,
  bool rewrite_root,
  int aggr_histo
) {
  _run = run;
  _pass = pass;
  _path = path;
  _hname = hname;
  _title = title;
  _opt_hist = opt_hist;
  _opt_plot = opt_plot;
  _str_sor = str_sor;
  _str_eor = str_eor;
  _plot_next = plot_next;
  _rewrite_root = rewrite_root;
  _aggr_histo = aggr_histo;
  run_moving_windows();
}
  