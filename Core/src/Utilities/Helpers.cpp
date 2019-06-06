// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Utilities/Helpers.hpp"

namespace FW {

namespace PlotHelpers {
  TH1F*
  bookHisto(const char*    histName,
            const char*    histTitle,
            const Binning& varBinning)
  {
    TH1F* hist = new TH1F(
        histName, histTitle, varBinning.nBins, varBinning.min, varBinning.max);
    hist->GetXaxis()->SetTitle(varBinning.title.c_str());
    hist->GetYaxis()->SetTitle("Entries");
    hist->Sumw2();
    return hist;
  }

  TH2F*
  bookHisto(const char*    histName,
            const char*    histTitle,
            const Binning& varXBinning,
            const Binning& varYBinning)
  {
    TH2F* hist = new TH2F(histName,
                          histTitle,
                          varXBinning.nBins,
                          varXBinning.min,
                          varXBinning.max,
                          varYBinning.nBins,
                          varYBinning.min,
                          varYBinning.max);
    hist->GetXaxis()->SetTitle(varXBinning.title.c_str());
    hist->GetYaxis()->SetTitle(varYBinning.title.c_str());
    hist->Sumw2();
    return hist;
  }

  void
  fillHisto(TH1F* hist, float value, float weight)
  {
    if (hist == nullptr) {
      std::bad_alloc();
    } else {
      hist->Fill(value, weight);
    }
  }

  void
  fillHisto(TH2F* hist, float xValue, float yValue, float weight)
  {
    if (hist == nullptr) {
      std::bad_alloc();
    } else {
      hist->Fill(xValue, yValue, weight);
    }
  }

  void
  anaHisto(TH1D* inputHist, int j, TH1F* meanHist, TH1F* widthHist)
  {
    // evaluate mean and width via the Gauss fit
    if (inputHist == nullptr) {
      std::bad_alloc();
    } else if (inputHist->GetEntries() > 0) {
      TFitResultPtr r = inputHist->Fit("gaus", "QS0");
      if (r.Get() and ((r->Status() % 1000) == 0)) {
        // fill the mean and width into 'j'th bin of the meanHist and widthHist,
        // respectively
        meanHist->SetBinContent(j, r->Parameter(1));
        meanHist->SetBinError(j, r->ParError(1));
        widthHist->SetBinContent(j, r->Parameter(2));
        widthHist->SetBinError(j, r->ParError(2));
      }
    }
  }

  TEfficiency*
  bookEff(const char* effName, const char* effTitle, const Binning& varBinning)
  {
    TEfficiency* efficiency = new TEfficiency(
        effName, effTitle, varBinning.nBins, varBinning.min, varBinning.max);
    return efficiency;
  }

  TEfficiency*
  bookEff(const char*    effName,
          const char*    effTitle,
          const Binning& varXBinning,
          const Binning& varYBinning)
  {
    TEfficiency* efficiency = new TEfficiency(effName,
                                              effTitle,
                                              varXBinning.nBins,
                                              varXBinning.min,
                                              varXBinning.max,
                                              varYBinning.nBins,
                                              varYBinning.min,
                                              varYBinning.max);
    return efficiency;
  }

  void
  fillEff(TEfficiency* efficiency, float value, bool status)
  {
    if (efficiency == nullptr) {
      std::bad_alloc();
    } else {
      efficiency->Fill(status, value);
    }
  }

  void
  fillEff(TEfficiency* efficiency, float valueX, float valueY, bool status)
  {
    if (efficiency == nullptr) {
      std::bad_alloc();
    } else {
      efficiency->Fill(status, valueX, valueY);
    }
  }
}  // namespace PlotHelpers

}  // namespace FW
