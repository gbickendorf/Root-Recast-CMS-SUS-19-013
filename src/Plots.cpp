#include "Plots.h"

void Plots::PlotPTMiss(vector<PassedEvent> events)
{
    TCanvas *c1 = new TCanvas("cPTMiss", "c1",800,700);
    gPad->SetTickx();
    gPad->SetTicky();
    c1->SetRightMargin(0.04);
    c1->SetLeftMargin(0.13);
    c1->SetBottomMargin(0.13);
    c1->cd();
    const double bins[7] = {300.0, 450.0, 600.0, 800.0, 1000.0, 1200.0, 2000.0};
    TH1 *histMetZ = new TH1F("ptmiss_Z", "ptmiss Z", 6, bins);
    TH1 *histMetW = new TH1F("ptmiss_W", "ptmiss W", 6, bins);
    TH1 *histMettt = new TH1F("ptmiss_tt", "ptmiss tt", 6, bins);
    for (size_t i = 0; i < events.size(); i++)
    {
        if (events[i].status == 0)
            histMetZ->Fill(events[i].ptmiss);
        if (events[i].status == 1)
            histMetW->Fill(events[i].ptmiss);
        if (events[i].status == 2)
            histMettt->Fill(events[i].ptmiss);
    }
    gPad->SetLogy();
    histMetZ->SetFillColor(kGreen+1);
    
    histMetZ->GetXaxis()->SetTitle("p^{miss}_{T} [GeV]");
    histMetZ->Draw("PLC");
    //c1->SaveAs("Plots/METZ.pdf");

    histMetW->SetFillColor(kBlue+1);
    histMetW->Draw();
    //c1->SaveAs("Plots/METW.pdf");

    histMettt->SetFillColor(kCyan+1);
    histMettt->Draw();
    //c1->SaveAs("Plots/METtt.pdf");

    THStack *hs = new THStack("hs", "Simulation 137 fb^{-1} (13 TeV)");
    histMettt->SetStats(0);
    histMetW->SetStats(0);
    histMetZ->SetStats(0);
    hs->Add(histMettt);
    hs->Add(histMetW);
    hs->Add(histMetZ);
    TLegend leg(.3, .8, .9, .9, "");
    leg.SetFillColor(0);
    leg.SetNColumns(2);
    leg.AddEntry(histMetZ, "Z + jets");
    leg.AddEntry(histMetW, "W + jets");
    leg.AddEntry(histMettt, "tt + jets");

    //hs->GetXaxis()->SetTitle();
    hs->Draw("");
    leg.Draw("");
    gPad->Modified();
    gPad->Update();
    hs->GetYaxis()->SetTitle("Events/bin");
    hs->GetYaxis()->SetTitleSize(0.06F);
    hs->GetXaxis()->SetTitle("p^{miss}_{T} [GeV]");
    hs->GetXaxis()->SetTitleSize(0.06F);

    hs->GetYaxis()->SetLimits(1, 6000); // "xmin", "xmax"
    hs->GetYaxis()->SetRangeUser(1,3000);
    hs->SetMinimum(2.5e0); // "ymin"
    hs->SetMaximum(3.8e3);
    c1->SaveAs("Plots/Stack.pdf");

    TFile *file = new TFile("Analysis.root", "UPDATE");
    histMetW->Write("histMetW", TObject::kOverwrite);
    histMetZ->Write("histMetZ", TObject::kOverwrite);
    histMettt->Write("histMettt", TObject::kOverwrite);
    hs->Write("hsMet", TObject::kOverwrite);
    file->Close();
    delete file;
}

void Plots::PlotMJ1(vector<PassedEvent> events, vector<double> params)
{
    TF1 *linFunc = (TF1 *)gROOT->GetFunction("pol1");
    linFunc->SetParameters(params[0], params[1]);
    linFunc->SetParError(0, params[2]);
    linFunc->SetParError(0, params[3]);
    TF1 *linFuncLeft = (TF1 *)linFunc->Clone();
    linFuncLeft->SetRange(40, 70);
    TF1 *linFuncMiddle = (TF1 *)linFunc->Clone();
    linFuncMiddle->SetRange(70, 100);
    linFuncMiddle->SetLineStyle(kDashed);
    TF1 *linFuncRight = (TF1 *)linFunc->Clone();
    linFuncRight->SetRange(100, 140);
    linFuncMiddle->SetLineColor(kBlue);
    linFuncLeft->SetLineColor(kBlue);
    linFuncRight->SetLineColor(kBlue);

    TCanvas *c1 = new TCanvas("cPlot2", "c1",800,700);
    gPad->SetTickx();
    gPad->SetTicky();
    c1->SetRightMargin(0.04);
    c1->SetLeftMargin(0.13);
    c1->SetBottomMargin(0.13);
    c1->cd();
    TH1F *leadMjZ = new TH1F("mj1Z", "", 20, 40.0, 140.0);
    TH1F *leadMjW = new TH1F("mj1W", "", 20, 40.0, 140.0);
    TH1F *leadMjtt = new TH1F("mj1tt", "", 20, 40.0, 140.0);
    TH2F *mj = new TH2F("mjs", "", 20, 40.0, 140.0, 20, 40.0, 140.0);
    for (size_t i = 0; i < events.size(); i++)
    {

        if (events[i].mj2 < 70.0 || events[i].mj2 > 100.0)
            continue;
        mj->Fill(events[i].mj1, events[i].mj2);
        if (events[i].status == 0)
            leadMjZ->Fill(events[i].mj1);
        if (events[i].status == 1)
            leadMjW->Fill(events[i].mj1);
        if (events[i].status == 2)
            leadMjtt->Fill(events[i].mj1);
    }
    leadMjZ->SetFillColor(kGreen+1);
    leadMjW->SetFillColor(kBlue+1);
    leadMjtt->SetFillColor(kCyan+1);
    leadMjZ->Draw();
    leadMjZ->GetYaxis()->SetRangeUser(0, 150);
    leadMjZ->GetListOfFunctions()->Add(linFuncLeft);
    leadMjZ->GetListOfFunctions()->Add(linFuncMiddle);
    leadMjZ->GetListOfFunctions()->Add(linFuncRight);
    THStack *hs = new THStack("hsLeadingJetM", "Simulation 137 fb^{-1} (13 TeV)");
    hs->Add(leadMjtt);
    hs->Add(leadMjW);
    //leadMjtt->Add(linFunc);
    hs->Add(leadMjZ);
    hs->Draw();
    hs->GetYaxis()->SetTitle("Events/(5 GeV)");
    hs->GetYaxis()->SetTitleSize(0.06F);
    hs->GetXaxis()->SetTitle("Leading jet m_{jet} [GeV]");
    hs->GetXaxis()->SetTitleSize(0.06F);
    hs->Draw();
    hs->GetYaxis()->SetRangeUser(0, 150);

    TLegend leg(.3, .7, .9, .9, "Subleading jet m_{jet} in Z signal window");
    leg.SetFillColor(0);
    leg.SetNColumns(2);
    leg.AddEntry(leadMjZ, "Z + jets");
    leg.AddEntry(leadMjW, "W + jets");
    leg.AddEntry(leadMjtt, "tt + jets");
    leg.AddEntry(linFuncLeft, "Linear fit with unc.");
    leg.Draw();
    linFuncLeft->SetFillColor(kBlue);
    linFuncLeft->SetFillStyle(3345);
    linFuncLeft->Draw("SAME E2");
    linFuncRight->SetFillColor(kBlue);
    linFuncRight->SetFillStyle(3345);
    linFuncRight->Draw("SAME E2");
    hs->SetMaximum(150);
    c1->SaveAs("Plots/MJ1.pdf");
}

void Plots::PlotSignalRegion(vector<PassedEvent> events, vector<double> NCRi, double transferFactor,double transferFactorErr)
{
    const double bins[7] = {300.0, 450.0, 600.0, 800.0, 1000.0, 1200.0, 2000.0};
    TH1 *histSR = new TH1F("SignalRegion", "", 6, bins);
    TH1 *histBG = new TH1F("BackGround", "", 6, bins);
    for (int i = 0; i < 6; i++)
    {
        histBG->SetBinContent(i + 1, NCRi[i] * transferFactor);
        histBG->SetBinError(i+1,sqrt(pow(NCRi[i] * transferFactorErr,2)+pow(sqrt(NCRi[i]) * transferFactor,2)));
        /* code */
    }
    TCanvas *c1 = new TCanvas("cSRPlot", "c1",800,600);
    c1->cd();
    gPad->SetLogy();
    for (size_t i = 0; i < events.size(); i++)
    {
        PassedEvent event = events[i];
        if (event.status >= 10)
            continue;
        if (event.mj1 > 70 && event.mj1 < 100 && event.mj2 > 70 && event.mj2 < 100)
            histSR->Fill(event.ptmiss);

        /* code */
    }

    
    histBG->SetLineColor(kRed);
    histBG->GetXaxis()->SetTitle("p^{miss}_{T} [GeV]");
    histBG->GetYaxis()->SetTitle("Events / bin");
    TH1 * histbgCopy = (TH1F*) histBG->Clone();
    histbgCopy->SetStats(0);
    histbgCopy->Draw("HIST");
    histBG->SetStats(0);
    //histBG->SetFillColor(kRed);
    histBG->SetFillStyle(3345);
    histBG->SetFillColorAlpha(kRed-9,0.3);
    histBG->Draw("same E2"); 
    histBG->SetStats(0);

    histSR->Draw("same E");
    histSR->GetYaxis()->SetRangeUser(0.01, 1000.0);
    histSR->SetStats(0);
    histSR->SetFillStyle(0);
    histSR->SetMarkerStyle(kFullCircle);
    histSR->SetLineColor(kBlack);

    TLegend leg(.3, .8, .9, .9, "");
    leg.SetFillColor(0);
    leg.SetNColumns(2);
    leg.AddEntry(histBG, "Pred. (stat. #oplus syst.) unc.");
    leg.AddEntry(histSR, "Observed data");
    leg.Draw();
    

    //histBG->Draw("same hist"); 
    //histBG->Draw("e2SAME");
    c1->SaveAs("Plots/SR.pdf");
}

void Plots::PlotPTShape(vector<PassedEvent> events)
{
    TCanvas *c1 = new TCanvas("cPTShape", "c1");
    c1->cd();
    gPad->SetLogy();
    const double bins[12] = {300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 1100.0, 1200.0, 1400.0, 2000.0};
    TH1 *histSRW = new TH1F("ShapeSignalRegionW", "SR", 11, bins);
    TH1 *histCRW = new TH1F("ShapeControlregionW", "CR", 11, bins);
    TH1 *histSRZ = new TH1F("ShapeSignalRegionZ", "SR", 11, bins);
    TH1 *histCRZ = new TH1F("ShapeControlregionZ", "CR", 11, bins);
    for (size_t i = 0; i < events.size(); i++)
    {
        PassedEvent evnt = events[i];
        if (evnt.mj1 > 70 && evnt.mj1 < 100 && evnt.mj2 > 70 && evnt.mj2 < 100)
        {
            if (evnt.status == 0)
                histSRZ->Fill(evnt.ptmiss);
            if (evnt.status == 1)
                histSRW->Fill(evnt.ptmiss);
        }
        if ((evnt.mj1 < 70 || evnt.mj1 > 100) && (evnt.mj2 < 70 || evnt.mj2 > 100))
        {
            if (evnt.status == 0)
                histCRZ->Fill(evnt.ptmiss);
            if (evnt.status == 1)
                histCRW->Fill(evnt.ptmiss);
        }
    }
    histCRZ->Scale(1.0 / histCRZ->Integral());
    histSRZ->Scale(1.0 / histSRZ->Integral());
    histCRW->Scale(1.0 / histCRW->Integral());
    histSRW->Scale(1.0 / histSRW->Integral());
    THStack *hs = new THStack("hsZShape", "Simulation 137 fb^{-1} (13 TeV)");
    hs->Add(histSRZ);
    hs->Add(histCRZ);
    hs->Draw("nostack");
    c1->SaveAs("Plots/zShape.pdf");
}

void Plots::PlotPhotonLeptonValidation(vector<PassedEvent> events)
{
    TCanvas *c1 = new TCanvas("cLepPhoVal", "c1");
    c1->cd();
    const double bins[7] = {200.0, 300.0, 450.0, 600.0, 800.0, 1000.0, 1400.0};
    TH1 *histPhotonValidationSR = new TH1F("PhotonValidationSR", "Photon SR", 6, bins);
    TH1 *histLeptonValidationSR = new TH1F("LeptonValidationSR", "Lepton SR", 6, bins);
    TH1 *histPhotonValidationCR = new TH1F("PhotonValidationCR", "Photon CR", 6, bins);
    TH1 *histLeptonValidationCR = new TH1F("LeptonValidationCR", "Lepton CR", 6, bins);
    double prev = 0;
    for (size_t i = 0; i < events.size(); i++)
    {
        PassedEvent evnt = events[i];
        if (prev == evnt.ptmiss)
            continue;
        prev = evnt.ptmiss;
        if (evnt.mj1 > 70 && evnt.mj1 < 100 && evnt.mj2 > 70 && evnt.mj2 < 100)
        {
            if (evnt.status == 11)
                histPhotonValidationSR->Fill(evnt.ptmiss);
            if (evnt.status == 10)
                histLeptonValidationSR->Fill(evnt.ptmiss);
        }
        if ((evnt.mj1 < 70 || evnt.mj1 > 100) && (evnt.mj2 < 70 || evnt.mj2 > 100))
        {
            if (evnt.status == 11)
                histPhotonValidationCR->Fill(evnt.ptmiss);
            if (evnt.status == 10)
                histLeptonValidationCR->Fill(evnt.ptmiss);
        }
    }
    cout << "Photon SR Val: " << histPhotonValidationSR->Integral() << endl;
    cout << "Photon CR Val: " << histPhotonValidationCR->Integral() << endl;
    cout << "Lepton SR Val: " << histLeptonValidationSR->Integral() << endl;
    cout << "Lepton CR Val: " << histLeptonValidationCR->Integral() << endl;
    double NormLepCR = histLeptonValidationCR->Integral();
    double NormLepSR = histLeptonValidationSR->Integral();
    double NormPhoCR = histPhotonValidationCR->Integral();
    double NormPhoSR = histPhotonValidationSR->Integral();

    TH1F *RatioPhoton = new TH1F("RatioPhoton", "RatioPhoton", 6, bins);
    TH1F *RatioLepton = new TH1F("RatioLepton", "RatioLepton", 6, bins);
    TH1F *RatioPhotonSRErr = new TH1F("RatioPhotonSRErr", "RatioPhotonSRErr", 6, bins);
    TH1F *RatioLeptonSRErr = new TH1F("RatioLeptonSRErr", "RatioLeptonSRErr", 6, bins);
    for (size_t i = 0; i < 7; i++)
    {
        if (histPhotonValidationCR->GetBinContent(i) != 0.0)
        {
            RatioPhoton->SetBinContent(i, histPhotonValidationSR->GetBinContent(i) / histPhotonValidationCR->GetBinContent(i));
            RatioPhoton->SetBinError(i, RatioPhoton->GetBinContent(i) / sqrt(histPhotonValidationCR->GetBinContent(i)));
        }
        if (histLeptonValidationCR->GetBinContent(i) != 0.0)
        {
            RatioLepton->SetBinContent(i, histLeptonValidationSR->GetBinContent(i) / histLeptonValidationCR->GetBinContent(i));
            RatioLepton->SetBinError(i, RatioLepton->GetBinContent(i) / sqrt(histLeptonValidationCR->GetBinContent(i)));
        }
    }
    RatioPhoton->Fit("pol0");
    RatioPhoton->Draw("E");

    for (size_t i = 0; i < 7; i++)
    {
        if (histPhotonValidationSR->GetBinContent(i) != 0.0)
        {
            RatioPhotonSRErr->SetBinContent(i, RatioPhoton->GetFunction("pol0")->GetParameter(0));
            cout << RatioPhoton->GetBinContent(i) / sqrt(histPhotonValidationSR->GetBinContent(i)) << endl;
            RatioPhotonSRErr->SetBinError(i, RatioPhoton->GetBinContent(i) / sqrt(histPhotonValidationSR->GetBinContent(i)));
        }
    }

    RatioPhotonSRErr->SetMarkerSize(0.);
    RatioPhotonSRErr->SetFillStyle(3003);
    RatioPhotonSRErr->SetMarkerColor(kBlack);
    RatioPhotonSRErr->SetFillColor(kBlack);
    RatioPhotonSRErr->Draw("SAME E2");
    RatioPhoton->GetYaxis()->SetRangeUser(0.0, 0.6);
    c1->SaveAs("Plots/RatioPhoton.pdf");

    RatioLepton->Fit("pol0");
    RatioLepton->Draw("E");

    for (size_t i = 0; i < 7; i++)
    {
        if (histLeptonValidationSR->GetBinContent(i) != 0.0)
        {
            RatioLeptonSRErr->SetBinContent(i, RatioLepton->GetFunction("pol0")->GetParameter(0));
            RatioLeptonSRErr->SetBinError(i, RatioLepton->GetBinContent(i) / sqrt(histLeptonValidationSR->GetBinContent(i)));
        }
    }
    RatioLeptonSRErr->SetMarkerSize(0.);
    RatioLeptonSRErr->SetFillStyle(3003);
    RatioLeptonSRErr->SetMarkerColor(kBlack);
    RatioLeptonSRErr->SetFillColor(kBlack);
    RatioLeptonSRErr->Draw("SAME E2");
    RatioLepton->GetYaxis()->SetRangeUser(0.0, 0.8);
    c1->SaveAs("Plots/RatioLepton.pdf");
    gPad->SetLogy();
    histPhotonValidationSR->Scale(1.0 / histPhotonValidationSR->Integral());
    histLeptonValidationSR->Scale(1.0 / histLeptonValidationSR->Integral());
    histPhotonValidationCR->Scale(1.0 / histPhotonValidationCR->Integral());
    histLeptonValidationCR->Scale(1.0 / histLeptonValidationCR->Integral());
    histPhotonValidationCR->SetLineColor(kBlue);
    histPhotonValidationSR->SetLineColor(kBlack);
    histLeptonValidationCR->SetLineColor(kBlue);
    histLeptonValidationSR->SetLineColor(kBlack);

    histLeptonValidationSR->SetMarkerStyle(kFullCircle);
    histPhotonValidationSR->SetMarkerStyle(kFullCircle);
    histLeptonValidationCR->SetMarkerStyle(kFullSquare);
    histPhotonValidationCR->SetMarkerStyle(kFullSquare);
    histLeptonValidationCR->SetMarkerColor(kBlue);
    histPhotonValidationCR->SetMarkerColor(kBlue);

    THStack *hsPho = new THStack("PhotonVal", "Simulation 137 fb^{-1} (13 TeV)");
    hsPho->Add(histPhotonValidationSR);
    hsPho->Add(histPhotonValidationCR);
    hsPho->Draw("nostack");
    //auto rpPho = new TRatioPlot(histPhotonValidationSR, histPhotonValidationCR);
    //rpPho->Draw();
    c1->SaveAs("Plots/PhotonVal.pdf");

    THStack *hsLep = new THStack("LeptonVal", "Simulation 137 fb^{-1} (13 TeV)");
    hsLep->Add(histLeptonValidationSR);
    hsLep->Add(histLeptonValidationCR);
    hsLep->Draw("nostack");
    //auto rpLep = new TRatioPlot(histLeptonValidationSR, histLeptonValidationCR);
    //rpPho->Draw();
    c1->SaveAs("Plots/LeptonVal.pdf");
}

void Plots::PlotAll()
{
    vector<PassedEvent> events;
    double transferfactor,transferfactorErr, bNorm,bNormErr;
    vector<double> NCRi, params;

    RootIO::ReadEvents("Normtest.root", events);
    Plots::PlotPTMiss(events);
    Analysis::FitCherbyshev(events, bNorm,bNormErr, params);
    Analysis::CalcTransferFactor(events, bNorm,bNormErr, NCRi, transferfactor,transferfactorErr);
    
    Plots::PlotMJ1(events, params);
    Plots::PlotSignalRegion(events, NCRi, transferfactor,transferfactorErr);
    Plots::PlotPTMiss(events);
    Plots::PlotPhotonLeptonValidation(events);
}