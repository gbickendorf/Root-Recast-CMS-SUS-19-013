#include "Plots.h"

void Plots::PlotPTMiss(vector<PassedEvent> events)
{
    TCanvas *c1 = new TCanvas("cPTMiss", "c1");
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
    histMetZ->SetFillColor(kGreen);
    histMetZ->GetYaxis()->SetRangeUser(0.1, 1e4);
    histMetZ->GetXaxis()->SetTitle("p^{miss}_{T} [GeV]");
    histMetZ->Draw("PLC");
    c1->SaveAs("Plots/METZ.eps");

    histMetW->SetFillColor(kBlue);
    histMetW->GetYaxis()->SetRangeUser(0.1, 1e4);
    histMetW->Draw();
    c1->SaveAs("Plots/METW.eps");

    histMettt->SetFillColor(kCyan);
    histMettt->GetYaxis()->SetRangeUser(0.1, 1e4);
    histMettt->Draw();
    c1->SaveAs("Plots/METtt.eps");

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

    //hs->GetXaxis()->SetTitle();
    hs->Draw("");
    hs->GetXaxis()->SetTitle("dASd");
    hs->Draw("same");
    leg.Draw("");
    gPad->Modified();
    gPad->Update();
    hs->GetYaxis()->SetTitle("Events/bin");
    hs->GetYaxis()->SetTitleSize(0.049F);
    hs->GetXaxis()->SetTitle("p^{miss}_{T} [GeV]");
    hs->GetXaxis()->SetTitleSize(0.039F);
    c1->SaveAs("Plots/Stack.eps");
    TFile *file = new TFile("Analysis.root", "UPDATE");
    histMetW->Write("histMetW", TObject::kOverwrite);
    histMetZ->Write("histMetZ", TObject::kOverwrite);
    hs->Write("hs", TObject::kOverwrite);
    file->Close();
    delete file;
}

void Plots::PlotMJ1(vector<PassedEvent> events, vector<double> params)
{
    TF1 *linFunc = (TF1 *)gROOT->GetFunction("pol1");
    linFunc->SetParameters(params[0], params[1]);
    TF1 *linFuncLeft = (TF1 *)linFunc->Clone();
    linFuncLeft->SetRange(40, 70);
    TF1 *linFuncMiddle = (TF1 *)linFunc->Clone();
    linFuncMiddle->SetRange(70, 100);
    linFuncMiddle->SetLineStyle(kDashed);
    TF1 *linFuncRight = (TF1 *)linFunc->Clone();
    linFuncRight->SetRange(100, 140);

    TCanvas *c1 = new TCanvas("cPlot2", "c1");
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
    leadMjZ->SetFillColor(kGreen);
    leadMjW->SetFillColor(kBlue);
    leadMjtt->SetFillColor(kCyan);
    leadMjZ->Draw();
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
    hs->GetYaxis()->SetTitleSize(0.049F);
    hs->GetXaxis()->SetTitle("Leading jet m_{jet} [GeV]");
    hs->GetXaxis()->SetTitleSize(0.039F);

    hs->Draw();
    //linFuncLeft->Draw("same");
    //linFunc->SetRange(100,140);
    //linFuncRight->Draw("same");
    c1->SaveAs("Plots/MJ1.eps");
}

void Plots::PlotSignalRegion(vector<PassedEvent> events, vector<double> &bi)
{
    const double bins[7] = {300.0, 450.0, 600.0, 800.0, 1000.0, 1200.0, 2000.0};
    TH1 *histSR = new TH1F("SignalRegion", "SR", 6, bins);
    TH1 *histBG = new TH1F("BackGround", "BG", 6, bins);
    for (int i = 0; i < 6; i++)
    {
        histBG->SetBinContent(i + 1, bi[i]);
        /* code */
    }
    TCanvas *c1 = new TCanvas("cSRPlot", "c1");
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
    histSR->Draw("E");
    histBG->Draw("HIST same");
    c1->SaveAs("Plots/SR.eps");
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
    c1->SaveAs("Plots/zShape.eps");
}

void Plots::PlotPhotonLeptonValidation(vector<PassedEvent> events)
{
    TCanvas *c1 = new TCanvas("cLepPhoVal", "c1");
    c1->cd();
    gPad->SetLogy();
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
    cout << "Photon SR Val: " <<histPhotonValidationSR->Integral() << endl;
    cout << "Photon CR Val: " <<histPhotonValidationCR->Integral() << endl;
    cout << "Lepton SR Val: " <<histLeptonValidationSR->Integral() << endl;
    cout << "Lepton CR Val: " <<histLeptonValidationCR->Integral() << endl;
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
    c1->SaveAs("Plots/PhotonVal.eps");

    THStack *hsLep = new THStack("LeptonVal", "Simulation 137 fb^{-1} (13 TeV)");
    hsLep->Add(histLeptonValidationSR);
    hsLep->Add(histLeptonValidationCR);
    hsLep->Draw("nostack");
    c1->SaveAs("Plots/LeptonVal.eps");
}

void Plots::PlotAll()
{
    vector<PassedEvent> events;
    double transferfactor, bNorm;
    vector<double> B_i, params;

    RootIO::ReadEvents("Normtest.root",events);
    Plots::PlotPTMiss(events);
    Analysis::FitCherbyshev(events, bNorm, params);
    Analysis::CalcTransferFactor(events, bNorm, B_i, transferfactor);
    Plots::PlotMJ1(events, params);
    Plots::PlotSignalRegion(events, B_i);
    Plots::PlotPTMiss(events);
    Plots::PlotPhotonLeptonValidation(events);
}