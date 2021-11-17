#include "Plots.h"

void Plots::Plot1(vector<PassedEvent> events)
{
    TCanvas *c1 = new TCanvas("c1", "c1");
    c1->cd();
    const double bins[7] = {300.0, 450.0, 600.0, 800.0, 1000.0, 1200.0, 2000.0};
    TH1 *histMetZ = new TH1F("ptmiss_Z", "ptmiss Z", 6, bins);
    TH1 *histMetW = new TH1F("ptmiss_W", "ptmiss W", 6, bins);
    for (size_t i = 0; i < events.size(); i++)
    {
        if (events[i].status == 0)
            histMetZ->Fill(events[i].ptmiss);
        if (events[i].status == 1)
            histMetW->Fill(events[i].ptmiss);
    }
    gPad->SetLogy();
    histMetZ->SetFillColor(kGreen);
    histMetZ->GetYaxis()->SetRangeUser(0.1, 1e4);
    histMetZ->GetXaxis()->SetTitle("p^{miss}_{T} [GeV]");
    //ax1->Lat
    histMetZ->Draw("PLC");
    c1->SaveAs("Plots/METZ.eps");
    histMetW->SetFillColor(kBlue);
    histMetW->GetYaxis()->SetRangeUser(0.1, 1e4);
    histMetW->Draw();
    c1->SaveAs("Plots/METW.eps");
    THStack *hs = new THStack("hs", "Simulation 137 fb^{-1} (13 TeV)");
    histMetW->SetStats(0);
    histMetZ->SetStats(0);
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
    gPad->Modified(); gPad->Update();
    hs->GetYaxis()->SetTitle("Events/bin");
    hs->GetYaxis()->SetTitleSize(0.049F);
    hs->GetXaxis()->SetTitle("p^{miss}_{T} [GeV]");
    hs->GetXaxis()->SetTitleSize(0.039F);
    c1->SaveAs("Plots/Stack.eps");
    TFile *file = new TFile("Analysis.root", "UPDATE");
    histMetW->Write("histMetW",TObject::kOverwrite);
    histMetZ->Write("histMetZ",TObject::kOverwrite);
    hs->Write("hs",TObject::kOverwrite);
    file->Close();
    delete file;
}

void Plots::Plot2(vector<PassedEvent> events)
{
    TCanvas *c1 = new TCanvas("c1", "c1");
    c1->cd();
    TH1F *leadMjZ = new TH1F("mj1Z","",20,40.0,140.0);
    TH1F *leadMjW = new TH1F("mj1W","",20,40.0,140.0);
    TH1F *leadMjtt = new TH1F("mj1tt","",20,40.0,140.0);
    for (size_t i = 0; i < events.size(); i++)
    {
        if(events[i].status==0)
            leadMjZ->Fill(events[i].mj1);
        if(events[i].status==1)
            leadMjW->Fill(events[i].mj1);
        if(events[i].status==2)
            leadMjtt->Fill(events[i].mj1);
    }
    leadMjZ->SetFillColor(kGreen);
    leadMjW->SetFillColor(kBlue);
    leadMjtt->SetFillColor(kCyan);
    leadMjZ->Draw();
    THStack *hs = new THStack("hsLeadingJetM", "Simulation 137 fb^{-1} (13 TeV)");
    hs->Add(leadMjtt);
    hs->Add(leadMjW);
    hs->Add(leadMjZ);
    
    hs->Draw();
    hs->GetYaxis()->SetTitle("Events/(5 GeV)");
    hs->GetYaxis()->SetTitleSize(0.049F);
    hs->GetXaxis()->SetTitle("Leading jet m_{jet} [GeV]");
    hs->GetXaxis()->SetTitleSize(0.039F);

    hs->Draw();
    c1->SaveAs("Plots/MJ1.eps");
}