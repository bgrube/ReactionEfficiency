// see https://root-forum.cern.ch/t/multiple-graphs-in-tgraph2d/20717/2

void
testGraph2D()
{
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(1111);

   TCanvas *c = new TCanvas("c","Graph2D example",0,0,800,800);

   TH3F *frame3d = new TH3F("frame3d","frame3d",10,0,10,10,1,10,10,0,10);
   frame3d->Draw();

   Double_t x1[4] = {1,1,2,2};
   Double_t y1[4] = {1,2,1,2};
   Double_t z1[4] = {1,2,3,4};
   TGraph2D *g2d1 = new TGraph2D(4,x1,y1,z1);
   g2d1->SetMarkerStyle(20);
   g2d1->Draw("PCOL SAME");

   Double_t x2[4] = {5,5,7,7};
   Double_t y2[4] = {5,7,5,7};
   Double_t z2[4] = {5,7,3,8};
   TGraph2D *g2d2 = new TGraph2D(4,x2,y2,z2);
   g2d2->SetMarkerStyle(21);

   g2d2->Draw("PCOL SAME");
}
