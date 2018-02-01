#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <string>
#include <vector>
#include <map>
#include <iomanip>
#include <sstream>
#include <time.h>
#include "TF1.h"
#include "TF2.h"
#include "TH1.h"
#include "TH2.h"
#include "TAxis.h"
#include "TMath.h"
#include "TEllipse.h"
#include "TStyle.h"
#include "TFile.h"
#include "TColor.h"
#include "TSpectrum.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TMinuit.h"
#include "TView3D.h"
#include "TView.h" 


// This code does the following :
//  - creates a 2-D function of an annulus. 
//  - draws the function as a surf2 plot. 
// This example can be executed via the interpreter or ACLIC
//   root > .x FitAnnulus.C
//   root > .x FitAnnulus.C++


// Functions : 
Double_t g3(Double_t *x, Double_t *par); 
Double_t fun3(Double_t *x, Double_t *par); 
void PixToMMzoom(TH2F * h1, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax); 

void Show(int iRun); 
// Global switches : 
bool canvas1 = true;
bool canvas2 = true; 
bool fithisto = true;
bool writehisto = true;  
bool savecanvas = true; 

//Global variables : 
int RunCount=-1; 
vector<Int_t> runs;
vector<Double_t>  ampFit;
vector<Double_t>  ampFitout;
vector<Double_t>  a_Fit; 
vector<Double_t>  b_Fit;
vector<Double_t>  thetaFit;
vector<Double_t>  kbT_Fit; 
vector<Double_t>  r_outFit; 
vector<Double_t>  r_inFit; 
vector<Double_t>  x0Fit;
vector<Double_t>  y0Fit; 



// Body of the code : 
 void annulus_fit(char* list){
    RunCount=0; 
    int RunNow;
    ifstream flist(list);
    while (!flist.eof() ) { 
        flist>> RunNow; 
        runs.push_back(RunNow); 
            for(int i=0; i<5; i++){
                cout << "*************************" << "\t" <<RunNow << "*************************" << endl; 
            }
            Show(RunNow); 
            RunCount++;
    }
    if(fithisto){
        ofstream out_run_list_file; // here I go... 
        ostringstream s;
        s << "list.txt";
        out_run_list_file.open(s.str().c_str());
        // Writing list of runs
        if(out_run_list_file.is_open()) {
            for(UInt_t i=0; i<(UInt_t)runs.size(); i++) {
                if(runs[i]!=0) {
                    out_run_list_file 
                    << runs[i] << "\t" 
                    << ampFit[i] << "\t" 
                    << ampFitout[i] << "\t" 
                    << a_Fit[i] << "\t" 
                    << b_Fit[i] << "\t" 
                    << thetaFit[i] << "\t" 
                    << kbT_Fit[i] << "\t" 
                    << r_outFit[i] << "\t" 
                    << r_inFit[i] << "\t" 
                    << x0Fit[i] << "\t" 
                    << y0Fit[i] << endl;
                }
            }
        }
    }
    flist.close();   
}

// Function that performs and shows the fit : 
void Show(int iRun) {

    // Set style : 
    gStyle->SetPalette(1);

    // Parameters for fit. 
    Double_t pi = 3.14159265359; 
    const Int_t npar = 10;

    // Coordinates I want to zoom in on : 
    Double_t x0 = 2.8; Double_t x1 = 13.2; 
    Double_t y0 = 3.2; Double_t y1 = 13.6; 
    /*Double_t x0 = 3.5; Double_t x1 = 9.; 
    Double_t y0 = 4.5; Double_t y1 = 10.;*/

    // Setup Canvas : 
    Double_t w = 700;
    Double_t h = 700;
    TCanvas * can_histo = new TCanvas("elliptical annulus", "study of the elliptical annulus ", w, h);       // name, title, width, height. 
    TCanvas * can_Fit = new TCanvas("data & fit", "data & fit", 2*w, h);
    can_histo->Divide(2,2);                                    // // create 6 pads (2 divisions along x, 1 along y). 
    can_Fit  ->Divide(2,1); 

    // Get Image : 
    char nameFile[30];
    sprintf(nameFile, "../root_files/raw_data/%d_annulus.root", iRun);
    TFile* file = new TFile(nameFile, "READ"); 
    TH2F * histo = (TH2F*)file->Get("ImgFull5THDetCMOS_0"); 

    // Convert pixels to mm in the original image and zoom in : 
    PixToMMzoom(histo, x0, x1, y0, y1, 0., 700.);  //void PixToMMzoom(TH2F * h1, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax); 

    // Prepare the fit function : 
    TF2 *f3 = new TF2("f3", fun3, x0, x1, y0, y1, npar); 
    Double_t paramsf3[npar] = {10000. , 10.5/1.1 , 9.0/1.1,   0.3        ,  20.   ,   109.    , 75.     , 575.    , 610.   , 10000. }; 
   
    //Double_t paramsf3[npar] =
    //{10000.,              2.4,   2.0,            0.3,         0.2        , 1.2        , 1.1      , 6.4    , 6.7    , 10000.     }; //2*pi/180*60.
    //amp_out[0]   , a[1] ,  b[2]     , theta[3]    , kbT[4]   , r_out[5], r_in[6], x0[7]  , y0[8], amp_in[9] 
    // Fuer die groesseren Ringe : {10000., 1.7, 1.05, 0.3, 0.2, 1.6, 1.5, 6.15, 7.5, 10000.}; //2*pi/180*60.
    //amp_out[0]   , a[1] ,b[2]     , theta[3]    , kbT[4]   , r_out[5], r_in[6], x0[7]  , y0[8], amp_in[9]   

    // Set the initial parameters : 
    f3->SetParameters(paramsf3);

    if(fithisto){
        f3->SetParNames("amp_out[0]", "a[1]", "b[2]", "theta[3]", "kbT[4]","r_out[5]", "r_in[6]", "x0[7]", "y0[8]", "amp_in[9]"); 

        // Fit histogram  with function f3 : 
        histo->Fit("f3", "0");

        // Retrieve the parameters from the fit function : 
        ampFit.push_back(f3 -> GetParameter(0));
        ampFitout.push_back(f3 -> GetParameter(9)); 
        a_Fit.push_back(f3->GetParameter(1)); 
        b_Fit.push_back(f3->GetParameter(2)); 
        thetaFit.push_back(f3->GetParameter(3)); 
        kbT_Fit.push_back(f3->GetParameter(4)); 
        r_outFit.push_back(f3->GetParameter(5)); 
        r_inFit.push_back(f3->GetParameter(6)); 
        x0Fit.push_back(f3->GetParameter(7)); 
        y0Fit.push_back(f3->GetParameter(8)); 
    }  
    // Set the function f3 with the parameters returned from the fit : 
    Double_t fitPass[npar] =
    {ampFit[RunCount], a_Fit[RunCount], b_Fit[RunCount], thetaFit[RunCount], kbT_Fit[RunCount],
    r_outFit[RunCount],r_inFit[RunCount], x0Fit[RunCount], y0Fit[RunCount], ampFitout[RunCount]};
    f3->SetParameters(fitPass);
        if(canvas1){
            can_histo->cd(1);  histo->Draw("contour2"); 
            can_histo->cd(2);  histo->Draw("surf2"); 
            gPad -> Modified(); gPad -> Update(); 
            gPad -> GetView() -> TopView(); 
        if(fithisto){
            can_histo->cd(3);  histo->Draw("contour2");
            f3->Draw("same"); 
            can_histo->cd(4);  histo->Draw("surf2");  
            f3->Draw("same");   
            gPad -> Modified(); gPad -> Update(); 
            gPad -> GetView() -> TopView(); 

        }
        if(savecanvas){
            char name_png[30];
            sprintf(name_png,"../root_files/rotated_coord_trans/png_files/%d_surf2_top_view_v2.png", iRun);
            can_histo -> Modified();   can_histo -> Update();
            can_histo -> SaveAs(name_png);
        }
        }
        if(canvas2){
            can_Fit -> cd(1); histo->Draw("surf2");  
            can_Fit -> cd(2);  
            f3->SetNpx(100);
            f3->SetNpy(100);
            f3->SetFillColorAlpha(5, 0.3); 
            histo->Draw("surf2");
            f3->Draw("surf2 same"); 
            if(savecanvas){
                char name1_png[30];
                sprintf(name1_png,"../root_files/rotated_coord_trans/png_files/%d_SurfacePlotDataAndFit.png", iRun);
                can_Fit -> Modified();   can_Fit -> Update();
                can_Fit -> SaveAs(name1_png);
            }
        }
        if(writehisto){
            char name[30];
            sprintf(name, "../root_files/rotated_coord_trans/%d_annulusFit.root", iRun);
            TFile* file1 = new TFile( name, "RECREATE");  
            histo -> Write(); 
            file1 -> Close(); 
            delete file1; 
        }
}


/*********************** Function or Method section **********************/

// Puck function : 
 Double_t g3(Double_t *x, Double_t *par) {                                             
  Double_t x1 = Double_t((x[0]-par[7])*cos(par[3])+(x[1]-par[8])*sin(par[3])); 
  Double_t x2 = Double_t(-(x[0]-par[7])*sin(par[3])+(x[1]-par[8])*cos(par[3])); 
  Double_t  a = Double_t(par[1]); 
  Double_t  b = Double_t(par[2]); 
  Double_t r1 = Double_t((x1)*(x1)/(a*a));
  Double_t r2 = Double_t((x2)*(x2)/(b*b));     
  Double_t A1 = Double_t((r1+r2-par[5])/par[4]);
  Double_t A2 = Double_t((r1+r2-par[6])/par[4]);
  return par[0]*(+1/(exp(A1)+1))-par[9]*(1/(exp(A2)+1));  
}

// Parameter assignment function for Puck function : 
 Double_t fun3(Double_t *x, Double_t *par) {                     // Both x and par are vectors. 
    Double_t *p1 = &par[0];                                      // p1 is a pointer to the address of parameter 0. 
    Double_t result = g3(x,p1); 
    return result;
 }

// Converts the histogram from pixels to mm and sets the range. 
  void PixToMMzoom(TH2F * h1, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax) {
    int sizeAxis = h1 -> GetXaxis() -> GetXmax();
    Double_t dimensions = ((sizeAxis-2)/43.5)*0.48;                  
    h1-> GetXaxis() -> SetLimits(0.,dimensions); 	
    h1-> GetYaxis() -> SetLimits(0.,dimensions); 
    h1-> GetXaxis()->SetTitle("mm");		
    h1-> GetYaxis()->SetTitle("mm");
    h1 -> GetZaxis() -> SetRangeUser(zmin, zmax);  
    h1 -> GetXaxis() -> SetRangeUser(xmin, xmax);                                       // Set zAxis to some value you like!!! 
    h1 -> GetYaxis() -> SetRangeUser(ymin, ymax);
    cout << "xmin: " << xmin << "\t" << "xmax: " << xmax << "\t" << (xmax-xmin) << endl; 
    cout << "ymin: " << ymin << "\t" << "ymax: " << ymax << "\t" << (ymax-ymin) << endl; 
    
  }
  
