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

//Global variables : 
int RunCount=0; 
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


// Functions : 
Double_t g3(Double_t *x, Double_t *par); 
Double_t fun3(Double_t *x, Double_t *par); 
void PixToMMzoom(TH2F * h1, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax); 




// Body of the code : 
 void annulus_fit_v4(int iRun){

    // Set style : 
    gStyle->SetPalette(1);

    // Parameters for fit. 
    Double_t pi = 3.14159265359; 
    const Int_t npar = 10;


    // Setup Canvas : 
    Double_t w = 1400;
    Double_t h = 700;
    TCanvas * can_histo = new TCanvas("elliptical annulus", "study of the elliptical annulus ", w, h);       // name, title, width, height. 
    can_histo -> Divide(2,1); 
    
    // Get Image : 
    char nameFile[60];
    sprintf(nameFile, "../root_files/raw_data/%d_annulus.root", iRun);
    TFile* file = new TFile(nameFile, "READ"); 
    //TH2F * histo = (TH2F*)file->Get("histo_gauss"); // How is the original histogram called? 
    TH2F * histo = (TH2F*)file->Get("ImgFull5THDetCMOS_0");


    // Declaration and definition of the parameters x0, x1, y0, y1 
    // Can I here limit the range of the fit function? 
    /*Double_t x0 = 400.; 
    Double_t x1 = 750.;
    Double_t y0 = 450.; 
    Double_t y1 = 750.; */

    Double_t x0 = 0.; 
    Double_t x1 = 1148.;
    Double_t y0 = 0.; 
    Double_t y1 = 1148.; 

    //histo -> GetXaxis() -> SetRangeUser(x0, x1);                                       // Set zAxis to some value you like!!! 
    //histo -> GetYaxis() -> SetRangeUser(y0, y1);



    TF2 *f3 = new TF2("f3", fun3, x0, x1, y0, y1, npar); 

    /*
    Double_t paramsf3[npar] = {10000. , 105. , 90., 0.3 , 10. , 109. , 75. , 575. ,610. , 10000. }; 
    //amp_out[0]   , a[1] ,  b[2]     , theta[3]    , kbT[4]   , r_out[5], r_in[6], x0[7]  , y0[8], amp_in[9] 
    f3 -> SetParameters(paramsf3); 
    f3->SetParNames("amp_out[0]", "a[1]", "b[2]", "theta[3]", "kbT[4]","r_out[5]", "r_in[6]", "x0[7]", "y0[8]", "amp_in[9]"); 
*/

    Double_t paramsf3[npar] = {10000. , 10.5/1.1 , 9.0/1.1,   0.3        ,  20.   ,   109.    , 75.     , 575.    , 610.   , 10000. }; 
                      //amp_out[0]    , a[1] , b[2], theta[3]    , kbT[4]   , r_out[5], r_in[6] , x0[7]   , y0[8]  , amp_in[9] 
    f3 -> SetParameters(paramsf3); 
    f3->SetParNames("amp_out[0]", "a[1]", "b[2]", "theta[3]", "kbT[4]","r_out[5]", "r_in[6]", "x0[7]", "y0[8]", "amp_in[9]"); 

    //Fit the histogram with the Puck function f3. 
    histo -> Fit("f3", "MR0"); 

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

    Double_t fitPass[npar] =
    {ampFit[RunCount], a_Fit[RunCount], b_Fit[RunCount], thetaFit[RunCount], kbT_Fit[RunCount],
    r_outFit[RunCount],r_inFit[RunCount], x0Fit[RunCount], y0Fit[RunCount], ampFitout[RunCount]};
    f3->SetParameters(fitPass);



    // Draw image: 
    can_histo -> cd(1); 
    histo -> Draw("surf2"); 
    gPad -> Modified(); gPad -> Update(); 
    gPad -> GetView() -> TopView();
    can_histo -> cd(2); 
    histo -> Draw("surf2"); 

    //Draw fit: 
    can_histo -> cd(1);
    f3->Draw("same");

    can_histo -> cd(2); 
    f3->SetNpx(100);
    f3->SetNpy(100);
    f3->SetFillColorAlpha(5, 0.3); 
    f3->Draw("surf2 same"); 

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
  
