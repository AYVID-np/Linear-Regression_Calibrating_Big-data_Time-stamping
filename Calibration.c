void Calibration()
{
TFile *f1 = new TFile("TDC_cal1.root");
FILE *f2 = fopen("calib1.dat","w");

Int_t period;

fprintf(f2,"INTERCEPT                     SLOPE                  Chisq  \n"); 

//******************************  Function definition *********************

TF1 *func = new TF1("func","pol1(0)",0,4096);   

//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx Scintillator Detectors in Plane 1   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

for(Int_t i=27; i<=297; i=i+18)
   {
   if(i!=153 && i!=171 && i!=189) 
  {
 TH1D *h  = new  TH1D(Form("T_%03d",i),Form("T_%03d",i),4096,0,4096);
 TH1D *hnew  = new  TH1D(Form("Tn_%03d",i),Form("Tn_%03d",i),4096,0,4096);
       f1->GetObject(Form("TOF8%03d",i),h);

TSpectrum *s=new TSpectrum(10,0.002);    
Int_t nfound = s->Search(h,1);           
Double_t *peakX;
peakX=s->GetPositionX();

//***************************** Channel number [PositionX()] Sorting algorithm ***********
 
 long long sizes = nfound;

  long long inds[10];  
  for(int k=0; k<nfound; k++) inds[k] = 0;     
  TMath::Sort(sizes,peakX,inds,kFALSE);
  

period=0;
cout << "sorted values " << endl;
for(int j=0; j<nfound; j++)
    {
      period=period+40;
      hnew->SetBinContent(peakX[ inds[j] ],period);           
    }


//*****************************  Linear Fit function **************************

h->Draw();  

func->SetParameters(1,1);
hnew->Fit("func");                                   
  
func->GetParameter(0);                              
func->GetParError(0);
func->GetParameter(1);

fprintf(f2,"Detector %s \n",(Form("TOF8%03d",i))); 
fprintf(f2,"%f+-%f     %f+-%f    %f \n", func->GetParameter(0),func->GetParError(0),func->GetParameter(1),func->GetParError(1),func->GetChisquare());  
}

}


//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx Scintillator Detectors in Plane 2  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

for(Int_t i=18; i<=342; i=i+18)
   {
   if(i!=180 && i!=306 && i!=324) 
  {
     TH1D *h = new TH1D(Form("T_7%03d",i),Form("T_7%03d",i),4096,0,4096);
     TH1D *hnew = new TH1D(Form("T_7%03d_new",i),Form("T_7%03d_new",i),4096,0,4096);
     
       f1->GetObject(Form("TOF7%03d",i),h);

TSpectrum *s=new TSpectrum(10,0.002);    
Int_t nfound = s->Search(h,1);   

Double_t *peakX;
peakX=s->GetPositionX();

//***************************** Channel number [PositionX()] Sorting algorithm ***********
 
 long long sizes = nfound;

  long long inds[10];  
  for(int k=0; k<nfound; k++) inds[k] = 0;
  TMath::Sort(sizes,peakX,inds,kFALSE);
  

period=0;
cout << "sorted values " << endl;
for(int j=0; j<nfound; j++)
    {

      period=period+40;
      hnew->SetBinContent(peakX[ inds[j] ],period);           
    }


//*****************************  Linear Regression  **************************
hnew->Draw();  

func->SetParameters(1,1);
hnew->Fit("func");                                    

fprintf(f2,"Detector %s \n",(Form("TOF7%03d",i))); 
fprintf(f2,"%f+-%f     %f+-%f    %f \n", func->GetParameter(0),func->GetParError(0),func->GetParameter(1),func->GetParError(1),func->GetChisquare());  
}

}

//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx  Proportional Detector-1  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

TH1D *h11 = new TH1D("h11","h11",4096,0,4096);      TH1D *h12 = new TH1D("h12","h12",4096,0,4096);     TH1D *h13 = new TH1D("h13","h13",4096,0,4096);   TH1D *h14 = new TH1D("h14","h14",4096,0,4096);

TH1D *h11n = new TH1D("h11n","h11n",4096,0,4096);     TH1D *h12n = new TH1D("h12n","h12n",4096,0,4096);     TH1D *h13n = new TH1D("h13n","h13n",4096,0,4096);   TH1D *h14n = new TH1D("h14n","h14n",4096,0,4096);

    
       f1->GetObject("MWPC1_XR",h11);            f1->GetObject("MWPC1_XL",h12);             f1->GetObject("MWPC1_YU",h13);              f1->GetObject("MWPC1_YD",h14);     


TSpectrum *s11=new TSpectrum(10,0.002);     TSpectrum *s12=new TSpectrum(10,0.002);      TSpectrum *s13=new TSpectrum(10,0.002);      TSpectrum *s14=new TSpectrum(10,0.002);    
Int_t nfound11 = s11->Search(h11,1);   Int_t nfound12 = s12->Search(h12,1);   Int_t nfound13 = s13->Search(h13,1);   Int_t nfound14 = s14->Search(h14,1);   


Double_t *peakX11, *peakX12,*peakX13,*peakX14;
peakX11=s11->GetPositionX();         peakX12=s12->GetPositionX();      peakX13=s13->GetPositionX();     peakX14=s14->GetPositionX();

//***************************** Channel number [PositionX()] Sorting algorithm ***********
 
 long long sizes11 = nfound11;      long long sizes12 = nfound12;     long long sizes13 = nfound13;    long long sizes14 = nfound14;
  long long inds11[10];          long long inds12[10];      long long inds13[10];       long long inds14[10];   
  
for(int k=0; k<nfound11; k++) inds11[k] = 0;
  TMath::Sort(sizes11,peakX11,inds11,kFALSE);

for(int l=0; l<nfound12; l++) inds12[l] = 0;
  TMath::Sort(sizes12,peakX12,inds12,kFALSE);

for(int m=0; m<nfound13; m++) inds13[m] = 0;
  TMath::Sort(sizes13,peakX13,inds13,kFALSE);

for(int n=0; n<nfound14; n++) inds14[n] = 0;
  TMath::Sort(sizes14,peakX14,inds14,kFALSE);

period=0;
for(int j=0; j<nfound11; j++)
    {
     period=period+40;
      h11n->SetBinContent(peakX11[ inds11[j] ],period);           
     }

period=0;
for(int j=0; j<nfound12; j++)
    {
     period=period+40;
      h12n->SetBinContent(peakX12[ inds12[j] ],period);           
     }

period=0;
for(int j=0; j<nfound13; j++)
    {
     period=period+40;
      h13n->SetBinContent(peakX13[ inds13[j] ],period);  
     }

period=0;
for(int j=0; j<nfound14; j++)
    {
     period=period+40;
      h14n->SetBinContent(peakX14[ inds14[j] ],period);           
     } 


//*****************************  Linear Regression **************************
  
func->SetParameters(1,1);

h11n->Fit("func");                                                                                                          

fprintf(f2,"  MWPC1-XR \n");
fprintf(f2,"%f+-%f     %f+-%f    %f \n", func->GetParameter(0),func->GetParError(0),func->GetParameter(1),func->GetParError(1),func->GetChisquare());  
//-----------------------------------------------------------------------

h12n->Fit("func");

fprintf(f2," MWPC1-XL \n");
fprintf(f2,"%f+-%f     %f+-%f    %f \n", func->GetParameter(0),func->GetParError(0),func->GetParameter(1),func->GetParError(1),func->GetChisquare());  

//-----------------------------------------------------------------------

h13n->Fit("func");

fprintf(f2,"  MWPC1-YU \n");
fprintf(f2,"%f+-%f     %f+-%f    %f \n", func->GetParameter(0),func->GetParError(0),func->GetParameter(1),func->GetParError(1),func->GetChisquare());  

//-----------------------------------------------------------------------

h14n->Fit("func");

fprintf(f2," MWPC1-YD \n");
fprintf(f2,"%f+-%f     %f+-%f    %f \n", func->GetParameter(0),func->GetParError(0),func->GetParameter(1),func->GetParError(1),func->GetChisquare());  




//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx  Proportional Detector-2  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


TH1D *h21 = new TH1D("h21","h21",4096,0,4096);      TH1D *h22 = new TH1D("h22","h22",4096,0,4096);     TH1D *h23 = new TH1D("h23","h23",4096,0,4096);   TH1D *h24 = new TH1D("h24","h24",4096,0,4096);

TH1D *h21n = new TH1D("h21n","h21n",4096,0,4096);     TH1D *h22n = new TH1D("h22n","h22n",4096,0,4096);     TH1D *h23n = new TH1D("h23n","h23n",4096,0,4096);   TH1D *h24n = new TH1D("h24n","h24n",4096,0,4096);

    
       f1->GetObject("MWPC2_XR",h21);            f1->GetObject("MWPC2_XL",h22);             f1->GetObject("MWPC2_YU",h23);              f1->GetObject("MWPC2_YD",h24);     


TSpectrum *s21=new TSpectrum(10,0.002);     TSpectrum *s22=new TSpectrum(10,0.002);      TSpectrum *s23=new TSpectrum(10,0.002);      TSpectrum *s24=new TSpectrum(10,0.002);    
Int_t nfound21 = s21->Search(h21,1);   Int_t nfound22 = s22->Search(h22,1);   Int_t nfound23 = s23->Search(h23,1);   Int_t nfound24 = s24->Search(h24,1);   


Double_t *peakX21, *peakX22,*peakX23,*peakX24;
peakX21=s21->GetPositionX();         peakX22=s22->GetPositionX();      peakX23=s23->GetPositionX();     peakX24=s24->GetPositionX();

//***************************** Channel number [PositionX()] Sorting algorithm ***********
 
 long long sizes21 = nfound21;      long long sizes22 = nfound22;     long long sizes23 = nfound23;    long long sizes24 = nfound24;

  long long inds21[10];          long long inds22[10];      long long inds23[10];       long long inds24[10];   
  
for(int k=0; k<nfound21; k++) inds21[k] = 0;
  TMath::Sort(sizes21,peakX21,inds21,kFALSE);

for(int l=0; l<nfound22; l++) inds22[l] = 0;
  TMath::Sort(sizes22,peakX22,inds22,kFALSE);

for(int m=0; m<nfound23; m++) inds23[m] = 0;
  TMath::Sort(sizes23,peakX23,inds23,kFALSE);

for(int n=0; n<nfound24; n++) inds24[n] = 0;
  TMath::Sort(sizes24,peakX24,inds24,kFALSE);
 
period=0;
for(int j=0; j<nfound21; j++)
    {
     period=period+40;
      h21n->SetBinContent(peakX21[ inds21[j] ],period);          
     }

period=0;
for(int j=0; j<nfound22; j++)
    {
     period=period+40;
      h22n->SetBinContent(peakX22[ inds22[j] ],period);           
     }

period=0;
for(int j=0; j<nfound23; j++)
    {
     period=period+40;
      h23n->SetBinContent(peakX23[ inds23[j] ],period);  
     }

period=0;
for(int j=0; j<nfound24; j++)
    {
     period=period+40;
      h24n->SetBinContent(peakX24[ inds24[j] ],period);           
     } 


//*****************************  Linear Regression **************************
  

func->SetParameters(1,1);

h21n->Fit("func");                                                                                                          

fprintf(f2,"  MWPC2-XR \n");
fprintf(f2,"%f+-%f     %f+-%f    %f \n", func->GetParameter(0),func->GetParError(0),func->GetParameter(1),func->GetParError(1),func->GetChisquare());  
//-----------------------------------------------------------------------

h22n->Fit("func");

fprintf(f2," MWPC2-XL \n");
fprintf(f2,"%f+-%f     %f+-%f    %f \n", func->GetParameter(0),func->GetParError(0),func->GetParameter(1),func->GetParError(1),func->GetChisquare());  

//-----------------------------------------------------------------------

h23n->Fit("func");

fprintf(f2,"  MWPC2-YU \n");
fprintf(f2,"%f+-%f     %f+-%f    %f \n", func->GetParameter(0),func->GetParError(0),func->GetParameter(1),func->GetParError(1),func->GetChisquare());  

//-----------------------------------------------------------------------

h24n->Fit("func");

fprintf(f2," MWPC2-YD \n");
fprintf(f2,"%f+-%f     %f+-%f    %f \n", func->GetParameter(0),func->GetParError(0),func->GetParameter(1),func->GetParError(1),func->GetChisquare());  

//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


fclose(f2);

}



