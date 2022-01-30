#include<string>
#include<time.h>

using namespace std;

void ReadData(
	 string filename, string &comments1, string &comments2,
	 int& Na, int& Ns, float& lambdaXray,
	 float* &DF1, float* &DF2, 
	 float* &FA1, float* &FA2, float* &FA3, float* &FA4, 
	 float* &FB1, float* &FB2, float* &FB3, float* &FB4, 
	 float* &FC,  float* &FE,
	 int &Nparc, int* &parc,
	 float& Smin, float& Smax, float& dS, float& Alpha, float& Drmax, float& delR, float& R0e, float& eps, float &CHFE,
	 int* &ps, float* &x, float* &y, float* &z);

float findRmax(const int Na, const float *x, const float *y, const float *z);

void interatomic(const int Na, const float *x, const float *y, const float *z, 
		const int Ns, const int Npairs, const int Nall, const float eps, const int *ps, const int *idxA,
		float* &Rij_0, float* &DispRij_0, int* &NRij_0);

void calcNrAvrDisperse(const int Nir, const int Npairs, const int Nreps, 
		const int *invidxSA1, const int *invidxSA2,
		int* &NRij_0, float* &Rij_0, float* &DispRij_0, 
		int* &SA1, int* &SA2,                            //возвращаемое зн.
		int* &NRij, float* &Rij, float* &DispRij);       //возвращаемое зн.


void calcFAR(
	const int Ns, 
	const int NSma, const float Smin, const float dS, const float Alpha,
	const float *DF1, const float *DF2,
	const float *FA1, const float *FA2, const float *FA3, const float *FA4, 
	const float *FB1, const float *FB2, const float *FB3, const float *FB4, 
	const float *FC,  const float *FE,
	float* &Gf2_expSAlpha,     //возвращаемое значение
	float  &GfactorSum1,       //возвращаемое значение
	float* &KFAR               //возвращаемое значение
	);

void calcINTS(
	const int Na, const int Ns, const int Nir, 
	const int *ps, const int *SA1, const int *SA2,  
	const int NSma, const float Smin, const float dS, const float CHFE,
	const int *NRij, const float *Rij, const float *DispRij,
	const float *Gf2_expSAlpha,     
	const float *KFAR,
	float* &Io,                  //возвращаемое значение
	float* &Inorm,				//возвращаемое значение
	float* &Hs					//возвращаемое значение
	);

void calcDr(
	const float Drmax, const float delR, 
	const int NSma, const float Smin, const float dS, 
	const float R0e, 
	const float GfactorSum1,  
	const float *Hs,  
	int NRma,          
	float* &Dr           //возвращаемое значение
	);


void Write_Rij_Nrij_Dispij(string fileName, 
						   const int Ns,
						   const int Nir,
						   const int *SA1, const int *SA2, 
						   const int *NRij, const float *Rij, const float *DispRij, 
						   const string comments1, const string comments2, const string comments3);

void Write_Is(string fileName, 
			  const int NSma, const float Smin, const float dS,
			  const float *Io,
			  const string comments1, const string comments2, const string comments3);

void Write_Isnorm(string fileName, 
				  const int NSma, const float Smin, const float dS,
				  const float *Inorm,
				  const string comments1, const string comments2, const string comments3);

void Write_Hs(string fileName, 
			  const int NSma, const float Smin, const float dS,
			  const float *Hs,
			  const string comments1, const string comments2, const string comments3);

void Write_Dr(string fileName, 
			  const int NRma, const float delR, const float *Dr, 
			  const string comments1, const string comments2, const string comments3);

void LogTime(ofstream *fs, string text, clock_t logt, bool isConsOutput);

string DateTimeCur();

