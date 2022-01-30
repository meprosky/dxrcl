#include<iostream>
#include<fstream>
#include<sstream>
#include<iomanip>
#include<time.h>
#include<string>
#include<locale>
#include<math.h>
#include<time.h>
#include<algorithm>
#include"func.h"

using namespace std;

int main(int argc, char *argv[])
{

	setlocale(LC_CTYPE, "");			//����� � ������� ���� ������� �����

	string comments1, comments2,	    //��� ������ ��������� ��. �����
		comments3, comments4;			//��� ������ ��� ���. �����

	int Na, Ns;							//����� ������, ����� ������
	float lambdaXray;					//����� ����� �����. ���������

	int Nparc;							//����� ����������� ������������� ��� ������, ��������� 
										//� ����� ��� "��". ���� nparcr = 0, �� �� ������� � ����. ������
										//������ parc ������ �������������

	int *parc;							//������ ��� ������

	float 
		Smin, Smax, dS, 
		Alpha, Drmax, delR, R0e, eps;  //Smin, Smax, dS, Al, Rmax, dR, R0e, eps

	float CHFE;							//����� ���������� ������ � ��������

	int   *ps;							//������ ������ ������

	float *x, *y, *z;					//���������� ������ x, y, z

	int Nir;							//����� ���������� �� R
	int NSma;							//����� ���������� �� S
	float *Rij_0, *DispRij_0;			//"�����" ������ ���� R, ���������
	int *NRij_0;						//���-��
	float *Rij, *DispRij;				//������� R, ���������
	int *NRij;							//���-�� �� R
	int *SA1, *SA2;						//�������. ������� ��������
	
	int Nreps, Npairs;					//���-�� ���������� �� eps(����������), ��� ������������

	float *DF1, *DF2;					 //������ � ������ ������������� ��������
	
	//9 ������������� �� ������ ���� ����� 
	//��� ������� ������� �������� ��������� 
	float 
		*FA1, *FA2, *FA3, *FA4, 
		*FB1, *FB2, *FB3, *FB4, *FC, *FE;

	cout << endl;
	cout << "dxrcl_p9.exe" << endl;
	cout << "������  �������������  ������������� ���������  �������������" << endl;
    cout << "����� �����������, ���������� �� ������������-���������������" << endl; 
	cout << "���������. ��������������  ���������������  ����������� 2013." << endl;
	//cout << "��������� ����������� � �������������� �� C++ ���������� �.�." << endl;
	//cout << "�� ������ ���� dxrcl44.for �� ����� Fortran �� ���������� �������� �.�." << endl;
	//cout << endl;
	cout << endl;

	clock_t t0 = clock();
	string curdt = DateTimeCur();
	ostringstream ostr;
	
	string inputfn;

	if(argc == 1)         //��������� ������ �������
	{
		std::cout << "���������� ��� �������� ����� ������..." << endl;
		getchar();
		exit(1);
		//inputfn = "pdxr.xyz";
	}
	else
	{
		inputfn = argv[1];
	}
	
	string bufstr = inputfn.substr( inputfn.find_last_of(".") + 1);									//���� ���������� .dxr
	transform(bufstr.begin(), bufstr.end(), bufstr.begin(), ::tolower);								//��������� ��� � ������ �������
	string outputfn = (bufstr == "dxr") ? inputfn.substr(0, inputfn.find_last_of(".")) : inputfn;   //���� �����, dxr ������������ � �������� ��� �����
	
	ofstream logf(outputfn + ".log");		//���� ������ ���� ������

	bufstr = DateTimeCur();                 //������� ����, ����� � �����. ������
	ostr.str("");
	ostr << "   ������ ������     - " + bufstr;
	LogTime(&logf, ostr.str(), 0, true); 

	//************************************************************************
	//
    //     ������ ������� ������
	//
	//************************************************************************
	ReadData(inputfn,
		comments1, comments2,
		Na, Ns, lambdaXray,
		DF1, DF2, FA1, FA2, FA3, FA4, FB1, FB2, FB3, FB4, FC, FE,
		Nparc, parc,
		Smin, Smax, dS, Alpha, Drmax, delR, R0e, eps, CHFE,
		ps, 
		x, y, z
		);
	
	comments3 = "���� � ����� - " + bufstr;
	LogTime(&logf, "������ " + inputfn + " ���������", clock() - t0, true); 
	
	Npairs = Ns * (Ns + 1) / 2;				// ����� ������������ ������ ��� ������
	
	//************************************************************************
    //     ��������������� ������� �������� ��� ������� � ��������� ��������
	//************************************************************************
	int Ns1 = Ns + 1;
	int *idxA = new int[Ns1*Ns1];
	int *invidxSA1 = new int[Npairs];
	int *invidxSA2 = new int[Npairs];
	int cc = 0;
	for(int i = 1; i <= Ns; ++i)         //��������� ������� �������� ��� ������� � ��������� ��������
	{
		for(int j = i; j <=Ns; ++j)
		{
			idxA[i*Ns1+j] = cc; 
			idxA[j*Ns1+i] = cc;
			invidxSA1[cc] = i;
			invidxSA2[cc] = j;

			++cc;
		}
	}

	//************************************************************************
    //     ���� ����. ��������� ��������
	//************************************************************************
	float Rmax = findRmax(Na, x, y, z);

	ostr.str("");
	ostr << Rmax;
	LogTime(&logf, "����. �����. � �������� Rmax = " + ostr.str(), clock() - t0, true); 

	Nreps = (int)(Rmax / eps + 1.1);    //���������� ���������� eps � ����. ���������
	
	int Nall = Npairs * Nreps;          //����� �������� � ������ ������������ ��� (1-2 = 2-1)
	                                    //�������� ��� 3 ������ ������ ����� 1-1 1-2 1-3 2-2 2-3 3-3 ����� 6 ���

	ostr.str("");
	ostr << "����� ���. ���������� ����������. ������ = " << Na; 
	LogTime(&logf, ostr.str(), clock() - t0, true); 
	
	//************************************************************************
    //		����� ��������� �� ������ ���������� (��� ������� �� �� N(N-1)/2 )
	//		���������� ������������ ���������� ���������� � ��������
    //		� ������������� ����� ���������� (� �������� EPS) ����������
    //		����� ����������� ������ ������ (A - B � B - A - ���������).
    //		���. ������� �����. ���������� � ��������� ����� ����������
    //		�� ������������ ���������� � �������� EPS.
    //************************************************************************
	Rij_0   =    new float[Nall];     //���������� ����������
	DispRij_0  = new float[Nall];     //��������� ����������
	NRij_0  =    new   int[Nall];     //���-�� ���������� � ��������� eps
	
	unsigned int sz_float = Nall * sizeof(float);
	unsigned int sz_int   = Nall * sizeof(int);
	
	memset(Rij_0,      0, sz_float);   //��������
	memset(DispRij_0,  0, sz_float);
	memset(NRij_0,     0, sz_int);
	
	interatomic(
		Na, x, y, z, 
		Ns, Npairs, Nall, eps, ps, idxA,
		Rij_0, DispRij_0, NRij_0);


	ostr.str("");
	ostr << "������. ���. ���������� ����������. ������ = " << Na*(Na-1)/2;
	LogTime(&logf, ostr.str(), clock() - t0, true); 

	//************************************************************************
	//	������������ ����� ��������� ��������� Nir
	//	�.�. ��� �����. ���������� ������� ����������� ���� �� ���� ���
	//	�������� ��� ��������� ������ ��� �������� ����� ��� ������
	//************************************************************************
	Nir = 0; 
	for(int i = 0; i < Nreps; ++i){      //��������� �� ���������� eps
		for(int j = 0; j < Npairs; ++j){ //��������� �� ����������� ��� ������ 1-1 1-2 1-3 2-2 2-3 � �.�.
			if(NRij_0[i * Npairs + j]){  //�� ����� 0
				++Nir;                  //������� ��������� ����������
			}
		}
	}
	

		
	Rij     = new float[Nir];     //���������� ����������
	DispRij = new float[Nir];     //��������� ����������
	NRij    = new   int[Nir];     //���-�� ���������� � ��������� eps
	SA1     = new   int[Nir];     //������ ���� ����� � ����
	SA2     = new   int[Nir];     //������ ���� ����� � ����

	//****************************************************************************
	//
	//   ���������� ������� ����������, �� ���������� � ���������, ��������� 
	//   ������� "c������" ������, �.�. ��������� �������� � NRij = 0 (��� ����������)
	//
	//****************************************************************************
	calcNrAvrDisperse(
		Nir, Npairs, Nreps, invidxSA1, invidxSA2, 
		NRij_0, Rij_0, DispRij_0,
		SA1, SA2, 
		NRij, Rij, DispRij);
		
	delete [] invidxSA1; //����������� ������
	delete [] invidxSA2;
	delete [] Rij_0;                                       
	delete [] NRij_0;
	delete [] DispRij_0;


	NSma = (int)((Smax - Smin) / dS + 1.1f);       //���-�� ���������� dS � ��������� Smin-Smax

	//****************************************************************************
	//
	//   ���������� �������� ��� � ����������� �� S � �� ��������� ��
	//   ����� ������ � ��������, �� �������������� �� �����������, ���������� � �.�.
	//
	//****************************************************************************
	float *Gf2_expSAlpha = new float[NSma];           //�������� �-�� ������������ G-������� � ���-�� �� S ������ ���������� � ������� ��� �������������
	
	//����-�� ���  � ���-�� �� S
	//������ ��� - ���� 3� ������ ������, ������ - �������� S 
	//� �����. ��� 2� ������ ������ ��� i � j �����
	float *KFAR = new float[NSma * Ns * Ns];
	
	float GfactorSum1;

	calcFAR(Ns, NSma, Smin, dS, Alpha,
		DF1, DF2, FA1, FA2, FA3, FA4, FB1, FB2, FB3, FB4, FC,  FE,
		Gf2_expSAlpha,     //������������ ��������
		GfactorSum1,       //������������ ��������
		KFAR               //������������ ��������
	);

	ostr.str("");
	ostr << "������. ���. ���. ������ ������� = " << NSma << "*" << Ns << "*" << Ns << "=" << NSma*Ns*Ns;
	LogTime(&logf, ostr.str(), clock() - t0, true); 

	//*****************************************************************************************
	//
	//	 �������� ����
	//   ���������� ������� I(s), H(s), D(r)
	//   
	//*****************************************************************************************
	float *Io    = new float[NSma];
	float *Inorm = new float[NSma];
	float *Hs    = new float[NSma];

	calcINTS(Na,  Ns,  Nir, 
		ps,  SA1,  SA2, 
		NSma, Smin,  dS,  CHFE, 
		NRij,  Rij,  DispRij, 
		Gf2_expSAlpha, KFAR,
		Io,	 Inorm, Hs);

	//*****************************************************************************************
	//      ������ D(r)
	//*****************************************************************************************
	int NRma = (int)(Drmax / delR + 1.1f); //���-�� ���������� ��� D(r)
	float *Dr = new float[NRma];       
	
	calcDr( Drmax, delR, NSma, Smin, dS, R0e, GfactorSum1, Hs, 
		NRma, Dr);

	ostr.str("");
	ostr << "������. ���. I(s) H(s) D(r)";
	LogTime(&logf, ostr.str(), clock() - t0, true); 

	ostr.str("");
	ostr << (clock() - t0)/CLOCKS_PER_SEC;
	
	comments3 += ". ����. ��������,���= " + ostr.str();

	//*****************************************************************************************
	//
	//      ������ �������� 
	//		������� R, ���������, ���-� R � ��������� eps,
	//      I(S), Inorm(S),H(S), D(r) � �����
	//
	//*****************************************************************************************
	Write_Rij_Nrij_Dispij(outputfn + ".ndr", 
						   Ns, Nir, SA1, SA2, 
						   NRij, Rij, DispRij, 
						   comments1, comments2, comments3);
	
	Write_Is(outputfn + "_Is.xy",
		NSma, Smin, dS, Io,
		comments1, comments2, comments3);

	Write_Isnorm(outputfn + "_Isnorm.xy",
		NSma, Smin, dS, Inorm,
		comments1, comments2, comments3);

	Write_Hs(outputfn + "_Hs.xy",
		NSma, Smin, dS, Hs,
		comments1, comments2, comments3);

	Write_Dr(outputfn + "_Dr.xy", 
		NRma, delR, Dr, 
		comments1, comments2, comments3);
		
	//����������� ������
	delete [] Rij;                                      
	delete [] NRij;
	delete [] DispRij;

	delete [] Io;
	delete [] Inorm;
	delete [] Hs;
	delete [] Dr;

	delete [] KFAR;
	delete [] Gf2_expSAlpha;
	delete [] idxA;
	delete [] SA1;
	delete [] SA2;


	delete [] DF1;  
	delete [] DF2;  
	delete [] FA1;  
	delete [] FA2; 
	delete [] FA3;
	delete [] FA4;
	delete [] FB1;
	delete [] FB2;
	delete [] FB3;
	delete [] FB4;
	delete [] FC;
	delete [] FE; 

	ostr.str("");
	ostr << "   ���������� ������ - " + DateTimeCur();
	LogTime(&logf, ostr.str(), clock() - t0, true); 
	//cout <<  "\n������� ����� ������� ��� ������...";
	//getchar();
	return 0;
}
