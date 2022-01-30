#include<iostream>
#include<fstream>
#include<sstream>
#include<time.h>
#include<string>
#include<iomanip>


using namespace std;

//*****************************************************************************************
//
//      «апись значений 
//		средних R, дисперсий, кол-в R в интервале eps,
//      Inorm(S),H(S), D(r) в файлы
//
//*****************************************************************************************
void Write_Rij_Nrij_Dispij(string fileName, 
						   const int Ns,
						   const int Nir,
						   const int *SA1, const int *SA2, 
						   const int *NRij, const float *Rij, const float *DispRij, 
						   const string comments1, const string comments2, const string comments3)
{
	ofstream outf(fileName);
	outf << comments1 << endl << comments2 << endl << comments3 << endl;
	outf << setw(7)  << "i-j";
	outf << setw(7)  << "єпп"; 
	outf << setw(7)  << "N";
	outf << setw(15) << "Rij";
	outf << setw(15) << "DispRij" << endl;
		
	ostringstream str("");	


	for(int i = 1; i <= Ns; ++i)
		for(int j = i; j <= Ns; ++j)
		{
			int l = 1;
			for(int k = 0; k < Nir; ++k)
			{
				if(SA1[k] == i &&  SA2[k] == j)  
				{
					str.str("");
					str << i << "-" << j;
					outf << setw(7)  << str.str();
					outf << setw(7)  << l;
					outf << setw(7)  << NRij[k];
					outf << setw(15) << setprecision(5) << fixed << Rij[k];
					outf << setw(15) << setprecision(5) << fixed << DispRij[k];
					outf << resetiosflags(ios::scientific | ios::fixed);
					outf << endl;
					
					++l;
				}
			}
			outf << endl;
		}

}

void Write_Is(string fileName, 
			  const int NSma, const float Smin, const float dS,
			  const float *Io,
			  const string comments1, const string comments2, const string comments3)
{
	ofstream fo(fileName);

	fo  << comments1 << endl
		<< comments2 << endl
		<< comments3 << endl
		<< "I(S)" << endl
		<< NSma << " 30 30 1 0" << endl; 

	float S = Smin;
	for(int k = 0; k < NSma; ++k)
	{
		//S = Smin + k * dS;
		fo << setw(8)  << setprecision(4) << fixed << S;
		fo << setw(15) << setprecision(4) << fixed << Io[k] << endl;
		S = S + dS;
	}
}



void Write_Isnorm(string fileName, 
				  const int NSma, const float Smin, const float dS,
				  const float *Inorm,
				  const string comments1, const string comments2, const string comments3)



{
	ofstream fo(fileName);

	fo  << comments1 << endl
		<< comments2 << endl
		<< comments3 << endl
		<< "Inorm(S)" << endl
		<< NSma << " 3 3 1 0" << endl; 

	float S = Smin;
	for(int k = 0; k < NSma; ++k)
	{
		//S = Smin + k * dS;
		fo << setw(8)  << setprecision(4) << fixed << S;
		fo << setw(15) << setprecision(4) << fixed << Inorm[k] << endl;
		S = S + dS;
	}
}

void Write_Hs(string fileName, 
			  const int NSma, const float Smin, const float dS,
			  const float *Hs,
			  const string comments1, const string comments2, const string comments3)
{
	ofstream fo(fileName);

	fo  << comments1 << endl
		<< comments2 << endl
		<< comments3 << endl
		<< "H(S)" << endl
		<< NSma << " 1 1 1 0" << endl; 

	float S = Smin;
	for(int k = 0; k < NSma; ++k)
	{
		//S = Smin + k * dS;
		fo << setw(8)  << setprecision(4) << fixed << S;
		fo << setw(15) << setprecision(4) << fixed << Hs[k] << endl;
		S = S + dS;
	}
}

void Write_Dr(string fileName, 
			  const int NRma, const float delR, const float *Dr, 
			  const string comments1, const string comments2, const string comments3)
{
	ofstream fo(fileName);

	fo  << comments1 << endl
		<< comments2 << endl
		<< comments3 << endl
		<< "D(r)" << endl
		<< NRma << " 15 15 1 0" << endl; 

	float R = 0;
	for(int k = 0; k < NRma; ++k)
	{
		fo << setw(8)  << setprecision(3) << fixed << R;
		fo << setw(15) << setprecision(3) << fixed << Dr[k] << endl;
		
		R = R + delR;
	}
}

//логирование в файл
void LogTime(ofstream *fs, string text, clock_t logt, bool isConsOutput)
{
	float ftime_sec = (float)logt / CLOCKS_PER_SEC;
	
	*fs << 
		setw(60) << left << text << 
		setw(10) << setprecision(3) << fixed << "t0+сек: "<< ftime_sec << endl;
	if(isConsOutput)
		cout << 
		setw(60) << left << text << right <<
		setw(10) << setprecision(3) << fixed << "t0+сек: "<< ftime_sec << endl; 
}

//текущее дата, врем€ в виде строки
string DateTimeCur()
{
	time_t tsec = time(NULL);						//врем€ в секундах с 1970 года
	tm *datetime = datetime = localtime(&tsec);		//поучаем структуру с пол€ми даты, времени и т.п.
	char buff[80];                                  //строка с датой и временем
	
	strftime(buff, 80, "%d.%m.%Y %H:%M:%S", datetime);	//формируем строку даты времени типа "07.05.2013 12:20:10"

	return buff;

}