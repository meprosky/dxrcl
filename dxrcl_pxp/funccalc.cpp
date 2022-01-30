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

#define PI 3.14159265f

using namespace std;

//************************************************************************
//
//     Чтение входных данных
//
//************************************************************************
void ReadData(
	 string filename, string &comments1, string &comments2,
	 int& Na, int& Ns, float& lambdaXray,
	 float* &DF1, float* &DF2, 
	 float* &FA1, float* &FA2, float* &FA3, float* &FA4, 
	 float* &FB1, float* &FB2, float* &FB3, float* &FB4, 
	 float* &FC, 
	 float* &FE,
	 int& Nparc, int* &parc,
	 float& Smin, float& Smax, float& dS, float& Alpha, float& Drmax, float& delR, float& R0e, float& eps, float &CHFE,
	 int* &ps, float* &x, float* &y, float* &z)
{
	ifstream inputDataFile(filename);
	stringstream strstream;
	

	if(!inputDataFile)
	{
		cout << "Ошибка чтения входного файла" << endl;
		getchar();
		exit(1);
	}

	strstream << inputDataFile.rdbuf();  //весь файл в строковый поток
	inputDataFile.close();               //вх. файл больше не нужен, закрываем 
	
	string inputdata = strstream.str();  //весь поток направляем в строку
	
	getline(strstream, comments1); //читаем строку 1
	getline(strstream, comments2); //читаем строку 2

	strstream >> Na >> Ns; //число атомов, число сортов

	strstream >> lambdaXray; // длина волны рентг. излучения

	DF1 = new float[Ns];  //первая дисперсионная поправка
	DF2 = new float[Ns];  //вторая дисперсионная поправка

	for(int i = 0; i < Ns; ++i) strstream >> DF1[i];
	for(int i = 0; i < Ns; ++i) strstream >> DF2[i];

	FA1 = new float[Ns];  //9 коэффициентов на каждый сорт атома 
	FA2 = new float[Ns];  //для расчета функций атомного рассеяния 
	FA3 = new float[Ns];
	FA4 = new float[Ns];
	FB1 = new float[Ns];
	FB2 = new float[Ns];
	FB3 = new float[Ns];
	FB4 = new float[Ns];
	FC  = new float[Ns];

	for(int i = 0; i < Ns; ++i) strstream >> FA1[i];
	for(int i = 0; i < Ns; ++i) strstream >> FA2[i];
	for(int i = 0; i < Ns; ++i) strstream >> FA3[i];
	for(int i = 0; i < Ns; ++i) strstream >> FA4[i];
	for(int i = 0; i < Ns; ++i) strstream >> FB1[i];
	for(int i = 0; i < Ns; ++i) strstream >> FB2[i];
	for(int i = 0; i < Ns; ++i) strstream >> FB3[i];
	for(int i = 0; i < Ns; ++i) strstream >> FB4[i];
	for(int i = 0; i < Ns; ++i) strstream >> FC[i];

	FE = new float[Ns];  //"формульная единица" или единица состава 
	                     // Например, Al2O3 - FE(1) = 2., FE(2) = 3. 
	                     //( Признак сорта атома PS -  Al - 1, O - 2) 
	
	for(int i = 0; i < Ns; ++i) strstream >> FE[i];

	strstream >> Nparc; //число парциальных распределений пар атомов, выводимых 
	                    //в файлы для "ху". Если nparcr = 0, то не выводим и след. строка
	                    //должна отсутствовать
	

	if(Nparc > 0)
	{
		parc = new int[Nparc];
		string s;
		getline(strstream, s);
		getline(strstream, s);
		
		replace(s.begin(), s.end(), ',', ' ');
		
		stringstream ss2(s);

		for(int i = 0; i < Nparc; ++i) ss2 >> parc[i];
		
	}
	else if(Nparc < 0)
	{
		cout << "Ошибка числа парциальных распределений Nparc";
	    exit(1);
	}
	else //0
	{
	}
	
	string s;
	getline(strstream, s);
	replace(s.begin(), s.end(), ',', ' ');
	
	stringstream ss2(s);
	ss2 >> Smin >> Smax >> dS >> Alpha >> Drmax >> delR >> R0e >> eps; //Smin, Smax, dS, Al, Rmax, dR, R0e, eps
	
	strstream >> CHFE; // число формульных единиц в кластере

	ps = new   int[Na];   //сорт атома
	x =  new float[Na];   //координаты атомов x, y, z
	y =  new float[Na];
	z =  new float[Na];

	for(int i = 0; i < Na; ++i) 	
		strstream >> ps[i] >> x[i] >> y[i] >> z[i];

}

//****************************************************************************
//
//   Вычисление значений ФАР в зависимости от S и не зависимых от
//   числа атомов в кластере, от распеределений по расстояниям, дисперсиям и т.п.
//
//****************************************************************************
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
	)
{

	
	float L4Pi = 1.0f / (4.0f * PI);
	float ST;
	
	int Ns2 = Ns * Ns;

	GfactorSum1 = 0.0f
		;
	for(int i = 0; i < Ns; ++i)
	{
		GfactorSum1 = GfactorSum1 + FE[i] * (FC[i] + FA1[i] + FA2[i] + FA3[i] + FA4[i]);
	}
	float GfactorSum1_inv = 1 / GfactorSum1;           //получаем обратную величину 
	                                                   //потом будем умножать (это быстрее чем деление)

	float S = Smin;
	for(int k = 0; k < NSma; ++k)
	{
		ST = S * L4Pi;
		ST = ST * ST;

		int _idxs = k * Ns2;                   //всп. индекс в 3х мерном массиве

		float GfactorSum2 = 0.0f;

		for (int i = 0; i < Ns; ++i)
		{
			float KFARFi = 
				FC[i] +
				FA1[i] * exp(-FB1[i] * ST) +			
				FA2[i] * exp(-FB2[i] * ST) +
				FA3[i] * exp(-FB3[i] * ST) +
				FA4[i] * exp(-FB4[i] * ST) + 
				DF1[i];

			float KP1 = KFARFi * KFARFi + DF2[i] * DF2[i];
			
			int _idx = i * Ns;             //всп. индекс в 3х мерном массиве

			KFAR[_idxs + _idx + i] = KP1;  //сорт i-i
			
			//сумма для расчета обостряющего G фактора
			KP1 = sqrtf(KP1);
			GfactorSum2 = GfactorSum2 + FE[i] * KP1;

			for (int j = i+1; j < Ns; ++j)
			{
				float KFARFj = 
					FC[j] + 
					FA1[j] * exp(-FB1[j] * ST) +			
					FA2[j] * exp(-FB2[j] * ST) +
					FA3[j] * exp(-FB3[j] * ST) +
					FA4[j] * exp(-FB4[j] * ST) + 
					DF1[j];

				float KP2 = 
					KFARFi * KFARFj + DF2[i] * DF2[j];
				
				KFAR[_idxs + _idx + j] = KP2;           //элемент матрицы - сорт i-j
				KFAR[_idxs + j*Ns + i] = KP2;           //элемент матрицы - сорт j-i (в принципе можно убрать)

			}
		}

		float Gf = GfactorSum2 * GfactorSum1_inv;
		float sAlpha = S * Alpha;
		
		Gf2_expSAlpha[k] = expf(- sAlpha * sAlpha) / (Gf * Gf);		// эксп деленная на квадрат обостряющего G фактор 
		                                                            //предвычисляем для дальнейшего расчета H(s)

																	//Hs[k] =            интерференционная функция
																	//  S * (Inorm[k] - IndScat / CHFE) * Gf2_expSAlpha[k]

		S = S + dS;
	}



}


//************************************************************************
//     Ищем макс. диагональ кластера
//************************************************************************
float findRmax(const int Na, const float *x, const float *y, const float *z)
{
	float minx = x[0], miny = y[0], minz = z[0];
	float maxx = x[0], maxy = y[0], maxz = z[0];
	
	for(int i = 0; i < Na; ++i) //ищем мин и макс значения координат
	{
		if(x[i] < minx) minx = x[i];
		if(y[i] < miny) miny = y[i];
		if(z[i] < minz)	minz = z[i];

		if(x[i] > maxx) maxx = x[i];
		if(y[i] > maxy) maxy = y[i];
		if(z[i] > maxz)	maxz = z[i];
	}
	
	float dx=maxx-minx, dy=maxy-miny, dz=maxz-minz;
	
	return sqrtf(dx*dx + dy*dy + dz*dz);            //макс. диагональ
}

void calcNrAvrDisperse(const int Nir, const int Npairs, const int Nreps, 
		const int *invidxSA1, const int *invidxSA2,
		int* &NRij_0, float* &Rij_0, float* &DispRij_0, 
		int* &SA1, int* &SA2,                         //возвращаемое зн.
		int* &NRij, float* &Rij, float* &DispRij    //возвращаемое зн.
		)
{
	int c = 0;
	float Rr;

	for(int i = 0; i < Nreps; ++i)							//пробегаем по интервалам eps
	{
		for(int j = 0; j < Npairs; ++j)						//пробегаем по комбинациям пар атомов 1-1 1-2 1-3 2-2 2-3 и т.д.
		{
			int k = i * Npairs + j;
			int Nr = NRij_0[k];
			
			if(Nr > 1)
			{
				NRij[c] = Nr;
				
				Rr = Rij[c] = Rij_0[k] / Nr;				//среднее расстояние в интервале n*Eps - (n+1)*Eps
				
				DispRij[c] =								//дисперсия
					sqrtf(
					fabsf(DispRij_0[k] - (float)Nr*Rr*Rr) /
					((float)Nr * (Nr - 1))
					); 
				
				SA1[c] = invidxSA1[j];						//сорт атома 1
				SA2[c] = invidxSA2[j];						//сорт атома 1 
				
				++c;

			}
			else if(Nr == 1)								//если расстояние в интервале n*Eps - (n+1)*Eps всего одно то дисперсия = 0
			{
				NRij[c] = 1;				
				Rij[c] = Rij_0[k];
				DispRij[c] = 0.0f;
				SA1[c] = invidxSA1[j];
				SA2[c] = invidxSA2[j];
				++c;
			}
		}
	}




	}

//************************************************************************
//     Вычисление совокупности межатомных расстояний в кластере
//     с суммированием числа одинаковых (в пределах EPS) расстояний
//     между одинаковыми парами атомов (A - B и B - A - одинаковы).
//     Выч. среднее межат. расстояние и дисперсия этого расстояния
//     по совокупности расстояний в пределах EPS.
//************************************************************
void interatomic(const int Na, const float *x, const float *y, const float *z, 
		const int Ns, const int Npairs, const int Nall, const float eps, const int *ps, const int *idxA,
		float* &Rij_0, float* &DispRij_0, int* &NRij_0)
{
	
	int Nam1 = Na - 1;
	int Ns1 = Ns + 1;
	int idxa, idxeps, idxG;
	
	//хранение необходимо по столбцам, чтобы были в порядке возрастания R
	//в строках интервалы по eps,  в столбцах пары атомов
	
	float eps_inv = 1.0f / eps;
	float dr, dr_p2, dx, dy, dz;
	

	//****************************************************************************
	//
	//   Вычисление межатомных расстояний, квадратов для дисперсий,
	//   сумм расстоняий
	//
	//****************************************************************************
	for(int i = 0; i < Nam1; ++i)
		{
			for(int j = i + 1; j < Na; ++j)
			{

			dx = x[j] - x[i];
			dy = y[j] - y[i];
			dz = z[j] - z[i];

			dr_p2 = dx*dx + dy*dy + dz*dz;     //квадрат
			dr = sqrtf(dr_p2);
			
			idxa = idxA[ ps[i] * Ns1 + ps[j] ];

			idxeps = (int)(dr * eps_inv);    //индекc чтоб попасть в интервал n*eps - (n+1)*eps

			idxG = idxeps * Npairs + idxa;   //глобальный индекс в массивах по Eps

			Rij_0[idxG] += dr;               //сумма расстояний, потом вычислим среднее
			
			DispRij_0[idxG] += dr_p2;        //для выч. дисперсии			
			
			NRij_0[idxG] += 1;               //подсчет кол-ва интервалов расстояний
			
			}
		}

	}

//*****************************************************************************************
//
//   Вычисление функций I(s), H(s), D(r)
//   
//*****************************************************************************************
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
	)
{

	//*****************************************************************************************
	//
	//   Блок вычислений функций I(s), H(s), D(r)
	//   
	//   Замечание! Поскольку сорта атомов начинаются с 1 (1,2,3...и т.д.) а массивы все с нуля,
	//   то при доступе к элементу массива зависящему от сорта атома, от сорта атома мы отнимаем 1
	//
	//*****************************************************************************************

	float S = Smin;
	int Ns2 = Ns * Ns;

	for(int k = 0; k < NSma; ++k)
	{
		int _idxs = k * Ns2;                                   //часть индекса для доступа в массив ФАР
		
		//*****************************************************************************************
		//     Расчет независимого рассеяния
        //		         Na
        //        I(s) = SUM fk(s)*fk(s)
        //               k=1
		//*****************************************************************************************

		float IndScat = 0.0f;									//слагаемое независимого рассеяния
		
		for(int i = 0; i < Na; ++i)
		{
			int _is = ps[i] - 1;								//см. замечание в заголовке
			int index = _idxs + _is*Ns + _is;					//индекс в массиве посчитанных коэффициентов ФАР
																//для сорта атома ps и значения S. См. замечание в заголовке

			//int index = _idxs + (ps[i] - 1) * (Ns + 1);         
			
			IndScat = IndScat + KFAR[index];                   //слагаемое независимого рассеяния
		}
		
		//*****************************************************************************************
        //        Расчет интерференционного слагаемого : 
        //                 Na-1 Na              *       *        
        //        I(s) =2* SUM SUM 1/2 *[fi(s)*fj(s) + fi(s)*fj(s)]*
        //                 i=1 j=i+1             
		//
        //        *SIN(s*rij)/(s*rij)*EXP(-0.5*(drij*s)**2) 
		//*****************************************************************************************

		float InterfScat = 0.0f;                                //интерференционное слагаемое

		for(int i = 0; i < Nir; ++i)
		{
			int _is = SA1[i] - 1;
			int _js = SA2[i] - 1;

			int index = _idxs + _is*Ns + _js;                 //индекс в массиве посчитанных коэффициентов ФАР
			                                                  //для сортов атомов SA1 и SA2 и значения S
			float d = S * DispRij[i];

			InterfScat = InterfScat +                          //вычисляем интерференционное слагаемое 
				  KFAR[index] *
				  NRij[i] * sinf(S * Rij[i]) /
				  (S * Rij[i]) *
				  exp(- 0.5f * d * d );
		}

		float SumScat = IndScat + InterfScat * 2;             //суммарная интенсивность I(s)

		Io[k] = SumScat;                                       //I(s)

		Inorm[k] = SumScat / CHFE;                             //I(s) нормированная на единицу состава CHFE

		Hs[k] =                                                //H(s) интерференционная функция
			S * (Inorm[k] - IndScat / CHFE) * 
			Gf2_expSAlpha[k];                                  // Gf2_expSAlpha[k] =
		                                                       //    exp(- (S*Alpha)^2) / (GfactorSum2 / GfactorSum1)^2 ;
			

		S = S + dS;                                            //увеличиваем S на интервал dS
	}



}

//*****************************************************************************************
//
//      Расчет D(r)
//
//*****************************************************************************************
void calcDr(
	const float Drmax, const float delR, 
	const int NSma, const float Smin, const float dS, 
	const float R0e, 
	const float GfactorSum1,  
	const float *Hs,  
	int NRma,          
	float* &Dr           //возвращаемое значение
	)

{

	float FSL = 2.0f * PI * PI * R0e * GfactorSum1;

	float S;

	float R = 0;
	for(int i = 0; i < NRma; ++i)
	{
		S = Smin;
		float Sum = 0.0f;
		for(int k = 0; k < NSma; ++k)
		{
			Sum = Sum + Hs[k] * sinf(S * R);
			S = S + dS;
		}
		Dr[i] = FSL * R + Sum * dS;
		R = R + delR;
	}

}







