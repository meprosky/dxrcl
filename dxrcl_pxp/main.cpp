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

	setlocale(LC_CTYPE, "");			//чтобы в консоли были русские буквы

	string comments1, comments2,	    //две строки заголовка вх. файла
		comments3, comments4;			//две строки для вых. файла

	int Na, Ns;							//число атомов, число сортов
	float lambdaXray;					//длина волны рентг. излучения

	int Nparc;							//число парциальных распределений пар атомов, выводимых 
										//в файлы для "ху". Если nparcr = 0, то не выводим и след. строка
										//массив parc должен отсутствовать

	int *parc;							//массив пар атомов

	float 
		Smin, Smax, dS, 
		Alpha, Drmax, delR, R0e, eps;  //Smin, Smax, dS, Al, Rmax, dR, R0e, eps

	float CHFE;							//число формульных единиц в кластере

	int   *ps;							//массив сортов атомов

	float *x, *y, *z;					//координаты атомов x, y, z

	int Nir;							//число интервалов по R
	int NSma;							//число интервалов по S
	float *Rij_0, *DispRij_0;			//"сырые" данные сумм R, квадартов
	int *NRij_0;						//кол-ва
	float *Rij, *DispRij;				//средние R, дисперсии
	int *NRij;							//кол-ва по R
	int *SA1, *SA2;						//вспомог. массивы индексов
	
	int Nreps, Npairs;					//кол-вл интервалов по eps(расстояний), пар перестановок

	float *DF1, *DF2;					 //первая и вторая дисперсионные поправки
	
	//9 коэффициентов на каждый сорт атома 
	//для расчета функций атомного рассеяния 
	float 
		*FA1, *FA2, *FA3, *FA4, 
		*FB1, *FB2, *FB3, *FB4, *FC, *FE;

	cout << endl;
	cout << "dxrcl_p9.exe" << endl;
	cout << "Расчет  распределения  интенсивности рассеяния  рентгеновских" << endl;
    cout << "лучей материалами, состоящими из беспорядочно-ориентированных" << endl; 
	cout << "кластеров. Петрозаводский  государственный  университет 2013." << endl;
	//cout << "Программа перепписана и оптимизирована на C++ Прохорским М.Е." << endl;
	//cout << "За основу взят dxrcl44.for на языке Fortran за авторством Фофанова А.Д." << endl;
	//cout << endl;
	cout << endl;

	clock_t t0 = clock();
	string curdt = DateTimeCur();
	ostringstream ostr;
	
	string inputfn;

	if(argc == 1)         //обработка строки запуска
	{
		std::cout << "Отсутсвует имя входного файла данных..." << endl;
		getchar();
		exit(1);
		//inputfn = "pdxr.xyz";
	}
	else
	{
		inputfn = argv[1];
	}
	
	string bufstr = inputfn.substr( inputfn.find_last_of(".") + 1);									//ищем расширение .dxr
	transform(bufstr.begin(), bufstr.end(), bufstr.begin(), ::tolower);								//переводим его в нижний регистр
	string outputfn = (bufstr == "dxr") ? inputfn.substr(0, inputfn.find_last_of(".")) : inputfn;   //если нашли, dxr отбрасываеми и получаем имя файла
	
	ofstream logf(outputfn + ".log");		//файл вывода лога работы

	bufstr = DateTimeCur();                 //текущие дата, время в текст. строку
	ostr.str("");
	ostr << "   Начало работы     - " + bufstr;
	LogTime(&logf, ostr.str(), 0, true); 

	//************************************************************************
	//
    //     Чтение входных данных
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
	
	comments3 = "Дата и время - " + bufstr;
	LogTime(&logf, "Данные " + inputfn + " прочитаны", clock() - t0, true); 
	
	Npairs = Ns * (Ns + 1) / 2;				// число перестановок сортов пар атомов
	
	//************************************************************************
    //     Вспомогательная матрица индексов для доступа к элементам массивов
	//************************************************************************
	int Ns1 = Ns + 1;
	int *idxA = new int[Ns1*Ns1];
	int *invidxSA1 = new int[Npairs];
	int *invidxSA2 = new int[Npairs];
	int cc = 0;
	for(int i = 1; i <= Ns; ++i)         //заполняем матрицу индексов для доступа к элементам массивов
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
    //     Ищем макс. диагональ кластера
	//************************************************************************
	float Rmax = findRmax(Na, x, y, z);

	ostr.str("");
	ostr << Rmax;
	LogTime(&logf, "Макс. расст. в кластере Rmax = " + ostr.str(), clock() - t0, true); 

	Nreps = (int)(Rmax / eps + 1.1);    //количество интервалов eps в макс. диагонали
	
	int Nall = Npairs * Nreps;          //длина массивов с учетом идентичности пар (1-2 = 2-1)
	                                    //например для 3 сортов атомов будет 1-1 1-2 1-3 2-2 2-3 3-3 всего 6 пар

	ostr.str("");
	ostr << "Старт выч. межатомных расстояний. Атомов = " << Na; 
	LogTime(&logf, ostr.str(), clock() - t0, true); 
	
	//************************************************************************
    //		САМЫЕ ЗАТРАТНЫЕ ПО ВРЕМНИ ВЫЧИСЛЕНИЯ (ибо зависят от от N(N-1)/2 )
	//		Вычисление совокупности межатомных расстояний в кластере
    //		с суммированием числа одинаковых (в пределах EPS) расстояний
    //		между одинаковыми парами атомов (A - B и B - A - одинаковы).
    //		Выч. среднее межат. расстояние и дисперсия этого расстояния
    //		по совокупности расстояний в пределах EPS.
    //************************************************************************
	Rij_0   =    new float[Nall];     //межатомные расстояния
	DispRij_0  = new float[Nall];     //дисперсии расстояний
	NRij_0  =    new   int[Nall];     //кол-во расстояний в интервале eps
	
	unsigned int sz_float = Nall * sizeof(float);
	unsigned int sz_int   = Nall * sizeof(int);
	
	memset(Rij_0,      0, sz_float);   //обнуляем
	memset(DispRij_0,  0, sz_float);
	memset(NRij_0,     0, sz_int);
	
	interatomic(
		Na, x, y, z, 
		Ns, Npairs, Nall, eps, ps, idxA,
		Rij_0, DispRij_0, NRij_0);


	ostr.str("");
	ostr << "Заверш. выч. межатомных расстояний. Циклов = " << Na*(Na-1)/2;
	LogTime(&logf, ostr.str(), clock() - t0, true); 

	//************************************************************************
	//	Подсчитываем число ненулевых элементов Nir
	//	т.е. тех межат. расстояния которых встречаются хотя бы один раз
	//	делается для выделения памяти для массивов точно под размер
	//************************************************************************
	Nir = 0; 
	for(int i = 0; i < Nreps; ++i){      //пробегаем по интервалам eps
		for(int j = 0; j < Npairs; ++j){ //пробегаем по комбинациям пар атомов 1-1 1-2 1-3 2-2 2-3 и т.д.
			if(NRij_0[i * Npairs + j]){  //не равно 0
				++Nir;                  //сколько ненулевых интервалов
			}
		}
	}
	

		
	Rij     = new float[Nir];     //межатомные расстояния
	DispRij = new float[Nir];     //дисперсии расстояний
	NRij    = new   int[Nir];     //кол-во расстояний в интервале eps
	SA1     = new   int[Nir];     //первый сорт атома в паре
	SA2     = new   int[Nir];     //второй сорт атома в паре

	//****************************************************************************
	//
	//   Вычисление средних расстояний, их количества в интервале, дисперсий 
	//   попутно "cжимаем" массив, т.е. исключаем элементы с NRij = 0 (нет расстояний)
	//
	//****************************************************************************
	calcNrAvrDisperse(
		Nir, Npairs, Nreps, invidxSA1, invidxSA2, 
		NRij_0, Rij_0, DispRij_0,
		SA1, SA2, 
		NRij, Rij, DispRij);
		
	delete [] invidxSA1; //освобождаем память
	delete [] invidxSA2;
	delete [] Rij_0;                                       
	delete [] NRij_0;
	delete [] DispRij_0;


	NSma = (int)((Smax - Smin) / dS + 1.1f);       //кол-во интервалов dS в диапазоне Smin-Smax

	//****************************************************************************
	//
	//   Вычисление значений ФАР в зависимости от S и не зависимых от
	//   числа атомов в кластере, от распеределений по расстояниям, дисперсиям и т.п.
	//
	//****************************************************************************
	float *Gf2_expSAlpha = new float[NSma];           //квадраты ф-ий обостряющего G-фактора в зав-ти от S память выделяется в функции где расчитывается
	
	//коэф-ты ФАР  в зав-ти от S
	//массив ФАР - суть 3х мерный массив, строки - значения S 
	//и соотв. ему 2х мерный массив ФАР i и j атома
	float *KFAR = new float[NSma * Ns * Ns];
	
	float GfactorSum1;

	calcFAR(Ns, NSma, Smin, dS, Alpha,
		DF1, DF2, FA1, FA2, FA3, FA4, FB1, FB2, FB3, FB4, FC,  FE,
		Gf2_expSAlpha,     //возвращаемое значение
		GfactorSum1,       //возвращаемое значение
		KFAR               //возвращаемое значение
	);

	ostr.str("");
	ostr << "Заверш. выч. ФАР. Размер массива = " << NSma << "*" << Ns << "*" << Ns << "=" << NSma*Ns*Ns;
	LogTime(&logf, ostr.str(), clock() - t0, true); 

	//*****************************************************************************************
	//
	//	 ОСНОВНОЙ БЛОК
	//   Вычисление функций I(s), H(s), D(r)
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
	//      Расчет D(r)
	//*****************************************************************************************
	int NRma = (int)(Drmax / delR + 1.1f); //кол-во интервалов для D(r)
	float *Dr = new float[NRma];       
	
	calcDr( Drmax, delR, NSma, Smin, dS, R0e, GfactorSum1, Hs, 
		NRma, Dr);

	ostr.str("");
	ostr << "Заверш. выч. I(s) H(s) D(r)";
	LogTime(&logf, ostr.str(), clock() - t0, true); 

	ostr.str("");
	ostr << (clock() - t0)/CLOCKS_PER_SEC;
	
	comments3 += ". Длит. расчетов,сек= " + ostr.str();

	//*****************************************************************************************
	//
	//      Запись значений 
	//		средних R, дисперсий, кол-в R в интервале eps,
	//      I(S), Inorm(S),H(S), D(r) в файлы
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
		
	//освобождаем память
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
	ostr << "   Завершение работы - " + DateTimeCur();
	LogTime(&logf, ostr.str(), clock() - t0, true); 
	//cout <<  "\nНажмите любую клавишу для выхода...";
	//getchar();
	return 0;
}
