
#ifndef PERIODICCRYSTAL1D
#define PERIODICCRYSTAL1D

#include <vector>
#include <Eigen>

class SpinWaveProblem1D {
public:
	SpinWaveProblem1D(int N, int baseSplit, int kStpes, int omegaSteps, double H1, double H2, double l1, double l2, double d, double omegaStart, double omegaEnd, double kStart, double kEnd, bool debug);
private:
	int N, baseSplit; //первое - количество гармоник в разложении, второе - мелкость разбиения при уточнении минимумов
	int kSteps, omegaSteps; //количество шагов по сетке ка-омега
	double omegaDelta; //мелкость по оси частот
	double l1, l2, d; //линейные размеры
	double gamma; //гиромагнитное отношения для циклической частоты
	double H1, H2, M4pi;  //параметры поля и намагниченности
	double omegaH1, omegaH2, omegaM; //всякие частоты
	double mu1(double omega); //коэффициенты в тензоре мю в разных полях
	double mu2(double omega);
	std::complex<double> M(int n, double omega); //для результатов разложения компоненты мю в ряд
	double b(int m); //блоховское волновое число

	void fillMatrix1(double k, double omega);
	Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces;
	void fixEigens();
	void refillMatrix1(double k);
	void goThroughGrid();
	std::complex<double> findDeterminant(double k, double omega);
	bool isMinimum(std::vector<std::complex<double> >& determinants, int i, int depth);

	int Debug = 0;
	double checkNull(double& k, double omega1, double omega2);
	double checkNull(double& k, double omega1, double omega2, double& startingDet, double& foundOmega);
	double omegaIdeal(double k, double omegaH, int mode);

	Eigen::MatrixXcd Matrix1;
	Eigen::MatrixXcd eigenVectors;
	Eigen::VectorXcd eigenValues;

	double omegaStart;
	double omegaEnd;
	double kStart;
	double kEnd;
	bool debug;
};

#endif // PERIODICCRYSTAL1D
