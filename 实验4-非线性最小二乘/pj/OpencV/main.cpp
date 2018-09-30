#include <opencv2/opencv.hpp>
#include "hw3_gn.h"
using namespace cv;

int main(int argc, char* argv[]) {
	double *Coefficient = new double[3];
	for (int i = 0; i < 3; i++) {
		Coefficient[i] = 1;
	}
	ResidualFunction* s = new EclipseFunction();
	GaussNewtonParams param = GaussNewtonParams();
	GaussNewtonReport report;
	Sover2193* mysolver = new Sover2193();
	mysolver->solve(s, Coefficient, param, &report);
	std::cout << "Report:" << std::endl;
	std::cout << "Num of iteration: " << report.n_iter << std::endl;
	char co[3] = { 'A','B','C' };
	for (int i = 0; i < 3; i++) {
		std::cout << co[i] << " " << Coefficient[i] << std::endl;
	}
	switch (report.stop_type) {
	case report.STOP_RESIDUAL_TOL:
		std::cout << "Stop Cause: 余项达到阈值" << std::endl; break;
	case report.STOP_GRAD_TOL:
		std::cout << "Stop Cause: 梯度达到阈值" << std::endl; break;
	case report.STOP_NO_CONVERGE:
		std::cout << "Stop Cause: 不收敛" << std::endl; break;
	default:
		std::cout << "Stop Cause: 其它数值错误" << std::endl; break;
	}
	getchar();
}