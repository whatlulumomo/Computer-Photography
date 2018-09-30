#include <opencv2/opencv.hpp>
#include <math.h>
using namespace cv;
// get shift value
void fftshift(const Mat &src, Mat &dst) {
	dst.create(src.size(), src.type());

	int rows = src.rows, cols = src.cols;
	Rect roiTopBand, roiBottomBand, roiLeftBand, roiRightBand;

	if (rows % 2 == 0) {
		roiTopBand = Rect(0, 0, cols, rows / 2);
		roiBottomBand = Rect(0, rows / 2, cols, rows / 2);
	}
	else {
		roiTopBand = Rect(0, 0, cols, rows / 2 + 1);
		roiBottomBand = Rect(0, rows / 2 + 1, cols, rows / 2);
	}

	if (cols % 2 == 0) {
		roiLeftBand = Rect(0, 0, cols / 2, rows);
		roiRightBand = Rect(cols / 2, 0, cols / 2, rows);
	}
	else {
		roiLeftBand = Rect(0, 0, cols / 2 + 1, rows);
		roiRightBand = Rect(cols / 2 + 1, 0, cols / 2, rows);
	}

	Mat srcTopBand = src(roiTopBand);
	Mat dstTopBand = dst(roiTopBand);
	Mat srcBottomBand = src(roiBottomBand);
	Mat dstBottomBand = dst(roiBottomBand);
	Mat srcLeftBand = src(roiLeftBand);
	Mat dstLeftBand = dst(roiLeftBand);
	Mat srcRightBand = src(roiRightBand);
	Mat dstRightBand = dst(roiRightBand);
	flip(srcTopBand, dstTopBand, 0);
	flip(srcBottomBand, dstBottomBand, 0);
	flip(dst, dst, 0);
	flip(srcLeftBand, dstLeftBand, 1);
	flip(srcRightBand, dstRightBand, 1);
	flip(dst, dst, 1);
}
int main() {
		// 傅里叶变换
		Mat I(512, 512, CV_32FC1);
		I(Rect(256 - 10, 256 - 30, 20, 60)) = 1.0;
	
		Mat J(I.size(), CV_32FC2);
		dft(I, J, DFT_COMPLEX_OUTPUT);
		fftshift(J, J);
	
		Mat Mag;
		vector<Mat> K;
		split(J, K);        // 将实数和虚数部分分解到 K[0] 和 K[1]
		pow(K[0], 2, K[0]); // 计算平方
		pow(K[1], 2, K[1]);
		Mag = K[0] + K[1];    // 两个分量的平方和
	
		Mat logMag;
		log(Mag + 1, logMag);
		normalize(logMag, logMag, 1.0, 0.0, CV_MINMAX);
		// ...
		imshow("Magnitude", logMag);
		waitKey(0); // 等待用户按下键盘
		destroyWindow("in"); // 销毁窗口 "hello"
		destroyWindow("out"); // 销毁窗口 "hello"
		return 0;
	
}