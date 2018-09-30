#include <opencv2/opencv.hpp>
#include <math.h>
using namespace cv;

// Compute gauss kernel
Mat GetGaussianKernel(const double sigma) {
	const double PI = acos(double(-1));
	int newsize = 2 * 5 * sigma + 1;

	Mat gaus = Mat(newsize, newsize, CV_64FC1); // 构造高斯核

	int center = newsize / 2;
	double sum = 0;

	// 根据正太分布式对kernel进行填充
	for (int i = 0; i < newsize; i++) {
		for (int j = 0; j < newsize; j++) {
			gaus.at<double>(i,j) = (1 / (2 * PI*sigma*sigma))*exp(-((i - center)*(i - center) + (j - center)*(j - center)) / (2 * sigma*sigma));
			sum += gaus.at<double>(i, j);
		}
	}

	// 进行归一化
	for (int i = 0; i < newsize; i++) {
		for (int j = 0; j < newsize; j++) {
			gaus.at<double>(i, j) /= sum;
		}
	}

	return gaus;
}

int main(int argc, char* argv[]) {
	if(argc != 4)
	{
		std::cout << "GaussianFilter <input-image> <output-image> <sigma>" << std::endl;
	}
	String inname = argv[1];
	String outname = argv[2];
	double sigma = atof(argv[3]);
	
	Mat image = imread(inname);
	Mat outimage = imread(inname);
	
	// compute filtely
	Mat gaus = GetGaussianKernel(sigma);
	filter2D(image,outimage,image.depth(),gaus);

	imshow("in", image); // show original image
	imshow("out", outimage); // show revised image
	imwrite(outname, outimage);

	waitKey(0); 
	destroyWindow("in"); // close window
	destroyWindow("out"); 

	return 0;
}