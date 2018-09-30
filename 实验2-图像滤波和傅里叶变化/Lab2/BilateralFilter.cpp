#include <opencv2/opencv.hpp>
#include <math.h>
using namespace cv;
#define Pi 3.1415926


// get gaussian value
double Gaussian(double sigma, double x)
{
	return exp(-pow(x, 2) / (2 * pow(sigma, 2))) / (sigma*pow(2 * Pi, 0.5));
}

// get Euclidean distance
double getDistance(int x1, int y1, int x2, int y2) {
	return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
}

// get pix according to pixel cooperation
Vec3b getPix(Mat image, int y, int x) {
	int ny = (y + image.rows) % image.rows;
	int nx = (x + image.cols) % image.cols;
	return image.at<Vec3b>(ny, nx);
}

// get wp
double GetWp(Mat image, int y, int x, double sigma1, double sigma2, int window, int k) {

	double wp = 0;
	int w = window / 2;

	for (int pi = y - w; pi < y + w + 1; pi++) {
		for (int pj = x - w; pj < x + w + 1; pj++) {
			double distance = getDistance(x, y, pj, pi);
			double lightdiff = abs(getPix(image, y, x)[k] - getPix(image, pi, pj)[k]);
			wp += Gaussian(sigma1, distance)*Gaussian(sigma2, lightdiff);
		}
	}

	return wp;
}

// get sum
double GetWpI(Mat image, int y, int x, double sigma1, double sigma2, int window, int k) {

	double wp = 0;
	int w = window / 2;

	for (int pi = y - w; pi < y + w + 1; pi++) {
		for (int pj = x - w; pj < x + w + 1; pj++) {
			double distance = getDistance(x, y, pj, pi);
			double lightdiff = abs(getPix(image, y, x)[k] - getPix(image, pi, pj)[k]);
			wp += Gaussian(sigma1, distance)*Gaussian(sigma2, lightdiff)*getPix(image, pi, pj)[k];
		}
	}

	return wp;
}



int main(int argc, char* argv[]) {

	if (argc != 5) {
		std::cout << "BilateralFilter <input-image> <output-image> <sigma-s> <sigma-r>" << std::endl;
		return 1;
	}

	String inname = argv[1];
	String outname = argv[2];
	Mat image = imread(inname);
	Mat outimage = imread(inname);

	int const window = 5;
	const int width = image.cols;
	const int height = image.rows;
	double sigma1 = atoi(argv[3]);
	double sigma2 = atoi(argv[4]);

	//cv::bilateralFilter(image, outimage, image.depth(), 20, 24);
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			for (int k = 0; k < 3; k++) {	// every pix need to process its neighboring info
				double wp = GetWp(image, i, j, sigma1, sigma2, window, k);
				double sum = GetWpI(image, i, j, sigma1, sigma2, window, k);
				outimage.at<Vec3b>(i, j)[k] = sum / wp;
			}
		}
	}

	//imwrite(outname, outimage);
	imshow("in", image); // show image
	imshow("out", outimage); 


	waitKey(0); 
	destroyWindow("in"); // close window
	destroyWindow("out");
	return 0;
}