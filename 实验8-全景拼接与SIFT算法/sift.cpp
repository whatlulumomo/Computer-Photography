#include <opencv2\opencv.hpp>
#include <opencv2\highgui\highgui.hpp>
#include <opencv2\nonfree\features2d.hpp>
#include <opencv2\legacy\legacy.hpp>
#include <math.h>
#include <vector>
#include <string>
using namespace cv;
using namespace std;

class CylindricalPanorama {
public:
	virtual bool makePanorama(vector<Mat>& img_vec, Mat& img_out, double f) = 0;
};

class Oanorama2193: public CylindricalPanorama {
	double r = 500;
public:
	vector<Mat> image;
	void readImage(vector<string> name, double f);
	bool makePanorama(vector<Mat>& img_vec, Mat& img_out, double f);
	void imagemerge(Mat& img, Mat& img2);
	Mat ScalaLarge(Mat img, int rows, int cols);
	void displayAllImage();
	Mat project2cylinder(Mat image, double r, double f);
	void checkborder(Mat &img);
};

void Oanorama2193::checkborder(Mat &img) {
	int right;
	for (int j = 0; j < img.cols; j++) {
		bool flag = false;
		for (int i = 0; i < img.rows; i++) {
			if (norm(img.at<cv::Vec3b>(j, i)) != 0) {
				continue;
			}
			flag = true;
		}
		right = j + 1;
		if (flag == true) {
			break;
		}
	}

	Mat newimg = Mat(img.rows, right, CV_8UC3);
	for (int i = 0; i < newimg.rows; i++) {
		for (int j = 0; j < newimg.cols; j++) {
			newimg.at<Vec3b>(j,i) = img.at<Vec3b>(j, i);
		}
	}
	newimg.copyTo(img);
}

void Oanorama2193::displayAllImage() {
	for (int i = 0; i < image.size(); i++) {
		imshow(to_string(i), image[i]);
	}
}

void Oanorama2193::readImage(vector<string> name, double f) {
	Mat t = imread(name[0]);
	int COLS = t.cols;
	int ROWS = t.rows;

	for (int i = 0; i < name.size(); i++) {
		Mat s = imread(name[i]);
		s = project2cylinder(s, r, f);
		s = ScalaLarge(s, ROWS, 4 * COLS);
		image.push_back(s);
	}
}

bool Oanorama2193::makePanorama(vector<Mat>& image, Mat& img_out, double f) {
	int num = image.size();
	if (num <= 0) {
		return false;
	}
	int mid = num / 2;
	Mat img = image[mid];
	for (int i = 1; i < num; i++) {
		int index;
		if (i % 2 == 1) {
			index = mid - (i + 1) / 2;
			cout << index;
			Mat img2;
			image[index].copyTo(img2);
			imagemerge(img2, img);
			img2.copyTo(img);
		}
		else {
			index = mid + (i + 1) / 2;
			cout << index;
			Mat img2;
			image[index].copyTo(img2);
			imagemerge(img, img2);
		}
	}
	//checkborder(img_out);
	img.copyTo(img_out);
	return true;
}



Mat Oanorama2193::project2cylinder(Mat image, double r, double f) {
	int rows = image.rows;
	int cols = image.cols;
	int max_x = (cols - 1) / 2;
	int max_y = (rows - 1) / 2;
	int max_x_ = r * atan(max_x / f);
	int max_y_ = r * max_y / f;
	int max_cols = max_x_ * 2 + 1;
	int max_rows = max_y_ * 2 + 1;

	Mat M(max_rows, max_cols, CV_8UC3, cv::Scalar::all(1));
	for (int i = 0; i < max_rows; i++) {
		for (int j = 0; j < max_cols; j++) {

			double x = f * tan((j - max_x_) / r);
			double y = (i - max_y_) / r * sqrt(x*x + f * f);
			x = x + max_x;
			y = y + max_y;

			if (x < 0 || y < 0 || x >= cols || y >= rows) {
				M.at<cv::Vec3b>(i, j) = cv::Vec3b(0,0,0);
			}
			else {
				M.at<cv::Vec3b>(i, j) = image.at<cv::Vec3b>(y, x);
			}
			
		}
	}
	return M;
}

void Oanorama2193::imagemerge(Mat& img, Mat& img2) {
	//siftÌØÕ÷¼ì²â  
	SiftFeatureDetector  siftdtc;
	vector<KeyPoint>kp1, kp2;

	siftdtc.detect(img, kp1);
	Mat outimg1;
	drawKeypoints(img, kp1, outimg1);
	//imshow("image1 keypoints", outimg1);


	siftdtc.detect(img2, kp2);
	Mat outimg2;
	drawKeypoints(img2, kp2, outimg2);
	//imshow("image2 keypoints", outimg2);


	SiftDescriptorExtractor extractor;
	Mat descriptor1, descriptor2;
	BruteForceMatcher<L2<float>> matcher;
	vector<DMatch> matches;
	Mat img_matches;
	extractor.compute(img, kp1, descriptor1);
	extractor.compute(img2, kp2, descriptor2);

	matcher.match(descriptor1, descriptor2, matches);

	// Keep 200 best matches only.
	// We sort distance between descriptor matches
	Mat index;
	int nbMatch = int(matches.size());
	Mat tab(nbMatch, 1, CV_32F);
	for (int i = 0; i < nbMatch; i++)
		tab.at<float>(i, 0) = matches[i].distance;
	sortIdx(tab, index, SORT_EVERY_COLUMN + SORT_ASCENDING);
	vector<DMatch> bestMatches;

	for (int i = 0; i < 200; i++)
		bestMatches.push_back(matches[index.at < int >(i, 0)]);


	drawMatches(img, kp1, img2, kp2, bestMatches, img_matches);
	//imshow("matches", img_matches);

	std::vector<Point2f> dst_pts;                   //1st
	std::vector<Point2f> source_pts;                //2nd

	for (vector<DMatch>::iterator it = bestMatches.begin(); it != bestMatches.end(); ++it) {
		//-- Get the keypoints from the good matches
		dst_pts.push_back(kp1[it->queryIdx].pt);
		source_pts.push_back(kp2[it->trainIdx].pt);
	}

	Mat H = findHomography(source_pts, dst_pts, CV_RANSAC);

	//cout << H << endl;
	Mat wim_2;
	warpPerspective(img2, wim_2, H, img.size());


	// We can do this since im_1 and wim_2 have the same size.
	for (int i = 0; i < img.cols; i++)
		for (int j = 0; j < img.rows; j++) {
			Vec3b color_im1 = img.at<Vec3b>(Point(i, j));
			Vec3b color_im2 = wim_2.at<Vec3b>(Point(i, j));
			// ÐÂ
			if (norm(color_im1) < norm(color_im2))
				img.at<Vec3b>(Point(i, j)) = color_im2;
			// ¾É
			//if (norm(color_im1) == 0)
			//	img.at<Vec3b>(Point(i, j)) = color_im2;

		}
}

Mat  Oanorama2193::ScalaLarge(Mat img, int rows, int cols) {
	Mat newimg(3*rows, cols, CV_8UC3);
	for (int i = 0; i < img.rows; i++) {
		for (int j = 0; j < img.cols; j++) {
			newimg.at<Vec3b>(Point(j, i+rows)) = img.at<Vec3b>(Point(j, i));
		}
	}
	return newimg;
}


vector<string> data1 = {
	"./image/data1/DSC01599.jpg",
	"./image/data1/DSC01600.jpg",
	"./image/data1/DSC01601.jpg",
	"./image/data1/DSC01602.jpg",
	"./image/data1/DSC01603.jpg",
	"./image/data1/DSC01604.jpg",
	"./image/data1/DSC01605.jpg",
	"./image/data1/DSC01606.jpg",
	"./image/data1/DSC01607.jpg",
	"./image/data1/DSC01608.jpg",
	"./image/data1/DSC01609.jpg",
	"./image/data1/DSC01610.jpg", 
	"./image/data1/DSC01611.jpg",
	"./image/data1/DSC01612.jpg",
	"./image/data1/DSC01613.jpg",
	"./image/data1/DSC01614.jpg",
	"./image/data1/DSC01615.jpg",
	"./image/data1/DSC01616.jpg",
	"./image/data1/DSC01617.jpg",
	"./image/data1/DSC01618.jpg"
};

vector<string> data2 = {
	"./image/data2/DSC01538.jpg",
	"./image/data2/DSC01539.jpg",
	"./image/data2/DSC01540.jpg",
	"./image/data2/DSC01541.jpg",
	"./image/data2/DSC01542.jpg",
	"./image/data2/DSC01543.jpg",
	"./image/data2/DSC01544.jpg",
	"./image/data2/DSC01545.jpg",
	"./image/data2/DSC01546.jpg",
	"./image/data2/DSC01547.jpg",
	"./image/data2/DSC01548.jpg",
	"./image/data2/DSC01549.jpg",
};

int main() {
	double f = 512.89;
	Mat img;
	Oanorama2193 myparanoma;
	//myparanoma.readImage(data2, f);
	myparanoma.readImage(data1, f);
	//myparanoma.displayAllImage();
	myparanoma.makePanorama(myparanoma.image, img, f);


	imshow("paranoma_data.jpg", img);
	imwrite("paranoma_data_1.jpg", img);

	waitKey(0);


	
}
