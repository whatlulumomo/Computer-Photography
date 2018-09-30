//#include <opencv2/opencv.hpp>
//#include <math.h>
//using namespace cv;
//#define Pi 3.1415926
//
//int main(int argc, char* argv[]) {
//	if (argc != 5) {
//		std::cout << "BilateralFilter <input-image> <output-image> <sigma-s> <sigma-r>" << std::endl;
//		return 1;
//	}
//
//
//	/*Mat image = imread("Emma-in.png");
//	Mat outimage = imread("Emma-in.png");*/
//	
//	String inname = argv[1];
//	String outname = argv[2];
//	Mat image = imread(inname);
//	Mat outimage = imread(outname); 
//
//	int w = 3;
//	int h = 3;
//	for (int i = 0; i < image.rows; i++) {
//		for (int j = 0; j < image.cols; j++) {
//			int sx = j - w + 1;
//			int sy = i - h + 1;
//			int ex = j + w - 1;
//			int ey = i + h - 1;
//			vector<int> testx;
//			vector<int> testy;
//			vector<int> testz;
//			Vec3b pix;
//			for (int ti = sy; ti <= ey; ti++) {
//				for (int tj = sx; tj <= ex; tj++) {
//					pix = image.at<Vec3b>((ti + image.rows)%image.rows, (tj+ image.cols)%image.cols);
//					testx.push_back(pix[0]);
//					testy.push_back(pix[1]);
//					testz.push_back(pix[2]);
//				}
//			}
//			sort(testx.begin(), testx.end());
//			sort(testy.begin(), testy.end());
//			sort(testz.begin(), testz.end());
//			int pos = testx.size() / 2;
//			int mx = testx[pos];
//			int my = testy[pos];
//			int mz = testz[pos];
//			outimage.at<Vec3b>(i, j) = Vec3b(mx,my,mz);
//		}
//	}
//
//	imwrite(outname, outimage);
//	imshow("in", image); // 在窗口 "hello" 中显示图片
//	imshow("out", outimage); // 在窗口 "hello" 中显示图片
//
//
//	waitKey(0); // 等待用户按下键盘
//	destroyWindow("in"); // 销毁窗口 "hello"
//	destroyWindow("out"); // 销毁窗口 "hello"
//	return 0;
//}