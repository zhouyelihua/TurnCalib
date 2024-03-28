/*
  作者：叶於平
  时间：2014/11/01
  作用：利用线性回归拟合平面计算转轴的方向和转轴的空间位置
*/
#include<iostream>
#include<highgui.h>
#include<cv.h>
#include<iomanip>
#include<fstream>
#include<cmath>
#include<vector>
#include<gsl/gsl_sf_bessel.h>
#include<gsl\gsl_multimin.h>
/*
GSL Configuration Link: 
	https://blog.csdn.net/zhouyelihua/article/details/48019347?locationNum=8&fps=1
*/
using namespace std;
/////v3=crossmul(v1,v2)
void CrossMultiple(vector<double>v1, vector<double>v2, vector<double>&v3);

double calculate_angle(vector<double>v1, vector<double>v2);

void calulate_planerotation_angle(vector<vector<double>>plane, vector<double>&rotation_angle);

double vector_norm2(vector<double>v);

void calculate_axis_orientation(vector<vector<double>>plane, vector<double>&axis);

void calulate_distance(vector<vector<double>>plane, vector<double>axis);

void optim(vector<double>&Circle_P);
double my_f(const gsl_vector *v, void *params);
vector<vector<double>>all_plane;
int num_plane;

int main()
{
	vector<double> CirclePoint(3, 0);
	fstream output;
	output.open("..\\PlaneDataSet\\para.txt", std::fstream::out);
	vector<vector<double>>distance;
	vector<double>axis(3, 0);
	vector<double>angle;
	int num_plane;
	std::cout << "平面的个数" << endl;
	std::cin >> num_plane;
	vector<double> X_vector;
	vector<double> Y_vector;
	vector<double> Z_vector;
	vector<double> plane(4, 0);
	double x_temp, y_temp, z_temp, r_temp, g_temp, b_temp;
	fstream in;
	int ii = 0;
	output << "拟合的各个平面方程：\n";
		for (; ii < num_plane;ii++)
	{ 
		char filename[30] = {0};
		sprintf(filename, "..\\PlaneDataSet\\PointCloud_%d_0.asc", ii);
		/*fstream in(filename);*/
		
		in.open(filename, fstream::in);
		in.seekg(0);
		int NumPoint;
		//////////////////////////
		std::string line;
		while (getline(in, line))
		{
			if (line[0] != '#')
			{
				std::istringstream iss(line);
				iss >> x_temp >> y_temp >> z_temp >> r_temp >> g_temp >> b_temp;
				X_vector.push_back(x_temp);
				Y_vector.push_back(y_temp);
				Z_vector.push_back(z_temp);
			}
		}
		////////////////////////////
		/*while (!in.eof())
		{
			in >> x_temp >> y_temp >> z_temp >> r_temp >> g_temp >> b_temp;
			X_vector.push_back(x_temp);
			Y_vector.push_back(y_temp);
			Z_vector.push_back(z_temp);
		}*/
		NumPoint = X_vector.size();
		CvMat* X = cvCreateMat(NumPoint, 3, CV_64FC1);
		CvMat*XT = cvCreateMat(3, NumPoint, CV_64FC1);
		for (int i = 0; i < NumPoint; i++)
		{
			CV_MAT_ELEM(*X, double, i, 0) = 1;
			CV_MAT_ELEM(*X, double, i, 1) = X_vector[i];
			CV_MAT_ELEM(*X, double, i, 2) = Y_vector[i];
			CV_MAT_ELEM(*XT, double, 0, i) = 1;
			CV_MAT_ELEM(*XT, double, 1, i) = X_vector[i];
			CV_MAT_ELEM(*XT, double, 2, i) = Y_vector[i];

		}
		CvMat*b = cvCreateMat(NumPoint, 1, CV_64FC1);
		for (int i = 0; i < NumPoint; i++)
		{
			CV_MAT_ELEM(*b, double, i, 0) = Z_vector[i];

		}
		CvMat*XT_mul_X = cvCreateMat(3, 3, CV_64FC1);
		CvMat*XT_mul_b = cvCreateMat(3, 1, CV_64FC1);
		cvMatMul(XT, X, XT_mul_X);
		cvMatMul(XT, b, XT_mul_b);
		CvMat*Inv_XT_mul_X = cvCreateMat(3, 3, CV_64FC1);
		cvInvert(XT_mul_X, Inv_XT_mul_X);
		CvMat*ans = cvCreateMat(3, 1, CV_64FC1);
		cvMatMul(Inv_XT_mul_X, XT_mul_b, ans);
		plane[0] = CV_MAT_ELEM(*ans, double, 1, 0);
		plane[1] = CV_MAT_ELEM(*ans, double, 2, 0);
		plane[2] = -1;
		plane[3] = CV_MAT_ELEM(*ans, double, 0, 0);
		std::cout << "\n";
		std::cout << "A:" << CV_MAT_ELEM(*ans, double, 1, 0) << endl;
		std::cout << "B:" << CV_MAT_ELEM(*ans, double, 2, 0) << endl;
		std::cout << "C:-1" << endl;
		std::cout << "D:" << CV_MAT_ELEM(*ans, double, 0, 0) << endl;

		output << setw(17) << setprecision(17) << CV_MAT_ELEM(*ans, double, 1, 0) << "x+" << CV_MAT_ELEM(*ans, double, 2, 0) << "\ty-z+\t" << CV_MAT_ELEM(*ans, double, 0, 0) << "=0\n";
		all_plane.push_back(plane);
		X_vector.clear();
		Y_vector.clear();
		Z_vector.clear();

		in.clear();
		in.close();
	} //得到各个平面的方程
		if (all_plane.size() > 1)
			calculate_axis_orientation(all_plane, axis);
		else
		{
			fprintf(stderr, "error\n");
			exit(0);
		}
			

		calulate_planerotation_angle(all_plane, angle);
		output << "角度\n";
		double pi=acos(-1.0);
		for (auto ita = angle.begin(); ita != angle.end();ita++)
			output << setw(17) << setprecision(17) << *ita/pi*180 << "\t" ;
		output << "\n";
		//output << setw(15) << setprecision(12) << axis[0] << "\t" << axis[1] << "\t" << axis[2] << "\n";
		double distemp;
		vector<double> crosstemp(3, 0);
		vector<double>distancetemp(4,0);
		for (auto it = all_plane.begin(); it != all_plane.end(); it++)
		{
			CrossMultiple(*it, axis, crosstemp);
			distemp = vector_norm2(axis) / vector_norm2(crosstemp);
			distancetemp[0] = (*it)[0] * distemp;
			distancetemp[1] = (*it)[1] * distemp;
			distancetemp[2] = (*it)[2] * distemp;
			distancetemp[3] = (*it)[3] * distemp;
			distance.push_back(distancetemp);
		}
		int dis_size = distance.size();
		CvMat* A = cvCreateMat(dis_size, 3, CV_64FC1);
		CvMat* b1 = cvCreateMat(dis_size, 1, CV_64FC1);
		CvMat* circleP = cvCreateMat(3, 1, CV_64FC1);
		for (int i = 0; i < dis_size - 1;i++)
		{
			CV_MAT_ELEM(*A, double, i, 0) = distance[i][0]-distance[i+1][0];
			CV_MAT_ELEM(*A, double, i, 1) = distance[i][1] - distance[i + 1][1];
			CV_MAT_ELEM(*A, double, i, 2) = distance[i][2] - distance[i + 1][2];
			CV_MAT_ELEM(*b1, double, i, 0) = -(distance[i][3] - distance[i + 1][3]);
		}
		CV_MAT_ELEM(*A, double, dis_size-1, 0) = axis[0];
		CV_MAT_ELEM(*A, double, dis_size - 1, 1) = axis[1];
		CV_MAT_ELEM(*A, double, dis_size - 1, 2) = axis[2];
		CV_MAT_ELEM(*b1, double, dis_size - 1, 0) = 0;

		cvSolve(A, b1, circleP, CV_SVD);
		//CV_MAT_ELEM(*circleP, double, 0, 0) = -(CV_MAT_ELEM(*circleP, double, 1, 0)*axis[1] + CV_MAT_ELEM(*circleP, double, 2, 0)*axis[2]) / axis[0];
		std::cout << "转轴的空间位置\n";
		std::cout << setiosflags(ios::fixed) << setiosflags(ios::right) << setw(17) << setprecision(17) << CV_MAT_ELEM(*circleP, double, 0, 0) << "\t" << CV_MAT_ELEM(*circleP, double, 1, 0) << "\t" << CV_MAT_ELEM(*circleP, double, 2, 0) << endl;
		output << "转轴的方向:\n";
		output << setw(17) << setprecision(17) << axis[0] << "\t" << axis[1] << "\t" << axis[2] << "\n";
		output << "转轴的位置:\n";
		output << setw(17) << setprecision(17) << CV_MAT_ELEM(*circleP, double, 0, 0) << "\t" << CV_MAT_ELEM(*circleP, double, 1, 0) << "\t" << CV_MAT_ELEM(*circleP, double, 2, 0) << endl;
		CirclePoint.at(0) = CV_MAT_ELEM(*circleP, double, 0, 0);
		CirclePoint.at(1) = CV_MAT_ELEM(*circleP, double, 1, 0);
		CirclePoint.at(2) = CV_MAT_ELEM(*circleP, double, 2, 0);
		calulate_distance(all_plane, CirclePoint);
		optim(CirclePoint);
		calulate_distance(all_plane, CirclePoint);
		output << "优化之后的转轴位置:\n";
		output << setw(17) << setprecision(17) << CirclePoint.at(0) << "\t" << CirclePoint.at(1) << "\t" << CirclePoint.at(2) << endl;
		//*/
	//
	//cout << "A:"<<CV_MAT_ELEM(*ans, double, 1,0) << endl;
	//cout << "B:"<<CV_MAT_ELEM(*ans, double, 2, 0) << endl;
	//cout << "C:-1" << endl;
	//cout << "D:" << CV_MAT_ELEM(*ans, double, 0, 0) << endl;

}
void CrossMultiple(vector<double>v1, vector<double>v2, vector<double>&v3)
{
	v3[0] = v1[1] * v2[2] - v2[1] * v1[2];
	v3[1] = v2[0] * v1[2] - v1[0] * v2[2];
	v3[2] = v1[0] * v2[1] - v2[0] * v1[1];

	//vector<double>a = { 1, 2, 3, 4 };
	//vector<double>b = { 1, 4, 5, 4 };
	//vector<double>c(3, 0);
	//CrossMultiple(a, b, c);
	//cout << c[0] << c[1] << c[2] << endl;
}
double calculate_angle(vector<double>v1, vector<double>v2)
{
	double angle;
	double acos_value = (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]) / sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v2[2]) / sqrt(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]);
	if (acos_value < 0)
		acos_value = -acos_value;
	angle = acos(acos_value);
	//angle = angle / 3.1415926 / 2 * 360;//转化为角度制
	cout << "\nangle:" << angle / 3.1415926535626897 / 2 * 360 << endl;
	return angle;
}
void calulate_planerotation_angle(vector<vector<double>>plane, vector<double>&rotation_angle)
{
	double angle_temp;
	for (auto it1 = plane.begin(); it1 != plane.end() - 1; it1++)
	{
		angle_temp = calculate_angle(*it1, *(it1 + 1));
		rotation_angle.push_back(angle_temp);


	}


}
double vector_norm2(vector<double>v)
{
	return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}
void calculate_axis_orientation(vector<vector<double>>plane, vector<double>&axis)
{

	vector<vector<double>>res;
	vector<double>temp(3, 0);
	double sum_axis_x = 0;
	double sum_axis_y = 0;
	double sum_axis_z = 0;
	double norm2temp = 0;

	//////计算全部的平面交线
	for (auto it1 = plane.begin(); it1 != plane.end(); it1++)
		for (auto it2 = it1 + 1; it2 != plane.end(); it2++)
		{
		CrossMultiple(*it1, *it2, temp);
		norm2temp = vector_norm2(temp);
		temp[0] /= norm2temp;
		temp[1] /= norm2temp;
		temp[2] /= norm2temp;
		cout << setiosflags(ios::fixed) << setiosflags(ios::right) << setw(15) << setprecision(12) << temp[0] << "\t" << temp[1] << "\t" << temp[2] << endl;
		res.push_back(temp);
		}
	/*for (auto it1 = plane.begin(); it1 != plane.end()-1; it1++)
	{
	CrossMultiple(*it1, *(it1+1), temp);
	res.push_back(temp);
	}
	CrossMultiple(*(plane.end()), *(plane.begin()), temp);
	res.push_back(temp);*/
	int n = res.size();
	for (auto it3 = res.begin(); it3 != res.end(); it3++)
	{
		if ((*it3)[1] < 0)
		{
			(*it3)[0] = -(*it3)[0];
			(*it3)[1] = -(*it3)[1];
			(*it3)[2] = -(*it3)[2];

		}
		sum_axis_x += (*it3)[0];

		sum_axis_y += (*it3)[1];
		sum_axis_z += (*it3)[2];
	}
	axis[0] = sum_axis_x / n;
	axis[1] = sum_axis_y / n;
	axis[2] = sum_axis_z / n;
	//归一化
	double normalization;
	normalization = sqrt(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]);
	axis[0] /= normalization;
	axis[1] /= normalization;
	axis[2] /= normalization;
	cout << "交线的平均值（转轴的方向）\n";

	cout << setiosflags(ios::fixed) << setiosflags(ios::right) << setw(15) << setprecision(12) << axis[0] << "\t" << axis[1] << "\t" << axis[2] << endl;
}
void calulate_distance(vector<vector<double>>plane, vector<double>axis)
{
	double distance_temp_y;
	double temp_temp;
	double A, B, C, D;
	cout << "distance:";
	for (auto it = plane.begin(); it < plane.end(); ++it)
	{
		A = (*it).at(0);
		B = (*it).at(1);
		C = (*it).at(2);
		D = (*it).at(3);
		//cout << A << "\t" << B << "\t" << C << "\t" << D << endl;
		temp_temp = axis.at(0)*A + axis.at(1)*B + axis.at(2)*C + D;
		if (temp_temp < 0)
			temp_temp = -temp_temp;
		distance_temp_y = temp_temp / vector_norm2(*it);
		cout << distance_temp_y << "\t";
	}
	cout << endl;
}
double my_f(const gsl_vector *v, void *params)
{
	//double x, y;
	//double *dp = (double *)params;
	double axis_x, axis_y, axis_z;
	axis_x = gsl_vector_get(v, 0);
	axis_y = gsl_vector_get(v, 1);
	axis_z = gsl_vector_get(v, 2);
	double err = 0;
	/////
	vector<double>distance_temp_yyp;
	double yyp_temp;
	double temp_temp;
	double A, B, C, D;
	//cout << "distance:";
	for (auto it = all_plane.begin(); it < all_plane.end(); ++it)
	{
		//cout << "ssssssssssssssssssssssssssd" << endl;
		A = (*it).at(0);
		B = (*it).at(1);
		C = (*it).at(2);
		D = (*it).at(3);
		//cout << A << "\t" << B << "\t" << C << "\t" << D << endl;
		temp_temp = axis_x*A + axis_y*B + axis_z*C + D;
		if (temp_temp < 0)
			temp_temp = -temp_temp;
		yyp_temp = temp_temp / vector_norm2(*it);
		//	cout << "yyp_temp=" << yyp_temp << endl;
		distance_temp_yyp.push_back(yyp_temp);
		//cout << distance_temp_y << "\t";
	}
	int vec_size = distance_temp_yyp.size();
	for (int i = 0; i < vec_size - 1; ++i)
	{
		err += (distance_temp_yyp.at(i) - distance_temp_yyp.at(i + 1))*(distance_temp_yyp.at(i) - distance_temp_yyp.at(i + 1));
	}
	//err += (distance_temp_yyp.at(vec_size - 1) - distance_temp_yyp.at(0))*(distance_temp_yyp.at(vec_size - 1) - distance_temp_yyp.at(0));
	//cout << endl;
	//////
	cout << "err=" << err << endl;
	return err;
}
void optim(vector<double>&Circle_P)
{

	cout << "waiting...\n\n";
	size_t np = 3;
	double par[2] = { 1, 2 };
	const gsl_multimin_fminimizer_type *T =
		gsl_multimin_fminimizer_nmsimplex;
	gsl_multimin_fminimizer *s = NULL;
	gsl_vector *ss, *x;
	gsl_multimin_function minex_func;

	size_t iter = 0, i;
	int status;
	double size;

	/* Initial vertex size vector */
	ss = gsl_vector_alloc(np);

	/* Set all step sizes to 1 */
	gsl_vector_set_all(ss, 0.0000001);

	/* Starting point */
	x = gsl_vector_alloc(np);

	gsl_vector_set(x, 0, Circle_P.at(0));
	gsl_vector_set(x, 1, Circle_P.at(1));
	gsl_vector_set(x, 2, Circle_P.at(2));

	/* Initialize method and iterate */
	minex_func.f = &my_f;
	minex_func.n = np;
	minex_func.params = (void *)&par;

	s = gsl_multimin_fminimizer_alloc(T, np);
	gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

	do
	{
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);

		if (status)
			break;

		size = gsl_multimin_fminimizer_size(s);
		status = gsl_multimin_test_size(size, 1e-12);

		if (status == GSL_SUCCESS)
		{
			printf("converged to minimum at\n");
			printf("%5d ", iter);
			printf("f() = %12.11f size = %.3f\n", s->fval, size);
		}

		//	printf("%5d ", iter);
		for (i = 0; i < np; i++)
		{
			;//printf("%10.3e ", gsl_vector_get(s->x, i));
		}
		//printf("f() = %12.11f size = %.3f\n", s->fval, size);
	} while (status == GSL_CONTINUE && iter < 1000);


	///

	Circle_P.at(0) = gsl_vector_get(s->x, 0);
	Circle_P.at(1) = gsl_vector_get(s->x, 1);
	Circle_P.at(2) = gsl_vector_get(s->x, 2);
	///
	//cout << "sssssssssssssssssssssssssssssssssssssssssss\n";
	//calulate_distance(all_plane, Circle_P);
	gsl_vector_free(x);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free(s);

	//return status;


}