#define _CRT_SECURE_NO_WARNINGS
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define PI 3.1415926

void main()
{
	while (1)
	{
		//正向转换XYZ->色貌属性
		double  X = 0;
		double  Y = 0;
		double  Z = 0;  //初始化色度坐标XYZ

		double  Xw = 94.81000;//测试光源下的白点的三刺激值
		double  Yw = 100.0;//在CIECAM02模型中参考光源的白点三刺激值均为100，这在模型中是固定的
		double  Zw = 107.330;  //即Xwr=Ywr=Zwr=100

		printf("请输入X、Y、Z:\n");
		scanf("%lf %lf %lf", &X, &Y, &Z);   //输入XYZ的值
		double LA = 20.0;//适应场的亮度LA(cd/m2)，通常取实验照明条件下完全漫反射平面亮度LW的20%

		double F = 0;               //明度对比因子：F
		double c = 0;               //环境响应参数：c
		double Nc = 0;            //色诱导因子：Nc
		F = 1.0; c = 0.69; Nc = 1.0;

		//--------------------- 正向模型计算-------------------------\\
					//1. 由X、Y、Z计算R、G、B-------------------------------------------
		double   xyz[3][1] = { { X },{ Y },{ Z } };// 定义XYZ数组，作为矩阵
		double Mcat[3][3] = { { 0.401288,0.650173, -0.051461 },{-0.250268,1.204414,0.045854 },{ -0.002079,0.048952,0.953127} };//输入转换矩阵M16
		double XwYwZw[3][1] = { { Xw },{ Yw },{ Zw } };
		int i, j, k;
		double temp[1] = { 0.0 };
		double RGB[3][1] = { 0, 0, 0 };
		for (i = 0; i < 3; i++)
		{
			for (k = 0; k < 1; k++)
				temp[k] = 0;
			for (j = 0; j < 3; j++)
			{//当前行的每个元素
				for (k = 0; k < 1; k++) {
					temp[k] += Mcat[i][j] * xyz[j][k];
				}
			}
			for (k = 0; k < 1; k++)
			{
				RGB[i][k] = temp[k];
				//printf("\n%.2lf\t", RGB[i][k]);
			}
			//printf("\n");
		}
		double R = RGB[0][0]; double G = RGB[1][0];   double B = RGB[2][0];
		//利用XwYwZw求RwGwBw
		int i1, j1, k1;
		double temp1[1] = { 0.0 };
		double RwGwBw[3][1] = { { 0.0 },{ 0.0 },{ 0.0 } };
		for (i1 = 0; i1 < 3; i1++)
		{
			for (k1 = 0; k1 < 1; k1++)
				temp1[k1] = 0;
			for (j1 = 0; j1 < 3; j1++)
			{//当前行的每个元素
				for (k1 = 0; k1 < 1; k1++) {
					temp1[k1] += Mcat[i1][j1] * XwYwZw[j1][k1];
				}
			}
			for (k1 = 0; k1 < 1; k1++)
			{
				RwGwBw[i1][k1] = temp1[k1];
				//printf("\n%lf\t", RwGwBw[i1][k1]);
			}
			//printf("\n");
		}

		double Rw = RwGwBw[0][0]; double Gw = RwGwBw[1][0];   double Bw = RwGwBw[2][0];

		//2.计算色适应因子D----------------------------------------------------------
		double   D = F*(1 - (1.0 / 3.6)*exp((-LA - 42) / 92));
		if (D > 1) D = 1;
		else if (D < 0)  D = 0;
		else;
		//3.由 R,G,B 计算在参考条件下样品对应色的三刺激值------------------------------
		double DR = D*(Yw / Rw) + 1 - D;
		double DG = D*(Yw / Gw) + 1 - D;
		double DB = D*(Yw / Bw) + 1 - D;
		//由R,G,B 计算在参考条件下样品对应色的三刺激值（进行了色适应转换）
		double  Rc = DR*R;
		double  Gc = DG*G;
		double  Bc = DB*B;
		//double RcGcBc = { { Rc }, { Gc }, { Bc } };
		double  Rwc = DR*Rw;
		double  Gwc = DG*Gw;
		double  Bwc = DB*Bw;
		//4.计算亮度水平适应因子FL---------------------------------------------------------
		double   kk = 1.0 / (5 * LA + 1);
		double FL = 0.2*pow(kk, 4)*(5 * LA) + 0.1*pow((1 - pow(kk, 4)), 2)*pow((5 * LA), (1.0 / 3.0));
		//5.计算诱导因子Nbb和Ncb---------------------------------------------------------
		double Lw = Yw;
		double  Yb = 20.0;       //Yb的值一般取20                     //double  Yb = 100.0 * Lb / Lw;
		double n = Yb / Yw;
		double z = 1.48 + sqrt(n);  // 计算非线性底数幂
		double Nbb = 0.725*pow((1.0 / n), 0.2);     double Ncb = Nbb;
		//6. 变换到HPE(Hunt-Pointer-Estevez)空间-----------------
		double Mcat_1[3][3] = { { 1.096124, -0.278869, 0.182745 },{ 0.454369, 0.473533, 0.072098 },{ -0.009628, -0.005698, 1.015326 } };
		double Mhpe[3][3] = { { 0.38971, 0.68898, -0.07868 },{ -0.22981, 1.18340, 0.04641 },{ 0.00000, 0.00000, 1.00000 } };
		double R_zhong = Mcat_1[0][0] * Rc + Mcat_1[0][1] * Gc + Mcat_1[0][2] * Bc;
		double G_zhong = Mcat_1[1][0] * Rc + Mcat_1[1][1] * Gc + Mcat_1[1][2] * Bc;
		double B_zhong = Mcat_1[2][0] * Rc + Mcat_1[2][1] * Gc + Mcat_1[2][2] * Bc;
		//(R' G' B') =MHPE *M-1CAT02(Rc Gc Bc)
		double R_ = Mhpe[0][0] * R_zhong + Mhpe[0][1] * G_zhong + Mhpe[0][2] * B_zhong;
		double G_ = Mhpe[1][0] * R_zhong + Mhpe[1][1] * G_zhong + Mhpe[1][2] * B_zhong;
		double B_ = Mhpe[2][0] * R_zhong + Mhpe[2][1] * G_zhong + Mhpe[2][2] * B_zhong;
		//在此用R_ G_ B_代表 R' G' B'

		double Rzhong = Mcat_1[0][0] * Rwc + Mcat_1[0][1] * Gwc + Mcat_1[0][2] * Bwc;
		double Gzhong = Mcat_1[1][0] * Rwc + Mcat_1[1][1] * Gwc + Mcat_1[1][2] * Bwc;
		double Bzhong = Mcat_1[2][0] * Rwc + Mcat_1[2][1] * Gwc + Mcat_1[2][2] * Bwc;
		//(Rw' Gw' Bw') =MHPE *M-1CAT02(Rwc Gwc Bwc)
		double Rw_ = Mhpe[0][0] * Rzhong + Mhpe[0][1] * Gzhong + Mhpe[0][2] * Bzhong;
		double Gw_ = Mhpe[1][0] * Rzhong + Mhpe[1][1] * Gzhong + Mhpe[1][2] * Bzhong;
		double Bw_ = Mhpe[2][0] * Rzhong + Mhpe[2][1] * Gzhong + Mhpe[2][2] * Bzhong;
		//在此用Rw_ Gw_ Bw_代表 Rw' Gw' Bw'
		double Raw_ = 400 * pow(((FL*Rw_) / 100.0), 0.42) / (pow(((FL*Rw_) / 100.0), 0.42) + 27.13) + 0.1;
		double Gaw_ = 400 * pow(((FL*Gw_) / 100.0), 0.42) / (pow(((FL*Gw_) / 100.0), 0.42) + 27.13) + 0.1;
		double Baw_ = 400 * pow(((FL*Bw_) / 100.0), 0.42) / (pow(((FL*Bw_) / 100.0), 0.42) + 27.13) + 0.1;

		double Aw = (2 * Raw_ + Gaw_ + Baw_ / 20.0 - 0.305)*Nbb;     //消色差响应Aw

																	 //7. 计算更接近于视锥体反应的基本值（非线性压缩）
		double Ra_ = 400 * pow(((FL*R_) / 100.0), 0.42) / (pow(((FL*R_) / 100.0), 0.42) + 27.13) + 0.1;
		double Ga_ = 400 * pow(((FL*G_) / 100.0), 0.42) / (pow(((FL*G_) / 100.0), 0.42) + 27.13) + 0.1;
		double Ba_ = 400 * pow(((FL*B_) / 100.0), 0.42) / (pow(((FL*B_) / 100.0), 0.42) + 27.13) + 0.1;
		if (Ra_ < 0)   Ra_ = -400 * pow(((-FL*R_) / 100.0), 0.42) / (pow(((-FL*R_) / 100.0), 0.42) + 27.13) + 0.1;
		else;
		if (Ga_ < 0)  Ga_ = -400 * pow(((-FL*G_) / 100.0), 0.42) / (pow(((-FL*G_) / 100.0), 0.42) + 27.13) + 0.1;
		else;
		if (Ba_ < 0)   Ba_ = -400 * pow(((-FL*B_) / 100.0), 0.42) / (pow(((-FL*B_) / 100.0), 0.42) + 27.13) + 0.1;
		else;
		//用 Ra_ Ga_ Ba_代表Ra' Ga' Ba' 


		//8. 计算红绿信息a、黄蓝信息b、色调角h-------------------------------------------
		double a = Ra_ - (12.0*Ga_) / 11.0 + Ba_ / 11.0;     //红绿信息a
		double b = (Ra_ + Ga_ - 2.0*Ba_) / 9.0;                  //黄蓝信息b
		double h = atan(b / a) * 180 / PI + 180;                                               //色调角h
		if (h > 360 || h < 0)   abort();
		else;

		//9.计算离心率 et 和色相构成H----------------------------------------------------------------
		//Unique hue data for calculation of hue quadrature:
		double h1 = 20.14; double e1 = 0.8; double H1 = 0.0;         //Red      i=1
		double h2 = 90.00; double e2 = 0.7; double H2 = 100.0;     //Yellow i=2
		double h3 = 164.25; double e3 = 1.0; double H3 = 200.0;   //Green   i=3
		double h4 = 237.53; double e4 = 1.2; double H4 = 300.0;   //Blue     i=4
		double h5 = 380.14; double e5 = 0.8; double H5 = 400.0;    //Red      i=5
		double h_;     double H;   //色相构成H
		if (h < h1)  h_ = h + 360;
		else h_ = h;
		double et = 0.25*(cos(h_*PI / 180.0 + 2) + 3.8);   //离心率 et
		if (h1 <= h_ &&h_ <= h2)         H = H1 + (100 * (h_ - h1) / e1) / ((h_ - h1) / e1 + (h2 - h_) / e2);
		else if (h2 <= h_&&h_ <= h3)  H = H2 + (100 * (h_ - h2) / e2) / ((h_ - h2) / e2 + (h3 - h_) / e3);
		else if (h3 <= h_&&h_ <= h4)  H = H3 + (100 * (h_ - h3) / e3) / ((h_ - h3) / e3 + (h4 - h_) / e4);
		else if (h4 <= h_ &&h_ <= h5)   H = H4 + (100 * (h_ - h4) / e4) / ((h_ - h4) / e4 + (h5 - h_) / e5);
		else printf("H溢出，出现错误\n");

		//10.计算消色差响应A和明度J及其他参数--------------------------------------------------------
		double A = (2 * Ra_ + Ga_ + Ba_ / 20.0 - 0.305)*Nbb;      //消色差响应A
		double J = 100 * pow((A / Aw), c*z);                                  //心理明度J
		double Q = (4.0 / c)*pow((J / 100.0), 0.5)*(Aw + 4.0)*pow(FL, 0.25);   //绝对亮度
		double t = (50000.0 / 13.0*Nc*Ncb)*et*pow((pow(a, 2) + pow(b, 2)), 0.5) / (Ra_ + Ga_ + 21.0 / 20.0*Ba_);
		double C = pow(t, 0.9)*pow((J / 100.0), 0.5)*pow((1.64 - pow(0.29, n)), 0.73);    //彩度C
		double M = C*pow(FL, 0.25);      //绝对彩度M
		double  s = 100 * pow((M / Q), 0.5);   //饱和度s

											   //	-----------------------------------------------------CIECAM02-UCS--------------------------------------------------------------------------------
											   // CIECAM02 色貌模型提出了色饱和度（satuation）、视彩度（colorfulness）和视亮度（brightness）属性，并在明度和视亮度部分加入了预测因子以避免出现负值，并且还考虑了由不同周边环境产生的影响。
											   /*利用上面的颜色属性可以组成一定形式的颜色空间，Luo等人在2006年使用三种色差数据：日光照明条件大色差和小色差数据，以及A
											   光源下的小色差数据建立了新的基于CIECAM02色貌模型的均匀颜色空间以及色差公式CAM02 - LCD、CAM02 - SCD、
											   CAM02 - UCS，分别适用于大色差数据、小色差数据以及大、小色差联合数据，在这些色差数据中，大色差的均值为 10 个CIELAB
											   色差单位，小色差的均值为 2.5 个CIELAB色差单位。此系列均匀颜色空间是以J'，a'，b'为三维坐标*/

		double KL = 1.00;
		double c1 = 0.007;
		double c2 = 0.0228;    //CIECAM02-UCS色差公式参数设置
		double J_ = (1 + 100 * c1)*J / (1 + c1*J);
		double M_ = (1 / c2)*log(1 + c2*M);
		double a_ = M_*cos(h*PI / 180);
		double b_ = M_*sin(h*PI / 180);   //J'a'b'空间的建立 


		/*printf("\n\nCIECAM02-UCS色貌模型参数：\n");
		printf("\n心理明度J‘:%lf\n", J_);
		printf("	a':%lf\n", a_);
		printf("	b':%lf\n", b_);
		printf("	C':%lf\n", C);
		printf("	h':%lf\n", h);

		char mystr[100] = { 0 };
		sprintf(mystr, "echo %lf   %lf    %lf   %lf >> E:\\1.txt", J_, a_, b_, C);
		system(mystr);*/
	}
	system("pause");
}