#include "TMath.h"
#include "TGraph.h"
#include "TF1.h"


Double_t koshixy(Double_t t_max = 75, Double_t vect_v0 = 1, Double_t pend_L = 5, Double_t coord_x0 = 3, Double_t coord_y0 = -4){

	Int_t numb_N = 7500;
	Double_t dev_g = 0.05;
	
	Double_t t_min = 0;
	Double_t time_h = (t_max - t_min)/numb_N;
//	Double_t epsilon = 1.E-8;
// max dev x and y at init cond for vect_v0 = 0 and t_max = 3
// N=300:  j=231  t=2.31    x=-3.07594 y=-3.9419
// N=600:  j=461  t=2.305   x=-3.038   y=-3.97122
// N=900:  j=691  t=2.30333 x=-3.02548 y=-3.98076
// N=1200: j=921  t=2.3025  x=-3.01924 y=-3.9855
// N=1500: j=1152 t=2.304   x=-3.0155  y=-3.98833
	Int_t check = -1;

	TGraph* graph_Xt = new TGraph(numb_N + 1);
	TGraph* graph_Yt = new TGraph(numb_N + 1);
	TGraph* graph_X2Y2 = new TGraph(numb_N + 1);	
	for (int i = 0; i <= numb_N; i++){
		graph_Xt->SetPoint(i, time_h*i, 0);
		graph_Yt->SetPoint(i, time_h*i, 0);
		graph_X2Y2->SetPoint(i, time_h*i, 0);
	}

	graph_Xt->SetPoint(0, 0, coord_x0);
	graph_Yt->SetPoint(0, 0, coord_y0);
	graph_X2Y2->SetPoint(0, 0, graph_Xt->GetPointY(0)*graph_Xt->GetPointY(0) + graph_Yt->GetPointY(0)*graph_Yt->GetPointY(0));
	graph_Xt->SetPointY(1, time_h*vect_v0*coord_y0/pend_L + coord_x0);
	graph_Yt->SetPointY(1, -time_h*vect_v0*coord_x0/pend_L + coord_y0);
	graph_X2Y2->SetPointY(1, graph_Xt->GetPointY(1)*graph_Xt->GetPointY(1) + graph_Yt->GetPointY(1)*graph_Yt->GetPointY(1));

	Double_t val_k = 0;
//	Double_t val_m = 0, val_l = 0;
	Double_t coef_Ax = 1, coef_Bx = -1, coef_Cx = 0;
	Double_t coef_Ay = 1, coef_By = -1, coef_Cy = 0;
	Double_t val_x0 = 0, val_y0 = 0, val_x1 = 0;
	Double_t dis_x = 0, dis_y = 0;

	TF1* func_Xt = new TF1("func Xt","[0]*x*x + [1]*x + [2]", -pend_L - 1, pend_L + 1);
	TF1* func_Yt = new TF1("func Yt","[0]*x*x + [1]*x + [2]", -pend_L - 1, pend_L + 1);
	func_Xt->SetParameters(coef_Ax, coef_Bx, coef_Cx);
	func_Yt->SetParameters(coef_Ay, coef_By, coef_Cy);

	for (int j = 0; j < numb_N - 1; j++){
		val_k = 2*graph_Yt->GetPointY(j + 1) - time_h*time_h*(9.81 + dev_g*sin(2*TMath::Pi()*time_h*j));
/*	
////////////////////////////////////////////////////////////////////////////////////
		val_m = 2*graph_Yt->GetPointY(j)*graph_Xt->GetPointY(j+1)/graph_Xt->GetPointY(j) - val_k;
		val_l = graph_Xt->GetPointY(j)*val_k/graph_Yt->GetPointY(j) - 2*graph_Xt->GetPointY(j+1);

		coef_Ax = graph_Yt->GetPointY(j)*graph_Yt->GetPointY(j)/graph_Xt->GetPointY(j)/graph_Xt->GetPointY(j) + 1;
		coef_Bx = -2*graph_Yt->GetPointY(j)*val_m/graph_Xt->GetPointY(j);
		coef_Cx = val_m*val_m - pend_L*pend_L;

		coef_Ay = graph_Xt->GetPointY(j)*graph_Xt->GetPointY(j)/graph_Yt->GetPointY(j)/graph_Yt->GetPointY(j) + 1;
		coef_By = -2*graph_Xt->GetPointY(j)*val_l/graph_Yt->GetPointY(j);
		coef_Cy = val_l*val_l - pend_L*pend_L;	
////////////////////////////////////////////////////////////////////////////////////
*/

////////////////////////////////////////////////////////////////////////////////////
		val_x0 = graph_Xt->GetPointY(j);
		val_x1 = graph_Xt->GetPointY(j + 1);
		val_y0 = graph_Yt->GetPointY(j);	
	
		coef_Ax = pend_L*pend_L;
		coef_Bx = 2*val_x0*val_y0*val_k - 4*val_y0*val_y0*val_x1;
		coef_Cx = 4*val_y0*val_y0*val_x1*val_x1 + val_x0*val_x0*val_k*val_k - 4*val_y0*val_k*val_x0*val_x1 - val_x0*val_x0*pend_L*pend_L;
	
		coef_Ay = pend_L*pend_L;
		coef_By = 4*val_x0*val_x1*val_y0 - 2*val_x0*val_x0*val_k;
		coef_Cy = val_x0*val_x0*val_k*val_k + 4*val_x1*val_x1*val_y0*val_y0 - 4*val_y0*val_k*val_x0*val_x1 - val_y0*val_y0*pend_L*pend_L;
/////////////////////////////////////////////////////////////////////////////////////
	
//		func_Xt->SetParameters(coef_Ax, coef_Bx, coef_Cx);
//		func_Yt->SetParameters(coef_Ay, coef_By, coef_Cy);

		if ( j==check ){
			func_Xt->SetParameters(coef_Ax, coef_Bx, coef_Cx);
			func_Yt->SetParameters(coef_Ay, coef_By, coef_Cy);
		
			func_Yt->Draw();
			func_Xt->Draw("same");
		}	

/*
//////////////////////////////////////////////////////////////////////////////////////	
		if (graph_Xt->GetPointY(j + 1) > 0){
			if(2*graph_Xt->GetPointY(j + 1) - graph_Xt->GetPointY(j) > 0){
				graph_Xt->SetPointY(j + 2, func_Xt->GetX(0, 0, pend_L + 0.1, epsilon));
				graph_Yt->SetPointY(j + 2, func_Yt->GetX(0, -pend_L - 0.1, 0, epsilon));
			} else {
				graph_Xt->SetPointY(j + 2,-graph_Xt->GetPointY(j + 1));
				graph_Yt->SetPointY(j + 2,graph_Yt->GetPointY(j + 1));
			}
		} else {
			graph_Xt->SetPointY(j + 2,func_Xt->GetX(0, -pend_L - 0.1, 0,  epsilon));
			graph_Yt->SetPointY(j + 2,func_Yt->GetX(0, -pend_L - 0.1, 0, epsilon));
		}
//////////////////////////////////////////////////////////////////////////////////////
*/

//////////////////////////////////////////////////////////////////////////////////////
		dis_x = coef_Bx*coef_Bx - 4*coef_Ax*coef_Cx;
		dis_y = coef_By*coef_By - 4*coef_Ay*coef_Cy;
		if (dis_x > 0){
			if(val_x0 > 0){
				graph_Xt->SetPointY(j + 2, (-coef_Bx + sqrt(dis_x))/coef_Ax/2);
			} else {
				graph_Xt->SetPointY(j + 2, (-coef_Bx - sqrt(dis_x))/coef_Ax/2);
			}
		} else {
			graph_Xt->SetPointY(j + 2, -coef_Bx/coef_Ax/2);
		}
		if (dis_y > 0){
			if(val_y0 > 0){
				graph_Yt->SetPointY(j + 2, (-coef_By + sqrt(dis_y))/coef_Ay/2);
			} else {
				graph_Yt->SetPointY(j + 2, (-coef_By - sqrt(dis_y))/coef_Ay/2);
			}
		} else {
			graph_Yt->SetPointY(j + 2, -coef_By/coef_Ay/2);
		}
		graph_X2Y2->SetPointY(j + 2, graph_Xt->GetPointY(j+2)*graph_Xt->GetPointY(j+2) + graph_Yt->GetPointY(j+2)*graph_Yt->GetPointY(j+2));
//////////////////////////////////////////////////////////////////////////////////////

		if (j == check){
//			cout << "k=" << val_k << " m=" << val_m << " l=" << val_l << endl;
			cout << "k=" << val_k << " x1=" << val_x1 << " y0=" << val_y0 << " x0=" << val_x0 << endl;
			cout << "Ax=" << coef_Ax << " Bx=" << coef_Bx << " Cx=" << coef_Cx << endl;
			cout << "Ay=" << coef_Ay << " By=" << coef_By << " Cy=" << coef_Cy << endl;
			cout << "disX=" << dis_x << " disY=" << dis_y << endl;
		
			cout << "j="<<j << " t="<< time_h*j << " x=" << graph_Xt->GetPointY(j) << " y="<< graph_Yt->GetPointY(j)<<endl;
			cout<<"j+1="<<j+1<<" t="<<time_h*(j+1)<<" x="<<graph_Xt->GetPointY(j+1)<<" y="<<graph_Yt->GetPointY(j+1)<<endl;
			cout<<"j+2="<<j+2<<" t="<<time_h*(j+2)<<" x="<<graph_Xt->GetPointY(j+2)<<" y="<<graph_Yt->GetPointY(j+2)<<endl;
		}	
		//cout << "j=" << j << " t=" << time_h*j << " x=" << graph_Xt->GetPointY(j) << " y=" << graph_Yt->GetPointY(j) << endl;
	}

	cout << "/////////////////" << endl;
	
	for (int j = 0; j <= numb_N; j++){
//		cout << "j=" << j << " t=" << time_h*j << " x=" << graph_Xt->GetPointY(j) <<" y=" << graph_Yt->GetPointY(j) << " x2y2="<< graph_X2Y2->GetPointY(j) << endl;
	}

//	graph_Xt->Draw();
	graph_Yt->Draw();
//	graph_X2Y2->Draw();

		return numb_N;

}
/* for N=10000 and tmax=100 crash programm: 
j=7565 t=75.65 x=-4.453 y=-2.36228 x2y2=25.4095
j=7566 t=75.66 x=-4.70871 y=-1.68168 x2y2=25
j=7567 t=75.67 x=-4.95137 y=-0.973887 x2y2=25.4645
j=7568 t=75.68 x=-4.99614 y=-0.196401 x2y2=25
j=7569 t=75.69 x=-5.01332 y=0.591676 x2y2=25.4834
j=7570 t=75.7 x=-4.80357 y=1.3877 x2y2=25
j=7571 t=75.71 x=-4.54511 y=2.25307 x2y2=25.7344
j=7572 t=75.72 x=-3.9785 y=3.02845 x2y2=25
j=7573 t=75.73 x=-3.34482 y=3.8511 x2y2=26.0188
j=7574 t=75.74 x=-2.36194 y=4.40695 x2y2=25
j=7575 t=75.75 x=-1.28072 y=4.97638 x2y2=26.4046
j=7576 t=75.76 x=0.0929565 y=4.99914 x2y2=25
j=7577 t=75.77 x=1.48015 y=4.91279 x2y2=26.3263
j=7578 t=75.78 x=2.85395 y=4.10548 x2y2=25
j=7579 t=75.79 x=4.47956 y=2.92203 x2y2=28.6048
j=7580 t=75.8 x=4.9979 y=0.144763 x2y2=25
j=7581 t=75.81 x=3.26279 y=-3.45698 x2y2=22.5965
j=7582 t=75.82 x=0.205591 y=-7.098 x2y2=50.4238
j=7583 t=75.83 x=-6.20879 y=-5.86002 x2y2=72.8889
j=7584 t=75.84 x=-7.46008 y=-0.744648 x2y2=56.2073
j=7585 t=75.85 x=-29.7222 y=6.49939 x2y2=925.65
j=7586 t=75.86 x=-4.20665 y=10.8578 x2y2=135.588
j=7587 t=75.87 x=51.0551 y=1394.63 x2y2=1.94761e+06
j=7588 t=75.88 x=626.571 y=2160.9 x2y2=5.0621e+06
j=7589 t=75.89 x=1.70371e+08 y=-14990.4 x2y2=2.90262e+16
j=7590 t=75.9 x=1.27291e+14 y=-5.87595e+08 x2y2=1.62029e+28
j=7591 t=75.91 x=4.3365e+21 y=1.09951e+10 x2y2=1.88052e+43
j=7592 t=75.92 x=2.39561e+38 y=0 x2y2=5.73894e+76
j=7593 t=75.93 x=4.63378e+57 y=0 x2y2=2.14719e+115
j=7594 t=75.94 x=2.39561e+38 y=-4.49286e+72 x2y2=2.01858e+145
j=7595 t=75.95 x=0 y=-inf x2y2=inf
j=7596 t=75.96 x=-inf y=-inf x2y2=inf
j=7597 t=75.97 x=-nan y=-nan x2y2=-nan
*/


