#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

/* Header files with a description of contents used */
#include <cvodes/cvodes.h>                  /* prototypes for CVODE fcts., consts.  */
#include <cvodes/cvodes_dense.h>
//#include <cvodes/cvodes_band.h>
//#include <cvode/cvode_diag.h>
//#include <cvodes/cvodes_spgmr.h>
///#include <cvode/cvode_sptfqmr.h>
//#include <cvode/cvode_spbcgs.h>
#include <nvector/nvector_serial.h>       /* serial N_Vector types, fcts., macros */

#include "APD.c"

#define Ith(v,i)    NV_Ith_S(v,i-1)  

#define Y0  RCONST(-80.50146)  //Volt
#define Y1  RCONST(0.004164108) // m
#define Y2  RCONST(0.6735613) // h
#define Y3  RCONST(0.6729362)  //j
#define Y4  RCONST(0.000002171081)  //d
#define Y5  RCONST(0.9999529) //f11
#define Y6  RCONST(0.9999529)  //f12
#define Y7  RCONST(9.913102e-01)  //Cainact
#define Y8  RCONST(0.002191519)  //r
#define Y9  RCONST(0.9842542)   //s
#define Y10 RCONST(0.6421196)  //sslow
#define Y11 RCONST(0.002907171)  //rss
#define Y12 RCONST(0.3142767) //sss
#define Y13 RCONST(0.003578708) //y
#define Y14 RCONST(10.73519) //Nai
#define Y15 RCONST(139.2751) //Ki
#define Y16 RCONST(0.00007901351) //Cai
#define Y17 RCONST(0.06600742) //CaNSR
#define Y18 RCONST(0.00008737212) // CaSS 
#define Y19 RCONST(0.06607948) // CaJSR
#define Y20 RCONST(0.6348229) //PC1
#define Y21 RCONST(0.0004327548) // Po1
#define Y22 RCONST(0.0000000006062540) //Po2
#define Y23 RCONST(0.3647471) //PC2
#define Y24 RCONST(0.005161900) //LTRPNCa  
#define Y25 RCONST(0.1394301) //HTRPNCa


#define NEQ 26
#define DT 0.1

double I_stim = 0.0;
int ctype=1;
static int update(realtype t, N_Vector y, N_Vector ydot, void *user_data);

//定义全局变量
const double Vss = 1.2e-09;//transfer pL to uL
const double VJSR = 5.6e-08;//transfer pL to uL
const double VNSR = 5.04e-07;//transfer pL to uL
const double Vmyo = 9.36e-06;//transfer pL to uL
const double Ko = 5.4;//mM
const double Nao = 140.0;//mM
const double Cao = 1.2;//mM
// Membrane capacitance
const double Cm = 100e-06; //μF
// Faraday constant
const double F = 96487.0; //C/mol
// Absolute temperature
const double T = 273 + 22;//295; //K
// Ideal gas constant
const double R = 8314.5; //mJ/molK
// Maximum conductance for ICaL
const double gCaL = 0.031;
// Maximum conductance for Iss
const double gss = 0.007;
// Maximum conductance for IK1
const double gK1 = 0.024;
// Maximum conductance for IBNa
const double gBNa = 8.015e-05;
// Maximum conductance for IBCa
const double gBCa = 3.24e-05;
// Maximum conductance for IBK
const double gBK = 13.8e-05;
// Maximum conductance for If
const double gf = 0.00145;
// Maximum INaK current
const double INaKmax = 0.08;//单位nA
// Half-maximum Na+ binding constant for INaK
const double kmNa = 10.0;// 单位mM
// Half-maximum Na+ binding constant for INaK
const double kmK = 1.5;// 单位mM
// Maximum ICaP current
const double ICaPmax = 0.004;//单位nA
// Scaling factor for INaCa
const double kNaCa = 0.9984e-05; //单位是(mM)-4
// Denominator constant for INaCa
const double dNaCa =  0.0001; //单位是(mM)-4
// Position of energy barrier controlling voltage dependence for INaCa
const double YNaCa = 0.5;
const double v1 = 1.8e03;
const double Kfb = 0.168e-03;
const double Krb = 3.29;
const double KSR = 1.0;
const double Nfb = 1.2;
const double Nrb = 1.0;
const double vmaxf = 0.4e-1;
const double vmaxr = 0.9;
const double tautr = 0.5747e-03;
const double tauxfer = 26.7e-3;
const double Kapos = 12.15e12;
const double Kaneg = 0.576e03;
const double Kbpos = 4.05e09;
const double Kbneg = 1.930e03;
const double Kcpos = 0.1e03;
const double Kcneg = 0.0008e03;
const double nRyR = 4;
const double mRyR = 3;// to distinguish with the original m
const double LTRPNtot = 70e-03;
const double HTRPNtot = 140e-03;
const double Khtrpnpos = 200e03;
const double Khtrpnneg = 66.0e-03;
const double Kltrpnpos = 40e03;
const double Kltrpnneg = 0.04e03;
const double CMDNtot = 50.0e-03;
const double CSQNtot = 15.0;
const double EGTAtot = 10.0;
const double KmCMDN = 2.38e-03;
const double KmCSQN = 0.8;
const double KmEGTA = 1.5e-04;
// inverse potential
const double ECaL = 65.0;
// constants used in equtations
const double fNa = 0.2;
const double fK  = 0.8; // fK = 1-fNa;
// time constant(tau) variables
const double tauCainact = 0.009;
const double tausss = 2.1;
// ********************** KATP Current, IKATP related **********************
const double Mgi = 0.3;//0.5;//mmol/L
//const double numKATP = 1. * 8500;
// ********** KATP configuration ***********
const double ATPi = 6800;//8000;//uM/L
const double ADPi = (7192.86 - 6800)*84/2200;//8000;//uM/L   //(7000 - ATPi)*0.05;

int main(int argc, char *argv[]){

	int i_beatsnum, d_t_total, i_pcl;
	double d_tmp, d_tmp1, d_tmp2 = 0.0;
	realtype t, d_t;	// realtype is double, typedefined in cvode  
	AP_infor AP; // 结构体，把AP相关的东西全部封在里面了
	void *p_cvode_mem = NULL; // this pointer points to the memory block used by the cvode
	N_Vector p = NULL;
	FILE *f_apd_measure_out, *f_current;
	clock_t time_begin, time_end;
	i_pcl = atoi(argv[1]);
	i_beatsnum = atoi(argv[2]);
	AP_information_initialise(&AP);// 把所有值初始化赋值0

	d_t_total = i_beatsnum * i_pcl;

	p = N_VNew_Serial(NEQ); // for solver
	// TODO can these lines be shortened using macros?
	Ith(p, 1) = Y0; 
	Ith(p, 2) = Y1;
	Ith(p, 3) = Y2;
	Ith(p, 4) = Y3;
	Ith(p, 5) = Y4;
	Ith(p, 6) = Y5;
	Ith(p, 7) = Y6;
	Ith(p, 8) = Y7;
	Ith(p, 9) = Y8;
	Ith(p, 10) = Y9;
	Ith(p, 11) = Y10;
	Ith(p, 12) = Y11;
	Ith(p, 13) = Y12;
	Ith(p, 14) = Y13;
	Ith(p, 15) = Y14;
	Ith(p, 16) = Y15;
	Ith(p, 17) = Y16;
	Ith(p, 18) = Y17;
	Ith(p, 19) = Y18;
	Ith(p, 20) = Y19;
	Ith(p, 21) = Y20;
	Ith(p, 22) = Y21;
	Ith(p, 23) = Y22;
	Ith(p, 24) = Y23;
	Ith(p, 25) = Y24;
	Ith(p, 26) = Y25;
	
	p_cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
	CVodeInit(p_cvode_mem, update, 0.0, p); 
	CVodeSStolerances(p_cvode_mem, 10e-6, 1e-6); //10e-6  //迭代控制相对误差
	CVDense(p_cvode_mem, NEQ);
	//CVBand(p_cvode_mem, NEQ,3,3);
	//CVDiag(p_cvode_mem);直接崩掉
	//CVSpgmr(p_cvode_mem,PREC_LEFT,0);//速度稍微快一点了 
	//CVSptfqmr(p_cvode_mem,PREC_NONE,0); 速度一般
	//CVSpbcg(p_cvode_mem,PREC_NONE,0); 不行
	//CVodeSetMaxNumSteps(p_cvode_mem,5); 
	CVodeSetMaxStep(p_cvode_mem, DT);


	// 写文件时用到的一些文件名一类
	char filename[30];
	sprintf(filename, "current_%d.dat", i_pcl);
	f_current = fopen(filename, "w");
	sprintf(filename, "APD_measure_out_%d.dat", i_pcl);
	f_apd_measure_out = fopen(filename, "w");
	time_begin = clock();
	int count = 0;

	// simulation begins from 0ms, but the result begins from dt, for example 0.1ms.
	for (d_t = DT; d_t <= d_t_total; d_t += DT,count++) // 这里属于风格问题，d_t代表我自己代码里的t或者time
	{
		// deal with the precision problem of double, if don't, the error will accumulate
		d_t = (long) (d_t * 1000000 + 0.5) / 1000000.0;

		// stimulate during the first 5ms of each pcl
		if (fmod(d_t, (double) i_pcl) >= 25 && fmod(d_t, (double) i_pcl) < 30)
			I_stim = -0.6;
		else
			I_stim = 0.0;

		CVode(p_cvode_mem, d_t, p, &t, CV_NORMAL);  //CV_NORMAL  选择步长

		// print a message every second
		modf(d_t / 1000, &d_tmp1); // d_tmp1取的是(d_t /1000)的整数部分 ，该函数的返回值是其小数部分
		if (d_tmp2 != d_tmp1) // 
			printf("%5ds simulation has completed.\n", (int) floor(d_t / 1000));// 以秒来计时，没甚么难的
		modf(d_t / 1000, &d_tmp2);

		// print data to files when the time is on each 0.1 second
		if (d_t >= 0 && modf(d_t*100, &d_tmp) < 1e-12 && count%2 ==0)
		{
			fprintf(f_current,"%20.10f%20.10f\n",d_t, Ith(p, 1));
		}

		APD90_DS_double_Value_measure(t, I_stim, i_pcl, DT, Ith(p, 1),	Ith(p, 18), Ith(p, 17), &AP, f_apd_measure_out); // f_apd_measure_out是文件名
		// APD90_measure(t, I_stim, i_pcl, dt, Ith(y, 39), &AP, f_apd_measure_out);
	} /* end of time loop */

	time_end = clock();
	//fclose(f_output_all);
	fclose(f_current);
	fclose(f_apd_measure_out);
	printf("CPU time = %f seconds\n",
			((double) (time_end - time_begin)) / CLOCKS_PER_SEC); // 统计程序运行时间
	// printf("%f %f %f\n", i_pcl, AP.APD_out_end[1], AP.APD_out_end[0] );
	printf("%d %f %f %.10f %.10f\n", i_pcl, AP.APD90, AP.APD50, // 统计各项指标
			AP.Systolic_value_L, AP.Systolic_value_H);

	return 0;
}



static int update(realtype t, N_Vector p, N_Vector pdot, void *user_data)
{
		double IK1,IB,INa,It,Iss,If,ICaL,INaCa,ICaP,INaK,IfNa,IfK,IBK,IBCa,IBNa,IKATP;
		double EK,ENa;
		double gNa,gt,a,b;
		double m_,h_  ,j_  ,d_  ,f11_,f12_,r_  ,s_  ,sslow_,rss_,sss_,Cainact_,yinf;
		double taum,tauh,tauj,taud,tauf11,tauf12,taur,taus,tausslow,taurss,tauy;
		double sigma,Jrel,Jup,Jtr,Jxfer,Jtrpn,fb,rb,betai,betass,betaJSR,temp;
		double KhMg,KhMgKo,delNa,delMg,KhNa0,KhNa,g0,gamma0,fMgATP,fNaATP,fTATP,KmATP,fATP,	hATP;

		realtype P[NEQ];
		realtype dP[NEQ];
		for(int i=0; i<NEQ; i++)
			P[i] = Ith(p, i+1);

	    // ******************* Compute Membrane currents, gate infinite & tau*************************
		// Na+ current
		
		if(ctype == 1) 
		{
			gNa = 0.8;
			gt = 0.035*1.073;
			a = 0.886;
			b = 0.114;
		}	
	
		if(ctype == 0) 
		{
			gNa = 1.33*0.8;
			gt = 0.4647*0.035*1.165;
			a = 0.583;
			b = 0.417;
		}
		ENa = (R*T/F)*log(Nao/P[14]);
		INa = gNa*pow(P[1],3)*P[2]*P[3]*(P[0]-ENa);	
		// Na+ gate infinite & tau
		m_ = 1.0/(1.0+exp((P[0]+45.0)/-6.5));
		h_ = 1.0/(1.0+exp((P[0]+76.1)/6.07));
		j_ = h_;
		taum = 0.00136/(((0.32*(P[0]+47.13))/(1.0-exp(-0.1*(P[0]+47.13))))+0.08*exp(-P[0]/11.0));
		dP[1] = (m_- P[1])/taum;

		//if(P[0] >= -40) 
		//{
			tauh=0.0004537*(1.0+exp(-(P[0]+10.66)/11.1));
			tauj=0.01163*(1.0+exp(-0.1*(P[0]+32.0)))/exp(-0.0000002535*P[0]);
	//	}
		if(P[0] < -40)
		{
			tauh=0.00349/(0.135*exp(-(P[0]+80.0)/6.8)+3.56*exp(0.079*P[0])+310000*exp(0.35*P[0]));
			temp = ((P[0]+37.78)/(1.0+exp(0.311*(P[0]+79.23))))
					*(-127140*exp(0.2444*P[0])-0.00003474*exp(-0.04391*P[0]))
					+0.1212*exp(-0.01052*P[0])/(1.0+exp(-0.1378*(P[0]+40.14)));
			tauj=0.00349/temp;
		}
		dP[2] = (h_- P[2])/tauh;
		dP[3] = (j_- P[3])/tauj;

		// L-type Ca2+ current		
		ICaL = gCaL*P[4]*((0.9+P[7]/10.0)*P[5]+(0.1-P[7]/10.0)*P[6])*(P[0]-65.0);
		d_ = 1.0/(1+exp((P[0]+15.3)/-5.0));
		f11_ = 1.0/(1+exp((P[0]+26.7)/5.4));
		f12_ = f11_;
		Cainact_ = 1.0/(1.0+P[18]/0.01);
		taud = 0.00305*exp(-0.0045*pow(P[0]+7.0,2))+0.00105*exp(-0.002*pow(P[0]-18.0,2))+0.00025;
		dP[4] = (d_- P[4])/taud;


		tauf11 = 0.105*exp(-pow((P[0]+45.0)/12.0,2))
				+0.04/(1.0+exp((-P[0]+25.0)/25.0))
				+0.015/(1.0+exp((P[0]+75.0)/25.0))
				+0.0017;
		tauf12 = 0.041*exp(-pow((P[0]+47.0)/12.0,2))
				+0.08/(1.0+exp((P[0]+55.0)/-5.0))
				+0.015/(1.0+exp((P[0]+75.0)/25.0))
				+0.0017;
		dP[5] = (f11_- P[5])/tauf11;
		dP[6] = (f12_- P[6])/tauf12;
		dP[7] = (Cainact_- P[7])/tauCainact;

		// Ca2+-independent transient outward K+ current, It
		EK = (R*T/F)*log(Ko/P[15]);
		It = gt*P[8]*(a*P[9]+b*P[10])*(P[0]-EK);
		// Ca2+-independent transient outward K+ gate infinite & tau
		r_ = 1.0/(1+exp((P[0]+10.6)/-11.42));
		s_ = 1.0/(1+exp((P[0]+45.3)/6.8841));
		sslow_ = s_;
		taur = 1.0/(45.16*exp(0.03577*(P[0]+50.0))+98.9*exp(-0.1*(P[0]+38.0)));
		dP[8] = (r_- P[8])/taur;

		//if(ctype == 0) 
		//{
			taus = 0.35*exp(-pow((P[0]+70.0)/15.0,2))+0.035;
			tausslow = 3.7*exp(-pow((P[0]+70.0)/30.0,2))+0.035;
		//}


		if(ctype == 1) 
		{
			taus = 0.55*exp(-pow((P[0]+70.0)/25.0,2))+0.049;
			tausslow = 3.3*exp(-pow((P[0]+70.0)/30.0,2))+0.049;
		}
		dP[9] = (s_- P[9])/taus;
		dP[10] = (sslow_- P[10])/tausslow;

		// Steady-state outward K+ current,Iss, gate and tau update
		Iss = gss*P[11]*P[12]*(P[0]-EK);
		rss_ = 1.0/(1+exp((P[0]+11.5)/-11.82));
		sss_ = 1.0/(1+exp((P[0]+87.5)/10.3));
		taurss = 10.0/(45.16*exp(0.03577*(P[0]+50.0))+98.9*exp(-0.1*(P[0]+38.0)));
		dP[11] = (rss_- P[11])/taurss;
		dP[12] = (sss_- P[12])/tausss;

		// Inward rectifier, IK1 
		IK1 = (48.0/(exp((P[0]+37.0)/25.0)+exp((P[0]+37.0)/-25.0))+10.0)*(0.001/(1.0+exp((P[0]-EK-76.77)/-17.0)))
				+gK1*(P[0]-EK-1.73)/( 1.0+exp(1.613*F*(P[0]-EK-1.73)/(R*T)) *(1.0+exp((Ko-0.9988)/-0.124)));

		// Hyperpolarization-activated current, If, gate and tau update
		IfNa = gf*P[13]*fNa*(P[0]-ENa);
		IfK = gf*P[13]*fK*(P[0]-EK);
		If = IfNa + IfK;
		yinf = 1.0/(1+exp((P[0]+138.6)/10.48));
		tauy = 1.0/(0.11885*exp((P[0]+80.00)/28.37)+0.56236*exp((P[0]+80.00)/-14.19));
		dP[13] = (yinf- P[13])/tauy;

		// Background currents
		IBNa = gBNa*(P[0]-ENa);
		IBK = gBK*(P[0]-EK);
		IBCa = gBCa*(P[0]-ECaL);
		IB = IBNa+IBCa+IBK;


		// Na+-K+ pump current,INaK
		sigma = (exp(Nao/67.3)-1.0)/7.0;
		INaK = INaKmax*(1.0/(1.0+0.1245*exp(-0.1*P[0]*F/(R*T))+0.0365*sigma*exp(-P[0]*F/(R*T))))*(Ko/(Ko+kmK))*(1.0/(1.0+pow(kmNa/P[14],1.5)));


		// Sarcolemmal Ca2+ pump current, ICaP
		ICaP = ICaPmax*(P[16]/(P[16]+0.0004));


		// Na+-Ca2+ ion exchanger current, INaCa
		INaCa = kNaCa*((P[14]*P[14]*P[14]*Cao*exp(0.03743*YNaCa*P[0])-Nao*Nao*Nao*P[16]*exp(0.03743*(YNaCa-1.0)*P[0]))
				/(1.0+dNaCa*(Nao*Nao*Nao*P[16]+P[14]*P[14]*P[14]*Cao)));



		// ATP sensitive K+ current, IKATP
		gamma0 = 1150000 * pow(1.0*Ko/5.4,0.24);//pS
		

		delMg = 0.32;
		KhMgKo = 0.65/pow(Ko+5.0,0.5);
		KhMg = KhMgKo * exp(-2.0*delMg*F*P[0]/(R*T));
		fMgATP = 1.0/(1.0+pow(1.0*Mgi/KhMg,1));

		delNa = 0.35;
		KhNa0 = 25.9; //mmol/L
		KhNa = KhNa0 * exp(-delNa*F*P[0]/(R*T));
		fNaATP = 1.0/(1.0+pow(P[14]/KhNa,2));

		fTATP = pow(1.3,(T - 309.0)/10.0);

		g0 = gamma0*fMgATP*fNaATP*fTATP;//pS

		
		// version 3.0 
		//if(ctype == 0) 
		//{
			KmATP = 0.4*0.7*(35.8 + 17.9*pow(ADPi,0.256)); // Shimokawa version
			hATP = (1.3 + 0.74*6*exp(-0.09*ADPi)); // Shimokawa version
		//}

		if(ctype == 1)
		{
			KmATP = 0.7*(35.8 + 17.9*pow(ADPi,0.256)); // Shimokawa version
			hATP = (1.3 + 0.74*6*exp(-0.09*ADPi)); // Shimokawa version
		}

		
		//hATP = 1.3 + 0.74*6*exp(-0.99*ADPi);
		fATP = 1.0/(1.0 + pow(ATPi/KmATP,hATP));
		IKATP = 1.0e-6*g0*fATP*(P[0] - EK);

		// ***************************** Ca2+ handling mechanisms ***********************************
		// Calcium release channel in sarcoplasmic reticulum
		Jrel = v1*(P[21]+P[22])*(P[19]-P[18]);
		
		// SERCA2a Ca2+ pump
		fb = pow(P[16]/Kfb,Nfb);
		rb = pow(P[17]/Krb,Nrb);
		Jup = KSR*(vmaxf*fb-vmaxr*rb)/(1.0+fb+rb);

		// Intracellular and sarcoplasmic reticulum Ca2+ fluxes
		Jtr = (P[17]-P[19])/tautr;
		Jxfer = (P[18]-P[16])/tauxfer;
		dP[25] = Khtrpnpos*P[16]*(HTRPNtot-P[25])-Khtrpnneg*P[25];	
		dP[24] = Kltrpnpos*P[16]*(LTRPNtot-P[24])-Kltrpnneg*P[24];
		Jtrpn = dP[25] + dP[24];

		// Calcium release channel in sarcoplasmic reticulum
		dP[20] = -Kapos*pow(P[18],nRyR)*P[20]+Kaneg*P[21];
		dP[21] = Kapos*pow(P[18],nRyR)*P[20]-Kaneg*P[21]-Kbpos*pow(P[18],mRyR)*P[21]+Kbneg*P[22]-Kcpos*P[21]+Kcneg*P[23];	
		dP[22] = -Kbneg*P[22] + Kbpos*pow(P[18],mRyR)*P[21];
		dP[23] = Kcpos*P[21]-Kcneg*P[23];


		// ********************** Update Intracellular ion concentrations ****************************
		dP[14] = -(INa+IBNa+3.0*INaCa+3.0*INaK+IfNa)*1.0/(Vmyo*F);
		dP[15] = -(IKATP+Iss+IBK+It+IK1+IfK-2*INaK)*1.0/(Vmyo*F);	
		betai = 1.0/(1.0+CMDNtot*KmCMDN/pow(KmCMDN+P[16],2)+EGTAtot*KmEGTA/pow(KmEGTA+P[16],2));
		// betai = 1.0/(1.0+CMDNtot*KmCMDN/pow(KmCMDN+P[16],2));
		betass = 1.0/(1.0+CMDNtot*KmCMDN/pow(KmCMDN+P[18],2));
		betaJSR = 1.0/(1.0+CSQNtot*KmCSQN/pow(KmCSQN+P[19],2));
		dP[16] = betai*(Jxfer-Jup-Jtrpn-(IBCa-2.0*INaCa+ICaP)*1.0/(2.0*Vmyo*F));
		dP[18] = betass*(Jrel*VJSR/Vss - Jxfer*Vmyo/Vss - ICaL*1.0/(2.0*Vss*F));	
		dP[19] = betaJSR*(Jtr-Jrel);	
		dP[17] = Jup*Vmyo/VNSR - Jtr*VJSR/VNSR;
		dP[0] = -(INa+ICaL+It+Iss+If+IK1+IB+INaK+INaCa+ICaP+I_stim+IKATP)/Cm;



		for (int i = 0; i < NEQ; i++)
			Ith(pdot, i + 1) = dP[i]*0.001; //0.001用来把计量单位换算成秒

		return 0;

}

