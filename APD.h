/*

    simple function to measure the APDs

    Haibo Ni

    qiangzi.ni@gmail.com

    update Feb 13. 2014
*/

#ifndef APD_H
#define APD_H

#include <stdio.h>

typedef struct
{
    double Vmax, Vmin, timeAPDstart, timeAPDend, dVdtmax, APD_out_end[2];
    int APD_switch, APD_count, Vmax_switch, APD_out_swhich;
    double APD90, APD50, APD20, APD75;
    double Vm_prev;
    double AMP, AMP_last;
    double APD90_prev;
    double t_since_up_stroke;
    double v90, v75, v50, v20;
    double ICaL_in, RyR_in;
    double Istim_prev;
    double Diastolic_value;
    double Systolic_value_H, Systolic_value_L;
    unsigned int N_stim;
    double Diastolic_value_2, Systolic_value_L_2, Systolic_value_H_2;
    FILE *out;
} AP_infor;

void AP_information_initialise(AP_infor *AP);
void APD90_measure(double t, double Istim, double BCL, double dt, double Vm, AP_infor *AP, FILE *out);
void APD90_DS_Value_measure(double t, double Istim, double BCL, double dt, double Vm, double Value, AP_infor *AP, FILE *out);
void APD90_DS_double_Value_measure(double t, double Istim, double BCL, double dt, double Vm, double Value_1, double Value_2, AP_infor *AP, FILE *out);
void APD90_DS_Value_measure_with_INa(double INa_threshold, double INa, double t, double Istim, double BCL, double dt, double Vm, double Value, AP_infor *AP, FILE *out);
void APD90_DS_Value_measure_with_StrokeTime_threshold(double StrokeTime, double t, double Istim, double BCL, double dt, double Vm, double Value, AP_infor *AP, FILE *out);
void calcium_accumulation(double t, double Istim, double BCL, double dt, double ICaL, double J_RyR, AP_infor *AP, FILE *out);
#endif
