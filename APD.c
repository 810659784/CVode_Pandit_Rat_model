#include <stdio.h>
#include <math.h>
#include "APD.h"

void AP_information_initialise(AP_infor *AP)
{
	AP->Vmax = 0.0;
	AP->Vmin = 0.0;
	AP->timeAPDstart = 0.0;
	AP->timeAPDend = 0.0;
	AP->dVdtmax = 0.0;
	AP->APD_out_end[0] = 0.0;
	AP->APD_out_end[1] = 0.0;

	AP->APD_switch = 0;
	AP->APD_count = 0;
	AP->Vmax_switch = 0;
	AP->APD_out_swhich = 0;
	AP->APD90 = AP->APD75 = AP->APD50 = AP->APD20 = 0.0;
	AP->Vm_prev = 0.0;
	AP->v90 = 0.0;
	AP->v75 = 0.0;
	AP->v50 = 0.0;
	AP->v20 = 0.0;
	AP->ICaL_in = 0.0;
	AP->RyR_in = 0.0;
	AP->Istim_prev = 0.0;
	AP->AMP = 0.0;
	AP->AMP_last = 0.0;
	;
	AP->APD90_prev = 0.0;
	AP->t_since_up_stroke = 0.0;

	AP->Diastolic_value = 0.0;
	AP->Systolic_value_L = 0.0;
	AP->Systolic_value_H = 0.0;
	AP->Diastolic_value_2 = 0.0;
	AP->Systolic_value_L_2 = 0.0;
	AP->Systolic_value_H_2 = 0.0;
	AP->N_stim = 0;
}

void APD90_measure(double t, double Istim, double BCL, double dt, double Vm,
		AP_infor *AP, FILE *out)
{

	if (((Istim >= 1e-3 || Istim <= -1e-3) && (fabs(AP->Istim_prev) < 1e-10)))
	{
		AP->APD_switch = 1;
		AP->Vmax = Vm;
		AP->timeAPDstart = t;
		AP->Vmin = Vm;
		AP->dVdtmax = 0.0;
		AP->N_stim++;
		AP->Vmax_switch = 0;

		// printf("ssssssssss\n");
	}

	if (AP->APD_switch == 1)
	{
		if (Vm > -30.0)
		{
			if (Vm >= AP->Vmax)
			{
				AP->Vmax = Vm;
			}
			else
				AP->Vmax_switch = 1;
			if (AP->Vmax_switch == 1)
			{
				AP->APD_switch = 2;
				AP->AMP_last = AP->AMP;
				AP->AMP = AP->Vmax - AP->Vmin;
				/*printf("%f %f\n", AP->Vmax, AP->Vmin);

				 printf("%f %f\n", AP->AMP, AP->AMP_last);*/
				AP->v90 = AP->Vmax - 0.9 * (AP->Vmax - AP->Vmin);
				AP->v75 = AP->Vmax - 0.75 * (AP->Vmax - AP->Vmin);
				AP->v50 = AP->Vmax - 0.50 * (AP->Vmax - AP->Vmin);
				AP->v20 = AP->Vmax - 0.20 * (AP->Vmax - AP->Vmin);
				// printf("asfde\n");
			}
		}
	}
	if (AP->APD_switch == 2)
	{

		if ((AP->Vm_prev >= AP->v20) && (Vm <= AP->v20))
		{
			AP->APD20 = t - AP->timeAPDstart;
			// printf("aaaa\n");

		}
		else if ((AP->Vm_prev >= AP->v50) && (Vm <= AP->v50))
		{
			AP->APD50 = t - AP->timeAPDstart;
			// printf("bbbbbbbbb\n");

		}
		else if (AP->Vm_prev >= AP->v75 && Vm <= AP->v75)
		{
			AP->APD75 = t - AP->timeAPDstart;
			// printf("ccccccccccc\n");

		}
		else if (AP->Vm_prev >= AP->v90 && Vm <= AP->v90)
		{
			AP->APD_switch = 0;
			AP->APD_count++;
			AP->Vmax_switch = 0;
			AP->APD90 = t - AP->timeAPDstart;
			// printf("dddddddddddd\n");
			fprintf(out, "%.2f %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n",
					AP->APD_count * BCL, AP->APD90, AP->APD75, AP->APD50,
					AP->APD20, AP->Vmax, AP->Vmin, AP->dVdtmax);
			if (AP->APD_out_swhich == 0)
			{
				AP->APD_out_end[0] = t - AP->timeAPDstart;
				AP->APD_out_swhich = 1;
			}
			else if (AP->APD_out_swhich == 1)
			{
				AP->APD_out_end[1] = t - AP->timeAPDstart;
				AP->APD_out_swhich = 0;
			}
		}
	}

	double dVdt = (Vm - AP->Vm_prev) / dt;
	if (dVdt > AP->dVdtmax)
		AP->dVdtmax = dVdt;

	AP->Vm_prev = Vm;
	AP->Istim_prev = Istim;
}

void APD90_DS_Value_measure(double t, double Istim, double BCL, double dt,
		double Vm, double Value, AP_infor *AP, FILE *out)
{
	if (((Istim >= 1e-8 || Istim <= -1e-8) && (fabs(AP->Istim_prev) < 1e-10)
			&& (AP->APD_switch == 0)))
	{
		if (AP->APD_count > 50)
		{
			// printf("dddddddddddd\n");
			// printf("%f %f %f\n", AP->Diastolic_value, AP->Systolic_value_L, AP->Systolic_value_H);

			fprintf(out,
					"%.2f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n",
					AP->APD_count * BCL, BCL, AP->APD90, AP->APD75, AP->APD50,
					AP->APD20, AP->Systolic_value_L, AP->Systolic_value_H,
					AP->Vmax, AP->Vmin, AP->dVdtmax);

		}
		AP->APD_switch = 1;
		AP->Vmax = -60;
		AP->timeAPDstart = t;
		AP->Vmin = Vm;
		AP->dVdtmax = 0.0;
		AP->Diastolic_value = Value;
		AP->Systolic_value_L = Value;
		AP->Systolic_value_H = Value;
		AP->N_stim++;
		AP->Vmax_switch = 0;

		// printf("ssssssssss\n");
	}
	if (AP->Systolic_value_L > Value)
	{
		AP->Systolic_value_L = Value;
	}
	if (AP->Systolic_value_H < Value)
	{
		AP->Systolic_value_H = Value;
	}
	if (AP->APD_switch == 1)
	{
		if (Vm > -80.0)
		{
			if (Vm >= AP->Vmax)
			{
				AP->Vmax = Vm;
			}
			else
				AP->Vmax_switch = 1;
			if (AP->Vmax_switch == 1)
			{
				AP->APD_switch = 2;
				AP->AMP_last = AP->AMP;
				AP->AMP = AP->Vmax - AP->Vmin;
				/*printf("%f %f\n", AP->Vmax, AP->Vmin);

				 printf("%f %f\n", AP->AMP, AP->AMP_last);*/
				AP->v90 = AP->Vmax - 0.9 * (AP->Vmax - AP->Vmin);
				AP->v75 = AP->Vmax - 0.75 * (AP->Vmax - AP->Vmin);
				AP->v50 = AP->Vmax - 0.50 * (AP->Vmax - AP->Vmin);
				AP->v20 = AP->Vmax - 0.20 * (AP->Vmax - AP->Vmin);
				// printf("asfde\n");
			}
		}
	}
	if (AP->APD_switch == 2)
	{

		if ((AP->Vm_prev >= AP->v20) && (Vm <= AP->v20))
		{
			AP->APD20 = t - AP->timeAPDstart;
			// printf("aaaa\n");

		}
		else if ((AP->Vm_prev >= AP->v50) && (Vm <= AP->v50))
		{
			AP->APD50 = t - AP->timeAPDstart;
			// printf("bbbbbbbbb\n");

		}
		else if (AP->Vm_prev >= AP->v75 && Vm <= AP->v75)
		{
			AP->APD75 = t - AP->timeAPDstart;
			// printf("ccccccccccc\n");

		}
		else if (AP->Vm_prev >= AP->v90 && Vm <= AP->v90)
		{
			AP->APD_switch = 0;
			AP->APD_count++;
			AP->Vmax_switch = 0;
			AP->APD90 = t - AP->timeAPDstart;
			// printf("dddddddddddd\n");
			// printf("%f %f %f\n", AP->Diastolic_value, AP->Systolic_value_L, AP->Systolic_value_H);

			// fprintf (out, "%.2f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n", AP->APD_count * BCL, AP->APD90, AP->APD75, AP->APD50, AP->APD20, AP->Vmax, AP->Vmin, AP->dVdtmax);
			if (AP->APD_out_swhich == 0)
			{
				AP->APD_out_end[0] = t - AP->timeAPDstart;
				AP->APD_out_swhich = 1;
			}
			else if (AP->APD_out_swhich == 1)
			{
				AP->APD_out_end[1] = t - AP->timeAPDstart;
				AP->APD_out_swhich = 0;
			}
		}
	}

	double dVdt = (Vm - AP->Vm_prev) / dt;
	if (dVdt > AP->dVdtmax)
		AP->dVdtmax = dVdt;
	AP->Istim_prev = Istim;
	AP->Vm_prev = Vm;
}

void APD90_DS_double_Value_measure(double t, double Istim, double BCL,
		double dt, double Vm, double Value_1, double Value_2, AP_infor *AP,
		FILE *out)
{
	if (((Istim >= 1e-8 || Istim <= -1e-8) && (fabs(AP->Istim_prev) < 1e-10)
			&& (AP->APD_switch == 0)))
	{
		if (AP->APD_count > 10)
		{
			// printf("dddddddddddd\n");
			// printf("%f %f %f\n", AP->Diastolic_value, AP->Systolic_value_L, AP->Systolic_value_H);

			fprintf(out,
					"%.2f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n",
					AP->APD_count * BCL, BCL, AP->APD90, AP->APD75, AP->APD50,
					AP->APD20, AP->Systolic_value_L, AP->Systolic_value_H,
					AP->Systolic_value_L_2, AP->Systolic_value_H_2, AP->Vmax,
					AP->Vmin, AP->dVdtmax);

		}
		AP->APD_switch = 1;
		AP->Vmax = -60;
		AP->timeAPDstart = t;
		AP->Vmin = Vm;
		AP->dVdtmax = 0.0;
		AP->Diastolic_value = Value_1;
		AP->Systolic_value_L = Value_1;
		AP->Systolic_value_H = Value_1;
		AP->Diastolic_value_2 = Value_2;
		AP->Systolic_value_L_2 = Value_2;
		AP->Systolic_value_H_2 = Value_2;
		AP->N_stim++;
		AP->Vmax_switch = 0;

		// printf("ssssssssss\n");
	}
	if (AP->Systolic_value_L > Value_1)
	{
		AP->Systolic_value_L = Value_1;
	}
	if (AP->Systolic_value_H < Value_1)
	{
		AP->Systolic_value_H = Value_1;
	}
	if (AP->Systolic_value_L_2 > Value_2)
	{
		AP->Systolic_value_L_2 = Value_2;
	}
	if (AP->Systolic_value_H_2 < Value_2)
	{
		AP->Systolic_value_H_2 = Value_2;
	}
	if (AP->APD_switch == 1)
	{
		if (Vm > -30.0)
		{
			if (Vm >= AP->Vmax)
			{
				AP->Vmax = Vm;
			}
			else
				AP->Vmax_switch = 1;
			if (AP->Vmax_switch == 1)
			{
				AP->APD_switch = 2;
				AP->AMP_last = AP->AMP;
				AP->AMP = AP->Vmax - AP->Vmin;
				/*printf("%f %f\n", AP->Vmax, AP->Vmin);

				 printf("%f %f\n", AP->AMP, AP->AMP_last);*/
				AP->v90 = AP->Vmax - 0.9 * (AP->Vmax - AP->Vmin);
				AP->v75 = AP->Vmax - 0.75 * (AP->Vmax - AP->Vmin);
				AP->v50 = AP->Vmax - 0.50 * (AP->Vmax - AP->Vmin);
				AP->v20 = AP->Vmax - 0.20 * (AP->Vmax - AP->Vmin);
				// printf("asfde\n");
			}
		}
	}
	if (AP->APD_switch == 2)
	{

		if ((AP->Vm_prev >= AP->v20) && (Vm <= AP->v20))
		{
			AP->APD20 = t - AP->timeAPDstart;
			// printf("aaaa\n");

		}
		else if ((AP->Vm_prev >= AP->v50) && (Vm <= AP->v50))
		{
			AP->APD50 = t - AP->timeAPDstart;
			// printf("bbbbbbbbb\n");

		}
		else if (AP->Vm_prev >= AP->v75 && Vm <= AP->v75)
		{
			AP->APD75 = t - AP->timeAPDstart;
			// printf("ccccccccccc\n");

		}
		else if (AP->Vm_prev >= AP->v90 && Vm <= AP->v90)
		{
			AP->APD_switch = 0;
			AP->APD_count++;
			AP->Vmax_switch = 0;
			AP->APD90 = t - AP->timeAPDstart;
			// printf("dddddddddddd\n");
			// printf("%f %f %f\n", AP->Diastolic_value, AP->Systolic_value_L, AP->Systolic_value_H);

			// fprintf (out, "%.2f %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n", AP->APD_count * BCL, AP->APD90, AP->APD75, AP->APD50, AP->APD20, AP->Vmax, AP->Vmin, AP->dVdtmax);
			if (AP->APD_out_swhich == 0)
			{
				AP->APD_out_end[0] = t - AP->timeAPDstart;
				AP->APD_out_swhich = 1;
			}
			else if (AP->APD_out_swhich == 1)
			{
				AP->APD_out_end[1] = t - AP->timeAPDstart;
				AP->APD_out_swhich = 0;
			}
		}
	}

	double dVdt = (Vm - AP->Vm_prev) / dt;
	if (dVdt > AP->dVdtmax)
		AP->dVdtmax = dVdt;
	AP->Istim_prev = Istim;
	AP->Vm_prev = Vm;
}

/*INa is used to measure the upstroke */
void APD90_DS_Value_measure_with_INa(double INa_threshold, double INa, double t,
		double Istim, double BCL, double dt, double Vm, double Value,
		AP_infor *AP, FILE *out)
{
	if (((Istim >= 1e-3 || Istim <= -1e-3) && (fabs(AP->Istim_prev) < 1e-10)))
	{
		if (AP->APD_count > 0)
		{
			// printf("dddddddddddd\n");
			// printf("%f %f %f\n", AP->Diastolic_value, AP->Systolic_value_L, AP->Systolic_value_H);

			fprintf(out,
					"%.2f %.3f %.3f %.3f %.3f %.3f %.3f  %.3f %.3f %.3f %.3f\n",
					AP->APD_count * BCL, BCL, AP->APD90, AP->APD75, AP->APD50,
					AP->APD20, AP->Systolic_value_L, AP->Systolic_value_H,
					AP->Vmax, AP->Vmin, AP->dVdtmax);

		}
		AP->APD_switch = 1;
		AP->Vmax = -60;
		AP->timeAPDstart = t;
		AP->Vmin = Vm;
		AP->dVdtmax = 0.0;
		AP->Diastolic_value = Value;
		AP->Systolic_value_L = Value;
		AP->Systolic_value_H = Value;
		AP->N_stim++;
		AP->Vmax_switch = 0;

		// printf("ssssssssss\n");
	}
	if (AP->Systolic_value_L > Value)
	{
		AP->Systolic_value_L = Value;
	}
	if (AP->Systolic_value_H < Value)
	{
		AP->Systolic_value_H = Value;
	}
	if (AP->APD_switch == 1)
	{
		if (Vm > -30.0)
		{
			if ((Vm > AP->Vmax) && (INa <= INa_threshold))
			{
				AP->Vmax = Vm;
			}
			else
				AP->Vmax_switch = 1;
			if (AP->Vmax_switch == 1)
			{
				AP->APD_switch = 2;
				AP->AMP_last = AP->AMP;
				AP->AMP = AP->Vmax - AP->Vmin;
				/*printf("%f %f\n", AP->AMP, AP->AMP_last);
				 printf("%f %f\n", AP->Vmax, AP->Vmin);*/

				AP->v90 = AP->Vmax - 0.9 * (AP->Vmax - AP->Vmin);
				AP->v75 = AP->Vmax - 0.75 * (AP->Vmax - AP->Vmin);
				AP->v50 = AP->Vmax - 0.50 * (AP->Vmax - AP->Vmin);
				AP->v20 = AP->Vmax - 0.20 * (AP->Vmax - AP->Vmin);
				// printf("asfde\n");
			}
		}
	}
	if (AP->APD_switch == 2)
	{

		if ((AP->Vm_prev >= AP->v20) && (Vm <= AP->v20))
		{
			AP->APD20 = t - AP->timeAPDstart;
			// printf("aaaa\n");

		}
		else if ((AP->Vm_prev >= AP->v50) && (Vm <= AP->v50))
		{
			AP->APD50 = t - AP->timeAPDstart;
			// printf("bbbbbbbbb\n");

		}
		else if (AP->Vm_prev >= AP->v75 && Vm <= AP->v75)
		{
			AP->APD75 = t - AP->timeAPDstart;
			// printf("ccccccccccc\n");

		}
		else if (AP->Vm_prev >= AP->v90 && Vm <= AP->v90)
		{
			AP->APD_switch = 0;
			AP->APD_count++;
			AP->Vmax_switch = 0;
			AP->APD90 = t - AP->timeAPDstart;
			if (AP->APD_out_swhich == 0)
			{
				AP->APD_out_end[0] = t - AP->timeAPDstart;
				AP->APD_out_swhich = 1;
			}
			else if (AP->APD_out_swhich == 1)
			{
				AP->APD_out_end[1] = t - AP->timeAPDstart;
				AP->APD_out_swhich = 0;
			}
		}
	}

	double dVdt = (Vm - AP->Vm_prev) / dt;
	if (dVdt > AP->dVdtmax)
		AP->dVdtmax = dVdt;
	AP->Istim_prev = Istim;
	AP->Vm_prev = Vm;
}

void APD90_DS_Value_measure_with_StrokeTime_threshold(double StrokeTime,
		double t, double Istim, double BCL, double dt, double Vm, double Value,
		AP_infor *AP, FILE *out)
{
	if (((Istim >= 1e-3 || Istim <= -1e-3) && (fabs(AP->Istim_prev) < 1e-10)))
	{
		AP->APD_switch = 1;
		AP->Vmax = -60;
		AP->timeAPDstart = t;
		AP->Vmin = Vm;
		AP->dVdtmax = 0.0;
		AP->Diastolic_value = Value;
		AP->Systolic_value_L = Value;
		AP->Systolic_value_H = Value;
		AP->N_stim++;
		AP->t_since_up_stroke = 0.0;
		AP->Vmax_switch = 0;

		// printf("ssssssssss\n");
	}
	if (AP->Systolic_value_L > Value)
	{
		AP->Systolic_value_L = Value;
	}
	if (AP->Systolic_value_H < Value)
	{
		AP->Systolic_value_H = Value;
	}
	AP->t_since_up_stroke += dt;

	if (AP->APD_switch == 1)
	{

		if (Vm > -10.0)
		{
			if ((Vm > AP->Vmax) && (AP->t_since_up_stroke <= StrokeTime))
			{
				AP->Vmax = Vm;
				// printf("%f\n", AP->t_since_up_stroke);
				// printf("%f\n", AP->Vmax);

			}
			else
				AP->Vmax_switch = 1;
			if (AP->Vmax_switch == 1)
			{
				AP->APD_switch = 2;
				AP->AMP_last = AP->AMP;
				AP->AMP = AP->Vmax - AP->Vmin;
				// printf("%f %f\n", AP->AMP, AP->AMP_last);
				// printf("%f %f\n", AP->Vmax, AP->Vmin);

				AP->v90 = AP->Vmax - 0.9 * (AP->Vmax - AP->Vmin);
				AP->v75 = AP->Vmax - 0.75 * (AP->Vmax - AP->Vmin);
				AP->v50 = AP->Vmax - 0.50 * (AP->Vmax - AP->Vmin);
				AP->v20 = AP->Vmax - 0.20 * (AP->Vmax - AP->Vmin);
				// printf("asfde\n");
			}
		}
	}
	if (AP->APD_switch == 2)
	{

		if ((AP->Vm_prev >= AP->v20) && (Vm <= AP->v20))
		{
			AP->APD20 = t - AP->timeAPDstart;
			// printf("aaaa\n");

		}
		else if ((AP->Vm_prev >= AP->v50) && (Vm <= AP->v50))
		{
			AP->APD50 = t - AP->timeAPDstart;
			// printf("bbbbbbbbb\n");

		}
		else if (AP->Vm_prev >= AP->v75 && Vm <= AP->v75)
		{
			AP->APD75 = t - AP->timeAPDstart;
			// printf("ccccccccccc\n");

		}
		else if (AP->Vm_prev >= AP->v90 && Vm <= AP->v90)
		{
			AP->APD_switch = 0;
			AP->APD_count++;
			AP->Vmax_switch = 0;
			AP->APD90 = t - AP->timeAPDstart;
			// printf("dddddddddddd\n");
			// printf("%f %f %f\n", AP->Diastolic_value, AP->Systolic_value_L, AP->Systolic_value_H);

			fprintf(out, "%.2f %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n",
					AP->APD_count * BCL, AP->APD90, AP->APD75, AP->APD50,
					AP->APD20, AP->Vmax, AP->Vmin, AP->dVdtmax);
			if (AP->APD_out_swhich == 0)
			{
				AP->APD_out_end[0] = t - AP->timeAPDstart;
				AP->APD_out_swhich = 1;
			}
			else if (AP->APD_out_swhich == 1)
			{
				AP->APD_out_end[1] = t - AP->timeAPDstart;
				AP->APD_out_swhich = 0;
			}
		}
	}

	double dVdt = (Vm - AP->Vm_prev) / dt;
	if (dVdt > AP->dVdtmax)
		AP->dVdtmax = dVdt;
	AP->Istim_prev = Istim;
	AP->Vm_prev = Vm;
}

/* measure the accumulations of calcium build up from the ICaL and RyR */
void calcium_accumulation(double t, double Istim, double BCL, double dt,
		double ICaL_con, double J_RyR, AP_infor *AP, FILE *out)
{
	if ((Istim >= 1e-3 || Istim <= -1e-3) && fabs(AP->Istim_prev) < 1e-3)
	{
		// printf("%s\n", );
		fprintf(out, "%.3f %f %f\n", t, AP->ICaL_in, AP->RyR_in);
		AP->timeAPDstart = t;
		AP->ICaL_in = 0.0;
		AP->RyR_in = 0.0;
	}

	AP->ICaL_in += ICaL_con * dt;
	AP->RyR_in += J_RyR * dt;
	AP->Istim_prev = Istim;
}
