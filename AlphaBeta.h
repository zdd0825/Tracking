#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 


///////////////////¸Ä//////////////////////////////////
typedef struct {
	float alpha; //alpha value (effects x, eg pos)
	float beta; //beta value (effects v, eg vel)
	float xk_1; //current x-estimate 
	float vk_1; //current v-estimate 
	float x_mea[100];
} AlphaBeta;

void InitializeAlphaBeta(float x_measured, float alpha, float beta, AlphaBeta* pab) {
	pab->xk_1 = x_measured;
	pab->vk_1 = 0;
	pab->alpha = alpha;
	pab->beta = beta;
}

float Cal_var(float *x_mea, int n) {
	float ave = 0.;
	float ret = 0.;
	for (int i = 0; i<n; i++)
		ave += x_mea[i];

	ave /= n;

	for (int i = 0; i<n; i++)
		ret += (x_mea[i] - ave)*(x_mea[i] - ave);

	ret /= n;
	return ret;
}


int AlphaBetaFilter(float x_measured, float dt, AlphaBeta* pab, int flag) {
	float xk_1 = pab->xk_1;
	float vk_1 = pab->vk_1;
	float alpha = pab->alpha;
	float beta = pab->beta;

	float xk; //current system state (ie: position)
	float vk; //derivative of system state (ie: velocity)
	float rk; //residual error 

	//update our (estimated) state 'x' from the system (ie pos = pos + vel (last).dt)
	xk = xk_1 + dt * vk_1;
	//printf("################################################################ %6.3f\n", xk);
	//update (estimated) velocity  
	vk = vk_1;
	//what is our residual error (mesured - estimated) 
	rk = x_measured - xk;
	//update our estimates given the residual error. 
	xk = xk + alpha * rk;
	vk = vk + beta / dt * rk;
	//finished! 

	//now all our "currents" become our "olds" for next time 
	pab->vk_1 = vk;
	pab->xk_1 = xk;
	//printf("################################################################ %6.3f\n", xk);

	pab->x_mea[flag % 100] = x_measured;
	float sigma_v = Cal_var(pab->x_mea, flag<100 ? (flag + 1) : 100);
	float sigma_w = 2.0 * dt;

	if (sigma_v != 0){
		float lamb = (sigma_w*dt*dt) / sigma_v;
		pab->alpha = 1 / 8.0*(-lamb*lamb - 8.0*lamb + (lamb + 4)*sqrt(lamb*lamb + 8.0*lamb));
		pab->beta = 1 / 4.0*(lamb*lamb + 4.0*lamb - lamb*sqrt(lamb*lamb + 8.0*lamb));
	}

	return flag+1;
}



// void main()
// {
// AlphaBeta ab_x;
// InitializeAlphaBeta(1, 0.85, 0.001, &ab_x); //x position
// //double error = (2 * 0.271*0.271 + 2 * 0.0285 + 0.271*0.0285) / 0.271*(4 - 2 * 0.271 - 0.0285) * 0.01;
// //printf("error = %f\n", error);
// float m[] = { 0, 20, 25, 16, 28, 21, 17, 19, 20, 23,26,30,25,16,20,22,14,24 };
// float dis = 0.;
// int flag = 0;
// for (int i = 0; i<18; i++)
// {
// AlphaBetaFilter(m[i], 1, &ab_x,flag);
// printf("%6.3f\n", ab_x.xk_1);
// //printf("alpha: %6.3f\n", ab_x.alpha);
// //printf("beta: %6.3f\n", ab_x.beta);
// dis += (20.0 - ab_x.xk_1) *  (20.0 - ab_x.xk_1);
// }
// printf("²î¾à£º%6.3f\n", dis);
// }