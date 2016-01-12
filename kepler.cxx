# include <iostream>
# include <cmath>
# include <fstream>

using namespace std;

void Euler(double* p1, double* p2, double* q1, double*q2, const double dt,const double N){
  for(int i=0;i<N-1;i++){
    p1[i+1] = p1[i] - dt*q1[i]*pow(pow(q1[i],2.0)+pow(q2[i],2.0),-1.5);
    p2[i+1] = p2[i] - dt*q2[i]*pow(pow(q1[i],2.0)+pow(q2[i],2.0),-1.5);
    q1[i+1] = q1[i] + dt*p1[i+1];
    q2[i+1] = q2[i] + dt*p2[i+1];   
  }}
void Hamilton(double* p1, double* p2, double* q1, double*q2,double* H,const double N){
  for(int i=0;i<N;i++){
   H[i] = 0.5 * (p1[i] * p1[i] + p2[i] * p2[i]) - 1.0/sqrt(q1[i]*q1[i] + q2[i] *q2[i] ); 
  }
}
void Speichern(double*p1,double*p2, double*q1, double*q2, double* H, const double dt,const double N){
  ofstream out ("Euler_dt_0.0005.txt");
  for(int i=0; i<N; i++){
    out << dt*i << "\t" << p1[i] << "\t" << p2[i] << "\t" << q1[i] << "\t" << q2[i] << "\t" << H[i] <<  endl;
  }
  out.close(); 
}

int main(void){
  const double e = 0.6;
  const double dt=0.0005;	
  const double tanf=0.0,tend=20*M_PI;
  const int N= ((tend-tanf)/dt) +1;
  double q1[N],q2[N],p1[N],p2[N],H[N];
  q1[0] = 1-e, q2[0] = 0.0, p1[0] = 0.0, p2[0] = sqrt((1+e)/(1-e));
  Euler(p1,p2,q1,q2,dt,N);
  Hamilton(p1,p2,q1,q2,H,N);
  Speichern(p1,p2,q1,q2,H,dt,N);
  
  return 0;
}

