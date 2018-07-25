
# include <stdio.h>

int main(){
// variable we are going to calculate on every 0.1 second 
double t = 0.0; //time 
double s  = 0.0;//height
double v = 0.0;//velocity 
double a = 0.0;//acceleration 
double F = 0.0;//total force 
double m = 0.0340 + 0.0242;//  mass = 0.0340 + 0.0242 = 0.0582 

double g = 9.80665;
double Fg = m * g;
//int i = 0; 

double Fd_body; 
double Fd_fins;
double Ft;
double dv;
double ds;
double cd_body = 0.45;
double rho = 1.293;
double A_body = 0.000506; 
double cd_fins = 0.01;

double A_fins = 0.00496; 



do{ 
  
  double dt = 0.1;
    t = t + dt;
      
  if (t == 0.1) 
     Ft = 6.0;
  if (t== 0.2) 
     Ft = 14.1;
  if (t == 0.3)
     Ft = 5.0;
  if (t >= 0.4 && t < 1.9)  
     Ft = 4.5;
  if (t >= 1.9) 
      Ft = 0;
      
  m = m - 0.0001644 * Ft; 
  Fg = m * g;
  Fd_body =  (cd_body * rho * A_body * v * v)/2; 
    
                              
  Fd_fins =  (cd_fins * rho * A_fins * v * v)/2; 
      

  F = Ft - ( Fd_body + Fd_fins + Fg);
  a = F/m;
  dv = a * dt;
  v = v + dv;
  ds = v * dt; 
  s = s + ds;
  //   t = t + dt;
  // updating variables 
 
      

  printf ("%s%f%s\n", "t =",   t, " sec");
  printf ("%s%f%s\n", "s =",   s, " meter");
  printf ("%s%f%s\n", "v=",    v, " meter per second");
  printf ("%s%f%s\n", "a =",   a, " meter per second squared");
  printf ("%s%f%s\n", "m =" ,  m, " kilogram");
 printf ("\n \n"  );
 }while (v > 0);
 }

