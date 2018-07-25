/* least_square_fit.c  tailorable code, provide your input and setup */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#undef  abs
#define abs(x) (((x)<0.0)?(-(x)):(x))
#undef  max
#define max(x,y) (((x)>(y))?((x)):(y))
#undef  min
#define min(x,y) (((x)<(y))?((x)):(y))

static void simeq(int n, double A[], double Y[], double X[]);
static FILE *inp = NULL;

/*
 * The purpose of this package is to provide a reliable and convenient
 * means for fitting existing data by a few coefficients. The companion
 * package check_fit provides the means to use the coefficients for
 * interpolation and limited extrapolation.
 *
 * This package implements the least square fit. 
 *
 * The problem is stated as follows :
 *   Given measured data for values of Y based on values of X1,X2 and X3. e.g.
 *
 *    Y_actual         X1      X2     X3  
 *    --------       -----   -----  -----
 *     32.5           1.0     2.5    3.7 
 *      7.2           2.0     2.5    3.6 
 *      6.9           3.0     2.7    3.5 
 *     22.4           2.2     2.1    3.1 
 *     10.4           1.5     2.0    2.6 
 *     11.3           1.6     2.0    3.1 
 *
 *  Find a, b and c such that   Y_approximate =  a * X1 + b * X2 + c * X3
 *  and such that the sum of (Y_actual - Y_approximate) squared is minimized.
 *
 * The method for determining the coefficients a, b and c follows directly
 * form the problem definition and mathematical analysis. (See more below)
 *
 * Y is called the dependent variable and X1 .. Xn the independent variables.
 * The procedures below implements a few special cases and the general case.
 *    The number of independent variables can vary.
 *    The approximation equation may use powers of the independent variables
 *    The user may create additional independent variables e.g. X2 = SIN(X1)
 *    with the restriction that the independent variables are linearly
 *    independent.  e.g.  Xi not equal  p Xj + q  for all i,j,p,q
 *
 *
 *
 * The mathematical derivation of the least square fit is as follows :
 *
 * Given data for the independent variable Y in terms of the dependent
 * variables S,T,U and V  consider that there exists a function F
 * such that     Y = F(S,T,U,V)
 * The problem is to find coefficients a,b,c and d such that
 *            Y_approximate = a * S + b * T + c * U + d * V
 * and such that the sum of ( Y - Y_approximate ) squared is minimized.
 *
 * Note: a, b, c, d are scalars. S, T, U, V, Y, Y_approximate are vectors.
 *
 * To find the minimum of  SUM( Y - Y_approximate ) ** 2
 * the derivatives must be taken with respect to a,b,c and d and
 * all must equal zero simultaneously. The steps follow :
 *
 *  SUM( Y - Y_approximate ) ** 2 = SUM( Y - a*S - b*T - c*U - d*V ) ** 2
 *
 * d/da =  -2 * S * SUM( Y - A*S - B*T - C*U - D*V )
 * d/db =  -2 * T * SUM( Y - A*S - B*T - C*U - D*V )
 * d/dc =  -2 * U * SUM( Y - A*S - B*T - C*U - D*V )
 * d/dd =  -2 * V * SUM( Y - A*S - B*T - C*U - D*V )
 *
 * Setting each of the above equal to zero and putting constant term on left
 *    the -2 is factored out,
 *    the independent variable is moved inside the summation
 *
 *  SUM( a*S*S + b*S*T + c*S*U + d*S*V = S*Y )
 *  SUM( a*T*S + b*T*T + c*T*U + d*T*V = T*Y )
 *  SUM( a*U*S + b*U*T + c*U*U + d*U*V = U*Y )
 *  SUM( a*V*S + b*V*T + c*V*U + d*V*V = V*Y )
 *
 * Distributing the SUM inside yields
 *
 *  a * SUM(S*S) + b * SUM(S*T) + c * SUM(S*U) + d * SUM(S*V) = SUM(S*Y)
 *  a * SUM(T*S) + b * SUM(T*T) + c * SUM(T*U) + d * SUM(T*V) = SUM(T*Y)
 *  a * SUM(U*S) + b * SUM(U*T) + c * SUM(U*U) + d * SUM(U*V) = SUM(U*Y)
 *  a * SUM(V*S) + b * SUM(V*T) + c * SUM(V*U) + d * SUM(V*V) = SUM(V*Y)
 *
 * To find the coefficients a,b,c and d solve the linear system of equations
 *
 *    | SUM(S*S)  SUM(S*T)  SUM(S*U)  SUM(S*V) |   | a |   | SUM(S*Y) |
 *    | SUM(T*S)  SUM(T*T)  SUM(T*U)  SUM(T*V) | x | b | = | SUM(T*Y) |
 *    | SUM(U*S)  SUM(U*T)  SUM(U*U)  SUM(U*V) |   | c |   | SUM(U*Y) |
 *    | SUM(V*S)  SUM(V*T)  SUM(V*U)  SUM(V*V) |   | d |   | SUM(V*Y) |
 *
 * Some observations :
 *     S,T,U and V must be linearly independent.
 *     There must be more data sets (Y, S, T, U, V) than variables.
 *     The analysis did not depend on the number of independent variables
 *     A polynomial fit results from the substitutions S=1, T=X, U=X**2, V=X**3
 *     The general case for any order polynomial follows, fit_pn.
 *     Any substitution such as three variables to various powers may be used.
 */


double time [20] = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9}; 
double thrust [20] = { 6.0, 14.1, 5.0, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 0.0}; 

static int data_set(double *y, double *x) /* returns 1 for data, 0 for end  */
{                                         /* sets value of y for value of x */
  //int k = 100; /* 100 */
  int k = 20; /* 20 */
  static i = 0;
  double t1 = 0.0;
  double dt1 = 0.1; /* approx Pi */
  double xx;
  double yy;

  i++;  
  if(i>k) {i=0; return 0;} /* ready for check */
  xx = t1 + (double)i * dt1;
xx = time[i];
  yy = thrust[i];
  *x = xx;
  *y = yy;
  return 1;
} /* end data_set */

static void fit_pn(int n, double A[], double Y[], double C[])
{                   /* n is number of coefficients, highest power-1 */
  int i, j, k;
  double x, y, t;
  double pwr[40]; /* at least n */
  
  for(i=0; i<n; i++)
  {
    for(j=0; j<n; j++)
    {
      A[i*n+j] = 0.0;
    }
    Y[i] = 0.0;
  }
  while(data_set(&y, &x))
  {
    pwr[0] = 1.0;
    for(i=1; i<n; i++) pwr[i] = pwr[i-1]*x;
    for(i=0; i<n; i++)
    {
      for(j=0; j<n; j++)
      {
        A[i*n+j] = A[i*n+j] + pwr[i]*pwr[j];
      }
      Y[i] = Y[i] + y*pwr[i];
    }
  }
  simeq(n, A, Y, C);
  for(i=0; i<n; i++) printf("C[%d]=%g \n", i, C[i]);
} /* end fit_pn */

static void check_pn(int n, double C[],
                     double *rms_err, double *avg_err, double *max_err)
{
  double x, y, ya, diff;
  double sumsq = 0.0;
  double sum = 0.0;
  double maxe = 0.0;
  double xmin, xmax, ymin, ymax, xbad, ybad;
  int i, k, imax;
  
  k = 0;
  while(data_set(&y, &x))
  {
    if(k==0)
    {
      xmin=x;
      xmax=x;
      ymin=y;
      ymax=y;
      imax=0;
      xbad=x;
      ybad=y;
    }
    if(x>xmax) xmax=x;
    if(x<xmin) xmin=x;
    if(y>ymax) ymax=y;
    if(y<ymin) ymin=y;
    k++;
    ya = C[n-1]*x;
    for(i=n-2; i>0; i--)
    {
      ya = (C[i]+ya)*x;
    }
    ya = ya + C[0];
    diff = abs(y-ya);
    if(diff>maxe)
    {
      maxe=diff;
      imax=k;
      xbad=x;
      ybad=y;
    }
    sum = sum + diff;
    sumsq = sumsq + diff*diff;
  }
  printf("check_pn k=%d, xmin=%g, xmax=%g, ymin=%g, ymax=%g \n",
          k, xmin, xmax, ymin, ymax);
  *max_err = maxe;
  *avg_err = sum/(double)k;
  *rms_err = sqrt(sumsq/(double)k);
  printf("max=%g at %d, xbad=%g, ybad=%g\n", maxe, imax, xbad, ybad);
} /* end check_pn */

static void fit_file(int n, double A[], double Y[], double C[])
{
  int i, j, k;
  double x, y, t;
  double pwr[30]; /* at least 2n */
  
  for(i=0; i<n; i++)
  {
    for(j=0; j<n; j++)
    {
      A[i*n+j] = 0.0;
    }
    Y[i] = 0.0;
  }
  while(!feof(inp))
  {
    fscanf(inp,"%lf %lf\n", &x, &y);
    pwr[0] = 1.0;
    for(i=1; i<n; i++) pwr[i] = pwr[i-1]*x;
    for(i=0; i<n; i++)
    {
      for(j=0; j<n; j++)
      {
        A[i*n+j] = A[i*n+j] + pwr[i]*pwr[j];
      }
      Y[i] = Y[i] + y*pwr[i];
    }
  }
  rewind(inp);
  simeq(n, A, Y, C);
  for(i=0; i<n; i++) printf("C[%d]=%g \n", i, C[i]);
} /* end fit_file */

static void check_file(int n, double C[],
                     double *rms_err, double *avg_err, double *max_err)
{
  double x, y, ya, diff;
  double sumsq = 0.0;
  double sum = 0.0;
  double maxe = 0.0;
  double xmin, xmax, ymin, ymax, xbad, ybad;
  int i, k, imax;
  
  k = 0;
  while(!feof(inp))
  {
    fscanf(inp,"%lf %lf\n", &x, &y);
    if(k==0)
    {
      xmin=x;
      xmax=x;
      ymin=y;
      ymax=y;
      imax=0;
      xbad=x;
      ybad=y;
    }
    if(x>xmax) xmax=x;
    if(x<xmin) xmin=x;
    if(y>ymax) ymax=y;
    if(y<ymin) ymin=y;
    k++;
    ya = C[n-1]*x;
    for(i=n-2; i>0; i--)
    {
      ya = (C[i]+ya)*x;
    }
    ya = ya + C[0];
    diff = abs(y-ya);
    if(diff>maxe)
    {
      maxe=diff;
      imax=k;
      xbad=x;
      ybad=y;
    }
    sum = sum + diff;
    sumsq = sumsq + diff*diff;
  }
  rewind(inp);
  printf("check_pn k=%d, xmin=%g, xmax=%g, ymin=%g, ymax=%g \n",
          k, xmin, xmax, ymin, ymax);
  *max_err = maxe;
  *avg_err = sum/(double)k;
  *rms_err = sqrt(sumsq/(double)k);
  printf("max=%g at %d, xbad=%g, ybad=%g\n", maxe, imax, xbad, ybad);
} /* end check_pn */

static void simeq(int n, double A[], double Y[], double X[])
{

/*      PURPOSE : SOLVE THE LINEAR SYSTEM OF EQUATIONS WITH REAL     */
/*                COEFFICIENTS   [A] * |X| = |Y|                     */
/*                                                                   */
/*      INPUT  : THE NUMBER OF EQUATIONS  n                          */
/*               THE REAL MATRIX  A   should be A[i][j] but A[i*n+j] */
/*               THE REAL VECTOR  Y                                  */
/*      OUTPUT : THE REAL VECTOR  X                                  */
/*                                                                   */
/*      METHOD : GAUSS-JORDAN ELIMINATION USING MAXIMUM ELEMENT      */
/*               FOR PIVOT.                                          */
/*                                                                   */
/*      USAGE  :     simeq(n,A,Y,X);                                 */
/*                                                                   */
/*                                                                   */
/*    WRITTEN BY : JON SQUIRE , 28 MAY 1983                          */
/*    ORIGINAL DEC 1959 for IBM 650, TRANSLATED TO OTHER LANGUAGES   */
/*    e.g. FORTRAN converted to Ada converted to C                   */

    double *B;           /* [n][n+1]  WORKING MATRIX */
    int *ROW;            /* ROW INTERCHANGE INDICES */
    int HOLD , I_PIVOT;  /* PIVOT INDICES */
    double PIVOT;        /* PIVOT ELEMENT VALUE */
    double ABS_PIVOT;
    int i,j,k,m;

    B = (double *)calloc((n+1)*(n+1), sizeof(double));
    ROW = (int *)calloc(n, sizeof(int));
    m = n+1;

    /* BUILD WORKING DATA STRUCTURE */
    for(i=0; i<n; i++){
      for(j=0; j<n; j++){
        B[i*m+j] = A[i*n+j];
      }
      B[i*m+n] = Y[i];
    }
    /* SET UP ROW  INTERCHANGE VECTORS */
    for(k=0; k<n; k++){
      ROW[k] = k;
    }

    /* BEGIN MAIN REDUCTION LOOP */
    for(k=0; k<n; k++){

      /* FIND LARGEST ELEMENT FOR PIVOT */
      PIVOT = B[ROW[k]*m+k];
      ABS_PIVOT = abs(PIVOT);
      I_PIVOT = k;
      for(i=k; i<n; i++){
        if( abs(B[ROW[i]*m+k]) > ABS_PIVOT){
          I_PIVOT = i;
          PIVOT = B[ROW[i]*m+k];
          ABS_PIVOT = abs ( PIVOT );
        }
      }

      /* HAVE PIVOT, INTERCHANGE ROW POINTERS */
      HOLD = ROW[k];
      ROW[k] = ROW[I_PIVOT];
      ROW[I_PIVOT] = HOLD;

      /* CHECK FOR NEAR SINGULAR */
      if( ABS_PIVOT < 1.0E-10 ){
        for(j=k+1; j<n+1; j++){
          B[ROW[k]*m+j] = 0.0;
        }
        printf("redundant row (singular) %d \n", ROW[k]);
      } /* singular, delete row */
      else{

        /* REDUCE ABOUT PIVOT */
        for(j=k+1; j<n+1; j++){
          B[ROW[k]*m+j] = B[ROW[k]*m+j] / B[ROW[k]*m+k];
        }

        /* INNER REDUCTION LOOP */
        for(i=0; i<n; i++){
          if( i != k){
            for(j=k+1; j<n+1; j++){
              B[ROW[i]*m+j] = B[ROW[i]*m+j] - B[ROW[i]*m+k] * B[ROW[k]*m+j];
            }
          }
        }
      }
      /* FINISHED INNER REDUCTION */
    }

    /* END OF MAIN REDUCTION LOOP */
    /* BUILD  X  FOR RETURN, UNSCRAMBLING ROWS */
    for(i=0; i<n; i++){
      X[i] = B[ROW[i]*m+n];
    }
    free(B);
    free(ROW);
} /* end simeq */

int main(int argc, char *argv[])
{
  int n;
  double A[400];
  double C[20];
  double Y[20];
  double rms_err, avg_err, max_err;
  
  printf("least_square_fit.c\n");
  if(argc>1)
  {
    printf("trying to open %s \n",argv[1]);
    inp = fopen(argv[1],"r");
  }
  
  /* sample polynomial least square fit, 3th power */
  n=3+1; /* need constant term and five powers 1,2,3,4,5 */
  printf("fitting upthrust vs time curve 20 points 0.0 to 1.9, %d degree polynomial\n", n-1);
  fit_pn(n, A, Y, C);
  check_pn(n, C, &rms_err, &avg_err, &max_err);
  printf("rms_err=%g, avg_err=%g, max_err=%g \n\n",
         rms_err, avg_err, max_err);
  n=4+1; /* need constant term and five powers 1,2,3,4,5 */
  printf("fitting upthrust vs time curve 20 points 0.0 to 1.9, %d degree polynomial\n", n-1);
  fit_pn(n, A, Y, C);
  check_pn(n, C, &rms_err, &avg_err, &max_err);
  printf("rms_err=%g, avg_err=%g, max_err=%g \n\n",
         rms_err, avg_err, max_err);
  n=5+1; /* need constant term and five powers 1,2,3,4,5 */
  printf("fitting upthrust vs time curve 20 points 0.0 to 1.9, %d degree polynomial\n", n-1);
  fit_pn(n, A, Y, C);
  check_pn(n, C, &rms_err, &avg_err, &max_err);
  printf("rms_err=%g, avg_err=%g, max_err=%g \n\n",
         rms_err, avg_err, max_err);
  n=6+1;
  printf("fitting uthrust vs time curve  20 points 0.0 to 1.0, %d degree polynomial\n", n-1);
  fit_pn(n, A, Y, C);
  check_pn(n, C, &rms_err, &avg_err, &max_err);
  printf("rms_err=%g, avg_err=%g, max_err=%g \n\n",
         rms_err, avg_err, max_err);
  n=7+1;
  printf("fit exp(x)*sin(x) 100 points 0.0 to Pi, %d degree polynomial\n", n-1);
  fit_pn(n, A, Y, C);
  check_pn(n, C, &rms_err, &avg_err, &max_err);
  printf("rms_err=%g, avg_err=%g, max_err=%g \n",
         rms_err, avg_err, max_err);
  if(inp!=NULL)
  {
    n = 8;
    printf("\n\nfitting file %s, n=%d \n", argv[1], n);
    fit_file(n, A, Y, C);
    check_file(n, C, &rms_err, &avg_err, &max_err);
    printf("rms_err=%g, avg_err=%g, max_err=%g \n",
           rms_err, avg_err, max_err);
    n = 9;
    printf("\n\nfitting file %s, n=%d \n", argv[1], n);
    fit_file(n, A, Y, C);
    check_file(n, C, &rms_err, &avg_err, &max_err);
    printf("rms_err=%g, avg_err=%g, max_err=%g \n",
           rms_err, avg_err, max_err);
    n = 10;
    printf("\n\nfitting file %s, n=%d \n", argv[1], n);
    fit_file(n, A, Y, C);
    check_file(n, C, &rms_err, &avg_err, &max_err);
    printf("rms_err=%g, avg_err=%g, max_err=%g \n",
           rms_err, avg_err, max_err);
    n = 11;
    printf("\n\nfitting file %s, n=%d \n", argv[1], n);
    fit_file(n, A, Y, C);
    check_file(n, C, &rms_err, &avg_err, &max_err);
    printf("rms_err=%g, avg_err=%g, max_err=%g \n",
           rms_err, avg_err, max_err);
    n = 12;
    printf("\n\nfitting file %s, n=%d \n", argv[1], n);
    fit_file(n, A, Y, C);
    check_file(n, C, &rms_err, &avg_err, &max_err);
    printf("rms_err=%g, avg_err=%g, max_err=%g \n",
           rms_err, avg_err, max_err);
    fclose(inp);
  }
  return 0;
} /* end main for least_square_fit.c */
