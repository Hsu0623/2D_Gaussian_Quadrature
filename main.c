#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <math.h>

#define N 9
#define K 10
#define pi acos(-1)

//integral range = [-1.0, 1.0]*[-1.0, 1.0].
double  a = -1.0, b = 1.0;

//actual intergration 
double actual_sum = 0.160429671298544;

// Sample points X or Y in [-1,1], 2-, 3- ,4-, ....., 10- points. 
double p[N][K] = { {-0.5773502691896257, 0.5773502691896257, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                   {0.0, -0.7745966692414834, 0.7745966692414834, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                   {-0.3399810435848563, 0.3399810435848563, -0.8611363115940526, 0.8611363115940526, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                   {0.0, -0.5384693101056831, 0.5384693101056831, -0.9061798459386640, 0.9061798459386640, 0.0, 0.0, 0.0, 0.0, 0.0},
                   {0.6612093864662645, -0.6612093864662645, -0.2386191860831969, 0.2386191860831969, -0.9324695142031521, 0.9324695142031521, 0.0, 0.0, 0.0, 0.0},
                   {0.0, 0.4058451513773972, -0.4058451513773972, -0.7415311855993945, 0.7415311855993945, -0.9491079123427585, 0.9491079123427585, 0.0, 0.0, 0.0},
                   {-0.1834346424956498, 0.1834346424956498, -0.5255324099163290, 0.5255324099163290, -0.7966664774136267, 0.7966664774136267, -0.9602898564975363, 0.9602898564975363, 0.0, 0.0},
                   {0.0, -0.8360311073266358, 0.8360311073266358, -0.9681602395076261, 0.9681602395076261, -0.3242534234038089, 0.3242534234038089, -0.6133714327005904, 0.6133714327005904, 0.0},
                   {-0.1488743389816312, 0.1488743389816312, -0.4333953941292472, 0.4333953941292472, -0.6794095682990244, 0.6794095682990244, -0.8650633666889845, 0.8650633666889845, -0.9739065285171717, 0.9739065285171717}};


//weights Aj, Ai, for 2-, 3-, ...., 10-points.
double wgt[N][K] = { {1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                     {0.8888888888888888, 0.5555555555555556, 0.5555555555555556, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                     {0.6521451548625461, 0.6521451548625461, 0.3478548451374538, 0.3478548451374538, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                     {0.5688888888888889, 0.4786286704993665, 0.4786286704993665,0.2369268850561891,0.2369268850561891, 0.0, 0.0, 0.0, 0.0, 0.0},
                     {0.3607615730481386, 0.3607615730481386, 0.4679139345726910, 0.4679139345726910, 0.1713244923791704, 0.1713244923791704, 0.0, 0.0, 0.0, 0.0},
                     {0.4179591836734694, 0.3818300505051189, 0.3818300505051189, 0.2797053914892766, 0.2797053914892766, 0.1294849661688697, 0.1294849661688697, 0.0, 0.0, 0.0},
                     {0.3626837833783620, 0.3626837833783620, 0.3137066458778873, 0.3137066458778873, 0.2223810344533745, 0.2223810344533745, 0.1012285362903763, 0.1012285362903763, 0.0, 0.0},
                     {0.3302393550012598, 0.1806481606948574, 0.1806481606948574, 0.0812743883615744, 0.0812743883615744, 0.3123470770400029, 0.3123470770400029, 0.2606106964029354, 0.2606106964029354, 0.0},
                     {0.2955242247147529, 0.2955242247147529, 0.2692667193099963, 0.2692667193099963, 0.2190863625159820, 0.2190863625159820, 0.1494513491505806, 0.1494513491505806, 0.0666713443086881, 0.0666713443086881}};


/*------------------------------------------------------------------------------------------------------
 * The function to be integrated: f(x, y) = (sin(2*pi*x)/2*pi*x) * (sin(3*pi*y)/3*pi*y).
 * The indefinite integral I(f(x)=F(x) =  arctan(x)
 */
double  my_func(double x, double y)
{
    double t0, t1, t2, t3;
    t0 = 2 * pi * x;
    t1 = 3 * pi * y;

    if (t0 == 0.0)
        t2 = 1.0;
    else
        t2 = sin(t0) / t0;

    if (t1 == 0.0)
        t3 = 1.0;
    else
        t3 = sin(t1) / t1;

    return t2 * t3;
}


/*-----------------------------------------------------------------------------------------------------
 * Conformation mapping from [-1, 1] to [a, b]
 *  Input:
 *     t: Gaussian point in [-1, 1].
 *     a, b: the 2 ends of the interval.
 *  Return:
 *     Gaussian point in [a, b].
 *  Called before evaluating function values.
 */
double  conform_map(double t, double a, double b)
{
    double x;
    x = (b - a) * t / 2.0 + (a + b) / 2.0;
    return (x);
}


/*---------------------------------------------------------------------------------------------------------
 * Gaussian quadrature integration subroutine.
 * Input:
 *   numInterval: numner of intervela(n*n)
 *   k: number of Gaussian points used(but actually k*k points),
 *   a_x, b_x: 2 ends of the X interval,
 *   a_y, b_y: 2 ends of the Y interval,
 *   f(): the funtion.
 * Return:
 *   The integral value.
 */
double  gauss_quadrature(int numInterval, int k, double a_x, double b_x, double a_y, double b_y, double f())
{
    int i, j;
    double sum = 0.0;
    double x, y, t_x, t_y;

    for (i = 0; i < k; i++) {
        t_x = p[k - 2][i];
        for (j = 0; j < k; j++) {
            t_y = p[k - 2][j];
            x = conform_map(t_x, a_x, b_x);
            y = conform_map(t_y, a_y, b_y);
            sum += wgt[k - 2][i] * wgt[k - 2][j] * f(x, y);
        }
    }
    numInterval = numInterval * numInterval;
    sum = sum / (double)numInterval;
    return sum;
}



/*-----------------------------------------------------------------------
 * computed relative_error;
 * Input:
 *  computed_sum: computed_sum.
 * Return:
 *  relative_error.
*/
double relative_error(double computed_sum) {
    double abs_err = fabs(computed_sum - actual_sum);
    return abs_err / actual_sum;
}


void Q5D()
{

    // numInterval is number of meshes
    int k, numInterval = 3, i, j;
    char filename[] = "Q4D_numinterva3.txt";
    //char filename[] = "Q4D_numinterval1.plt";
    
    double gaussSum, localSum, relative_err;
    double h;
    FILE* fp;

    fp = fopen(filename, "w");
    //using different numbers of intervals(numInterval*numInterval)
    fprintf(fp, "No. of interval = %d*%d\n", numInterval, numInterval);
    h = (b - a) / numInterval;
    // using different number of Gaussian points(k*k).
    for (k = 2; k <= 10; k++) {
        fprintf(fp, "   No. of Gaussian points(k*k) = %d*%d\n", k, k);
        gaussSum = 0.0;
        for (i = 1; i <= numInterval; i++) {
            for (j = 1; j <= numInterval; j++) {
                localSum = gauss_quadrature(numInterval, k, a + (j - 1) * h, a + j * h, a + (i - 1) * h, a + i * h, my_func);
                gaussSum += localSum;
            }
        }
        relative_err = relative_error(gaussSum);
        fprintf(fp, "       Integral= %15.15f\n", gaussSum);
        fprintf(fp, "       relative_err= %15.15f\n", relative_err);
        //fprintf(fp, "%d %15.15f\n", k, relative_err);
    }
    fclose(fp);

}


void Q5C()
{
    int numInterval, i, j;
    double gaussSum, localSum, relative_err;
    double h;
    FILE* fp;

    char filename[] = "Q4C_node1.plt";
    //char filename[] = "Q4C_node2.txt";
    fp = fopen(filename, "w");
    // k is number of sample node(k*k)
    int k = 1;

    // using different number of Gaussian points(k*k).
    //fprintf(fp, "No. of Gaussian points(k*k) = %d*%d\n", k, k);
    // using different numbers of intervals(numInterval*numInterval)
    for (numInterval = 1; numInterval <= 8; numInterval++) {
        gaussSum = 0.0;
        //fprintf(fp, "   No. of interval = %d*%d\n", numInterval, numInterval);
        h = (b - a) / numInterval;
        for (i = 1; i <= numInterval; i++) {
            for (j = 1; j <= numInterval; j++) {
                localSum = gauss_quadrature(numInterval, k, a + (j - 1) * h, a + j * h, a + (i - 1) * h, a + i * h, my_func);
                gaussSum += localSum;
            }
        }
        relative_err = relative_error(gaussSum);
        //fprintf(fp, "       Integral= %15.15f\n", gaussSum);
        //fprintf(fp, "       relative_err= %15.15f\n", relative_err);
        fprintf(fp, "%d %15.15f\n", numInterval, relative_err);
    }

    fclose(fp);
}
void Q4()
{
    int k, numInterval, i, j;
    char filename[] = "2Dgauss.txt";
    double gaussSum, localSum, relative_err;
    double h;
    FILE* fp;

    fp = fopen(filename, "w");

    // using different number of Gaussian points(k*k).
    for (k = 2; k <= 4; k++) {
        fprintf(fp, "No. of Gaussian points(k*k) = %d*%d\n", k, k);
        // using different numbers of intervals(numInterval*numInterval)
        for (numInterval = 1; numInterval <= 4; numInterval++) {
            if (numInterval == 3) continue;
            gaussSum = 0.0;
            fprintf(fp, "   No. of interval = %d*%d\n", numInterval, numInterval);
            h = (b - a) / numInterval;
            for (i = 1; i <= numInterval; i++) {
                for (j = 1; j <= numInterval; j++) {
                    localSum = gauss_quadrature(numInterval, k, a + (j - 1) * h, a + j * h, a + (i - 1) * h, a + i * h, my_func);
                    gaussSum += localSum;
                }
            }
            fprintf(fp, "       Integral= %15.15f\n", gaussSum);
        }
    }
    fclose(fp);

}
void Q5B()
{
    int k, numInterval, i, j;
    char filename[] = "2Dgauss&relative_err.txt";
    double gaussSum, localSum, relative_err;
    double h;
    FILE* fp;

    fp = fopen(filename, "w");

    // using different number of Gaussian points(k*k).
    for (k = 2; k <= 4; k++) {
        fprintf(fp, "No. of Gaussian points(k*k) = %d*%d\n", k, k);
        // using different numbers of intervals(numInterval*numInterval)
        for (numInterval = 1; numInterval <= 4; numInterval++) {
            if (numInterval == 3) continue;
            gaussSum = 0.0;
            fprintf(fp, "   No. of interval = %d*%d\n", numInterval, numInterval);
            h = (b - a) / numInterval;
            for (i = 1; i <= numInterval; i++) {
                for (j = 1; j <= numInterval; j++) {
                    localSum = gauss_quadrature(numInterval, k, a + (j - 1) * h, a + j * h, a + (i - 1) * h, a + i * h, my_func);
                    gaussSum += localSum;
                }
            }
            relative_err = relative_error(gaussSum);
            fprintf(fp, "       Integral= %15.15f\n", gaussSum);
            fprintf(fp, "       relative_err= %15.15f\n", relative_err);
        }
    }
    fclose(fp);

}

int main() {
    //Q4();
    //Q5B();
    //Q5C();
    Q5D();
    return 0;
}