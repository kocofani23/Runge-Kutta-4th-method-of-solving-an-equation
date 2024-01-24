#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MAX 10

//struct to create and initialize values of equation
struct diff_equation
{
    int variable_no;
    int x_count;
    int y_count;

    double x0;
    double y0;
    double x_target;
    double h;
    double exact_value;

    double **x_coeffs;
    double y_coeff;
    double y_dev_coeff;

} equation;

//function to create equation based on the inputs of the user
void create_equation(int count, struct diff_equation *equation)
{
    int i;

    printf("Function looks like this: ");
    printf("[y_dev_coeff * y' + x_coef * x ^ x_power + y_coef * y = 0]\n");
    printf("Total number of variables 'x' and 'y' should be : %d\n", equation->variable_no);
    printf("Total number of x variables should be %d\n", equation->variable_no - 1);
    printf("Enter coefficient of first order derivative: ");
    scanf("%lf", &equation->y_dev_coeff);

    equation->x_count = equation->variable_no - 1;
    printf("\nEnter number of x variables: ");
    scanf("%d", &equation->x_count);
    equation->x_coeffs = (double**)calloc(2, sizeof(double*));
    for(int i=0; i<2; i++){
        equation->x_coeffs[i] = (double *)calloc(equation->x_count, sizeof(double));            //allocate memory dinamically
    }

    printf("Enter coefficients of the 'x' variables and their respective powers: ");
    for(i=0; i<equation->x_count; i++){
        scanf("%lf %lf", &equation->x_coeffs[0][i], &equation->x_coeffs[1][i]);         //store coefficients and powers of x values
    }                                                                                   //coefficients are at row no.0 and powers are at row no.1

    equation->y_count = 1;

    printf("\nEnter the coefficient of the 'y' variable: ");
    scanf("%lf", &equation->y_coeff);


}

//function to evaluate equation based on the inputs and numbers used in runge kutta function
double return_equation(double x, double y, struct diff_equation equation)
{
    double value = 0.0;
    int i;
    int x_count = equation.x_count;
    int y_count = equation.y_count;

    for(i=0; i<x_count; i++){
        value += equation.x_coeffs[0][i] * pow(x, equation.x_coeffs[1][i]);
    }
    value += y * equation.y_coeff;

    return value / equation.y_dev_coeff;                        //divide by the coefficient of the derivative if it is != 1

}

//function to evaluate the approximate value of the ODE based on the runge kutta recursive method
void runge_kutta_function(struct diff_equation equation)
{
    double c1, c2, c3, c4;
    double x = equation.x0;
    double h = equation.h;
    double y = equation.y0;
    double exact_value = equation.exact_value;
    double target = equation.x_target;
    int count = equation.variable_no;

    while(x < target){                      //while x is smaller than target iterate

        printf("x0 = %lf\tx = %lf\ty = %lf\tExact value = %lf\tepsilon = %lf\n", equation.x0, x, y, exact_value, (exact_value-y));

        c1 = return_equation(x, y, equation);
        c2 = return_equation((x + h/2), (y + c1*h/2), equation);
        c3 = return_equation((x + h/2), (y + c2*h/2), equation);
        c4 = return_equation((x + h), (y + c3*h), equation);

        y = y + ((c1 + 2*c2 + 2*c3 + c4) * h / 6);                  //finding y value
        x += h;                                                     //increasing x value by h
    }

}

int main()
{

    struct diff_equation equation;
    int i;
    printf("****************\tRUNGE-KUTTA METHOD\t****************\n");

    printf("Enter number of variables of equation(max number of variables accepted is 10): ");
    scanf("%d", &equation.variable_no);
    if(equation.variable_no > 10){
        printf("Number of variables in equation should be less than 10...");
        return 1;
    }

    create_equation(equation.variable_no, &equation);

    printf("Enter initial value of x: ");
    scanf("%lf", &equation.x0);

    printf("Enter initial value of y [y(x0) = y0]: ");
    scanf("%lf", &equation.y0);

    printf("Enter target value of x: ");
    scanf("%lf", &equation.x_target);

    printf("Enter step size: ");
    scanf("%lf", &equation.h);

    printf("Enter exact value of the input equation: ");
    scanf("%lf", &equation.exact_value);

    runge_kutta_function(equation);
    for(i=0; i<equation.x_count; i++){
        free(equation.x_coeffs[i]);
    }
    free(equation.x_coeffs);


    return 0;
}
