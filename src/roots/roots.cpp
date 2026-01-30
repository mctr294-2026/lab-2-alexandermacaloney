#include "roots.hpp"
#include <cmath>      // for std::abs
#include <algorithm>  // for std::swap

bool bisection(std::function<double(double)> f, 
// Returns true if the root has been found. 
// Returns false otherwise (i.e.: if unable to finmd a root or if there are invalid inputs)
// std::function< > holds the callable function
// double(double) takes an input of type double (a type of precise float number) and spits out a double
               double a, double b, // This line represents the two endpoints of the interval [a,b]
               double* root) // Pointer to root. Answer stored here. Allows function to write into variable that exists outside of the function
{
    if (root == nullptr) {     //Check pointer. Null or invalid pointer returns 'false'
        return false;
    }

    if (a>b){ //Check if a is larger than b
        std::swap(a, b); // If so, then interchange the two numbers, so a is the smaller value
    }

    double fa = f(a); //Evalates f(a) and assigns the result to the variable fa
    double fb = f(b); //Same but for f(b)

    if (!std::isfinite(fa) || !std::isfinite(fb)) {  // Is f(a) not a finite number? Or f(b)? If so, return false
        return false;
    }

    const double tolerance = 1e-6; // Set the required tolerance to 1x10^-6. Creates a const variable of type double. Const variables can't be accidentally changed.
    const int max_iteration = 1000000;  // Set the max iteration as a constant varriable
    
    if (std::abs(fa) <= tolerance){ // Check if the absolute value of a is less than the toleranmce value
        *root = a; // If so, store the value of a in the root pointer
        return true; // return true to the function since the root has been found
    }
    if (std::abs(fb) <= tolerance){ // Check if the absolute value of b is less than the toleranmce value
        *root = b; // If so, store the value of b in the root pointer
        return true; // return true to the function since the root has been found
    }

    if (fa*fb > 0.0) { // Check if the value is positive. The value must be negative since a positive # x a negative # is negative
        return false;
    }

    for (int iteration = 0; iteration < max_iteration; ++iteration) { // Starts at interation 0 and increments byh 1 upwards while its less than the max iteration value
        const double mid = 0.5 * (a + b);   // Find midpoint between a and b and assign it toa new constant variable
        const double fmid = f(mid); // evaluate the function at the midpoint and assign it to a new constant variable

        if (!std::isfinite(fmid)) { // Is f(mid) not a finite number? If so, return false
            return false;
        }

        // Pass
        if (std::abs(fmid) <= tolerance || 0.5 * (b - a) <= tolerance) { // If the absolute value of the midpoint is less than the tolerance OR if the width of the interval itself is less than the tolerance,
            *root = mid; // then the root is the midpoint. Assign it the root pointer
            return true;
        }

        // Shrink interval: keep the half that still has the sign change
        if (fa * fmid < 0.0) { // root in interval [a, m]. a is the negative value, so fmid must be positive for this to be true
            b = mid; // replace b with the midpoint value
            fb = fmid; // replace fb with fmid
        } else { // root in interval [m, b]. b is the positive number. so fmid must be negative
            a = mid;// replace a with the midpoint value
            fa = fmid; // replace fa with fmid
        }
    }

    // If max_iteratioin is reached without converging
    return false; // End the function and return the bool value False
}

bool regula_falsi(std::function<double(double)> f, 
// Returns true if the root has been found. 
// Returns false otherwise (i.e.: if unable to finmd a root or if there are invalid inputs)
// std::function< > holds the callable function.
// double(double) takes an input of type double (a type of precise float number) and spits out a double
               double a, double b, // This line represents the two endpoints of the interval [a,b]
               double* root) // Pointer to root. Answer stored here. Allows function to write into variable that exists outside of the function
{
    if (root == nullptr) {     //Check pointer. Null or invalid pointer returns 'false'
        return false;
    }

    if (a>b){ //Check if a is larger than b
        std::swap(a, b); // If so, then interchange the two numbers, so a is the smaller value
    }

    double fa = f(a); //Evalates f(a) and assigns the result to the variable fa
    double fb = f(b); //Same but for f(b)

    if (!std::isfinite(fa) || !std::isfinite(fb)) {  // Is f(a) not a finite number? Or f(b)? If so, return false
        return false;
    }

    const double tolerance = 1e-6; // Set the required tolerance to 1x10^-6. Creates a const variable of type double. Const variables can't be accidentally changed.
    const int max_iteration = 1000000;  // Set the max iteration as a constant varriable
    
    if (std::abs(fa) <= tolerance){ // Check if the absolute value of a is less than the toleranmce value
        *root = a; // If so, store the value of a in the root pointer
        return true; // return true to the function since the root has been found
    }
    if (std::abs(fb) <= tolerance){ // Check if the absolute value of b is less than the toleranmce value
        *root = b; // If so, store the value of b in the root pointer
        return true; // return true to the function since the root has been found
    }

    if (fa*fb > 0.0) { // Check if the value is positive. The value must be negative since a positive # x a negative # is negative
        return false;
    }

    for (int iteration = 0; iteration < max_iteration; ++iteration) { // Starts at interation 0 and increments byh 1 upwards while its less than the max iteration value
        const double denominator = (fb-fa);   // Find the value of the denominator and store it in a constant variable
        if (denominator == 0.0){
            return false; //Not allowed to divide by zero
        }

        const double c = (a*fb-b*fa) / denominator; // Formula for midpoint calculation in regula falsi
        const double fc = f(c); // //Evalates f(c) and assigns the result to the variable fc

        if (!std::isfinite(fc)) { // Is f(c) not a finite number? If so, return false
            return false;
        }

        // Pass
        if (std::abs(fc) <= tolerance) { // If the absolute value of the midpoint is less than the tolerance,
            *root = c; // then the root is the midpoint. Assign it the root pointer
            return true;
        }

        // Shrink interval: keep the half that still has the sign change
        if (fa * fc < 0.0) { // root in interval [a, c]. a is the negative value, so fc must be positive for this to be true
            b = c; // replace b with the midpoint value
            fb = fc; // replace fb with fc
        } else { // root in interval [c, b]. b is the positive number. so fc must be negative
            a = c;// replace a with the midpoint value
            fa = fc; // replace fa with fc
        }
    }

    // If max_iteratioin is reached without converging
    return false; // End the function and return the bool value False
}

bool newton_raphson(std::function<double(double)> f,  // Two callable functions needed for newton_raphson. This one is f(x)
                    std::function<double(double)> g, //This one is g(x)
                    double a, double b, double c, // Third double is the intial guess x
                    double* root) // Pointer to root. Answer stored here. Allows function to write into variable that exists outside of the function
{
    if (root == nullptr) {     //Check pointer. Null or invalid pointer returns 'false'
        return false;
    }

    if (a>b){ //Check if a is larger than b
        std::swap(a, b); // If so, then interchange the two numbers, so a is the smaller value
    }

    const double tolerance = 1e-6; // Set the required tolerance to 1x10^-6. Creates a const variable of type double. Const variables can't be accidentally changed.
    const int max_iteration = 1000000;  // Set the max iteration as a constant varriable
    
    double x = c; //Initial guess assigned to variable x

    if (x < a || x > b){ //Check to ensure that the initial guess x is between the a and b. False otherwwise
        return false;
    }

    for (int iteration = 0; iteration < max_iteration; ++iteration) { // Starts at interation 0 and increments byh 1 upwards while its less than the max iteration value
        const double fx = f(x);   //Evalates f(x) and assigns the result to the variable fx
        const double gx = g(x); //Same but for g(x)

        if (!std::isfinite(fx) || !std::isfinite(gx)) { // Is f(x) OR g(x) not a finite number? If so, return false
            return false;
        }

        // Pass
        if (std::abs(fx) <= tolerance) { // If the absolute value of f(x) is less than the tolerance
            *root = x; // then the root is the initial guess x. Assign it the root pointer
            return true;
        }
        
        if (gx == 0.0){ // Divide by zero is not allowed
            return false;
        }

        const double x_new = x - (fx / gx);  // generate a new guess for x by taking the orignal guess and subtraccting the quotient of the two functions

        if (x_new < a || x_new > b){ //Check to ensure that the new guess x is between the a and b. False otherwwise
            return false;
        }

        x = x_new; //Assign the new guess as x before looping
    }

    // If max_iteratioin is reached without converging
    return false; // End the function and return the bool value False
}

bool secant(std::function<double(double)> f, // One callable function f(x)
            double a, double b, double c, // range and the intial guess x
            double* root) // Pointer to root. Answer stored here. Allows function to write into variable that exists outside of the function
{
    if (root == nullptr) {     //Check pointer. Null or invalid pointer returns 'false'
        return false;
    }

    if (a>b){ //Check if a is larger than b
        std::swap(a, b); // If so, then interchange the two numbers, so a is the smaller value
    }

    const double tolerance = 1e-6; // Set the required tolerance to 1x10^-6. Creates a const variable of type double. Const variables can't be accidentally changed.
    const int max_iteration = 1000000;  // Set the max iteration as a constant varriable
    
    double x0 = c; //Initial guess assigned to variable x0
    double x1 = b; //One of the endpoints

    if (x0 < a || x0 > b){ //Check to ensure that the initial guess x0 is between a and b. False otherwwise
        return false;
    }

    double f0 = f(x0); //Evalates f(x0) and assigns the result to the variable f0
    double f1 = f(x1); // Same but for f(x1)

    if (!std::isfinite(f0) || !std::isfinite(f1)) { // Is f(0) OR f(1) not a finite number? If so, return false
        return false;
    }

    if (std::abs(f0) <= tolerance){ // Check if the endpoints or initial guess happened to be a root
        *root = x0; // If so, store the value of x0 in the root pointer
        return true;
    }
    
    if (std::abs(f1) <= tolerance){ // Check if the endpoints or initial guess happened to be a root
        *root = x1; // If so, store the value of x1 in the root pointer
        return true;
    }

    for (int iteration = 0; iteration < max_iteration; ++iteration) { // Starts at interation 0 and increments byh 1 upwards while its less than the max iteration value
        const double denominator = (f1 - f0);   //Evalates value of denominator
        if (denominator == 0.0){
            return false; // Can't divide by zero
        }

        const double x2 = x1 - (f1 * (x1-x0)/denominator); // Generates a new guess to update the secant line based on recursive logic

        if (x2 < a || x2 > b){ //Check to ensure that the generated guess x2 is between a and b. False otherwwise
            return false;
        }

        const double f2 = f(x2); //Evalates f(x2) and assigns the result to the variable f2
        if (!std::isfinite(f2)){ // Is f(2) not a finite number? If so, return false
            return false;
        }

        // Pass
        if (std::abs(f2) <= tolerance) { // If the absolute value of f(2) is less than the tolerance
            *root = x2; // then the root is the  guess x2. Assign it the root pointer
            return true;
        }
        
         //Adjust the number in prep for next loop
        x0 = x1;
        f0 = f1;
        x1 = x2;
        f1 = f2;
        
        
    }

    // If max_iteratioin is reached without converging
    return false; // End the function and return the bool value False
}

// Some of the above functions were partially written by AI