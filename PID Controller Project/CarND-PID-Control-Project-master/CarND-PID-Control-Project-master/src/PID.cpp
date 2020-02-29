#include "PID.h"
#include <iostream>
#include <cmath>
#include <limits>

using namespace std;

/*
* TODO: Complete the PID class.
*/

PID::PID() {}

PID::~PID() {}

void PID::Init(double Kp, double Ki, double Kd) {

 // Initialize the tunable parameters
 PID::Kp = Kp ;
 PID::Ki = Ki ;
 PID::Kd = Kd ;
 
 //Initialize the error variables
 p_error = 0.0f;
 d_error = 0.0f;
 i_error = 0.0f;

 // Implementing Twiddle
 twiddle_toggle = false;
 dp = {0.1*Kp,0.01*Ki,0.1*Kd};
 step = 1;
 param_idx = 0;
 settle_steps = 100;
 eval_steps = 1000;
 error_total = 0;
 error_max = std::numeric_limits<double>::max();
 add = false;
 subtract = false;

 
}

void PID::UpdateError(double cte)
{

d_error = cte - p_error; // Since p_error represents previous cte
p_error = cte;
i_error += cte;

// Updating Total Error
if(step %(settle_steps+eval_steps) > settle_steps)
{
 error_total = error_total + cte*cte;
}

//Twiddle
if(twiddle_toggle && step%(eval_steps+settle_steps)==0)
{
 std::cout<<"Step:"<<step<<endl;
 std::cout<<"Total Error:"<<error_total<<endl;
 std::cout<<"Best Error:"<<error_max<<endl;

 if(error_total<error_max)
 {
  error_max = error_total;
  if(step != settle_steps + eval_steps)
  {
   dp[param_idx]=dp[param_idx]*1.1;
  }
  //Moving on to the next parameter
  param_idx=(param_idx+1)%3;
  add = false;
  subtract = false;
 }

 if(!add && !subtract)
 {
  TwiddleMod(param_idx,dp[param_idx]);
  add = true;
 }
 else  if(add && !subtract)
 {
  TwiddleMod(param_idx,-2*dp[param_idx]);
  subtract = true;
 }
 else
 {
  TwiddleMod(param_idx,dp[param_idx]);
  dp[param_idx]=dp[param_idx]*0.9;
  // Next parameter
  param_idx = (param_idx+1)%3;
  add = false;
  subtract = false;
 }
 error_total = 0;
 std::cout<<"New parameters"<<endl;
 std::cout<<"P :"<<Kp<<", I :"<<Ki<<", D :"<<Kd<<endl;
 std::cout<<"Param Index: "<<param_idx<<endl;
}
step++;

}

double PID::TotalError() {

double output;
output = -Kp*p_error - Ki*i_error - Kd*d_error;

return output;

}

void PID::TwiddleMod(int index,double amount)
{
 if(index == 0)
 {
  Kp = Kp + amount;
 }
 else if(index == 1)
 {
  Ki = Ki + amount;
 }
 else if(index == 2)
 {
  Kd = Kd + amount;
 }
 else
 {
  std::cout<<"Index out of bounds"<<endl;
 }

}
