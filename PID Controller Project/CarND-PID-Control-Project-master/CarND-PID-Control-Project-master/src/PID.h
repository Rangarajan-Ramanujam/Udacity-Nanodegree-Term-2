#ifndef PID_H
#define PID_H
#include <vector>

class PID {
public:
  /*
  * Errors
  */
  double p_error;
  double i_error;
  double d_error;

  /*
  * Coefficients
  */ 
  double Kp;
  double Ki;
  double Kd;
  
  //Twiddle parameters
  std::vector<double> dp;
  int step, param_idx;
  int settle_steps, eval_steps;
  double error_total, error_max;
  bool twiddle_toggle,add,subtract;
 
  /*
  * Constructor
  */
  PID();

  /*
  * Destructor.
  */
  virtual ~PID();

  /*
  * Initialize PID.
  */
  void Init(double Kp, double Ki, double Kd);

  /*
  * Update the PID error variables given cross track error.
  */
  void UpdateError(double cte);

  /*
  * Calculate the total PID error.
  */
  double TotalError();

  // Adaptive Gain
  void TwiddleMod(int index , double amount);
};

#endif /* PID_H */
