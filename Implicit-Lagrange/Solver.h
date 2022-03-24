#ifndef ELEMENTS_SOLVER_H
#define ELEMENTS_SOLVER_H  

class Solver{

public:
  Solver();
  ~Solver();
  
  virtual void setup() {}

  virtual void run(int argc, char *argv[]) = 0;

  virtual void finalize() {}

  virtual void exit_solver(int status);

  int setup_flag, finalize_flag;

};

#endif // end Header Guard
