#include <ceres_mpc/mpc_mm_ceres.h>

using namespace Eigen;
using namespace std;
using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;
using ceres::CostFunction;
using ceres::SizedCostFunction;

bool SolveMyOptimizationProblem(ceres::Problem& problem) {
  //CHECK(problem != NULL);

  // Run the solver!
  ceres::Solver::Options options;
  options.max_num_iterations = 1000;
  options.linear_solver_type = ceres::DENSE_QR;
  options.minimizer_progress_to_stdout = false;
  ceres::Solver::Summary summary;
  ceres::Solve(options, &problem, &summary);

  //std::cout << summary.FullReport() << '\n';
  std::cout << summary.BriefReport() << '\n';
std::cout << summary.IsSolutionUsable() << '\n';
  return summary.IsSolutionUsable();
}

int main(int argc, char** argv) {
  google::InitGoogleLogging(argv[0]);



MatrixXd model_A_ = MatrixXd::Random(8,8);
MatrixXd model_B_ = MatrixXd::Random(model_A_.rows(),4);
MatrixXd model_Bd_ = MatrixXd::Random(model_A_.rows(),6);
MatrixXd Q = MatrixXd::Random(model_A_.rows(), model_A_.rows());
MatrixXd Q_final = MatrixXd::Random(model_A_.rows(),model_A_.rows());
MatrixXd R = MatrixXd::Random( model_B_.cols(), model_B_.cols());
MatrixXd R_delta = MatrixXd::Random(model_B_.cols(),model_B_.cols());
MatrixXd estimated_disturbances_ = MatrixXd::Random(model_Bd_.cols(),1);


  double x[] = {1.0,1.0, 1.0,1.0,1.0,1.0, 1.0,1.0,1.0,1.0, 1.0,1.0,1.0,1.0, 1.0,1.0,1.0,1.0, 1.0,1.0,1.0,1.0, 1.0,1.0,
  1.0,1.0, 1.0,1.0,1.0,1.0, 1.0,1.0,1.0,1.0, 1.0,1.0,1.0,1.0, 1.0,1.0,1.0,1.0, 1.0,1.0,1.0,1.0, 1.0,1.0};

 
  MPC_cost*  cost1 = new MPC_cost(  model_A_,  model_B_,  model_Bd_,  Q, Q_final,
                              R,  R_delta, estimated_disturbances_, 4, 12);

  ceres::Problem problem;



  problem.AddResidualBlock(cost1, NULL, x);

for(int i = 0; i< 48; i++){
    problem.SetParameterLowerBound(x,i, -10);
  problem.SetParameterUpperBound(x,i, -5);
  }
for(int i = 0; i< 4; i++){

/*
problem.SetParameterLowerBound(x,1, lower_bounds(1,0)+i*10);
  problem.SetParameterUpperBound(x,1, upper_bounds(1,0)+i*10);
*/

if(!SolveMyOptimizationProblem(problem))
    cout << "The solve was not successful, exiting." << endl;
cout << x[0] <<"  " << x[1] << endl;
cout << "zavrsio sam " << endl;
}
  // ceres::Solver::Options options;
  // options.max_num_iterations = 100;
  // options.linear_solver_type = ceres::DENSE_QR;
  // options.minimizer_progress_to_stdout = true;
  // ceres::Solver::Summary summary;
  // ceres::Solve(options, &problem, &summary);

  // std::cout << summary.FullReport() << '\n';
  
 return 0;
}
