#include <ceres_mpc/mpc_mm_ceres.h>

MPC_cost::MPC_cost( MatrixXd A, MatrixXd B, MatrixXd Bd, MatrixXd Q,
                       MatrixXd Q_final, MatrixXd R, MatrixXd R_delta,
                        MatrixXd disturbance , int num_params, int horizon) : num_params_(num_params),horizon(horizon) ,A_(A), B_(B), Bd_(Bd), Q_(Q), Q_final_(Q_final), R_(R),R_delta_(R_delta), disturbance_(disturbance)
            {   
                this->insecure_ = this->Bd_ * disturbance;
                x_states = MatrixXd::Zero(A_.rows(), horizon+1);
                deriv_wrt_u =  MatrixXd::Zero(B_.cols(), horizon);
                u = MatrixXd::Zero(B_.cols(), 1);
                u_past = MatrixXd::Zero(B_.cols(), 1);
                u_current = MatrixXd::Zero(B_.cols(), 1);
                lambdas_x = MatrixXd::Zero(A_.rows(), horizon);
                lambdas_u = MatrixXd::Zero(B_.cols(), horizon);
                lambdas_u_ref = MatrixXd::Zero(B_.cols(), horizon);
                u_horizon = MatrixXd::Zero(B_.cols(), horizon);
                lambdas_u_ref = MatrixXd::Zero(B_.cols(), horizon);
                x0_ =  MatrixXd::Zero(x_states.rows(), 1);
                x_ss_ =  MatrixXd::Zero(x_states.rows(), 1);
            }

//When you add the const keyword to a method the this pointer will essentially become a pointer to const object, and you cannot therefore change any member data. (This is not completely true, because you can mark a member as mutable and a const method can then change it. It's mostly used for internal counters and stuff.).

bool MPC_cost::Evaluate(double const* const* x,
                                   double* residuals,
                                   double** jacobians) const {
 



               
        x_states.block(0,0,x0_.rows(), x0_.cols()) = x0_;
        
        u_past.block(0,0,u_current.rows(), u_current.cols())= u_current;
       //bilo je horizon -1
        for (int i = 0; i< horizon; i++){
		
		for(int j = 0; j < num_params_; j++){
		    u(j,0) = x[0][j*horizon + (i)],           
		       
		}
            x_states.block(0,i+1,x0_.rows(), x0_.cols()) = A_ * x_states.block(0,i,x0_.rows(), x0_.cols()) + B_*u + insecure_;
            lambdas_x.block(0,i,x0_.rows(), x0_.cols())  = -1*x_ss_ + x_states.block(0,i,x0_.rows(), x0_.cols());
            lambdas_u.block(0,i,u.rows(), u.cols()) = u - u_past;
            lambdas_u_ref.block(0,i,u.rows(), u.cols()) = u - u_current;
            deriv_wrt_u.block(0,i,u.rows(), u.cols()) =(2*u.transpose()*R_).transpose() + (4*u.transpose()*R_delta_).transpose()
                                         + (-1*R_delta_*u_past) + (-1*u_past.transpose()*R_delta_).transpose();
            if(i > 0){ deriv_wrt_u.block(0,i-1,u.rows(), u.cols()) =(deriv_wrt_u.block(0,i-1,u.rows(), u.cols()) + (-1*u.transpose()*R_delta_).transpose() -1*R_delta_*u).eval();}
            
            u_past.block(0,0,u.rows(), u.cols()) = u;
        }
        
        residuum = ((lambdas_x.cwiseProduct(Q_*lambdas_x)).sum() + (lambdas_u_ref.cwiseProduct(R_*lambdas_u_ref)).sum()                 + (lambdas_u.cwiseProduct(R_delta_*lambdas_u)).sum()     +((-x_ss_ + x_states.block(0,horizon,x0_.rows(),x0_.cols())).transpose()*Q_final_*(-x_ss_ + x_states.block(0,horizon,x0_.rows(),x0_.cols()))).sum() );

        
        

        residuals[0] = residuum;

        if (jacobians != NULL) {
                if (jacobians[0] != NULL) {
                    for(int i = 0; i<num_params_; i++){
                        for(int j = 0; j< horizon; j++){
                            jacobians[0][i*horizon + j] = deriv_wrt_u(i,j);
                            }
                        } 
                }
            
        }

        return true;
}

 
