void ITSMC_REF_force_control(void)
{

    int controller = 1;
    int axes = 0;

    double time_sliding = 0;

    double pos[] = {0, 0, 0, 0, 0, 0};
    double vel[] = {0, 0, 0, 0, 0, 0};
    double torques[] = {0, 0, 0, 0, 0, 0};
    double tau_f[] = {0, 0, 0, 0, 0, 0};
    double tor[] = {0, 0, 0, 0, 0, 0};
    double c_robot[] = {0, 0, 0, 0, 0, 0};
    double g_robot[] = {0, 0, 0, 0, 0, 0};
    double tau_f_c_robot[] = {0, 0, 0, 0, 0, 0};
    double tau_f_v_robot[] = {0, 0, 0, 0, 0, 0};
    double tau_f_b_robot[] = {0, 0, 0, 0, 0, 0};
    double d_robot[] = {0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0};

    double F[] = {0, 0, 0, 0, 0, 0};
    double Fd[] = {0, 0, 15, 0, 0, 0};
    double JF[] = {0, 0, 0, 0, 0, 0};
    double J[] = {0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0};
    double qd[] = {0, 0, 0, 0, 0, 0};          
    double dqd[] = {0, 0, 0, 0, 0, 0};  
    double ddqd[] = {0, 0, 0, 0, 0, 0};  
    double Ddqdd[] = {0, 0, 0, 0, 0, 0};  

    double Xd=0;
    double dXd=0;
    double ddXd=0;
    double Yd=0;
    double dYd=0;
    double ddYd=0;
    double Zd=0;
    double dZd=0;
    double ddZd=0;
    double Zp=0;
    double dZp=0;
    double ddZp=0;

    // Eigen::MatrixXd d_matri(6,6);
    Eigen::MatrixXd inv_d(6,6);
    Eigen::MatrixXd tran_d(6,6);
    // Eigen::MatrixXd J_matri(6,6);
    Eigen::MatrixXd tran_J(6,6);

    Eigen::MatrixXd sliding_matri(6,1);
    Eigen::MatrixXd sliding_mode_matri(6,1);
    Eigen::VectorXd tau_f_0(6);

    //滑模的积分结果
    Eigen::MatrixXd integral_sliding_mode(6,1);

    //滑模最终结果
    Eigen::MatrixXd result_integral_sliding_mode(6,1);

    for(int k = 0; k < 6; k++)
    {
        integral_sliding_mode(k) = 0;
    }

    // FTC+SM
    double k1[] = {30, 50, 20, 20, 20, 5};
    double k2[] = {25, 30, 20, 20, 15, 4};
    double k3[] = {3, 3, 2, 2, 1, 1};

    // PD+SM
    double pd_k1[] = {80, 150, 40, 80, 50, 40};
    double pd_k2[] = {85, 50, 20, 20, 10, 10};

    double alpha_1 = 0, alpha_2 = 0;

    double e1[] = {0, 0, 0, 0, 0, 0};
    double e2[] = {0, 0, 0, 0, 0, 0};

    double deno = 0;

    int d_mark = 0;
    int j_mark = 0;

    double k_1 = 0, k_2 = 0, k_3 = 0;

    sleep(5);
    power_on();

    alpha_1 = 0.35;
    alpha_2 = 2 * alpha_1 / (1 + alpha_1);
    time_sliding = 0;

    struct timespec st;
    double t0 = 0.0;
    double t1 = 0.0;
    double dt = 0.0;
    int count = 0;

    // 在循环外初始化起始时间
    clock_gettime(CLOCK_REALTIME, &st);
    double t_start = st.tv_sec + (st.tv_nsec / 1e9);

    while (true)
    {

        clock_gettime(CLOCK_REALTIME, &st);
        t0 = st.tv_sec + (st.tv_nsec / 1e9);
        double t = t0 - t_start;

        for (int k = 0; k < 6; k++)
        {
            pos[k] = get_pos(k) * PI / 180;
        }

        double* jacobian = yakeb(pos);
        for (int i = 0; i < 36; ++i) {
            J[i] = jacobian[i];
        }
        Xd=0.35+0.05*cos((pi/5)*t);
        Yd=0+0.05*sin((pi/5)*t);
        Zp=0.4;

        dXd=-0.01*pi*sin((pi/5)*t);
        dYd=0.01*pi*cos((pi/5)*t);
        ddXd=-0.002*pi*pi*cos((pi/5)*t);
        ddYd=-0.002*pi*pi*sin((pi/5)*t);

        //力传感器读取力 
        double* sensor_values1 = read_sensor();
        for (int k = 0; k < 6; ++k) {
            F[k] = sensor_values1[k];
        } 
        impedance_control(Zp,dZp,ddZp,F,Fd,Zd,dZd,ddZd);   
        inv_kin(Xd,dXd,ddXd,Yd,dYd,ddYd,Zd,dZd,ddZd,J,qd,dqd,ddqd);

        for (int k = 0; k < 6; k++)
        {
            if (3 == k or 1 == k or 5 == k)
            {
                vel[k] = - get_vel(k);
            }
            else
            {
                vel[k] = get_vel(k);
            }
        }

        dyn_mat_ER3A(pos, vel, d_robot, c_robot, g_robot, tau_f_v_robot, tau_f_c_robot,tau_f_b_robot);

        Eigen::Map<Eigen::MatrixXd> d_matri(d_robot, 6, 6);  // 直接映射double数组为MatrixXd

        //计算误差
        for(int k = 0; k < 6; k++)
        {
            e1[k] = pos[k] - qd[k];
            e2[k] = vel[k] - dqd[k];
        }

        for (int k = 0; k < 6; k++)
        {
            torques[k] =  k1[k] * sig(e1[k], alpha_1) + k2[k] * sig(e2[k], alpha_2);
        }

        inv_d = d_matri.inverse();
        tran_d = inv_d.transpose();

        for(int k = 0; k < 6; k++)
        {
            sliding_matri(k) = torques[k] + c_robot[k]*e2[k] + tau_f_v_robot[k]*e2[k];
        }

        //滑模的积分部分 得到被积函数矩阵
        sliding_mode_matri = inv_d * sliding_matri;
        integral_sliding_mode =  integral_sliding_mode  + sliding_mode_matri * dt;

        //滑模最终结果
        for (int k = 0; k < 6; k++)
        {
            result_integral_sliding_mode(k) =  integral_sliding_mode(k)  + e2[k];
        }

        //D^{-T}S
        result_integral_sliding_mode = tran_d * result_integral_sliding_mode;

        Eigen::Map<Eigen::MatrixXd> J_matri(J, 6, 6);  // 直接映射double数组为MatrixXd
        tran_J=J_matri.transpose();

        double* sensor_values2 = read_sensor();
        for (int k = 0; k < 6; ++k) {
            F[k] = sensor_values2[k];
        } 

        //有限时间控制器
        Ddqdd=d_matri*ddqd;
        JF=tran_J*F;
        torques[0] =  - torques[0] - k3[0]*sat(result_integral_sliding_mode(0),0.8) + g_robot[0] + Ddqdd[0] + tau_f_c_robot[0]*sgn(e2[0]) + tau_f_b_robot[0] - JF[0];
        torques[1] =  - torques[1] - k3[1]*sat(result_integral_sliding_mode(1),0.8)  + g_robot[1] + Ddqdd[1] + tau_f_c_robot[1]*sgn(e2[1]) + tau_f_b_robot[1] - JF[1];
        torques[2] =  - torques[2] - k3[2]*sat(result_integral_sliding_mode(2),0.8)  + g_robot[2] + Ddqdd[2] + tau_f_c_robot[2]*sgn(e2[2]) + tau_f_b_robot[2] - JF[2];
        torques[3] =  - torques[3] - k3[3]*sat(result_integral_sliding_mode(3),0.8) + g_robot[3] + Ddqdd[3] + tau_f_c_robot[3]*sgn(e2[3]) + tau_f_b_robot[3] - JF[3];
        torques[4] =  - torques[4] - k3[4]*sat(result_integral_sliding_mode(4),0.8)  + g_robot[4] + Ddqdd[4] + tau_f_c_robot[4]*sgn(e2[4]) + tau_f_b_robot[4] - JF[4];
        torques[5] =  - torques[5] - k3[5]*sat(result_integral_sliding_mode(5),0.8)  + g_robot[5] + Ddqdd[5] + tau_f_c_robot[5]*sgn(e2[5]) + tau_f_b_robot[5] - JF[5];

        Eigen::Map<const Eigen::VectorXd> pos_eigen(pos, 6);
        Eigen::Map<const Eigen::VectorXd> vel_eigen(vel, 6);
        Eigen::Map<const Eigen::VectorXd> qd_eigen(qd, 6);
        Eigen::Map<const Eigen::VectorXd> dqd_eigen(dqd, 6);
        Eigen::Map<const Eigen::VectorXd> ddqd_eigen(ddqd, 6);

        tau_f_0 = rbf_tau(pos_eigen,vel_eigen,qd_eigen,dqd_eigen,ddqd_eigen);
        for (int i = 0; i < 6; ++i) {
            tau_f[i] = tau_f_0(i);
        }
        torques[0] = torques[0] - tau_f[0];
        torques[1] = torques[1] - tau_f[1];
        torques[2] = torques[2] - tau_f[2];
        torques[3] = torques[3] - tau_f[3];
        torques[4] = torques[4] - tau_f[4];
        torques[5] = torques[5] - tau_f[5];
        
        for (int k = 0; k < 6; k++) 
        {
            
            if (3 == k or 1 == k or 5 == k)
            {
                torques[k] = -torques[k];
            }
            else
            {
                torques[k] = torques[k];
            }
            set_torque(k, torques[k]);
        }

        clock_gettime(CLOCK_REALTIME, &st);
        t1 = st.tv_sec + (st.tv_nsec / 1e9);
        dt = t1 - t0;

        time_sliding = dt;
        count++;
    }
}

