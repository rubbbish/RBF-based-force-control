VectorXd rbf_tau(const VectorXd& q, const VectorXd& dq, const VectorXd& qd, const VectorXd& dqd, const VectorXd& ddqd)
{

    static bool initialized = false;
    double s_int[6] = {0};
    double smc2[6] = {0};
    double s[6] = {0};
    static double dt = 0.001;
    double D[36] = {0};
    double J[36] = {0};
    double C[36] = {0};
    double G[6] = {0};
    double x1[6] = {0};
    double x2[6] = {0};
    double f[6] = {0};
    double fc[6] = {0};
    double fv[6] = {0};
    double fb[6] = {0};
    double tau[6] = {0};
    double taud[6] = {0};
    double ddq[6] = {0};
    double ddq_0[6] = {0};
    int j_mark = 0;
    struct timespec st;
    double t0 = 0.0;
    double t1 = 0.0;
    Eigen::MatrixXd J_matri(6,6);
    Eigen::MatrixXd tran_J(6,6);

    static VectorXd w1(5), w1_1(5), w1_2(5);
    static VectorXd w2(5), w2_1(5), w2_2(5);
    static VectorXd w3(5), w3_1(5), w3_2(5);
    static VectorXd w4(5), w4_1(5), w4_2(5);
    static VectorXd w5(5), w5_1(5), w5_2(5);
    static VectorXd w6(5), w6_1(5), w6_2(5);
    static Eigen::VectorXd prev_tau_f;    

    double q_prev[6] = {0};
    double q_pprev[6] = {0};

    if (!initialized) {
        w1.setZero(); w1_1.setZero(); w1_2.setZero();
        w2.setZero(); w2_1.setZero(); w2_2.setZero();
        w3.setZero(); w3_1.setZero(); w3_2.setZero();
        w4.setZero(); w4_1.setZero(); w4_2.setZero();
        w5.setZero(); w5_1.setZero(); w5_2.setZero();
        w6.setZero(); w6_1.setZero(); w6_2.setZero();
        prev_tau_f = Eigen::VectorXd::Zero(6); 
        initialized = true;
    }

    clock_gettime(CLOCK_REALTIME, &st);
    t0 = st.tv_sec + (st.tv_nsec / 1e9);
    dyn_mat_ER3B(q.data(), dq.data(), D, C, G, fv, fc, fb);

    for (int i = 0; i < 6; ++i) {
        x1[i] = q(i) - qd(i);
        x2[i] = dq(i) - dqd(i);
        f[i] = fc[i] * sign(dq(i)) + fb[i] + fv[i] * dq(i);
    }


    double e1[6] = {0};
    double e2[6] = {0};
    double alpha1 = 0.35;
    double alpha2 = 2 * alpha1 / (1 + alpha1);
    for (int i = 0; i < 6; ++i) {
        e1[i] = sig(x1[i], alpha1);
        e2[i] = sig(x2[i], alpha2);
    }

    Eigen::Map<MatrixXd> D_matri(D, 6, 6);
    Eigen::Map<VectorXd> e1_eigen(e1, 6);
    Eigen::Map<VectorXd> e2_eigen(e2, 6);
    Eigen::Map<VectorXd> x2_eigen(x2, 6);
    

    double k1_vals[] = {30, 50, 20, 20, 20, 5};
    double k2_vals[] = {25, 30, 20, 20, 15, 4};
    double k3_vals[] = {6, 6, 4, 4, 2, 2};
    Eigen::Map<Eigen::VectorXd> k1(k1_vals, 6);
    Eigen::Map<Eigen::VectorXd> k2(k2_vals, 6);
    Eigen::Map<Eigen::VectorXd> k3(k3_vals, 6);
    
    Eigen::Map<MatrixXd> C_matri(C, 6, 6);
    Eigen::Map<VectorXd> fv_eigen(fv, 6);
    VectorXd s_eigen = D_matri.inverse() * (k1.cwiseProduct(e1_eigen) + k2.cwiseProduct(e2_eigen) + C_matri * x2_eigen + fv_eigen.cwiseProduct(x2_eigen));
    
    for (int i = 0; i < 6; ++i) {
        s[i] = s_eigen(i);
        s_int[i] += s[i] * dt;
        smc2[i] = s_int[i] + x2[i];
    }


    double* jacobian = yakeb(q.data());
    memcpy(J, jacobian, sizeof(double) * 36);
    
    j_mark = 0;
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            J_matri(i,j) = J[j_mark++];
        }
    }
    tran_J = J_matri.transpose();

    VectorXd F = read_sensor();


    Eigen::Map<VectorXd> tau_eigen(tau, 6);
    Eigen::Map<VectorXd> G_eigen(G, 6);
    Eigen::Map<VectorXd> fc_eigen(fc, 6);
    Eigen::Map<VectorXd> fb_eigen(fb, 6);
    

    VectorXd sat_eigen(6);
    VectorXd smc2_eigen(6);
    for (int i = 0; i < 6; ++i) {
        smc2_eigen(i) = smc2[i];
    }
    VectorXd D_inv_T_smc2 = D_matri.inverse().transpose() * smc2_eigen;
    for (int i = 0; i < 6; ++i) {
        sat_eigen(i) = sat(D_inv_T_smc2(i), 0.8);
    }
    

    VectorXd sign_eigen(6);
    for (int i = 0; i < 6; ++i) {
        sign_eigen(i) = sign(e2[i]);
    }
    
    tau_eigen = -k1.cwiseProduct(e1_eigen) - k2.cwiseProduct(e2_eigen) - k3.cwiseProduct(sat_eigen) 
                + D_matri * ddqd + G_eigen + fc_eigen.cwiseProduct(sign_eigen) + fb_eigen - prev_tau_f - tran_J * F;

    for (int i = 0; i < 6; ++i) {
        tau[i] = tau_eigen(i);
    }

    ddq_0 = acc_fcn(q.data(), q_prev, q_pprev, dt);
    ddq = filter(ddq_0);


    for (int i = 0; i < 6; ++i) {
        taud[i] = D[i*6 + i] * ddq[i] + C[i*6 + i]* dq[i] + G[i] + f[i] - tau[i];
    }


    double bj[6] = {1.92, 1.53, 1.35, 1.43, 1.32, 0.65};
    double xite = 0.5;  
    double alfa = 0.05; 


    MatrixXd c(2, 5);
    c << -2.5, -1.5, -0.5, 0.5, 1.5,
         -0.4, -0.2, 0.0, 0.2, 0.4;
    c *= 6.6;


    VectorXd xi1(2);
    xi1 << q(0), dq(0);
    VectorXd h1 = VectorXd::Zero(5);
    
    for (int j = 0; j < 5; ++j) {
        h1(j) = std::exp(-std::pow((xi1 - c.col(j)).norm(), 2) / (2 * std::pow(bj[0], 2)));
    }
    
    double t1_0 = w1.dot(h1);
    VectorXd d_w1 = xite * (taud[0] - t1_0) * h1;
    w1 = w1_1 + d_w1 + alfa * (w1_1 - w1_2);
    w1_2 = w1_1;
    w1_1 = w1;


    VectorXd xi2(2);
    xi2 << q(1), dq(1);
    VectorXd h2 = VectorXd::Zero(5);
    
    for (int j = 0; j < 5; ++j) {
        h2(j) = std::exp(-std::pow((xi2 - c.col(j)).norm(), 2) / (2 * std::pow(bj[1], 2)));
    }
    
    double t2_0 = w2.dot(h2);
    VectorXd d_w2 = xite * (taud[1] - t2_0) * h2;
    w2 = w2_1 + d_w2 + alfa * (w2_1 - w2_2);
    w2_2 = w2_1;
    w2_1 = w2;

    VectorXd xi3(2);
    xi3 << q(2), dq(2);
    VectorXd h3 = VectorXd::Zero(5);
    
    for (int j = 0; j < 5; ++j) {
        h3(j) = std::exp(-std::pow((xi3 - c.col(j)).norm(), 2) / (2 * std::pow(bj[2], 2)));
    }
    
    double t3_0 = w3.dot(h3);
    VectorXd d_w3 = xite * (taud[2] - t3_0) * h3;
    w3 = w3_1 + d_w3 + alfa * (w3_1 - w3_2);
    w3_2 = w3_1;
    w3_1 = w3;


    VectorXd xi4(2);
    xi4 << q(3), dq(3);
    VectorXd h4 = VectorXd::Zero(5);
    
    for (int j = 0; j < 5; ++j) {
        h4(j) = std::exp(-std::pow((xi4 - c.col(j)).norm(), 2) / (2 * std::pow(bj[3], 2)));
    }
    
    double t4_0 = w4.dot(h4);
    VectorXd d_w4 = xite * (taud[3] - t4_0) * h4;
    w4 = w4_1 + d_w4 + alfa * (w4_1 - w4_2);
    w4_2 = w4_1;
    w4_1 = w4;


    VectorXd xi5(2);
    xi5 << q(4), dq(4);
    VectorXd h5 = VectorXd::Zero(5);
    
    for (int j = 0; j < 5; ++j) {
        h5(j) = std::exp(-std::pow((xi5 - c.col(j)).norm(), 2) / (2 * std::pow(bj[4], 2)));
    }
    
    double t5_0 = w5.dot(h5);
    VectorXd d_w5 = xite * (taud[4] - t5_0) * h5;
    w5 = w5_1 + d_w5 + alfa * (w5_1 - w5_2);
    w5_2 = w5_1;
    w5_1 = w5;

    VectorXd xi6(2);
    xi6 << q(5), dq(5);
    VectorXd h6 = VectorXd::Zero(5);
    
    for (int j = 0; j < 5; ++j) {
        h6(j) = std::exp(-std::pow((xi6 - c.col(j)).norm(), 2) / (2 * std::pow(bj[5], 2)));
    }
    
    double t6_0 = w6.dot(h6);
    VectorXd d_w6 = xite * (taud[5] - t6_0) * h6;
    w6 = w6_1 + d_w6 + alfa * (w6_1 - w6_2);
    w6_2 = w6_1;
    w6_1 = w6;

    VectorXd tau_f(6);
    tau_f << t1_0, t2_0, t3_0, t4_0, t5_0, t6_0;
    prev_tau_f = tau_f;

    memcpy(q_pprev, q_prev, sizeof(q_prev));
    for (int i = 0; i < 6; ++i) {
        q_prev[i] = q(i);
        q_pprev[i] = q_prev[i];
    }

    clock_gettime(CLOCK_REALTIME, &st);
    t1 = st.tv_sec + (st.tv_nsec / 1e9);
    dt = t1 - t0;
    return tau_f;

}