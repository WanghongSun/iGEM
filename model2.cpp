#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <cstdlib>
using namespace std;
struct Model2Params {
    // 细菌生长参数
    double k = 0.131;     // 细菌生长速率常数
    double K = 100000;     // 环境容纳量
    // 二肽分泌与消耗参数
    double m = 0.04;    // 分泌比例常数
    double kx = 2;    // 二肽消耗速率常数
};
// N: 大肠杆菌数量
// X: 胞内二肽浓度
void computeDerivatives(
    double N, double X,
    const Model2Params& params,
    double& dNdt, double& dXdt
) {
    // 细菌生长方程：dN/dt = k*N*(1 - N/K)
    dNdt = params.k * N * (1.0 - N / params.K);
    // 二肽分泌-消耗方程：dX/dt = m*N - kx*X
    dXdt = params.m * N - params.kx * X;
    N = std::min(N, params.K); 
    if (dNdt < 0 && N <= 0) dNdt = 0;
    if (dXdt < 0 && X <= 0) dXdt = 0;
}
void rk4Step(
    double& N, double& X,
    double dt, const Model2Params& params
) {
    // 计算k1
    double k1_N, k1_X;
    computeDerivatives(N, X, params, k1_N, k1_X);
    // 计算k2
    double k2_N, k2_X;
    computeDerivatives(
        N + 0.5 * dt * k1_N,
        X + 0.5 * dt * k1_X,
        params, k2_N, k2_X
    );
    // 计算k3
    double k3_N, k3_X;
    computeDerivatives(
        N + 0.5 * dt * k2_N,
        X + 0.5 * dt * k2_X,
        params, k3_N, k3_X
    );
    // 计算k4
    double k4_N, k4_X;
    computeDerivatives(
        N + dt * k3_N,
        X + dt * k3_X,
        params, k4_N, k4_X
    );
    N += dt * (k1_N + 2 * k2_N + 2 * k3_N + k4_N) / 6.0;
    X += dt * (k1_X + 2 * k2_X + 2 * k3_X + k4_X) / 6.0;
    N = std::max(0.0, N);
    X = std::max(0.0, X);
}
int main() {
    Model2Params params;
    double N = 5000;     // 初始细菌数量
    double X = 0.0;     // 初始二肽浓度
    // 模拟参数（覆盖完整生长和分泌周期）
    double total_time = 100.0;  // 总模拟时间
    double dt = 0.05;             // 时间步长
    int steps = static_cast<int>(total_time / dt);  // 总迭代步数
    // 输出文件
    std::ofstream outFile("model2_igem_whu_results.txt");
    if (!outFile.is_open()) {
        std::cerr << "错误：无法打开输出文件！" << std::endl;
        return 1;
    }
    outFile << std::fixed << std::setprecision(6);
    outFile << "Time\tN\tX\n";  // 列名：时间、细菌数量、二肽浓度
    for (int i = 0; i < steps; ++i) {
        double t = i * dt;
        outFile << t << "\t" << N << "\t" << X << "\n";
        rk4Step(N, X, dt, params);
    }
    outFile.close();
    system("gnuplot plot_N.gp");
    system("gnuplot plot_X.gp");
    std::cout << "图像生成完成，已保存为N_vs_time.png和X_vs_time.png" << std::endl;
    return 0;
}
