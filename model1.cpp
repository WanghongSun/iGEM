#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <cstdlib>
using namespace std;
// 模型参数结构体
struct ModelParams {
    // 生长相关参数
    double k_c1 = 0.56;    // 捕食者(C1)最大生长速率 
    double k_c2 = 0.5;     // 猎物(C2)最大生长速率 
    double c_max = 100.0;  // 最大种群浓度（OD600归一化）
    // 信号分子相关参数
    double k_d1 = 0.10;    // C1分泌3OC12HSL(Ae1)的速率常数 
    double k_d2 = 0.12;    // C2分泌3OC6HSL(Ae2)的速率常数 
    double d1 = 0.08;      // 3OC12HSL(Ae1)的降解速率 
    double d2 = 0.09;      // 3OC6HSL(Ae2)的降解速率 
    // Hill系数参数（调控死亡率）
    double K1 = 5.0;       // C1死亡率的Hill常数（3OC6HSL的结合系数）
    double K2 = 4.0;       // C2死亡率的Hill常数（3OC12HSL的结合系数）
    double beta = 2.0;     // Hill指数
    // 额外调节参数
    double mort_base = 0.05; // 基础死亡率 
    double d = 0.1125;       // 稀释率 
    double iptg = 5.0;       // IPTG浓度（诱导表达效率）
    double dc1 = 20.0;       // CcdB毒性效应系数
};
void computeDerivatives(
    double c1, double c2,       
    double Ae1, double Ae2,     
    const ModelParams& params,
    double& dc1_dt, double& dc2_dt,  
    double& dAe1_dt, double& dAe2_dt 
) {
    double csum = c1 + c2;
    double competition = csum / params.c_max;  
    // 1. 种群生长速率（含稀释率影响）
    double growth_c1 = params.k_c1 * c1 * (1 - competition) - params.d * c1;
    double growth_c2 = params.k_c2 * c2 * (1 - competition) - params.d * c2;
    // 2. 死亡率函数（受IPTG调控）
    double mortality_c1 = (params.mort_base * params.dc1) / 
                        (1 + pow(Ae2 / params.K1, params.beta)) * 
                        (1.0 / (1.0 + params.iptg / 10.0));                  
    double mortality_c2 = params.mort_base * 
                        (1 + pow(Ae1 / params.K2, params.beta)) * 
                        (params.iptg / (5.0 + params.iptg));
    // 3. 种群净增长（生长-死亡）
    dc1_dt = growth_c1 - mortality_c1 * c1;
    dc2_dt = growth_c2 - mortality_c2 * c2;
    // 4. 信号分子动力学（含稀释率影响）
    dAe1_dt = params.k_d1 * c1 - params.d1 * Ae1 - params.d * Ae1;
    dAe2_dt = params.k_d2 * c2 - params.d2 * Ae2 - params.d * Ae2;
    // 确保分母不会导致变量负向溢出
    if (dc1_dt < 0 && c1 <= 0) dc1_dt = 0;
    if (dc2_dt < 0 && c2 <= 0) dc2_dt = 0;
    if (dAe1_dt < 0 && Ae1 <= 0) dAe1_dt = 0;
    if (dAe2_dt < 0 && Ae2 <= 0) dAe2_dt = 0;
}
void rk4Step(
    double& c1, double& c2, double& Ae1, double& Ae2,
    double dt, const ModelParams& params
) {
    double k1_c1, k1_c2, k1_Ae1, k1_Ae2;
    computeDerivatives(c1, c2, Ae1, Ae2, params, k1_c1, k1_c2, k1_Ae1, k1_Ae2);
    double k2_c1, k2_c2, k2_Ae1, k2_Ae2;
    computeDerivatives(
        c1 + 0.5 * dt * k1_c1,
        c2 + 0.5 * dt * k1_c2,
        Ae1 + 0.5 * dt * k1_Ae1,
        Ae2 + 0.5 * dt * k1_Ae2,
        params, k2_c1, k2_c2, k2_Ae1, k2_Ae2
    );
    double k3_c1, k3_c2, k3_Ae1, k3_Ae2;
    computeDerivatives(
        c1 + 0.5 * dt * k2_c1,
        c2 + 0.5 * dt * k2_c2,
        Ae1 + 0.5 * dt * k2_Ae1,
        Ae2 + 0.5 * dt * k2_Ae2,
        params, k3_c1, k3_c2, k3_Ae1, k3_Ae2
    );
    double k4_c1, k4_c2, k4_Ae1, k4_Ae2;
    computeDerivatives(
        c1 + dt * k3_c1,
        c2 + dt * k3_c2,
        Ae1 + dt * k3_Ae1,
        Ae2 + dt * k3_Ae2,
        params, k4_c1, k4_c2, k4_Ae1, k4_Ae2
    );
    c1 += dt * (k1_c1 + 2*k2_c1 + 2*k3_c1 + k4_c1) / 6.0;
    c2 += dt * (k1_c2 + 2*k2_c2 + 2*k3_c2 + k4_c2) / 6.0;
    Ae1 += dt * (k1_Ae1 + 2*k2_Ae1 + 2*k3_Ae1 + k4_Ae1) / 6.0;
    Ae2 += dt * (k1_Ae2 + 2*k2_Ae2 + 2*k3_Ae2 + k4_Ae2) / 6.0;
    // 确保变量非负
    c1 = std::max(0.0, c1);
    c2 = std::max(0.0, c2);
    Ae1 = std::max(0.0, Ae1);
    Ae2 = std::max(0.0, Ae2);
}
int main() {
    ModelParams params;  
    double c1 = 20.0;   // 初始捕食者浓度
    double c2 = 20.0;   // 初始猎物浓度
    double Ae1 = 0.0;   // 初始A1信号浓度
    double Ae2 = 0.0;   // 初始A2信号浓度
    double total_time = 1000.0;  // 总模拟时间
    double dt = 0.1;              // 时间步长
    int steps = static_cast<int>(total_time / dt);  
    std::ofstream outFile("qs_model_results.txt");
    if (!outFile.is_open()) {
        std::cerr << "无法打开输出文件" << std::endl;
        return 1;
    }
    outFile << std::fixed << std::setprecision(6);
    outFile << "Time\tc1\tc2\tAe1\tAe2\n";  
    for (int i = 0; i < steps; ++i) {
        double t = i * dt;  
        outFile << t << "\t" << c1 << "\t" << c2 << "\t" << Ae1 << "\t" << Ae2 << "\n";
        rk4Step(c1, c2, Ae1, Ae2, dt, params);
    }
    outFile.close();
    system("gnuplot plot_c1_c2.gp"); 
    std::cout << "模拟完成，结果保存至qs_model_results.txt，图像已绘制为c1_c2_dynamics.png" << std::endl;
    return 0;
}
