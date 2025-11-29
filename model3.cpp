#include <SFML/Graphics.hpp>
#include <vector>
#include <cmath>
#include <random>
#include <iostream>
#include <map>
#include <string>
#include <algorithm>

// --- 基础参数 ---
const int WIDTH = 1500;
const int HEIGHT = 1000;
const float CELL_RADIUS = 4.0f;
const int FPS = 60;

// --- 颜色定义 ---
const sf::Color COLOR_POS(255, 50, 50);     // 阳性: 鲜红
const sf::Color COLOR_NEG(50, 150, 255);    // 阴性: 鲜蓝
const sf::Color COLOR_NEU(200, 200, 200);   // 中性: 灰白
const sf::Color BG_COLOR(10, 10, 10);

// --- 数据结构 ---
enum ParticleType { POS, NEG, NEU };

struct Particle {
    ParticleType type;
    float x, y;
    bool bonded_state;      // 结合状态：false=未结合，true=已结合
    int bonded_partners;    // 已结合的伙伴数量
    float bond_distance;    // 结合距离阈值
    std::vector<int> bonded_indices;  // 已结合的粒子索引
};


struct ScenarioParams {
    float wall_bias_neu;
    float wall_bias_charged;
    float wall_dist_limit;
    float base_repulsion;
    float pair_attraction;
    float same_repulsion;
    float cross_repulsion;
    int required_bonds_positive;  // 阳性粒子需要结合的伙伴数阈值
    int required_bonds_negative;  // 阴性粒子需要结合的伙伴数阈值
};

struct Scenario {
    std::wstring name; // 使用宽字符支持中文
    float radius_factor;
    std::pair<int, int> ratio;
    ScenarioParams params;
};

// --- 场景配置 (硬编码以匹配 Python 字典) ---
std::map<int, Scenario> SCENARIOS = {
    {1, {L"1. 资源掠夺 (趋壁效应)", 240.0f, {1, 1}, {2.5f, 0.8f, 200.0f, 500.0f, 10.0f, 5.0f, 0.0f, 1, 1}}},
    {2, {L"2. 族群分离 (带隙效应)", 30.0f, {1, 1}, {0.0f, 0.0f, 0.0f, 500.0f, 60.0f, 15.0f, 80.0f, 1, 1}}},
    {3, {L"3. 自由扩散 (小感知半径)", 7.5f, {1, 1}, {0.0f, 0.0f, 0.0f, 500.0f, 5.0f, 2.0f, 0.0f, 1, 1}}},
    {4, {L"4. 比例 1:1 (紧密配对)", 30.0f, {1, 1}, {1.0f, 1.0f, 150.0f, 500.0f, 150.0f, 20.0f, 40.0f, 1, 1}}},
    {5, {L"5. 比例 1:2 (1阴被2阳包围)", 30.0f, {2, 1}, {1.0f, 1.0f, 150.0f, 500.0f, 150.0f, 8.0f, 40.0f, 1, 2}}},
    {6, {L"6. 比例 1:3 (1阴被3阳包围)", 30.0f, {3, 1}, {1.0f, 1.0f, 150.0f, 500.0f, 180.0f, 15.0f, 40.0f, 1, 3}}}
};

// --- 随机数生成器 ---
std::random_device rd;
std::mt19937 gen(rd());

float random_float(float min, float max) {
    std::uniform_real_distribution<float> dis(min, max);
    return dis(gen);
}

int random_int(int min, int max) {
    std::uniform_int_distribution<> dis(min, max);
    return dis(gen);
}

// --- 核心逻辑 ---

float get_distance_sq(const Particle& p1, const Particle& p2) {
    return (p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y);
}

void update_bonding_states(std::vector<Particle>& particles, int scenario_id) {
    // 重置所有粒子的结合伙伴计数和索引，但保留结合状态
    for (auto& p : particles) {
        p.bonded_partners = 0;
        p.bonded_indices.clear();
    }
    
    const ScenarioParams& params = SCENARIOS[scenario_id].params;
    float bond_distance_sq = std::pow(CELL_RADIUS * 2.5f, 2); // 结合距离阈值的平方
    
    // 检查所有粒子对之间的距离
    for (size_t i = 0; i < particles.size(); ++i) {
        Particle& p1 = particles[i];
        if (p1.type == NEU) continue; // 中性粒子不参与结合
        
        for (size_t j = i + 1; j < particles.size(); ++j) {
            Particle& p2 = particles[j];
            if (p2.type == NEU) continue; // 中性粒子不参与结合
            
            // 只有异性电荷粒子才能结合
            if (p1.type != p2.type) {
                float dist_sq = get_distance_sq(p1, p2);
                
                // 如果粒子距离小于结合阈值
                if (dist_sq < bond_distance_sq) {
                    // 更新结合伙伴计数
                    p1.bonded_partners++;
                    p2.bonded_partners++;
                    
                    // 记录结合的粒子索引
                    p1.bonded_indices.push_back(j);
                    p2.bonded_indices.push_back(i);
                }
            }
        }
    }
    
    // 根据结合伙伴数量更新结合状态
    for (auto& p : particles) {
        if (p.type == NEU) continue; // 中性粒子没有结合状态
        
        int required_bonds = 0;
        if (p.type == POS) {
            required_bonds = params.required_bonds_positive;
        } else if (p.type == NEG) {
            required_bonds = params.required_bonds_negative;
        }
        
        // 如果结合伙伴数量达到或超过阈值，设置为已结合状态
        if (p.bonded_partners >= required_bonds) {
            p.bonded_state = true;
        } else {
            p.bonded_state = false;
        }
    }
}

float calculate_fitness(const Particle& particle, const std::vector<Particle>& particles, int p_idx, float sensing_radius, int scenario_id) {
    float score = 0.0f;
    float sensing_radius_sq = sensing_radius * sensing_radius;
    const ScenarioParams& params = SCENARIOS[scenario_id].params;

    // 基础距离定义
    float R = CELL_RADIUS;
    float hard_shell_sq = std::pow(2.1f * R, 2);
    float repulsion_range_sq = std::pow(12.0f * R, 2);

    // --- 1. 墙壁效应 ---
    if (params.wall_bias_neu > 0 || params.wall_bias_charged > 0) {
        float d_left = particle.x;
        float d_right = WIDTH - particle.x;
        float d_top = particle.y;
        float d_bottom = HEIGHT - particle.y;
        float min_wall_dist = std::min({d_left, d_right, d_top, d_bottom});

        float limit = params.wall_dist_limit;
        if (min_wall_dist < limit) {
            float factor = 0;
            if (particle.type == NEU) factor = params.wall_bias_neu;
            else factor = params.wall_bias_charged;
            score += (limit - min_wall_dist) * factor;
        }
    }

    // --- 2. 粒子间相互作用 ---
    for (size_t i = 0; i < particles.size(); ++i) {
        if (i == p_idx) continue; // 跳过自己

        const Particle& other = particles[i];
        float dist_sq = get_distance_sq(particle, other);

        // 对于场景1，特殊处理：保留阴阳粒子结合，跳过其他距离限制
        if (scenario_id == 1) {
            // 只保留阴阳粒子间的吸引力
            bool is_p_pos = (particle.type == POS);
            bool is_p_neg = (particle.type == NEG);
            bool is_o_pos = (other.type == POS);
            bool is_o_neg = (other.type == NEG);
            
            // 只处理阳性和阴性粒子之间的吸引力
            if ((is_p_pos && is_o_neg) || (is_p_neg && is_o_pos)) {
                if (params.pair_attraction > 0 && !particle.bonded_state && !other.bonded_state) {
                    // 进一步扩大作用范围到24.0f * R，显著增加寻路范围
                    if (dist_sq < std::pow(50.0f * R, 2)) {
                        // 进一步增加吸引力强度的系数，从2000.0f增加到4000.0f
                        score += params.pair_attraction * (1.0f / (dist_sq + 1.0f)) * 8000.0f;
                    }
                }
            }
            // 跳过所有其他距离限制检查
            continue;
        }

        // A. 绝对刚体排斥 - 增强版本
        if (dist_sq < hard_shell_sq) {
            // 使用随距离减小而指数增长的排斥力，防止粒子重合卡住
            // 当距离趋近于0时，排斥力趋近于无穷大
            float min_distance_sq = std::pow(2.0f * R, 2); // 粒子刚好接触时的距离平方
            if (dist_sq < min_distance_sq) {
                // 当粒子重叠时，使用更强的排斥力
                float overlap_factor = 1.0f - (dist_sq / min_distance_sq);
                // 指数增长的排斥力，系数可以根据需要调整
                score -= params.base_repulsion * (1.0f + overlap_factor * 1000.0f);
            } else {
                // 当粒子接近但未重叠时，使用基础排斥力
                score -= params.base_repulsion * 10.0f;
            }
            continue;
        }

        if (dist_sq > sensing_radius_sq) continue;

        bool is_p_charged = (particle.type == POS || particle.type == NEG);
        bool is_o_charged = (other.type == POS || other.type == NEG);

        // B. 带电 vs 中性
        if (is_p_charged != is_o_charged) {
            if (params.cross_repulsion > 0) {
                float expanded_repulsion_range_sq = std::pow(20.0f * R, 2); // 12.0f * R * 5 = 60.0f * R
                if (dist_sq < expanded_repulsion_range_sq) {
                    score -= params.cross_repulsion * (1.0f - dist_sq / expanded_repulsion_range_sq) * 20.0f;
                }
            }
        }
        // C. 带电 vs 带电
        else if (is_p_charged && is_o_charged) {
            if (particle.type == other.type) {
                // 同性排斥
                if (params.same_repulsion > 0) {
                    if (dist_sq < repulsion_range_sq) {
                        score -= params.same_repulsion * (1.0f - dist_sq / repulsion_range_sq) * 80.0f;
                    }
                }
            } else {
                // 异性吸引
                if (params.pair_attraction > 0) {
                    // 检查粒子是否已结合 - 如果任意一个已结合，则不计算吸引力
                    if (!particle.bonded_state && !other.bonded_state) {
                        // 进一步扩大作用范围到24.0f * R，显著增加寻路范围
                        if (dist_sq < std::pow(50.0f * R, 2)) {
                            // 进一步增加吸引力强度的系数，从2000.0f增加到4000.0f
                            score += params.pair_attraction * (1.0f / (dist_sq + 1.0f)) * 8000.0f;
                        }
                    }
                }
            }
        }
        // D. 中性 vs 中性
        else if (!is_p_charged && !is_o_charged) {
            // 中性粒子之间的排斥力，保持适当距离
            float neu_repulsion_range_sq = std::pow(10.0f * R, 2);
            if (dist_sq < neu_repulsion_range_sq) {
                // 使用适当的排斥力系数，与同性排斥类似但强度可能需要调整
                score -= params.same_repulsion * (1.0f - dist_sq / neu_repulsion_range_sq) * 60.0f;
            }
        }
    }
    return score;
}

void init_particles(int scenario_id, std::vector<Particle>& particles, float& sensing_radius) {
    particles.clear();
    const Scenario& config = SCENARIOS[scenario_id];
    
    int total_p, num_neu;
    if (scenario_id == 1) { total_p = 600; num_neu = 300; }
    else if (scenario_id == 2) { total_p = 400; num_neu = 200; }
    // 为配对模型（4,5,6）增加中性粒子数量
    else if (scenario_id == 4 || scenario_id == 5 || scenario_id == 6) { 
        total_p = 850; // 增加总粒子数
        num_neu = 300; // 显著增加中性粒子数量，从50个增加到200个
    }
    else { total_p = 350; num_neu = 50; } // 其他模型保持不变

    int remaining = total_p - num_neu;
    float r_pos = (float)config.ratio.first;
    float r_neg = (float)config.ratio.second;
    float total_ratio = r_pos + r_neg;

    // 优化粒子数量比例，确保异性粒子能够正确配对
    // 对于非1:1比例的场景，保证数量较少的类型有足够的粒子
    int num_pos, num_neg;
    if (scenario_id == 4 || scenario_id == 5 || scenario_id == 6) { // 这些场景专为配对设计
        // 确保数量较少的类型有足够的粒子与数量较多的类型配对
        if (r_pos > r_neg) { // 阳多阴少
            num_neg = std::max(30, (int)(remaining / (r_pos + r_neg) * r_neg)); // 至少30个阴性粒子
            num_pos = remaining - num_neg;
        } else if (r_neg > r_pos) { // 阴多阳少
            num_pos = std::max(30, (int)(remaining / (r_pos + r_neg) * r_pos)); // 至少30个阳性粒子
            num_neg = remaining - num_pos;
        } else { // 1:1比例
            num_pos = (int)(remaining * 0.5f);
            num_neg = remaining - num_pos;
        }
    } else {
        // 其他场景保持原有计算方式
        num_pos = (int)(remaining * (r_pos / total_ratio));
        num_neg = remaining - num_pos;
    }
    
    // 确保至少有10个阳性和10个阴性粒子，保证配对可能性
    num_pos = std::max(10, num_pos);
    num_neg = std::max(10, num_neg);
    
    // 重新平衡总数，确保不超过remaining
    int total_charged = num_pos + num_neg;
    if (total_charged > remaining) {
        float scale_factor = (float)remaining / (float)total_charged;
        num_pos = (int)(num_pos * scale_factor);
        num_neg = remaining - num_pos;
    }

    auto add_particles = [&](int count, ParticleType type) {
        for (int i = 0; i < count; ++i) {
            Particle p;
            p.type = type;
            p.x = random_float(50, WIDTH - 50);
            p.y = random_float(50, HEIGHT - 50);
            p.bonded_state = false;      // 初始未结合
            p.bonded_partners = 0;       // 初始结合伙伴数为0
            p.bond_distance = CELL_RADIUS * 2.5f;  // 结合距离阈值
            // 初始化结合索引向量
            particles.push_back(p);
        }
    };

    add_particles(num_pos, POS);
    add_particles(num_neg, NEG);
    add_particles(num_neu, NEU);

    // 计算场景对角线长度，确保感知半径足够大以覆盖整个场景
    float diagonal = std::sqrt(WIDTH * WIDTH + HEIGHT * HEIGHT);
    sensing_radius = diagonal / 2.0f; // 设置感知半径为对角线的一半，确保覆盖整个场景
}

int main() {
    // 启用抗锯齿
    sf::ContextSettings settings;
    settings.antialiasingLevel = 8;
    sf::RenderWindow window(sf::VideoMode(WIDTH, HEIGHT), "Model3 Particle Animation (C++)", sf::Style::Default, settings);
    window.setFramerateLimit(FPS);

    // 加载字体 (Windows 默认路径)
    sf::Font font;
    if (!font.loadFromFile("C:/Windows/Fonts/simhei.ttf")) {
        // 如果找不到黑体，尝试 Arial
        if (!font.loadFromFile("C:/Windows/Fonts/arial.ttf")) {
            std::cerr << "Error loading font" << std::endl;
            return -1;
        }
    }

    int current_scenario = 1;
    std::vector<Particle> particles;
    float sensing_radius;
    init_particles(current_scenario, particles, sensing_radius);

    float T = 20.0f;
    long long iteration = 0;

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
            
            if (event.type == sf::Event::KeyPressed) {
                if (event.key.code >= sf::Keyboard::Num1 && event.key.code <= sf::Keyboard::Num6) {
                    int idx = event.key.code - sf::Keyboard::Num0;
                    if (SCENARIOS.count(idx)) {
                        current_scenario = idx;
                        init_particles(current_scenario, particles, sensing_radius);
                        T = 20.0f;
                        iteration = 0;
                    }
                }
            }
        }

        // --- 蒙特卡洛迭代 ---
        int steps_per_frame = 1000;
        std::uniform_int_distribution<> p_dist(0, particles.size() - 1);

        for (int k = 0; k < steps_per_frame; ++k) {
            int idx = p_dist(gen);
            Particle& p = particles[idx];

            float old_x = p.x;
            float old_y = p.y;
            float old_score = calculate_fitness(p, particles, idx, sensing_radius, current_scenario);

            float step = 5.0f;
            float dx = random_float(-step, step);
            float dy = random_float(-step, step);

            p.x += dx;
            p.y += dy;

            // 边界限制
            p.x = std::max(CELL_RADIUS, std::min((float)WIDTH - CELL_RADIUS, p.x));
            p.y = std::max(CELL_RADIUS, std::min((float)HEIGHT - CELL_RADIUS, p.y));

            float new_score = calculate_fitness(p, particles, idx, sensing_radius, current_scenario);
            float delta_S = new_score - old_score;

            bool accept = false;
            if (delta_S > 0) {
                accept = true;
            } else {
                if (T > 0.01f) {
                    if (random_float(0.0f, 1.0f) < std::exp(delta_S / T)) {
                        accept = true;
                    }
                }
            }

            if (!accept) {
                p.x = old_x;
                p.y = old_y;
            }
        }

        // 温度衰减
        if (T > 0.1f) T *= 0.999f;
        iteration += steps_per_frame;
        
        // 更新粒子结合状态
        update_bonding_states(particles, current_scenario);

        // --- 绘制 ---
        window.clear(BG_COLOR);

        // 绘制粒子
        sf::CircleShape circle(CELL_RADIUS);
        circle.setOrigin(CELL_RADIUS, CELL_RADIUS); // 设置圆心为原点

        for (const auto& p : particles) {
            circle.setPosition(p.x, p.y);
            // 保持原有颜色，不随结合状态改变
            if (p.type == POS) circle.setFillColor(COLOR_POS);
            else if (p.type == NEG) circle.setFillColor(COLOR_NEG);
            else circle.setFillColor(COLOR_NEU);
            
            // 移除结合状态对视觉效果的影响，所有粒子使用相同的描边设置
            circle.setOutlineThickness(-1.0f); // 向内描边
            circle.setOutlineColor(circle.getFillColor()); // 保持同色
            
            window.draw(circle);
        }

        // UI 显示
        const Scenario& info = SCENARIOS[current_scenario];
        
        sf::Text text;
        text.setFont(font);
        text.setCharacterSize(24);
        text.setFillColor(sf::Color::White);

        // 第一行
        text.setString(info.name);
        text.setPosition(20, 20);
        window.draw(text);

        // 第二行
        text.setString(L"感知半径: " + std::to_wstring((int)info.radius_factor) + L" * R");
        text.setFillColor(sf::Color(200, 200, 100));
        text.setPosition(20, 50);
        window.draw(text);

        // 第三行
        // std::to_string 对于浮点数格式化不如 printf 灵活，这里简单处理
        char buffer[100];
        sprintf_s(buffer, "迭代: %lld | 温度 T: %.2f", iteration, T);
        // 将 char* 转换为 wstring 以支持中文环境下的统一显示
        std::string s_buff(buffer);
        std::wstring ws_buff(s_buff.begin(), s_buff.end());
        
        text.setString(ws_buff);
        text.setFillColor(sf::Color(150, 150, 150));
        text.setPosition(20, 80);
        window.draw(text);

        // 第四行
        text.setString(L"按键 1-6 切换场景");
        text.setFillColor(sf::Color(100, 100, 100));
        text.setPosition(20, 110);
        window.draw(text);

        window.display();
    }

    return 0;
}
