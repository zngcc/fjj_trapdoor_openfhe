#include <iostream>
#include <vector>
#include <string>
#include <chrono>

// 性能测量工具类
class Timer {
private:
    std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
public:
    void start() { start_time = std::chrono::high_resolution_clock::now(); }
    double stop() {
        auto end_time = std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double>(end_time - start_time).count();
    }
};


// 打印向量
template<typename T>
void print_vector(const std::vector<T>& vec, const std::string& name = "") {
    if (!name.empty()) {
        std::cout << name << ": ";
    }
    for (const auto& val : vec) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
}

// 计算向量误差
double calculate_max_error(const std::vector<double>& result, 
                          const std::vector<double>& expected) {
    // 这里需要实现误差计算逻辑
    // 例如：
    double max_error = 0.0;
    for (size_t i = 0; i < result.size() && i < expected.size(); ++i) {
        double error = std::abs(result[i] - expected[i]);
        if (error > max_error) {
            max_error = error;
        }
    }
    return max_error;
}
