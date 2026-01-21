#include <iostream>
#include <iomanip>
#include <vector>
#include <thread>
#include <mutex>
#include <chrono>
#include <algorithm>
#include "TSP.h"

std::mutex ioMutex;

struct Result {
    double Ts;
    double Te;
    double a;

    std::vector<double> distances;
    double mean;
    double median;
    double timeSeconds;
};

double mean(const std::vector<double>& v) {
    double sum = 0.0;
    for (double x : v) sum += x;
    return sum / v.size();
}

double median(std::vector<double> v) {
    std::sort(v.begin(), v.end());
    size_t n = v.size();
    if (n % 2 == 0)
        return (v[n/2 - 1] + v[n/2]) / 2.0;
    return v[n/2];
}


int main() {


    double Ts = 15000.0;
    double Te = 0.0001;
    double time  = 180.0;
    const int runs = 20;

    std::cout << std::fixed << std::setprecision(15);

    //generateMap(250,10000);
    downloadMap("Mapa.txt");

    std::vector<City> basePath = greedyAlgorithm();

    std::cout << "Greedy path:\n";
    displayPath(basePath);
    std::cout << "Greedy distance: "
              << calculateDistance(basePath) << "\n";

    std::cout << "-------------------------------------------------\n";

    std::vector<std::thread> threads;
    std::vector<std::vector<City>> results(runs);

    auto worker = [&](int id) {
        auto localPath = basePath; // KOPIA
        results[id] = simulatedAnnealing(
            localPath,
                Ts,
                Te,
                time
        );

        std::lock_guard<std::mutex> lock(ioMutex);
        std::cout << "[Thread " << id << "] finished\n";
    };

    for (int i = 0; i < runs; ++i) {
        threads.emplace_back(worker, i);
    }

    for (auto& t : threads) {
        t.join();
    }

    std::cout << "\n========== RESULTS ==========\n";

    for (int i = 0; i < runs; ++i) {
        std::cout << "\nRun " << i << ":\n";
        std::cout << "Distance: "
                  << calculateDistance(results[i]) << "\n";
    }
  /*
  std::cout << std::fixed << std::setprecision(6);

    downloadMap("Mapa.txt");
    std::vector<City> basePath = greedyAlgorithm();

    const int runs = 10;

    std::vector<double> Ts_values = {
        1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000,
        15000, 20000, 25000, 30000, 40000, 50000, 60000, 70000, 80000, 900000, 100000, 110000, 120000
    };
    std::vector<double> Te_values = { 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12 };
    std::vector<double> a_values  = { 0.995, 0.997, 0.999, 0.9995, 0.9997, 0.9999, 0.99995, 0.99997, 0.99999 };

    std::vector<Result> results;

    size_t totalCombinations = Ts_values.size() * Te_values.size() * a_values.size();
    size_t completedCombinations = 0;
    std::mutex progressMutex;

    for (double Ts : Ts_values)
        for (double Te : Te_values)
            for (double a  : a_values) {

                Result res;
                res.Ts = Ts;
                res.Te = Te;
                res.a  = a;
                res.distances.resize(runs);

                auto start = std::chrono::high_resolution_clock::now();

                std::vector<std::thread> threads;

                for (int i = 0; i < runs; ++i) {
                    threads.emplace_back([&, i]() {
                        auto path = basePath;
                        auto finalPath = simulatedAnnealing(path, Ts, Te, a);
                        res.distances[i] = calculateDistance(finalPath);
                    });
                }

                for (auto& t : threads)
                    t.join();

                auto end = std::chrono::high_resolution_clock::now();

                res.mean   = mean(res.distances);
                res.median = median(res.distances);
                res.timeSeconds =
                    std::chrono::duration<double>(end - start).count();

                results.push_back(res);

                // 🔹 wyświetlenie postępu
                {
                    std::lock_guard<std::mutex> lock(progressMutex);
                    completedCombinations++;
                    std::cout << "[COMPLETED] Ts=" << Ts
                              << " Te=" << Te
                              << " a="  << a
                              << " | remaining: "
                              << (totalCombinations - completedCombinations)
                              << " / " << totalCombinations
                              << "\n";
                }
            }

    // sortowanie wyników po średniej rosnąco
    std::sort(results.begin(), results.end(),
        [](const Result& a, const Result& b) {
            return a.mean < b.mean;
        });

    // wypisanie wyników
    std::cout << "\n=========== BENCHMARK RESULTS (SORTED) ===========\n";

    for (const auto& r : results) {
        std::cout << "\nTs=" << r.Ts
                  << " Te=" << r.Te
                  << " a="  << r.a << "\n";

        std::cout << "Mean:    " << r.mean << "\n";
        std::cout << "Median:  " << r.median << "\n";
        std::cout << "Time[s]: " << r.timeSeconds << "\n";

        std::cout << "Distances:\n";
        for (double d : r.distances)
            std::cout << "  " << d << "\n";
    }
    */
    return 0;
}
