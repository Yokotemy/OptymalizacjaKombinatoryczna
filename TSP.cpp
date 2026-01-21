#include "TSP.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <algorithm>
#include <chrono>
#include <set>
#include <vector>

std::vector<City> cities;

void generateMap(int number, int maxDistance, const std::string& fileName = "Mapa.txt") {
    std::ofstream file(fileName);
    if (!file) {
        std::cerr << "Nie udało się otworzyć pliku!\n";
        exit(1);
    }

    std::mt19937 gen(std::random_device{}());
    std::uniform_int_distribution<int> dist(0, maxDistance);

    for (int i = 0; i < number; i++) {
        file << i + 1 << " " << dist(gen) << " " << dist(gen) << "\n";
    }
    file.close();
    std::cout << "Wygenerowano " << number << " miast do pliku: " << fileName << "\n";
}

void downloadMap(const std::string& fileName) {
    cities.clear();
    std::ifstream file(fileName);
    if (!file) {
        std::cerr << "Nie udało się otworzyć pliku: " << fileName << "\n";
        return;
    }
    int id;
    double x, y;
    while (file >> id >> x >> y) {
        cities.push_back({id, x, y});
    }
    file.close();
}

double distance(const City& a, const City& b) {
    return std::sqrt(std::pow(a.x - b.x, 2) + std::pow(a.y - b.y, 2));
}

double calculateDistance(const std::vector<City>& path) {
    double dist = 0.0;
    for (size_t i = 0; i + 1 < path.size(); ++i)
        dist += distance(path[i], path[i+1]);
    return dist;
}

std::vector<City> greedyAlgorithm() {
    std::vector<City> greedy;
    if (cities.empty()) return greedy;

    City current = cities[0];
    greedy.push_back(current);
    std::vector<bool> visited(cities.size(), false);
    visited[0] = true;

    for (size_t step = 1; step < cities.size(); step++) {
        double bestDist = 1e23;
        int bestIndex = -1;

        for (size_t i = 0; i < cities.size(); i++) {
            if (!visited[i]) {
                double d = distance(current, cities[i]);
                if (d < bestDist) {
                    bestDist = d;
                    bestIndex = i;
                }
            }
        }
        if (bestIndex == -1) break;
        visited[bestIndex] = true;
        current = cities[bestIndex];
        greedy.push_back(current);
    }
    if (!greedy.empty())
        greedy.push_back(greedy[0]);

    return greedy;
}

void displayPath(const std::vector<City>& path) {
    for (size_t i = 0; i < path.size(); ++i) {
        std::cout << path[i].id;
        if (i != path.size() - 1) std::cout << "->";
    }
    std::cout << "\n";
}

void twoOpt(std::vector<City>& path, int i, int j) {
    std::reverse(path.begin() + i, path.begin() + j + 1);
}

// Szybka kalkulacja delty - klucz do 3-minutowej optymalizacji
inline double getDelta(const std::vector<City>& path, int i, int j) {
    const City& A = path[i - 1];
    const City& B = path[i];
    const City& C = path[j];
    const City& D = path[j + 1];

    double distBefore = distance(A, B) + distance(C, D);
    double distAfter  = distance(A, C) + distance(B, D);
    return distAfter - distBefore;
}

std::vector<City> simulatedAnnealing(const std::vector<City>& startPath,
    double T_start,
    double T_end,
    double limitSeconds) {

    std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<double> rand01(0.0, 1.0);

    int n = startPath.size();
    std::vector<City> current = startPath;
    std::vector<City> best = startPath;
    double currentDist = calculateDistance(current);
    double bestDist = currentDist;

    auto startTime = std::chrono::steady_clock::now();
    std::uniform_int_distribution<int> distIdx(1, n - 2);

    long long iterations = 0;
    double T = T_start;

    while (true) {
        // 1. Aktualizacja parametrów co 65536 iteracji
        if ((iterations & 65535) == 0) {
            auto now = std::chrono::steady_clock::now();
            double elapsed = std::chrono::duration<double>(now - startTime).count();

            if (elapsed >= limitSeconds) break;

            double progress = elapsed / limitSeconds;

            if (n >= 400) {
                // --- MECHANIZM REHEAT DLA DUŻYCH ZBIORÓW ---
                // Jeśli jesteśmy między 80% a 85% czasu, podbijamy temperaturę (faza "ratunkowa")
                if (progress > 0.80 && progress < 0.85) {
                    // Reheat: liniowo wracamy na chwilę do wyższej temperatury
                    double reheatProgress = (progress - 0.80) / 0.05;
                    T = T_start * 0.1 * (1.0 - reheatProgress); // Skok do 10% T_start i powolny spadek
                } else {
                    // Standardowe chłodzenie potęgowe (skalowane do n)
                    double power = (n >= 900) ? 6.0 : 4.0;
                    T = T_start * std::pow(1.0 - progress, power);
                }
            } else {
                // Dla małych zbiorów - czyste chłodzenie wykładnicze
                T = T_start * std::pow(T_end / T_start, progress);
            }
        }

        // 2. Rdzeń algorytmu (2-opt)
        int i = distIdx(gen);
        int j = distIdx(gen);

        if (i != j) {
            if (i > j) std::swap(i, j);

            // Szybka delta O(1)
            double delta = getDelta(current, i, j);

            if (delta < 0 || (T > 0 && rand01(gen) < std::exp(-delta / T))) {
                std::reverse(current.begin() + i, current.begin() + j + 1);
                currentDist += delta;

                if (currentDist < bestDist) {
                    bestDist = currentDist;
                    best = current;
                }
            }
        }
        iterations++;
    }
    return best;
}
