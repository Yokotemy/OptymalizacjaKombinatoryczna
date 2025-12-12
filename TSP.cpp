#include "TSP.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <algorithm>
#include <set>

std::vector<City> cities;
void generateMap(int number, int maxDistance) {
    std::ofstream file("citiesMap.txt");
    if (!file) { std::cerr << "Nie udało się otworzyć pliku!\n"; exit(1); }

    std::mt19937 gen(std::random_device{}());
    std::uniform_int_distribution<int> dist(0, maxDistance);

    for (int i = 0; i < number; i++) {
        file << i + 1 << " " << dist(gen) << " " << dist(gen) << "\n";
    }
    file.close();
    std::cout << "Wygenerowano miasta.\n";
}

void downloadMap(const std::string& fileName) {
    cities.clear();
    std::ifstream file(fileName);
    if (!file) { std::cerr << "Nie udało się otworzyć pliku!\n"; return; }

    int id, x, y;
    while (file >> id >> x >> y) {
        cities.push_back({id, x, y});
    }
    file.close();
}
void displayMap() {
    for (const auto& d : cities)
        std::cout << d.id << " " << d.x << " " << d.y << "\n";
}

double distance(const City& a, const City& b) {
    return std::sqrt(std::pow(a.x - b.x, 2) + std::pow(a.y - b.y, 2));
}
double calculateDistance(const std::vector<City>& path) {
    double dist = 0.0;
    for (size_t i = 0; i < path.size() - 1; ++i)
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

        visited[bestIndex] = true;
        current = cities[bestIndex];
        greedy.push_back(current);
    }

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
bool hasDuplicate(const std::vector<City>& path) {
    std::set<int> ids;
    for (size_t i = 0; i < path.size() - 1; i++) {
        if (!ids.insert(path[i].id).second)
            return true;
    }
    return false;
}
std::vector<City> simulatedAnnealing(
    const std::vector<City>& startPath,
    double T_start,
    double T_end,
    double alpha
) {
    std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<double> rand01(0.0, 1.0);

    int n = startPath.size();
    std::vector<City> best = startPath;
    std::vector<City> current = startPath;

    double bestDist = calculateDistance(best);
    double currentDist = bestDist;
    double T = T_start;

    while (T > T_end) {
        std::uniform_int_distribution<int> distI(1, n - 3);
        int i = distI(gen);
        std::uniform_int_distribution<int> distJ(i + 1, n - 2);
        int j = distJ(gen);

        std::vector<City> candidate = current;
        twoOpt(candidate, i, j);

        if (hasDuplicate(candidate)) {
            std::cerr << "Błąd: duplikaty po 2-opt!\n";
            exit(1);
        }

        double candidateDist = calculateDistance(candidate);
        double delta = candidateDist - currentDist;

        if (delta < 0 || rand01(gen) < std::exp(-delta / T)) {
            current = candidate;
            currentDist = candidateDist;

            if (currentDist < bestDist) {
                best = current;
                bestDist = currentDist;
            }
        }

        T *= alpha;
    }

    return best;
}
