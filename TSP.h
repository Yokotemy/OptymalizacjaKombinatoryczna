#ifndef LOCALSEARCH_TSP_H
#define LOCALSEARCH_TSP_H

#include <vector>
#include <string>



struct City {
    int id;
    int x;
    int y;
};

void generateMap(int num, int dist);
void downloadMap(const std::string& fileName);
void displayMap();

std::vector<City> greedyAlgorithm();
void displayPath(const std::vector<City>& path);
double calculateDistance(const std::vector<City>& path);

void twoOpt(std::vector<City>& path, int i, int j);
bool hasDuplicate(const std::vector<City>& path);

std::vector<City> simulatedAnnealing(
    const std::vector<City>& path,
    double T_start = 10000.0,
    double T_end   = 0.001,
    double alpha   = 0.999
);

#endif
