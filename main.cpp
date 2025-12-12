#include <iostream>
#include <ostream>
#include<string>
#include<iomanip>
#include "TSP.h"


int main() {
    std::cout << std::fixed << std::setprecision(15);

    //generateMap(1000, 10000);
    downloadMap("citiesMap.txt");
    std::vector<City> Path = greedyAlgorithm();
    displayPath(Path);
    std::cout << "Droga algorytmu zachlannego wynosi: " << calculateDistance(Path) << std::endl;

    std::cout << "-----------------------------------------------------------------" << std::endl;

    std::vector<City> bestPath = simulatedAnnealing(Path, 10000.0, 0.000001, 0.9999);
    displayPath(bestPath);
    std::cout << "Droga algorytmu metaheurystycznego wynosi: " << calculateDistance(bestPath) << std::endl;

    return 0;
}