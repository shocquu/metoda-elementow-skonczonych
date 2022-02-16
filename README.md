
# Project Title

Metoda Elementów Skończonych - rozkład ciepła

## Description

Zbiór algorytmów i metod umożliwiających obliczenie temperatury węzłów w elemencie skończonym. Program korzysta z dwupunktowego schematu całkowania do numerycznego wyznaczenia wartości.

## Getting Started

### Zależności

Do graficznej reprezentacji danych wymagana jest biblioteka ![Matplotlib for C++](https://matplotlib-cpp.readthedocs.io/en/latest/)

### Uruchomienie

Program można uruchomić tworząc siatkę o ustalonych wymiarach i liczbę węzłów oraz przekazując parametry typu SimulationData. Mogą być one nadpisane lub pozostawione bez zmian; w tym przypadku symulacja odbędzie się dla domyślnych wartości. 
```cpp
SimulationData dataset;
Grid(0.1, 0.1, 4, 4).start(dataset);
```
Możliwi jest także odczyt z pliku. Wtedy nie ma konieczności przekazywania parametrów symulacji.
```cpp
Grid("data/Test4_31_31_trapez.txt").start()
```
Ponadto możliwe jest wyświetlenie rozkładu ciepła stosując metodę plotHeatMap().
```cpp
Grid("data/MES_31_31_v2.txt").heatMap();
```

### Podgląd
![Console](https://github.com/shocquu/metoda-elementow-skonczonych/blob/master/output/results.png?raw=true)

![Plot](https://github.com/shocquu/metoda-elementow-skonczonych/blob/master/output/heatmap.png?raw=true)
