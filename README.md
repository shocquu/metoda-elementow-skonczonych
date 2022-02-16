# Metoda Elementów Skończonych - rozkład ciepła

## Description

Zbiór algorytmów i metod umożliwiających obliczenie temperatury węzłów w elemencie skończonym. Program korzysta z dwupunktowego schematu całkowania do numerycznego wyznaczenia wartości.



## Getting Started

### Zależności

Do graficznej reprezentacji danych wymagana jest biblioteka [Matplotlib for C++](https://matplotlib-cpp.readthedocs.io/en/latest/)



### Uruchomienie

Program można uruchomić poprzez utworzenie siatki o ustalonych wymiarach i liczbie węzłów oraz przekazaniu parametrów typu SimulationData do metody *start()*. Argumenty mogą być nadpisane lub pozostawione bez zmian; w tym przypadku symulacja odbędzie się dla domyślnych wartości.
```cpp
SimulationData dataset;
Grid(0.1, 0.1, 4, 4).start(dataset);
```
Inną opcją wczytania danych jest odczyt z pliku. Ścieżkę do pliku tekstowego należy podać bezpośrednio w konstruktorze lub jako argument dla metody *launch()* wywołanej na obiekcie struktury *Grid*. Przy takim uruchomieniu nie ma konieczności przekazywania parametrów symulacji.
```cpp
Grid().launch("data/Test4_31_31_trapez.txt");
Grid("data/Test4_31_31_trapez.txt").start()
```
Ponadto możliwe jest wyświetlenie rozkładu ciepła stosując metodę plotHeatMap().
```cpp
Grid("data/MES_31_31_v2.txt").heatMap();
```



### Podgląd
![Console](https://github.com/shocquu/metoda-elementow-skonczonych/blob/master/output/results.png?raw=true)

![Plot](https://github.com/shocquu/metoda-elementow-skonczonych/blob/master/output/heatmap.png?raw=true)
