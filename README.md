# Тестовое задание в VK Maps

Выполнил **Камлев Виталий**

Файл `BFS.cpp` содержит реализацию алгоритма обхода графа, переданного в заданном формате файлом `.txt`, в ширину (*BFS*) с использованием шаблонов и паттерна *«посетитель»* для вычисления кратчайших путей. Ниже приведено подробное описание каждой составляющей.

---

## 1. Структура данных для ребра

### Шаблон `Edge`
- **Описание:**  
  Шаблонная структура `Edge` принимает тип `Vertex` и представляет ребро графа, хранящее пару вершин (начало и конец).
  
- **Методы:**
  - `GetStart()` – возвращает исходную вершину.
  - `GetTarget()` – возвращает целевую вершину.
  
- **Внутренняя реализация:**  
  Ребро хранится в виде `std::pair<Vertex, Vertex>`, где первая вершина — начало, а вторая — цель.

---

## 2. Класс графа

### Шаблон `Graph`
- **Описание:**  
  Класс `Graph` построен на шаблоне, по умолчанию тип вершины – `int`. Он принимает список ребер и формирует структуру смежности для хранения исходящих ребер для каждой вершины.
  
- **Конструктор:**  
  Принимает вектор ребер, затем для каждого ребра добавляет его в список исходящих ребер вершины-источника, используя `std::unordered_map`.

- **Методы:**
  - `GetOutgoingEdges(const Vertex& vertex)` – возвращает вектор ребер, исходящих из данной вершины.

---

## 3. Обход графа (BFS)

### Функция `BFS`
- **Описание:**  
  Шаблонная функция для обхода графа в ширину. Принимает граф, начальную вершину и объект-посетитель, который собирает результаты обхода.
  
- **Основные шаги:**
  1. Инициализация расстояния до начальной вершины как 0.
  2. Использование очереди для обхода графа. Каждая вершина помечается расстоянием (число шагов от начальной точки).
  3. Для каждого ребра, исходящего из текущей вершины, если найден путь короче ранее зафиксированного, обновляется расстояние и вершина добавляется в очередь.
  4. После завершения обхода вызывается метод `Fill` посетителя, чтобы передать вычисленные расстояния.

---

## 4. Паттерн Посетитель

### Абстрактный класс `AbstractBFSVisitor`
- **Описание:**  
  Определяет интерфейс для посетителей, которые обрабатывают найденные ребра во время BFS. Объявляет виртуальную функцию `DiscoverEdge` для обработки каждого ребра.
  
- **Свойства:**  
  Хранит вектор `answer`, предназначенный для результатов (например, расстояния до каждой вершины).

### Класс `ShortestPathsVisitor`
- **Наследование:**  
  Наследуется от `AbstractBFSVisitor`.
  
- **Цель:**  
  Определяет посетителя, который вычисляет кратчайшие пути от заданной начальной вершины до всех остальных.
  
- **Основные методы:**
  - `DiscoverEdge` – обновляет предков и минимальное расстояние для целевой вершины ребра.
  - `CreateAns` – инициализирует вектор `answer`, заполняя его значениями `std::numeric_limits<int>::max()` и устанавливая начальное расстояние равным 0.
  - `Fill` – обновляет вектор `answer` на основе результатов обхода (полученных в `BFS`).
  - `GetAns` – возвращает расстояние до вершины с заданным индексом.

---

## 5. Вызов обхода и вывод результата

### Функция `CallBFS`
- **Описание:**  
  Функция, объединяющая создание графа, вызов обхода BFS и получение посетителя с вычисленными кратчайшими путями.
  
- **Параметры:**  
  - `all_edges` – вектор всех ребер графа.
  - `vertexes` – количество вершин.
  - `start` – начальная вершина обхода.

- **Возвращаемое значение:**  
  Возвращает объект `ShortestPathsVisitor`, содержащий вычисленные кратчайшие расстояния.

### Функция `Print`
- **Описание:**  
  Выводит кратчайшие расстояния для каждой вершины в консоль.

---

## 6. Чтение данных из файла и обработка ошибок

### Функция `GetData`
- **Описание:**  
  Функция отвечает за считывание входных данных из файла `graph.txt`.
  
- **Процесс:**
  1. Открытие файла и проверка успешности открытия.
  2. Считывание количества вершин и ребер.
  3. Чтение ребер графа. Каждое ребро добавляется в список дважды (для неориентированного графа — в обоих направлениях).
  4. Считывание номера начальной вершины.
  5. Запуск обхода BFS через функцию `CallBFS` и вывод результатов с помощью `Print`.
  6. Обработка исключений и вывод сообщений об ошибках при сбоях чтения или открытия файла.

---

## 7. Функция `main`
- **Описание:**  
  Точка входа в программу, которая вызывает функцию `GetData` для инициализации работы с графом и обработки данных.

---

## Заключение

Данный код демонстрирует:
- Создание обобщённых структур данных для представления графа.
- Реализацию алгоритма BFS для поиска кратчайших путей.
- Применение паттерна «посетитель» для обработки результатов обхода.
- Чтение входных данных из файла с обработкой исключений.

Такой подход позволяет легко расширять функциональность и повторно использовать компоненты при работе с различными типами графов и алгоритмами обхода.
