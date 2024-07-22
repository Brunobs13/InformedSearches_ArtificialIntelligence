# Informed Searches in Artificial Intelligence

## Description

This repository contains the solution for a project in the course "Introdução à Inteligência Artificial". The goal is to minimize the movement cost of families within specific territories by strategically placing stations. The solution implements informed search algorithms, specifically the A* algorithm, to find the optimal placement of these stations.

## Project Structure

This project includes the following files:
- `mainef.py`: Initial commit for EfolioB2201529BrunoFerreiraIIA.
- `mapas.py`: Contains the maps and coordinates for the project.

## Instructions for Running the Code

### Prerequisites

Ensure you have Python installed on your system. The project was developed using Python 3.x.

### Running the Python Code

To run the Python code, navigate to the directory containing the files and execute the following command:

```bash
python3 mainf.py
```
You will be prompted to choose the algorithm by name and whether to use specific coordinates or generate random positions.

## Implementation Details

### Problem Analysis

The objective is to position stations in a way that minimizes the movement cost for families in specific territories. The territories are represented by matrices of various dimensions, with families located in zones divided by squares. The minimum distance calculation from a zone to a station is crucial, and for this problem, the Chebyshev distance is used.

### Heuristic Function

An admissible heuristic function penalizes the cost of moving families to a station, subtracting the weight of the family set within a distance of 1 (where the movement cost is zero). The weight assigned to the number of families in this zone is significant, as high population density within this radius could reduce the number of required stations.

### Code Structure

The code is defined in two files:

- **Static Coordinates**: Contains maps and coordinates initialized statically for testing.
- **Random Position Generator**: Generates random positions based on the number of stations defined in the main file.

### Key Functions

- **chebyshev_distance**: Calculates the minimum distance considering horizontal, vertical, and diagonal movements.
- **travel_cost(distance)**: Returns the cost based on the specified distance.
- **Heuristic function**: Calculates the penalty of the movement cost minus the weight of families within a distance of 1 from the current station, avoiding negative numbers by using a 30% weight.
- **calculate_travel_cost(map, stations)**: Calculates and returns the average travel cost.
- **minimize_cost(map, stations)**: Formula to minimize the cost for the problem.
- **a_star_algorithm(map, initial_stations)**: Implements the A* algorithm, working in conjunction with **expand_and_generate_states(map, current_stations, visited, frontier, counter)** to apply the state expansion logic of A*.

### Results Table

The results of the A* algorithm applied to various instances are summarized in the table below:

| Algorithm/Configuration | Instance 1 | Instance 2 | ... | Instance 20 |
|-------------------------|------------|------------|-----|-------------|
| **A* Expansions**       | 1          | 1          | ... | 17          |
| **Generations**         | 1          | 1          | ... | 6473        |
| **Evaluations**         | 1          | 1          | ... | 6473        |
| **Cost**                | 1264       | 1134       | ... | 4295        |
| **Time (msec)**         | 0.09       | 0.05       | ... | 1342.47     |
| **Best result**         | 1114       | 1171       | ... | 3273        |


### Conclusion

This project demonstrates the effective application of the A* search algorithm to minimize the movement cost for families in various territories. The heuristic function and state expansion logic are crucial in achieving optimal results.


