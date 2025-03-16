# Team 01 - The Pear Project

## Introduction

As part of the Project in Mathematical Engineering, we were tasked to predict internal oxygen and carbon dioxyde concentrations of a pear. To do so, we are using the Finite Element Method (FEM) to implement a respiration-diffusion system. The implemented system can then be tested depending on the storage conditions of the pear. Those are linked to a specific number as given below:


| Number | Storage Condition |
| :---:  | --------          |
|  1     | Orchard           |
|  2     | Shelf life        |
|  3     | Refrigerator      |
|  4     | Precooling        |
|  5     | Disorder inducing |
|  6     | Optimal CA        |

The rest of the parameters, as well as the complete project description, can be found in the Literature folder.


## Running the code on your own

As the Eigen library is already included in this project, there are no specific prerequisites to be installed. However, the Matlab part of the code will require installing Matlab or using the online Matlab editor to retrieve a graphical interface.

After cloning this repository, the code can be run as follows:
```
g++ -o test main.cpp
./test n
```
with *n* the number between 1 and 6 referring to the storage conditions in the table above.

## Contents of the different files

In this section, we analyse in a bit more details how exactly each file in the `src/` folder works.
* **`parameters.cpp`** defines the initial parameters of the problem. It also initializes the parameters that depend on the storage conditions.
* **`read_file.hpp`** creates coherent mesh matrices based on the files with the mesh nodes as input.
* **`jacobian.hpp`** creates the vectors dependent on the concentrations. It also computes the necessary for the Jacobian, that is the matrices containing the derivations of those vectors. 
* **`matrix.cpp`** computes the system matrices that were described in the mathematical computations.
* **`newton_raphson.cpp`** contains the Newton-Raphson method for solving a non-linear problem.
* **`main.cpp`** contains the solution of both the linear and the non-linear problem, calling the different functions in the previous files to create a solver for the given problem. It also contains a timer implemented to retrieve the time spent on the solution of the problem.




## Results

After running the lines above, different information will be displayed in the terminal, like the residue after each iteration of the Newton-Raphson algorithm or the total time after which a sufficient tolerance is reached.
![Cu_0_6 (1)](https://github.com/user-attachments/assets/42ba591e-7644-4b4b-bbe4-c6f8540e4e09)
![Cv_0_6 (1)](https://github.com/user-attachments/assets/ed85e8b7-e0f7-4418-9afd-6c730cf164a0)
![Cu_Final_6 (1)](https://github.com/user-attachments/assets/621338f5-36de-42a9-9dcd-354b2236ea87)
![Cu_0_6 (1)](https://github.com/user-attachments/assets/cdb3c9d3-44aa-4693-a998-57eff08011ab)




In the file **`src/result_c0.txt`**, the initial vector is displayed. This vector refers to the solution of the linearised problem, and is used as starting value for the non-linear problem.

The files **`src/result_cu.txt`** and **`src/result_cv.txt`** contain the final result of the Newton-Raphson method and can be plotted in the Matlab environment to get a graphical result of the solution. The terms cu and cv correspond to the concentrations in oxygen and carbon dioxyde respectively. A tolerance of 10<sup>-10</sup> is allowed in the algorithm.
