# Simplex Method with Big M - LLP Solver

## Introduction
This Python script implements the Simplex method with the Big M technique to solve Linear Programming Problems (LLP). It takes input in the form of linear inequalities and equations and determines the optimal solution. The code can handle various conditions such as unbounded solutions, infeasible solutions, and multiple optimal solutions.

## Usage
1. **Objective Function:** Enter the coefficients of the objective function in the format (e.g., `3x1 + 4x2 - 3x3`).
2. **Constraint Equations:** Enter each constraint equation one by one. Use the format `ax1 + bx2 - cx3 >= d`, `ax1 + bx2 - cx3 == d`, or `ax1 + bx2 - cx3 <= d`. You can use spaces for clarity (e.g., `2x1 + 3x2 - 4x3 >= 5`).

## Code Explanation
The script is divided into several functions to enhance readability and modularity:

### 1. Variable Assignment Functions
- `assign_variables(n)`: Assign variables x1, x2, ..., xn.
- `assign_slag_variables(n)`: Assign slack variables s1, s2, ..., sn.
- `assign_artificial_variables(n)`: Assign artificial variables a1, a2, ..., an.

### 2. Matrix and Array Operations Functions
- `create_identity(n)`: Create an identity matrix of size n.
- `min_positive_index(arr)`: Find the index of the minimum positive element in an array.
- `merge_matrix(A, B)`: Merge two matrices element-wise.
- `find_ratio(A, B)`: Find the ratio of elements in arrays A and B.
- `find_sigmaAB(A, B)`: Calculate the dot product of arrays A and B.
- `get_col(A, n)`: Extract a column from a matrix.
- `sub_arr(A, B)`: Subtract array B from array A.

### 3. Input Parsing Functions
- `parse_input()`: Parse user input for LLP, including the objective function and constraint equations.
- `parse_equation(eqn, num_vars)`: Parse a single equation and return its coefficients.

### 4. Simplex Algorithm Functions
- `all_greater_equal_zero(z)`: Check if all elements in the array are greater than or equal to zero.
- `find_cefficient(var, c, no_of_variables)`: Find the coefficient of a variable in the objective function.
- `print_table(...)`: Display the simplex tableau.

### 5. Main Solver Functions
- `simplex(c, val, b, no_of_variables, bv=None, bv_with_var=None, c_with_var=None)`: Implement the simplex algorithm.
- `big_M(c, val, b, no_of_variables, extra_coefficient_matrix)`: Integrate Big M into the simplex algorithm.

## Running the Code
1. Ensure you have Python installed on your machine.
2. Copy and paste the script into a Python environment or save it as a `.py` file.
3. Run the script.
4. Follow the prompts to input the LLP problem details.

## Example
```bash
python llp_solver.py
```

## Notes
- Ensure all input coefficients are integers or decimals.
- The script uses the Big M method to handle artificial variables and ensure feasibility.
- The output includes the optimal solution, variable values, and information about the solution status.

Feel free to modify the script or integrate it into your project as needed.
