import numpy as np
from math import floor
import random


def two_phase_simplex(tableau):
    original_objective = tableau[-1, :].copy()
    tableau[-1, :] = 0
    artificial_vars = []

    for i in range(tableau.shape[0] - 1):
        if tableau[i, -1] < 0:
            tableau[i, :] = -tableau[i, :]
        if tableau[i, -1] > 0:
            artificial_var = np.zeros((tableau.shape[0], 1))
            artificial_var[i, 0] = 1
            artificial_var[-1, 0] = 1
            tableau = np.hstack((tableau[:, :-1], artificial_var, tableau[:, -1:]))
            artificial_vars.append(tableau.shape[1] - 2)

    tableau[-1, artificial_vars] = 1

    while True:
        entering_col = bland_rule_entering_col(tableau)

        if entering_col is None:
            break

        leaving_row = bland_rule_leaving_row(tableau, entering_col)
        if leaving_row is None:
            raise ValueError("Unbounded problem")

        pivot(tableau, leaving_row, entering_col)

    # Check if any artificial variables are left in the basis
    for i in range(tableau.shape[0] - 1):
        if tableau[i, -1] != 0 and i in artificial_vars:
            raise ValueError("Infeasible problem")

    tableau = tableau[:, [i for i in range(tableau.shape[1]) if i not in artificial_vars]]
    tableau[-1, :] = original_objective

    return tableau

# def is_integer(value):
#     return abs(round(value) - value) < 1e-9

def is_integer(value):
    return np.isclose(round(value), value, rtol=1e-9)

def find_fractional_row(tableau):
    for i in range(tableau.shape[0] - 1):
        integ = np.floor(tableau[i, -1])
        frac = tableau[i, -1] - integ
        if (frac < 0.00000000001 or frac > 1 - 0.00000000001):
            continue
        else:
            return i
    return None


def bland_rule_exiting_row_dual(tableau):
    exiting_row = None
    most_negative_value = 0
    for i in range(tableau.shape[0] - 1):
        if tableau[i, -1] < 0 and (exiting_row is None or tableau[i, -1] < most_negative_value):
            most_negative_value = tableau[i, -1]
            exiting_row = i
    return exiting_row

def bland_rule_entering_col_dual(tableau, exiting_row):
    entering_col = None
    min_ratio = float('inf')
    for i in range(tableau.shape[1] - 2):
        if tableau[exiting_row, i] < 0:
            ratio = -tableau[-1, i] / tableau[exiting_row, i]
            #if ratio < min_ratio and not np.isclose(ratio, 0):
            if ratio < min_ratio:
                min_ratio = ratio
                entering_col = i
    return entering_col


def bland_rule_entering_col(tableau):
    entering_col = None
    most_negative_value = 0
    for i in range(tableau.shape[1] - 1):
        if tableau[-1, i] < 0 and (entering_col is None or tableau[-1, i] < most_negative_value):
            most_negative_value = tableau[-1, i]
            entering_col = i
    return entering_col

def bland_rule_leaving_row(tableau, entering_col):
    leaving_row = None
    min_ratio = float('inf')
    for i in range(tableau.shape[0] - 1):
        if tableau[i, entering_col] > 0:
            ratio = tableau[i, -1] / tableau[i, entering_col]
            #if ratio < min_ratio or (np.isclose(ratio, min_ratio) and i < leaving_row):
            if ratio < min_ratio:
                min_ratio = ratio
                leaving_row = i
    return leaving_row

def pivot_dual(tableau, row, col):
    pivot_value = tableau[row, col]
    tableau[row, :] /= pivot_value

    for i in range(tableau.shape[0]):
        if i != row:
            multiplier = tableau[i, col] / tableau[row, col]
            tableau[i, :] -= multiplier * tableau[row, :]
            
def pivot(tableau, row, col):
    tableau[row, :] /= tableau[row, col]
    for i in range(tableau.shape[0]):
        if i != row:
            tableau[i, :] -= tableau[i, col] * tableau[row, :]

def gomory(filename):
    with open(filename, 'r') as f:
        n, m = map(int, f.readline().split())
        b = list(map(int, f.readline().split()))
        c = list(map(int, f.readline().split()))
        A = [list(map(int, f.readline().split())) for i in range(m)]

    tableau = np.zeros((m + 1, n + m + 1))
    tableau[:-1, :n] = A
    tableau[:-1, n:n + m] = np.eye(m)
    tableau[:-1, -1] = b
    tableau[-1, :n] = -np.array(c)
    tableau[-1, -1] = 0

    if any(value < 0 for value in b) or any(value < 0 for value in c):
        tableau = two_phase_simplex(tableau)

    while True:
        
        entering_col = bland_rule_entering_col(tableau)

        if entering_col is None:
            break

        leaving_row = bland_rule_leaving_row(tableau, entering_col)
        if leaving_row is None:
            raise ValueError("Unbounded problem")

        pivot(tableau, leaving_row, entering_col)

    while True:
        integer_solution = True
        for row in range(tableau.shape[0]):
            #Last column ke values
            if not is_integer(tableau[row, -1]):
                integer_solution = False
                break
        # If all variables have integer values, break out of the loop
        if integer_solution:
            break

        frac_row = find_fractional_row(tableau)
      

        if frac_row is not None:
            #integer_part = np.floor(tableau[frac_row, :-1])
            integer_part = np.floor(tableau[frac_row])
            
            fractional_part = tableau[frac_row] - integer_part
     
            for i in range(len(fractional_part)):
                if fractional_part[i] < 0.00000000001 or fractional_part[i] > 1 - 0.00000000001:
                    fractional_part[i] = 0


            #An array of zeros
            gomory_cut = np.zeros(tableau.shape[1])

            gomory_cut = -fractional_part
        
            # Add a slack variable
            gomory_cut = np.insert(gomory_cut, -1, 1)
            # Add a new column for the slack variable
            new_slack_col = np.zeros((tableau.shape[0], 1))
            tableau = np.hstack((tableau[:, :-1], new_slack_col, tableau[:, -1:]))
            # Add the Gomory cut constraint to the tableau
            tableau = np.vstack((tableau, gomory_cut))
            
            rows = tableau.shape[0]
            temp1 = tableau[rows-1].copy()
            temp2 = tableau[rows-2].copy()
            tableau[rows-1] = temp2
            tableau[rows-2] = temp1
      
            # Update basic variables
            # basic_vars = np.where(np.isclose(tableau[:-1, :-1].sum(axis=0), 1))[0]
            # for row, var in enumerate(basic_vars):
            #     if row != frac_row:
            #         tableau[row, -2] = 0
            # tableau[-2, -2] = 1

            basic_vars = np.where(np.isclose(tableau[:-1, :-1].sum(axis=0), 1))[0]
            for i, var in enumerate(basic_vars):
                if var < tableau.shape[0] - 1 and i != frac_row:
                    tableau[var, -2] = 0
            tableau[-2, -2] = 1

            m += 1  # Update the number of constraints
            
        while True:
   
            exiting_row = bland_rule_exiting_row_dual(tableau)
            if exiting_row is None:
                break

            enter_col = bland_rule_entering_col_dual(tableau, exiting_row)
            if enter_col is None:
                raise ValueError("Unbounded problem")
            pivot_dual(tableau, exiting_row, enter_col)
    

    x = np.zeros(n)
    for j in range(n):
        pivot_rows = np.where(np.isclose(tableau[:, j], 1))[0]
        if len(pivot_rows) == 1 and pivot_rows[0] < m:
            x[j] = tableau[pivot_rows[0], -1]

    x = np.round(x)
    #objective_value = np.dot(c, x)

    return x

if __name__ == '__main__':
# Read input from file
    filename = input("Enter the input file name: ")
    # Call the gomory() function to solve the ILP
    x  = gomory(filename)
    #print("Optimal solution:", x)
    
