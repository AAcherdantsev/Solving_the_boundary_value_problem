"""
funcs[0] = A(x)
funcs[1] = B(x)
funcs[2] = C(x)
Остальные имена переменных сохранены
"""

import sympy as sp
import matplotlib.pyplot as plt
from numpy import linspace

def input_data(file_name): # ввод данных
    with open(file_name) as file:
        for _ in range(5):
            file.readline()  # пропускаем первые 5 строк
        funcs = [sp.simplify(file.readline().split("=")[1]) for _ in range(3)]
        F = [float(sp.simplify(file.readline().split("=")[1])) for _ in range(2)]
        D = [float(sp.simplify(file.readline().split("=")[1])) for _ in range(2)]
        E = [float(sp.simplify(file.readline().split("=")[1])) for _ in range(2)]
        a = float(sp.simplify(file.readline().split("=")[1]))
        b = float(sp.simplify(file.readline().split("=")[1]))
        n = int(file.readline().split("=")[1])
    return funcs, F, D, E, a, b, n


def run_through_method(funcs, F, D, E, a, b, n, is_first_order_of_approximation = True): #  метод прогонки
    h = (b - a) / n   # шаг
    x = [a + i * h for i in range(n + 1)] # сетка   
    matrix = [[0 for __ in range(n + 1)] for _ in range(n + 1)] 
    vector_d = [0 for _ in range(len(matrix))]
    # Дальше пошло применение формул из учебника

    if is_first_order_of_approximation == True:
        matrix[0][0] = F[0]*h - D[0]
        matrix[0][1] = D[0]
        matrix[-1][-1] = D[1]
        matrix[-1][-2] = F[1]*h - D[1]
        for i in range(1, n):
            matrix[i][i - 1] = 2 - funcs[0].subs("x", x[i]) * h
            matrix[i][i] = (-4 + funcs[1].subs("x", x[i])*2*h**2)
            matrix[i][i + 1] = 2 + funcs[0].subs("x", x[i])*h
        vector_d[0] = E[0] * h
        vector_d[-1] = E[1] * h
        for i in range(1, n):
            vector_d[i] = funcs[2].subs("x", x[i])*2*(h**2)

    if is_first_order_of_approximation == False:
        matrix[0][0] = -(-F[0]*h + D[0] + D[0]*(funcs[0].subs("x", a)\
                       - funcs[1].subs("x", a) * h)*(h/2))
        matrix[0][1] = funcs[0].subs("x", a)*D[0]*(h/2) + D[0]
        matrix[-1][-1] = -(-F[1]*h - D[1] + D[1]*(funcs[0].subs("x", b) \
                         + funcs[1].subs("x", b)* h)*(h/2))
        matrix[-1][-2] = funcs[0].subs("x", b)*D[1]*(h/2) - D[1]

        for i in range(1, n):
            matrix[i][i - 1] = 1 - funcs[0].subs("x", x[i])*(h/2)
            matrix[i][i] = -(2 - funcs[1].subs("x", x[i])*(h**2))
            matrix[i][i + 1] = 1 + funcs[0].subs("x", x[i])*(h/2)
        
        vector_d[0] = E[0]*h + funcs[2].subs("x", a)*D[0]*((h**2)/2)
        vector_d[-1] = E[1]*h - funcs[2].subs("x", b)*D[1]*((h**2)/2)
        for i in range(1, n):
            vector_d[i] = funcs[2].subs("x", x[i])*(h**2)

    n = n + 1
    y_solve = [0 for _ in range(n)]
    # прямой ход прогонки 
    m = 1;
    for i in range(1,n):
        m = matrix[i][i - 1]/matrix[i-1][i-1]
        matrix[i][i] = matrix[i][i] - m*matrix[i-1][i] 
        vector_d[i] = vector_d[i] - m*vector_d[i-1] 
    #Обратный ход
    y_solve[n-1] = vector_d[n-1]/matrix[n-1][n-1];
    for i in range(n - 2, -1, -1):
        y_solve[i]=(vector_d[i] - matrix[i][i + 1]* y_solve [i+1]) / matrix[i][i]
    return x, list(map(float, y_solve))


def show_graphics(x, y, x2, y2, my_func, show_dots = True):
    x_for_func = linspace(x[0], x[-1])
    y_for_func = [my_func.subs("x", x_for_func[i]) for i in range(len(x_for_func))]
    plt.plot(x, y, label = "Первая аппроксимация")  # построение графика
    plt.plot(x_for_func, y_for_func,  label = "Решение из таблицы")  # построение графика
    plt.plot(x2, y2, label = "Вторая аппроксимация")  # построение графика
    if show_dots == True:
        for i in range(len(x)):
            plt.scatter(x[i], y[i], s = 10, c = "blue")
        for i in range(len(x2)):
            plt.scatter(x2[i], y2[i], s = 10, c = "red")
    plt.title("График") # заголовок
    plt.xlabel("x")     # ось абсцисс
    plt.ylabel("u(x)")  # ось ординат
    plt.grid()          # включение отображение сетки
    plt.legend()
    plt.show() 


def print_table(x, y):
    for _ in range(32):
        print("-", end = "")
    print()
    print("{:15}".format("        x"), end = "|")
    print("{:15}".format("      u(x)"), end = "|")
    print()
    for _ in range(32):
        print("-", end = "")
    print()
    for i in range(len(x)):
        print("{:15}".format(round(x[i], 8)), end = "|")
        print("{:15}".format(round(y[i], 8)), end = "|")
        print()
        for _ in range(32):
            print("-", end = "")
        print()
    print()
    print()


def do_fourth_task(funcs, F, D, E, a, b, my_func):
    n_for_test = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 100, 200]
    h_for_test = [(b - a)/n_i for n_i in n_for_test]
    list_for_first_approximation = []
    list_for_second_approximation = []
    for i in range(len(n_for_test)):
        x1, y1 = run_through_method(funcs, F, D, E, a, b, n_for_test[i], True)
        x2, y2 = run_through_method(funcs, F, D, E, a, b, n_for_test[i], False)
        y_correct = [my_func.subs("x", x1[i]) for i in range(len(x1))]
        list_for_first_approximation.append(max([abs(y1[i] - y_correct[i]) for i in range(len(x1))]))
        list_for_second_approximation.append(max([abs(y2[i] - y_correct[i]) for i in range(len(x1))]))

    plt.plot(h_for_test, list_for_first_approximation, label = "Для первого порядка аппроксимации", color = "blue")  # построение графика
    plt.title("График норм") # заголовок
    plt.xlabel("h")     # ось абсцисс
    plt.ylabel("y")     # ось ординат
    plt.grid()          # включение отображение сетки
    plt.legend()
    plt.show()    
    #########################
    plt.plot(h_for_test, list_for_second_approximation,  label = "Для второго порядка аппроксимации", color = "green")  # построение графика
    plt.title("График норм") # заголовок
    plt.xlabel("h")     # ось абсцисс
    plt.ylabel("y")     # ось ординат
    plt.grid()          # включение отображение сетки
    plt.legend()
    plt.show() 
   

if __name__ == "__main__":
    file_name = "input.txt"
    my_func = sp.simplify("(x / 2) + (2*cos(2*x) - sin(2*x)) / 20")   # та функция из таблицы (по вариантам которая)
    funcs, F, D, E, a, b, n = input_data(file_name)
    x_first_app, y_first_app = run_through_method(funcs, F, D, E, a, b, n, True)
    x_second_app, y_second_app = run_through_method(funcs, F, D, E, a, b, n, False)
    # print("Первый порядок аппроксимации")
    # print_table(x_first_app, y_first_app)
    # print("Второй порядок аппроксимации")
    # print_table(x_second_app, y_second_app)
    show_graphics(x_first_app, y_first_app,x_second_app, y_second_app, my_func)  
    do_fourth_task(funcs, F, D, E, a, b, my_func)  # Выполнение 4-го пункта