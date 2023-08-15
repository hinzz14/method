import tkinter as tk
import math 
from scipy.interpolate import lagrange
from scipy.interpolate import barycentric_interpolate
from sympy import Symbol, lambdify, sympify, symbols, diff
from scipy.integrate import quad

def trapezoidal_rule(f, a, b, n):
    h = (b - a) / n
    x = a
    sum_result = f(a)
    for i in range(1, n):
        x += h
        sum_result += 2 * f(x)
    sum_result += f(b)
    return (h / 2) * sum_result

def simpson_rule(f, a, b, n):
    h = (b - a) / n
    x = a
    sum_result = f(a)
    for i in range(1, n):
        x += h
        if i % 2 == 0:
            sum_result += 2 * f(x)
        else:
            sum_result += 4 * f(x)
    sum_result += f(b)
    return (h / 3) * sum_result

def bisection_method(a, b, tolerance):
    # Tạo biến ký hiệu x
    x = symbols('x')

    # Chuyển đổi biểu thức từ chuỗi thành biểu thức toán học
    expr = sympify(entry_function.get())

    # Kiểm tra hợp lệ của biểu thức
    if expr is None:
        return None

    # Kiểm tra tính chất đầu mút của khoảng
    if expr.subs(x, a) * expr.subs(x, b) >= 0:
        return None

    # Áp dụng phương pháp chia đôi
    while (b - a) >= tolerance:
        c = (a + b) / 2
        if expr.subs(x, c) == 0:
            return c
        if expr.subs(x, c) * expr.subs(x, a) < 0:
            b = c
        else:
            a = c

    return (a + b) / 2

def calculate_dathucnoisuy():
    x_values = [float(x) for x in x_dathuc_entry.get().split()]
    y_values = [float(y) for y in y_dathuc_entry.get().split()]
    interp_value = float(noisuy_dathuc_entry.get())
    
    # Tính đa thức nội suy bằng phương pháp Lagrange
    lagrange_result = lagrange(x_values, y_values)
    result = lagrange_result(interp_value)
    lagrange_label.config(text=f"Lagrange: {result}")
    lagrange_label.pack()  # Hiển thị dòng kết quả Lagrange

    #tính đa thức nội suy bằng phương pháp newton
    newton_result = barycentric_interpolate(x_values, y_values, interp_value)
    newton_label.config(text=f"Newton: {newton_result}")
    newton_label.pack()

    x_dathuc_label.pack_forget()
    x_dathuc_entry.pack_forget()
    y_dathuc_label.pack_forget()
    y_dathuc_entry.pack_forget()
    noisuy_dathuc_label.pack_forget()
    noisuy_dathuc_entry.pack_forget()
    button_tinhdathuc.pack_forget()
    button_back.pack_forget()
    button_backdathuc.pack()

def calculate_tinhgandungtichphan():
    function_str = function_entry.get()
    lower_limit = float(lower_entry.get())
    upper_limit = float(upper_entry.get())
    num_segments = int(segments_entry.get())
    
    x = Symbol('x')
    f = lambdify(x, function_str)
    
    trapezoidal_result = trapezoidal_rule(f, lower_limit, upper_limit, num_segments)
    simpson_result = simpson_rule(f, lower_limit, upper_limit, num_segments)
    
    result_gandungtichphan_label.config(text=f"Công thức hình thang: {trapezoidal_result}\n"
                                            f"Công thức Simpson: {simpson_result}")
def calculate_chiadoi():
    a = float(entry_a.get())
    b = float(entry_b.get())
    tolerance = float(entry_tolerance.get())

    result = bisection_method(a, b, tolerance)
    if result is None:
        result_chiadoi_label.configure(text="Không tìm thấy nghiệm trong khoảng này.")
    else:
        result_chiadoi_label.configure(text="Nghiệm gần đúng: {:.6f}".format(result))
    result_chiadoi_label.pack()
    button_backphuongtrinh.pack()

def calculate_lapdon():
    # Lấy giá trị từ các ô nhập liệu
    function_lapdon_text = function_lapdon_entry.get()
    initial_lapdon_value = float(initial_lapdon_entry.get())
    error_lapdon = float(error_lapdon_entry.get())
    
    # Phương pháp lặp đơn
    def g(x):
        return eval(function_lapdon_text)
    
    # Khởi tạo giá trị ban đầu
    x0 = initial_lapdon_value
    
    # Áp dụng phương pháp lặp đơn
    while True:
        x1 = g(x0)
        
        # Kiểm tra điều kiện dừng
        if abs(x1 - x0) < error_lapdon:
            break
        
        x0 = x1
    
    # Hiển thị kết quả
    result_lapdon_label.config(text=f"Kết quả: {x1}")

def secant_method(a, b, tolerance_daycung):
    # Tạo biến ký hiệu x
    x = symbols('x')

    # Chuyển đổi biểu thức từ chuỗi thành biểu thức toán học
    expr_daycung = sympify(entry_daycung_function.get())

    # Kiểm tra hợp lệ của biểu thức
    if expr_daycung is None:
        return None

    # Tính đạo hàm của biểu thức
    f_prime = diff(expr_daycung, x)

    # Áp dụng phương pháp dây cung
    while abs(b - a) >= tolerance_daycung:
        f_a = expr_daycung.subs(x, a)
        f_b = expr_daycung.subs(x, b)

        if f_a == f_b:
            return None

        a, b = b, b - (f_b * (b - a)) / (f_b - f_a)

    return b

def calculate_daycung():
    a = float(entry_daycung_a.get())
    b = float(entry_daycung_b.get())
    tolerance_daycung = float(entry_tolerance_daycung.get())

    result_daycung = secant_method(a, b, tolerance_daycung)
    if result_daycung is None:
        result_daycung_label.configure(text="Không tìm thấy nghiệm trong khoảng này.")
    else:
        result_daycung_label.configure(text="Nghiệm gần đúng: {:.6f}".format(result_daycung))

def jacobi_method(matrix_jacobi_hpt, b_jacobi_hpt, tolerance_jacobi_hpt):
    n_jacobi_hpt = len(matrix_jacobi_hpt)
    x_jacobi_hpt = [0] * n_jacobi_hpt

    iterations_jacobi_hpt = 0
    while True:
        next_x_jacobi_hpt = x_jacobi_hpt.copy()
        for i in range(n_jacobi_hpt):
            sum_jacobi_hpt = 0
            for j in range(n_jacobi_hpt):
                if j != i:
                    sum_jacobi_hpt += matrix_jacobi_hpt[i][j] * x_jacobi_hpt[j]
            next_x_jacobi_hpt[i] = (b_jacobi_hpt[i] - sum_jacobi_hpt) / matrix_jacobi_hpt[i][i]

        # Kiểm tra điều kiện dừng
        error_jacobi_hpt = max(abs(next_x_jacobi_hpt[i] - x_jacobi_hpt[i]) for i in range(n_jacobi_hpt))
        if error_jacobi_hpt < tolerance_jacobi_hpt:
            return next_x_jacobi_hpt, iterations_jacobi_hpt

        x_jacobi_hpt = next_x_jacobi_hpt
        iterations_jacobi_hpt += 1

def gauss_seidel_method(matrix_seidel_hpt, b_seidel_hpt, tolerance_seidel_hpt):
    n_seidel_hpt = len(matrix_seidel_hpt)
    x_seidel_hpt = [0] * n_seidel_hpt

    iterations_seidel_hpt = 0
    while True:
        next_x_seidel_hpt = x_seidel_hpt.copy()
        for i in range(n_seidel_hpt):
            sum_seidel_hpt = 0
            for j in range(n_seidel_hpt):
                if j != i:
                    sum_seidel_hpt += matrix_seidel_hpt[i][j] * next_x_seidel_hpt[j]
            next_x_seidel_hpt[i] = (b_seidel_hpt[i] - sum_seidel_hpt) / matrix_seidel_hpt[i][i]

        # Kiểm tra điều kiện dừng
        error_seidel_hpt = max(abs(next_x_seidel_hpt[i] - x_seidel_hpt[i]) for i in range(n_seidel_hpt))
        if error_seidel_hpt < tolerance_seidel_hpt:
            return next_x_seidel_hpt, iterations_seidel_hpt

        x_seidel_hpt = next_x_seidel_hpt
        iterations_seidel_hpt += 1

def solve_equations_hpt():
    rows_hpt = int(entry_rows_hpt.get())
    cols_hpt = int(entry_cols_hpt.get())

    matrix_hpt = []
    b_hpt = []

    for i in range(rows_hpt):
        row_hpt = []

        for j in range(cols_hpt):
            value_hpt = float(entry_matrix_hpt[i][j].get())
            row_hpt.append(value_hpt)

        matrix_hpt.append(row_hpt)

        b_value_hpt = float(entry_b_hpt[i].get())
        b_hpt.append(b_value_hpt)

    tolerance_hpt = float(entry_tolerance_hpt.get())

    result_jacobi_hpt, iterations_jacobi_hpt = jacobi_method(matrix_hpt, b_hpt, tolerance_hpt)
    result_seidel_hpt, iterations_seidel_hpt = gauss_seidel_method(matrix_hpt, b_hpt, tolerance_hpt)

    if result_jacobi_hpt is None:
        result_label_hpt.configure(text="Phương pháp Jacobi không hội tụ sau số lần lặp cho trước.")
    else:
        result_text_hpt = "Phương pháp Jacobi:\n"
        for i in range(len(result_jacobi_hpt)):
            result_text_hpt += "x{} = {:.6f}\n".format(i+1, result_jacobi_hpt[i])
        result_text_hpt += "Số lần lặp: {}".format(iterations_jacobi_hpt)
        result_text_hpt += "\n\n"

    if result_seidel_hpt is None:
        result_label_hpt.configure(text="Phương pháp Gauss-Seidel không hội tụ sau số lần lặp cho trước.")
    else:
        result_text_hpt += "Phương pháp Gauss-Seidel:\n"
        for i in range(len(result_seidel_hpt)):
            result_text_hpt += "x{} = {:.6f}\n".format(i+1, result_seidel_hpt[i])
        result_text_hpt += "Số lần lặp: {}".format(iterations_seidel_hpt)

    result_label_hpt.configure(text=result_text_hpt)

def create_matrix_entries_hpt():
    rows_hpt = int(entry_rows_hpt.get())
    cols_hpt = int(entry_cols_hpt.get())

    for i in range(rows_hpt):
        row_hpt = []

        for j in range(cols_hpt):
            entry_hpt = tk.Entry(matrix_frame_hpt)
            entry_hpt.grid(row=i, column=j, padx=5, pady=5)
            row_hpt.append(entry_hpt)

        entry_matrix_hpt.append(row_hpt)

        label_hpt = tk.Label(matrix_frame_hpt, text="b{}:".format(i+1))
        label_hpt.grid(row=i, column=cols_hpt, padx=5, pady=5)

        entry_b_i_hpt = tk.Entry(matrix_frame_hpt)
        entry_b_i_hpt.grid(row=i, column=cols_hpt+1, padx=5, pady=5)
        entry_b_hpt.append(entry_b_i_hpt)

def gauss_elimination_khugauss(matrix_khugauss):
    n_khugauss = len(matrix_khugauss)
    
    for i_khugauss in range(n_khugauss):
        if matrix_khugauss[i_khugauss][i_khugauss] == 0:
            return None
            
        for j_khugauss in range(i_khugauss+1, n_khugauss):
            ratio_khugauss = matrix_khugauss[j_khugauss][i_khugauss] / matrix_khugauss[i_khugauss][i_khugauss]
            
            for k_khugauss in range(n_khugauss+1):
                matrix_khugauss[j_khugauss][k_khugauss] -= ratio_khugauss * matrix_khugauss[i_khugauss][k_khugauss]
    
    result_khugauss = [0] * n_khugauss
    
    for i_khugauss in range(n_khugauss-1, -1, -1):
        result_khugauss[i_khugauss] = matrix_khugauss[i_khugauss][n_khugauss] / matrix_khugauss[i_khugauss][i_khugauss]
        
        for j_khugauss in range(i_khugauss-1, -1, -1):
            matrix_khugauss[j_khugauss][n_khugauss] -= matrix_khugauss[j_khugauss][i_khugauss] * result_khugauss[i_khugauss]
    
    return result_khugauss

def solve_equations_khugauss():
    rows_khugauss = int(entry_rows_khugauss.get())
    cols_khugauss = int(entry_cols_khugauss.get())

    matrix_khugauss = []
    
    for i_khugauss in range(rows_khugauss):
        row_khugauss = []
        
        for j_khugauss in range(cols_khugauss):
            value_khugauss = float(entry_matrix_khugauss[i_khugauss][j_khugauss].get())
            row_khugauss.append(value_khugauss)
        
        matrix_khugauss.append(row_khugauss)
    
    result_khugauss = gauss_elimination_khugauss(matrix_khugauss)
    
    if result_khugauss is None:
        result_label_khugauss.configure(text="Hệ phương trình không có nghiệm duy nhất.")
    else:
        result_text_khugauss = "Nghiệm:\n"
        for i_khugauss in range(len(result_khugauss)):
            result_text_khugauss += "x{} = {:.6f}\n".format(i_khugauss+1, result_khugauss[i_khugauss])
        result_label_khugauss.configure(text=result_text_khugauss)

def create_matrix_entries_khugauss():
    rows_khugauss = int(entry_rows_khugauss.get())
    cols_khugauss = int(entry_cols_khugauss.get())
    
    for i_khugauss in range(rows_khugauss):
        row_khugauss = []
        
        for j_khugauss in range(cols_khugauss):
            entry_khugauss = tk.Entry(matrix_frame_khugauss)
            entry_khugauss.grid(row=i_khugauss, column=j_khugauss, padx=5, pady=5)
            row_khugauss.append(entry_khugauss)
        
        entry_matrix_khugauss.append(row_khugauss)

def show_buttons_dathuc():
    x_dathuc_label.pack()
    x_dathuc_entry.pack()
    y_dathuc_label.pack()
    y_dathuc_entry.pack()
    noisuy_dathuc_label.pack()
    noisuy_dathuc_entry.pack()
    button_giaiphuongtrinh.pack_forget()
    button_dathuc.pack_forget()
    button_tinhdathuc.pack()
    button_back.pack()
    button_backdathuc.pack_forget()
    lagrange_label.pack_forget()
    newton_label.pack_forget()
    button_tinhgandungtichphan.pack_forget()
    button_giaihephuongtrinh.pack_forget()

def show_buttons_gandungtichphan():
    button_dathuc.pack_forget()
    button_tinhgandungtichphan.pack_forget()
    button_giaiphuongtrinh.pack_forget()
    function_label.pack()
    function_entry.pack()
    lower_label.pack()
    lower_entry.pack()
    upper_label.pack()
    upper_entry.pack()
    segments_label.pack()
    segments_entry.pack()
    calculate_button.pack()
    result_gandungtichphan_label.pack()
    button_back.pack()
    button_giaihephuongtrinh.pack_forget()

def show_buttons_giaiphuongtrinh():
    function_lapdon_label.pack_forget()
    function_lapdon_entry.pack_forget()
    initial_lapdon_label.pack_forget()
    initial_lapdon_entry.pack_forget()
    error_lapdon_label.pack_forget()
    error_lapdon_entry.pack_forget()
    calculate_lapdon_button.pack_forget()
    result_lapdon_label.pack_forget()
    button_dathuc.pack_forget()
    button_tinhgandungtichphan.pack_forget()
    button_giaiphuongtrinh.pack_forget()
    button_chiadoi.pack()
    button_backphuongtrinh.pack_forget()
    label_function.pack_forget()
    entry_function.pack_forget()
    label_a.pack_forget()
    entry_a.pack_forget()
    label_b.pack_forget()
    entry_b.pack_forget()
    label_tolerance.pack_forget()
    entry_tolerance.pack_forget()
    calculate_chiadoi_button.pack_forget()
    result_chiadoi_label.pack_forget()
    button_lapdon.pack()
    button_daycung.pack()
    button_back.pack()
    label_daycung_function.pack_forget()
    entry_daycung_function.pack_forget()
    label_daycung_a.pack_forget()
    entry_daycung_a.pack_forget()
    label_daycung_b.pack_forget()
    entry_daycung_b.pack_forget()
    label_tolerance_daycung.pack_forget()
    entry_tolerance_daycung.pack_forget()
    calculate_daycung_button.pack_forget()
    result_daycung_label.pack_forget()
    button_giaihephuongtrinh.pack_forget()

def show_buttons_chiadoi():
    button_dathuc.pack_forget()
    button_tinhgandungtichphan.pack_forget()
    button_giaiphuongtrinh.pack_forget()
    label_function.pack()
    entry_function.pack()
    label_a.pack()
    entry_a.pack()
    label_b.pack()
    entry_b.pack()
    label_tolerance.pack()
    entry_tolerance.pack()
    calculate_chiadoi_button.pack()
    button_chiadoi.pack_forget()
    result_chiadoi_label.pack_forget()
    button_backphuongtrinh.pack()
    button_back.pack_forget()
    button_lapdon.pack_forget()
    button_daycung.pack_forget()

def show_buttons_lapdon():
    button_chiadoi.pack_forget()
    button_dathuc.pack_forget()
    button_tinhgandungtichphan.pack_forget()
    button_giaiphuongtrinh.pack_forget()
    function_lapdon_label.pack()
    function_lapdon_entry.pack()
    initial_lapdon_label.pack()
    initial_lapdon_entry.pack()
    error_lapdon_label.pack()
    error_lapdon_entry.pack()
    calculate_lapdon_button.pack()
    result_lapdon_label.pack()
    button_backphuongtrinh.pack()
    button_back.pack_forget()
    button_lapdon.pack_forget()
    button_daycung.pack_forget()

def show_buttons_daycung():
    button_chiadoi.pack_forget()
    button_dathuc.pack_forget()
    button_tinhgandungtichphan.pack_forget()
    button_giaiphuongtrinh.pack_forget()
    label_daycung_function.pack()
    entry_daycung_function.pack()
    label_daycung_a.pack()
    entry_daycung_a.pack()
    label_daycung_b.pack()
    entry_daycung_b.pack()
    label_tolerance_daycung.pack()
    entry_tolerance_daycung.pack()
    calculate_daycung_button.pack()
    result_daycung_label.pack()
    button_backphuongtrinh.pack()
    button_back.pack_forget()
    button_lapdon.pack_forget()
    button_daycung.pack_forget()

def show_buttons_giaihephuongtrinh():
    button_dathuc.pack_forget()
    button_tinhgandungtichphan.pack_forget()
    button_giaiphuongtrinh.pack_forget()
    button_giaihephuongtrinh.pack_forget()
    button_tinhgandungtichphan.pack_forget()
    button_giaiphuongtrinh.pack_forget()
    button_chiadoi.pack_forget()
    button_backphuongtrinh.pack_forget()
    button_khugauss.pack()
    button_ptlapdon.pack()
    button_back.pack()

def show_buttons_backhephuongtrinh():
    button_dathuc.pack_forget()
    button_tinhgandungtichphan.pack_forget()
    button_giaiphuongtrinh.pack_forget()
    button_giaihephuongtrinh.pack_forget()
    button_tinhgandungtichphan.pack_forget()
    button_giaiphuongtrinh.pack_forget()
    button_chiadoi.pack_forget()
    button_backphuongtrinh.pack_forget()
    button_khugauss.pack()
    button_ptlapdon.pack()
    button_back.pack()
    label_rows_hpt.pack_forget()
    entry_rows_hpt.pack_forget()
    label_cols_hpt.pack_forget()
    entry_cols_hpt.pack_forget()
    matrix_frame_hpt.pack_forget()
    create_button_hpt.pack_forget()
    label_tolerance_hpt.pack_forget()
    entry_tolerance_hpt.pack_forget()
    solve_button_hpt.pack_forget()
    result_label_hpt.pack_forget()
    button_backhephuongtrinh.pack_forget()
    label_rows_khugauss.pack_forget()
    entry_rows_khugauss.pack_forget()
    label_cols_khugauss.pack_forget()
    entry_cols_khugauss.pack_forget()
    matrix_frame_khugauss.pack_forget()
    create_button_khugauss.pack_forget()
    solve_button_khugauss.pack_forget()
    result_label_khugauss.pack_forget()

def show_buttons_khugauss():
    label_rows_khugauss.pack()
    entry_rows_khugauss.pack()
    label_cols_khugauss.pack()
    entry_cols_khugauss.pack()
    matrix_frame_khugauss.pack()
    entry_matrix_khugauss = []
    create_button_khugauss.pack()
    solve_button_khugauss.pack()
    result_label_khugauss.pack()
    button_dathuc.pack_forget()
    button_tinhgandungtichphan.pack_forget()
    button_giaiphuongtrinh.pack_forget()
    button_giaihephuongtrinh.pack_forget()
    button_backhephuongtrinh.pack()
    button_khugauss.pack_forget()
    button_ptlapdon.pack_forget()
    button_back.pack_forget()

def show_buttons_ptlapdon():
    label_rows_hpt.pack()
    entry_rows_hpt.pack()
    label_cols_hpt.pack()
    entry_cols_hpt.pack()
    matrix_frame_hpt.pack()
    entry_matrix_hpt = []
    entry_b_hpt = []
    create_button_hpt.pack()
    label_tolerance_hpt.pack()
    entry_tolerance_hpt.pack()
    solve_button_hpt.pack()
    result_label_hpt.pack()
    button_dathuc.pack_forget()
    button_tinhgandungtichphan.pack_forget()
    button_giaiphuongtrinh.pack_forget()
    button_giaihephuongtrinh.pack_forget()
    button_backhephuongtrinh.pack()
    button_khugauss.pack_forget()
    button_ptlapdon.pack_forget()
    button_back.pack_forget()

def show_initial_buttons():
    button_dathuc.pack()
    button_tinhgandungtichphan.pack()
    button_giaiphuongtrinh.pack()
    button_giaihephuongtrinh.pack()
    button_tinhdathuc.pack_forget()
    button_back.pack_forget()
    x_dathuc_label.pack_forget()
    x_dathuc_entry.pack_forget()
    y_dathuc_label.pack_forget()
    y_dathuc_entry.pack_forget()
    noisuy_dathuc_label.pack_forget()
    noisuy_dathuc_entry.pack_forget()
    lagrange_label.pack_forget()
    function_label.pack_forget()
    function_entry.pack_forget()
    lower_label.pack_forget()
    lower_entry.pack_forget()
    upper_label.pack_forget()
    upper_entry.pack_forget()
    segments_label.pack_forget()
    segments_entry.pack_forget()
    calculate_chiadoi_button.pack_forget()
    result_gandungtichphan_label.pack_forget()
    button_chiadoi.pack_forget()
    label_function.pack_forget()
    entry_function.pack_forget()
    label_a.pack_forget()
    entry_a.pack_forget()
    label_b.pack_forget()
    entry_b.pack_forget()
    label_tolerance.pack_forget()
    entry_tolerance.pack_forget()
    calculate_button.pack_forget()
    button_lapdon.pack_forget()
    button_daycung.pack_forget()
    label_rows_hpt.pack_forget()
    entry_rows_hpt.pack_forget()
    label_cols_hpt.pack_forget()
    entry_cols_hpt.pack_forget()
    matrix_frame_hpt.pack_forget()
    create_button_hpt.pack_forget()
    label_tolerance_hpt.pack_forget()
    entry_tolerance_hpt.pack_forget()
    solve_button_hpt.pack_forget()
    result_label_hpt.pack_forget()
    button_khugauss.pack_forget()
    button_ptlapdon.pack_forget()

def go_back():
    show_initial_buttons()

root = tk.Tk()
root.title("Phạm Duy Lợi - 222601120")  
root.geometry("600x500")

button_dathuc = tk.Button(root, text="Đa thức nội suy", command=show_buttons_dathuc)
button_tinhdathuc = tk.Button(root, text="Tính", command=calculate_dathucnoisuy)
button_back = tk.Button(root, text="Trở về Menu",command=go_back)
x_dathuc_label = tk.Label(root, text="Nhập giá trị x (cách nhau bởi dấu cách):")
x_dathuc_entry = tk.Entry(root)
y_dathuc_label = tk.Label(root, text="Nhập giá trị y (cách nhau bởi dấu cách):")
y_dathuc_entry = tk.Entry(root)
noisuy_dathuc_label = tk.Label(root, text="Nhập giá trị cần nội suy:")
noisuy_dathuc_entry = tk.Entry(root)
lagrange_label=tk.Label(root, text ="Lagrange")
newton_label = tk.Label(root, text = "Newton")
button_backdathuc=tk.Button(root, text="back",command=show_buttons_dathuc)
button_tinhgandungtichphan = tk.Button(root, text = "Tính gần đúng tích phân",command=show_buttons_gandungtichphan)
function_label = tk.Label(root, text="Nhập hàm f:")
function_entry = tk.Entry(root)
lower_label = tk.Label(root, text="Giới hạn dưới:")
lower_entry = tk.Entry(root)
upper_label = tk.Label(root, text="Giới hạn trên:")
upper_entry = tk.Entry(root)
segments_label = tk.Label(root, text="số đoạn:")
segments_entry = tk.Entry(root)
calculate_button = tk.Button(root, text="Tính", command=calculate_tinhgandungtichphan)
result_gandungtichphan_label = tk.Label(root, text="")
button_giaiphuongtrinh = tk.Button(root, text = "Giải phương trình",command = show_buttons_giaiphuongtrinh)
button_giaihephuongtrinh = tk.Button(root, text = "Giải hệ phương trình tuyến tính",command = show_buttons_giaihephuongtrinh)
button_khugauss = tk.Button(root, text = "Phương pháp khử Gauss",command = show_buttons_khugauss)
button_ptlapdon = tk.Button(root,text = "Phương pháp lặp", command = show_buttons_ptlapdon)
button_backhephuongtrinh = tk.Button(root, text = "Back",command = show_buttons_backhephuongtrinh)
button_chiadoi = tk.Button(root, text = "Phương pháp chia đôi",command = show_buttons_chiadoi)
button_lapdon = tk.Button(root, text = "Phương pháp lặp đơn",command =  show_buttons_lapdon)
button_daycung = tk.Button(root, text = "Phương pháp dây cung",command = show_buttons_daycung)
label_function = tk.Label(root, text="Hàm số:")
entry_function = tk.Entry(root)
label_a = tk.Label(root, text="Giá trị a:")
entry_a = tk.Entry(root)
label_b = tk.Label(root, text="Giá trị b:")
entry_b = tk.Entry(root)
label_tolerance = tk.Label(root, text="Sai số:")
entry_tolerance = tk.Entry(root)
calculate_chiadoi_button = tk.Button(root, text="Tính toán", command=calculate_chiadoi)
result_chiadoi_label = tk.Label(root, text="")
button_backphuongtrinh = tk.Button(root, text = "back",command = show_buttons_giaiphuongtrinh)
function_lapdon_label = tk.Label(root, text="Hàm f(x) đã biến đổi:")
function_lapdon_entry = tk.Entry(root)
initial_lapdon_label = tk.Label(root, text="Giá trị ban đầu:")
initial_lapdon_entry = tk.Entry(root)
error_lapdon_label = tk.Label(root, text="Sai số:")
error_lapdon_entry = tk.Entry(root)
calculate_lapdon_button = tk.Button(root, text="Tính", command=calculate_lapdon)
result_lapdon_label = tk.Label(root)
label_daycung_function = tk.Label(root, text="Hàm số:")
entry_daycung_function = tk.Entry(root)
label_daycung_a = tk.Label(root, text="Giá trị dưới:")
entry_daycung_a = tk.Entry(root)
label_daycung_b = tk.Label(root, text="Giá trị trên:")
entry_daycung_b = tk.Entry(root)
label_tolerance_daycung = tk.Label(root, text="Sai số:")
entry_tolerance_daycung = tk.Entry(root)
calculate_daycung_button = tk.Button(root, text="Tính toán", command=calculate_daycung)
result_daycung_label = tk.Label(root, text="")
label_rows_hpt = tk.Label(root, text="Số hàng:")
entry_rows_hpt = tk.Entry(root)
label_cols_hpt = tk.Label(root, text="Số cột:")
entry_cols_hpt = tk.Entry(root)
matrix_frame_hpt = tk.Frame(root)
entry_matrix_hpt = []
entry_b_hpt = []
create_button_hpt = tk.Button(root, text="Tạo ma trận", command=create_matrix_entries_hpt)
label_tolerance_hpt = tk.Label(root, text="Sai số :")
entry_tolerance_hpt = tk.Entry(root)
solve_button_hpt = tk.Button(root, text="Giải hệ phương trình", command=solve_equations_hpt)
result_label_hpt = tk.Label(root, text="")
label_rows_khugauss = tk.Label(root, text="Số hàng:")
entry_rows_khugauss = tk.Entry(root)
label_cols_khugauss = tk.Label(root, text="Số cột:")
entry_cols_khugauss = tk.Entry(root)
matrix_frame_khugauss = tk.Frame(root)
entry_matrix_khugauss = []
create_button_khugauss = tk.Button(root, text="Tạo ma trận", command=create_matrix_entries_khugauss)
solve_button_khugauss = tk.Button(root, text="Giải hệ phương trình", command=solve_equations_khugauss)
result_label_khugauss = tk.Label(root, text="")

show_initial_buttons()

root.mainloop()