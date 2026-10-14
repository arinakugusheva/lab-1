"""
Лабораторная работа № 1. Вариант 13.
Итерационные методы вычисления sqrt(a), cbrt(b), ln(c).
a=31, b=100, c=0.2, eps=1e-5, критерий остановки: по разности (Р), x0 = a/2, доп. требование Д1.
"""
import math

A, B, C = 31, 100, 0.2
EPS = 1e-5
N_MAX = 100


def sqrt_newton(a, eps=EPS, n_max=N_MAX, trace=None):
    """sqrt(a) методом Ньютона (Герон). Остановка по разности |x_new - x| < eps."""
    if a < 0:
        raise ValueError("Подкоренное выражение должно быть неотрицательным")
    if a == 0:
        return 0.0, 0
    x = a / 2                              # x0 = a/2
    for n in range(1, n_max + 1):
        x_new = 0.5 * (x + a / x)
        if trace is not None:
            trace.append((n, x_new, abs(x_new - x)))
        if abs(x_new - x) < eps:
            return x_new, n
        x = x_new
    raise RuntimeError("Точность не достигнута за N_MAX итераций")


def cbrt_newton(b, eps=EPS, n_max=N_MAX, trace=None):
    """cbrt(b) методом Ньютона. Остановка по разности |x_new - x| < eps."""
    if b == 0:
        return 0.0, 0
    sign = -1 if b < 0 else 1
    b = abs(b)
    x = b / 2                              # x0 = b/2 (аналог x0 = a/2)
    for n in range(1, n_max + 1):
        x_new = (2 * x + b / (x * x)) / 3
        if trace is not None:
            trace.append((n, x_new, abs(x_new - x)))
        if abs(x_new - x) < eps:
            return sign * x_new, n
        x = x_new
    raise RuntimeError("Точность не достигнута за N_MAX итераций")


def ln_series(c, eps=EPS, n_max=10_000, trace=None):
    """ln(c) = 2*(y + y^3/3 + y^5/5 + ...), y=(c-1)/(c+1). Остановка по |член| < eps."""
    if c <= 0:
        raise ValueError("Аргумент логарифма должен быть положительным")
    y = (c - 1) / (c + 1)
    y2 = y * y
    power, total = y, 0.0
    for k in range(n_max):
        term = 2 * power / (2 * k + 1)
        if trace is not None:
            trace.append((k + 1, total + term, abs(term)))
        total += term                      # последний (малый) член тоже добавляем
        if abs(term) < eps:
            return total, k + 1
        power *= y2
    raise RuntimeError("Точность не достигнута за n_max итераций")


def correct_digits(x, ref):
    """Число верных десятичных знаков после запятой."""
    err = abs(x - ref)
    return 16 if err == 0 else max(0, int(-math.log10(err)))


def print_trace(title, trace, ref):
    print(f"\n{title}")
    print(f"{'n':>3}{'x_n':>18}{'|x_n - x_(n-1)|':>20}{'верных знаков':>16}")
    for n, x, d in trace:
        print(f"{n:>3}{x:>18.10f}{d:>20.3e}{correct_digits(x, ref):>16}")


if __name__ == "__main__":
    tables = {}
    rows = []
    for name, f, arg, ref in [("sqrt", sqrt_newton, A, math.sqrt(A)),
                              ("cbrt", cbrt_newton, B, B ** (1 / 3)),
                              ("ln", ln_series, C, math.log(C))]:
        tr = []
        val, n = f(arg, trace=tr)
        rows.append((name, arg, val, ref, abs(val - ref), n))
        tables[name] = (tr, ref)

    print(f"{'Функция':<8}{'Аргумент':>9}{'Результат':>16}{'Эталон':>16}"
          f"{'Погрешность':>14}{'Итераций':>10}")
    for name, arg, val, ref, err, n in rows:
        print(f"{name:<8}{arg:>9}{val:>16.10f}{ref:>16.10f}{err:>14.2e}{n:>10}")

    print_trace("Д1: sqrt(31)", *tables["sqrt"])
    print_trace("Д1: cbrt(100)", *tables["cbrt"])
    print_trace("Д1: ln(0.2) (x_n = частичная сумма, шаг = |член ряда|)", *tables["ln"])
