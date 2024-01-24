# Code made for SageMath version 9.3
# Initialize all the variables
r = var('r')
T0, T1, T2 = var('T0 T1 T2')
T = [T0, T1, T2]
r_one = [1, r]

f0, f1, f2, f3, f4, f5 = var('f0 f1 f2 f3 f4 f5')
f = [f0, f1, f2, f3, f4, f5]


# Returns f_n(x)
def n_angle(n, x):
    m = n // 2
    result = x^n
    for k in range(1, m + 1):
        p = binomial(n, 2 * k) * (x*x - 1)^k * x^(n - 2*k)
        result += p
    return result


# Returns \cos(\widehat{s1 s2 s3})
def get_angle(s1, s2, s3):
    return ((s1 + s3)^2 - (s1 + s2)^2 - (s2 + s3)^2) / (-2 * (s1 + s2) * (s2 + s3))


# Calculates T_{step} \in a + b
def calculate_polynomial(a, b, step):
    g = T[step]^2 - 2 * a * b * T[step] + a*a + b*b - 1
    return g


# Sorts the list_radii lexicographically
def sort_radii(list_radii):
    to_sort = [(str(item), item) for item in list_radii]
    to_sort.sort()
    return [item for str_item, item in to_sort]


# Returns all the angles in the around-corona list_radii
def get_angles(list_radii, around):
    n = len(list_radii)
    result = {}
    for i in range(len(list_radii)):
        a = [list_radii[i], list_radii[(i + 1) % n]]
        a = sort_radii(a)
        a.insert(1, around)
        lexi_angle = tuple(a)
        if lexi_angle in result:
            result[lexi_angle] += 1
        else:
            result[lexi_angle] = 1
    return result


# Returns the total angle by substituting r=sol into the angle formulas
def get_total_angle(angles, sol):
    num_sol = sol.numerical_approx()
    # print(sol, num_sol)
    cur_angle = 0
    for angle in angles:
        formula = get_angle(*angle)
        cos_f = formula.subs({r: sol})
        a = arccos(cos_f)
        cur_angle += a * angles[angle]
    num_approx = cur_angle.numerical_approx()
    print(cur_angle, num_approx)
    return num_approx


# Calculates the equation for the x-coordinate of the sum of all the elements in list_vars
def add_polynomials(list_vars):
    # No T's means the equation for the x-coordinate should be f0 == 1
    if len(list_vars) == 1:
        return list_vars[0] - 1

    prev_xs = []
    prev_var = None
    step = 0
    for p in list_vars:
        if prev_var is None:
            prev_var = p
        else:
            result = calculate_polynomial(prev_var, p, step)
            prev_xs.append(result)
            prev_var = T[step]
            step += 1
    # De gevallen \delta = 1 en \delta \geq 2
    if len(prev_xs) >= 2:
        prev_xs[-2] = prev_xs[-2].subs({T[step - 2]: list_vars[-1]})
        del prev_xs[-1]
    elif len(prev_xs) == 1:
        prev_xs[0] = prev_xs[0].subs({T[0]: 1})

    # Het geval \delta \geq 3
    while len(prev_xs) >= 2:
        step -= 1
        prev_xs[-2] = prev_xs[-2].resultant(prev_xs[-1], T[step - 2])
        del prev_xs[-1]
    return prev_xs[0]


# Returns the polynomial for each angle, which is f_n(\cos(\widehat{angle})) where n is the amount of times
# the angle appears in the r-corona
def get_angle_formulae(angles):
    angle_formulae = []
    for angle in angles:
        angle_formulae.append(n_angle(angles[angle], get_angle(*angle)))
    return angle_formulae


margin = 10 ^ (-15)


# Checks which solutions are correct
def verify_solutions(solutions, angles):
    verified_sols = []
    for sol in solutions:
        num_approx = get_total_angle(angles, sol)
        if abs(num_approx - pi * 2) < margin:
            verified_sols.append(sol)
    return verified_sols


# Moves each element of cor i places to the left
def permutate_cyclically(cor, i):
    new_cor = []
    for j in range(len(cor)):
        new_cor.append(cor[(i + j) % len(cor)])
    return new_cor


# Checks if there exists a corona for a given list of angle multiplicities (the vector \vec{h})
def check_angle_coding(ks):
    return ks[1] % 2 == 0 and (ks[0] == 0 or ks[1] != 0 or ks[1] == ks[2] == 0) and \
                              (ks[2] == 0 or ks[1] != 0 or ks[0] == ks[1] == 0)


# Gets all the coronas where needed is the amount of 1's and r's we want
def get_coronas_recursive(needed, cur_list, permutations, current):
    if len(cur_list) == sum(needed):
        cyclic_permutations = []
        for i in range(len(cur_list)):
            cyclic_permutations.append(permutate_cyclically(cur_list, i))
            cyclic_permutations.append([])
            for j in range(len(cur_list)):
                cyclic_permutations[-1].append(cur_list[(i - j) % len(cur_list)])
        cyclic_permutations.sort()
        if cyclic_permutations[0] == cur_list:
            permutations.append(cur_list)
    else:
        for i in range(len(needed)):
            if current[i] < needed[i]:
                new_list = cur_list.copy()
                new_list.append(i)
                new_current = current.copy()
                new_current[i] += 1
                get_coronas_recursive(needed, new_list, permutations, new_current)
    return permutations


# Gets the index for \vec{h}
def get_index(a, b):
    if a == 0 and b == 0:
        return 0
    elif (a == 1 and b == 0) or (a == 0 and b == 1):
        return 1
    elif a == 1 and b == 1:
        return 2


# Gets the multiplicity vector \vec{h} for a permutation
def get_coding(perm):
    cur_k = [0] * 3
    for i in range(len(perm)):
        cur_k[get_index(perm[i], perm[(i+1) % len(perm)])] += 1
    return cur_k


# Gets all the coronas from a given multiplicity vector ks
def get_coronas_from_coding(ks):
    ones = ks[0] + ks[1] // 2
    rrs = ks[1] // 2 + ks[2]
    perms = get_coronas_recursive([ones, rrs], [], [], [0, 0])
    good_perms = []
    for perm in perms:
        cur_k = get_coding(perm)
        if cur_k == list(ks):
            good_perms.append([r_one[i] for i in perm])
    return good_perms


# Gets all the coronas for all multipicity vectors kss
def get_coronas_from_checked_coding(kss):
    output = []
    for ks in kss:
        output.extend(get_coronas_from_coding(ks))
    return output


margin2 = 10 ^ (-5)


# Gets all the possible multiplicity vectors for a 1-corona, given a solution
def get_angles_around_one(solution):
    _r = solution
    angles = [arccos(get_angle(1, 1, 1)).numerical_approx(),
          arccos(get_angle(1, 1, _r)).numerical_approx(),
          arccos(get_angle(_r, 1, _r)).numerical_approx()]
    kss = []
    for k1 in range(int(2*pi / angles[0]) + 2):
        sub = k1*angles[0]
        for k2 in range(int((2*pi - sub) / angles[1]) + 2):
            sub = k1*angles[0] + k2*angles[1]
            for k3 in range(int((2*pi - sub) / angles[2]) + 2):
                sub = k1*angles[0] + k2*angles[1] + k3*angles[2]
                sub2 = sub - 2*pi
                sub2 = sub2.numerical_approx()
                if abs(sub2) < margin2 and check_angle_coding((k1, k2, k3)):
                    kss.append((k1, k2, k3))
    return kss


# Gets all the 1-coronas for a solution
def get_coronas_around_one(solution):
    kss = get_angles_around_one(solution)
    print('kss:', kss)
    return get_coronas_from_checked_coding(kss)


# Get all the possible one-coronas and possibilities for r for a given r-corona
def get_possible_radii(list_radii):
    print('Neem de $r$-corona', list_radii, '.')
    angles = get_angles(list_radii, r)

    add_formula = add_polynomials(f[:len(angles)])
    i = 0
    angle_formulae = get_angle_formulae(angles)
    print('Add formula:', add_formula)
    print('De vergelijkingen voor de hoeken: ', angle_formulae)
    for formula in angle_formulae:
        add_formula = add_formula.subs({f[i]: formula})
        i += 1
    num = add_formula.numerator()
    num = num.factor()
    factors = num._factor_list()

    good_factors = []
    for factor, count in factors:
        if r in factor.arguments() and len(factor.coefficients()) > 1:
            good_factors.append(factor)
    print('Factoren: ', good_factors)
    solutions = []
    for factor in good_factors:
        try:
            solutions.extend(factor.roots(r, ring=AA))
        except Exception:
            try:
                solutions.extend(factor.roots(r, ring=RR))
            except Exception:
                pass
    good_solutions = []
    for solution, count in solutions:
        if 1 > solution > 0:
            good_solutions.append(solution)

    print('Oplossingen: ', good_solutions)
    verified_sols = verify_solutions(good_solutions, angles)
    print('Goede oplossingen: ', verified_sols)
    if len(verified_sols) > 0:
        sol = verified_sols[0]
        print('Mogelijke $1$-coronas: \\\\')
        possible_ones = get_coronas_around_one(sol)
        print(possible_ones)

    print('------------------')
    return verified_sols


# Make all permutations of length length which contain possible_vars, modulo cyclic permutations
# Also checks if such a permutation is a valid r-corona, which means that there must be a r so that the angles
# add up to 2pi
def make_permutation(cur_list, length, all_lists, possible_vars):
    if len(cur_list) >= length:
        if cur_list.count(possible_vars[-1]) == length:
            return
        cyclic_permutations = []
        for i in range(length):
            cyclic_permutations.append([])
            for j in range(length):
                cyclic_permutations[i].append(cur_list[(i + j) % length])
        cyclic_permutations = sort_radii(cyclic_permutations)
        if cyclic_permutations[0] == cur_list and is_valid(cur_list):
            all_lists.append(cur_list)
    else:
        for p in possible_vars:
            new_list = cur_list.copy()
            new_list.append(p)
            make_permutation(new_list, length, all_lists, possible_vars)
        new_list = cur_list.copy()
        new_list.append(1)
        make_permutation(new_list, length, all_lists, possible_vars)


# Makes all the possible r-coronas in the case of 2 radii
def make_all_permutations(possible_vars):
    all_permutations = []
    for i in range(3, 6):
        make_permutation([], i, all_permutations, possible_vars)
    return all_permutations


# Do everything for the case N=2
def do_everything():
    possible_radii = []
    for permutation in make_all_permutations([r]):
        possible_radii.append((permutation, get_possible_radii(permutation)))
    print(possible_radii)


do_everything()
