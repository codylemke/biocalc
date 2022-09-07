import numpy as np


def neighbors(x, param_space, delta=1):  # somewhat based on the find_neighbors function from mlrose
    ns = list()
    for i in range(len(x)):
        center = x[i]
        before = center - delta
        after = center + delta
        if before >= param_space[1]:
            n = x.copy()
            n[i] = before
            ns.append(n)
        if after <= param_space[2]:
            n = x.copy()
            n[i] = after
            ns.append(n)

    return ns


def random_neighbor(x, param_space, delta=1):
    ns = neighbors(x, param_space, delta)
    return ns[np.random.randint(0, len(ns))]


def genetic_algorithm(f, param_space, population_size=50, keep_frac=0.5, max_iters=np.inf, influx_rate=0.2, delta=1,
                      mutation_rate=0.2, curve=True, rounds_to_convergence=20):
    """
        f is a function taking an array as input and returning a float as output
        param_space should be a 3-tuple (dimensions, min, max)
        This algorithm will start with a population of random start sites of size population_size.

        This implementation uses Truncate selection:
            Every iteration, the most fit top half of the population is kept, and the least fit half is destroyed.
            A number of randomly generated new individuals is created proportional to the influx_rate.
            The rest of the destroyed population is replaced by averaging the parameters of random individuals from the
            kept set.
    """

    maxv = float("-inf")
    argmax = None
    fitness_curve = list()

    if int(delta) == delta:  # delta is an integer, so we assume this is a discrete problem
        individuals = np.random.randint(param_space[1], param_space[2] + 1, size=(population_size, param_space[0]))
    else:  # delta is not an integer, so it's ok to start at a non-integer starting point
        # individuals = np.random.randint(param_space[1], param_space[2]+1, size=(population_size,param_space[0]))
        individuals = (param_space[2] - param_space[1]) * np.random.random_sample((population_size, param_space[0])) + \
                      param_space[1]
    # individuals = np.random.randint(param_space[1], param_space[2]+1, size=(population_size,param_space[0]))
    y = [f(x) for x in individuals]

    population = list([y[x], individuals[x]] for x in range(len(y)))
    population.sort(key=lambda x: x[0], reverse=True)
    influx_individuals = int(influx_rate * population_size)
    fitness_curve.append(population[0][0])
    i = 0
    tries = 0
    while (i < max_iters) and (tries < rounds_to_convergence):
        i += 1
        tries += 1
        start = int(len(population) * keep_frac)

        for j in range(0, influx_individuals):
            if int(delta) == delta:  # delta is an integer, so we assume this is a discrete problem
                new_individual = np.random.randint(param_space[1], param_space[2] + 1, size=(param_space[0]))
            else:  # delta is not an integer, so it's ok to start at a non-integer starting point
                # individuals = np.random.randint(param_space[1],
                #   param_space[2]+1, size=(population_size,param_space[0]))
                new_individual = (param_space[2] - param_space[1]) * np.random.random_sample(param_space[0]) + \
                                 param_space[1]
            new_y = f(new_individual)
            population[start + j] = [new_y, new_individual]

        start = start + influx_individuals
        end = int(len(population) * keep_frac)
        for j in range(start, len(population)):
            if param_space[0] <= 1:
                crossover = 0
            else:
                crossover = np.random.randint(0, param_space[0] - 1)
            parent_1 = population[np.random.randint(0, end)][1]
            parent_2 = population[np.random.randint(0, end)][1]

            # print(parent_1)

            offspring = parent_1.copy()
            for k in range(crossover + 1):
                offspring[k] = parent_2[k]

            offspring_y = f(offspring)
            population[j] = [offspring_y, offspring]
        # mutate
        mutate = np.random.random_sample(population_size) < mutation_rate
        for (j, v) in enumerate(mutate):
            if v:
                # print(population[j])
                population[j][1] = random_neighbor(population[j][1], param_space, delta)
                population[j][0] = f(population[j][1])
        population.sort(key=lambda x: x[0], reverse=True)
        if population[0][0] > maxv:
            maxv = population[0][0]
            argmax = population[0][1]
            tries = 0

        fitness_curve.append(maxv)

    if curve:
        return argmax, maxv, np.array(fitness_curve)
    return argmax, maxv
