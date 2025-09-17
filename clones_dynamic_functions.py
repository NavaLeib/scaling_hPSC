from scipy.optimize import curve_fit
from scipy.optimize import fsolve
from scipy.optimize import least_squares
import math
import numpy as np


class Model:
    """
    Simulation (modelling) of the clones dynamics.

    :argument:
    number_of_species = number of species / clones (int)
    mutant_percent  # percentage of mutation initially (float in [0,1])
    number_of_mutations # number of mutations
    self.net_growth_rates_array = net_growth_rates_array  # array (dimension : 2*number of species )
                                  of birth and death rates [rP1, rM1, rP2, rM2...] =
                                        [birth_normal, death_normal, birth_mut1, death_mut1, birth_mut2, death_mut2,...]
    self.interaction = interaction  # scalar of the interaction strength
    self.plating_efficiency = plating_efficiency  # scalar (>=1) captures the plating efficiency of mutations

    Methods:
    time_propagation(self, number_of_steps=10**6 , passaging = 4)
    """

    def __init__(self, number_of_species, mutant_percent, number_of_mutations, net_growth_rates_array, interaction,
                 plating_efficiency, initial_mean_clone_size):
        self.number_of_species = number_of_species  # Number of Species (clones)
        self.mutant_percent = mutant_percent  # percentage of mutation initially
        self.number_of_mutations = number_of_mutations  # number of mutations
        self.net_growth_rates_array = net_growth_rates_array  # array of birth and death rates [rP1, rM1, rP2, rM2...]
        self.interaction = interaction  # scalar of the interaction strength
        self.plating_efficiency = plating_efficiency  # scalar (>=1) captures the plating efficiency of mutations
        self.initial_mean_clone_size = initial_mean_clone_size


    def time_propagation_pass_every_time_EveryPassKeepForAnalysisRemainCells(self, number_of_steps=10 ** 6, passaging=4,
                                                                             percent_of_pass=0):
        """
        Propagate the model with time; simulate and saved the clones compositions in every generation
        When passaing occurs - save only the cells that DID NOT pass

        :param self: the parameters of the model
        :param number_of_steps: the total number of steps/reactions (int)
        :param passaging: the number of days before sampling and passaging
        :param percent_of_pass: the % of passaging

        :return: clones compositions in every generation
        """

        print('*************start************************')

        rP = np.zeros(self.number_of_species)
        rM = np.zeros(self.number_of_species)
        rates_array = self.net_growth_rates_array

        # growth rates:
        # rP[:int(self.number_of_species * self.mutant_percent)] = rates_array[0] * self.plating_efficiency
        rP[int(self.number_of_species * self.mutant_percent):] = rates_array[0]
        #
        # print(rates_array)

        for i in range(self.number_of_mutations):
            rP[int(i * self.number_of_species * self.mutant_percent / self.number_of_mutations)
               :int((i + 1) * self.number_of_species * self.mutant_percent / self.number_of_mutations)] = rates_array[
                2 * i + 2]

        # death rate:
        rM[int(self.number_of_species * self.mutant_percent / self.number_of_mutations):] = rates_array[1]
        for i in range(self.number_of_mutations):
            rM[int(i * self.number_of_species * self.mutant_percent / self.number_of_mutations)
               :int((i + 1) * self.number_of_species * self.mutant_percent / self.number_of_mutations)] = rates_array[
                2 * i + 3]

        X = np.ones(self.number_of_species)  # starting with one individual

        # print(rP, rM)

        # X = np.array([int(i) for i  in np.random.exponential(self.initial_mean_clone_size, self.number_of_species)])  #starting with m individuals ; m expo. distrubted
        # X = np.array([int(i) for i in np.random.poisson(self.initial_mean_clone_size,
        #                                                     self.number_of_species)])  # starting with m individuals ; m poisson. distrubted

        initial_number_of_cells = sum(X)
        print('initial number of cells', initial_number_of_cells, ', initial mutation ratio = ', 1 - sum(
            X[int(self.number_of_species * self.mutant_percent / self.number_of_mutations):]) / initial_number_of_cells)

        # for interaction (affect only on normal cells - suppressed by mutants
        normal_indicator = np.zeros(self.number_of_species)
        normal_indicator[int(self.number_of_species * self.mutant_percent):] = 1

        interaction = self.interaction

        X_X = np.zeros((100, self.number_of_species))
        X_X[0, :] = X
        time_gen = np.zeros(100)
        t = 0
        flag = 0
        rng = np.random.default_rng()

        for n in range(number_of_steps):

            RxPlus = np.array(rP) * np.array(X)
            RxMinus = np.array(rM) * np.array(X) + \
                      interaction * (np.array(X[:int(self.number_of_species * self.mutant_percent)]).sum()) * \
                      np.array(normal_indicator) * np.array(X)
            Rtotal = sum(RxPlus + RxMinus)

            rand = np.random.uniform()
            dt = -math.log(rand) / Rtotal
            t = t + dt

            r = np.random.uniform()

            RR = (np.concatenate((np.array([0]), (np.cumsum(np.concatenate((RxPlus, RxMinus))))))) / Rtotal

            # i = list(abs(RR - r)).index(min(abs(RR - r)))
            temp1 = np.where(r >= np.array(RR))
            temp2 = np.where(r < np.array(RR))
            s = temp1[0][-1]

            if temp1[0][-1] != (temp2[0][0] - 1):
                print('error')

            # print(s)
            X[s % self.number_of_species] = X[s % self.number_of_species] \
                                            + (self.number_of_species > s) - (self.number_of_species <= s)

            # print(s)

            # for s in range(2*NumberOfSpecies):
            #    if r>RR[s] and r<RR[s+1]:
            #        X[s % NumberOfSpecies] = X[s % NumberOfSpecies] + (NumberOfSpecies >= s) - (NumberOfSpecies < s)
            #        break

            # print((math.log2(sum(X) / initial_number_of_cells)) % (1) )

            ## the following line it save the state of the END each day!!
            # more accurate - saves the last moment before the new day, e.g. the last minute before t<24 will be saved
            # as in day 1. the last event before t<48 will be the saved as day 2. etc.
            # If propagation occur it will over-write that day with the new one
            X_X[math.ceil(t / 24), :] = X
            # print(t, int(t/24))
            # if t>24:
            #     print('days = ', int(t/24),
            #        X[:int(self.number_of_species * self.mutant_percent)].sum(axis=0) / (X.sum(axis=0)))

            if t > 24 * passaging * (flag + 1):
                print('passing every ' + str(passaging) + ' days \t t=', t / 24)
                print('mutation ratio before pass',
                      X[:int(self.number_of_species * self.mutant_percent)].sum(axis=0) / (X.sum(axis=0)))
                print('number of cells before pass' + str(X.sum(axis=0)))
                num_to_pass = (X.sum()) * percent_of_pass / 100
                # X = rng.multinomial(num_to_pass,
                #                     X / (X.sum()))  # THAT'S NOT EXACT SAMPLING _ CAN SAMPLE MORE THEN HAVE

                if (t / 24) < 16:
                        X = sample_cells_from_the_pool(int(num_to_pass), pool=X)

                        print('mutation ratio after pass',
                              X[:int(self.number_of_species * self.mutant_percent)].sum(axis=0) / (X.sum(axis=0)))

                        X_X[math.floor(t / 24), :] = X_X[math.floor(t / 24),
                                                     :] - X  ## for passaging time - keep the cells that did not pass - that are the cells
                        ## that are going to be sequanced

                        X_X[
                            X_X < 0] = 0  # sunce it is not actual sampling (can sample more than have - need to verify it doesn't give negative avlues)
                        print('number of cells after pass' + str(X.sum(axis=0)))
                else:
                        print('day 16 - no passaging (keeping all cells)split_batches')
                flag = flag + 1
                print(int(t / 24))
            if int(t / 24) > 16:  # no more than 16 days
                break

            if int(t / (passaging * 24)) == 10:  # no more than 10 passaging to keep

                break

        return X_X

    def split_batches(self, number_of_batches=10, number_of_steps=10 ** 7,
                      passaging=4, percent_of_pass=10):
        """
        at t=0 split into *number_of_batches* batches and save the clones compositions in every generation

        :param self: the parameters of the model
        :param number_of_batches: the number of batches (int)

        :return: clones compositions in every generatio

        """

        X_old = np.ones(self.number_of_species)  # starting with one individual
        X_new = np.zeros(self.number_of_species)  # starting with one individual

        # X = np.array([int(i) for i  in np.random.exponential(self.initial_mean_clone_size, self.number_of_species)])  #starting with m individuals ; m expo. distrubted

        initial_number_of_cells = sum(X_old)
        XX = []

        for i in range(number_of_batches):
            print('batch number', i)
            print(np.sum(X_old), np.sum(X_new))
            X_new = sample_cells_from_the_pool(int(self.number_of_species/number_of_batches), pool=X_old)

            X = X_new.copy()


            rP = np.zeros(self.number_of_species)
            rM = np.zeros(self.number_of_species)
            rates_array = self.net_growth_rates_array

            # growth rates:
            # rP[:int(self.number_of_species * self.mutant_percent)] = rates_array[0] * self.plating_efficiency
            rP[int(self.number_of_species * self.mutant_percent):] = rates_array[0]
            #
            # print(rates_array)

            for i in range(self.number_of_mutations):
                rP[int(i * self.number_of_species * self.mutant_percent / self.number_of_mutations)
                   :int((i + 1) * self.number_of_species * self.mutant_percent / self.number_of_mutations)] = \
                rates_array[2 * i + 2]

            # death rate:
            rM[int(self.number_of_species * self.mutant_percent / self.number_of_mutations):] = rates_array[1]
            for i in range(self.number_of_mutations):
                rM[int(i * self.number_of_species * self.mutant_percent / self.number_of_mutations)
                   :int((i + 1) * self.number_of_species * self.mutant_percent / self.number_of_mutations)] = \
                rates_array[2 * i + 3]

            initial_number_of_cells = sum(X)
            print('initial number of cells', initial_number_of_cells,
                  ', initial mutation ratio = ', 1 - sum(X[int(self.number_of_species * self.mutant_percent / self.number_of_mutations):]) / initial_number_of_cells)

            # for interaction (affect only on normal cells - suppressed by mutants
            normal_indicator = np.zeros(self.number_of_species)
            normal_indicator[int(self.number_of_species * self.mutant_percent):] = 1

            interaction = self.interaction

            X_X = np.zeros((30, self.number_of_species)).astype(int)
            X_X[0, :] = X
            time_gen = np.zeros(30)
            t = 0
            flag = 0
            rng = np.random.default_rng()

            for n in range(number_of_steps):

                RxPlus = np.array(rP) * np.array(X)
                RxMinus = np.array(rM) * np.array(X) + \
                          interaction * (np.array(X[:int(self.number_of_species * self.mutant_percent)]).sum()) * \
                          np.array(normal_indicator) * np.array(X)
                Rtotal = sum(RxPlus + RxMinus)

                rand = np.random.uniform()
                dt = -math.log(rand) / Rtotal
                t = t + dt

                r = np.random.uniform()

                RR = (np.concatenate((np.array([0]), (np.cumsum(np.concatenate((RxPlus, RxMinus))))))) / Rtotal

                # i = list(abs(RR - r)).index(min(abs(RR - r)))
                temp1 = np.where(r >= np.array(RR))
                temp2 = np.where(r < np.array(RR))
                s = temp1[0][-1]

                if temp1[0][-1] != (temp2[0][0] - 1):
                    print('error')

                X[s % self.number_of_species] = X[s % self.number_of_species] \
                                                + (self.number_of_species > s) - (self.number_of_species <= s)

                X_X[math.ceil(t / 24), :] = X

                if t > 24 * passaging * (flag + 1):
                    print('passing every ' + str(passaging) + ' days \t t=', t / 24)
                    print('mutation ratio before pass',
                          X[:int(self.number_of_species * self.mutant_percent)].sum(axis=0) / (X.sum(axis=0)))
                    print('number of cells before pass' + str(X.sum(axis=0)))
                    num_to_pass = (X.sum()) * percent_of_pass / 100
                    # X = rng.multinomial(num_to_pass,
                    #                     X / (X.sum()))  # THAT'S NOT EXACT SAMPLING _ CAN SAMPLE MORE THEN HAVE

                    if (t / 24) < 16:
                        X = sample_cells_from_the_pool(int(num_to_pass), pool=X)

                        print('passaging with plating efficiency', self.plating_efficiency)
                        X[int(self.number_of_species * self.mutant_percent):] = [round(self.plating_efficiency * item) for item in  X[int(self.number_of_species * self.mutant_percent):]]

                        print('mutation ratio after pass',
                              X[:int(self.number_of_species * self.mutant_percent)].sum(axis=0) / (X.sum(axis=0)))

                        X_X[math.floor(t / 24), :] = X_X[math.floor(t / 24),
                                                     :] - X  ## for passaging time - keep the cells that did not pass - that are the cells
                        ## that are going to be sequanced

                        X_X[
                            X_X < 0] = 0  # sunce it is not actual sampling (can sample more than have - need to verify it doesn't give negative avlues)
                        print('number of cells after pass' + str(X.sum(axis=0)))
                    else:
                        print('day 16 - no passaging (keeping all cells)split_batches')
                    flag = flag + 1

                    print(int(t / 24))

                if int(t / 24) > 16:  # no more than 12 days
                    # print(np.sum(X_old), np.sum(X_new))
                    break



            XX.append(X_X)


        return XX

    def split_batches_emergence_of_variants(self, number_of_batches=10, number_of_steps=10 ** 7,
                      passaging=4, percent_of_pass=10, var_emerge_rate = 0.001):
        """
        at t=0 split into *number_of_batches* batches and save the clones compositions in every generation

        :param self: the parameters of the model
        :param number_of_batches: the number of batches (int)

        :return: (clones compositions once a day, growth_rate_of_clones)

        """

        X_old = np.ones(self.number_of_species)  # starting with one individual
        X_new = np.zeros(self.number_of_species)  # starting with one individual

        # X = np.array([int(i) for i  in np.random.exponential(self.initial_mean_clone_size, self.number_of_species)])  #starting with m individuals ; m expo. distrubted

        initial_number_of_cells = sum(X_old)
        XX = []

        for i in range(number_of_batches):
            print('batch number', i)
            print(np.sum(X_old), np.sum(X_new))
            if number_of_batches > 1:
                X_new = sample_cells_from_the_pool(int(self.number_of_species / number_of_batches), pool=X_old)
            else:
                X_new = X_old.copy()

            X = X_new.copy()

            rP = np.zeros(self.number_of_species)
            rM = np.zeros(self.number_of_species)
            rates_array = self.net_growth_rates_array

            # growth rates:
            # rP[:int(self.number_of_species * self.mutant_percent)] = rates_array[0] * self.plating_efficiency
            rP[int(self.number_of_species * self.mutant_percent):] = rates_array[0]
            #
            # print(rates_array)

            for i in range(self.number_of_mutations):
                rP[int(i * self.number_of_species * self.mutant_percent / self.number_of_mutations)
                   :int((i + 1) * self.number_of_species * self.mutant_percent / self.number_of_mutations)] = \
                    rates_array[2 * i + 2]

            # death rate:
            rM[int(self.number_of_species * self.mutant_percent / self.number_of_mutations):] = rates_array[1]
            for i in range(self.number_of_mutations):
                rM[int(i * self.number_of_species * self.mutant_percent / self.number_of_mutations)
                   :int((i + 1) * self.number_of_species * self.mutant_percent / self.number_of_mutations)] = \
                    rates_array[2 * i + 3]

            initial_number_of_cells = sum(X)
            print('initial number of cells', initial_number_of_cells,
                  ', initial mutation ratio = ', 1 - sum(X[
                                                         int(self.number_of_species * self.mutant_percent / self.number_of_mutations):]) / initial_number_of_cells)

            # for interaction (affect only on normal cells - suppressed by mutants
            normal_indicator = np.zeros(self.number_of_species)
            normal_indicator[int(self.number_of_species * self.mutant_percent):] = 1

            interaction = self.interaction

            X_X = np.zeros((30, self.number_of_species)).astype(int)
            X_X[0, :] = X
            time_gen = np.zeros(30)
            t = 0
            flag = 0
            rng = np.random.default_rng()

            for n in range(number_of_steps):

                RxPlus = np.array(rP) * np.array(X)
                RxMinus = np.array(rM) * np.array(X) + \
                          interaction * (np.array(X[:int(self.number_of_species * self.mutant_percent)]).sum()) * \
                          np.array(normal_indicator) * np.array(X)

                RxVar = var_emerge_rate * (np.array(X) * np.array(normal_indicator)).sum()

                Rtotal = (RxPlus).sum() + (RxMinus).sum() + (RxVar)

                rand = np.random.uniform()
                dt = -math.log(rand) / Rtotal
                t = t + dt

                r = np.random.uniform()

                # normal_indicator, (rP, rM), X
                R = np.concatenate((np.array([0]), RxPlus, RxMinus, np.array([RxVar])))

                RR = np.cumsum(R) / Rtotal

                # i = list(abs(RR - r)).index(min(abs(RR - r)))
                temp1 = np.where(r >= np.array(RR))
                temp2 = np.where(r < np.array(RR))
                s = temp1[0][-1]

                # if temp1[0][-1] != (temp2[0][0] - 1):
                #     print('error')

                if s < len(RxPlus):
                    X[s] += 1
                elif (s >= len(RxPlus) and s < (len(RxPlus) + len(RxMinus))):
                    X[s] -=1
                elif s >= (len(RxPlus) + len(RxMinus)):
                    normal_indicator, (rP, rM), X = mutation_emergence(normal_indicator,rates_array,X)
                    X_X = np.c_[ X_X, np.zeros(len(time_gen))] ## adding a col to the X_X - so it can incorporate additional clone
                X_X[math.ceil(t / 24), :] = X

                if t > 24 * passaging * (flag + 1):
                    print('passing every ' + str(passaging) + ' days \t t=', t / 24)
                    print('mutation ratio before pass',
                          1-(np.array(X)*np.array(normal_indicator)).sum(axis=0) / (X.sum(axis=0)))
                    print('number of cells before pass' + str(X.sum(axis=0)))
                    num_to_pass = (X.sum()) * percent_of_pass / 100
                    # X = rng.multinomial(num_to_pass,
                    #                     X / (X.sum()))  # THAT'S NOT EXACT SAMPLING _ CAN SAMPLE MORE THEN HAVE

                    if (t / 24) < 16:
                        X = sample_cells_from_the_pool(int(num_to_pass), pool=X)

                        print('passaging with plating efficiency', self.plating_efficiency)
                        X[int(self.number_of_species * self.mutant_percent):] = [round(self.plating_efficiency * item)
                                                                                 for item in X[
                                                                                             int(self.number_of_species * self.mutant_percent):]]

                        print('mutation ratio after pass',
                          1-(np.array(X)*np.array(normal_indicator)).sum(axis=0) / (X.sum(axis=0)))

                        X_X[math.floor(t / 24), :] = X_X[math.floor(t / 24),
                                                     :] - X  ## for passaging time - keep the cells that did not pass - that are the cells
                        ## that are going to be sequanced

                        X_X[
                            X_X < 0] = 0  # sunce it is not actual sampling (can sample more than have - need to verify it doesn't give negative avlues)
                        print('number of cells after pass' + str(X.sum(axis=0)))
                    else:
                        print('day 16 - no passaging (keeping all cells)split_batches')
                    flag = flag + 1

                    print(int(t / 24))

                if int(t / 24) > 16:  # no more than 12 days
                    # print(np.sum(X_old), np.sum(X_new))
                    break

            XX.append(X_X)

        return XX

    def split_batches_tau_leaping(self, number_of_batches=10, number_of_steps=10 ** 7,
                      passaging=4, percent_of_pass=10):
        """
        at t=0 split into *number_of_batches* batches and save the clones compositions in every generation

        :param self: the parameters of the model
        :param number_of_batches: the number of batches (int)

        :return: clones compositions in every generatio

        """

        X_old = np.ones(self.number_of_species)  # starting with one individual
        X_new = np.zeros(self.number_of_species)  # starting with one individual

        # X = np.array([int(i) for i  in np.random.exponential(self.initial_mean_clone_size, self.number_of_species)])  #starting with m individuals ; m expo. distrubted

        initial_number_of_cells = sum(X_old)
        XX = []

        for i in range(number_of_batches):
            print('batch number', i)
            print(np.sum(X_old), np.sum(X_new))
            if number_of_batches > 1:
                X_new = sample_cells_from_the_pool(int(self.number_of_species/number_of_batches), pool=X_old)
            else:
                X_new = X_old
            X = X_new.copy()


            rP = np.zeros(self.number_of_species)
            rM = np.zeros(self.number_of_species)
            rates_array = self.net_growth_rates_array

            # growth rates:
            # rP[:int(self.number_of_species * self.mutant_percent)] = rates_array[0] * self.plating_efficiency
            rP[int(self.number_of_species * self.mutant_percent):] = rates_array[0]
            #
            # print(rates_array)

            for i in range(self.number_of_mutations):
                rP[int(i * self.number_of_species * self.mutant_percent / self.number_of_mutations)
                   :int((i + 1) * self.number_of_species * self.mutant_percent / self.number_of_mutations)] = \
                rates_array[2 * i + 2]

            # death rate:
            rM[int(self.number_of_species * self.mutant_percent / self.number_of_mutations):] = rates_array[1]
            for i in range(self.number_of_mutations):
                rM[int(i * self.number_of_species * self.mutant_percent / self.number_of_mutations)
                   :int((i + 1) * self.number_of_species * self.mutant_percent / self.number_of_mutations)] = \
                rates_array[2 * i + 3]

            initial_number_of_cells = sum(X)
            print('initial number of cells', initial_number_of_cells,
                  ', initial mutation ratio = ', 1 - sum(X[int(self.number_of_species * self.mutant_percent / self.number_of_mutations):]) / initial_number_of_cells)

            # for interaction (affect only on normal cells - suppressed by mutants
            normal_indicator = np.zeros(self.number_of_species)
            normal_indicator[int(self.number_of_species * self.mutant_percent):] = 1

            interaction = self.interaction

            X_X = np.zeros((30, self.number_of_species)).astype(int)
            X_X[0, :] = X
            time_gen = np.zeros(30)
            t = 0
            flag = 0
            rng = np.random.default_rng()

            for n in range(number_of_steps):

                RxPlus = np.array(rP) * np.array(X)
                RxMinus = np.array(rM) * np.array(X) + \
                          interaction * (np.array(X[:int(self.number_of_species * self.mutant_percent)]).sum()) * \
                          np.array(normal_indicator) * np.array(X)
                Rtotal = sum(RxPlus + RxMinus)

                # rand = np.random.uniform()
                mu_i = np.array(RxPlus)-np.array(RxMinus)
                epsilon = 0.003
                m = np.fmax(epsilon*np.array(X[np.nonzero(X)]),np.ones(shape=len(np.array(X[np.nonzero(X)]))))
                tau = np.nanmin(m/mu_i[np.nonzero(X)]) * 0.1 ## artifically decrease tau by 0.1 to ensure better results
                # tau1 = np.min(X[np.nonzero(X)])/Rtotal
                t = t + tau

                number_of_divisions = np.random.poisson(RxPlus * tau)

                added_number_in_a_leap = sum(number_of_divisions)
                X += number_of_divisions

                X_X[math.ceil(t / 24), :] = X

                if t > 24 * passaging * (flag + 1):
                    print('passing every ' + str(passaging) + ' days \t t=', t / 24)
                    print('mutation ratio before pass',
                          X[:int(self.number_of_species * self.mutant_percent)].sum(axis=0) / (X.sum(axis=0)))
                    print('number of cells before pass' + str(X.sum(axis=0)))
                    num_to_pass = (X.sum()) * percent_of_pass / 100
                    # X = rng.multinomial(num_to_pass,
                    #                     X / (X.sum()))  # THAT'S NOT EXACT SAMPLING _ CAN SAMPLE MORE THEN HAVE

                    if (t / 24) < 16:
                        X = sample_cells_from_the_pool(int(num_to_pass), pool=X)

                        if self.plating_efficiency < 1:
                            print('passaging with plating efficiency', self.plating_efficiency)
                            num_to_pass = X[int(self.number_of_species * self.mutant_percent):].sum(axis=0)*self.plating_efficiency
                            X[int(self.number_of_species * self.mutant_percent):] = sample_cells_from_the_pool(int(num_to_pass), pool=X[int(self.number_of_species * self.mutant_percent):])
                        print('mutation ratio after pass',
                              X[:int(self.number_of_species * self.mutant_percent)].sum(axis=0) / (X.sum(axis=0)))

                        X_X[math.floor(t / 24), :] = X_X[math.floor(t / 24),
                                                     :] - X  ## for passaging time - keep the cells that did not pass - that are the cells
                        ## that are going to be sequanced

                        X_X[
                            X_X < 0] = 0  # sunce it is not actual sampling (can sample more than have - need to verify it doesn't give negative avlues)
                        print('number of cells after pass' + str(X.sum(axis=0)))
                    else:
                        print('day 16 - no passaging (keeping all cells)split_batches')
                    flag = flag + 1

                    print(int(t / 24))

                if int(t / 24) > 16:  # no more than 12 days
                    # print(np.sum(X_old), np.sum(X_new))
                    break



            XX.append(X_X)


        return XX

    def time_propagation_pass_every_time(self, number_of_steps=10 ** 6, passaging=4, percent_of_pass=10):
        """
        Propagate the model with time; simulate and saved the clones compositions in every generation

        :param self: the parameters of the model
        :param number_of_steps: the total number of steps/reactions (int)
        :param passaging: the number of days before sampling and passaging
        :param percent_of_pass: the % of passaging

        :return: clones compositions in every generation
        """

        rP = np.zeros(self.number_of_species)
        rM = np.zeros(self.number_of_species)
        rates_array = self.net_growth_rates_array

        # growth rates:
        # rP[:int(self.number_of_species * self.mutant_percent)] = rates_array[0] * self.plating_efficiency
        rP[int(self.number_of_species * self.mutant_percent):] = rates_array[0]
        #
        # print(rates_array)

        for i in range(self.number_of_mutations):
            rP[int(i * self.number_of_species * self.mutant_percent / self.number_of_mutations)
               :int((i + 1) * self.number_of_species * self.mutant_percent / self.number_of_mutations)] = rates_array[
                2 * i + 2]

        # death rate:
        rM[int(self.number_of_species * self.mutant_percent / self.number_of_mutations):] = rates_array[1]
        for i in range(self.number_of_mutations):
            rM[int(i * self.number_of_species * self.mutant_percent / self.number_of_mutations)
               :int((i + 1) * self.number_of_species * self.mutant_percent / self.number_of_mutations)] = rates_array[
                2 * i + 3]

        X = np.ones(self.number_of_species)  # starting with one individual

        # X = np.array([int(i) for i  in np.random.exponential(self.initial_mean_clone_size, self.number_of_species)])  #starting with m individuals ; m expo. distrubted

        initial_number_of_cells = sum(X)
        print('initial number of cells', initial_number_of_cells)

        # for interaction (affect only on normal cells - suppressed by mutants
        normal_indicator = np.zeros(self.number_of_species)
        normal_indicator[int(self.number_of_species * self.mutant_percent):] = 1

        interaction = self.interaction

        X_X = np.zeros((100, self.number_of_species))
        X_X[0, :] = X
        time_gen = np.zeros(100)
        t = 0
        flag = 0
        rng = np.random.default_rng()

        for n in range(number_of_steps):

            RxPlus = np.array(rP) * np.array(X)
            RxMinus = np.array(rM) * np.array(X) + \
                      interaction * (np.array(X[:int(self.number_of_species * self.mutant_percent)]).sum()) * \
                      np.array(normal_indicator) * np.array(X)
            Rtotal = sum(RxPlus + RxMinus)

            rand = np.random.uniform()
            dt = -math.log(rand) / Rtotal
            t = t + dt

            r = np.random.uniform()

            RR = (np.concatenate((np.array([0]), (np.cumsum(np.concatenate((RxPlus, RxMinus))))))) / Rtotal

            # i = list(abs(RR - r)).index(min(abs(RR - r)))
            temp1 = np.where(r >= np.array(RR))
            temp2 = np.where(r < np.array(RR))
            s = temp1[0][-1]

            if temp1[0][-1] != (temp2[0][0] - 1):
                print('error')

            # print(s)
            X[s % self.number_of_species] = X[s % self.number_of_species] \
                                            + (self.number_of_species >= s) - (self.number_of_species < s)

            # for s in range(2*NumberOfSpecies):
            #    if r>RR[s] and r<RR[s+1]:
            #        X[s % NumberOfSpecies] = X[s % NumberOfSpecies] + (NumberOfSpecies >= s) - (NumberOfSpecies < s)
            #        break

            # print((math.log2(sum(X) / initial_number_of_cells)) % (1) )

            X_X[int(t / 24), :] = X
            # if t>24:
            #     print('days = ', int(t/24),
            #        X[:int(self.number_of_species * self.mutant_percent)].sum(axis=0) / (X.sum(axis=0)))

            if t > 24 * passaging * (flag + 1):
                print('passing every ' + str(passaging) + ' days \t t=', t / 24)
                print('mutation ratio',
                      X[:int(self.number_of_species * self.mutant_percent)].sum(axis=0) / (X.sum(axis=0)))
                print('number of cells before pass' + str(X.sum(axis=0)))
                num_to_pass = (X.sum()) * percent_of_pass / 100
                # X = rng.multinomial(num_to_pass,
                #                     X / (X.sum()))  # THAT'S NOT EXACT SAMPLING _ CAN SAMPLE MORE THEN HAVE

                X = sample_cells_from_the_pool(num_to_pass, pool=X)
                print('number of cells after pass' + str(X.sum(axis=0)))
                flag = flag + 1

            if int(t / 24) == 16:  # no more than 100 generations to keep
                break

            if int(t / (passaging * 24)) == 10:  # no more than 10 passaging to keep

                break

        print('saving data')
        np.save("sample.npy", X_X)

        return X_X




def sample_cells_from_the_pool(num_to_pass, pool):
    rng = np.random.default_rng()

    X_before = pool

    sampled_X = np.zeros(shape=np.shape(pool))

    for i in range(num_to_pass):
        prob_to_be_sampled = X_before / (X_before.sum())

        sampled_cell_id = rng.choice(len(X_before), size=1, p=prob_to_be_sampled)

        sampled_X[sampled_cell_id] += 1

        X_before[sampled_cell_id] -= 1

    return sampled_X

def mutation_emergence (normal_indicator, rates_array, X):
    """

    this function is called with a rate of mutation (predefined rate that a wt cell with be mut cell)

    Args:
        normal_indicator: an array with the size of num_of_clones - indicates if it's wt (1) or var (0)
        rates_array:
        X: state of clones

    Returns:

    """
    # rng = np.random.default_rng()
    # prob_to_be_sampled = normal_indicator / (normal_indicator.sum())
    # sampled_cell_id = rng.choice(len(normal_indicator), size=1, p=prob_to_be_sampled)
    # normal_indicator[sampled_cell_id] = 1 ## changes one wt


    normal_indicator = np.array(list(normal_indicator)+[0]) ## this add var indicator to the end of the list
    rP = rates_array[0] * np.array(normal_indicator) + rates_array[2] * (1-np.array(normal_indicator))
    rM = rates_array[1] * np.array(normal_indicator) + rates_array[3] * (1-np.array(normal_indicator))
    X = np.array(list(X)+[1])

    return normal_indicator, (rP,rM), X