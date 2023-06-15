import variables
from strain import strain
from Droplet_R import droplets_R
from variables import *
from decimal import *
import pandas as pd
import numpy as np
import warnings
import operator

##ignore warnings
warnings.filterwarnings("ignore")
getcontext().prec = 50

import time
start_time = time.time()


class DropTest(object):

       ## fct to calculate survival +  N(t) for a set nr of repeats
    def calc_survival_prob_total_nr_bact_diff_part(self, partmin, partmax, step, spec_time, total_sim):
        ## create list with part fact as column names
        total_prob = pd.DataFrame()
        df_total_mass = pd.DataFrame()

        for nr_sim in range(0, total_sim, step):
            prob_diff_part = []
            part_fact = []
            ps_array = np.empty([partmax*partmax,2])
            print('Simulation nr:', nr_sim)

            ## simulate droplets for different partitioning factors
            n=-1
            for i in range(partmin, partmax, step):  # change?? for k in (5 ** i * 2 ** j for i in range(partmin, partmax, step) for j in range(partmin, partmax, step)): OR for k in (5 ** ((n - partmin) * step) * 2 ** ((m - partmin) * step) for n in range(partmin, partmax, step) for m in range(partmin, partmax, step)):

                for j in range(partmin, partmax, step):
                    n=n+1
                    new_i = 5 ** i * 2 ** j
                    new_nr_drops_total_mass = new_i
                    new_volume = variables.volume * new_nr_drops_total_mass
                    total_drop_nr = round(variables.total_drop_nr / new_nr_drops_total_mass)

                    if total_drop_nr == 0:
                        print("Error in partitioning; the partitioning factor is so small that you're trying to simulate 0 droplets")
                    else:

                        part_fct = total_drop_nr
                        print('m (nr of subvolumes)', part_fct)

                        part_fact.append(part_fct)

                        ## simulate growth with new params
                        strain_R = strain(new_nr_drops_total_mass)
                        Droplet_exp = droplets_R(total_drop_nr, strain_R, AB_conc, new_volume, new_nr_drops_total_mass)  # 0.5, 300
                        Droplet_exp.run(loading, growth)

                        ## calculate the total nr of bacteria in all droplets, if any survived, prob survival = 1
                        Droplet_exp.countTotalMass(growth)
                        #Droplet_exp.calc_tau_det(Droplet_exp.N_r_array)
                        #Droplet_exp.plots(growth)

                        rho_T, N_T, ps, bigPs = Droplet_exp.calc_theo_survival_prob(Droplet_exp.N_r_array)
                        print('bigPs', bigPs)
                        #rho_T, N_T, ps, bigPs = Droplet_exp.calc_theo_survival_prob(Droplet_exp.N_r_array)
                        #ps_array[n,0]=ps
                        #ps_array[n, 1] = bigPs

                        nr_bact_each_ts = Droplet_exp.total_mass
                        ## append the nr of bacteria to dataframe with N(t) vs part factor
                        df_total_mass['{} {}'.format(part_fct, nr_sim)] = nr_bact_each_ts
                        index = int(spec_time / variables.dt)

                        if nr_bact_each_ts[index] == 0:
                            prob_survival = 0
                        else:
                            prob_survival = 1

                        ## append probs to a list with diff probs for diff part factors
                        prob_diff_part.append(prob_survival)

                ## append row at the end of each simulation for each partition factor
            prob_part_per_iteration = pd.DataFrame(data=[prob_diff_part])
            total_prob = pd.concat([total_prob, prob_part_per_iteration])

        total_prob.columns = part_fact
        #print('total prob before sorting ', total_prob.columns.tolist())
        #print('part fact before sorting', part_fact)
        part_fact = sorted(part_fact)

        #print('df columns', df_total_mass.columns.tolist())
        columns = list(df_total_mass.columns)

        columns_split = []
        for word in columns:
            tmp = word.split(' ')
            #print(tmp)
            tmp = [int(x) for x in tmp]
            columns_split.append(tmp)
        #print('split columns', columns_split)

        columns_sorted = sorted(columns_split, key=operator.itemgetter(0))
        #print('sorted columns', columns_sorted)
        columns_sorted_joined = []
        for word in columns_sorted:
            #print(word)
            #print(len(word))
            word = str(word[0]) +' '+ str(word[1])
            #print(word)

            columns_sorted_joined.append(word)
        #print('columns sorted joined', columns_sorted_joined)
        df_total_mass = df_total_mass.reindex(columns_sorted_joined, axis=1)
        #print('new df columns list', df_total_mass.columns.tolist())
        #print('new total prob columns list', sorted(total_prob.columns))
        total_prob = total_prob.reindex(sorted(total_prob.columns, reverse=True), axis = 1)
        ## calculate the survival fraction and add it to the last row of the dataframe
        total_prob.loc['surv_Frac'] = total_prob.sum() / (len(total_prob.axes[0]))
        #print('total prob after sorting ', total_prob.columns.tolist())
        #print('part fact after sorting', part_fact)
        ## make new dataframe with the last row of the df as survival fraction
        surv_df = total_prob.tail(1)
        #surv_df = surv_df.reindex(sorted(surv_df.columns), axis=1)
        #print(surv_df.head(10))

        ##start building the average N(t) over simulations
        total_nr_bact = np.zeros((df_total_mass.shape[0], len(part_fact)))

        ## for each partition factor, calculate the sum over the simulation of N(t)
        for i in range(0, len(part_fact), 1):
           for j in range(0, total_sim, 1):
               k = i * total_sim + j
               total_nr_bact[:, i] += df_total_mass.iloc[:, k]
           #print('stop')

        ## at this stage we have a np array with a sum of N(t) across all iterations -> we need to divide it by the nr of sim
        nr_simu = np.array(total_sim)
        avg_nr_bact = np.divide(total_nr_bact, nr_simu)
        avg_nr_bact = pd.DataFrame(avg_nr_bact, columns = part_fact)


       ## save all df
        # Get path of current folder
        curr_path = os.getcwd()
        # print(curr_path)

        # Check whether the specified path exists or not

        folder_name = 'output/dropnr_{}_loading_{}_growth_{}_initialN_{}_abconc_{}'.format(
            variables.total_drop_nr, variables.loading, variables.growth, variables.initialN, AB_conc)
        path = os.path.join(curr_path, folder_name)
        isExist = os.path.exists(path)

        if not isExist:
            # Create a new directory because it does not exist
            os.makedirs(path, exist_ok=True)
            print("The new directory is created!")

        os.chdir(path)

        pd.DataFrame(surv_df).to_csv('survival_fraction_growth_{}_loading_{}_ABconc{}.csv'.format(growth, loading, AB_conc), index=None)
        pd.DataFrame(avg_nr_bact).to_csv('average_df_growth_{}_loading_{}_ABconc{}.csv'.format(growth, loading, AB_conc), index=None)
        pd.DataFrame(df_total_mass).to_csv('df_growth_{}_starting_nr_drops_{}_ABconc{}.csv'.format(growth, variables.total_drop_nr, AB_conc), index=None)
        pd.DataFrame(ps_array).to_csv('theory_growth_{}_starting_nr_drops_{}_ABconc{}.csv'.format(growth, variables.total_drop_nr, AB_conc), index=None)
        np.savetxt('part_fact.txt', part_fact, delimiter=',')  # X is an array

        # print("--- %s seconds ---" % (time.time() - start_time))

    ## Test the survival fraction of droplets for different factors, ab conc and calculate the total mass
    ## function to sumulate growth and make df used for plotting Nf and Ni versus partitioning factor
    def count_total_mass(self, partmin, partmax, step):
        df_total_mass = pd.DataFrame()
        # loop for different partitioning factors:
        # loop for partitioning factors
        part_fact = []
        for i in range(partmin, partmax, step):
            for j in range(partmin, partmax, step):
                new_i = 5 ** i * 2 ** j
                new_nr_drops_total_mass = new_i
                new_volume = variables.volume * new_nr_drops_total_mass
                total_drop_nr = round(variables.total_drop_nr / new_nr_drops_total_mass)
                #print('total_droplets', total_drop_nr)
                if total_drop_nr == 0:
                    print("Error in partitioning; the partitioning factor is so small that you're trying to simulate 0 droplets")
                else:

                    strain_R = strain(new_nr_drops_total_mass)
                    Droplet_exp = droplets_R(total_drop_nr, strain_R, AB_conc, new_volume)  # 0.5, 300
                    Droplet_exp.run(loading, growth)
                    Droplet_exp.countTotalMass(growth)

                    ##this is the total number of bacteria at each timestep
                    nr_bact_each_ts = Droplet_exp.total_mass
                    part_fct = 1/total_drop_nr
                    #print('part fct', part_fct)
                    df_total_mass['{}'.format(part_fct)] = nr_bact_each_ts
                    ## append all parition factors, next step will transform this into partition factors
                    part_fact.append(part_fct)

        ### insert the time array into the dataframe
        df_total_mass.insert(loc=0, column='Time', value = np.linspace(t_start, t_end, num=round(t_end/dt)))
        ## See below how the dataframe should look like:
        ## Header Time      x0.1 x0.5 x1
        ##          0        100   100   100 ### Nr of bacteria for t_start
        ##          .         .     .     .
        ##          .         .     .     .
        ##          .         .     .     .
        ##         300        0     2     0


        pd.DataFrame(df_total_mass).to_csv('output/df_growth_{}_starting_nr_drops_{}.csv'.format(growth, variables.total_drop_nr), index = None)
        np.savetxt('output/part_fact.txt', part_fact, delimiter=',')  # X is an array

        print("--- %s seconds ---" % (time.time() - start_time))

        ## function to sumulate growth and make df used for plotting Nf and Ni versus partitioning factor FOR DIFFERENT AB CONC
    def count_total_mass_diff_ab(self, partmin, partmax, abmin, abmax, step):
        df_total_mass = pd.DataFrame()
        ab_list = []
        for ab in range(abmin, abmax, step):
            # loop for different partitioning factors:
            # loop for partitioning factors
            part_fact = []
            AB_conc = ab
            ab_list.append(AB_conc)
            print(ab)
            for i in range(partmin, partmax, step):
                for j in range(partmin, partmax, step):
                    new_i = 5 ** i * 2 ** j
                    new_nr_drops_total_mass = new_i
                    new_volume = variables.volume * new_nr_drops_total_mass
                    total_drop_nr = round(variables.total_drop_nr / new_nr_drops_total_mass)
                    print('total_droplets', total_drop_nr)
                    if total_drop_nr == 0:
                        print(
                            "Error in partitioning; the partitioning factor is so small that you're trying to simulate 0 droplets")
                    else:

                        strain_R = strain(new_nr_drops_total_mass)
                        Droplet_exp = droplets_R(total_drop_nr, strain_R, AB_conc, new_volume)  # 0.5, 300
                        Droplet_exp.run(loading, growth)
                        Droplet_exp.countTotalMass(growth)

                        ##this is the total number of bacteria at each timestep
                        nr_bact_each_ts = Droplet_exp.total_mass
                        part_fct = 1 / total_drop_nr
                        print('part fct', part_fct)
                        df_total_mass['p_{}_ab_{}'.format(part_fct, ab)] = nr_bact_each_ts
                        ## append all parition factors, next step will transform this into partition factors
                        part_fact.append(part_fct)

        print(df_total_mass)
        ### insert the time array into the dataframe
        df_total_mass.insert(loc=0, column='Time', value=np.linspace(t_start, t_end, num=round(t_end / dt)))
        ## See below how the dataframe should look like:
        ## Header Time      x0.1_ab _0 x0.5_ab_0 x1_ab_0 x1_ab_1 ...
        ##          0        100   100   100 ### Nr of bacteria for t_start
        ##          .         .     .     .
        ##          .         .     .     .
        ##          .         .     .     .
        ##         300        0     2     0

        pd.DataFrame(df_total_mass).to_csv('output/df_growth_{}_starting_nr_drops_{}.csv'.format(growth, variables.total_drop_nr), index=None)
        np.savetxt('output/part_fact.txt', part_fact, delimiter=',')  # X is an array
        np.savetxt('output/ab_conc.txt', ab_list, delimiter=',')  # X is an array

        print("--- %s seconds ---" % (time.time() - start_time))


simulate = DropTest()
#simulate.count_total_mass(part_min, part_max, step)
simulate.calc_survival_prob_total_nr_bact_diff_part(part_min, part_max, step, spec_time, total_sim)

