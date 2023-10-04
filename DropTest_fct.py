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

       ## fct to RUN SIMULATION and calculate survival +  N(t) for a set number of repeats
    def calc_survival_prob_total_nr_bact_diff_part(self, step, spec_time, total_sim):
        ## create list with part fact as column names
        total_prob = pd.DataFrame()
        df_total_mass = pd.DataFrame()

        for nr_sim in range(0, total_sim, step):
            prob_diff_part = []
            part_fact = []
            print('Simulation nr:', nr_sim)

            ## simulate droplets for different partitioning factors
            for i in range(len(droplet_list)):
                total_drop_nr = droplet_list[i]
                print('total droplets', total_drop_nr)
                m_fact = droplet_list[-1]/total_drop_nr
                ## volume for one small droplet is 1e-7, but we're starting the simulation from the one big droplet and further getting to smaller droplets
                new_volume = variables.volume * droplet_list[-1]/ droplet_list[i]

                if total_drop_nr == 0:
                    print("Error in partitioning; the partitioning factor is so small that you're trying to simulate 0 droplets")
                else:
                    part_fact.append(total_drop_nr)

                    ## simulate growth with new params
                    strain_R = strain(m_fact)
                    Droplet_exp = droplets_R(total_drop_nr, strain_R, variables.AB_conc, new_volume, m_fact)
                    Droplet_exp.run(loading, growth)
                    #print(i)
                    ## calculate the total nr of bacteria in all droplets, if any survived, prob survival = 1
                    Droplet_exp.countTotalMass(growth)
                    #Droplet_exp.plots(growth)
                    Droplet_exp.save('initialN{}_growthrate{}_MIC{}_totaldropnr{}_ABconc{}_'
                                     'dt{}_loading{}_growth{}.csv'.format(initialN, growthrate, MIC, total_drop_nr,variables.AB_conc, dt, loading, growth),
                                     'ABconc{}_loading{}_growth{}.csv'.format(variables.AB_conc, loading, growth), 'Time.csv',  variables.AB_conc)

                    nr_bact_each_ts = Droplet_exp.total_mass
                    ## append the nr of bacteria to dataframe with N(t) vs part factor
                    df_total_mass['{} {}'.format(nr_sim, total_drop_nr)] = nr_bact_each_ts
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
        columns = list(df_total_mass.columns)

        columns_split = []
        for word in columns:
            tmp = word.split(' ')
            #print(tmp)
            tmp = [int(x) for x in tmp]
            columns_split.append(tmp)
        #print('split columns', columns_split)

        columns_sorted = sorted(columns_split, key=operator.itemgetter(1))
        columns_sorted_joined = []
        for word in columns_sorted:
            word = str(word[0]) +' '+ str(word[1])

            columns_sorted_joined.append(word)
        df_total_mass = df_total_mass.reindex(columns_sorted_joined, axis=1)


        ## calculate the survival fraction and add it to the last row of the dataframe
        total_prob.loc['surv_Frac'] = total_prob.sum() / (len(total_prob.axes[0]))
        ## make new dataframe with the last row of the df as survival fraction
        surv_df = total_prob.tail(1)

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

        folder_name = 'output/dropnr_{}_loading_{}_growth_{}_initialN_{}_abconc_{}'.format(variables.droplet_list[-1], variables.loading, variables.growth, variables.initialN, variables.AB_conc)
        path = os.path.join(curr_path, folder_name)
        isExist = os.path.exists(path)

        if not isExist:
            # Create a new directory because it does not exist
            os.makedirs(path, exist_ok=True)
            print("The new directory is created!")

        os.chdir(path)

        pd.DataFrame(surv_df).to_csv('survival_fraction_growth_{}_loading_{}_ABconc{}.csv'.format(growth, loading, variables.AB_conc), index=None)
        pd.DataFrame(avg_nr_bact).to_csv('average_df_growth_{}_loading_{}_ABconc{}.csv'.format(growth, loading, variables.AB_conc), index=None)
        pd.DataFrame(df_total_mass).to_csv('df_growth_{}_starting_nr_drops_{}_ABconc{}.csv'.format(growth, variables.droplet_list[-1], variables.AB_conc), index=None)
        np.savetxt('part_fact.txt', part_fact, delimiter=',')  # X is an array

        # print("--- %s seconds ---" % (time.time() - start_time))

#simulate = DropTest()
#simulate.calc_survival_prob_total_nr_bact_diff_part(step, spec_time, total_sim)

