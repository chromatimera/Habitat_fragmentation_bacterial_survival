import variables
from strain import strain
from Droplet_R import droplets_R
from variables import *
from decimal import *
import pandas as pd
import numpy as np
import warnings

##ignore warnings
warnings.filterwarnings("ignore")
getcontext().prec = 50


import time
start_time = time.time()


class DropTest(object):
    def run(self): ## simulation for one AB concentration; used in troubleshooting

        strain_R = strain(nr_drops_total_mass=1)
        Droplet_exp = droplets_R(total_drop_nr, strain_R, AB_conc, volume)
        Droplet_exp.run(loading, growth)
        Droplet_exp.save('Ni{}'
                         '_MIC{}_totaldropnr{}_ABconc{}_'
                         '_loading{}_growth{}.csv'.format(initialN, growthrate, MIC, total_drop_nr, AB_conc, dt,
                                                              loading, growth), 'ABconc{}_loading{}_growth{}.csv'.format(AB_conc, loading, growth),'Time_list.csv', AB_conc)

        #print("--- %s seconds ---" % (time.time() - start_time))
        Droplet_exp.plots(growth)
        Droplet_exp.countSurvival(growth)

    def test_dt(self, dtmin, dtmax, stepdt): ## will be removed; fct to test difference in AB degradation between different dts -
        # can estimate error in linear approximation; uncomment Deg_list in the Exp_R file
        new = pd.DataFrame()

        # loop for variables
        for i in range(dtmin, dtmax, stepdt):
            variables.dt = 1 / (10 ** i)
            print(nr_timesteps)
            print(variables.dt)
            variables.t_end = round(variables.dt * nr_timesteps, i)
            print(t_end)
            print('dt = ', variables.dt, 't_end = ', variables.t_end)
            strain_R = strain(nr_drops_total_mass=1)
            Droplet_exp = droplets_R(strain_R, AB_conc)
            Droplet_exp.run(loading, growth)
            additional = pd.DataFrame({"" + str(variables.dt) + "": Droplet_exp.deg_list})
            new = pd.concat([new, additional], axis=1)
            Droplet_exp.save(
                'initialN{}_growthrate{}_MIC{}_totaldropnr{}_ABconc{}_dt{}_loading{}_growth{}.csv'.format(initialN,
                                                                                                          growthrate,
                                                                                                          MIC,
                                                                                                          total_drop_nr,
                                                                                                          AB_conc,
                                                                                                          variables.dt,
                                                                                                          loading,
                                                                                                          growth),
                'ABconc{}_loading{}_growth{}.csv'.format(AB_conc, loading, growth),'Time{}.csv'.format(variables.dt), AB_conc)

        # Get path of current folder
        curr_path = os.getcwd()
        # Check whether the specified path exists or not
        folder_name = 'output/df_test_dt_nr_timesteps_{}.csv'.format(nr_timesteps)
        path = os.path.join(curr_path, folder_name)
        isExist = os.path.exists(path)

        if not isExist:
            # Create a new directory because it does not exist
            os.makedirs(path, exist_ok=True)
            print("The new directory is created!")

        os.chdir(path)
        pd.DataFrame(new).to_csv('df_degradation_{}_nr_timesteps_{}.csv'.format(degradation, nr_timesteps), index = None)

        os.chdir(curr_path)

    ## fct to calculate survival for a set nr of repeats
    def calc_survival_prob_diff_part(self, partmin, partmax, step, spec_time, total_sim):

        ## create list with part fact as column names
        total_prob = pd.DataFrame()
        survival_fraction = []
        for nr_sim in range(0, total_sim, step):
            prob_diff_part = []
            part_fact = []
            ## simulate droplets for different partitioning factors
            for i in range(partmin, partmax, step):
                 for j in range(partmin, partmax, step):
                    ## set new values for diff part factor
                    new_i = 5**i * 2**j
                    new_nr_drops_total_mass = new_i
                    new_volume = variables.volume * new_nr_drops_total_mass
                    total_drop_nr = round(variables.total_drop_nr /new_nr_drops_total_mass)
                    part_fct = 1 / total_drop_nr
                    part_fact.append(part_fct)
                    ## simulate growth with new params
                    strain_R = strain(new_nr_drops_total_mass)
                    Droplet_exp = droplets_R(total_drop_nr, strain_R, AB_conc, new_volume)  # 0.5, 300
                    Droplet_exp.run(loading, growth)


                    ## calculate the total nr of bacteria in all droplets, if any survived, prob survival = 1
                    Droplet_exp.countTotalMass(growth)
                    nr_bact_each_ts = Droplet_exp.total_mass
                    index = int(spec_time/variables.dt)
                    if nr_bact_each_ts[index] == 0:
                        prob_survival = 0
                    else:
                        prob_survival = 1
                    ## append probs to a list with diff probs for diff part factors
                    prob_diff_part.append(prob_survival)
            ## append row at the end of each simulation for each partition factor
            additional = pd.DataFrame(data=[prob_diff_part])
            total_prob = pd.concat([total_prob, additional])
        total_prob.columns = part_fact
        print(total_prob)
        ## calculate the survival fraction and make new dataframe
        total_prob.loc['surv_Frac'] = total_prob.sum()/ (len(total_prob.axes[0]))
        print(total_prob)

        surv_df = total_prob.tail(1)
        pd.DataFrame(surv_df).to_csv('output/survival_fraction_growth_{}_loading_{}_ABconc{}.csv'.format(growth, loading, AB_conc), index = None)
        # print("--- %s seconds ---" % (time.time() - start_time))

    ## fct to calculate survival +  N(t) for a set nr of repeats
    def calc_survival_prob_total_nr_bact_diff_part(self, partmin, partmax, step, spec_time, total_sim):
        ## create list with part fact as column names
        total_prob = pd.DataFrame()
        df_total_mass = pd.DataFrame()

        for nr_sim in range(0, total_sim, step):
            prob_diff_part = []
            part_fact = []
            print('Simulation nr:', nr_sim)
            ## simulate droplets for different partitioning factors
            for i in range(partmin, partmax, step):
                #print('i', i)

                for j in range(partmin, partmax, step):
                    #print('j', j)

                    ## set new values for diff part factor
                    new_i = 5 ** i * 2 ** j
                    #print('new_i', new_i)

                    new_nr_drops_total_mass = new_i
                    new_volume = variables.volume * new_nr_drops_total_mass
                    total_drop_nr = round(variables.total_drop_nr / new_nr_drops_total_mass)

                    print('total_droplets', total_drop_nr)
                    if total_drop_nr == 0:
                        print("Error in partitioning; the partitioning factor is so small that you're trying to simulate 0 droplets")
                    else:

                        part_fct = 1 / total_drop_nr
                        #print('part_fact', part_fct)

                        part_fact.append(part_fct)

                        ## simulate growth with new params
                        strain_R = strain(new_nr_drops_total_mass)
                        Droplet_exp = droplets_R(total_drop_nr, strain_R, AB_conc, new_volume)  # 0.5, 300
                        Droplet_exp.run(loading, growth)

                        ## calculate the total nr of bacteria in all droplets, if any survived, prob survival = 1
                        Droplet_exp.countTotalMass(growth)
                        nr_bact_each_ts = Droplet_exp.total_mass
                        ## append the nr of bacteria to dataframe with N(t) vs part factor
                        df_total_mass['{}_{}'.format(part_fct, nr_sim)] = nr_bact_each_ts

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
        ## calculate the survival fraction and add it to the last row of the dataframe
        total_prob.loc['surv_Frac'] = total_prob.sum() / (len(total_prob.axes[0]))
        #print(total_prob)

        ## make new dataframe with the last row of the df as survival fraction
        surv_df = total_prob.tail(1)

        ##start building the average N(t) over simulations
        total_nr_bact = np.zeros((df_total_mass.shape[0], len(part_fact)))

        ## for each partition factor, calculate the sum over the simulation of N(t)
        for i in range(0, len(part_fact), 1):
           for j in range(0, total_sim, 1):
               k = i + j * len(part_fact)
               total_nr_bact[:, i] += df_total_mass.iloc[:, k]

        ## at this stage we have a np array with a sum of N(t) across all iterations -> we need to divide it by the nr of sim
        nr_simu = np.array(total_sim)
        avg_nr_bact = np.divide(total_nr_bact, nr_simu)
        avg_nr_bact = pd.DataFrame(avg_nr_bact, columns = part_fact)

       ## save all df
        pd.DataFrame(surv_df).to_csv('output/survival_fraction_growth_{}_loading_{}_ABconc{}.csv'.format(growth, loading, AB_conc), index=None)
        pd.DataFrame(avg_nr_bact).to_csv('output/average_df_growth_{}_loading_{}_ABconc{}.csv'.format(growth, loading, AB_conc), index=None)
        pd.DataFrame(df_total_mass).to_csv('output/df_growth_{}_starting_nr_drops_{}.csv'.format(growth, variables.total_drop_nr), index=None)
        np.savetxt('output/part_fact.txt', part_fact, delimiter=',')  # X is an array

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
                print('total_droplets', total_drop_nr)
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
                    print('part fct', part_fct)
                    df_total_mass['{}'.format(part_fct)] = nr_bact_each_ts
                    ## append all parition factors, next step will transform this into partition factors
                    part_fact.append(part_fct)


        #
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
#simulate.run()
#simulate.test_dt(0, 10, 1)
#simulate.count_total_mass(part_min, part_max, step)
#simulate.test_surv_frac_diff_ab_conc(abmin, abmax, step)
simulate.calc_survival_prob_total_nr_bact_diff_part(part_min, part_max, step, spec_time, total_sim)
#simulate.count_total_mass_diff_ab(part_min, part_max, abmin, abmax, step)
#simulate.count_total_mass(abmin, abmax, step)
