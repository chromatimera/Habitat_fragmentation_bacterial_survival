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

    def test_surv_frac_diff_ab_conc(self, abmin, abmax, step):
        new = pd.DataFrame()

        # loop for var
        for i in range(abmin, abmax, step):
            new_ab_conc = i
            print('antib_conc', new_ab_conc)
            strain_R = strain(nr_drops_total_mass= 1 )
            Droplet_exp = droplets_R(total_drop_nr, strain_R, new_ab_conc, variables.volume)
            Droplet_exp.run(loading, growth)
            Droplet_exp.countSurvival(growth)
            #surv_frac = Droplet_exp.Res_survival_fraction
            #print(surv_frac)
            additional = pd.DataFrame({"" + str(new_ab_conc) + "": [Droplet_exp.Res_survival_fraction]})
            # print(additional)
            new = pd.concat([new, additional], axis=1)
             #Droplet_exp.save(
            #     'initialN{}_growthrate{}_MIC{}_totaldropnr{}_ABconc{}_dt{}_loading{}_growth{}.csv'.format(initialN,
            #                                                                                               growthrate,
            #                                                                                               MIC,
            #                                                                                               total_drop_nr,
            #                                                                                               new_ab_conc,
            #                                                                                               dt,
            #                                                                                               loading,
            #                                                                                               growth),
            #     'ABconc{}_loading{}_growth{}.csv'.format(new_ab_conc, loading, growth),'Time_list.csv', Nsat, total_drop_nr, loading,
            #     growth, initialN, new_ab_conc, growthrate, dt)

            print("--- %s seconds ---" % (time.time() - start_time))
            # Droplet_exp.normed_histograms(50)
            #print(new)
            pd.DataFrame(new).to_csv('output/df_growth_{}_count_survival_nr_drop_{}_ab_range_{}_{}.csv'.format(growth,total_drop_nr,abmin,abmax), index= None)
        return new
    ### i think for now survival is not binary; it;s a ratio of all droplets - check with Nia
    def test_surv_frac_diff_partitioning(self, partmin, partmax, step): ## one AB conc
        new = pd.DataFrame()

        # loop for partitioning factors
        for i in range(partmin, partmax, step):
             for j in range(0,2,1):
                new_i = 5**i * 2**j
                new_nr_drops_total_mass = new_i
                new_volume = variables.volume * new_nr_drops_total_mass
                total_drop_nr = round(variables.total_drop_nr /new_nr_drops_total_mass)
                strain_R = strain(new_nr_drops_total_mass)
                Droplet_exp = droplets_R(total_drop_nr, strain_R, AB_conc, new_volume)  # 0.5, 300
                Droplet_exp.run(loading, growth)
                Droplet_exp.countSurvival(growth)
                additional = pd.DataFrame({"" + str(new_nr_drops_total_mass) + "": [Droplet_exp.Res_survival_fraction]})
                # print(additional)
                new = pd.concat([new, additional], axis=1)
                #Droplet_exp.plots(growth)
                #Droplet_exp.save('initialN{}_totaldropnr{}_vol{}_loading{}_growth{}.csv'.format(initialN,new_nr_drops_total_mass,new_volume,loading,growth),
                #                 'ABconc{}_loading{}_growth{}_vol{}.csv'.format(AB_conc, loading, growth, new_volume),'Time_list.csv', AB_conc) #, Nsat, total_drop_nr, loading,growth, initialN, new_volume, growthrate, dt)
        #print("--- %s seconds ---" % (time.time() - start_time))
        pd.DataFrame(new).to_csv('output/df_growth_{}_loading_{}_ABconc{}.csv'.format(growth, loading, AB_conc), index = None)
        #print('total_drop_nr', total_drop_nr)


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
                 for j in range(0, 2, 1):
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

    ## fct to calculate survival for a set nr of repeats
    def calc_survival_prob_total_nr_bact_diff_part(self, partmin, partmax, step, spec_time, total_sim):

        ## create list with part fact as column names
        total_prob = pd.DataFrame()
        survival_fraction = []
        for nr_sim in range(0, total_sim, step):
            df_total_mass = pd.DataFrame()
            prob_diff_part = []
            part_fact = []
            ## simulate droplets for different partitioning factors
            for i in range(partmin, partmax, step):
                for j in range(0, 2, 1):
                    ## set new values for diff part factor
                    new_i = 5 ** i * 2 ** j
                    new_nr_drops_total_mass = new_i
                    new_volume = variables.volume * new_nr_drops_total_mass
                    total_drop_nr = round(variables.total_drop_nr / new_nr_drops_total_mass)

                    print('total_droplets', total_drop_nr)
                    if total_drop_nr == 0:
                        print(
                            "Error in partitioning; the partitioning factor is so small that you're trying to simulate 0 droplets")
                    else:

                        part_fct = 1 / total_drop_nr
                        part_fact.append(part_fct)

                        ## simulate growth with new params
                        strain_R = strain(new_nr_drops_total_mass)
                        Droplet_exp = droplets_R(total_drop_nr, strain_R, AB_conc, new_volume)  # 0.5, 300
                        Droplet_exp.run(loading, growth)

                        ## calculate the total nr of bacteria in all droplets, if any survived, prob survival = 1
                        Droplet_exp.countTotalMass(growth)
                        nr_bact_each_ts = Droplet_exp.total_mass
                        ## append the nr of bacteria to dataframe with N(t) vs part factor
                        df_total_mass['{}'.format(part_fct)] = nr_bact_each_ts

                        index = int(spec_time / variables.dt)

                        if nr_bact_each_ts[index] == 0:
                            prob_survival = 0
                        else:
                            prob_survival = 1

                        ## append probs to a list with diff probs for diff part factors
                        prob_diff_part.append(prob_survival)
                ## append row at the end of each simulation for each partition factor
                additional = pd.DataFrame(data=[prob_diff_part])
                total_prob = pd.concat([total_prob, additional])
            print('total prob', total_prob)
            print('part fac', part_fact)

            total_prob.columns = part_fact
            ## calculate the survival fraction and make new dataframe
            total_prob.loc['surv_Frac'] = total_prob.sum() / (len(total_prob.axes[0]))
            print(total_prob)

            surv_df = total_prob.tail(1)
            pd.DataFrame(surv_df).to_csv(
                'output/survival_fraction_growth_{}_loading_{}_ABconc{}.csv'.format(growth, loading, AB_conc),
                index=None)
            pd.DataFrame(df_total_mass).to_csv(
                'output/df_growth_{}_starting_nr_drops_{}.csv'.format(growth, variables.total_drop_nr), index=None)
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
#simulate.test_surv_frac_diff_partitioning(0, 5, 1)
#simulate.count_total_mass(part_min, part_max, step)
#simulate.test_surv_frac_diff_ab_conc(abmin, abmax, step)
simulate.calc_survival_prob_total_nr_bact_diff_part(part_min, part_max, step, spec_time, total_sim)
#simulate.count_total_mass_diff_ab(part_min, part_max, abmin, abmax, step)
#simulate.count_total_mass(abmin, abmax, step)
