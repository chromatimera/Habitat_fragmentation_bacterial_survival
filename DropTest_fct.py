import variables
from strain import strain
from Droplet_R import droplets_R
from variables import *
from decimal import *
import pandas as pd


getcontext().prec = 50


import time
start_time = time.time()


class DropTest(object):
    def run(self): ## simulation for one AB concentration; used in troubleshooting

        strain_R = strain(nr_drops_total_mass=1) #prob_ext=None, time_ext=None)
        #2 Just resistant strain exp: def __init__(self, total_drop_number, strain_r, AB_conc, dt, t_end):
        Droplet_exp = droplets_R(strain_R, AB_conc)  #0.5, 300
        #3
        Droplet_exp.run(loading, growth)
        #4

        Droplet_exp.save('initialN{}'
                         '_growthrate{}_MIC{}_totaldropnr{}_ABconc{}_'
                         'dt{}_loading{}_growth{}.csv'.format(initialN, growthrate, MIC, total_drop_nr, AB_conc, dt,
                                                              loading, growth), 'ABconc{}_loading{}_growth{}.csv'.format(AB_conc, loading, growth),'Time_list.csv', AB_conc)



        print("--- %s seconds ---" % (time.time() - start_time))
        Droplet_exp.plots(growth)

        Droplet_exp.countSurvival(growth)
        #Droplet_exp.normed_histograms(50)
    def test_dt(self, dtmin, dtmax, stepdt): ## will be removed; fct to test difference in AB degradation between different dts - can estimate error in linear approximation
        new = pd.DataFrame()

        # loop for variables
        for i in range(dtmin, dtmax, stepdt):
            new_dt = 1 / (10 ** i)
            new_t_end = new_dt * nr_timesteps
            # new_t_end = 300
            print('dt = ', new_dt, 't_end = ', new_t_end)
            #print(new_t_end)
            strain_R = strain(initialN, growthrate, deathrate, MIC, new_dt, new_t_end, Nsat)
            Droplet_exp = droplets_R(total_drop_nr, strain_R, AB_conc, new_dt, new_t_end)  # 0.5, 300
            Droplet_exp.run(loading, growth)
            additional = pd.DataFrame({"" + str(new_dt) + "": Droplet_exp.deg_list})
            new = pd.concat([new, additional], axis=1)
            Droplet_exp.save(
                'initialN{}_growthrate{}_MIC{}_totaldropnr{}_ABconc{}_dt{}_loading{}_growth{}.csv'.format(initialN,
                                                                                                          growthrate,
                                                                                                          MIC,
                                                                                                          total_drop_nr,
                                                                                                          AB_conc,
                                                                                                          new_dt,
                                                                                                          loading,
                                                                                                          growth),
                'ABconc{}_loading{}_growth{}.csv'.format(AB_conc, loading, growth), Nsat, total_drop_nr, loading,
                growth, initialN, AB_conc, growthrate, new_dt)



            # print("--- %s seconds ---" % (time.time() - start_time))
            Droplet_exp.plots(growth)

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

        # loop for variables
        for i in range(abmin, abmax, step):
            new_ab_conc = i
            print('antib_conc', new_ab_conc)
            strain_R = strain(initialN, growthrate, deathrate, MIC, dt, t_end, Nsat)
            Droplet_exp = droplets_R(total_drop_nr, strain_R, new_ab_conc, dt, t_end)  # 0.5, 300
            Droplet_exp.run(loading, growth)
            Droplet_exp.countSurvival(growth)
            #surv_frac = Droplet_exp.Res_survival_fraction
            #print(surv_frac)
            additional = pd.DataFrame({"" + str(new_ab_conc) + "": [Droplet_exp.Res_survival_fraction]})
            # print(additional)
            new = pd.concat([new, additional], axis=1)
            # #Droplet_exp.save(
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
            #Droplet_exp.plots(growth)
            # Droplet_exp.normed_histograms(50)
            #print(new)
            pd.DataFrame(new).to_csv('output/df_growth_{}_count_survival_nr_drop_{}_ab_range_{}_{}.csv'.format(growth,total_drop_nr,abmin,abmax), index= None)
        return new

    ## Test the survival fraction of droplets for different factors, ab conc and calculate the total mass
    def test_survival_big_droplet_diff_ab_conc(self, abmin, abmax, step, n_drops_min, n_drops_max, step_drops):

        # loop for diff antibiotic conc
        for fakedrop in range(n_drops_min, n_drops_max, step_drops):
            print(fakedrop)
            new_total_drop_nr = int(variables.total_drop_nr/fakedrop)
            print('total_drop_nr', new_total_drop_nr)
            nr_drops_total_mass = fakedrop
            print('new_total-drop_nr',new_total_drop_nr)
            variables.Nsat = variables.Nsat * nr_drops_total_mass
            variables.volume = variables.volume * nr_drops_total_mass
            variables.initialN = variables.initialN * nr_drops_total_mass

            df_total_mass = pd.DataFrame()
            surv_frac_diff_ab = pd.DataFrame()

            for i in range(abmin, abmax, step):
                new_ab_conc = i
                #print('antib_conc', new_ab_conc)
                strain_R = strain(variables.initialN, growthrate, deathrate, MIC, dt, t_end, variables.Nsat, nr_drops_total_mass)
                Droplet_exp = droplets_R(strain_R, new_ab_conc, dt, t_end)  # 0.5, 300
                Droplet_exp.run(loading, growth)
                Droplet_exp.countSurvival(growth)
                additional = pd.DataFrame({"" + str(new_ab_conc) + "": [Droplet_exp.Res_survival_fraction]})
                # print(additional)
                surv_frac_diff_ab = pd.concat([surv_frac_diff_ab, additional], axis=1)
                pd.DataFrame(surv_frac_diff_ab).to_csv('output/df_growth_{}_count_survival_nr_drop_{}_ab_range_{}_{}.csv'.format(growth,new_total_drop_nr,abmin,abmax), index= None)

                Droplet_exp.countTotalMass(growth)
                additional_ab = Droplet_exp.total_mass
                df_total_mass['Total_mass, ab conc {}'.format(new_ab_conc)] = additional_ab

            df_total_mass.insert(loc=0, column='Time', value=Droplet_exp.time_list_array)

            # pd.DataFrame(new).to_csv('output/df_growth_{}_count_survival.csv'.format(growth), index=None)
            pd.DataFrame(df_total_mass).to_csv('output/df_growth_{}_nr_drops_{}_big_droplet_nr_{}.csv'.format(growth, total_drop_nr,nr_drops_total_mass), index=None)


    def count_total_mass(self, abmin, abmax, step):
        df_total_mass = pd.DataFrame()
        # loop for variables
        for i in range(abmin, abmax, step):
            new_ab_conc = i
            print('new ab conc', new_ab_conc)
            strain_R = strain(initialN, growthrate, deathrate, MIC, dt, t_end, Nsat)
            Droplet_exp = droplets_R(total_drop_nr, strain_R, new_ab_conc, dt, t_end)  # 0.5, 300
            Droplet_exp.run(loading, growth)


            Droplet_exp.countTotalMass(growth)
            additional = Droplet_exp.total_mass
            df_total_mass ['Total_mass, ab conc {}'.format(new_ab_conc)] = additional

        df_total_mass.insert(loc=0, column='Time', value=Droplet_exp.time_list_array)

        # pd.DataFrame(new).to_csv('output/df_growth_{}_count_survival.csv'.format(growth), index=None)
        pd.DataFrame(df_total_mass).to_csv('output/df_growth_{}_nr_drops_{}_big_droplet_nr_{}.csv'.format(growth, total_drop_nr, nr_drops_total_mass), index= None)
        print("--- %s seconds ---" % (time.time() - start_time))

simulate = DropTest()
simulate.run()
#simulate.test_dt(1, 2, 1)
#simulate.test_surv_frac_diff_ab_conc(abmin, abmax, step)
#simulate.test_survival_big_droplet_diff_ab_conc(abmin,abmax,step,nr_drop_min,nr_drop_max,step_drop)
#simulate.count_total_mass(abmin, abmax, step)
