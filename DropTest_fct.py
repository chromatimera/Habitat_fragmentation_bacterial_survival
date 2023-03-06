from strain import strain
from Droplet_R import droplets_R
import matplotlib.pyplot as plt
import numpy as np
from variables import *
from decimal import *
import pandas as pd


getcontext().prec = 50


import time
start_time = time.time()


class DropTest(object):
    def run(self):

        strain_R = strain(initialN, growthrate, deathrate, MIC, dt, t_end, Nsat) #prob_ext=None, time_ext=None)
        #2 Just resistant strain exp: def __init__(self, total_drop_number, strain_r, AB_conc, dt, t_end):
        Droplet_exp = droplets_R(total_drop_nr, strain_R, AB_conc, dt, t_end)  #0.5, 300
        #3
        Droplet_exp.run(loading, growth)
        #4

        Droplet_exp.save('initialN{}'
                         '_growthrate{}_MIC{}_totaldropnr{}_ABconc{}_'
                         'dt{}_loading{}_growth{}.csv'.format(initialN, growthrate, MIC, total_drop_nr, AB_conc, dt,
                                                              loading, growth), 'ABconc{}_loading{}_growth{}.csv'.format(AB_conc, loading, growth),'Time_list.csv', Nsat,
                                                              total_drop_nr, loading, growth, initialN, AB_conc, growthrate, dt)



        print("--- %s seconds ---" % (time.time() - start_time))
        Droplet_exp.plots(growth)

        Droplet_exp.countSurvival(growth)
        #Droplet_exp.normed_histograms(50)
    def test_dt(self, dtmin, dtmax, stepdt):
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
            # Droplet_exp.countSurvival()
            # Droplet_exp.normed_histograms(50)

        # Get path of current folder
        curr_path = os.getcwd()
        # print(curr_path)

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

            # print(strain_R.initialN)
            # print(strain_R.MIC)
            # print(Droplet_exp.AB_conc_molecules)

            print("--- %s seconds ---" % (time.time() - start_time))
            #Droplet_exp.plots(growth)
            # Droplet_exp.normed_histograms(50)
            #print(new)
            pd.DataFrame(new).to_csv('output/df_growth_{}_count_survival_nr_drop_{}_ab_range_{}_{}.csv'.format(growth,total_drop_nr,abmin,abmax), index= None)
        return new

    def count_total_mass(self, abmin, abmax, step):
        df_total_mass = pd.DataFrame()
        # loop for variables
        for i in range(abmin, abmax, step):
            new_ab_conc = i
            print('new ab conc', new_ab_conc)
            strain_R = strain(initialN, growthrate, deathrate, MIC, dt, t_end, Nsat)
            Droplet_exp = droplets_R(total_drop_nr, strain_R, new_ab_conc, dt, t_end)  # 0.5, 300
            Droplet_exp.run(loading, growth)
            # Droplet_exp.save(
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

            # Droplet_exp.plots(growth)
            ###count survival if necessary

            # Droplet_exp.countSurvival(growth)

            # additional = pd.DataFrame({"" + str(new_ab_conc) + "": [Droplet_exp.total_mass]})
            # print(additional)
            # new = pd.concat([new, additional], axis=1)


            Droplet_exp.countTotalMass(growth)
            additional = Droplet_exp.total_mass
            df_total_mass ['Total_mass, ab conc {}'.format(new_ab_conc)] = additional

        df_total_mass.insert(loc=0, column='Time', value=Droplet_exp.time_list_array)

        # pd.DataFrame(new).to_csv('output/df_growth_{}_count_survival.csv'.format(growth), index=None)
        pd.DataFrame(df_total_mass).to_csv('output/df_growth_{}_nr_drops_{}_big_droplet_nr_{}.csv'.format(growth, total_drop_nr, nr_drops_total_mass), index= None)
        print("--- %s seconds ---" % (time.time() - start_time))

    def test_surv_frac_convergence(self, dropmin, dropmax, nr_points):
        new = pd.DataFrame()
        drop_list = np.linspace(dropmin, dropmax, nr_points)
        # loop for variables
        for new_total_drop_nr in drop_list:
            new_total_drop_nr = round(new_total_drop_nr)
            print(new_total_drop_nr)
            strain_R = strain(initialN, growthrate, deathrate, MIC, dt, t_end, Nsat)
            Droplet_exp = droplets_R(new_total_drop_nr, strain_R, AB_conc, dt, t_end)  # 0.5, 300
            Droplet_exp.run(loading, growth)
            Droplet_exp.countSurvival(growth)
            #surv_frac = Droplet_exp.Res_survival_fraction
            #print(surv_frac)
            additional = pd.DataFrame({"" + str(new_total_drop_nr) + "": [Droplet_exp.Res_survival_fraction]})
            # print(additional)
            new = pd.concat([new, additional], axis=1)
            Droplet_exp.save(
                'initialN{}_growthrate{}_MIC{}_totaldropnr{}_ABconc{}_dt{}_loading{}_growth{}.csv'.format(initialN,
                                                                                                          growthrate,
                                                                                                          MIC,
                                                                                                          new_total_drop_nr,
                                                                                                          AB_conc,
                                                                                                          dt,
                                                                                                          loading,
                                                                                                          growth),
                'ABconc{}_loading{}_growth{}.csv'.format(AB_conc, loading, growth), 'Time_list.csv', Nsat, new_total_drop_nr, loading,
                growth, initialN, AB_conc, growthrate, dt)

            # print(strain_R.initialN)
            # print(strain_R.MIC)
            # print(Droplet_exp.AB_conc_molecules)

            # print("--- %s seconds ---" % (time.time() - start_time))
            #Droplet_exp.plots(growth)
            # Droplet_exp.normed_histograms(50)
            print(new)
            pd.DataFrame(new).to_csv('df_growth_{}_total_drop_nr_changing.csv'.format(growth), index= None)
        return new

    def test_ab(self):
        for i in range(abmin, abmax, step):
            strain_R = strain(initialN, growthrate, deathrate, MIC, dt, t_end, Nsat)

            Droplet_exp = droplets_R(total_drop_nr, strain_R, i, dt, t_end)  # 0.5, 300
            Droplet_exp.run(loading, growth)

            Droplet_exp.save(
                'initialN{}_growthrate{}_MIC{}_totaldropnr{}_ABconc{}_dt{}_loading{}_growth{}.csv'.format(initialN,
                                                                                                          growthrate,
                                                                                                          MIC,
                                                                                                          total_drop_nr,
                                                                                                          i, dt,
                                                                                                          loading,
                                                                                                          growth),
                'ABconc{}_loading{}_growth{}.csv'.format(i, loading, growth), 'Time_list.csv', Nsat, total_drop_nr, loading, growth,
                initialN, i, growthrate, dt)

            # print(strain_R.initialN)
            # print(strain_R.MIC)
            # print(Droplet_exp.AB_conc_molecules)

            # print("--- %s seconds ---" % (time.time() - start_time))
            #Droplet_exp.plots(growth)
            # Droplet_exp.countSurvival()
            # Droplet_exp.normed_histograms(50)

simulate = DropTest()
simulate.run()
#simulate.test_dt(1, 2, 1)
#simulate.test_surv_frac_diff_ab_conc(abmin, abmax, step)
#simulate.count_total_mass(abmin, abmax, step)
