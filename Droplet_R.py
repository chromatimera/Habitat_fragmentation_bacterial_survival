import variables
from strain import strain
from Experiment_R import Experiment_R
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os
import pandas as pd
import shutil
from variables import *
from strain import strain_script_name, strain_script_path
from Experiment_R import experiment_script_name, experiment_script_path
from decimal import *

getcontext().prec = 50


class droplets_R():

    def __init__(self, total_drop_nr, strain_r, AB_conc, volume, nr_drops_total_mass):
        self.total_drop_number = total_drop_nr
        self.strain_r = strain_r
        self.AB_conc = AB_conc
        self.dt = variables.dt
        self.t_end = variables.t_end
        self.timesteps = round(self.t_end / variables.dt)
        self.volume = volume
        self.nr_drops_total_mass = nr_drops_total_mass

        ## if doing deterministic growth
        self.N_r_array = np.empty((self.total_drop_number, self.timesteps))#.astype(int) # initialize empty array
        self.AB_conc_array = np.empty((self.total_drop_number, self.timesteps))  # initialize empty array
        self.time_list_array = np.empty((self.total_drop_number, self.timesteps))

        #if doing stochastic growth
        self.N_list_gillespie = [] # N_r_array for gillespie since the nr of rows are different between droplets, and cannot concatenate
        self.time_list_gillespie = []
        self.AB_conc_array_gillespie = []

        plt.rc('font', size=14)  # controls default text size
        plt.rcParams['axes.spines.right'] = False
        plt.rcParams['axes.spines.top'] = False

    def run(self, init_type, grow_meth):
        #run identical experiments in each droplet
        Exp = Experiment_R(self.strain_r, self.AB_conc, self.nr_drops_total_mass)

        for k in range(0, self.total_drop_number):
            #print('drop. nr:', k)
            if (grow_meth != "binary"):
                Exp.run(init_type, grow_meth)
                ## not necessarily gillespie, but the point is that the N, AB_conc and Time list have variable lengths
                self.N_list_gillespie.append(Exp.N_array)
                self.time_list_gillespie.append(Exp.ts)
                self.AB_conc_array_gillespie.append(Exp.AB_conc_array)
            else:
                Exp.run(init_type, grow_meth)
                self.N_r_array[k] = Exp.N_array
                self.time_list_array[k] = Exp.ts
                self.AB_conc_array[k] = Exp.AB_conc_array


    def plots(self, grow_meth):

        # Get path of current folder
        curr_path = os.getcwd()
        #print(curr_path)

        # Check whether the specified path exists or not


        folder_name = 'output/' \
                      'dropnr_{}_loading_{}_growth_{}_initialN_{}_abconc_{}_gr_{}_dt_{}_Nsat_{}'.format(variables.total_drop_nr,loading,growth, variables.initialN, self.AB_conc, growthrate, self.dt,variables.Nsat)
        path = os.path.join(curr_path, folder_name)
        isExist = os.path.exists(path)

        if not isExist:
            # Create a new directory because it does not exist
            os.makedirs(path, exist_ok=True)
            print("The new directory is created!")

        os.chdir(path)

        # create x axis of time in minutes:
        if (grow_meth != "binary"):

            plt.rcParams.update({'font.size': 14})

            fig, ax = plt.subplots()
            for i in range(0, self.total_drop_number):
                plt.plot(self.time_list_gillespie[i], self.N_list_gillespie[i])
            plt.grid(True)
            # plt.title('Growth of resistant strain')
            plt.title('{}_growth'.format(growth))
            plt.ylabel('Number of bacteria')
            plt.xlabel('Time (mins)')
            plt.xlim(0, self.t_end)
            plt.xlim(0, self.t_end)
            # tick_spacing = 2
            # ax.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
            plt.ylim(bottom=0)
            plt.savefig('plot_Nbact_loading_{}_growth_{}ab_conc_{}.png'.format(loading, growth, str(AB_conc)))


            plt.figure(2)
            for i in range(0, self.total_drop_number):
                plt.plot(self.time_list_gillespie[i], self.AB_conc_array_gillespie[i])
            plt.grid(True)
            # plt.title('Concentration of antibiotic over time in each droplet.')
            plt.ylabel('Antibiotic (ug/nL)')
            plt.xlabel('Time (mins)')
            plt.xlim(0, self.t_end)
            plt.ylim(bottom=0)
            plt.title('{}_growth'.format(growth))
            plt.savefig('plot_ABconc_{}_loading_{}_growth_{}.png'.format(str(AB_conc), loading, growth))
            plt.show()

            plt.figure(3)
            for i in range(0, self.total_drop_number):
                plt.plot(self.time_list_gillespie[i], self.N_list_gillespie[i]/self.volume)
            plt.grid(True)
            plt.ylabel('N/V')
            plt.xlabel('Time (min)')
            plt.xlim(0, self.t_end)
            plt.ylim(bottom=0)
            plt.title('{}_growth'.format(growth))
            plt.show()

        else:
            X = np.arange(0, self.t_end, self.dt).tolist()
            XX= np.tile(X, (self.total_drop_number, 1))

            plt.rcParams.update({'font.size': 14})

            fig, ax= plt.subplots()
            #plt.grid(True)
            plt.plot (XX.T, self.N_r_array.T)
            plt.grid(True)
            #plt.title('Growth of resistant strain')
            plt.title ('{}_growth'.format(growth))
            plt.ylabel('Number of bacteria')
            plt.xlabel('Time (mins)')
            plt.xlim(0,self.t_end)
            plt.xlim(0,self.t_end)
            #tick_spacing = 2
            #ax.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
            plt.ylim(bottom=0)
            plt.savefig('plot_Nbact_loading_{}_growth_{}ab_conc_{}'.format(loading, growth, AB_conc))


            plt.figure(2)
            #plt.grid(True)
            plt.plot(XX.T, self.AB_conc_array.T)
            plt.grid(True)
            #plt.title('Concentration of antibiotic over time in each droplet.')
            plt.ylabel('Antibiotic (ug/mL)')
            plt.xlabel('Time (mins)')
            plt.xlim(0, self.t_end)
            plt.ylim(bottom=0)
            plt.title ('{}_growth'.format(growth))
            plt.savefig('plot_ABconc_{}_loading_{}_growth_{}'.format(self.AB_conc, loading, growth))
            plt.show()

            plt.figure(3)
            plt.plot(XX.T, self.N_r_array.T / self.volume)
            plt.grid(True)
            plt.ylabel('N/V')
            plt.xlabel('Time (min)')
            plt.xlim(0, self.t_end)
            plt.ylim(bottom=0)
            plt.title('{}_growth'.format(growth))
            plt.show()

        os.chdir(curr_path)

    def countSurvival(self, grow_meth):
        self.last_N_list = []
        self.first_N_list = []

        if (grow_meth != "binary"):
            for i in range(0, len(self.N_list_gillespie)):
                self.last_N_list.append(self.N_list_gillespie[i][-1])
                self.first_N_list.append(self.N_list_gillespie[i][0])
            #print('N_list', self.last_N_list)
            Res_living = np.count_nonzero(self.last_N_list)
            if np.count_nonzero(self.first_N_list) == 0:
                self.Res_survival_fraction = 0
            else:
                self.Res_survival_fraction = Res_living / np.count_nonzero(self.first_N_list)
                print("Fraction of droplets where bacteria survived=  " + str(self.Res_survival_fraction))

        else:
            Res_living = np.count_nonzero(self.N_r_array[:, self.timesteps - 1])
            if np.count_nonzero(self.N_r_array[:, 0]) == 0:
                self.Res_survival_fraction = 0
            else:
                self.Res_survival_fraction = Res_living / np.count_nonzero(self.N_r_array[:, 0])
            print("Fraction of droplets where bacteria survived=  " + str(self.Res_survival_fraction))

    def countTotalMass(self, grow_meth):
        self.total_mass = np.empty(self.timesteps)

        if (grow_meth != "binary"):

            ## for all times with increment dt from 0 until t_end
            for t in range(0, self.timesteps):
                self.mass_droplet = []
                ## for each droplet
                for i in range(0, len(self.N_list_gillespie)):
                    time_list = sorted(j for j in self.time_list_gillespie[i] if j <= t * self.dt)
                    self.mass_droplet.append(int(self.N_list_gillespie[i][len(time_list)-1]))
                self.total_mass[t] = np.sum(self.mass_droplet)
        else:
            for i in range(0, self.timesteps):
                if np.count_nonzero(self.N_r_array[:, 0]) == 0:
                    self.total_mass[i] = 0
                else:
                    self.total_mass[i] = np.sum(self.N_r_array[:, i])
                #print(self.total_mass[i])

            #print('total_mass', self.total_mass)
            #print("Total mass in droplets where bacteria survived=  " + str(self.total_mass))

    def calc_tau_det(self, N_array):  # calculates tau (paper eq 7) for the det case
        b = 1
        F = (self.AB_conc - self.strain_r.MIC) + Km * np.log(self.AB_conc / self.strain_r.MIC)
        B = self.strain_r.deathrate * self.volume / (N_array[0,0] * Vmax * b)  # do we need to add a value for b in the simulations
        C = np.log(1 - B * F)
        #print(F,B,C)
        tau = (-1 / self.strain_r.deathrate) * C
        print('Tau= ', tau, 'min')
        return tau



    def save(self,NRfilename, ABfilename,Timefilename, AB_conc):

        #Get path of current folder
        curr_path = os.getcwd()
        #print(curr_path)

        # Check whether the specified path exists or not

        folder_name = 'output/dropnr_{}_loading_{}_growth_{}_initialN_{}_abconc_{}_gr_{}_dt_{}_Nsat_{}'.format(variables.total_drop_nr, variables.loading, variables.growth, variables.initialN, AB_conc, variables.growthrate, self.dt,variables.Nsat)
        path = os.path.join(curr_path, folder_name)
        isExist = os.path.exists(path)

        if not isExist:
            # Create a new directory because it does not exist
            os.makedirs(path, exist_ok = True)
            print("The new directory is created!")

        os.chdir(path)

        if (growth != "binary"):
            new_fn=os.path.join(path+NRfilename)   ##os.path.join(path+"N_R.csv")

            pd.DataFrame(self.AB_conc_array_gillespie).to_csv(ABfilename)
            pd.DataFrame(self.time_list_gillespie).to_csv(Timefilename)
            pd.DataFrame(self.N_list_gillespie).to_csv(NRfilename)
        else:
            pd.DataFrame(self.N_r_array).to_csv(NRfilename)
            pd.DataFrame(self.AB_conc_array).to_csv(ABfilename)


        # make dir
        scripts_dir = path + '/scripts'  # data_dir is the path to dir in which data is in
        if not os.path.exists(scripts_dir):
            os.mkdir(scripts_dir)

        # get script path
        current_script_path = __file__
        #print(current_script_path)
        current_script_name = os.path.basename(__file__)
        #print(current_script_name)

        shutil.copyfile(current_script_path, scripts_dir + '/' + current_script_name)
        shutil.copyfile(variables_script_path, scripts_dir+ '/' + variables_script_name)
        shutil.copyfile(strain_script_path, scripts_dir+ '/' + strain_script_name)
        shutil.copyfile(experiment_script_path, scripts_dir + '/' + experiment_script_name)

        os.chdir(curr_path)

    def save_dt(self, dt_filename, growth, initialN, AB_conc):

        #Get path of current folder
        curr_path = os.getcwd()
        #print(curr_path)

        # Check whether the specified path exists or not

        folder_name = 'output/df_degradation_{}_nr_timesteps_{}.csv'.format(degradation, nr_timesteps)
        path = os.path.join(curr_path, folder_name)
        isExist = os.path.exists(path)

        if not isExist:
            # Create a new directory because it does not exist
            os.makedirs(path, exist_ok = True)
            print("The new directory is created!")

        os.chdir(path)
        pd.DataFrame(self.AB_conc_array).to_csv(ABfilename)



        # make dir
        scripts_dir = path + '/scripts'  # data_dir is the path to dir in which data is in
        if not os.path.exists(scripts_dir):
            os.mkdir(scripts_dir)

        # get script path
        current_script_path = __file__
        #print(current_script_path)
        current_script_name = os.path.basename(__file__)
        #print(current_script_name)

        shutil.copyfile(current_script_path, scripts_dir + '/' + current_script_name)
        shutil.copyfile(variables_script_path, scripts_dir+ '/' + variables_script_name)
        shutil.copyfile(strain_script_path, scripts_dir+ '/' + strain_script_name)
        shutil.copyfile(experiment_script_path, scripts_dir + '/' + experiment_script_name)

        os.chdir(curr_path)
