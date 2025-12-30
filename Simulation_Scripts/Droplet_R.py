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
import scipy.special as sc
getcontext().prec = 50
from mpmath import *
import matplotlib.pyplot as plt

BIGGER_SIZE = 22

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('legend', fontsize=BIGGER_SIZE)    # legend fontsize


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

        plt.rcParams['axes.spines.right'] = False
        plt.rcParams['axes.spines.top'] = False

    def run(self, init_type, grow_meth):
        #run identical experiments in each droplet
        Exp = Experiment_R(self.strain_r, self.AB_conc, self.nr_drops_total_mass)
        for k in range(0, self.total_drop_number):
            #print('drop. nr:', k)
            # if growth method is different than binary or balanced
            if (grow_meth != "binary") and (grow_meth != "balanced"):
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
        # Get current working directory
        curr_path = os.getcwd()

        # Construct folder name based on parameters
        folder_name = 'output/dropnr_{}_loading_{}_growth_{}_initialN_{}_abconc_{}_gr_{}_dt_{}_Nsat_{}'.format(
            variables.droplet_list[-1],
            loading,
            growth,
            variables.initialN,
            self.AB_conc,
            growthrate,
            self.dt,
            variables.Nsat
        )
        path = os.path.join(curr_path, folder_name)
        isExist = os.path.exists(path)

        if not isExist:
            # Create directory if it doesn't exist
            os.makedirs(path, exist_ok=True)
            print("The new directory is created!")

        os.chdir(path)

        # Plotting for non-binary growth methods
        if grow_meth != "binary" and grow_meth != "balanced":
            # Plot Number of Bacteria Over Time
            plt.figure(figsize=(7.5, 5))
            for i in range(self.total_drop_number):
                plt.plot(self.time_list_gillespie[i], self.N_list_gillespie[i], label=f'Droplet {i + 1}')
            plt.grid(False)
            plt.ylabel('Number of Bacteria')
            plt.xlabel('Time (min)')
            plt.xlim(0, self.t_end)
            # Removed: plt.ylim(bottom=0) to allow automatic y-axis scaling
            plt.legend()
            plt.savefig(f'plot_Nbact_loading_{loading}_growth_{growth}_ab_conc_{self.AB_conc}.png')
            plt.close()

            # Plot Antibiotic Concentration Over Time
            plt.figure(figsize=(7.5, 5))
            for i in range(self.total_drop_number):
                plt.plot(self.time_list_gillespie[i], self.AB_conc_array_gillespie[i], label=f'Droplet {i + 1}')
            plt.grid(False)
            plt.ylabel(r'\bf{N(t)}')
            plt.xlabel(r'\bf{Time (min)}')
            plt.xlim(0, self.t_end)
            # Removed: plt.ylim(bottom=0) to allow automatic y-axis scaling
            plt.title(f'{growth}_growth')
            plt.legend()
            plt.savefig(f'plot_ABconc_{self.AB_conc}_loading_{loading}_growth_{growth}.png')
            plt.close()

            # Plot Number of Bacteria per Volume Over Time
            plt.figure(figsize=(7.5, 5))
            for i in range(self.total_drop_number):
                plt.plot(self.time_list_gillespie[i], self.N_list_gillespie[i] / self.volume, label=f'Droplet {i + 1}')
            plt.grid(False)
            plt.ylabel(f'Number of Bacteria per {self.volume} ml Volume')
            plt.xlabel('Time (min)')
            plt.xlim(0, self.t_end)
            # Removed: plt.ylim(bottom=0) to allow automatic y-axis scaling
            plt.title(f'{growth}_growth')
            plt.legend()
            plt.savefig(f'plot_Nbact_per_volume_loading_{loading}_growth_{growth}.png')
            plt.close()

        # Plotting for binary growth method
        else:
            # Create x-axis of time in minutes
            X = np.arange(0, self.t_end, self.dt)
            XX = np.tile(X, (self.total_drop_number, 1))

            # Plot Number of Bacteria Over Time
            plt.figure(figsize=(9, 7))
            for i in range(self.total_drop_number):
                plt.plot(X, self.N_r_array[i, :], label=f'Droplet {i + 1}')
            plt.grid(False)
            plt.ylabel(r'$N(t)$ (cells)')
            plt.xlabel(r'$t$ (min)')
            plt.xlim(0, self.t_end)
            # Removed: plt.ylim(bottom=0) to allow automatic y-axis scaling

            # Commented out shading
            """
            # Find positions in N_array where population survives (6) or dies (5)
            y_survival = np.argwhere(self.N_r_array.T[0, :] == 6)
            y_death = np.argwhere(self.N_r_array.T[0, :] == 5)
            print('y_surv', y_survival)
            print('y_death', y_death)

            # Initialize variables
            y_surv = []
            y_deth = []
            col_ind_y_survival = None
            col_ind_y_death = None

            # Handle y_survival
            if y_survival.size > 0:
                col_ind_y_survival = int(y_survival[0].item())
                y_surv = list(self.N_r_array.T[:, col_ind_y_survival])
            else:
                print("No survival cases found.")
                y_surv = [0] * self.timesteps  # Example default

            # Handle y_death
            if y_death.size > 0:
                col_ind_y_death = int(y_death[0].item())
                y_deth = list(self.N_r_array.T[:, col_ind_y_death])
            else:
                print("No death cases found.")
                y_deth = [0] * self.timesteps  # Example default

            # Proceed with plotting fill_between only if both indices are found
            if col_ind_y_survival is not None and col_ind_y_death is not None:
                print('y_surv', y_surv)
                print('y_deth', y_deth)
                plt.rcParams['hatch.color'] = 'g'
                plt.yticks(range(1, 11))
                plt.fill_between(X, y_surv, 100, color='green', alpha=0.2, hatch='.', interpolate=True)
                plt.fill_between(X, 0, y_deth, color='red', alpha=0.2, interpolate=True)
            else:
                print("Skipping fill_between due to missing survival or death data.")
            """

            plt.legend()
            plt.savefig(f'plot_Nbact_loading_{loading}_growth_{growth}ab_conc_{AB_conc}_m_{self.total_drop_number}.png')
            plt.savefig(f'plot_Nbact_loading_{loading}_growth_{growth}ab_conc_{AB_conc}_m_{self.total_drop_number}.svg', format="svg", dpi=600)
            plt.close()

            # Plot Antibiotic Concentration Over Time
            plt.figure(figsize=(9, 7))
            for i in range(self.total_drop_number):
                plt.plot(X, self.AB_conc_array[i, :], label=f'Droplet {i + 1}')
            plt.grid(False)
            plt.ylabel(r'$a(t)$ ($\mu$g/ml)')
            plt.xlabel(r'$t$ (min)')
            plt.xlim(0, self.t_end)
            # Removed: plt.ylim(bottom=0) to allow automatic y-axis scaling

            # Commented out shading
            """
            # Reuse col_ind_y_survival and col_ind_y_death for AB concentration
            if col_ind_y_survival is not None:
                y_surv = list(self.AB_conc_array.T[:, col_ind_y_survival])
            else:
                y_surv = [0] * self.timesteps  # Example default

            if col_ind_y_death is not None:
                y_deth = list(self.AB_conc_array.T[:, col_ind_y_death])
            else:
                y_deth = [0] * self.timesteps  # Example default

            print('y_surv', y_surv)
            print('y_deth', y_deth)
            print('XX.T', XX.T[:, 1])

            plt.rcParams['hatch.color'] = 'g'
            plt.fill_between(X, y_deth, 100, color='red', alpha=0.2, interpolate=True)
            plt.fill_between(X, 0, y_surv, color='green', alpha=0.2, hatch='.', interpolate=True)
            """

            plt.legend()
            plt.savefig(f'plot_ABconc_{self.AB_conc}_loading_{loading}_growth_{growth}_m_{self.total_drop_number}.png')
            plt.savefig(f'plot_ABconc_{self.AB_conc}_loading_{loading}_growth_{growth}_m_{self.total_drop_number}.svg', format="svg", dpi=600)
            plt.close()

            # Plot Number of Bacteria per Volume Over Time
            plt.figure(figsize=(8, 6))
            for i in range(self.total_drop_number):
                plt.plot(X, self.N_r_array[i, :] / self.volume, label=f'Droplet {i + 1}')
            plt.grid(False)
            plt.ylabel(f'Number of Bacteria per {self.volume} ml Volume')
            plt.xlabel('Time (min)')
            plt.xlim(0, self.t_end)
            # Removed: plt.ylim(bottom=0) to allow automatic y-axis scaling
            plt.title(f'{growth}_growth')
            plt.legend()
            plt.savefig(f'plot_Nbact_per_volume_loading_{loading}_growth_{growth}_m_{self.total_drop_number}.png')
            plt.savefig(f'plot_Nbact_per_volume_loading_{loading}_growth_{growth}_m_{self.total_drop_number}.svg', format="svg", dpi=600)

            plt.close()

        os.chdir(curr_path)

    def countSurvival(self, grow_meth):
        self.last_N_list = []
        self.first_N_list = []

        if (grow_meth != "binary" and grow_meth != "balanced"):
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

        if (grow_meth != "binary" and grow_meth != "balanced"):

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

    def calc_theo_survival_prob(self, N_array):
        ##calculate rho_threshold; eq (9) from paper
        b = 1
        F1 = self.strain_r.deathrate/(b * Vmax)
        F = (self.AB_conc - self.strain_r.MIC) + Km * np.log(self.AB_conc / self.strain_r.MIC) #eq4
        rho_T = F1 * F ## **units: cell/vol if Vmax is ug/min

        ## calculate N_T

        N_T = np.floor(rho_T * self.volume)
        print('N_T',N_T)

        ## calculate the theoretical survival probability; eq. (10) from paper
        ## rhobulk*v form the paper is initialN  (deterministic) from the simulations and the first value from the N_array
        #rho_bulk = variables.initialN * variables.total_drop_nr/variables.volume * variables.total_drop_nr
        rho_bulk = variables.initialN / variables.volume # constant in det
        #print('rho bulk', rho_bulk)
        #print('vol',self.volume)

        exp_fact = math.exp(-rho_bulk*self.volume)
        ps = exp_fact * nsum(lambda n: (rho_bulk*self.volume)**n/fac(n), [N_T+1, inf])
        bigPs= 1- (1-ps)**self.total_drop_number

        return rho_T, N_T, ps, bigPs



    def save(self,NRfilename, ABfilename,Timefilename, AB_conc):

        #Get path of current folder
        curr_path = os.getcwd()
        #print(curr_path)

        # Check whether the specified path exists or not

        folder_name = 'output/dropnr_{}_loading_{}_growth_{}_initialN_{}_abconc_{}_gr_{}'.format(variables.droplet_list[-1], variables.loading, variables.growth, variables.initialN, variables.AB_conc, variables.growthrate)
        path = os.path.join(curr_path, folder_name)
        isExist = os.path.exists(path)

        if not isExist:
            # Create a new directory because it does not exist
            os.makedirs(path, exist_ok = True)
            print("The new directory is created!")

        os.chdir(path)

        if (growth != "binary" and growth != "balanced"):
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
