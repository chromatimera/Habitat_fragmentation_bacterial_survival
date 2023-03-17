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
