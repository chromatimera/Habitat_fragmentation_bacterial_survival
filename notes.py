from variables import *
for i in range(part_min, part_max, step):
    for j in range(part_min, part_max, step):
        new_i = 5 ** i * 2 ** j
        new_nr_drops_total_mass = new_i
        new_volume = volume * new_nr_drops_total_mass
        new_total_drop_nr = round(total_drop_nr / new_nr_drops_total_mass)
        print('new_nr_drops', new_nr_drops_total_mass)
        print('new_volume', new_volume)
        print('drop nr', new_total_drop_nr)

for i in range(len(droplet_list)):
    total_drop_nr = droplet_list[i]
    new_nr_drops_total_mass = total_drop_nr
    new_volume = volume * droplet_list[-1] / droplet_list[i]
    print('new_nr_drop', new_nr_drops_total_mass)


