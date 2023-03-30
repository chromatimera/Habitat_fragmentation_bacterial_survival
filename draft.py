import math
### nr of droplets and the power of i and j in the partitioning loop are related i.e. for 100 droplets, the loop goes from 0 to 3.
total_drop_nr = 1000
## don't change the 2 lines below
part_min = 0
part_max = math.floor(math.log(total_drop_nr, 10)) + 1
print(part_max)
step = 1
for i in range(part_min, part_max, step):
    # print('i', i)
    for j in range(part_min, part_max, step):
        # print('j', j)

        ## set new values for diff part factor
        new_i = 5 ** i * 2 ** j
        print('new_i', new_i)


part_max = int(math.log10(1000) + 1)
print(part_max)