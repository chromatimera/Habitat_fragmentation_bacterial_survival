import pandas as pd
import os
import matplotlib.pyplot as plt
from variables import nr_timesteps

print(os.getcwd())
os.chdir("./output/df_test_dt_nr_timesteps_{}.csv".format(nr_timesteps))
print(os.getcwd())

expo_df = pd.read_csv('df_degradation_exponential_nr_timesteps_{}.csv'.format(nr_timesteps))#, usecols=range(0, 11))
print(expo_df)

linear_df = pd.read_csv('df_degradation_linear_nr_timesteps_{}.csv'.format(nr_timesteps))#, usecols=range(1, 11))
print(linear_df)

difference_df = pd.DataFrame()
print(difference_df)

for i in range(0,len(expo_df.axes[1]),1):
    t = 1/10**(i)
    #print(t)
    difference_df[''+str(t)+''] = expo_df.iloc[:, i] - linear_df.iloc[:, i]
    # difference_df[''+str(t)+''] = expo_df.iloc[:,i] - linear_df.iloc[:,i]

column_names = list(difference_df.columns.values)
print(column_names)
print("\nDifference of expo - linear :\n", difference_df)


pd.DataFrame(difference_df).to_csv('difference_df.csv', index=None)


plt.figure(1)
#plt.grid(True)
for i in range(0,len(expo_df.axes[1]),1):
    plt.plot(difference_df.iloc[:,i], label=column_names[i])

plt.grid(True)
#plt.title('Concentration of antibiotic over time in each droplet.')
plt.ylabel('Error in rate (log scale)')
plt.xlabel('Nr of dts (timesteps)')
plt.legend()
plt.yscale('log')
#plt.ylim(bottom=0)
#plt.show()
plt.savefig('Difference between linear and exponential degradation {} steps.svg'.format(nr_timesteps))
plt.savefig('Difference between linear and exponential degradation {} steps.png'.format(nr_timesteps))