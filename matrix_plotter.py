import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from populate_transition_matrix import transmission_matrix_generator

df = pd.read_csv('C:\Uganda\Malariatherapy_2017_06_29.csv', index_col=0)
df.Asexual = pd.to_numeric(df.Asexual, errors='coerce')
df = df[df.Treated==0]





mean_max = np.nanmean(pd.to_numeric(df.gam_max,errors='coerce'))

fig,axarr = plt.subplots(1,5,sharey=True)

[TM1,mask1] = transmission_matrix_generator(df,t_start=22,t_end= 200,interval= 10)
[TM2,mask2] = transmission_matrix_generator(df,t_start=22,t_end= 200,interval= 2)
[TM3,mask3] = transmission_matrix_generator(df,t_start=22,t_end= 200,interval= 3)
[TM4,mask4] = transmission_matrix_generator(df,t_start=22,t_end= 200,interval= 5)
[TM5,mask5] = transmission_matrix_generator(df,t_start=22,t_end= 200,interval= 10)
[TM6,mask6] = transmission_matrix_generator(df,t_start=22,t_end= 200,interval= 15)
[TM7,mask7] = transmission_matrix_generator(df,t_start=22,t_end= 200,interval= 25)

TM1.to_csv(r'C:\Uganda\Transition matrices\Updated\TM_gam_22_200_dt1.csv')
TM2.to_csv(r'C:\Uganda\Transition matrices\Updated\TM_gam_22_200_dt2.csv')
TM3.to_csv(r'C:\Uganda\Transition matrices\Updated\TM_gam_22_200_dt3.csv')
TM4.to_csv(r'C:\Uganda\Transition matrices\Updated\TM_gam_22_200_dt5.csv')
TM5.to_csv(r'C:\Uganda\Transition matrices\Updated\TM_gam_22_200_dt10.csv')
TM6.to_csv(r'C:\Uganda\Transition matrices\Updated\TM_gam_22_200_dt15.csv')
TM7.to_csv(r'C:\Uganda\Transition matrices\Updated\TM_gam_22_200_dt25.csv')




sns.heatmap(TM1,mask =mask1,ax = axarr[0],cmap='YlOrRd',cbar=False)
sns.heatmap(TM2,mask =mask2,ax = axarr[1],cmap='YlOrRd',cbar=False)
sns.heatmap(TM3,mask =mask3,ax = axarr[2],cmap='YlOrRd',cbar=False)
sns.heatmap(TM4,mask =mask3,ax = axarr[3],cmap='YlOrRd',cbar=False)
sns.heatmap(TM5,mask =mask4,ax = axarr[4],cmap='YlOrRd',cbar=False)


# TM5 = sns.heatmap(transmission_matrix_generator(df[df.Treated==1], 120, 1),ax=axarr[1,0],cmap = 'YlOrRd',cbar=False)
# TM6 = sns.heatmap(transmission_matrix_generator(df[df.Treated==1], 120, 5),ax=axarr[1,1],cmap = 'YlOrRd',cbar=False)
# TM7 = sns.heatmap(transmission_matrix_generator(df[df.Treated==1], 120, 10),ax=axarr[1,2],cmap = 'YlOrRd',cbar=False)
# TM8 = sns.heatmap(transmission_matrix_generator(df[df.Treated==1], 120, 15),ax=axarr[1,3],cmap = 'YlOrRd')

for ax in axarr:
    for tick in ax.get_xticklabels():
        tick.set_horizontalalignment('right')
        tick.set_rotation(45)


for ax in axarr:
    for tick in ax.get_yticklabels():
        tick.set_verticalalignment('top')
        tick.set_rotation(45)

(axarr[0]).plot([0,7],[0,7],color = 'k',linestyle ='--',linewidth=2)
(axarr[1]).plot([0,7],[0,7],color = 'k',linestyle ='--',linewidth=2)
(axarr[2]).plot([0,7],[0,7],color = 'k',linestyle ='--',linewidth=2)
(axarr[3]).plot([0,7],[0,7],color = 'k',linestyle ='--',linewidth=2)
(axarr[4]).plot([0,7],[0,7],color = 'k',linestyle ='--',linewidth=2)


fig.subplots_adjust(wspace = 0.05)
plt.show()
