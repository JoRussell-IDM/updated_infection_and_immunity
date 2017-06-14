import numpy as np
import pandas as pd

#populate a matrix of stratified density (orders of mag from <10^1 to >10^6 ) at t0 to density at t0+dt
#explore t, dt, and denisty binning

#malariatherapy
#conssiderations:
#infections all start from one definitive infection time
#micoscopy negative infections are a hidden state (sub detection threshold)

df = pd.read_csv('C:\Uganda\Malariatherapy_2017_05_31.csv', index_col=0)
df.Asexual = pd.to_numeric(df.Asexual,errors='coerce')

def transmission_matrix_generator(df,t_start,t_end,interval):
    #define the time t
    t_cutoff_pre = t_start
    t_cutoff_post = t_end
    #define the transition time dt
    dt = interval

    #define the bins for stratification of parasite density
    density_strata = ['truezero','submicroscopic','10^1-10^2','10^2-10^3','10^3-10^4','10^4-10^5','10^5-10^6']
    transition_matrix = pd.DataFrame(np.zeros((7,7)),index = density_strata, columns = density_strata)


    for patient, measurements in df.groupby('id'):
        df_test = measurements.copy().reset_index(drop=True)

        for i in df_test.index.values:
            current_strata = []
            future_strata = []
            if np.logical_and((df_test['Day'][i] < t_cutoff_post),(df_test['Day'][i] > t_cutoff_pre)):
                # current_strata = []
                # end = None
                if np.nansum(df_test.Asexual[i:]) == 0:
                    current_strata = 'truezero'
                elif df_test.Asexual[i] == 0:
                    current_strata = 'submicroscopic'
                elif df_test.Asexual[i] <= 100:
                    current_strata = '10^1-10^2'
                elif df_test.Asexual[i]<= 1000:
                    current_strata = '10^2-10^3'
                elif df_test.Asexual[i] <= 10000:
                    current_strata = '10^3-10^4'
                elif df_test.Asexual[i]<= 100000:
                    current_strata = '10^4-10^5'
                elif df_test.Asexual[i] <= 1000000:
                    current_strata = '10^5-10^6'

            if df_test['Day'][i] < (max(df_test['Day'])- dt):
                future_strata = []
                future_density = df_test.Asexual[i+dt]

                if np.nansum(df_test.Asexual[i+dt:]) == 0:
                    future_strata = 'truezero'
                elif future_density == 0:
                    future_strata = 'submicroscopic'
                elif future_density <= 100:
                    future_strata = '10^1-10^2'
                elif future_density <= 1000:
                    future_strata = '10^2-10^3'
                elif future_density <= 10000:
                    future_strata = '10^3-10^4'
                elif future_density <= 100000:
                    future_strata = '10^4-10^5'
                elif future_density <= 1000000:
                    future_strata = '10^5-10^6'

            if (np.isnan(bool(future_strata)) | np.isnan(bool(current_strata))):
                pass
            else:
                transition_matrix.loc[future_strata,current_strata] += 1


    for element in density_strata:
        transition_matrix[element] = transition_matrix[element]/(sum(transition_matrix[element])).astype(float)

    transition_matrix[transition_matrix==0] = 'NaN'

    mask = transition_matrix.isnull()

    return [transition_matrix, mask]

if __name__ == '__main__':
    df = pd.read_csv('C:\Uganda\Malariatherapy_2017_05_31.csv', index_col=0)
    df.Asexual = pd.to_numeric(df.Asexual, errors='coerce')
    # t = 120
    # interval = 20
    # TM = transmission_matrix_generator(df,t,interval)
    # TM.to_csv(r'C:\Uganda\TM_t'+str(t)+'_dt'+str(interval)+'.csv')