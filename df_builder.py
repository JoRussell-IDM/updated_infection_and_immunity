import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import psycopg2
import psycopg2.extensions
import numpy as np
from datetime import datetime
from inspect_individuals import plot_events_malariatherapy
from inspect_individuals import plot_events_Garki

data_folder = 'C:\Uganda'

# load the DF for Uganda

# df_PRISM = pd.read_csv(os.path.join(data_folder, 'Tororo_analyzed_2017_05_08.csv'), index_col=0)

# load the DF for Malariatherapy

df_therapy = pd.read_csv('Malariatherapy.txt',sep='\t')

# classify the patient histories by whether they received treatment at some point (Treated = 1)

for id,measurements in df_therapy.groupby('patient'):
    durations = pd.to_numeric(measurements.loc[:, 'Asexual'], errors='coerce')
    Temps = pd.to_numeric(measurements.loc[:,'Temp'], errors='coerce')
    Treats = pd.to_numeric(measurements.loc[:, 'Tretmentfree?'], errors='coerce')

    measurements.Asexual = durations
    asexual_infection_duration_index = max(i for i in measurements.Asexual.index if measurements.Asexual[i] > 0)
    gam_inf_duration_index = max(i for i in measurements.Gametocytes.index if measurements.Gametocytes[i] > 0)
    asexual_infection_duration = measurements.Day[asexual_infection_duration_index]
    gametocyte_infection_duration = measurements.Day[gam_inf_duration_index]
    if np.sum(pd.to_numeric(measurements.Gametocytes, errors='coerce'))>0:
        max_gam_index = measurements[pd.to_numeric(measurements.Gametocytes, errors='coerce') == max(pd.to_numeric(measurements.Gametocytes, errors='coerce'))].index
        max_gam_day = measurements[pd.to_numeric(measurements.Gametocytes, errors='coerce') == max(pd.to_numeric(measurements.Gametocytes, errors='coerce'))].Day
        df_therapy.loc[max_gam_index, 'gam_max'] = max_gam_day

    df_therapy.loc[df_therapy.patient == id,'asexual_infection_duration'] = asexual_infection_duration
    df_therapy.loc[df_therapy.patient == id, 'gametocyte_infection_duration'] = gametocyte_infection_duration
    df_therapy.loc[df_therapy.patient == id,'fever'] = Temps>37
    df_therapy.loc[df_therapy.patient == id,'Treats'] = Treats == 0

    if sum(df_therapy.loc[df_therapy.patient == id,'Treats']) >0:
        df_therapy.loc[df_therapy.patient == id, 'Treated'] = 1
    elif sum(df_therapy.loc[df_therapy.patient == id,'Treats']) ==0:
        df_therapy.loc[df_therapy.patient == id, 'Treated'] = 0
    else:
        df_therapy.loc[df_therapy.patient == id, 'Treated'] = 2

#sort infections by duration of positive parasite count
df_therapy.sort_values(by=['asexual_infection_duration','patient'], inplace=True)
new_columns = df_therapy.columns.values
new_columns[0] = 'id'
df_therapy.columns = new_columns

df_therapy= df_therapy[df_therapy.Treated == 0]
path = r'C:\Uganda'
#
# save malariatherapy df for future use
df_therapy.to_csv(os.path.join(path, r'Malariatherapy_untreated_2017_07_10.csv'))

#Parse Garki dataset from server as a DF
# try:
#     con = psycopg2.connect(host='ivlabsdssql01.na.corp.intven.com', port=5432, dbname="idm_db")
# except psycopg2.Error:
#     raise Exception(
#         "Failed connection to %s.  Is it possible that other VPN connections (e.g. Sabalcore) are interfering?" % 'ivlabsdssql01.na.corp.intven.com')
#
# query = "SELECT t.id,p.sex,t.datecoll, t.survey, t.village, v.vname, v.area, t.compound, p.dateb, t.fever, t.pfa, t.pfg, t.exam\
#                FROM malaria.parademo t \
#                INNER JOIN malaria.people p ON p.id = t.id \
#                INNER JOIN malaria.village v ON v.village = t.village \
#                WHERE p.dateb <> '0000-00-00' AND t.datecoll <> '0000-00-00' AND t.exam > 0 \
#                ORDER BY t.id, t.survey"
# cursor = con.cursor()
# cursor.execute(query)
# df_Garki = pd.DataFrame(cursor.fetchall())
# df_Garki.columns = [desc[0] for desc in cursor.description]
# cursor.close()
# vnames = df_Garki[df_Garki.survey>16].vname.unique()
# df_Garki = df_Garki[df_Garki.vname.isin(vnames)]
#
# df_Garki.datecoll = pd.to_datetime(df_Garki.datecoll, format= '%Y-%m-%d')
#
# #manually edit an errant datetime in the dataset
# df_Garki.dateb.replace(to_replace = '1973-09-31', value = '1973-09-01', inplace = True)
#
# #set new fields for DF (densities, age)
#
# df_Garki.dateb = pd.to_datetime(df_Garki.dateb, format= '%Y-%m-%d')
# df_Garki['age']= (df_Garki.datecoll - df_Garki.dateb )/ np.timedelta64(1, 'Y')
# df_Garki['IRS_status'] = ['pre_IRS' if date < datetime(1973,04,01) else 'during_IRS' if date < datetime(1974,10,31) else 'post_IRS' for date in df_Garki['datecoll']]
# df_Garki['asexual_density'] = [-416*np.log(1-float(df_Garki.loc[row, 'pfa'])/(df_Garki.loc[row,'exam'])) if df_Garki.loc[row,'pfa'] != df_Garki.loc[row,'exam'] else 'maxed' for row in df_Garki.index.values] # 416 times the value per frame to get per uL, 200 frames == 0.48uL
# df_Garki['gametocyte_density'] = [-416*np.log(1-float(df_Garki.loc[row, 'pfg'])/(df_Garki.loc[row,'exam'])) if df_Garki.loc[row,'pfg'] != df_Garki.loc[row,'exam'] else 'maxed' for row in df_Garki.index.values] # 416 times the value per frame to get per uL, 200 frames == 0.48uL
#
# for id, measurements in df_Garki.groupby('id'):
#      age_enroll = min(measurements.age)
#      df_Garki.loc[df_Garki.id == id,'age_at_enroll'] = age_enroll
#
# _ = df_Garki.sort_values(by=['age_at_enroll','id'], inplace=True)
#
# df_Garki.to_csv(os.path.join(path, r'DF_Garki.csv'))
# df_Garki = pd.read_csv('C:\Uganda\DF_Garki.csv', index_col=0,parse_dates=['datecoll'])
#
# fig, (ax1,ax2,ax3) = plt.subplots(nrows = 3, sharex = True)
#
# plot_events_Garki(df_Garki[(df_Garki.vname.isin(['Ajura','Nabanawa'])) & (df_Garki.sex ==1)], subtitle = 'Control Sites',fig =fig, ax =ax1)
# plot_events_Garki(df_Garki[(df_Garki.vname.isin(['Kukar Maikiya','Rafin Marke','Kawari'])) & (df_Garki.sex ==1)], subtitle ='West sites',fig =fig, ax =ax2)
# plot_events_Garki(df_Garki[(df_Garki.vname.isin(['Kargo Kudu','Nasakar','Bakan Sabara'])) & (df_Garki.sex ==1)], subtitle ='East Sites',fig =fig, ax = ax3)
#
# plt.show()
fig, (ax1,ax2) = plt.subplots(2, 1, sharex=True, num='events', figsize=(20, 10))
plot_events_malariatherapy(df_therapy[df_therapy.Treated == 1], 'Malariatherapy', 'Malariatherapy treated infections', fig, ax1)
plot_events_malariatherapy(df_therapy[df_therapy.Treated == 0], 'Malariatherapy', 'Malariatherapy untreated infections', fig, ax2)
plt.show()

print('oh wow')

### Notes
# x = age, y = incidence, hue = pre intervention (high EIR), post intervention (low EIR)
# x =  time until first infection, y = time until subsequent infection, hue = symptomatic
# x = age corrected number of infections, y = fraction symptomatic
# x = average frequency of lifetime infections, y = fraction symptomatic
# x = time since last infection, y = fraction symptomatic, hue = symptomatic

# how to efficiently calculate incidence by age bracket?
# cumulatively calculate new cases of disease and divide by the cumulative time in person years of observation
