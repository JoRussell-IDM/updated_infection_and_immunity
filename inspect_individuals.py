import os
import logging

import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style('white')
sns.set_context('paper')

log = logging.getLogger(__name__)

# data_folder = os.path.join('~', 'Dropbox', 'uganda_prism')
data_folder = "Q:/Malaria/Uganda/Final PRISM cohort study databases through June 2016"


def read_data():
    """
    Read complete data file and transform a few string channels to booleans.
    Other interesting channels:
    - itnlastnight
    - temperature, hemoglobin, sevmalcat
    - fever, fduration (fatigue, abdominalpain, anorexia, vomiting, diarrhea, cough, headache, jointpains, muscleaches, seizure, jaundice)
    - dx1code, med1code
    - careoutside
      - wcaregiven1, cdate1, diagnosis1, bsdone1
      - osantimalarial
        - antimalarial1, dose1, lgdate1
    - traveloutsc
      - datefirstnight1, datelastnight1, district1, parish1, sleptunderitn1, reasontravel1
    """

    df = pd.read_csv(os.path.join(data_folder, 'all_data_w_lat_long.csv'),
                     index_col=0, low_memory=False, parse_dates=['date'])
    df['fever'] = df.febrile.apply(lambda x: True if x == 'Yes' else False if x == 'No' else x)  # febrile = reported fever or measured high temperature?
    df['LAMP'] = df.LAMP.apply(lambda x: True if x == 'Positive' else False if x == 'Negative' else x)
    df['id'] = df.id.astype(int)


    sero = pd.read_csv(os.path.join(data_folder, 'serology_data.csv'),
                       usecols=['id', 'unique_id', 'Conc2', 'antigen']
                       )
    sero['date'] = pd.to_datetime(sero.unique_id.map(lambda x: '-'.join(x.split('-')[1:])))
    concentrations = sero.groupby(['id', 'date', 'antigen']).Conc2.mean().unstack('antigen')
    # print(concentrations.head())

    df = df.set_index(['id', 'date']).join(concentrations, how='left').reset_index()

    return df

def extract_events(df):
    """
    Extract relevant info to categorize patient state
    Routine visits (Enrollment, Routine):
    - regardless of fever --> microscopy --> if micro neg: LAMP

    Non-routine visits (T-cell count, non-Routine):
    - fever --> if fever pos: microscopy; else: usually no test (occasionally microscopy)
    """

    df = df[['date','id','siteid','fever','agecat','parasitedensity', 'gametocytes', 'LAMP', 'visittype', 'age', 'gender','hhid','Longitude','Latitude','hh_biting_propensity_wt','anymalaria']].copy()

    df_events = df ## Keeping the nans here so as to not improperly call them as true.
    df_events.loc[pd.notnull(df_events.parasitedensity), 'micro'] = (df_events.parasitedensity>0)
    df_events['routine'] = df_events.visittype.isin(['Routine', 'Enrollment'])
    events = df_events.groupby(['routine', 'fever', 'micro', 'LAMP']).id.count().rename('count')
    log.debug(events)

    return df_events

def extract_events_malariatherapy(df):
    """
    Extract relevant info to categorize patient state
    Routine visits (Enrollment, Routine):
    - regardless of fever --> microscopy --> if micro neg: LAMP

    Non-routine visits (T-cell count, non-Routine):
    - fever --> if fever pos: microscopy; else: usually no test (occasionally microscopy)
    """

    df_ = df[['Day','patient','Inst','Temp','Asexual', 'Gametocytes', 'gender','route']].copy()

    df_events = df_ ## Keeping the nans here so as to not improperly call them as true.
    df_events.loc[pd.notnull(df_events.Asexual), 'micro'] = (df_events.Asexual>0) or (df_events.Gametocytes>0)
    return df_events


def extract_antibodies(df):
    return df.loc[:, 'AMA1':].copy()


def human_format(num):
    magnitude = 0
    while abs(num) >= 1000:
        magnitude += 1
        num /= 1000.0
    # add more suffixes if you need them
    return '%.3g%s' % (num, ['', 'k', 'M', 'G'][magnitude])


def plot_events(df, suptitle='', fig=None, ax=None, detailed=False):

    if not fig:
        fig, ax = plt.subplots(1, 1, num='events', figsize=(20, 10))

    cmap = plt.cm.YlOrRd
    lamp_color = cmap(0.1)

    # Indexing of unique individuals sorted by age (for position in scatter)
    df.sort_values(by='age', inplace=True)
    unique_ids = df.drop_duplicates(subset=['id'], keep='first').copy()
    unique_ids['idx'] = range(len(unique_ids))
    df = df.merge(unique_ids[['id', 'idx']], on='id', how='outer')

    # Marker for visit (red if febrile)
    ax.scatter(df.date.values, df.idx.values, c='darkgray', s=5, alpha=1, lw=0)
    fever = df[df.fever]
    ax.scatter(fever.date.values, fever.idx.values, c='firebrick', s=5, alpha=1, lw=0)

    sero = df[pd.notnull(df.AMA1)]
    if detailed:
        ax.scatter(sero.date.values, sero.idx.values + 0.75, c='k', s=100, zorder=-10, marker=u'$\u2193$', alpha=0.5)

    # LAMP done only for microscopy negatives
    # Fill with low-density coloring and small marker if positive
    # If negative, show open gray marker to represent upper bound on parasitemia that is excluded by negative LAMP
    lamp = df[pd.notnull(df.LAMP)].copy()
    lamp['LAMP'] = lamp.LAMP.astype(bool)
    lamp_neg = lamp[~lamp.LAMP]
    ax.scatter(lamp_neg.date.values, lamp_neg.idx.values, s=30, alpha=0.5, c='None', edgecolors='darkgray', lw=0.5)
    if detailed:
        for i in range(len(lamp_neg)):
            ax.annotate('l', (lamp_neg.date.values[i],  lamp_neg.idx.values[i] - 0.5),
                        fontsize=8, alpha=0.5, ha='center')
    lamp_positive = lamp[lamp.LAMP]
    ax.scatter(lamp_positive.date.values, lamp_positive.idx.values, s=60, alpha=0.25, c=lamp_color, edgecolors='darkgray', lw=0.5)
    if detailed:
        for i in range(len(lamp_positive)):
            ax.annotate('L', (lamp_positive.date.values[i],  lamp_positive.idx.values[i] - 0.5),
                        fontsize=9, fontweight='bold', alpha=0.8, ha='center')

    # Microscopy done for all febrile and routine visits
    # For individuals with no downstream LAMP test,
    #   show open gray marker to represent upper bound on parasitemia that is excluded by negative smear
    micro_only = df[(pd.isnull(df.LAMP)) & (pd.notnull(df.parasitedensity))].copy()
    micro_neg = micro_only[micro_only.parasitedensity == 0]
    ax.scatter(micro_neg.date.values, micro_neg.idx.values, s=60, alpha=0.5, c='None', edgecolors='darkgray', lw=0.5)
    if detailed:
        for i in range(len(micro_neg)):
            ax.annotate('m', (micro_neg.date.values[i],  micro_neg.idx.values[i] - 0.5),
                        fontsize=8, alpha=0.5, ha='center')

    # For individuals with positive smear,
    #   size and color of markers represent measured parasitemia
    parasite_positive = micro_only[micro_only.parasitedensity > 0]
    marker_size = 50 * np.log10(parasite_positive.parasitedensity)
    marker_size[marker_size < 60] = 60
    marker_color = parasite_positive.parasitedensity.copy()
    # marker_color[parasite_positive.fever] = 1e5
    ax.scatter(parasite_positive.date.values, parasite_positive.idx.values,
               c=marker_color, s=marker_size, cmap=cmap, norm=LogNorm(10, 3e5),
               alpha=0.25, lw=0.5, edgecolors='darkgray')

    parasite_positive_with_gam = parasite_positive[parasite_positive.gametocytes > 0]
    marker_size = 50 * np.log10(parasite_positive_with_gam.parasitedensity)
    marker_size[marker_size < 60] = 60
    ax.scatter(parasite_positive_with_gam.date.values, parasite_positive_with_gam.idx.values,
               s=marker_size, alpha=0.5, c='None', edgecolors='olive', lw=1.5)

    if detailed:
        for i, v in enumerate(parasite_positive.parasitedensity.values):
            txt = human_format(v)
            ax.annotate(txt, (parasite_positive.date.values[i], parasite_positive.idx.values[i] + 0.25),
                        fontsize=9, fontweight='bold', alpha=0.8, ha='center', color='firebrick' if parasite_positive.fever.values[i] else 'k')

    # Label y-axis with age and gender (subsampled to fit)
    ax.get_yaxis().get_major_formatter().set_useOffset(False)
    subsample = 1 + len(unique_ids) / 100
    ax.get_yaxis().set_ticks(unique_ids.idx[::subsample])
    ax.get_yaxis().set_ticklabels(unique_ids.apply(lambda x: '%d %s' % (x.age, x.gender[0].upper()), axis=1)[::subsample])
    ax.set(ylim=(df.idx.min() - 1, df.idx.max() + 1))
    fig.set_tight_layout(True)
    if suptitle:
        fig.suptitle(suptitle, y=1)

    import datetime
    for i in range(len(unique_ids)):
        ax.plot([datetime.date(2011, 8, 20), datetime.date(2016, 7, 5)], [i, i], color='gray', alpha=0.1, lw=0.25)


def plot_events_malariatherapy(df, suptitle='', subtitle = '', fig=None, ax=None, detailed=False):

    if not fig:
        fig, ax = plt.subplots(1, 1, num='events', figsize=(20, 10))

    cmap = plt.cm.YlOrRd
    cmap_gam = plt.cm.GnBu
    lamp_color = cmap(0.1)

    # Indexing of unique individuals sorted by age (for position in scatter)

    unique_ids = df.drop_duplicates(subset=['id'], keep='first').copy()
    unique_ids['idx'] = range(len(unique_ids))
    df = df.merge(unique_ids[['id', 'idx']], on='id', how='outer')

    # Marker for visit (red if febrile)
    ax.scatter(df.Day.values, df.idx.values, c='darkgray', s=5, alpha=1, lw=0)

    Treatment = df[(df['Tretmentfree?'] == '0') | (df['Tretmentfree?'] == False)]
    #ax.scatter(fever.Day.values, fever.idx.values, c='firebrick', s=8, alpha=1, lw=0)


    # Microscopy done for all febrile and routine visits
    # For individuals with no downstream LAMP test,
    #   show open gray marker to represent upper bound on parasitemia that is excluded by negative smear
    micro_only = df[pd.notnull(df.Asexual)].copy()
    micro_neg = micro_only[micro_only.Asexual == 0]
    ax.scatter(micro_neg.Day.values, micro_neg.idx.values, s=60, alpha=0.5, c='None', edgecolors='darkgray', lw=0.5)

    # For individuals with positive smear,
    #   size and color of markers represent measured parasitemia
    #treated = micro_only.loc[micro_only['Tretmentfree?'].isin(['0', False])]
    #untreated = micro_only.loc[micro_only['Tretmentfree?'].isin(['1', True])]
    fever = micro_only[micro_only.fever]
    parasite_positive = micro_only[(pd.to_numeric(micro_only.Asexual, errors = 'coerce') > 0) | (pd.to_numeric(micro_only.Gametocytes, errors = 'coerce') >0)]


    Asexual_marker_size = 50*np.log10(pd.to_numeric(parasite_positive.Gametocytes,errors='coerce')) #change for Asexual!!!!!!


    Gam_marker_size = 50 * np.log10(pd.to_numeric(parasite_positive.Gametocytes, errors='coerce'))
    Asexual_marker_size[Asexual_marker_size < 60] = 60
    Gam_marker_size[Gam_marker_size < 60] = 60
    Asexual_marker_color = [pd.to_numeric(parasite_positive.Asexual,errors='coerce')]
    Gam_marker_color = [pd.to_numeric(parasite_positive.Gametocytes,errors='coerce')]
    # marker_color[parasite_positive.fever] = 1e5
    ax.scatter(parasite_positive.Day.values, parasite_positive.idx.values,
               c=Asexual_marker_color, s=Asexual_marker_size, cmap=cmap_gam, norm=LogNorm(10, 3e5),
               alpha=0.35, lw=0.5, edgecolors='darkgray')
    #ax.scatter(parasite_positive.Day.values, parasite_positive.idx.values,
    #           c=Gam_marker_color, s=Gam_marker_size, cmap=cmap_gam, norm=LogNorm(10, 3e5),
    #           alpha=0.25, lw=0.5, edgecolors='darkgray')

    ax.scatter(fever.Day.values, fever.idx.values, c='firebrick', s=6, alpha=1, lw=0)
   # ax.scatter(Treatment.Day.values, Treatment.idx.values, c='forestgreen', s=6, alpha=1, lw=0)

    # parasite_positive_with_gam = parasite_positive[parasite_positive.Gametocytes > 0]
    # marker_size = 50 * np.log10(parasite_positive_with_gam.Gametocytes)
    # marker_size[marker_size < 60] = 60
    # ax.scatter(parasite_positive_with_gam.Day.values, parasite_positive_with_gam.idx.values,
    #            s=marker_size, alpha=0.5, c='None', edgecolors='olive', lw=1.5)


    # Label y-axis with age and gender (subsampled to fit)
    ax.get_yaxis().get_major_formatter().set_useOffset(False)
    subsample = 1 + len(unique_ids) / 100
    ax.get_yaxis().set_ticks(unique_ids.idx[::subsample])
    ax.get_yaxis().set_ticklabels(unique_ids.apply(lambda x: '%s %s' % (x.id, str(x.Sex).upper()), axis=1)[::subsample])
    ax.set(ylim=(df.idx.min() - 1, df.idx.max() + 1))
    ax.set(xlim=(0,300))
    ax.set_xlabel('Days (since infection)')
    ax.set_ylabel('Individual (id, gender)')
    fig.set_tight_layout(True)
    if suptitle:
        fig.suptitle(suptitle, y=1)
    if subtitle:
        ax.set_title(subtitle)
    import datetime
    for i in range(len(unique_ids)):
         ax.plot([1, max(df.Day)], [i, i], color='gray', alpha=0.1, lw=0.25)


def plot_events_Garki(df, suptitle='', subtitle = '', fig=None, ax=None, detailed=False):

    if not fig:
        fig, ax = plt.subplots(1, 1, num='events', figsize=(20, 10))

    cmap = plt.cm.YlOrRd
    cmap_gam = plt.cm.GnBu
    lamp_color = cmap(0.1)

    # Indexing of unique individuals sorted by age (for position in scatter)

    unique_ids = df.drop_duplicates(subset=['id'], keep='first').copy()
    unique_ids['idx'] = range(len(unique_ids))
    df = df.merge(unique_ids[['id', 'idx']], on='id', how='outer')
    df.datecoll = df.datecoll.apply(lambda x: np.datetime64(x))
    df.asexual_density[df.asexual_density == 'maxed'] = float(2205)
    # Marker for visit (red if febrile)
    ax.scatter(df.datecoll.values, df.idx.values, c='darkgray', s=5, alpha=1, lw=0)

    # Treatment = df[(df['Tretmentfree?'] == '0') | (df['Tretmentfree?'] == False)]
    #ax.scatter(fever.Day.values, fever.idx.values, c='firebrick', s=8, alpha=1, lw=0)


    # Microscopy done for all febrile and routine visits
    # For individuals with no downstream LAMP test,
    #   show open gray marker to represent upper bound on parasitemia that is excluded by negative smear
    micro_only = df[pd.notnull(df.asexual_density)].copy()
    micro_neg = micro_only[micro_only.asexual_density == 0]
    ax.scatter(micro_neg.datecoll.values, micro_neg.idx.values, s=60, alpha=0.5, c='None', edgecolors='darkgray', lw=0.5)

    # For individuals with positive smear,
    #   size and color of markers represent measured parasitemia
    #treated = micro_only.loc[micro_only['Tretmentfree?'].isin(['0', False])]
    #untreated = micro_only.loc[micro_only['Tretmentfree?'].isin(['1', True])]
    fever = micro_only[micro_only.fever == 1]
    parasite_positive = micro_only[(pd.to_numeric(micro_only.asexual_density, errors = 'coerce') > 0)]# | (pd.to_numeric(micro_only.Gametocytes, errors = 'coerce') >0)]


    Asexual_marker_size = 50*np.log10(pd.to_numeric(parasite_positive.asexual_density,errors='coerce')) #change for Asexual!!!!!!


   # Gam_marker_size = 50 * np.log10(pd.to_numeric(parasite_positive.Gametocytes, errors='coerce'))
    Asexual_marker_size[Asexual_marker_size < 60] = 60
    #Gam_marker_size[Gam_marker_size < 60] = 60
    Asexual_marker_color = [pd.to_numeric(parasite_positive.asexual_density,errors='coerce')]
    #Gam_marker_color = [pd.to_numeric(parasite_positive.Gametocytes,errors='coerce')]
    # marker_color[parasite_positive.fever] = 1e5
    ax.scatter(parasite_positive.datecoll.values, parasite_positive.idx.values,
               c=Asexual_marker_color, s=Asexual_marker_size, cmap=cmap, norm=LogNorm(10, 3e5),
               alpha=0.35, lw=0.5, edgecolors='darkgray')
    #ax.scatter(parasite_positive.Day.values, parasite_positive.idx.values,
    #           c=Gam_marker_color, s=Gam_marker_size, cmap=cmap_gam, norm=LogNorm(10, 3e5),
    #           alpha=0.25, lw=0.5, edgecolors='darkgray')

    ax.scatter(fever.datecoll.values, fever.idx.values, c='firebrick', s=6, alpha=1, lw=0)
   # ax.scatter(Treatment.Day.values, Treatment.idx.values, c='forestgreen', s=6, alpha=1, lw=0)

    # parasite_positive_with_gam = parasite_positive[parasite_positive.Gametocytes > 0]
    # marker_size = 50 * np.log10(parasite_positive_with_gam.Gametocytes)
    # marker_size[marker_size < 60] = 60
    # ax.scatter(parasite_positive_with_gam.Day.values, parasite_positive_with_gam.idx.values,
    #            s=marker_size, alpha=0.5, c='None', edgecolors='olive', lw=1.5)


    # Label y-axis with age and gender (subsampled to fit)
    ax.get_yaxis().get_major_formatter().set_useOffset(False)
    subsample = 1 + len(unique_ids) / 100
    ax.get_yaxis().set_ticks(unique_ids.idx[::subsample])
    ax.get_yaxis().set_ticklabels(unique_ids.apply(lambda x: '%s %s' % (int(x.age_at_enroll), str(x.sex).upper()), axis=1)[::subsample])
    ax.set(ylim=(df.idx.min() - 1, df.idx.max() + 1))
    #ax.set(xlim=(0,24))
    ax.set_xlabel('Date')
    ax.set_ylabel('Individual (age, gender)')
    fig.set_tight_layout(True)
    if suptitle:
        fig.suptitle(suptitle, y=1)
    if subtitle:
        ax.set_title(subtitle)
    import datetime
  #  for i in range(len(unique_ids)):
   #     ax.plot([1, max(df.datecoll)], [i, i], color='gray', alpha=0.1, lw=0.25)



def plot_individual(ind_df, site='', fig=None, axs=None):

    if not fig:
        fig, axs = plt.subplots(7, 1, num='individual', figsize=(16, 8), sharex=True)

    plot_events(ind_df, ': '.join([site, str(ind_df.id.iloc[0])]), fig, axs[0], detailed=True)

    sero_df = ind_df.dropna(subset=['AMA1']).set_index('date').loc[:, 'AMA1':]

    for i, ab in enumerate(['GEXP1', 'Etramp5.Ag1', 'MSP1.19', 'GLURP', 'MSP2.CH150.9', 'AMA1']):
        ax = axs[i + 1]
        ax.plot(sero_df.index.values, sero_df[ab], '-o')
        ax.set(yscale='log')
        ax.set_title(ab, y=0.7)

    fig.set_tight_layout(True)


if __name__ == '__main__':

    logging.basicConfig(format='%(message)s')
    log.setLevel(logging.INFO)

    ind_df = read_data()

    # id_range = (3001, 3020)
    # ind = ind_df[ind_df.id.isin(range(*id_range))]

    log.debug(ind_df.siteid.value_counts())

    site = 'Tororo'  # Tororo, Kanungu, Jinja

    ind = ind_df[ind_df.siteid == site]
    events_df = extract_events(ind)

    antibodies_df = extract_antibodies(ind)
    events_df = events_df.join(antibodies_df)
    sampled_ids = events_df.dropna(subset=['AMA1']).id.unique()

    # print(events_df.visittype.value_counts())
    # events_df = events_df[events_df.visittype == 'Not routine']
    # events_df = events_df[events_df.visittype != 'Not routine']  # enrollment + routine

    events_df = events_df[events_df.id.isin(sampled_ids)]
    # for i, (ind, ind_df) in enumerate(events_df.groupby('id')):
    #     plot_individual(ind_df, site)
    #     if i > 10:
    #         break
    #     plt.show()

