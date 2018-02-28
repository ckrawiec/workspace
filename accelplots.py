import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import datetime as dt

def extvsintwind():
    env1 = pd.read_csv('/Users/Christina/DES/accel/envdata_windspeed_042114_051514.csv')
    env2 = pd.read_csv('/Users/Christina/DES/accel/envdata_windspeed_091714_102114.csv')

    ane1 = pd.read_csv('/Users/Christina/DES/accel/anemom042114_051514.csv')
    ane2 = pd.read_csv('/Users/Christina/DES/accel/anemom_09172014_10212014.csv')

    df_env = pd.concat([env1, env2])
    df_ane = pd.concat([ane1, ane2])

    x_env, x_ane = [], []
    #some data doesn't have decimal seconds
    for d in df_env['time_recorded']:
        try:
            x_env.append(dt.datetime.strptime(d,'%Y-%m-%d %H:%M:%S.%f+00:00').date())
        except ValueError:
            x_env.append(dt.datetime.strptime(d,'%Y-%m-%d %H:%M:%S+00:00').date())

    ane_mask = [np.isnan(d) for d in df_ane['time_recorded']]
    for d in df_ane['time_recorded'][~ane_mask]:
        print d, d[1]
        try:
            x_ane.append(dt.datetime.strptime(d[1],'%Y-%m-%d %H:%M:%S.%f+00:00').date())
        except ValueError:
            x_ane.append(dt.datetime.strptime(d[1],'%Y-%m-%d %H:%M:%S+00:00').date())

    plt.scatter(x_env, df_env['wind_speed'])
    plt.scatter(x_ane, df_ane['magnitude'][ane_mask])

    plt.show()


def main():
    extvsintwind()

if __name__=="__main__":
    main()
