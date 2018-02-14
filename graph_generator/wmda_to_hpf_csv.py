import pandas as pd
wmda_orig = pd.read_csv('data/wmda/freqs.txt', sep=';', names = ['hap', 'freq'])
wmda_orig['pop'] = 'CAU'
wmda = wmda_orig[['hap', 'pop', 'freq']]
wmda.to_csv('data/wmda/hpf.csv',index=False)
