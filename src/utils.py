import pandas as pd

from sklearn.preprocessing import StandardScaler

def shuffle_df(df : pd.DataFrame) -> pd.DataFrame:
    return df.sample(frac=1).reset_index(drop=True)

def create_shuffled_test_df(only_attack_df : pd.DataFrame, only_norm_df : pd.DataFrame) -> pd.DataFrame:
    only_attack_df['label'] = 'attack'
    only_norm_df['label']   = 'normal'
    test_ds = only_attack_df.append(only_norm_df)
    return shuffle_df(test_ds)

def normalize_df(df, columns, scalar = StandardScaler()):
    return pd.DataFrame(scalar.fit_transform(df), columns=columns)