from pathlib import Path
from matplotlib import pyplot as plt

from helper import get_mts_df

if __name__ == "__main__":
    mts_df = get_mts_df(Path("data/mts"))
    gb = mts_df.groupby('ID')
    mts_dfs = [gb.get_group(x) for x in gb.groups]

    for df in mts_dfs:
        title = df['ID'][0]
        df.plot(
            x='Axial Displacement',
            y='Axial Force',
            xlabel='Axial Displacement [mm]',
            ylabel='Axial Force [N]',
            legend=False,
            kind='line',
            title=title
        )
        plt.show()
