import pathlib
import pandas as pd
import re
import itertools
import plotly.express as px
import plotly.graph_objects as go


def main():
    df_it = pd.read_csv(pathlib.Path('data/IT_series/IT_iChiroPro.csv'))
    df_it["ITSensor"] = None

    for file in itertools.chain(
            pathlib.Path('data/IT_series/element').iterdir(),
            pathlib.Path('data/IT_series/nevo').iterdir()
    ):
        if len(re.findall("T[0-9]+_[0-9]{2}\.csv", file.name)) == 0:
            continue

        df_sample = pd.read_csv(file)
        df_it.loc[df_it["SampleNo"] == file.stem, "ITSensor"] = df_sample["Channel 1"].max() * 10

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=df_it["SampleNo"], y=df_it["ITSensor"], name="ITSensor"))
    fig.add_trace(go.Scatter(x=df_it["SampleNo"], y=df_it["ITChiro"], name="ITChiro"))
    fig.show()
    print()


if __name__ == '__main__':
    main()
