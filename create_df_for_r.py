from __future__ import annotations

from pathlib import Path

import pandas as pd

from helper import (
    get_bone_morpho_df,
    get_disp_at_max_df,
    get_it_df,
    get_mts_df,
    get_sample_df,
    get_stiffness_df,
    get_uf_df,
    strip_versioning,
)

if __name__ == "__main__":
    strip_versioning(Path(__file__).parent / "data")

    df_morpho = get_bone_morpho_df(Path("data/bone_morpho"))
    df_samples = get_sample_df(Path("data"))
    df_it = get_it_df(Path("data/it_ichiro"))
    df_mst = get_mts_df(Path("data/mts"))
    df_uf = get_uf_df(Path("data/mts"))

    df_st = get_stiffness_df(Path("data/mts"))
    df_dm = get_disp_at_max_df(Path("data/mts"))

    df = pd.merge(df_samples, df_morpho, on="MeasNoHighRes", how="outer")
    df = pd.merge(df, df_it, on="ID", how="outer")
    df = pd.merge(df, df_uf, on="ID", how="outer")
    df = pd.merge(df, df_st, on="ID", how="outer")
    df = pd.merge(df, df_dm, on="ID", how="outer")

    # remove unused samples
    df = df[df["selected"] == "yes"]

    df = df.sort_values(by=["group", "BVTV"])

    df_dm = df.copy(deep=True)
    df_dm = df_dm[["BVTV", "ID", "group", "DM"]]
    df_dm.to_csv(Path("data_generated/df_dm.csv"), index=False)

    df_st = df.copy(deep=True)
    df_st = df_st[["BVTV", "ID", "group", "ST"]]
    df_st.to_csv(Path("data_generated/df_st.csv"), index=False)

    df_uf = df.copy(deep=True)
    df_uf = df_uf[["BVTV", "ID", "group", "UF"]]
    df_uf.to_csv(Path("data_generated/df_uf.csv"), index=False)

    df_it = df.copy(deep=True)
    df_it = df_it[["BVTV", "ID", "group", "IT"]]
    # remove all the measurements where we have hit the maximal insertion torque of ichiro.
    df_it = df_it[df_it["IT"] < 70.5]
    df_it.to_csv(Path("data_generated/df_it.csv"), index=False)

    df.to_csv(Path("data_generated/df_all.csv"), index=False)
