from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


def strip_versioning(file: Path) -> Path | None:
    if not file.exists():
        return None

    if file.is_dir():
        for inner_file in file.iterdir():
            strip_versioning(inner_file)

    suffix_split = file.suffix.split(";")
    if len(suffix_split) == 1:
        return file

    name = file.stem + "_" + suffix_split[1]
    if not file.is_dir():
        name += suffix_split[0]

    file = file.rename(file.parent / name)

    return file


def read_results_bone_morpho(file_name: Path, version: int) -> dict:
    with file_name.open("r") as f:
        names = f.readline()
        for _ in range(version):
            values = f.readline()

    names = names.split("\t")
    values = values.split("\t")

    val = {a.strip(): b.strip() for (a, b) in zip(names, values)}

    return val


def results_bone_morpho_in_dir(data_dir: Path, version: int = 1) -> dict | None:
    if not data_dir.exists():
        return None
    results_bone_morpho = list(data_dir.glob("*3DRESULTS_BONE_MORPHO_1.TXT"))
    if len(results_bone_morpho) != 1:
        return None

    return read_results_bone_morpho(results_bone_morpho[0], version)


def get_bone_morpho_df(dir_name: Path, version: int = 2) -> pd.DataFrame:
    vals = list()
    for file_name in dir_name.iterdir():
        if "3DRESULTS_BONE_MORPHO_1.TXT" not in file_name.name:
            continue
        val = read_results_bone_morpho(file_name, version)
        vals.append(val)

    df_morpho = pd.DataFrame(vals, columns=["MeasNo", "VOX-BV/TV"])
    df_morpho.rename(
        columns={"VOX-BV/TV": "BVTV", "MeasNo": "MeasNoHighRes"}, inplace=True
    )
    df_morpho["BVTV"] = df_morpho["BVTV"].astype(float)
    df_morpho["MeasNoHighRes"] = df_morpho["MeasNoHighRes"].astype(int)
    df_morpho = df_morpho.sort_values(by=["MeasNoHighRes"])

    return df_morpho


def get_sample_df(dir_name: Path) -> pd.DataFrame:
    df_samples = pd.read_excel(dir_name / "samples.xlsx")
    df_samples = df_samples[["MeasNoHighRes", "ID", "group", "selected"]]

    return df_samples


def get_it_df(dir_name: Path) -> pd.DataFrame:
    dfs_it = list()
    for file_name in dir_name.glob("*.csv"):
        dfs_it.append(pd.read_csv(file_name, delimiter=";"))

    df_it = pd.concat(dfs_it)
    df_it.rename(
        columns={"Referenz": "ID", "Max. erreichtes Drehmoment": "IT"}, inplace=True
    )
    df_it = df_it.groupby("ID", as_index=False).first()

    df_it = df_it[["IT", "ID"]]

    return df_it


def get_mts_df(dir_name: Path) -> pd.DataFrame:
    dfs_mts = list()
    for file_name in dir_name.glob("*.csv"):
        df_mts = pd.read_csv(file_name, header=1)
        df_mts = df_mts[1:]  # Remove the first row to reindex the DataFrame
        df_mts.reset_index(drop=True, inplace=True)  # Reset index if desired
        df_mts["ID"] = file_name.stem
        dfs_mts.append(df_mts)

    df_mts = pd.concat(dfs_mts)
    df_mts.columns = [f.strip() for f in df_mts.columns]

    df_mts["Time"] = df_mts["Time"].astype(float)
    df_mts["Axial Count"] = df_mts["Axial Count"].astype(float)
    df_mts["Axial Force"] = -df_mts["Axial Force"].astype(float)
    df_mts["Axial Displacement"] = -df_mts["Axial Displacement"].astype(float)

    return df_mts


def get_uf_df(dir_name: Path) -> pd.DataFrame:
    df_mts = get_mts_df(dir_name)
    df_uf = df_mts.groupby("ID")["Axial Force"].max()
    df_uf.rename("UF", inplace=True)

    return df_uf


def get_disp_at_max_df(dir_name: Path) -> pd.DataFrame:
    df_mts = get_mts_df(dir_name).reset_index(drop=True)
    idx = df_mts.groupby("ID")["Axial Force"].idxmax()
    df_dm = df_mts.iloc[idx]
    df_dm = df_dm[["ID", "Axial Displacement"]].set_index("ID").squeeze()
    df_dm.rename("DM", inplace=True)

    return df_dm


def get_stiffness_df(dir_name: Path) -> pd.DataFrame:
    df_mts = get_mts_df(dir_name)
    gb = df_mts.groupby("ID")
    id = list()
    stiffness = list()

    for x in gb.groups:
        df = gb.get_group(x)
        axial_change = df["Axial Displacement"].rolling(window=21).mean().to_numpy()
        axial_change = (axial_change[1:] > axial_change[:-1]).astype(int)
        axial_change = np.insert(axial_change, 0, [0, 0])
        axial_change = pd.Series(axial_change)
        axial_change = (
            axial_change.rolling(window=21)
            .mean()
            .fillna(0)
            .round(0)
            .astype(int)
            .to_numpy()
        )
        axial_change = (axial_change[1:] != axial_change[:-1]).astype(int)
        axial_change = pd.Series(axial_change)
        axial_change = axial_change.cumsum()

        df = df[np.logical_and(axial_change > 7, axial_change < 13)].reset_index(
            drop=True
        )

        x = df["Axial Displacement"]
        y = df["Axial Force"]
        res = stats.linregress(x, y)

        stiffness.append(res.slope)
        id.append(df["ID"][0])

        # title = df['ID'][0]
        # plt.plot(x, y, 'o', label='data' )
        # plt.plot(x, res.intercept + res.slope * x, 'r', label='regression')
        # plt.title = title
        # plt.xlabel = 'Axial Displacement [mm]'
        # plt.ylabel = 'Axial Force [N]'
        # plt.legend()
        # plt.show()

    df_out = pd.DataFrame({"ID": id, "ST": stiffness})
    return df_out
