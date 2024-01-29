from __future__ import annotations

from pathlib import Path
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

implant_type = [
    "Nevo",
    "Element",
    "Nevo",
    "Element",
    "Nevo",
    "Element",
    "Nevo",
    "Element",
    "Nevo",
    "Element",
    "Nevo",
    "Element",
    "Nevo",
    "Element",
    "Nevo",
    "Element",
    "Nevo",
    "Element",
    "Nevo",
    "Element",
    "Nevo",
    "Element",
    "Nevo",
    "Element",
    "Nevo",
    "Element",
    "Nevo",
    "Element",
    "Nevo",
    "Element",
]

homogenous = [
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
]

sample_no = [
    "T2_10",
    "T2_05",
    "T2_09",
    "T2_04",
    "T2_08",
    "T2_03",
    "T2_07",
    "T2_02",
    "T2_06",
    "T2_01",
    "T1_06",
    "T1_05",
    "T1_04",
    "T1_03",
    "T1_02",
    "T1_01",
    "T7_01",
    "T7_02",
    "T7_03",
    "T7_04",
    "T3_01",
    "T3_02",
    "T5_01",
    "T5_02",
    "T5_03",
    "T5_04",
    "T5_05",
    "T5_06",
    "T5_07",
    "F1_02",
]


def strip_versioning(file: Path) -> Path | None:
    if not file.exists():
        return None

    if file.is_dir():
        for inner_file in file.iterdir():
            strip_versioning(inner_file)

    name_split = file.name.split(";")
    file = file.rename(file.parent / name_split[0])
    name_split = file.name.split(".DIR")
    file = file.rename(file.parent / name_split[0])

    return file


def read_results_bone_morpho(file_name: Path) -> dict:
    with file_name.open("r") as f:
        names = f.readline()
        values = f.readline()

    names = names.split("\t")
    values = values.split("\t")

    val = {a.strip(): b.strip() for (a, b) in zip(names, values)}

    return val


def results_bone_morpho(data_dir: Path) -> dict | None:
    if not data_dir.exists():
        return None
    results_bone_morpho = list(data_dir.glob("*3DRESULTS_BONE_MORPHO.TXT"))
    if len(results_bone_morpho) != 1:
        return None

    return read_results_bone_morpho(results_bone_morpho[0])


def main():
    strip_versioning(Path(__file__).parent / "data")

    vals = list()
    bv_tv = dict()
    for data_dir in Path("data").iterdir():
        val = results_bone_morpho(data_dir)
        vals.append(val)

    df = pd.DataFrame(vals, columns=["MeasNo", "TRI-BV/TV"])
    df["SampleNo"] = sample_no
    df["Homogenous"] = homogenous
    df["TRI-BV/TV"] = df["TRI-BV/TV"].astype(float)
    df = df.sort_values(by=["TRI-BV/TV", "Homogenous"])
    df = df.sort_values(by=["Homogenous", "TRI-BV/TV"])
    df["ImplantType"] = implant_type
    print(df)

    cmap = mpl.colors.ListedColormap([(1, 0, 0, 1), (0, 1, 0, 1)])
    fig, ax = plt.subplots(figsize=(4, 6))
    sc = ax.scatter(x=df["ImplantType"], y=df["TRI-BV/TV"], c=df["Homogenous"], cmap=cmap)

    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_ticks([0.25, 0.75])
    cbar.set_ticklabels(['not homogeneous', 'homogeneous'])

    ax.margins(x=0.4)
    plt.tight_layout()
    plt.show()

    print("stop")


if __name__ == "__main__":
    main()
