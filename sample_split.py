from __future__ import annotations

from pathlib import Path
import pandas as pd
pd.options.plotting.backend = "plotly"

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
    # "Nevo",
    # "Element",
]

homogenous = [
    3,
    3,
    3,
    2,
    3,
    3,
    2,
    2,
    2,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    2,
    1,
    2,
    2,
    3,
    2,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
]

# homogenous = [
#     "no",
#     "no",
#     "no",
#     "no",
#     "no",
#     "no",
#     "no",
#     "no",
#     "no",
#     "no",
#     "yes",
#     "yes",
#     "yes",
#     "yes",
#     "yes",
#     "yes",
#     "yes",
#     "yes",
#     "yes",
#     "yes",
#     "yes",
#     "yes",
#     "yes",
#     "yes",
#     "yes",
#     "yes",
#     "yes",
#     "yes",
#     "yes",
#     "yes",
# ]

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


def results_bone_morpho(data_dir: Path, version: int = 1) -> dict | None:
    if not data_dir.exists():
        return None
    results_bone_morpho = list(data_dir.glob(f"*3DRESULTS_BONE_MORPHO_1.TXT"))
    if len(results_bone_morpho) != 1:
        return None

    return read_results_bone_morpho(results_bone_morpho[0], version)


def main():
    strip_versioning(Path(__file__).parent / "data")

    vals = list()
    bv_tv = dict()
    for data_dir in Path("data").iterdir():
        val = results_bone_morpho(data_dir, 2)
        vals.append(val)

    df = pd.DataFrame(vals, columns=["MeasNo", "VOX-BV/TV"])
    df["VOX-BV/TV"] = df["VOX-BV/TV"].astype(float)
    df["MeasNo"] = df["MeasNo"].astype(int)
    df = df.sort_values(by=["MeasNo"])

    df["SampleNo"] = sample_no
    df["Homogenous"] = homogenous

    df["Homogenous"] = df["Homogenous"] != 3
    df["Homogenous"] = df["Homogenous"].astype(str)

    # remove unused samples:
    df.drop(df[df["SampleNo"] == "F1_02"].index, inplace=True)
    df.drop(df[df["SampleNo"] == "T7_03"].index, inplace=True)

    df = df.sort_values(by=["VOX-BV/TV"])
    df["ImplantType"] = implant_type

    fig = df.plot.scatter(x="ImplantType", y="VOX-BV/TV", color="Homogenous", text="SampleNo")
    fig.update_xaxes(type='category')
    fig.show()

    df = df.sort_values(by=["ImplantType", "VOX-BV/TV"])
    print(df)

    print("stop")


if __name__ == "__main__":
    main()
