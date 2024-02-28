from __future__ import annotations

from pathlib import Path
import pandas as pd


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
    df_morpho.rename(columns={"VOX-BV/TV": "BVTV", "MeasNo": "MeasNoHighRes"}, inplace=True)
    df_morpho["BVTV"] = df_morpho["BVTV"].astype(float)
    df_morpho["MeasNoHighRes"] = df_morpho["MeasNoHighRes"].astype(int)
    df_morpho = df_morpho.sort_values(by=["MeasNoHighRes"])

    return df_morpho


def get_sample_df(dir_name: Path) -> pd.DataFrame:
    df_samples = pd.read_excel(dir_name / "samples.xlsx")
    df_samples = df_samples[['MeasNoHighRes', 'ID', "group", "selected"]]

    return df_samples


def get_it_df(dir_name: Path) -> pd.DataFrame:
    dfs_it = list()
    for file_name in dir_name.glob("*.csv"):
        dfs_it.append(pd.read_csv(file_name, delimiter=";"))

    df_it = pd.concat(dfs_it)
    df_it.rename(columns={"Referenz": "ID", "Max. erreichtes Drehmoment": "IT"}, inplace=True)
    df_it = df_it.groupby("ID", as_index=False).first()

    df_it = df_it[['IT', 'ID']]

    return df_it


def get_uf_df(dir_name: Path) -> pd.DataFrame:
    dfs_uf = list()
    for file_name in dir_name.glob("*.csv"):
        df_uf = pd.read_csv(file_name, header=1)
        df_uf = df_uf[1:]  # Remove the first row to reindex the DataFrame
        df_uf.reset_index(drop=True, inplace=True)  # Reset index if desired
        df_uf["ID"] = file_name.stem
        dfs_uf.append(df_uf)

    df_uf = pd.concat(dfs_uf)
    df_uf.columns = [f.strip() for f in df_uf.columns]

    df_uf["Time"] = df_uf["Time"].astype(float)
    df_uf["Axial Count"] = df_uf["Axial Count"].astype(float)
    df_uf["Axial Force"] = -df_uf["Axial Force"].astype(float)
    df_uf["Axial Displacement"] = df_uf["Axial Displacement"].astype(float)

    df_uf = df_uf.groupby("ID")["Axial Force"].max()
    df_uf.rename("UF", inplace=True)

    return df_uf


def main():
    strip_versioning(Path(__file__).parent / "data")

    df_morpho = get_bone_morpho_df(Path("data/bone_morpho"))
    df_samples = get_sample_df(Path("data"))
    df_it = get_it_df(Path("data/it_ichiro"))
    df_uf = get_uf_df(Path("data/uf"))

    df = pd.merge(df_samples, df_morpho, on='MeasNoHighRes', how='outer')
    df = pd.merge(df, df_it, on='ID', how='outer')
    df = pd.merge(df, df_uf, on='ID', how='outer')

    # remove unused samples todo:
    df = df[df["selected"] == "yes"]

    df = df.sort_values(by=["group", "BVTV"])

    df.to_csv(Path("data/df_selected.csv"), index=False)

    print("stop")


if __name__ == "__main__":
    main()
