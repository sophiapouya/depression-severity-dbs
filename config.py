import json
from pathlib import Path

CONFIG_PATH = Path(__file__).parent / "config.json"

with open(CONFIG_PATH) as f:
    _config = json.load(f)

PATHS = _config["paths"]
BASE_DIR_SEEG = Path(PATHS["base_dir_seeg"])
BASE_DIR_DBS = Path(PATHS["base_dir_dbs"])
FEATURES_DBS = Path(PATHS["dbs_features_csv"])
FEATURES_SEEG = Path(PATHS["seeg_features_csv"])
FEATURES_SEEG_GM = Path(PATHS["seeg_features_gm_csv"])
ROOT_DIR = Path(PATHS["root_dir"])
CATDI_ELECTRODES = Path(PATHS["catdi_electrodes"])
CATDI_SCORES = Path(PATHS["catdi_scores"])