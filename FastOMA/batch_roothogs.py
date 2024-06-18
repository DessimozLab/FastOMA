
import shutil
from pathlib import Path
from ._wrappers import logger
from . import __version__ as fastoma_version

big_rhog_filesize_thresh = 400 * 1000
sum_list_rhogs_filesize_thresh = 2 * 1e6


"""

fastoma-batch-roothogs --input-roothogs omamer_rhogs --out-big rhogs_big  --out-rest rhogs_rest -vv

"""

class BatchBuilder:
    def __init__(self, outdir: Path, max_size: int):
        self.outdir = outdir
        self.max_size = max_size

    def __enter__(self):
        self.cur_batch = []
        self.cur_size = 0
        self.counter = 0
        self.outdir.mkdir(parents=True, exist_ok=True)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if len(self.cur_batch) > 0:
            self._flush()

    def add_hog(self, hog_file: Path):
        self.cur_batch.append(hog_file)
        self.cur_size += hog_file.stat().st_size
        logger.debug("adding %s with size %d to batch %d", hog_file, hog_file.stat().st_size, self.counter)
        if self.cur_size > self.max_size:
            self._flush()
            self.counter += 1

    def _flush(self):
        batch_dir = self.outdir / str(self.counter)
        batch_dir.mkdir()
        for fn in self.cur_batch:
            shutil.copy(fn, batch_dir)
        logger.debug("creating batch %s with %d families; total size of files is %d",
                     batch_dir, len(self.cur_batch), self.cur_size)
        self.cur_size = 0
        self.cur_batch = []


def folder_1h_rhog(roothog_path: Path, output_folder_big: Path, output_folder_rest: Path):
    # create a list of hogs in descending filesize order
    hog_size_tuples = sorted([(f, f.stat().st_size) for f in roothog_path.rglob("*.fa")], key=lambda x: -x[1])
    with BatchBuilder(output_folder_big, 1) as big_hogs, \
            BatchBuilder(output_folder_rest, sum_list_rhogs_filesize_thresh) as rest_hogs:
        for hog, fsize in hog_size_tuples:
            if fsize > big_rhog_filesize_thresh:
                big_hogs.add_hog(hog)
            else:
                rest_hogs.add_hog(hog)


def fastoma_batch_roothogs():
    import argparse
    parser = argparse.ArgumentParser(description="Analyse roothog families and create batches for analysis")
    parser.add_argument("--version", action="version", version="FastOMA v"+fastoma_version)
    parser.add_argument('--input-roothogs', required=True, help="folder where input roothogs are stored")
    parser.add_argument('--out-big', required=True, help="folder where the big single family hogs should be stored")
    parser.add_argument('--out-rest', required=True, help="folder where the remaining families should be stored in"
                                                          "batch subfolder structure.")
    parser.add_argument('-v', default=0, action="count", help="incrase verbosity")
    conf_batch_roothogs = parser.parse_args()
    logger.setLevel(level=30 - 10 * min(conf_batch_roothogs.v, 2))
    logger.debug("Arguments: %s", conf_batch_roothogs)

    folder_1h_rhog(Path(conf_batch_roothogs.input_roothogs), Path(conf_batch_roothogs.out_big), Path(conf_batch_roothogs.out_rest))

