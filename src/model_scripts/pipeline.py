"""Functions to create pipelines for different model scenarios."""
import datetime
import logging
from pathlib import Path

from src.neural.trainer import get_output_path


def create_run_name(name):
    # create run name by appending time and date
    run_name = name + datetime.datetime.now().strftime("_%Y-%m-%d_%H-%M-%S")
    return run_name


def create_logger(run_name):
    # create filepath for log
    log_file = get_output_path(
        base_name=run_name, file_name=Path(run_name).with_suffix(".log")
    )
    # create file logger
    log_fmt = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    logging.basicConfig(filename=log_file, level=logging.INFO, format=log_fmt)
    # apply settings to root logger, so that loggers in modules can inherit both the file and console logger
    logger = logging.getLogger()
    # add console logger
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    console.setFormatter(logging.Formatter(log_fmt))
    logger.addHandler(console)
