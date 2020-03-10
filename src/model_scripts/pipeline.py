"""Functions to create pipelines for different model scenarios."""
import datetime
import logging

from src.config import LOG_DIR


def create_run_name(name):
    # create run name by appending time and date
    run_name = name + datetime.datetime.now().strftime("_%Y-%m-%d_%H-%M-%S")
    return run_name


def create_logger(run_name):
    # create filepath for log
    log_file = LOG_DIR / run_name
    log_file = log_file.with_suffix(".log")
    log_file.parent.mkdir(parents=True, exist_ok=True)
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
