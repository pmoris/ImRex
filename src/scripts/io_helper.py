"""Functions to create pipelines for different model scenarios."""
import datetime
import logging
import os
from pathlib import Path

from src.neural.trainer import get_output_path


def create_run_name(name):
    # create run name by appending time and date
    run_name = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S_") + name
    return run_name


def create_evaluate_path(run_name, evaluate_model):
    """ Create output directory for evaluation in same directory as the supplied model file. """
    model_path = Path(evaluate_model).absolute().parent
    evaluate_dir = model_path / ("evaluate_" + run_name)
    evaluate_dir.mkdir(parents=True, exist_ok=False)
    return evaluate_dir


def create_logger(run_name, evaluate_dir=None, log_to_file=True, level=logging.INFO):
    if log_to_file:
        # create filepath for log
        if evaluate_dir:
            # for an evaluation run, the output directory is created separately in the scenario script
            log_file = (evaluate_dir / run_name).with_suffix(".log")
        else:
            log_file = get_output_path(
                base_name=run_name, file_name=Path(run_name).with_suffix(".log")
            )
        # create file logger
        log_fmt = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        logging.basicConfig(filename=log_file, level=level, format=log_fmt)
        # apply settings to root logger, so that loggers in modules can inherit both the file and console logger
        logger = logging.getLogger()
        # add console logger
        console = logging.StreamHandler()
        console.setLevel(level)
        console.setFormatter(logging.Formatter(log_fmt))
        logger.addHandler(console)
    else:
        log_fmt = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        logging.basicConfig(level=level, format=log_fmt)

    # suppress tf logging
    os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2"  # ERROR
    logging.getLogger("tensorflow").setLevel(logging.ERROR)
