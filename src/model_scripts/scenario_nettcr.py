""" Scenario for neural network. """
import src.bacli as bacli
from src.config import PROJECT_ROOT
from src.data.vdjdb_source import VdjdbSource
from src.models.model_nettcr import ModelNetTCR
from src.neural.trainer import Trainer
from src.processing.kfolds import (
    epitope_stratified_fold_splitter,
    fold_iterator,
    random_fold_splitter,
)
from src.processing.net_tcr_batch_generator import nettcr_batch_generator
from src.processing.splitter import splitter

bacli.set_description(__doc__)


@bacli.command
def run(
    batch_size: int = 128,
    val_split: float = None,
    epochs: int = 40,
    neg_ratio: float = 0.5,
    min_group: int = 32,
    name: str = "",
    n_folds: int = 5,
    early_stop=False,
    data_path=PROJECT_ROOT
    / "data/interim/vdjdb-2019-08-08/vdjdb-human-tra-trb-no10x.csv",
    stratified: bool = False,
):

    data_source = VdjdbSource(filepath=data_path)

    trainer = Trainer(epochs, include_early_stop=early_stop)
    model = ModelNetTCR(name_suffix=name)

    if val_split is not None:
        train, val = splitter(data_source, ratio=val_split)
        iterations = [(train, val)]
    else:
        fold_splitter = (
            epitope_grouped_fold_splitter if stratified else random_fold_splitter
        )
        folds = fold_splitter(data_source, n_folds)
        iterations = fold_iterator(folds)
    for index, (train, val) in enumerate(iterations):
        print("Iteration:", index)
        print("train set", len(train))
        print("val set", len(val))
        print("batch size", batch_size)

        train_stream = nettcr_batch_generator(train, neg_ratio, batch_size, min_group)
        val_stream = nettcr_batch_generator(val, neg_ratio, batch_size, min_group)

        trainer.train(model, train_stream, val_stream, iteration=index)
