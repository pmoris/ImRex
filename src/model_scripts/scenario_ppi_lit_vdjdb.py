""" Scenario for neural network. """
import src.bacli as bacli
from src.config import PROJECT_ROOT
from src.data.vdjdb_source import VdjdbSource
from src.models.model_ppi_lit_vdjdb import ModelPPILitVDJdb
from src.neural.trainer import Trainer
from src.processing.kfolds import fold_iterator, random_fold_splitter
from src.processing.ppi_lit_generator import ppi_lit_generator  # , PPILitGenerator2
from src.processing.splitter import splitter


bacli.set_description(__doc__)


@bacli.command
def run(
    batch_size: int = 512,
    val_split: float = None,
    epochs: int = 40,
    neg_ratio: float = 0.5,
    min1: int = 10,
    max1: int = 20,
    min2: int = 8,
    max2: int = 13,
    name: str = "",
    n_folds: int = 3,
    early_stop=False,
    data_path=PROJECT_ROOT / "data/interim/vdjdb-human-no10x.csv",
):

    ppi_source_pos = VdjdbSource(data_path)

    pep1_range = (min1, max1)
    pep2_range = (min2, max2)

    trainer = Trainer(epochs, include_early_stop=early_stop)
    model = ModelPPILitVDJdb(max1, max2, nameSuffix=name)

    if val_split is not None:
        train, val = splitter(ppi_source_pos, ratio=val_split)
        iterations = [(train, val)]
    else:
        folds = random_fold_splitter(ppi_source_pos, n_folds)
        iterations = fold_iterator(folds)

    for index, (train, val) in enumerate(iterations):
        print("Iteration:", index)
        print("batch size", batch_size)

        print("train set", len(train))
        print("val set", len(val))

        train_stream = ppi_lit_generator(
            train, neg_ratio, batch_size, pep1_range, pep2_range, symmetric=False
        )
        val_stream = ppi_lit_generator(
            val, neg_ratio, batch_size, pep1_range, pep2_range, symmetric=False
        )

        trainer.train(model, train_stream, val_stream, iteration=index)
