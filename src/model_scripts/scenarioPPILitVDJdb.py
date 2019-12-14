""" Scenario for neural network. """
import src.bacli as bacli
from src.config import PROJECT_ROOT
from src.data.vdjdbSource import VdjdbSource
from src.models.modelPPILitVDJdb import ModelPPILitVDJdb
from src.neural.trainer import Trainer
from src.processing.kfolds import FoldIterator, RandomFoldSplitter
from src.processing.ppi_lit_generator import PPILitGenerator  # , PPILitGenerator2
from src.processing.splitter import Splitter


bacli.setDescription(__doc__)


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
    nrFolds: int = 3,
    early_stop=False,
    data_path=PROJECT_ROOT / "data/interim/vdjdb-human-no10x.csv",
):

    ppiSourcePos = VdjdbSource(data_path)

    pep1Range = (min1, max1)
    pep2Range = (min2, max2)

    trainer = Trainer(epochs, includeEarlyStop=early_stop)
    model = ModelPPILitVDJdb(max1, max2, nameSuffix=name)

    if val_split is not None:
        train, val = Splitter(ppiSourcePos, ratio=val_split)
        iterations = [(train, val)]
    else:
        folds = RandomFoldSplitter(ppiSourcePos, nrFolds)
        iterations = FoldIterator(folds)

    for index, (train, val) in enumerate(iterations):
        print("Iteration:", index)
        print("batch size", batch_size)

        print("train set", len(train))
        print("val set", len(val))

        trainStream = PPILitGenerator(
            train, neg_ratio, batch_size, pep1Range, pep2Range, symmetric=False
        )
        valStream = PPILitGenerator(
            val, neg_ratio, batch_size, pep1Range, pep2Range, symmetric=False
        )

        trainer.train(model, trainStream, valStream, iteration=index)
