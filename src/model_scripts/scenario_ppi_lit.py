""" Scenario for neural network. """
import src.bacli as bacli
from src.config import PROJECT_ROOT
from src.data.ppi_source import PpiSource, SequencesMap
from src.models.model_ppi_lit import ModelPPILit
from src.neural.trainer import Trainer
from src.processing.kfolds import fold_iterator, random_fold_splitter
from src.processing.ppi_lit_generator import ppi_lit_generator, ppi_lit_generator2
from src.processing.splitter import splitter

bacli.set_description(__doc__)


@bacli.command
def run(
    batch_size: int = 512,
    val_split: float = None,
    epochs: int = 40,
    neg_ratio: float = 0.5,
    min_len: int = 0,
    max_len: int = 1200,
    name: str = "",
    n_folds: int = 3,
    early_stop=False,
    data_path=PROJECT_ROOT / "data/raw/ppi/PPI_positive.csv",
    negative_path=PROJECT_ROOT / "data/raw/ppi/data/PPI_negative.csv",
    sequences_path=PROJECT_ROOT / "data/raw/ppi/PPI_sequences.csv",
    swap=False,
):

    sequences_map = SequencesMap(sequences_path)
    ppi_source_pos = PpiSource(data_path, sequences_map, label=1)
    pep_range = (min_len, max_len)

    print("swap:", swap)

    trainer = Trainer(epochs, include_early_stop=early_stop)
    model = ModelPPILit(max_len, max_len, nameSuffix=name)

    if val_split is not None:
        train, val = splitter(ppi_source_pos, ratio=val_split)
        if negative_path:
            ppi_source_neg = PpiSource(negative_path, sequences_map, label=0)
            neg_train, neg_val = splitter(ppi_source_neg, ratio=val_split)
            iterations = [((train, neg_train), (val, neg_val))]
        else:
            iterations = [(train, val)]
    else:
        folds = random_fold_splitter(ppi_source_pos, n_folds)
        iterations = fold_iterator(folds)
        if negative_path:
            ppi_source_neg = PpiSource(negative_path, sequences_map, label=0)
            neg_folds = random_fold_splitter(ppi_source_neg, n_folds)
            neg_iterations = fold_iterator(neg_folds)
            iterations = [
                ((train, neg_train), (val, neg_val))
                for (train, val), (neg_train, neg_val) in zip(
                    iterations, neg_iterations
                )
            ]

    for index, (train, val) in enumerate(iterations):
        print("Iteration:", index)
        print("batch size", batch_size)
        if negative_path:
            train, neg_train = train
            val, neg_val = val
            print("train set", len(train))
            print("val set", len(val))

            print("neg train set", len(neg_train))
            print("neg val set", len(neg_val))

            train_stream = ppi_lit_generator2(
                train, neg_train, neg_ratio, batch_size, pep_range, pep_range, swap=swap
            )
            val_stream = ppi_lit_generator2(
                val, neg_val, neg_ratio, batch_size, pep_range, pep_range, swap=swap
            )

        else:
            print("train set", len(train))
            print("val set", len(val))

            train_stream = ppi_lit_generator(
                train, neg_ratio, batch_size, pep_range, pep_range
            )
            val_stream = ppi_lit_generator(
                val, neg_ratio, batch_size, pep_range, pep_range
            )

        trainer.train(model, train_stream, val_stream, iteration=index)


# if __name__ == "__main__":
#     run(val_split=0.2, batch_size=4, data_path="../data/PPI_positive_SUBSET.csv")
#     exit()
