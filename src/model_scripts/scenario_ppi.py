""" Scenario for neural network. """
import src.bacli as bacli
from src.bio.feature_builder import CombinedPeptideFeatureBuilder
from src.bio.peptide_feature import parse_features, parse_operator
from src.config import PROJECT_ROOT
from src.data.ppi_source import PpiSource, SequencesMap
from src.models.model_ppi_padded import ModelPPIPadded
from src.neural.trainer import Trainer
from src.processing.kfolds import fold_iterator, random_fold_splitter
from src.processing.padded_batch_generator import (
    padded_batch_generator,
    padded_batch_generator2,
)
from src.processing.splitter import splitter


bacli.set_description(__doc__)


@bacli.command
def run(
    batch_size: int = 16,
    val_split: float = None,
    epochs: int = 40,
    neg_ratio: float = 0.5,
    min_len: int = 0,
    max_len: int = 600,
    name: str = "",
    n_folds: int = 3,
    features: str = "hydrophob,isoelectric,mass,hydrophil,charge",
    operator: str = "best",  # can be: prod, diff, layer or best
    early_stop=False,
    data_path=PROJECT_ROOT / "data/raw/ppi/PPI_positive.csv",
    negative_path=PROJECT_ROOT / "data/raw/ppi/data/PPI_negative.csv",
    sequences_path=PROJECT_ROOT / "data/raw/ppi/PPI_sequences.csv",
    swap=False,
):

    features_list = parse_features(features)
    operator = parse_operator(operator)
    feature_builder = CombinedPeptideFeatureBuilder(features_list, operator)

    sequences_map = SequencesMap(sequences_path)
    ppi_source_pos = PpiSource(data_path, sequences_map, label=1)

    print("features:", features_list)
    print("operator:", operator)
    print("swap:", swap)

    pep_range = (min_len, max_len)

    trainer = Trainer(epochs, include_early_stop=early_stop)
    model = ModelPPIPadded(
        max_len, max_len, nameSuffix=name, channels=feature_builder.get_number_layers()
    )

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

            train_stream = padded_batch_generator2(
                train,
                neg_train,
                feature_builder,
                neg_ratio,
                batch_size,
                pep_range,
                pep_range,
                swap=swap,
            )
            val_stream = padded_batch_generator2(
                val,
                neg_val,
                feature_builder,
                neg_ratio,
                batch_size,
                pep_range,
                pep_range,
                swap=swap,
            )

        else:
            print("train set", len(train))
            print("val set", len(val))

            train_stream = padded_batch_generator(
                train,
                feature_builder,
                neg_ratio,
                batch_size,
                pep_range,
                pep_range,
                cache_images=False,
                swap=swap,
            )
            val_stream = padded_batch_generator(
                val,
                feature_builder,
                neg_ratio,
                batch_size,
                pep_range,
                pep_range,
                cache_images=False,
                swap=swap,
            )

        trainer.train(model, train_stream, val_stream, iteration=index)


# if __name__ == "__main__":
#     run(val_split=0.2, batch_size=128, features="mass,charge,hydrophob", operator="best", data_path="../data/PPI_positive_SUBSET.csv")
#     exit()
